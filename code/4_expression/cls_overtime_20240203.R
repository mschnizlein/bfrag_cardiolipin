library(ggplot2)
library(reshape2)
library(growthcurver)
library(ggtext)
library(gridtext)
library(data.table)
# library(glue)
library(mdthemes)
library(matrixStats)
library(readxl)
# library(data.table)
# library(genefilter)


substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
substrLeft <- function(x, n){
  substr(x, 1, n)
}

'%!in%'<-function(x,y)!('%in%'(x,y))

my_theme<-theme_bw() + theme(axis.line=element_line(color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(), text=element_text(size = 10))
my_theme2<-theme_bw() + theme(axis.line=element_line(color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(), text=element_text(size = 10),axis.text.x = element_text(angle=30,hjust = 1))

insert_minor <- function(major_labs, n_minor){
  labs <- c(sapply(major_labs, function(x) c(x, rep("", n_minor))))
  labs[1:(length(labs)-n_minor)]
}

setwd("C:/Users/mksch/Dropbox/matt_data/crosson_lab_msu_schnizlein/bfrag_bile_acid_stress/cls_expression/rtqpcr/data/")

# Read in data
rt.data<-read_xlsx(path="cls_rtqpcr_data_20240203_overtime.xlsx",sheet="Results",col_names = TRUE)
colnames(rt.data)<-gsub(" ","_",as.character(rt.data[47,]))
rt.data<-rt.data[48:nrow(rt.data),]

rt.data$rowID<-substrLeft(rt.data$Well_Position,1)
rt.data$colID<-gsub("[A-Z]","",rt.data$Well_Position)

rt.data.sub<-rt.data[rt.data$rowID %in% c("A","B","C","D","E","F") & rt.data$colID %in% c(1,2,3,4,5,6,7,8,9,10,11,12),]

rt.data.sub$gene<-NA
rt.data.sub[rt.data.sub$colID %in% c(1,2,3,4,5,6) & rt.data.sub$rowID %in% c("A","B","C"),]$gene<-"g2600"
rt.data.sub[rt.data.sub$colID %in% c(7,8,9,10,11,12) & rt.data.sub$rowID %in% c("A","B","C"),]$gene<-"clsA"
rt.data.sub[rt.data.sub$colID %in% c(7,8,9,10,11,12) &  rt.data.sub$rowID %in% c("D","E","F"),]$gene<-"clsB"

rt.data.sub[rt.data.sub$colID %in% c(1,2) &  rt.data.sub$rowID %in% c("D","E","F"),]$gene<-"clsB"
rt.data.sub[rt.data.sub$colID %in% c(3,4) &  rt.data.sub$rowID %in% c("D","E","F"),]$gene<-"clsA"
rt.data.sub[rt.data.sub$colID %in% c(5,6) &  rt.data.sub$rowID %in% c("D","E","F"),]$gene<-"g2600"

rt.data.sub$replicate<-NA
rt.data.sub[rt.data.sub$rowID %in% c("A","D"),]$replicate<-"WT1"
rt.data.sub[rt.data.sub$rowID %in% c("B","E"),]$replicate<-"WT2"
rt.data.sub[rt.data.sub$rowID %in% c("C","F"),]$replicate<-"WT3"

rt.data.sub$condition<-NA
rt.data.sub[rt.data.sub$rowID %in% c("A","B","C") & rt.data.sub$colID %in% c(1,2,3),]$condition<-"T1"
rt.data.sub[rt.data.sub$rowID %in% c("A","B","C") & rt.data.sub$colID %in% c(4,5,6),]$condition<-"T3"
rt.data.sub[rt.data.sub$colID %in% c(7,8,9),]$condition<-"T1"
rt.data.sub[rt.data.sub$colID %in% c(10,11,12),]$condition<-"T3"

rt.data.sub[rt.data.sub$rowID %in% c("D","E","F") & rt.data.sub$colID %in% c(1,2,3,4,5,6),]$condition<-"T3_nort"

rt.data.f<-rt.data.sub[,c("gene","Well_Position","colID","rowID","CT","condition","replicate")]

rt.data.f$CT<-as.numeric(rt.data.f$CT)

# rt.data.f[rt.data.f$CT >=20,]
# rt.data.f<-rt.data.f[rt.data.f$Well_Position %!in% c("C1","D1","D2"),]
# rt.data.f<-rt.data.f[rt.data.f$CT <=22,]
ggplot(rt.data.f,aes(x=condition,y=CT,color=gene))+geom_point()
ggplot(rt.data.f,aes(x=condition,y=CT,color=gene))+geom_point()

# delta delta Ct = delta Ct (treated sample) - delta Ct (untreated sample)
# delta Ct = Ct (gene of interest) - Ct (housekeeping gene)

rt.data.f.agg<-aggregate(CT~gene+condition+replicate,data=rt.data.f,FUN=mean)
rt.data.f.agg.nort<-rt.data.f.agg[rt.data.f.agg$condition %in% "T3_nort",]
rt.data.f.agg<-rt.data.f.agg[rt.data.f.agg$condition %!in% "T3_nort",]
rt.data.f.agg$cond_rep<-paste(rt.data.f.agg$condition,rt.data.f.agg$replicate,sep="_")

norm_data<-rt.data.f.agg[rt.data.f.agg$gene %in% "g2600",]

colnames(norm_data)<-paste0(colnames(norm_data),"_norm")

rt.data_norm<-merge(rt.data.f.agg[rt.data.f.agg$gene %!in% "g2600",],norm_data,by.x="cond_rep",by.y="cond_rep_norm",all.x=TRUE)

# data_f_norm_cond<-dcast(data=as.data.table(rt.data_norm),value.var =  c("CT","CT_norm"),formula= gene+condition~gene_normcondition)

rt.data_norm$delta_Ct<-rt.data_norm$CT - rt.data_norm$CT_norm # CT - control gene

# rt.data_norm$dc_logical<-substrRight(rt.data_norm$condition,2)
#rt.data_norm$sample<-substrLeft(rt.data_norm$condition,4)
rt.data_norm$sample_gene<-paste(rt.data_norm$replicate,rt.data_norm$gene,sep="_")
rt.data_norm.t3<-rt.data_norm[rt.data_norm$condition %in% "T3",c("sample_gene","delta_Ct")]
colnames(rt.data_norm.t3)<-c("sample_gene","delta_Ct_t3")

rt.data_norm.t1<-rt.data_norm[rt.data_norm$condition %in% "T1",]
rt.data_norm.f<-merge(rt.data_norm.t1,rt.data_norm.t3,by.all="sample_gene")
rt.data_norm.f$delta_delta_Ct<-rt.data_norm.f$delta_Ct - rt.data_norm.f$delta_Ct_t3 # calculating delta delta Ct

# Calculate fold gene expression 2^-(ddCt)
rt.data_norm.f$expression_fold_change<-2^-(rt.data_norm.f$delta_delta_Ct)
# rt.data_norm.f$expression_fold_change_norm<-2^-(rt.data_norm.f$delta_Ct_norm)


rt.data_norm.f$gene<-factor(rt.data_norm.f$gene, levels=c("clsA","clsB"))

ddCt.p<-ggplot(rt.data_norm.f,aes(x=gene,y=delta_delta_Ct))+geom_point()+
  as_md_theme(my_theme)+xlab("Gene ID")+ylab("\u0394\u0394 Ct")
plot(ddCt.p)

fc_gene.p<-ggplot(rt.data_norm.f,aes(x=gene,y=log2(expression_fold_change),color=gene))+geom_point(position = position_dodge(width=0.2))+
  scale_y_continuous(limits=c(-3,4),breaks=c(-3,-2,-1,0,1,2,3,4))+
  scale_color_manual(name="Gene Amplified",values=c("#006198","#299764"),limits=c("clsA","clsB"),labels=c("*clsA*","*clsB*"))+
  as_md_theme(my_theme)+theme(legend.position = "none")+
  xlab("Gene Amplified")+ylab("log2(FC) Early Log/Stationary")+ggtitle("")
plot(fc_gene.p)


ggsave(filename="../figures/rtqpcr_cls_overtime_20240203.pdf",plot=fc_gene.p,width=1.5,height=2,units="in",useDingbats=FALSE)

t.test(rt.data_norm.f[rt.data_norm.f$gene %in% "clsA",]$expression_fold_change,rt.data_norm.f[rt.data_norm.f$gene %in% "clsB",]$expression_fold_change)
# p = 0.007221

fc_gene_bar.p<-ggplot(rt.data_norm.f,aes(x=gene,y=log2(expression_fold_change),fill=gene))+geom_bar(stat="summary")+
  scale_y_continuous(limits=c(-2,4),breaks=c(-3,-2,-1,0,1,2,3,4))+
  scale_fill_manual(name="Gene Amplified",values=c("#006198","#299764"),limits=c("clsA","clsB"),labels=c("*clsA*","*clsB*"))+geom_errorbar(stat='summary', width=0.5, linewidth=1)+
  as_md_theme(my_theme)+theme(legend.position = "none")+
  xlab("Gene Amplified")+ylab("log2(FC) Early Log/Stationary")+ggtitle("")
plot(fc_gene_bar.p)

ggsave(filename="../figures/rtqpcr_cls_overtime_20240203_bar.pdf",plot=fc_gene_bar.p,width=1.5,height=2,units="in",device = cairo_pdf)

t.test(rt.data_norm.f[rt.data_norm.f$gene %in% "clsA" & rt.data_norm.f$condition %in% "T1",]$delta_Ct,rt.data_norm.f[rt.data_norm.f$gene %in% "clsA" & rt.data_norm.f$condition %in% "T1",]$delta_Ct_t3)
# p = 0.012
t.test(rt.data_norm.f[rt.data_norm.f$gene %in% "clsB" & rt.data_norm.f$condition %in% "T1",]$delta_Ct,rt.data_norm.f[rt.data_norm.f$gene %in% "clsB" & rt.data_norm.f$condition %in% "T1",]$delta_Ct_t3)
# p = 0.025


clsA.lm<-glm(CT~log2(dilution),rt.data.f[rt.data.f$gene %in% "clsA",],family="gaussian")
summary(clsA.lm)

clsB.lm<-glm(CT~log2(dilution),rt.data.f[rt.data.f$gene %in% "clsB",],family="gaussian")
summary(clsB.lm)

g2600.lm<-glm(CT~log2(dilution),rt.data.f[rt.data.f$gene %in% "g2600",],family="gaussian")
summary(clsA.lm)
