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
rt.data.str<-read_xlsx(path="p207_rtqpcr_t3cls_20240220.xlsx",sheet="Results",col_names = TRUE)
colnames(rt.data.str)<-gsub(" ","_",as.character(rt.data.str[47,]))
rt.data.str<-rt.data.str[48:nrow(rt.data.str),]

rt.data.str$rowID<-substrLeft(rt.data.str$Well_Position,1)
rt.data.str$colID<-gsub("[A-Z]","",rt.data.str$Well_Position)

# rt.data.str<-rt.data[rt.data$rowID %in% c("A","B","C","D","E","F","G","H") & rt.data$colID %in% c(1,2,3,4,5,6,7,8,9,10,11,12),]

# 2600 primers had 1 uL of RNA, the rest had 0.5 uL of RNA

rt.data.str$gene_prim<-NA
rt.data.str[rt.data.str$colID %in% c(1,2,3,4,5,6,7,8,9) & rt.data.str$rowID %in% c("A","B","C"),]$gene_prim<-"g2600"
rt.data.str[rt.data.str$colID %in% c(1,2,3,4,5,6,7,8,9) & rt.data.str$rowID %in% c("D","E","F"),]$gene_prim<-"clsA"
rt.data.str[rt.data.str$colID %in% c(10,11,12),]$gene_prim<-"clsB"
rt.data.str[rt.data.str$colID %in% c(7,8,9) & rt.data.str$rowID %in% c("G"),]$gene_prim<-"clsB"
rt.data.str[rt.data.str$colID %in% c(1,2,3) & rt.data.str$rowID %in% c("G","H"),]$gene_prim<-"g2600"
rt.data.str[rt.data.str$colID %in% c(4,5,6) & rt.data.str$rowID %in% c("G"),]$gene_prim<-"g2600"
rt.data.str[rt.data.str$colID %in% c(4,7) & rt.data.str$rowID %in% c("H"),]$gene_prim<-"clsA"
rt.data.str[rt.data.str$colID %in% c(5,8) & rt.data.str$rowID %in% c("H"),]$gene_prim<-"clsB"
rt.data.str[rt.data.str$colID %in% c(6,9) & rt.data.str$rowID %in% c("H"),]$gene_prim<-"g2600"

rt.data.str$replicate<-NA
rt.data.str[rt.data.str$colID %in% c(1,4,7,10) & rt.data.str$rowID %in% c("A","B","C","D","E","F","G","H"),]$replicate<-"1"
rt.data.str[rt.data.str$colID %in% c(2,5,8,11) & rt.data.str$rowID %in% c("A","B","C","D","E","F","G","H"),]$replicate<-"2"
rt.data.str[rt.data.str$colID %in% c(3,6,9,12) & rt.data.str$rowID %in% c("A","B","C","D","E","F","G","H"),]$replicate<-"3"

rt.data.str$strain<-NA
rt.data.str[rt.data.str$colID %in% c(1,2,3) & rt.data.str$rowID %in% c("A","B","C","D","E","F","G"),]$strain<-"WT"
rt.data.str[rt.data.str$colID %in% c(10,11,12) & rt.data.str$rowID %in% c("A","B","C"),]$strain<-"WT"

rt.data.str[rt.data.str$colID %in% c(4,5,6) & rt.data.str$rowID %in% c("A","B","C","D","E","F"),]$strain<-"clsA"
rt.data.str[rt.data.str$colID %in% c(10,11,12) & rt.data.str$rowID %in% c("D","E","F"),]$strain<-"clsA"
rt.data.str[rt.data.str$colID %in% c(1,2,3) & rt.data.str$rowID %in% "H",]$strain<-"clsA"
rt.data.str[rt.data.str$colID %in% c(7,8,9) & rt.data.str$rowID %in% c("A","B","C","D","E","F","G"),]$strain<-"clsB"
rt.data.str[rt.data.str$colID %in% c(4,5,6) & rt.data.str$rowID %in% "G",]$strain<-"clsB"
rt.data.str[rt.data.str$colID %in% c(10,11,12) & rt.data.str$rowID %in% c("G","H"),]$strain<-"clsB"

rt.data.str[rt.data.str$colID %in% c(4,5,6,7,8,9) & rt.data.str$rowID %in% "H",]$strain<-"neg"

rt.data.str$condition<-NA
rt.data.str[rt.data.str$rowID %in% c("A","B","C","D","E","F"),]$condition<-"T3"
rt.data.str[rt.data.str$rowID %in% c("G","H") & rt.data.str$colID %in% c(7,8,9,10,11,12),]$condition<-"T3"
rt.data.str[rt.data.str$rowID %in% c("G","H") & rt.data.str$colID %in% c(1,2,3,4,5,6),]$condition<-"T3_nort"

rt.data.str$rna_vol<-0.5
rt.data.str[rt.data.str$rowID %in% c("A","B","C") & rt.data.str$colID %in% c(1,2,3,4,5,6,7,8,9),]$rna_vol<-1

rt.data.f<-rt.data.str[,c("gene_prim","Well_Position","colID","rowID","CT","strain","replicate","condition","rna_vol")]

rt.data.f$CT<-as.numeric(rt.data.f$CT)
rt.data.f$gene_prim<-factor(rt.data.f$gene_prim)
rt.data.f$condition<-factor(rt.data.f$condition, levels=c("T3","T3_nort"))

# Since twice as much RNA was used in g2600 conditions, I will add 0.76 to the Ct value since that is the slope of a line that runs through 2-fold dilutions with that primer pair
rt.data.f[rt.data.f$colID %in% c(1,2,3,4,5,6,7,8,9) & rt.data.f$rowID %in% c("A","B","C"),]$CT<-rt.data.f[rt.data.f$colID %in% c(1,2,3,4,5,6,7,8,9) & rt.data.f$rowID %in% c("A","B","C"),]$CT+0.76

# rt.data.f[rt.data.f$CT >=20,]
# rt.data.f<-rt.data.f[rt.data.f$Well_Position %!in% c("C1","D1","D2"),]
# rt.data.f<-rt.data.f[rt.data.f$CT <=22,]
ggplot(rt.data.f,aes(x=strain,y=CT,color=strain,group=condition))+geom_jitter(position = position_jitterdodge(dodge.width=1,jitter.height=0, jitter.width = 0.2))+facet_wrap(~gene_prim)
ggplot(rt.data.f[rt.data.f$gene_prim %in% "g2600",],aes(x=strain,y=CT,color=condition))+geom_jitter(position = position_jitterdodge(dodge.width=1,jitter.height=0, jitter.width = 0.2))

# delta delta Ct = delta Ct (treated sample) - delta Ct (untreated sample)
# delta Ct = Ct (gene of interest) - Ct (housekeeping gene)

rt.data.f.agg<-aggregate(CT~gene_prim+condition+replicate+strain,data=rt.data.f,FUN=mean)
rt.data.f.agg.nort<-rt.data.f.agg[rt.data.f.agg$condition %in% "T3_nort",]
rt.data.f.agg<-rt.data.f.agg[rt.data.f.agg$condition %!in% "T3_nort",]
rt.data.f.agg$strain_rep<-paste(rt.data.f.agg$strain,rt.data.f.agg$replicate,sep="_")


norm_data<-rt.data.f.agg[rt.data.f.agg$gene_prim %in% "g2600",]

colnames(norm_data)<-paste0(colnames(norm_data),"_norm")


rt.data_norm<-merge(rt.data.f.agg[rt.data.f.agg$gene_prim %!in% "g2600",],norm_data,by.x="strain_rep",by.y="strain_rep_norm",all.x=TRUE)

# data_f_norm_cond<-dcast(data=as.data.table(rt.data_norm),value.var =  c("CT","CT_norm"),formula= gene+condition~gene_normcondition)

rt.data_norm$delta_Ct<-rt.data_norm$CT - rt.data_norm$CT_norm # CT - control gene

# rt.data_norm$dc_logical<-substrRight(rt.data_norm$condition,2)
#rt.data_norm$sample<-substrLeft(rt.data_norm$condition,4)

rt.data_norm$gene_rep<-paste(rt.data_norm$gene_prim,rt.data_norm$replicate,sep="_")
rt.data_norm.clsb<-rt.data_norm[rt.data_norm$strain %in% "clsB",c("gene_rep","delta_Ct")]
colnames(rt.data_norm.clsb)<-c("gene_rep","delta_Ct_clsB")
rt.data_norm.clsa<-rt.data_norm[rt.data_norm$strain %in% "clsA",c("gene_rep","delta_Ct")]
colnames(rt.data_norm.clsa)<-c("gene_rep","delta_Ct_clsA")
rt.data_norm.wt<-rt.data_norm[rt.data_norm$strain %in% "WT",]


rt.data_norm.f<-merge(rt.data_norm.wt,rt.data_norm.clsa,by.all="gene_rep",all.x=TRUE,all.y=TRUE)
rt.data_norm.f<-merge(rt.data_norm.f,rt.data_norm.clsb,by.all="gene_rep",all.x=TRUE,all.y=TRUE)

rt.data_norm.f$delta_delta_Ct_clsa.wt<-rt.data_norm.f$delta_Ct_clsA - rt.data_norm.f$delta_Ct # calculating delta delta Ct for clsa strain
rt.data_norm.f$delta_delta_Ct_clsb.wt<-rt.data_norm.f$delta_Ct_clsB - rt.data_norm.f$delta_Ct # calculating delta delta Ct for clsb strain
rt.data_norm.f$delta_delta_Ct_clsb.a<-rt.data_norm.f$delta_Ct_clsB - rt.data_norm.f$delta_Ct_clsA # calculating delta delta Ct for between clsa and clsb strains

# Calculate fold gene expression 2^-(ddCt)
rt.data_norm.f$expression_fold_change_clsa.wt<-2^-(rt.data_norm.f$delta_delta_Ct_clsa.wt)
rt.data_norm.f$expression_fold_change_clsb.wt<-2^-(rt.data_norm.f$delta_delta_Ct_clsb.wt)
rt.data_norm.f$expression_fold_change_clsb.clsa<-2^-(rt.data_norm.f$delta_delta_Ct_clsb.a)
# rt.data_norm.f$expression_fold_change_norm<-2^-(rt.data_norm.f$delta_Ct_norm)


rt.data_norm.f$gene_prim<-factor(rt.data_norm.f$gene_prim, levels=c("clsA","clsB"))

rt.data_norm.f.m<-melt.data.table(as.data.table(rt.data_norm.f),id.vars = c("gene_rep","strain_rep","gene_prim","condition","replicate","strain","CT","gene_prim_norm","condition_norm","replicate_norm","strain_norm","CT_norm","delta_Ct","delta_Ct_clsA","delta_Ct_clsB"))

rt.data_norm.f.m$strain_comp<-"none"
rt.data_norm.f.m[rt.data_norm.f.m$variable %in% c("delta_delta_Ct_clsa.wt","expression_fold_change_clsa.wt"),]$strain_comp<-"dclsA"
rt.data_norm.f.m[rt.data_norm.f.m$variable %in% c("delta_delta_Ct_clsb.wt","expression_fold_change_clsb.wt"),]$strain_comp<-"dclsB"

ddCt.wt.clsa.p<-ggplot(rt.data_norm.f[!is.na(rt.data_norm.f$delta_delta_Ct_clsa.wt),],aes(x=gene_prim,y=delta_delta_Ct_clsa.wt))+geom_point()+
  as_md_theme(my_theme)+xlab("Gene Amplified")+ylab("\u0394\u0394 Ct (\u0394clsA Ct - WT Ct)")
plot(ddCt.wt.clsa.p)

#ddCt.wt.clsa.small.p<-ggplot(rt.data_norm.f[!is.na(rt.data_norm.f$delta_delta_Ct_wt.clsa),],aes(x=gene_prim,y=delta_delta_Ct_wt.clsa))+geom_point()+
#  ylim(0,5)+
#  as_md_theme(my_theme)+xlab("Gene ID")+ylab("\u0394\u0394 Ct (\u0394clsA Ct - WT Ct)")
# plot(ddCt.wt.clsa.small.p)

ddCt.wt.clsb.p<-ggplot(rt.data_norm.f[!is.na(rt.data_norm.f$delta_delta_Ct_clsb.wt),],aes(x=gene_prim,y=delta_delta_Ct_clsb.wt))+geom_point()+
  as_md_theme(my_theme)+xlab("Gene Amplified")+ylab("\u0394\u0394 Ct (\u0394clsB Ct - WT Ct)")
plot(ddCt.wt.clsb.p)

# ddCt.wt.clsb.small.p<-ggplot(rt.data_norm.f[!is.na(rt.data_norm.f$delta_delta_Ct_wt.clsb),],aes(x=gene_prim,y=delta_delta_Ct_wt.clsb))+geom_point()+
#  ylim(0,5)+
#  as_md_theme(my_theme)+xlab("Gene ID")+ylab("\u0394\u0394 Ct (\u0394clsA Ct - WT Ct)")
# plot(ddCt.wt.clsb.small.p)

fc_gene_stat.p<-ggplot(rt.data_norm.f.m[rt.data_norm.f.m$variable %in% c("expression_fold_change_clsa.wt","expression_fold_change_clsb.wt"),])+
  geom_point(aes(x=strain_comp,y=log2(value),color=gene_prim))+
  geom_point(aes(x=strain_comp,y=log2(value),color=gene_prim))+
  scale_y_continuous(limits=c(-3,4),breaks=c(-3,-2,-1,0,1,2,3,4))+
  scale_x_discrete(limits=c("dclsA","dclsB"),labels=c("\u0394*clsA*","\u0394*clsB*"))+
  as_md_theme(my_theme)+theme(legend.position = "none")+
  xlab("Strain")+ylab("log2(fold-change) (*cls* KO)/WT")+ggtitle("")+
  scale_color_manual(name="Gene Amplified",values=c("#299764","#006198"),limits=c("clsB","clsA"),labels=c("*clsB*","*clsA*"))
plot(fc_gene_stat.p)

ggsave(plot=fc_gene_stat.p,filename="../figures/rtqpcr_cls_straincomp_stationary_20240226.pdf",width=1.5,height=2,units="in",device=cairo_pdf)

fc_gene_stat_bar.p<-ggplot(rt.data_norm.f.m[rt.data_norm.f.m$variable %in% c("expression_fold_change_clsa.wt","expression_fold_change_clsb.wt"),], aes(x=strain_comp,y=log2(value),fill=gene_prim))+
  geom_bar(stat="summary")+geom_errorbar(stat='summary', width=0.5, linewidth=1)+
  scale_y_continuous(limits=c(-2,4),breaks=c(-2,-1,0,1,2,3,4))+
  scale_x_discrete(limits=c("dclsA","dclsB"),labels=c("\u0394*clsA*","\u0394*clsB*"))+
  as_md_theme(my_theme)+theme(legend.position = "none")+
  xlab("Strain")+ylab("log2(fold-change) (*cls*/WT)")+ggtitle("")+
  scale_fill_manual(name="Gene Amplified",values=c("#299764","#006198"),limits=c("clsB","clsA"),labels=c("*clsB*","*clsA*"))
plot(fc_gene_stat_bar.p)

ggsave(plot=fc_gene_stat_bar.p,filename="../figures/rtqpcr_cls_straincomp_stationary_20240226_bar.pdf",width=1.5,height=2,units="in",device=cairo_pdf)


t.test(rt.data_norm.f[rt.data_norm.f$gene_prim %in% "clsA" & rt.data_norm.f$condition %in% "T3",]$delta_Ct,rt.data_norm.f[rt.data_norm.f$gene_prim %in% "clsA" & rt.data_norm.f$condition %in% "T3",]$delta_Ct_clsB)
# p = 0.358
t.test(rt.data_norm.f[rt.data_norm.f$gene_prim %in% "clsB" & rt.data_norm.f$condition %in% "T3",]$delta_Ct,rt.data_norm.f[rt.data_norm.f$gene_prim %in% "clsB" & rt.data_norm.f$condition %in% "T3",]$delta_Ct_clsA)
# p = 0.005









fc_gene_full.p<-ggplot(rt.data_norm.f.m[rt.data_norm.f.m$variable %in% c("expression_fold_change_clsa.wt","expression_fold_change_clsb.wt"),])+
  geom_point(aes(x=strain_comp,y=log2(value),color=gene_prim))+
  geom_point(aes(x=strain_comp,y=log2(value),color=gene_prim))+
  scale_y_continuous()+
  as_md_theme(my_theme)+xlab("Strain")+ylab("log2(fold-change) (*cls* KO)/WT")+ggtitle("Stationary Phase")
plot(fc_gene_full.p)


ggsave(filename="figures/rtqpcr_cls_overtime_20240203.pdf",width=8,height=6,units="in",useDingbats=FALSE)

t.test(log2(rt.data_norm.f.m[rt.data_norm.f.m$variable %in% "expression_fold_change_clsa.wt" & rt.data_norm.f.m$gene_prim %in% "clsB",]$value),log2(rt.data_norm.f.m[rt.data_norm.f.m$variable %in% "expression_fold_change_clsb.wt" & rt.data_norm.f.m$gene_prim %in% "clsA",]$value))
# p = 0.0107



clsA.lm<-glm(CT~log2(dilution),rt.data.f[rt.data.f$gene %in% "clsA",],family="gaussian")
summary(clsA.lm)

clsB.lm<-glm(CT~log2(dilution),rt.data.f[rt.data.f$gene %in% "clsB",],family="gaussian")
summary(clsB.lm)

g2600.lm<-glm(CT~log2(dilution),rt.data.f[rt.data.f$gene %in% "g2600",],family="gaussian")
summary(clsA.lm)
