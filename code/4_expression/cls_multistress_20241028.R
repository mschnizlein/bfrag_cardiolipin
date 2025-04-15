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
library(data.table)
library(ggpubr)
library(dplyr)
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
rtdata_p1<-as.data.table(read_xls(path="p207_rtqpcr_20241028_p1.xls",sheet="Results",col_names = TRUE))
rtdata_p2<-as.data.table(read_xls(path="p207_rtqpcr_20241028_p2.xls",sheet="Results",col_names = TRUE))
rtdata_p3<-as.data.table(read_xls(path="p207_rtqpcr_20241028_p3.xls",sheet="Results",col_names = TRUE))
rtdata.list<-list(rtdata_p1,rtdata_p2,rtdata_p3)
names(rtdata.list)<-c("p1","p2","p3")

for (i in 1:3){
  colnames(rtdata.list[[i]])<-gsub(" ","_",as.character(rtdata.list[[i]][47,]))
  rtdata.list[[i]]<-rtdata.list[[i]][48:nrow(rtdata.list[[i]]),]
  rtdata.list[[i]]<-cbind(rtdata.list[[i]],as.data.table(matrix(data=names(rtdata.list[i]),ncol=1,nrow=nrow(rtdata.list[[i]]),dimnames=list(1:nrow(rtdata.list[[i]]),c("plateID")))))
  rtdata.list[[i]]<-rtdata.list[[i]][,-"MTP"]
}

rtdata<-rbindlist(rtdata.list)
rtdata$rowID<-substrLeft(rtdata$Well_Position,1)
rtdata$colID<-gsub("[A-Z]","",rtdata$Well_Position)

# All had 1 uL of RNA

rtdata$gene_prim<-"NA"
rtdata[rtdata$rowID %in% c("C","F","H"),]$gene_prim<-"g2600"
rtdata[rtdata$colID %in% c(1:6) & rtdata$rowID %in% c("G"),]$gene_prim<-"g2600"
rtdata[rtdata$colID %in% c(1:5) & rtdata$rowID %in% c("H"),]$gene_prim<-"g2600"
rtdata[rtdata$rowID %in% c("A","D"),]$gene_prim<-"clsA"
rtdata[rtdata$rowID %in% c("G") & rtdata$colID %in% 7:9,]$gene_prim<-"clsA"
rtdata[rtdata$rowID %in% c("B","E"),]$gene_prim<-"clsB"
rtdata[rtdata$rowID %in% "G" & rtdata$colID %in% 10:12,]$gene_prim<-"clsB"


rtdata[rtdata$plateID %in% "p3" & rtdata$colID %in% 7:12 & rtdata$rowID %in% c("D","E","F","G","H"),]$gene_prim<-"empty"
rtdata[rtdata$plateID %in% "p3" & rtdata$colID %in% 4:6 & rtdata$rowID %in% "G",]$gene_prim<-"empty"
rtdata[rtdata$plateID %in% "p3" & rtdata$colID %in% 1:5 & rtdata$rowID %in% "H",]$gene_prim<-"empty"
rtdata[rtdata$Well_Position %in% c("G6","H5"),]$gene_prim<-"empty"
rtdata<-rtdata[rtdata$gene_prim %!in% "empty",]

# Define negative control primers
rtdata[rtdata$rowID %in% "H" & rtdata$colID %in% 6:7,]$gene_prim<-"clsA"
rtdata[rtdata$rowID %in% "H" & rtdata$colID %in% 6:7,]$gene_prim<-"clsB"

rtdata$replicate<-"NA"
# rtdata[rtdata$colID %in% c(1,4,7,10),]$replicate<-"1"
# rtdata[rtdata$colID %in% c(2,5,8,11),]$replicate<-"2"
# rtdata[rtdata$colID %in% c(3,6,9,12),]$replicate<-"3"


rtdata[rtdata$colID %in% c(1:3,10:12) & rtdata$rowID %in% c("A","B","C"),]$replicate<-"1"
rtdata[rtdata$colID %in% c(7:9) & rtdata$rowID %in% c("D","E","F"),]$replicate<-"1"
rtdata[rtdata$colID %in% c(4:6) & rtdata$rowID %in% c("A","B","C"),]$replicate<-"2"
rtdata[rtdata$colID %in% c(1:3,10:12) & rtdata$rowID %in% c("D","E","F"),]$replicate<-"2"
rtdata[rtdata$colID %in% c(7:9) & rtdata$rowID %in% c("A","B","C"),]$replicate<-"3"
rtdata[rtdata$colID %in% c(4:6) & rtdata$rowID %in% c("D","E","F"),]$replicate<-"3"
rtdata[rtdata$colID %in% c(9:12) & rtdata$rowID %in% c("G","H"),]$replicate<-"3"


rtdata$condition<-"NA"
rtdata[rtdata$plateID %in% c("p1","p2","p3") & rtdata$colID %in% 1:9 & rtdata$rowID %in% c("A","B","C"),]$condition<-"bhis"

# P1 conditions
rtdata[rtdata$plateID %in% "p1" & rtdata$colID %in% 10:12 & rtdata$rowID %in% c("A","B","C"),]$condition<-"K300"
rtdata[rtdata$plateID %in% "p1" & rtdata$colID %in% 1:6 & rtdata$rowID %in% c("D","E","F"),]$condition<-"K300"
rtdata[rtdata$plateID %in% "p1" & rtdata$colID %in% 7:12 & rtdata$rowID %in% c("D","E","F","G"),]$condition<-"SDS"
rtdata[rtdata$plateID %in% "p1" & rtdata$colID %in% 10:12 & rtdata$rowID %in% "H",]$condition<-"SDS"

# P2 conditions
rtdata[rtdata$plateID %in% "p2" & rtdata$colID %in% 10:12 & rtdata$rowID %in% c("A","B","C"),]$condition<-"Monensin"
rtdata[rtdata$plateID %in% "p2" & rtdata$colID %in% 1:6 & rtdata$rowID %in% c("D","E","F"),]$condition<-"Monensin"
rtdata[rtdata$plateID %in% "p2" & rtdata$colID %in% 7:12 & rtdata$rowID %in% c("D","E","F","G"),]$condition<-"Nigericin"
rtdata[rtdata$plateID %in% "p2" & rtdata$colID %in% 10:12 & rtdata$rowID %in% "H",]$condition<-"Nigericin"

# P3 conditions
rtdata[rtdata$plateID %in% "p3" & rtdata$colID %in% 10:12 & rtdata$rowID %in% c("A","B","C"),]$condition<-"deoxycholate"
rtdata[rtdata$plateID %in% "p3" & rtdata$colID %in% 1:6 & rtdata$rowID %in% c("D","E","F"),]$condition<-"deoxycholate"

rtdata$rt_status<-"rt"
rtdata[rtdata$colID %in% 1:6 & rtdata$rowID %in% "G",]$rt_status<-"no_rt"
rtdata[rtdata$colID %in% 1:5 & rtdata$rowID %in% "H",]$rt_status<-"no_rt"

rt.data.f<-rtdata[,c("gene_prim","Well_Position","colID","rowID","CT","replicate","condition","rt_status","plateID")]

rt.data.f$CT<-as.numeric(rt.data.f$CT)
rt.data.f$gene_prim<-factor(rt.data.f$gene_prim)
rt.data.f$condition<-factor(rt.data.f$condition, levels=c("bhis","K300","SDS","Monensin","Nigericin","deoxycholate"))

ggplot(rt.data.f,aes(x=gene_prim,y=CT,color=condition,group=condition))+geom_jitter(position = position_jitterdodge(dodge.width=1,jitter.height=0, jitter.width = 0.2))+facet_wrap(~condition)

ggplot(rt.data.f[rt.data.f$gene_prim %in% "g2600",],aes(x=rt_status,y=CT,color=condition))+geom_jitter(position = position_jitterdodge(dodge.width=1,jitter.height=0, jitter.width = 0.2))

# delta delta Ct = delta Ct (treated sample) - delta Ct (untreated sample)
# delta Ct = Ct (gene of interest) - Ct (housekeeping gene)

# reframe(across(.cols = "CT",.fns=mean), .by=c("gene_prim","condition","replicate","strain"), .data=rt.data.f)
rt.data.f.agg<-reframe(across(.cols = "CT",.fns=mean), .by=c("gene_prim","condition","replicate","rt_status","plateID"), .data=rt.data.f)
rt.data.f.agg.nort<-rt.data.f.agg[rt.data.f.agg$rt_status %in% "no_rt",]
rt.data.f.agg<-rt.data.f.agg[rt.data.f.agg$rt_status %!in% "no_rt",]
rt.data.f.agg$plate_cond_rep<-paste(rt.data.f.agg$condition,rt.data.f.agg$replicate,rt.data.f.agg$plateID,sep="_")

norm_data<-rt.data.f.agg[rt.data.f.agg$gene_prim %in% "g2600",]

colnames(norm_data)<-paste0(colnames(norm_data),"_norm")

rt.data_norm<-merge(rt.data.f.agg[rt.data.f.agg$gene_prim %!in% "g2600",],norm_data,by.x="plate_cond_rep",by.y="plate_cond_rep_norm",all.x=TRUE)
rt.data_norm<-rt.data_norm[!is.na(rt.data_norm$condition),]

# data_f_norm_cond<-dcast(data=as.data.table(rt.data_norm),value.var =  c("CT","CT_norm"),formula= gene+condition~gene_normcondition)

rt.data_norm$delta_Ct<-rt.data_norm$CT - rt.data_norm$CT_norm # CT - control gene

# rt.data_norm$dc_logical<-substrRight(rt.data_norm$condition,2)
#rt.data_norm$sample<-substrLeft(rt.data_norm$condition,4)

rt.data_norm$gene_rep<-paste(rt.data_norm$gene_prim,rt.data_norm$replicate,sep="_")
rt.data_norm.cast<-dcast(data=as.data.table(rt.data_norm), gene_prim+replicate+plateID~condition,value.var="delta_Ct")
rt.data_norm.f<-melt.data.table(as.data.table(rt.data_norm.cast),id.vars = c("gene_prim","replicate","plateID","bhis"),measure.vars = c("K300","SDS","Monensin","Nigericin","deoxycholate"))
rt.data_norm.f<-rt.data_norm.f[!is.na(rt.data_norm.f$value),]

rt.data_norm.f$delta_delta_Ct<-rt.data_norm.f$value - rt.data_norm.f$bhis # calculating delta delta Ct for all conditions

# Calculate fold gene expression 2^-(ddCt)
rt.data_norm.f$expression_fold_change<-2^-(rt.data_norm.f$delta_delta_Ct)

rt.data_norm.f$gene_prim<-factor(rt.data_norm.f$gene_prim, levels=c("clsA","clsB"))

# Plotting
ddCt.allstress.p<-ggplot(rt.data_norm.f,aes(x=gene_prim,y=delta_delta_Ct,color=replicate))+geom_point()+
  as_md_theme(my_theme)+xlab("Gene Amplified")+ylab("\u0394\u0394 Ct (pH 5.5 Ct - pH 7 Ct)")+
  ylim(-4,4)+facet_wrap(~variable,nrow=1)
plot(ddCt.allstress.p)

fc_gene_stat_bar.p<-ggplot(rt.data_norm.f, aes(x=gene_prim,y=log2(expression_fold_change),fill=gene_prim))+
  geom_bar(stat="summary",fun="mean")+stat_summary(geom="errorbar",fun.data=mean_sd, width=0.5, linewidth=0.5,position=position_dodge(width=0.8))+
  scale_y_continuous(name="log2(treated/untreated)",limits=c(-2,4),breaks=c(-2,-1,0,1,2,3,4))+
  scale_x_discrete(name="Gene Amplified",limits=c("clsA","clsB"),labels=c("*clsA*","*clsB*"))+
  as_md_theme(my_theme)+theme(legend.position = "none")+
  ggtitle("")+
  scale_fill_manual(name="Gene Amplified",values=c("#299764","#006198"),limits=c("clsB","clsA"),labels=c("*clsB*","*clsA*"))+facet_wrap(~variable,nrow=1)
plot(fc_gene_stat_bar.p)

ggsave(plot=ddCt.allstress.p,filename="../figures/rtqpcr_cls_KSDSMonNigDc_20241028_ddct.pdf",width=4,height=2,units="in",useDingbats=FALSE)
ggsave(plot=fc_gene_stat_bar.p,filename="../figures/rtqpcr_cls_KSDSMonNigDc_20241028_log2fc.pdf",width=6,height=2,units="in",device=cairo_pdf)

# Between log2(fold-changes) of each gene in a condition
t.test(log2(rt.data_norm.f[rt.data_norm.f$variable %in% "K300" & rt.data_norm.f$gene_prim %in% "clsB",]$expression_fold_change),log2(rt.data_norm.f[rt.data_norm.f$variable %in% "K300" & rt.data_norm.f$gene_prim %in% "clsA",]$expression_fold_change))
# p = 0.00232
t.test(log2(rt.data_norm.f[rt.data_norm.f$variable %in% "SDS" & rt.data_norm.f$gene_prim %in% "clsB",]$expression_fold_change),log2(rt.data_norm.f[rt.data_norm.f$variable %in% "SDS" & rt.data_norm.f$gene_prim %in% "clsA",]$expression_fold_change))
# p = ns
t.test(log2(rt.data_norm.f[rt.data_norm.f$variable %in% "Monensin" & rt.data_norm.f$gene_prim %in% "clsB",]$expression_fold_change),log2(rt.data_norm.f[rt.data_norm.f$variable %in% "Monensin" & rt.data_norm.f$gene_prim %in% "clsA",]$expression_fold_change))
# p = 0.005292
t.test(log2(rt.data_norm.f[rt.data_norm.f$variable %in% "Nigericin" & rt.data_norm.f$gene_prim %in% "clsB",]$expression_fold_change),log2(rt.data_norm.f[rt.data_norm.f$variable %in% "Nigericin" & rt.data_norm.f$gene_prim %in% "clsA",]$expression_fold_change))
# p = 0.0297
t.test(log2(rt.data_norm.f[rt.data_norm.f$variable %in% "deoxycholate" & rt.data_norm.f$gene_prim %in% "clsB",]$expression_fold_change),log2(rt.data_norm.f[rt.data_norm.f$variable %in% "deoxycholate" & rt.data_norm.f$gene_prim %in% "clsA",]$expression_fold_change))
# p = 0.0011

## Between delta Cts of each gene compared to no treatment
rt.data_norm.f_stats<-melt.data.table(as.data.table(rt.data_norm.cast),id.vars = c("gene_prim","replicate","plateID"),measure.vars = c("K300","SDS","Monensin","Nigericin","deoxycholate","bhis"))
rt.data_norm.f_stats<-rt.data_norm.f_stats[!is.na(rt.data_norm.f_stats$value),]
rt.data_norm.f_stats$variable<-factor(rt.data_norm.f_stats$variable,levels=c("bhis","K300","SDS","Monensin","Nigericin","deoxycholate"))
rt.data_norm.f_stats$gene_prim<-factor(rt.data_norm.f_stats$gene_prim, levels=c("clsA","clsB"))

p1_clsA.lm<-glm(formula = value ~ variable,data = rt.data_norm.f_stats[rt.data_norm.f_stats$gene_prim %in% "clsA" & rt.data_norm.f_stats$plateID %in% "p1",])
summary(p1_clsA.lm)
# K300 0.00626
p1_clsB.lm<-glm(formula = value ~ variable,data = rt.data_norm.f_stats[rt.data_norm.f_stats$gene_prim %in% "clsB" & rt.data_norm.f_stats$plateID %in% "p1",])
summary(p1_clsB.lm)
# K300 0.000287

p2_clsA.lm<-glm(formula = value ~ variable,data = rt.data_norm.f_stats[rt.data_norm.f_stats$gene_prim %in% "clsA" & rt.data_norm.f_stats$plateID %in% "p2",])
summary(p2_clsA.lm)
# Monensin 0.031
p2_clsB.lm<-glm(formula = value ~ variable,data = rt.data_norm.f_stats[rt.data_norm.f_stats$gene_prim %in% "clsB" & rt.data_norm.f_stats$plateID %in% "p2",])
summary(p2_clsB.lm)
# Monensin 0.0183

p3dc_clsA.lm<-glm(formula = value ~ variable,data = rt.data_norm.f_stats[rt.data_norm.f_stats$gene_prim %in% "clsA" & rt.data_norm.f_stats$plateID %!in% c("p1","p2"),])
summary(p3dc_clsA.lm)
# Dc 7.88E-5
p3dc_clsB.lm<-glm(formula = value ~ variable,data = rt.data_norm.f_stats[rt.data_norm.f_stats$gene_prim %in% "clsB" & rt.data_norm.f_stats$plateID %!in% c("p1","p2"),])
summary(p3dc_clsB.lm)
# Dc 0.00861
