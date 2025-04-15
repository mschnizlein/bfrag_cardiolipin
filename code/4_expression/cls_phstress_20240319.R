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
rtdata_p1<-read_xlsx(path="p207_rtqpcr_20240318.xlsx",sheet="Results",col_names = TRUE)
rtdata_p2<-read_xls(path="p207_rtqpcr_ph_20241106.xls",sheet="Results",col_names = TRUE)

rtdata.list<-list(rtdata_p1,rtdata_p2)
names(rtdata.list)<-c("p1","p2")

for (i in 1:2){
  colnames(rtdata.list[[i]])<-gsub(" ","_",as.character(rtdata.list[[i]][47,]))
  rtdata.list[[i]]<-rtdata.list[[i]][48:nrow(rtdata.list[[i]]),]
  rtdata.list[[i]]<-cbind(rtdata.list[[i]],as.data.table(matrix(data=names(rtdata.list[i]),ncol=1,nrow=nrow(rtdata.list[[i]]),dimnames=list(1:nrow(rtdata.list[[i]]),c("plateID")))))
}


rtdata<-rbindlist(rtdata.list)
rtdata$rowID<-substrLeft(rtdata$Well_Position,1)
rtdata$colID<-gsub("[A-Z]","",rtdata$Well_Position)

# for P1 2600 primers had 1 uL of RNA, the rest had 0.5 uL of RNA
# for P2 all had the same amount of RNA

# P1 Primer Conditions
rtdata$gene_prim<-"NA"
rtdata[rtdata$colID %in% c(1,2,3,4,5,6) & rtdata$rowID %in% c("A","B","C") & rtdata$plateID %in% "p1",]$gene_prim<-"g2600"
rtdata[rtdata$colID %in% c(7,8,9,10,11,12) & rtdata$rowID %in% c("D","E") & rtdata$plateID %in% "p1",]$gene_prim<-"g2600"
rtdata[rtdata$colID %in% c(7,8,9,10,11,12) & rtdata$rowID %in% c("A","B","C") & rtdata$plateID %in% "p1",]$gene_prim<-"clsA"
rtdata[rtdata$colID %in% c(1,2,3,4,5,6) & rtdata$rowID %in% c("D","E","F") & rtdata$plateID %in% "p1",]$gene_prim<-"clsB"
rtdata[rtdata$colID %in% c(1,2,3,4,5,6) & rtdata$rowID %in% c("G","H") & rtdata$plateID %in% "p1",]$gene_prim<-"g1852"
rtdata[rtdata$colID %in% c(7,8,9,10,11,12) & rtdata$rowID %in% "G" & rtdata$plateID %in% "p1",]$gene_prim<-"g1852"

rtdata$replicate<-"NA"
rtdata[rtdata$colID %in% c(1,3,5,7,9,11) & rtdata$plateID %in% "p1",]$replicate<-"1"
rtdata[rtdata$colID %in% c(2,4,6,8,10,12) & rtdata$plateID %in% "p1",]$replicate<-"2"

rtdata$strain<-"WT"
rtdata[rtdata$colID %in% c(7,8,9,10,11,12) & rtdata$rowID %in% "F" & rtdata$plateID %in% "p1",]$strain<-"neg"

# P2
rtdata[rtdata$rowID %in% c("C","F","H") & rtdata$plateID %in% "p2",]$gene_prim<-"g2600"
rtdata[rtdata$colID %in% c(1:6) & rtdata$rowID %in% c("G") & rtdata$plateID %in% "p2",]$gene_prim<-"g2600"
rtdata[rtdata$colID %in% c(1:5) & rtdata$rowID %in% c("H") & rtdata$plateID %in% "p2",]$gene_prim<-"g2600"
rtdata[rtdata$rowID %in% c("A","D") & rtdata$plateID %in% "p2",]$gene_prim<-"clsA"
rtdata[rtdata$rowID %in% c("G") & rtdata$colID %in% 7:9 & rtdata$plateID %in% "p2",]$gene_prim<-"clsA"
rtdata[rtdata$rowID %in% c("B","E") & rtdata$plateID %in% "p2",]$gene_prim<-"clsB"
rtdata[rtdata$rowID %in% "G" & rtdata$colID %in% 10:12 & rtdata$plateID %in% "p2",]$gene_prim<-"clsB"

# Define negative control primers
rtdata[rtdata$rowID %in% "H" & rtdata$colID %in% 4:5 & rtdata$plateID %in% "p2",]$gene_prim<-"clsA"
rtdata[rtdata$rowID %in% "H" & rtdata$colID %in% 6:7 & rtdata$plateID %in% "p2",]$gene_prim<-"clsB"
rtdata[rtdata$rowID %in% "H" & rtdata$colID %in% 8:9 & rtdata$plateID %in% "p2",]$gene_prim<-"g2600"

rtdata[rtdata$colID %in% c(1:3,10:12) & rtdata$rowID %in% c("A","B","C") & rtdata$plateID %in% "p2",]$replicate<-"3"
rtdata[rtdata$colID %in% c(7:9) & rtdata$rowID %in% c("D","E","F") & rtdata$plateID %in% "p2",]$replicate<-"3"
rtdata[rtdata$colID %in% c(4:6) & rtdata$rowID %in% c("A","B","C") & rtdata$plateID %in% "p2",]$replicate<-"4"
rtdata[rtdata$colID %in% c(1:3,10:12) & rtdata$rowID %in% c("D","E","F") & rtdata$plateID %in% "p2",]$replicate<-"4"
rtdata[rtdata$colID %in% c(7:9) & rtdata$rowID %in% c("A","B","C","G") & rtdata$plateID %in% "p2",]$replicate<-"5"
rtdata[rtdata$colID %in% c(4:6) & rtdata$rowID %in% c("D","E","F") & rtdata$plateID %in% "p2",]$replicate<-"5"
rtdata[rtdata$colID %in% c(10:12) & rtdata$rowID %in% c("G","H") & rtdata$plateID %in% "p2",]$replicate<-"5"

rtdata$condition<-"NA"
rtdata[rtdata$plateID %in% c("p1","p2","p3") & rtdata$colID %in% 1:9 & rtdata$rowID %in% c("A","B","C"),]$condition<-"bhis"

# P1 conditions
rtdata$condition<-"NA"
rtdata[rtdata$colID %in% c(1,2,7,8) & rtdata$plateID %in% "p1",]$condition<-"pH5.5"
rtdata[rtdata$colID %in% c(3,4,9,10) & rtdata$plateID %in% "p1",]$condition<-"pH7"
rtdata[rtdata$colID %in% c(5,6,11,12) & rtdata$plateID %in% "p1",]$condition<-"pH9"

# P2 conditions
rtdata[rtdata$plateID %in% "p2" & rtdata$colID %in% 1:9 & rtdata$rowID %in% c("A","B","C"),]$condition<-"pH7"
rtdata[rtdata$plateID %in% "p2" & rtdata$colID %in% 10:12 & rtdata$rowID %in% c("A","B","C"),]$condition<-"pH5.5"
rtdata[rtdata$plateID %in% "p2" & rtdata$colID %in% 1:6 & rtdata$rowID %in% c("D","E","F"),]$condition<-"pH5.5"
rtdata[rtdata$plateID %in% "p2" & rtdata$colID %in% 7:12 & rtdata$rowID %in% c("D","E","F","G"),]$condition<-"pH9"
rtdata[rtdata$plateID %in% "p2" & rtdata$colID %in% 10:12 & rtdata$rowID %in% "H",]$condition<-"pH9"

rtdata$rt_status<-"rt"
rtdata[rtdata$plateID %in% "p1" & rtdata$colID %in% 7:12 & rtdata$rowID %in% c("D","E"),]$rt_status<-"no_rt"
rtdata[rtdata$plateID %in% "p2" & rtdata$colID %in% 1:6 & rtdata$rowID %in% "G",]$rt_status<-"no_rt"
rtdata[rtdata$plateID %in% "p2" & rtdata$colID %in% 1:5 & rtdata$rowID %in% "H",]$rt_status<-"no_rt"

rt.data.f<-rtdata[,c("gene_prim","Well_Position","colID","rowID","CT","replicate","condition","rt_status","plateID")]

rt.data.f$CT<-as.numeric(rt.data.f$CT)
rt.data.f$gene_prim<-factor(rt.data.f$gene_prim)
rt.data.f$condition<-factor(rt.data.f$condition, levels=c("bhis","pH5.5","pH7","pH9"))

ggplot(rt.data.f,aes(x=gene_prim,y=CT,color=condition,group=condition))+geom_jitter(position = position_jitterdodge(dodge.width=1,jitter.height=0, jitter.width = 0.2))+facet_wrap(~condition)

ggplot(rt.data.f[rt.data.f$gene_prim %in% "g2600",],aes(x=rt_status,y=CT,color=condition))+geom_jitter(position = position_jitterdodge(dodge.width=1,jitter.height=0, jitter.width = 0.2))

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
rt.data_norm.f<-melt.data.table(as.data.table(rt.data_norm.cast),id.vars = c("gene_prim","replicate","plateID","pH7"),measure.vars = c("pH5.5","pH9"))
rt.data_norm.f<-rt.data_norm.f[!is.na(rt.data_norm.f$value),]

rt.data_norm.f.pH5comp<-melt.data.table(as.data.table(rt.data_norm.cast),id.vars = c("gene_prim","replicate","plateID","pH9"),measure.vars = c("pH7","pH5.5"))
rt.data_norm.f.pH5comp<-rt.data_norm.f.pH5comp[!is.na(rt.data_norm.f.pH5comp$value),]

rt.data_norm.f$delta_delta_Ct<-rt.data_norm.f$value - rt.data_norm.f$pH7 # calculating delta delta Ct for all conditions
rt.data_norm.f.pH5comp$delta_delta_Ct<-rt.data_norm.f.pH5comp$value - rt.data_norm.f.pH5comp$pH9

# Calculate fold gene expression 2^-(ddCt)
rt.data_norm.f$expression_fold_change<-2^-(rt.data_norm.f$delta_delta_Ct)
rt.data_norm.f.pH5comp$expression_fold_change<-2^-(rt.data_norm.f.pH5comp$delta_delta_Ct)

rt.data_norm.f$gene_prim<-factor(rt.data_norm.f$gene_prim, levels=c("clsA","clsB","g1852"))
rt.data_norm.f.pH5comp$gene_prim<-factor(rt.data_norm.f.pH5comp$gene_prim, levels=c("clsA","clsB","g1852"))

rt.data_norm.f<-rt.data_norm.f[rt.data_norm.f$plateID %in% "p2",]

# Plotting
ddCt.allstress.p<-ggplot(rt.data_norm.f,aes(x=gene_prim,y=delta_delta_Ct,color=replicate))+geom_point()+
  as_md_theme(my_theme)+xlab("Gene Amplified")+ylab("\u0394\u0394 Ct (pH 5.5 Ct - pH 7 Ct)")+
  ylim(-4,4)+facet_wrap(~variable,nrow=1)
plot(ddCt.allstress.p)

fc_gene_stat_bar.p<-ggplot(rt.data_norm.f[rt.data_norm.f$gene_prim %in% c("clsA","clsB") & rt.data_norm.f$plateID %in% "p2",], aes(x=gene_prim,y=log2(expression_fold_change),fill=gene_prim))+
  geom_bar(stat="summary",fun="mean")+stat_summary(geom="errorbar",fun.data=mean_sd, width=0.5, linewidth=0.5,position=position_dodge(width=0.8))+
  scale_y_continuous(name="log2(treated/untreated)",limits=c(-2,4),breaks=c(-2,-1,0,1,2,3,4))+
  scale_x_discrete(name="Gene Amplified",limits=c("clsA","clsB"),labels=c("*clsA*","*clsB*"))+
  as_md_theme(my_theme)+theme(legend.position = "none")+
  ggtitle("")+
  scale_fill_manual(name="Gene Amplified",values=c("#299764","#006198"),limits=c("clsB","clsA"),labels=c("*clsB*","*clsA*"))+facet_wrap(~variable,nrow=1)
plot(fc_gene_stat_bar.p)

fc_gene_stat_bar_pH5.5comp.p<-ggplot(rt.data_norm.f.pH5comp, aes(x=gene_prim,y=log2(expression_fold_change),fill=gene_prim))+
  geom_bar(stat="summary",fun="mean")+stat_summary(geom="errorbar",fun.data=mean_sd, width=0.5, linewidth=0.5,position=position_dodge(width=0.8))+
  scale_y_continuous(name="log2(treated/untreated)",limits=c(-2,4),breaks=c(-2,-1,0,1,2,3,4))+
  scale_x_discrete(name="Gene Amplified",limits=c("clsA","clsB"),labels=c("*clsA*","*clsB*"))+
  as_md_theme(my_theme)+theme(legend.position = "none")+
  ggtitle("")+
  scale_fill_manual(name="Gene Amplified",values=c("#299764","#006198"),limits=c("clsB","clsA"),labels=c("*clsB*","*clsA*"))+facet_wrap(~variable,nrow=1)
plot(fc_gene_stat_bar_pH5.5comp.p)

ggsave(plot=ddCt.allstress.p,filename="../figures/rtqpcr_cls_pH_20241106_ddct.pdf",width=4,height=2,units="in",useDingbats=FALSE)
ggsave(plot=fc_gene_stat_bar.p,filename="../figures/rtqpcr_cls_pH_20241106_log2fc.pdf",width=6,height=2,units="in",device=cairo_pdf)
ggsave(plot=fc_gene_stat_bar_pH5.5comp.p,filename="../figures/rtqpcr_cls_pH_20241106_log2fc_pH5comp.pdf",width=6,height=2,units="in",device=cairo_pdf)

# Between log2(fold-changes) of each gene in a condition
t.test(log2(rt.data_norm.f[rt.data_norm.f$variable %in% "pH5.5" & rt.data_norm.f$gene_prim %in% "clsB",]$expression_fold_change),log2(rt.data_norm.f[rt.data_norm.f$variable %in% "pH5.5" & rt.data_norm.f$gene_prim %in% "clsA",]$expression_fold_change))
# p = ns
t.test(log2(rt.data_norm.f[rt.data_norm.f$variable %in% "pH9" & rt.data_norm.f$gene_prim %in% "clsB",]$expression_fold_change),log2(rt.data_norm.f[rt.data_norm.f$variable %in% "pH9" & rt.data_norm.f$gene_prim %in% "clsA",]$expression_fold_change))
# p = 0.099

## Between delta Cts of each gene compared to no treatment
rt.data_norm.f_stats<-melt.data.table(as.data.table(rt.data_norm.cast),id.vars = c("gene_prim","replicate","plateID"),measure.vars = c("pH5.5","pH7","pH9"))
rt.data_norm.f_stats<-rt.data_norm.f_stats[!is.na(rt.data_norm.f_stats$value),]
rt.data_norm.f_stats$variable<-factor(rt.data_norm.f_stats$variable,levels=c("pH7","pH5.5","pH9"))
rt.data_norm.f_stats$gene_prim<-factor(rt.data_norm.f_stats$gene_prim, levels=c("clsA","clsB"))
rt.data_norm.f_stats<-rt.data_norm.f_stats[rt.data_norm.f_stats$plateID %in% "p2",]

clsA.lm<-glm(formula = value ~ variable,data = rt.data_norm.f_stats[rt.data_norm.f_stats$gene_prim %in% "clsA",])
summary(clsA.lm)
# pH 9 0.0173
clsB.lm<-glm(formula = value ~ variable,data = rt.data_norm.f_stats[rt.data_norm.f_stats$gene_prim %in% "clsB",])
summary(clsB.lm)
# pH9 6.85E-5



###########
rt.data.f.subset<-dcast(as.data.table(rt.data.f[rt.data.f$gene_prim %in% c("g2600","g1852") & rt.data.f$condition %!in% c("pH5.5_nort","pH7_nort","pH9_nort"),c("replicate","condition","gene_prim","CT")]),formula = replicate+condition~gene_prim,value.var = "CT",fun.aggregate = mean)

ggplot(rt.data.f.subset,aes(x=g2600,y=g1852,color=replicate))+geom_point()+
  ylim(6,14)+
  xlim(6,14)+
  geom_smooth(method="lm",color="black")

g2600.1852.lm<-glm(data=rt.data.f.subset,formula=g2600~g1852, family="gaussian")
summary(g2600.1852.lm)

# delta delta Ct = delta Ct (treated sample) - delta Ct (untreated sample)
# delta Ct = Ct (gene of interest) - Ct (housekeeping gene)

rt.data.f.agg<-aggregate(CT~gene_prim+condition+replicate,data=rt.data.f,FUN=mean)
rt.data.f.agg.nort<-rt.data.f.agg[rt.data.f.agg$condition %in% c("pH5.5_nort","pH7_nort","pH9_nort"),]
rt.data.f.agg<-rt.data.f.agg[rt.data.f.agg$condition %!in% c("pH5.5_nort","pH7_nort","pH9_nort"),]
rt.data.f.agg$cond_rep<-paste(rt.data.f.agg$condition,rt.data.f.agg$replicate,sep="_")


norm_data<-rt.data.f.agg[rt.data.f.agg$gene_prim %in% "g2600",]

colnames(norm_data)<-paste0(colnames(norm_data),"_norm")


rt.data_norm<-merge(rt.data.f.agg[rt.data.f.agg$gene_prim %!in% "g2600",],norm_data,by.x="cond_rep",by.y="cond_rep_norm",all.x=TRUE)

# data_f_norm_cond<-dcast(data=as.data.table(rt.data_norm),value.var =  c("CT","CT_norm"),formula= gene+condition~gene_normcondition)

rt.data_norm$delta_Ct<-rt.data_norm$CT - rt.data_norm$CT_norm # CT - control gene

# rt.data_norm$dc_logical<-substrRight(rt.data_norm$condition,2)
#rt.data_norm$sample<-substrLeft(rt.data_norm$condition,4)

rt.data_norm$gene_rep<-paste(rt.data_norm$gene_prim,rt.data_norm$replicate,sep="_")
rt.data_norm.ph5<-rt.data_norm[rt.data_norm$condition %in% "pH5.5",c("gene_rep","delta_Ct")]
colnames(rt.data_norm.ph5)<-c("gene_rep","delta_Ct_ph5")
rt.data_norm.ph9<-rt.data_norm[rt.data_norm$condition %in% "pH9",c("gene_rep","delta_Ct")]
colnames(rt.data_norm.ph9)<-c("gene_rep","delta_Ct_ph9")
rt.data_norm.ph7<-rt.data_norm[rt.data_norm$condition %in% "pH7",]


rt.data_norm.f<-merge(rt.data_norm.ph7,rt.data_norm.ph5,by.all="gene_rep",all.x=TRUE,all.y=TRUE)
rt.data_norm.f<-merge(rt.data_norm.f,rt.data_norm.ph9,by.all="gene_rep",all.x=TRUE,all.y=TRUE)

rt.data_norm.f$delta_delta_Ct_ph5.7<-rt.data_norm.f$delta_Ct_ph5 - rt.data_norm.f$delta_Ct # calculating delta delta Ct for clsa strain
rt.data_norm.f$delta_delta_Ct_ph9.7<-rt.data_norm.f$delta_Ct_ph9 - rt.data_norm.f$delta_Ct # calculating delta delta Ct for clsb strain
rt.data_norm.f$delta_delta_Ct_ph5.9<-rt.data_norm.f$delta_Ct_ph5 - rt.data_norm.f$delta_Ct_ph9 # calculating delta delta Ct for between clsa and clsb strains

# Calculate fold gene expression 2^-(ddCt)
rt.data_norm.f$expression_fold_change_ph5.7<-2^-(rt.data_norm.f$delta_delta_Ct_ph5.7)
rt.data_norm.f$expression_fold_change_ph9.7<-2^-(rt.data_norm.f$delta_delta_Ct_ph9.7)
rt.data_norm.f$expression_fold_change_ph5.9<-2^-(rt.data_norm.f$delta_delta_Ct_ph5.9)
# rt.data_norm.f$expression_fold_change_norm<-2^-(rt.data_norm.f$delta_Ct_norm)


rt.data_norm.f$gene_prim<-factor(rt.data_norm.f$gene_prim, levels=c("clsA","clsB","g1852"))

rt.data_norm.f.m<-melt.data.table(as.data.table(rt.data_norm.f),id.vars = c("gene_rep","cond_rep","gene_prim","condition","replicate","CT","gene_prim_norm","condition_norm","replicate_norm","CT_norm","delta_Ct","delta_Ct_ph5","delta_Ct_ph9"))

rt.data_norm.f.m$cond_comp<-"none"
rt.data_norm.f.m[rt.data_norm.f.m$variable %in% c("delta_delta_Ct_ph5.7","expression_fold_change_ph5.7"),]$cond_comp<-"pH7"
rt.data_norm.f.m[rt.data_norm.f.m$variable %in% c("delta_delta_Ct_ph9.7","expression_fold_change_ph9.7"),]$cond_comp<-"pH7"
rt.data_norm.f.m[rt.data_norm.f.m$variable %in% c("delta_delta_Ct_ph5.9","expression_fold_change_ph5.9"),]$cond_comp<-"pH9"

ddCt.wt.ph5.7.p<-ggplot(rt.data_norm.f[!is.na(rt.data_norm.f$delta_delta_Ct_ph5.7),],aes(x=gene_prim,y=delta_delta_Ct_ph5.7,color=replicate))+geom_point()+
  as_md_theme(my_theme)+xlab("Gene Amplified")+ylab("\u0394\u0394 Ct (pH 5.5 Ct - pH 7 Ct)")+
  ylim(-4,4)
plot(ddCt.wt.ph5.7.p)

#ddCt.wt.clsa.small.p<-ggplot(rt.data_norm.f[!is.na(rt.data_norm.f$delta_delta_Ct_wt.clsa),],aes(x=gene_prim,y=delta_delta_Ct_wt.clsa))+geom_point()+
#  ylim(0,5)+
#  as_md_theme(my_theme)+xlab("Gene ID")+ylab("\u0394\u0394 Ct (\u0394clsA Ct - WT Ct)")
# plot(ddCt.wt.clsa.small.p)

ddCt.wt.ph9.7.p<-ggplot(rt.data_norm.f[!is.na(rt.data_norm.f$delta_delta_Ct_ph9.7),],aes(x=gene_prim,y=delta_delta_Ct_ph9.7,color=replicate))+geom_point()+
  as_md_theme(my_theme)+xlab("Gene Amplified")+ylab("\u0394\u0394 Ct (pH 9 Ct - pH 7 Ct)")+
  ylim(-4,4)
plot(ddCt.wt.ph9.7.p)

ddCt.wt.ph5.9.p<-ggplot(rt.data_norm.f[!is.na(rt.data_norm.f$delta_delta_Ct_ph5.9),],aes(x=gene_prim,y=delta_delta_Ct_ph5.9,color=replicate))+geom_point()+
  as_md_theme(my_theme)+xlab("Gene Amplified")+ylab("\u0394\u0394 Ct (pH 9 Ct - pH 7 Ct)")+
  ylim(-4,4)
plot(ddCt.wt.ph5.9.p)

# ddCt.wt.clsb.small.p<-ggplot(rt.data_norm.f[!is.na(rt.data_norm.f$delta_delta_Ct_wt.clsb),],aes(x=gene_prim,y=delta_delta_Ct_wt.clsb))+geom_point()+
#  ylim(0,5)+
#  as_md_theme(my_theme)+xlab("Gene ID")+ylab("\u0394\u0394 Ct (\u0394clsA Ct - WT Ct)")
# plot(ddCt.wt.clsb.small.p)

fc_gene_stat.p<-ggplot(rt.data_norm.f.m[rt.data_norm.f.m$variable %in% c("expression_fold_change_ph5.7","expression_fold_change_ph9.7","expression_fold_change_ph5.9") & rt.data_norm.f.m$gene_prim %in% c("clsA","clsB"),],aes(x=variable,y=log2(value),color=gene_prim))+
  geom_point(position=position_dodge(width=0.2))+
  scale_y_continuous(limits=c(-4,4),breaks=c(-4,-3,-2,-1,0,1,2,3,4))+
  scale_x_discrete(limits=c("expression_fold_change_ph5.7","expression_fold_change_ph9.7","expression_fold_change_ph5.9"),labels=c("pH 5.5/pH 7","pH 9/pH 7","pH 5.5/pH 9"))+
  as_md_theme(my_theme)+xlab("Comparison")+ylab("log2(fold-change)")+ggtitle("")+
  scale_color_manual(name="Gene Amplified",values=c("#006198","#299764"),limits=c("clsA","clsB"),labels=c("*clsA*","*clsB*"))
plot(fc_gene_stat.p)


fc_gene_stat_ph5.9.p<-ggplot(rt.data_norm.f.m[rt.data_norm.f.m$variable %in% "expression_fold_change_ph5.9" & rt.data_norm.f.m$gene_prim %in% c("clsA","clsB"),],aes(x=gene_prim,y=log2(value),color=gene_prim))+
  geom_point(position=position_dodge(width=0.2))+
  scale_y_continuous(limits=c(-4,4),breaks=c(-4,-3,-2,-1,0,1,2,3,4))+
  scale_x_discrete()+
  as_md_theme(my_theme)+theme(legend.position = "none")+
  xlab("Gene Amplified")+ylab("log2(FC) pH 5.5/pH 9")+ggtitle("")+
  scale_color_manual(name="Gene Amplified",values=c("#006198","#299764"),limits=c("clsA","clsB"),labels=c("*clsA*","*clsB*"))
plot(fc_gene_stat_ph5.9.p)

fc_gene_stat_ph5.9_bar.p<-ggplot(rt.data_norm.f.m[rt.data_norm.f.m$variable %in% "expression_fold_change_ph5.9" & rt.data_norm.f.m$gene_prim %in% c("clsA","clsB"),], aes(x=gene_prim,y=log2(value),fill=gene_prim))+
  geom_bar(stat="summary")+stat_summary(geom="errorbar",fun.data=mean_sd, width=0.5, linewidth=0.5,position=position_dodge(width=0.8))+
  scale_y_continuous(limits=c(-2,4),breaks=c(-2,-1,0,1,2,3,4))+
  scale_x_discrete(limits=c("clsA","clsB"),labels=c("\u0394*clsA*","\u0394*clsB*"))+
  as_md_theme(my_theme)+theme(legend.position = "none")+
  xlab("Strain")+ylab("log2(FC) (pH 5.5/pH 9)")+ggtitle("")+
  scale_fill_manual(name="Gene Amplified",values=c("#299764","#006198"),limits=c("clsB","clsA"),labels=c("*clsB*","*clsA*"))
plot(fc_gene_stat_ph5.9_bar.p)


ggsave(plot=fc_gene_stat_ph5.9.p,filename="../figures/rtqpcr_cls_pH_20240319.pdf",width=1.5,height=2,units="in",useDingbats=FALSE)
ggsave(plot=fc_gene_stat_ph5.9_bar.p,filename="../figures/rtqpcr_cls_pH_20240319_bar.pdf",width=1.5,height=2,units="in",device=cairo_pdf)

t.test(log2(rt.data_norm.f.m[rt.data_norm.f.m$variable %in% "expression_fold_change_ph5.9" & rt.data_norm.f.m$gene_prim %in% "clsB",]$value),log2(rt.data_norm.f.m[rt.data_norm.f.m$variable %in% "expression_fold_change_ph5.9" & rt.data_norm.f.m$gene_prim %in% "clsA",]$value))
# p = 0.0182


t.test(rt.data_norm.f[rt.data_norm.f$gene_prim %in% "clsA",]$delta_Ct_ph5,rt.data_norm.f[rt.data_norm.f$gene_prim %in% "clsA",]$delta_Ct)
# p = 0.4732
t.test(rt.data_norm.f[rt.data_norm.f$gene_prim %in% "clsB",]$delta_Ct_ph9,rt.data_norm.f[rt.data_norm.f$gene_prim %in% "clsB",]$delta_Ct)
# p = 0.1069



clsA.lm<-glm(CT~log2(dilution),rt.data.f[rt.data.f$gene %in% "clsA",],family="gaussian")
summary(clsA.lm)

clsB.lm<-glm(CT~log2(dilution),rt.data.f[rt.data.f$gene %in% "clsB",],family="gaussian")
summary(clsB.lm)

g2600.lm<-glm(CT~log2(dilution),rt.data.f[rt.data.f$gene %in% "g2600",],family="gaussian")
summary(clsA.lm)
