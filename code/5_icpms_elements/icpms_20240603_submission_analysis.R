# ICP-MS Submission on 06/03/2024
# Includes experiments performed on 04/24/24, 04/25/24, and 05/31/2024

library(ggplot2)
library(ggtext)
library(gridtext)
# library(glue)
library(mdthemes)
library(matrixStats)
library(data.table)
library(magrittr)
library(readxl)
library(scales)
library(colorBlindness)
library(RColorBrewer)
library(tidyr)


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


breaks_strain<-c("wt","clsA","clsA_comp","clsB","clsB_comp","clsAB")
labels_strain<-c("WT",
                 "\u0394 *clsA*","\u0394 *clsA*  att::*clsA*",
                 "\u0394 *clsB*","\u0394 *clsB* att::*clsB*",
                 "\u0394 *clsAclsB*")
colors_strain<-c("#000000","#006198","#56B4E9","#0f682c","#299764","#e38200")

setwd("C:/Users/mksch/Dropbox/matt_data/crosson_lab_msu_schnizlein/bfrag_bile_acid_stress/icp-ms/")

cfu_0424<-read_xlsx(path="icp_20240424_gad_cls/icp_gad_calc_20240424_cfutiter.xlsx")
cfu_0424$cond<-"bhis"
cfu_0424$exp<-"exp0424"
cfu_0425<-read_xlsx(path="icp_20240425_gad_cls/icp_gad_calc_20240425_cfutiter.xlsx")
cfu_0425$cond<-"bhis"
cfu_0425$exp<-"exp0425"
cfu_0531<-read_xlsx(path="icp_20240531_ionophore_dc/icp_20240531_ionophore_cfu.xlsx")
cfu_0531$exp<-"exp0531"
cfu_data<-rbind(cfu_0424,cfu_0425,cfu_0531)

submission_data<-read_xlsx(path="icp_submission_20240603.xlsx")

####
cfu_data.m<-melt.data.table(as.data.table(cfu_data),id.vars = c("tubeID","strain","biological_replicate","sampling_replicate","cond","exp"))
cfu_data.m$colonies<-gsub("_.*","",cfu_data.m$value) %>% as.numeric(.)
cfu_data.m$dilution<-gsub(".*_","",cfu_data.m$value) %>% as.numeric(.)
cfu_data.m[cfu_data.m$exp %in% c("exp0424","exp0425"),]$dilution<- -5
cfu_data.m$dilution<-paste0("1e",cfu_data.m$dilution) %>% as.numeric(.)

cfu_data.m$spot_vol_ul<-5
cfu_data.m[cfu_data.m$exp %in% "exp0531",]$spot_vol_ul<-8

# term 1: colonies/uL in spot; term 2: convert to colonies/mL in spot; term 3: convert to colonies/mL in tube
cfu_data.m$cfu_mL<-(cfu_data.m$colonies/cfu_data.m$spot_vol_ul)*(1000/1)*(1/cfu_data.m$dilution)

cfu_data.m$sampling_replicate<-as.factor(cfu_data.m$sampling_replicate)

ggplot(cfu_data.m,aes(x=tubeID,y=cfu_mL,color=sampling_replicate))+geom_jitter()+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)), limits=c(1e8, 1e10))+
  facet_wrap(~exp)

cfu_data.m.avg<-NULL
cfu_data.m.avg<-aggregate(cfu_mL~tubeID+strain+biological_replicate+sampling_replicate+cond+exp, data=cfu_data.m[,c("tubeID","strain","biological_replicate","sampling_replicate","variable","cond","exp","cfu_mL")],FUN=mean)
cfu_data.m.avg.f<-aggregate(cfu_mL~tubeID+strain+biological_replicate+cond+exp, data=cfu_data.m.avg,FUN=mean)


cfu_data.final<-cfu_data.m.avg.f[,c("tubeID","cfu_mL")]

subdata.f<-merge(submission_data,cfu_data.final,by.all="tubeID",all.x = TRUE)
subdata.f$cfu_total<-subdata.f$cfu_mL*5

ggplot(subdata.f,aes(x=optical_density,y=cfu_total,color=strain))+geom_jitter()+geom_smooth(method="lm",color="black")

hist(subdata.f$cfu_total)
hist(subdata.f$optical_density)

cfu.od.lm<-lm(log10(cfu_total)~optical_density,data=subdata.f)
summary(cfu.od.lm)
# write.csv(subdata.f,file = "icp_submission_20240603_final.csv",row.names = FALSE)
cfu_0918<-read.csv(file="icp_20230918_gad_cls_vtpK/icp_20230918_cfu_averaged.csv")
cfu_0918<-cfu_0918[,c("tubeID","strain","growth_replicate","cfu_ml")]
colnames(cfu_0918)<-c("tubeID","strain","replicate","cfu_mL")
cfu_0918$experiment<-"20230918"
cfu_0918$condition<-"bhis"
cfu_0918$strain<-gsub("cAcB","clsAB",cfu_0918$strain) %>% gsub("cA","clsA",.) %>% gsub("cB","clsB",.)

subdata_comb<-rbind(subdata.f[subdata.f$sample_type %in% "cell_pellet",colnames(cfu_0918)],cfu_0918[cfu_0918$strain %!in% "vtpK",])
subdata_comb$strain<-factor(subdata_comb$strain,levels=breaks_strain)


ggplot(subdata.f[subdata.f$condition %in% "bhis",],aes(x=cfu_mL,y=optical_density,color=strain))+geom_point()+

  scale_color_discrete()+
  as_md_theme(my_theme)

# scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x,n=4), labels = trans_format("log10", label_math(10^.x)),limits=c(1E8,1E9),name="log10(CFU/mL)")+
  annotation_logticks(base=10,sides = "b",scaled=TRUE,short = unit(0.05, "cm"),mid = unit(0.1, "cm"),long = unit(0.2, "cm"))

strain.cfu.p<-ggplot(subdata_comb[subdata_comb$condition %in% "bhis",],aes(x=strain,y=cfu_mL,fill=strain))+geom_boxplot()+
  scale_y_continuous(name="CFU/mL",limits=c(1e8,8e8),breaks=c(2e8,4e8,6e8,8e8))+
  scale_x_discrete(name="",breaks=breaks_strain,labels=labels_strain)+
  scale_fill_manual(breaks=breaks_strain,labels=labels_strain,values=colors_strain)+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(strain.cfu.p)

strain.od.p<-ggplot(subdata.f[subdata.f$condition %in% "bhis",],aes(x=strain,y=optical_density,fill=strain))+geom_boxplot()+
  scale_y_continuous()+
  scale_x_discrete()+

  as_md_theme(my_theme)+theme(legend.position = "none")
plot(strain.od.p)

subdata.f$strain<-factor(subdata.f$strain,levels=c("wt","clsA","clsA_comp","clsB","clsB_comp","clsAB"))
subdata.f$strain2<-factor(subdata.f$strain,levels=c("clsAB","wt","clsA","clsA_comp","clsB","clsB_comp"))

subdata.f$cfu_odnorm<-subdata.f$optical_density/(subdata.f$cfu_mL/1E9)
subdata.f$cfu_odnorm_wtfc<-subdata.f$cfu_odnorm/mean(subdata.f[!is.na(subdata.f$cfu_odnorm) & subdata.f$strain %in% "wt" & subdata.f$condition %in% "bhis",]$cfu_odnorm)

strain.od.p<-ggplot(subdata.f[!is.na(subdata.f$optical_density) & subdata.f$condition %in% "bhis",],aes(x=strain,y=optical_density,fill=strain))+geom_boxplot()+
  scale_y_continuous(name="OD600")+
  scale_x_discrete()+
  scale_fill_manual(breaks=breaks_strain,labels=labels_strain,values=colors_strain)+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(strain.od.p)

strain.cfu.odnorm.p<-ggplot(subdata.f[!is.na(subdata.f$optical_density) & subdata.f$condition %in% "bhis",],aes(x=strain,y=cfu_odnorm,fill=strain))+geom_boxplot()+
  scale_y_continuous(name="OD600/CFU per mL")+
  scale_x_discrete()+
  scale_fill_manual(breaks=breaks_strain,labels=labels_strain,values=colors_strain)+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(strain.cfu.odnorm.p)

od.strain.lm<-lm(optical_density~strain,data=subdata.f[!is.na(subdata.f$cfu_odnorm) & subdata.f$condition %in% "bhis",])
summary(od.strain.lm)
cfu.odnorm.strain.lm<-lm(cfu_odnorm~strain,data=subdata.f[!is.na(subdata.f$cfu_odnorm) & subdata.f$condition %in% "bhis",])
summary(cfu.odnorm.strain.lm)
cfu.odnorm.strain2.lm<-lm(cfu_odnorm~strain2,data=subdata.f[!is.na(subdata.f$cfu_odnorm) & subdata.f$condition %in% "bhis",])
summary(cfu.odnorm.strain2.lm)

cfu.strain.lm<-lm(log10(cfu_mL)~strain,data=subdata_comb[subdata_comb$condition %in% "bhis",])
summary(cfu.strain.lm)
# None significant

treat.cfu.p<-ggplot(subdata.f[subdata.f$experiment %in% "20240531",],aes(x=condition,y=cfu_mL,color=strain))+geom_point(color="black")+
  scale_y_continuous(name="CFU/mL",limits=c(1E8,8E8))+
  scale_x_discrete(name="",breaks=c("bhis","dc","monensin","nigericin"),labels=c("BHIS","0.01% Dc","0.8 ug/mL Monensin A","0.1 ug/mL Nigericin"))+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(treat.cfu.p)

cfu.treat.lm<-lm(log10(cfu_mL)~condition,data=subdata_comb[subdata_comb$strain %in% "wt",])
summary(cfu.treat.lm)

ggsave(filename = "figures/revised_202407/cfu_treat.pdf",width=2,height=2,units="in",device=cairo_pdf)
ggsave(filename = "figures/revised_202407/cfu_strain.pdf",width=2,height=2,units="in",device=cairo_pdf)
ggsave(filename = "figures/revised_202407/OD_cfunorm_strain.pdf",plot = strain.cfu.odnorm.p,width=2,height=2,units="in",device=cairo_pdf)
