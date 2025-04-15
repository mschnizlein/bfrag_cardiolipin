# Summary statistics for all growth curves related to the cls manuscript
# 06/14/24

library(ggplot2)
# library(reshape2)
library(growthcurver)
library(ggtext)
library(gridtext)
# library(glue)
library(mdthemes)
library(matrixStats)
library(data.table)
library(magrittr)
library(ggpubr)
# library(data.table)
# library(genefilter)

setwd("C:/Users/mksch/Dropbox/matt_data/crosson_lab_msu_schnizlein/bfrag_bile_acid_stress/genetics_bfrag/growth_curves_general/growth_curve_analysis_pooled")

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

# These files are output of the ipolygrowth package

gc_nahk<-read.csv(file="summarizedgrowth_NaHK.csv",header = TRUE)
gc_ionoph<-read.csv(file="summarizedgrowth_ionophore.csv",header=TRUE)
gc_dcsds<-read.csv(file="summarizedgrowth_dcsds.csv",header=TRUE)

gc_combined<-rbind(gc_nahk,gc_ionoph,gc_dcsds) %>% as.data.table(.)

strain_factors<-c("wt","cA","cA_comp","cB","cB_comp","cAB","cAB_acomp","cAB_bcomp")
cond_factors<-c("bhis","bhis_dc_01","sds_2","na300","k300","pH5","pH9","nigericin_0.1","monensin_0.8")
cond_labels<-c("BHIS","0.01% Dc","0.002% SDS","300 mM Na","300 mM K","pH 5","pH 9","0.1 ug/mL Nigericin","0.8 ug/mL Monensin")
strain_labels<-c("WT",
            "\u0394 *clsA* att:EV","\u0394 *clsA*  att::*clsA*",
            "\u0394 *clsB* att:EV","\u0394 *clsB* att::*clsB*",
            "\u0394 *clsAclsB* att:EV A","\u0394 *clsAclsB* att::*clsA*"
            ,"\u0394 *clsAclsB* att::*clsB*")
strain_colors<-c("black","#006198","#56B4E9","#0f682c","#299764","#e38200","#7D4564","#CC79A7")

gc_combined.m<-melt.data.table(gc_combined,id.vars=c("strain","media_type2","str_med2","str_med2_rep","mintime","exp_rep"),measure.vars = c("peak.growth.rate","peak.growth.time","doubling.time","lag.time","max.y","max.y.time"))
gc_combined.m[gc_combined.m$strain %in% "WT",]$strain<-"wt"
gc_combined.m[gc_combined.m$media_type2  %in% "bhis_0",]$media_type2<-"bhis"
gc_combined.m$strain<-factor(gc_combined.m$strain, levels=strain_factors)
gc_combined.m$media_type2<-factor(gc_combined.m$media_type2, levels=cond_factors)

# number of replicates
gc_combined[gc_combined$media_type2 %in% c("bhis","bhis_0") & gc_combined$strain %in% c("wt","WT"),] # 9 bhis replicates
gc_combined[gc_combined$media_type2 %in% c("bhis_dc_01") & gc_combined$strain %in% c("cA"),] # 4 dc replicates (only two for WT)
gc_combined[gc_combined$media_type2 %in% c("sds_2") & gc_combined$strain %in% c("wt","WT"),] # 3 sds replicates
gc_combined[gc_combined$media_type2 %in% c("na300") & gc_combined$strain %in% c("cB"),] # 3 Na replicates (only 2 for cA do it failing to fit a polynomial)
gc_combined[gc_combined$media_type2 %in% c("k300") & gc_combined$strain %in% c("wt","WT"),] # 3 K replicates
gc_combined[gc_combined$media_type2 %in% c("pH5") & gc_combined$strain %in% c("wt","WT"),] # 3 pH replicates
gc_combined[gc_combined$media_type2 %in% c("monensin_0.8") & gc_combined$strain %in% c("wt","WT"),] # 3 monensin replicates
gc_combined[gc_combined$media_type2 %in% c("nigericin_0.1") & gc_combined$strain %in% c("wt","WT"),] # 3 nigericin replicates

summ_stats_bhis.p<-ggplot(gc_combined.m[gc_combined.m$media_type2 %in% "bhis",],aes(x=strain,y=value,color=variable))+geom_point()+facet_wrap(~variable)+as_md_theme(my_theme2)
plot(summ_stats_bhis.p)

summ_stats_conditions_maxy.p<-ggplot(gc_combined.m[gc_combined.m$variable %in% c("max.y") & gc_combined.m$media_type2 %in% c("bhis","bhis_dc_01","sds_2","k300","na300","monensin_0.8","nigericin_0.1"),],aes(x=media_type2,y=value,color=strain,group=strain))+
  stat_summary(fun=mean,fun.min=mean,fun.max=mean,geom="crossbar",width=0.5,linewidth=0.5,position=position_dodge(width=0.8)) +stat_summary(geom="errorbar",fun.data=mean_sd, width=0.5, linewidth=0.5,position=position_dodge(width=0.8))+
  scale_y_continuous(limits=c(0,1),name = "Max OD600")+scale_x_discrete(breaks=cond_factors,labels=cond_labels, name= "")+
  scale_fill_manual(breaks=strain_factors,labels=strain_labels,values=strain_colors)+
  scale_color_manual(breaks=strain_factors,labels=strain_labels,values=strain_colors)+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(summ_stats_conditions_maxy.p)

gc_combined.m$value2<-gc_combined.m$value
gc_combined.m[gc_combined.m$variable %in% "doubling.time" & gc_combined.m$value >110,]$value2<-80

summ_stats_conditions_doubletime.p<-ggplot(gc_combined.m[gc_combined.m$variable %in% c("doubling.time") & gc_combined.m$media_type2 %in% c("bhis","bhis_dc_01","sds_2","na300","k300","nigericin_0.1","monensin_0.8"),],aes(x=media_type2,y=value2,color=strain,group=strain))+
  stat_summary(fun=mean,fun.min=mean,fun.max=mean,geom="crossbar",width=0.5,linewidth=0.5,position=position_dodge(width=0.8))+stat_summary(geom="errorbar",fun.data=mean_sd, width=0.5, linewidth=0.5,position=position_dodge(width=0.9))+
  scale_y_continuous(name="Doubling Time (hrs)",breaks=c(0,25,50,75,100))+scale_x_discrete(breaks=cond_factors,labels=cond_labels, name= "")+
  scale_fill_manual(breaks=strain_factors,labels=strain_labels,values=strain_colors)+
  scale_color_manual(breaks=strain_factors,labels=strain_labels,values=strain_colors)+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(summ_stats_conditions_doubletime.p)

ggsave(plot=summ_stats_conditions_maxy.p,filename="figures/gc_stats_maxy.pdf",width=4,height=1.5,units="in",device=cairo_pdf)
ggsave(plot=summ_stats_conditions_doubletime.p,filename="figures/gc_stats_doubletime.pdf",width=4,height=1.5,units="in",device=cairo_pdf)

##
# calculating fold change of doubling time and max OD from WT

gc_combined_wt<-aggregate(value~strain+media_type2+variable, data=gc_combined.m[gc_combined.m$strain %in% "wt",],FUN=mean)

gc_combined.m$value_fc<-0

metric_vector<-unique(gc_combined.m$variable) %>% as.character(.)
media_vector<-unique(gc_combined.m$media_type2) %>% as.character(.)

for (i in 1:6){
  for (j in 1:9){
    gc_combined.m[gc_combined.m$variable %in% metric_vector[i] & gc_combined.m$media_type2 %in% media_vector[j],]$value_fc<-gc_combined.m[gc_combined.m$variable %in% metric_vector[i] & gc_combined.m$media_type2 %in% media_vector[j],]$value/gc_combined_wt[gc_combined_wt$variable %in% metric_vector[i] & gc_combined_wt$media_type2 %in% media_vector[j],]$value
  }
}

summ_stats_conditions_maxy_fc.p<-ggplot(gc_combined.m[gc_combined.m$variable %in% c("max.y") & gc_combined.m$media_type2 %!in% c("sds_2","pH9"),],aes(x=media_type2,y=log2(value_fc),color=strain,group=strain))+
  stat_summary(fun=mean,fun.min=mean,fun.max=mean,geom="crossbar",width=0.5,linewidth=0.5,position=position_dodge(width=0.8)) + stat_summary(geom="errorbar",fun.data=mean_sd, width=0.5, linewidth=0.5,position=position_dodge(width=0.8))+
  scale_y_continuous(name = "log2(FC) Max OD600 (Strain/WT)",limits=c(-3,1))+scale_x_discrete(breaks=cond_factors,labels=cond_labels, name= "")+
  scale_fill_manual(breaks=strain_factors,labels=strain_labels,values=strain_colors)+
  scale_color_manual(breaks=strain_factors,labels=strain_labels,values=strain_colors)+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(summ_stats_conditions_maxy_fc.p)

summ_stats_conditions_maxy_fc_sds.p<-ggplot(gc_combined.m[gc_combined.m$variable %in% c("max.y") & gc_combined.m$media_type2 %in% c("sds_2"),],aes(x=media_type2,y=log2(value_fc),color=strain,group=strain))+
  stat_summary(fun=mean,fun.min=mean,fun.max=mean,geom="crossbar",width=0.5,linewidth=0.5,position=position_dodge(width=0.8)) +stat_summary(geom="errorbar",fun.data=mean_sd, width=0.5, linewidth=0.5,position=position_dodge(width=0.8))+
  scale_y_continuous(name = "log2(FC) Max OD600 (Strain/WT)")+scale_x_discrete(breaks=cond_factors,labels=cond_labels, name= "")+
  scale_fill_manual(breaks=strain_factors,labels=strain_labels,values=strain_colors)+
  scale_color_manual(breaks=strain_factors,labels=strain_labels,values=strain_colors)+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(summ_stats_conditions_maxy_fc_sds.p)


summ_stats_conditions_doubletime_fc.p<-ggplot(gc_combined.m[gc_combined.m$variable %in% c("doubling.time") & gc_combined.m$media_type2 %!in% c("sds_2","pH9"),],aes(x=media_type2,y=log2(value_fc),color=strain,group=strain))+
  stat_summary(fun=mean,fun.min=mean,fun.max=mean,geom="crossbar",width=0.5,linewidth=0.5,position=position_dodge(width=0.8)) + stat_summary(geom="errorbar",fun.data=mean_sd, width=0.5, linewidth=0.5,position=position_dodge(width=0.8))+
  scale_y_continuous(name="log2(FC) Doubling Time (hrs) (Strain/mean(WT))")+scale_x_discrete(breaks=cond_factors,labels=cond_labels, name= "")+
  scale_fill_manual(breaks=strain_factors,labels=strain_labels,values=strain_colors)+
  scale_color_manual(breaks=strain_factors,labels=strain_labels,values=strain_colors)+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(summ_stats_conditions_doubletime_fc.p)

summ_stats_conditions_doubletime_fc_sds.p<-ggplot(gc_combined.m[gc_combined.m$variable %in% c("doubling.time") & gc_combined.m$media_type2 %in% c("sds_2"),],aes(x=media_type2,y=log2(value_fc),color=strain,group=strain))+
  stat_summary(fun=mean,fun.min=mean,fun.max=mean,geom="crossbar",width=0.5,linewidth=0.5,position=position_dodge(width=0.8)) + stat_summary(geom="errorbar",fun.data=mean_sd, width=0.5, linewidth=0.5,position=position_dodge(width=0.8))+
  scale_y_continuous(name="log2(FC) Doubling Time (hrs) (Strain/mean(WT))")+scale_x_discrete(breaks=cond_factors,labels=cond_labels, name= "")+
  scale_fill_manual(breaks=strain_factors,labels=strain_labels,values=strain_colors)+
  scale_color_manual(breaks=strain_factors,labels=strain_labels,values=strain_colors)+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(summ_stats_conditions_doubletime_fc_sds.p)

ggsave(plot=summ_stats_conditions_maxy_fc.p,filename="figures/gc_stats_maxy_fc.pdf",width=4.5,height=3,units="in",device=cairo_pdf)
ggsave(plot=summ_stats_conditions_doubletime_fc.p,filename="figures/gc_stats_doubletime_fc.pdf",width=4.5,height=3,units="in",device=cairo_pdf)
ggsave(plot=summ_stats_conditions_maxy_fc_sds.p,filename="figures/gc_stats_maxy_fc_sds.pdf",width=1.5,height=3,units="in",device=cairo_pdf)
ggsave(plot=summ_stats_conditions_doubletime_fc_sds.p,filename="figures/gc_stats_doubletime_fc_sds.pdf",width=1.5,height=3,units="in",device=cairo_pdf)


maxy_bhis.lm<-glm(log2(value_fc)~strain,data=gc_combined.m[gc_combined.m$media_type2 %in% "bhis" & gc_combined.m$variable %in% "max.y",])
summary(maxy_bhis.lm) # none
maxy_dc.lm<-glm(log2(value_fc)~strain,data=gc_combined.m[gc_combined.m$media_type2 %in% "bhis_dc_01" & gc_combined.m$variable %in% "max.y",])
summary(maxy_dc.lm) # cAB 0.00825, cAB a 0.03009
maxy_sds.lm<-glm(log2(value_fc)~strain,data=gc_combined.m[gc_combined.m$media_type2 %in% "sds_2" & gc_combined.m$variable %in% "max.y",])
summary(maxy_sds.lm) # cB 0.0036, cAB 0.00459
maxy_na.lm<-glm(log2(value_fc)~strain,data=gc_combined.m[gc_combined.m$media_type2 %in% "na300" & gc_combined.m$variable %in% "max.y",])
summary(maxy_na.lm) # ns
maxy_k.lm<-glm(log2(value_fc)~strain,data=gc_combined.m[gc_combined.m$media_type2 %in% "k300" & gc_combined.m$variable %in% "max.y",])
summary(maxy_k.lm) # cA 0.0267, cAB 0.00533
maxy_ph5.lm<-glm(log2(value_fc)~strain,data=gc_combined.m[gc_combined.m$media_type2 %in% "pH5" & gc_combined.m$variable %in% "max.y",])
summary(maxy_ph5.lm) # none
maxy_nig.lm<-glm(log2(value_fc)~strain,data=gc_combined.m[gc_combined.m$media_type2 %in% "nigericin_0.1" & gc_combined.m$variable %in% "max.y",])
summary(maxy_nig.lm) # none (non-transformed cA 0.0326, cAB_bcomp 0.00752)
maxy_mon.lm<-glm(log2(value_fc)~strain,data=gc_combined.m[gc_combined.m$media_type2 %in% "monensin_0.8" & gc_combined.m$variable %in% "max.y",])
summary(maxy_mon.lm) # cA 0.0337, cB 0.004647, cAB 0.0251, cAB_acomp 0.000627, cAB_bcomp 0.043116

dt_bhis.lm<-glm(log2(value_fc)~strain,data=gc_combined.m[gc_combined.m$media_type2 %in% "bhis" & gc_combined.m$variable %in% "doubling.time",])
summary(dt_bhis.lm) # none
dt_dc.lm<-glm(log2(value_fc)~strain,data=gc_combined.m[gc_combined.m$media_type2 %in% "bhis_dc_01" & gc_combined.m$variable %in% "doubling.time",])
summary(dt_dc.lm) # cAB 0.00789
dt_sds.lm<-glm(log2(value_fc)~strain,data=gc_combined.m[gc_combined.m$media_type2 %in% "sds_2" & gc_combined.m$variable %in% "doubling.time",])
summary(dt_sds.lm) # cB 0.0013, cAB 0.000705
dt_na.lm<-glm(log2(value_fc)~strain,data=gc_combined.m[gc_combined.m$media_type2 %in% "na300" & gc_combined.m$variable %in% "doubling.time",])
summary(dt_na.lm) # na
dt_k.lm<-glm(log2(value_fc)~strain,data=gc_combined.m[gc_combined.m$media_type2 %in% "k300" & gc_combined.m$variable %in% "doubling.time",])
summary(dt_k.lm) # cAB 0.0496
dt_ph5.lm<-glm(log2(value_fc)~strain,data=gc_combined.m[gc_combined.m$media_type2 %in% "pH5" & gc_combined.m$variable %in% "doubling.time",])
summary(dt_ph5.lm) # cB 0.0341
dt_nig.lm<-glm(log2(value_fc)~strain,data=gc_combined.m[gc_combined.m$media_type2 %in% "nigericin_0.1" & gc_combined.m$variable %in% "doubling.time",])
summary(dt_nig.lm) # none sig
dt_mon.lm<-glm(log2(value_fc)~strain,data=gc_combined.m[gc_combined.m$media_type2 %in% "monensin_0.8" & gc_combined.m$variable %in% "doubling.time",])
summary(dt_mon.lm) # cA 0.000919, cB 7.09E-6, cAB 0.0177, cAB_acomp 0.000116, cAB_bcomp 0.000815
