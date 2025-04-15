library(ggplot2)
# library(reshape2)
# library(growthcurver)
library(ipolygrowth)
library(ggtext)
library(gridtext)
# library(glue)
library(mdthemes)
library(matrixStats)
library(data.table)
library(magrittr)
library(growthcurver)
library(ggpubr)
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


#####

# In Magellan, format the Tecan output file (ASC-II) with the column/well at the top of the file. Include temperature and time data
# Previously, I manually edited the .asc file to fit formatting requirements; currently modified code to fit the input data

setwd("C:/Users/mksch/Dropbox/matt_data/crosson_lab_msu_schnizlein/bfrag_bile_acid_stress/genetics_bfrag/growth_curves_general/growth_curve_analysis_pooled/")

# Prior to import, needed to remove heading and specify names for first two columns ("Time" and "Temp") and also replace all degree symbols with nothing
data_sds_0329<-read.csv(file="../sds/p207_cls_sds_20240329_updatedconc/p207_sds_20240329_formatted.csv")
data_sds_0523<-read.csv(file="../sds/p207_cls_sds_20240523/p207_sds_20240523_formatted.csv")
# data_dc_0214<-read.csv(file="../dc_stress/p207_clscomp_dcstress_20240214_atccomp_THIS/p207_cls_dc_20240214_comp_formatted.csv")
data_dc_1023<-read.csv(file="../dc_stress/p207_clscomp_dc_20241023/p207_dc_20241023_formatted.csv",header=TRUE)
data_dc_1030<-read.csv(file="../dc_stress/p207_clscomp_dc_20241030/p207_dc_20241030_formatted.csv",header=TRUE)
data_dc_1031<-read.csv(file="../dc_stress/p207_clscomp_dc_20241031/p207_dc_20241031_formatted.csv",header=TRUE)
data_dc_1107<-read.csv(file="../k/p207_cls_k_dc_20241107/p207_k_dc_20241107_formatted.csv",header=TRUE)
# data_dc_0822<-data_dc_0822[data_dc_0822$media_type %in% c("bhis","bhis_dc"),]

data_sds_0329$experiment<-"exp0329_sds"
data_sds_0329$replicate<-1
data_sds_0523$experiment<-"exp0523_sds"
data_dc_1023$experiment<-"exp1023_dc"
data_dc_1023$replicate<-1
data_dc_1030$experiment<-"exp1030_dc"
data_dc_1030$replicate<-2
data_dc_1031$experiment<-"exp1031_dc"
data_dc_1031$replicate<-3
data_dc_1107$experiment<-"exp1107_dc"
data_dc_1107$replicate<-4

data_comb<-rbind(data_sds_0329,data_sds_0523,data_dc_1023,data_dc_1030,data_dc_1031,data_dc_1107)

data_comb$media_type2<-data_comb$media_type
# data_comb[data_comb$media_type %in% c("k030","bhis","na138","pH7"),]$media_type2<-"bhis"
data_comb.f<-data_comb[data_comb$Time <=20,]
data_comb.f$Time2<-round(data_comb.f$Time, digits=1)
data_comb.f$strain_class<-gsub("_.*","",data_comb.f$strain) %>% gsub("cAcB","cAB",.)
data_comb.f<-data_comb.f[!is.na(data_comb.f$strain),]
data_comb.f[data_comb.f$strain %in% "wt_a",]$strain<-"wt"
data_comb.f[data_comb.f$strain %in% c("cAB_a","cAB_b","cAcB"),]$strain<-"cAB"
data_comb.f[data_comb.f$media_type %in% "dc_01",]$media_type2<-"bhis_dc_01"
data_comb.f$str_med2<-paste0(data_comb.f$strain,"_",data_comb.f$media_type2)

data_comb.f<-data_comb.f[data_comb.f$media_type %!in% c("k300","bhis_dc_075","bhis_dc_01_noatc","sds_3","sds_1","dc_075"),]

str_vector<-unique(data_comb.f$strain)

data_comb.f.avg<-aggregate(od~strain+experiment+media_type2+Time2+strain_class+str_med2+replicate,data=data_comb.f,FUN=mean)

### Growth groups
## P207: B fragilis WT P207 from glycerol into BHIS broth grown overnight, backdiluted culture 20x and then diluted to OD 0.025 in the plate
## P207_wt_a: wt w/EVa
## P207_clsA: deletion of clsA w/EV
## P207_clsB: deletion of clsB w/EV
## P207_clsAB_a: deletion of clsA and B w/EV
## P207_clsAcomp: complementation
## P207_clsBcomp: complementation
## P207_clsAB_Acomp: complementation
## P207_clsAB_Bcomp: complementation

## negcont: uninoculated media
### Media groups
## "bhis" =  BHIS
## "bhis_sds_1" = BHIS + 0.001% SDS
## "bhis_sds_2" = BHIS + 0.002% SDS
## "bhis_sds_3" = BHIS + 0.003% SDS
## "bhis_sds_4" = BHIS + 0.004% SDS
## "bhis_dc_01" = BHIS + 0.01% Deoxycholate

#####

breaks_plots<-c("wt_bhis","cA_bhis","cA_comp_bhis","cB_bhis","cB_comp_bhis","cAB_bhis","cAB_acomp_bhis","cAB_bcomp_bhis",
                "wt_sds_2","cA_sds_2","cA_comp_sds_2","cB_sds_2","cB_comp_sds_2","cAB_sds_2","cAB_acomp_sds_2","cAB_bcomp_sds_2",
                "wt_bhis_dc_01","cA_bhis_dc_01","cA_comp_bhis_dc_01","cB_bhis_dc_01","cB_comp_bhis_dc_01","cAB_bhis_dc_01","cAB_acomp_bhis_dc_01","cAB_bcomp_bhis_dc_01")
labels_plots<-c("WT","\u0394 *clsA*","\u0394 *clsA*  att::*clsA*","\u0394 *clsB*","\u0394 *clsB*  att::*clsB*","\u0394 *clsAclsB*","\u0394 *clsAclsB* att::*clsA*","\u0394 *clsAclsB* att::*clsB*",
                "WT (0.002% SDS)","\u0394 *clsA* (0.002% SDS)","\u0394 *clsA*  att::*clsA* (0.002% SDS)","\u0394 *clsB* (0.002% SDS)","\u0394 *clsB*  att::*clsB* (0.002% SDS)","\u0394 *clsAclsB* (0.002% SDS)","\u0394 *clsAclsB* att::*clsA* (0.002% SDS)","\u0394 *clsAclsB* att::*clsB* (0.002% SDS)",
                "WT (0.005% Dc)","\u0394 *clsA* (0.005% Dc)","\u0394 *clsA*  att::*clsA* (0.005% Dc)","\u0394 *clsB* (0.005% Dc)","\u0394 *clsB*  att::*clsB* (0.005% Dc)","\u0394 *clsAclsB* (0.005% Dc)","\u0394 *clsAclsB* att::*clsB* (0.005% Dc)","\u0394 *clsAclsB* att::*clsB* (0.005% Dc)")

# colors_plots<-c("#000000","#424242","#999999","#424242","#999999","#333333","#636363","#999999",
#                "#c35600","#006198","#56B4E9","#0f682c","#299764","#e38200","#46c8a1","#69fae2",
#                "#c35600","#006198","#56B4E9","#0f682c","#299764","#e38200","#46c8a1","#69fae2")

colors_plots<-c("#000000","#006198","#56B4E9","#0f682c","#299764","#e38200","#7D4564","#CC79A7",
                "#000000","#006198","#56B4E9","#0f682c","#299764","#e38200","#7D4564","#CC79A7",
                "#000000","#006198","#56B4E9","#0f682c","#299764","#e38200","#7D4564","#CC79A7")



breaks_1<-c("wt","cA","cA_comp","cB","cB_comp","cAB","cAB_acomp","cAB_bcomp")
# colors_1<-c("#000000","#303030","#464646","#5e5e5e","#767676","#909090","#c35600","#006198","#56B4E9","#0f682c","#299764","#46c8a1")
labels_1<-c("WT",
            "\u0394 *clsA* att:EV","\u0394 *clsA*  att::*clsA*",
            "\u0394 *clsB* att:EV","\u0394 *clsB* att::*clsB*",
            "\u0394 *clsAclsB* att:EV A","\u0394 *clsAclsB* att::*clsA*"
            ,"\u0394 *clsAclsB* att::*clsB*")

colors_1<-c("#c35600","#006198","#56B4E9","#0f682c","#299764","#e38200","#7D4564","#CC79A7")


breaks_2<-c("sds_4","sds_3","sds_2","sds_1")
labels_2<-c("0.004% SDS","0.003% SDS","0.002% SDS","0.001% SDS")
colors_2<-c("#c35600","#56B4E9","#0f682c","#e38200")

# Dc (excluding 10/22 because of increased overall toxicity in that sample set)

mediaall_wt_dc.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type %in% c("bhis","bhis_dc_01") & data_comb.f.avg$strain %in% "wt",],aes(x=Time2,y=od,color=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_sd,geom="ribbon",color="transparent",alpha=0.3)+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_y_continuous(limits=c(-0.1,1))+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_wt_dc.p)

mediaall_ca_dc.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis","bhis_dc_01") & data_comb.f.avg$strain_class %in% c("wt","cA") & data_comb.f.avg$replicate %in% 1:2,],aes(x=Time2,y=od,color=str_med2,group=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_sd,geom="ribbon",color="transparent",alpha=0.3)+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_y_continuous(limits=c(-0.1,1))+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_ca_dc.p)

mediaall_cb_dc.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis","bhis_dc_01") & data_comb.f.avg$strain_class %in% c("wt","cB") & data_comb.f.avg$replicate %in% 1:2,],aes(x=Time2,y=od,color=str_med2,group=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_sd,geom="ribbon",color="transparent",alpha=0.3)+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_y_continuous(limits=c(-0.1,1))+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_cb_dc.p)

mediaall_cab_dc.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis","bhis_dc_01") & data_comb.f.avg$strain_class %in% c("wt","cAB") & data_comb.f.avg$replicate %in% 1:2,],aes(x=Time2,y=od,color=str_med2,group=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_sd,geom="ribbon",color="transparent",alpha=0.3)+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_y_continuous(limits=c(-0.1,1))+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_cab_dc.p)

ggsave(filename = "figures/growthcurve_dc_combined_wt.pdf",plot=mediaall_wt_dc.p,width=2,height = 2,units="in",device=cairo_pdf)
ggsave(filename = "figures/growthcurve_dc_combined_clsa.pdf",plot=mediaall_ca_dc.p,width=2,height = 2,units="in",device=cairo_pdf)
ggsave(filename = "figures/growthcurve_dc_combined_clsb.pdf",plot=mediaall_cb_dc.p,width=2,height = 2,units="in",device=cairo_pdf)
ggsave(filename = "figures/growthcurve_dc_combined_clsaclsb.pdf",plot=mediaall_cab_dc.p,width=2,height = 2,units="in",device=cairo_pdf)

# SDS

mediaall_wt_sds.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type %in% c("bhis","sds_2") & data_comb.f.avg$strain %in% "wt",],aes(x=Time2,y=od,color=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_y_continuous(limits=c(0,1))+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_wt_sds.p)

mediaall_ca_sds.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis","sds_2") & data_comb.f.avg$strain_class %in% c("wt","cA"),],aes(x=Time2,y=od,color=str_med2,group=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_y_continuous(limits=c(0,1))+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_ca_sds.p)

mediaall_cb_sds.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis","sds_2") & data_comb.f.avg$strain_class %in% c("wt","cB"),],aes(x=Time2,y=od,color=str_med2,group=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_y_continuous(limits=c(0,1))+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_cb_sds.p)

mediaall_cab_sds.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis","sds_2") & data_comb.f.avg$strain_class %in% c("wt","cAB"),],aes(x=Time2,y=od,color=str_med2,group=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  as_md_theme(my_theme)+theme(legend.position="none")+
  scale_y_continuous(limits=c(0,1))+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_cab_sds.p)

ggsave(filename = "figures/growthcurve_sds_combined_wt.pdf",plot=mediaall_wt_sds.p,width=2,height = 2,units="in",device=cairo_pdf)
ggsave(filename = "figures/growthcurve_sds_combined_clsa.pdf",plot=mediaall_ca_sds.p,width=2,height = 2,units="in",device=cairo_pdf)
ggsave(filename = "figures/growthcurve_sds_combined_clsb.pdf",plot=mediaall_cb_sds.p,width=2,height = 2,units="in",device=cairo_pdf)
ggsave(filename = "figures/growthcurve_sds_combined_clsaclsb.pdf",plot=mediaall_cab_sds.p,width=2,height = 2,units="in",device=cairo_pdf)

# Summary statistics (ipolygrowth)
# https://github.com/cran/ipolygrowth/tree/master
# data_comb.singleavg<-aggregate(od~strain+media_type2+Time2+strain_class+str_med2,data=data_comb.f.avg,FUN=sum)

# max(data_comb.singleavg[data_comb.singleavg$strain %in% "wt" & data_comb.singleavg$media_type2 %in% "sds_1",]$od)

# library(growthrates)
# testtest<-all_splines(od ~ time | str_med2_rep+strain,data = data_comb.f.avg, spar= 0.5)
# summary(testtest)
# coef(testtest)
# plot(testtest)

data_comb.f.avg$str_med2_rep<-paste0(data_comb.f.avg$str_med2,"_",data_comb.f.avg$replicate)
data_growthcurve<-data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis","sds_2","bhis_dc_01") & data_comb.f.avg$strain_class %in% c("wt","cA","cB","cAB") & data_comb.f.avg$str_med2_rep %!in% c("cAB_bhis_dc_01_1","cAB_bhis_dc_01_3"),]

data_growthcurve.bhis<-data_growthcurve[data_growthcurve$media_type2 %in% "bhis" & data_growthcurve$Time2 >=0.5 & data_growthcurve$Time2 <= 14,]
data_growthcurve.lagged<-data_growthcurve[data_growthcurve$media_type2 %!in% "bhis" & data_growthcurve$Time2 >=0.5,]
data_growthcurve.wtsds<-data_growthcurve[data_growthcurve$media_type2 %in% "sds_2" & data_growthcurve$strain_class %in% "wt" & data_growthcurve$Time2 >=5 & data_growthcurve$Time2 <=15,]

data_growthcurve[data_growthcurve$od <0,]$od<-0

out.multi.f.bhis<-ipg_multisample(data=data_growthcurve.bhis,id="str_med2_rep",time.name="Time2",y.name="od",epsilon = 0.2/100)
out.multi.f.bhis$estimates

out.multi.f.lagged<-ipg_multisample(data=data_growthcurve.lagged,id="str_med2_rep",time.name="Time2",y.name="od",epsilon = 0.2/100)
out.multi.f.lagged$estimates

out.multi.f.wtsds<-ipg_multisample(data=data_growthcurve.wtsds,id="str_med2_rep",time.name="Time2",y.name="od",epsilon = 0.2/100)
out.multi.f.wtsds$estimates

out.multi.f.test<-ipg_multisample(data=data_growthcurve_test,id="str_med2_rep",time.name="Time3",y.name="od",epsilon = 0.2/100)
out.multi.f.test$estimates


ipg.plot_sds.all<-ggplot()+
  geom_point(data =data_growthcurve.wtsds, aes(x = Time2, y = od))+
  geom_line(data = out.multi.f.wtsds$fitted, aes(x = time, y = fit))+
  facet_wrap(~ str_med2_rep)+
  labs(color = "replicate")+
  theme_bw()
plot(ipg.plot_sds.all)

# Trimming growth curves for better fit
growth_curve_metrics<-rbind(out.multi.f.bhis$estimates,out.multi.f.lagged$estimates)
growth_curve_metrics$`peak growth time2`<-round(growth_curve_metrics$`peak growth time`,0)

growthcurves_nogrowth<-unique(growth_curve_metrics[growth_curve_metrics$`peak growth time` <= 2,]$str_med2_rep)

growthcurve_data.f<-as.data.table(matrix(ncol=9,dimnames = list(c(1),c(colnames(data_growthcurve)))))
growthcurve_data.f$mintime<-0
growthcurve_data.f$Time3<-0
strmedrep_breaks<-unique(data_growthcurve$str_med2_rep)[unique(data_growthcurve$str_med2_rep) %!in% growthcurves_nogrowth]

# growthcurvedata_temp<-NULL

for (i in 1:length(strmedrep_breaks)){
  growthcurvedata_temp<-data_growthcurve[data_growthcurve$str_med2_rep %in% strmedrep_breaks[i],]
  growthcurvedata_temp<-growthcurvedata_temp[growthcurvedata_temp$Time2 < (growth_curve_metrics[growth_curve_metrics$str_med2_rep %in% strmedrep_breaks[i],]$`peak growth time` + 6) & growthcurvedata_temp$Time2 > (growth_curve_metrics[growth_curve_metrics$str_med2_rep %in% strmedrep_breaks[i],]$`peak growth time` - 6),]
  growthcurvedata_temp$mintime<-min(growthcurvedata_temp$Time2)
  growthcurvedata_temp$Time3<-growthcurvedata_temp$Time2 - growthcurvedata_temp$mintime
  growthcurve_data.f<-rbind(growthcurve_data.f,growthcurvedata_temp)
}

growthcurve_data.f<-growthcurve_data.f[!is.na(growthcurve_data.f$Time2),]
growthcurve_data.f[growthcurve_data.f$od <= 0,]$od<-0

growthcurve_data.f$exp_rep<-paste0(growthcurve_data.f$experiment,growthcurve_data.f$str_med2_rep)

out.multi.f<-ipg_multisample(data=growthcurve_data.f,id="exp_rep",time.name="Time3",y.name="od",epsilon = 0.2/100)
out.multi.f$estimates

data_growthcurve.unique<-growthcurve_data.f[!duplicated(growthcurve_data.f$exp_rep),c("exp_rep","media_type2","strain","str_med2","str_med2_rep","mintime")]


data_growthcurve.f<-merge(data_growthcurve.unique,out.multi.f$estimates,by="exp_rep",all.x=TRUE)
data_growthcurve.f[,c("peak growth time","lag time","max y time")]<-data_growthcurve.f[,c("peak growth time","lag time","max y time")]+data_growthcurve.f$mintime # adjusting for time taking out earlier to obtain better fit

#  peak growth rate = calculate growth rate (slope) at peak growth times (real number only)
# doubling time = doubling time at peak growth = ln(2)/growth rate
# lag time = the time when the line of intercept from the polynomial model intersects with the linear line of peak growth at peak growth time


# For checking fit of polynomials
ipg.plot.all<-ggplot()+
  geom_point(data = growthcurve_data.f, aes(x = Time3, y = od))+
  geom_line(data = out.multi.f$fitted, aes(x = time, y = fit))+
  facet_wrap(~ exp_rep)+
  labs(color = "replicate")+
  theme_bw()
plot(ipg.plot.all)


ggsave(filename = "figures/ipolygrowth_mapped_dcsds_all.pdf",plot=ipg.plot.all,width=10, height=8,units="in",device=cairo_pdf)

write.csv(x=data_growthcurve.f,file="summarizedgrowth_dcsds.csv",row.names = FALSE)
