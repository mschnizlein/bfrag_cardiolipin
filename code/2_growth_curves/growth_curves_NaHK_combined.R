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
data_na0412<-read.csv(file="../na/p207_cls_na_20240412/p207_cls_na_20240412_formatted.csv",fill=NA,header=TRUE)
data_na0417<-read.csv(file="../na/p207_cls_na_20240417/p207_cls_na_20240417_formatted.csv",fill=NA,header=TRUE)
data_na0927<-read.csv(file="../na/p207_cls_na_dc_20240927_dc5/p207_cls_na_dc_20240927_formatted.csv",fill=NA,header=TRUE)
data_k0418<-read.csv(file="../k/p207_cls_k_20240418/p207_cls_k_20240418_formatted.csv",fill=NA,header=TRUE)
data_k1105<-read.csv(file="../k/p207_cls_k_20241105/p207_k_20241105_formatted.csv",fill=NA,header=TRUE)
data_k1107<-read.csv(file="../k/p207_cls_k_dc_20241107/p207_k_dc_20241107_formatted.csv",fill=NA,header=TRUE)
data_h0229<-read.csv(file="../ph/p207_clscomp_phstress_20240229/p207_clscomp_phstress_20240229_formatted.csv")
data_h0822<-read.csv(file="../dc_stress/p207_clscomp_dcstress_20240822_lowtoxicity/p207_ph_dc_20240822_formatted.csv")
data_h0822<-data_h0822[data_h0822$media_type %in% c("bhis","bhis_ph"),]

# Add Na 9/27 (1 replicate); K 9/26 (2 replicates)

data_na0412$experiment<-"exp0412_na"
data_na0417$experiment<-"exp0417_na"
data_na0927$experiment<-"exp0927_na"
data_k0418$experiment<-"exp0418_k"
data_k1105$experiment<-"exp1105_k"
data_k1107$experiment<-"exp1107_k"
data_h0229$experiment<-"exp0229_h"
data_h0822$experiment<-"exp0822_h"

data_na0412$replicate<-1
data_na0417$replicate<-2
data_na0927$replicate<-3
data_k0418$replicate<-1
data_k1105$replicate<-2
data_k1107$replicate<-3
data_h0229$replicate<-1
data_h0822$replicate<-data_h0822$replicate+1


data_comb<-rbind(data_na0412,data_na0417,data_na0927,data_k0418,data_k1105,data_k1107,data_h0229,data_h0822)

# Plotting multiple experiments next to each other (Na 04/12, Na 04/17, K 04/18)
### Growth groups
## P207: B fragilis WT P207 from glycerol into BHIS broth grown overnight, backdiluted culture 20x and then diluted to OD 0.025 in the plate
## P207_wt_a: wt w/EV
## P207_clsA: deletion of clsA w/EV
## P207_clsB: deletion of clsB w/EV
## P207_clsAB_a: deletion of clsA and B w/EV
## P207_clsAcomp: complementation
## P207_clsBcomp: complementation
## P207_clsAB_Acomp: complementation
## P207_clsAB_Bcomp: complementation

## negcont: uninoculated media
### Media groups
## "bhis" =  BHIS (30 mM K, 138 mM Na, pH 7)
## "k300" = BHIS (300 mM K)
## "k450" = BHIS (450 mM K)
## "k600" = BHIS (600 mM K)
## "na300" = BHIS (300 mM Na)
## "na450" = BHIS (450 mM Na)
## "na600" = BHIS (600 mM Na)
## "pH 5" = BHIS adjusted to pH 5 with HCl
## "pH 9" = BHIS adjusted to pH 9 with NaOH


data_comb$media_type2<-data_comb$media_type
data_comb[data_comb$media_type %in% c("k030","bhis","na138","pH7"),]$media_type2<-"bhis"
data_comb.f<-data_comb[data_comb$Time <=20,]
data_comb.f$Time2<-round(data_comb.f$Time, digits=1)
data_comb.f$strain_class<-substrLeft(data_comb.f$strain,3) %>% gsub("_","",.)
data_comb.f<-data_comb.f[!is.na(data_comb.f$strain),]
data_comb.f[data_comb.f$strain %in% c("WT_a","wt"),]$strain<-"WT"
data_comb.f[data_comb.f$strain_class %in% c("wt"),]$strain_class<-"WT"
data_comb.f[data_comb.f$strain %in% c("cAB_a","cAB_b"),]$strain<-"cAB"
data_comb.f[data_comb.f$media_type %in% "bhis_ph",]$media_type2<-"pH5"
data_comb.f$str_med2<-paste0(data_comb.f$strain,"_",data_comb.f$media_type2)

str_vector<-unique(data_comb.f$strain)

data_comb.f.avg<-aggregate(od~strain+experiment+media_type2+Time2+strain_class+str_med2+replicate,data=data_comb.f,FUN=mean)

#####
breaks_wt_ion<-c("bhis","k300","na300")

breaks_plots<-c("WT_bhis","cA_bhis","cA_comp_bhis","cB_bhis","cB_comp_bhis","cAB_bhis","cAB_acomp_bhis","cAB_bcomp_bhis",
                "WT_na300","cA_na300","cA_comp_na300","cB_na300","cB_comp_na300","cAB_na300","cAB_acomp_na300","cAB_bcomp_na300",
                "WT_k300","cA_k300","cA_comp_k300","cB_k300","cB_comp_k300","cAB_k300","cAB_acomp_k300","cAB_bcomp_k300",
                "WT_pH5","cA_pH5","cA_comp_pH5","cB_pH5","cB_comp_pH5","cAB_pH5","cAB_acomp_pH5","cAB_bcomp_pH5",
                "WT_pH9","cA_pH9","cA_comp_pH9","cB_pH9","cB_comp_pH9","cAB_pH9","cAB_acomp_pH9","cAB_bcomp_pH9")
labels_plots<-c("WT","\u0394 *clsA*","\u0394 *clsA*  att::*clsA*","\u0394 *clsB*","\u0394 *clsB*  att::*clsB*","\u0394 *clsAclsB*","\u0394 *clsAclsB* att::*clsA*","\u0394 *clsAclsB* att::*clsB*",
                "WT (300 mM Na+)","\u0394 *clsA* (300 mM Na+)","\u0394 *clsA*  att::*clsA* (300 mM Na+)","\u0394 *clsB* (300 mM Na+)","\u0394 *clsB*  att::*clsB* (300 mM Na+)","\u0394 *clsAclsB* (300 mM Na+)","\u0394 *clsAclsB* att::*clsA* (300 mM Na+)","\u0394 *clsAclsB* att::*clsB* (300 mM Na+)",
                "WT (300 mM K+)","\u0394 *clsA* (300 mM K+)","\u0394 *clsA*  att::*clsA* (300 mM K+)","\u0394 *clsB* (300 mM K+)","\u0394 *clsB*  att::*clsB* (300 mM K+)","\u0394 *clsAclsB* (300 mM K+)","\u0394 *clsAclsB* att::*clsB* (300 mM K+)","\u0394 *clsAclsB* att::*clsB* (300 mM K+)",
                "WT (pH 5)","\u0394 *clsA* (pH 5)","\u0394 *clsA*  att::*clsA* (pH 5)","\u0394 *clsB* (pH 5)","\u0394 *clsB*  att::*clsB* (pH 5)","\u0394 *clsAclsB* (pH 5)","\u0394 *clsAclsB* att::*clsA* (pH 5)","\u0394 *clsAclsB* att::*clsB* (pH 5)",
                "WT (pH 9)","\u0394 *clsA* (pH 9)","\u0394 *clsA*  att::*clsA* (pH 9)","\u0394 *clsB* (pH 9)","\u0394 *clsB*  att::*clsB* (pH 9)","\u0394 *clsAclsB* (pH 9)","\u0394 *clsAclsB* att::*clsA* (pH 9)","\u0394 *clsAclsB* att::*clsB* (pH 9)")

colors_plots<-c("#000000","#006198","#56B4E9","#0f682c","#299764","#e38200","#7D4564","#CC79A7",
                "#000000","#006198","#56B4E9","#0f682c","#299764","#e38200","#7D4564","#CC79A7",
                "#000000","#006198","#56B4E9","#0f682c","#299764","#e38200","#7D4564","#CC79A7",
                "#000000","#006198","#56B4E9","#0f682c","#299764","#e38200","#7D4564","#CC79A7",
                "#000000","#006198","#56B4E9","#0f682c","#299764","#e38200","#7D4564","#CC79A7")


breaks_na<-c("WT_a","cA","cA_comp","cB","cB_comp","cAB_a","cAB_acomp","cAB_bcomp")

labels_1<-c("WT",
            "\u0394 *clsA* att:EV","\u0394 *clsA*  att::*clsA*",
            "\u0394 *clsB* att:EV","\u0394 *clsB* att::*clsB*",
            "\u0394 *clsAclsB* att:EV A","\u0394 *clsAclsB* att::*clsA*"
            ,"\u0394 *clsAclsB* att::*clsB*")

colors_1<-c("#c35600","#006198","#56B4E9","#0f682c","#299764","#e38200","#7D4564","#CC79A7")

breaks_2<-c("k030","k300","k450","k600")
labels_2<-c("BHIS","300 mM K","0.002% SDS","0.001% SDS")
colors_2<-c("#c35600","#56B4E9","#0f682c","#e38200")

mediaall_wt.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis","na300") & data_comb.f.avg$strain %in% "WT",],aes(x=Time2,y=od,color=str_med2,fill=str_med2, linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  scale_y_continuous(limits=c(-0.1,1))+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_wt.p)

mediaall_ca.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis","na300") & data_comb.f.avg$strain_class %in% c("WT","cA"),],aes(x=Time2,y=od,color=str_med2,group=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_y_continuous(limits=c(-0.1,1))+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_ca.p)

mediaall_cb.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis","na300") & data_comb.f.avg$strain_class %in% c("WT","cB"),],aes(x=Time2,y=od,color=str_med2,group=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_y_continuous(limits=c(-0.1,1))+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_cb.p)

mediaall_cab.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis","na300") & data_comb.f.avg$strain_class %in% c("WT","cAB"),],aes(x=Time2,y=od,color=str_med2,group=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_y_continuous(limits=c(-0.1,1))+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_cab.p)

ggsave(filename = "figures/growthcurve_na_combined_wt.pdf",plot=mediaall_wt.p,width=2,height = 2,units="in",device=cairo_pdf)
ggsave(filename = "figures/growthcurve_na_combined_clsa.pdf",plot=mediaall_ca.p,width=2,height = 2,units="in",device=cairo_pdf)
ggsave(filename = "figures/growthcurve_na_combined_clsb.pdf",plot=mediaall_cb.p,width=2,height = 2,units="in",device=cairo_pdf)
ggsave(filename = "figures/growthcurve_na_combined_clsaclsb.pdf",plot=mediaall_cab.p,width=2,height = 2,units="in",device=cairo_pdf)

# K
mediaall_wt_k.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis","k300") & data_comb.f.avg$strain %in% "WT",],aes(x=Time2,y=od,color=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_y_continuous(limits=c(-0.1,1))+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_wt_k.p)

mediaall_ca_k.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis","k300") & data_comb.f.avg$strain_class %in% c("WT","cA"),],aes(x=Time2,y=od,color=str_med2,group=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_y_continuous(limits=c(-0.1,1))+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_ca_k.p)

mediaall_cb_k.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis","k300") & data_comb.f.avg$strain_class %in% c("WT","cB"),],aes(x=Time2,y=od,color=str_med2,group=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_y_continuous(limits=c(-0.1,1))+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_cb_k.p)

mediaall_cab_k.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis","k300") & data_comb.f.avg$strain_class %in% c("WT","cAB"),],aes(x=Time2,y=od,color=str_med2,group=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_y_continuous(limits=c(-0.1,1))+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_cab_k.p)

ggsave(filename = "figures/growthcurve_k_combined_wt.pdf",plot=mediaall_wt_k.p,width=2,height = 2,units="in",device=cairo_pdf)
ggsave(filename = "figures/growthcurve_k_combined_clsa.pdf",plot=mediaall_ca_k.p,width=2,height = 2,units="in",device=cairo_pdf)
ggsave(filename = "figures/growthcurve_k_combined_clsb.pdf",plot=mediaall_cb_k.p,width=2,height = 2,units="in",device=cairo_pdf)
ggsave(filename = "figures/growthcurve_k_combined_clsaclsb.pdf",plot=mediaall_cab_k.p,width=2,height = 2,units="in",device=cairo_pdf)

# H (pH 5 and 9)
mediaall_wt_h5.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis","pH5") & data_comb.f.avg$strain %in% "WT",],aes(x=Time2,y=od,color=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_y_continuous(limits=c(-0.1,1))+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_wt_h5.p)

mediaall_wt_h9.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis","pH9") & data_comb.f.avg$strain %in% "WT",],aes(x=Time2,y=od,color=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_y_continuous(limits=c(-0.1,1))+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_wt_h9.p)

mediaall_ca_h5.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis","pH5") & data_comb.f.avg$strain_class %in% c("WT","cA"),],aes(x=Time2,y=od,color=str_med2,group=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_y_continuous(limits=c(-0.1,1))+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_ca_h5.p)

mediaall_ca_h9.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis","pH9") & data_comb.f.avg$strain_class %in% c("WT","cA"),],aes(x=Time2,y=od,color=str_med2,group=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_y_continuous(limits=c(-0.1,1))+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_ca_h9.p)

mediaall_cb_h5.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis","pH5") & data_comb.f.avg$strain_class %in% c("WT","cB"),],aes(x=Time2,y=od,color=str_med2,group=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_y_continuous(limits=c(-0.1,1))+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_cb_h5.p)

mediaall_cb_h9.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis","pH9") & data_comb.f.avg$strain_class %in% c("WT","cB"),],aes(x=Time2,y=od,color=str_med2,group=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_y_continuous(limits=c(-0.1,1))+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_cb_h9.p)

mediaall_cab_h5.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis","pH5") & data_comb.f.avg$strain_class %in% c("WT","cAB"),],aes(x=Time2,y=od,color=str_med2,group=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_y_continuous(limits=c(-0.1,1))+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_cab_h5.p)

mediaall_cab_h9.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis","pH9") & data_comb.f.avg$strain_class %in% c("WT","cAB"),],aes(x=Time2,y=od,color=str_med2,group=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_y_continuous(limits=c(-0.1,1))+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_cab_h9.p)

ggsave(filename = "figures/growthcurve_ph5_combined_wt.pdf",plot=mediaall_wt_h5.p,width=2,height = 2,units="in",device=cairo_pdf)
ggsave(filename = "figures/growthcurve_ph5_combined_clsa.pdf",plot=mediaall_ca_h5.p,width=2,height = 2,units="in",device=cairo_pdf)
ggsave(filename = "figures/growthcurve_ph5_combined_clsb.pdf",plot=mediaall_cb_h5.p,width=2,height = 2,units="in",device=cairo_pdf)
ggsave(filename = "figures/growthcurve_ph5_combined_clsaclsb.pdf",plot=mediaall_cab_h5.p,width=2,height = 2,units="in",device=cairo_pdf)

ggsave(filename = "figures/growthcurve_ph9_combined_wt.pdf",plot=mediaall_wt_h9.p,width=2,height = 2,units="in",device=cairo_pdf)
ggsave(filename = "figures/growthcurve_ph9_combined_clsa.pdf",plot=mediaall_ca_h9.p,width=2,height = 2,units="in",device=cairo_pdf)
ggsave(filename = "figures/growthcurve_ph9_combined_clsb.pdf",plot=mediaall_cb_h9.p,width=2,height = 2,units="in",device=cairo_pdf)
ggsave(filename = "figures/growthcurve_ph9_combined_clsaclsb.pdf",plot=mediaall_cab_h9.p,width=2,height = 2,units="in",device=cairo_pdf)


# Growth curve summary statistics

data_comb.f.avg$str_med2_rep<-paste0(data_comb.f.avg$str_med2,"_",data_comb.f.avg$replicate)
data_growthcurve<-data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis","k300","na300","pH5"),]

data_growthcurve[data_growthcurve$od <0,]$od<-0

data_growthcurve<-data_growthcurve[data_growthcurve$str_med2_rep %!in% c("cA_na300_3"),] # This replicate triggers an error in the fit

out.multi.f.1stpass<-ipg_multisample(data=data_growthcurve,id="str_med2_rep",time.name="Time2",y.name="od",epsilon = 0.2/100)
out.multi.f.1stpass$estimates

# Trimming growth curves for better fit
growth_curve_metrics<-out.multi.f.1stpass$estimates
growth_curve_metrics$`peak growth time2`<-round(growth_curve_metrics$`peak growth time`,0)

growthcurves_nogrowth<-unique(growth_curve_metrics[growth_curve_metrics$`peak growth time` <= 2,]$str_med2_rep)

growthcurve_data.f<-as.data.table(matrix(ncol=ncol(data_growthcurve),dimnames = list(c(1),c(colnames(data_growthcurve)))))
growthcurve_data.f$mintime<-0
growthcurve_data.f$Time3<-0
strmedrep_breaks<-unique(data_growthcurve$str_med2_rep)[unique(data_growthcurve$str_med2_rep) %!in% growthcurves_nogrowth]


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

data_growthcurve.unique<-growthcurve_data.f[!duplicated(growthcurve_data.f$str_med2_rep),c("exp_rep","media_type2","strain","str_med2","str_med2_rep","mintime")]


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


ggsave(filename = "figures/ipolygrowth_mapped_nahk_all.pdf",plot=ipg.plot.all,width=10, height=8,units="in",device=cairo_pdf)

write.csv(x=data_growthcurve.f,file="summarizedgrowth_nahk.csv",row.names = FALSE)

