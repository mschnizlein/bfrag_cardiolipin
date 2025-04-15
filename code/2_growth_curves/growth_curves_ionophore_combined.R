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
ion_0508<-read.csv(file="../ionophores/p207_ionophore_20240508/p207_ionophore_20240508_formatted.csv")
ion_0514<-read.csv(file="../ionophores/p207_ionophore_20240514/p207_ionophore_20240514_formatted.csv")
ion_0515<-read.csv(file="../ionophores/p207_ionophore_20240515/p207_ionophore_20240515_formatted.csv")
ion_0522<-read.csv(file="../ionophores/p207_ionophore_20240522/p207_ionophore_20240522_formatted.csv")

ion_0508$experiment<-"exp0508_ion"
ion_0508$replicate<-1
ion_0514$experiment<-"exp0514_ion"
ion_0514<-ion_0514[ion_0514$abx %!in% "monensin",]
ion_0514$replicate<-2
ion_0515$experiment<-"exp0515_ion"
ion_0515<-ion_0515[ion_0515$abx %!in% "monensin",]
ion_0515$replicate<-3
ion_0522$experiment<-"exp0522_ion"


data_comb<-rbind(ion_0508,ion_0514,ion_0515,ion_0522)

data_comb$media_type2<-data_comb$abx_cond
# data_comb[data_comb$media_type %in% c("k030","bhis","na138","pH7"),]$media_type2<-"bhis"
data_comb.f<-data_comb[data_comb$Time <=20,]
data_comb.f$Time2<-round(data_comb.f$Time, digits=1)
data_comb.f$strain_class<-substrLeft(data_comb.f$strain,3) %>% gsub("_","",.)
data_comb.f<-data_comb.f[!is.na(data_comb.f$strain),]
# data_comb.f[data_comb.f$strain %in% "wt_a",]$strain<-"wt"
# data_comb.f[data_comb.f$strain %in% c("cAB_a","cAB_b"),]$strain<-"cAB"
data_comb.f$str_med2<-paste0(data_comb.f$strain,"_",data_comb.f$media_type2)

data_comb.f<-data_comb.f[data_comb.f$strain %!in% "none",]

str_vector<-unique(data_comb.f$strain)

data_comb.f.avg<-aggregate(od~abx+abx_conc+strain+abx_cond+str_abx_cond+experiment+media_type2+Time2+strain_class+str_med2+replicate,data=data_comb.f,FUN=mean)

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
## "bhis_0" =  BHIS + 1% vehicle (EtOH)
## "monensin_0.4" = BHIS + 0.4 ug/mL monensin A
## "monensin_0.8" = BHIS + 0.8 ug/mL monensin A
## "nigericin_0.1" = BHIS + 0.1 ug/mL nigericin

#####

breaks_plots<-c("wt_bhis_0","cA_bhis_0","cA_comp_bhis_0","cB_bhis_0","cB_comp_bhis_0","cAB_bhis_0","cAB_acomp_bhis_0","cAB_bcomp_bhis_0",
                "wt_monensin_0.4","cA_monensin_0.4","cA_comp_monensin_0.4","cB_monensin_0.4","cB_comp_monensin_0.4","cAB_monensin_0.4","cAB_acomp_monensin_0.4","cAB_bcomp_monensin_0.4",
                "wt_monensin_0.8","cA_monensin_0.8","cA_comp_monensin_0.8","cB_monensin_0.8","cB_comp_monensin_0.8","cAB_monensin_0.8","cAB_acomp_monensin_0.8","cAB_bcomp_monensin_0.8",
                "wt_nigericin_0.1","cA_nigericin_0.1","cA_comp_nigericin_0.1","cB_nigericin_0.1","cB_comp_nigericin_0.1","cAB_nigericin_0.1","cAB_acomp_nigericin_0.1","cAB_bcomp_nigericin_0.1")
labels_plots<-c("WT","\u0394 *clsA*","\u0394 *clsA*  att::*clsA*","\u0394 *clsB*","\u0394 *clsB*  att::*clsB*","\u0394 *clsAclsB*","\u0394 *clsAclsB* att::*clsA*","\u0394 *clsAclsB* att::*clsB*",
                "WT (0.4 ug/mL Monensin)","\u0394 *clsA* (0.4 ug/mL Monensin)","\u0394 *clsA*  att::*clsA* (0.4 ug/mL Monensin)","\u0394 *clsB* (0.4 ug/mL Monensin)","\u0394 *clsB*  att::*clsB* (0.4 ug/mL Monensin)","\u0394 *clsAclsB* (0.4 ug/mL Monensin)","\u0394 *clsAclsB* att::*clsA* (0.4 ug/mL Monensin)","\u0394 *clsAclsB* att::*clsB* (0.4 ug/mL Monensin)",
                "WT (0.8 ug/mL Monensin)","\u0394 *clsA* (0.8 ug/mL Monensin)","\u0394 *clsA*  att::*clsA* (0.8 ug/mL Monensin)","\u0394 *clsB* (0.8 ug/mL Monensin)","\u0394 *clsB*  att::*clsB* (0.8 ug/mL Monensin)","\u0394 *clsAclsB* (0.8 ug/mL Monensin)","\u0394 *clsAclsB* att::*clsB* (0.8 ug/mL Monensin)","\u0394 *clsAclsB* att::*clsB* (0.8 ug/mL Monensin)",
                "WT (0.1 ug/mL Nigericin)","\u0394 *clsA* (0.1 ug/mL Nigericin)","\u0394 *clsA*  att::*clsA* (0.1 ug/mL Nigericin)","\u0394 *clsB* (0.1 ug/mL Nigericin)","\u0394 *clsB*  att::*clsB* (0.1 ug/mL Nigericin)","\u0394 *clsAclsB* (0.1 ug/mL Nigericin)","\u0394 *clsAclsB* att::*clsB* (0.1 ug/mL Nigericin)","\u0394 *clsAclsB* att::*clsB* (0.1 ug/mL Nigericin)")

# colors_plots<-c("#000000","#424242","#999999","#424242","#999999","#333333","#636363","#999999",
#                "#c35600","#006198","#56B4E9","#0f682c","#299764","#e38200","#46c8a1","#69fae2",
#                "#c35600","#006198","#56B4E9","#0f682c","#299764","#e38200","#46c8a1","#69fae2",
#                "#c35600","#006198","#56B4E9","#0f682c","#299764","#e38200","#46c8a1","#69fae2")

colors_plots<-c("#000000","#006198","#56B4E9","#0f682c","#299764","#e38200","#7D4564","#CC79A7",
                "#000000","#006198","#56B4E9","#0f682c","#299764","#e38200","#7D4564","#CC79A7",
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

# 0.4 ug/mL Monensin A
mediaall_wt_mon4.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type %in% c("bhis_0","monensin_0.4") & data_comb.f.avg$strain %in% "wt",],aes(x=Time2,y=od,color=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_y_continuous(limits=c(0,1))+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_wt_mon4.p)

mediaall_ca_mon4.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis_0","monensin_0.4") & data_comb.f.avg$strain_class %in% c("wt","cA"),],aes(x=Time2,y=od,color=str_med2,group=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_y_continuous(limits=c(0,1))+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_ca_mon4.p)

mediaall_ca_mon4_legend.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis_0","monensin_0.4") & data_comb.f.avg$strain_class %in% c("wt","cA"),],aes(x=Time2,y=od,color=str_med2,group=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_y_continuous(limits=c(0,1))+
  as_md_theme(my_theme)+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_ca_mon4_legend.p)

mediaall_cb_mon4.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis_0","monensin_0.4") & data_comb.f.avg$strain_class %in% c("wt","cB"),],aes(x=Time2,y=od,color=str_med2,group=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_y_continuous(limits=c(0,1))+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_cb_mon4.p)

mediaall_cb_mon4.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis_0","monensin_0.4") & data_comb.f.avg$strain_class %in% c("wt","cB"),],aes(x=Time2,y=od,color=str_med2,group=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_y_continuous(limits=c(0,1))+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_cb_mon4.p)

mediaall_cab_mon4.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis_0","monensin_0.4") & data_comb.f.avg$strain_class %in% c("wt","cAB"),],aes(x=Time2,y=od,color=str_med2,group=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_y_continuous(limits=c(0,1))+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_cab_mon4.p)

mediaall_cab_mon4.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis_0","monensin_0.4") & data_comb.f.avg$strain_class %in% c("wt","cAB"),],aes(x=Time2,y=od,color=str_med2,group=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_y_continuous(limits=c(0,1))+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_cab_mon4.p)

ggsave(filename = "figures/growthcurve_monensin4_combined_wt.pdf",plot=mediaall_wt_mon4.p,width=2,height = 2,units="in",device=cairo_pdf)
ggsave(filename = "figures/growthcurve_monensin4_combined_clsa.pdf",plot=mediaall_ca_mon4.p,width=2,height = 2,units="in",device=cairo_pdf)
ggsave(filename = "figures/growthcurve_monensin4_combined_clsb.pdf",plot=mediaall_cb_mon4.p,width=2,height = 2,units="in",device=cairo_pdf)
ggsave(filename = "figures/growthcurve_monensin4_combined_clsaclsb.pdf",plot=mediaall_cab_mon4.p,width=2,height = 2,units="in",device=cairo_pdf)

# 0.8 ug/mL monensin A
mediaall_wt_mon8.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type %in% c("bhis_0","monensin_0.8") & data_comb.f.avg$strain %in% "wt",],aes(x=Time2,y=od,color=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_y_continuous(limits=c(0,1))+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_wt_mon8.p)

mediaall_ca_mon8.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis_0","monensin_0.8") & data_comb.f.avg$strain_class %in% c("wt","cA"),],aes(x=Time2,y=od,color=str_med2,group=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_y_continuous(limits=c(0,1))+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_ca_mon8.p)

mediaall_cb_mon8.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis_0","monensin_0.8") & data_comb.f.avg$strain_class %in% c("wt","cB"),],aes(x=Time2,y=od,color=str_med2,group=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_y_continuous(limits=c(0,1))+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_cb_mon8.p)

mediaall_cab_mon8.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis_0","monensin_0.8") & data_comb.f.avg$strain_class %in% c("wt","cAB"),],aes(x=Time2,y=od,color=str_med2,group=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_y_continuous(limits=c(0,1))+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_cab_mon8.p)

ggsave(filename = "figures/growthcurve_monensin8_combined_wt.pdf",plot=mediaall_wt_mon8.p,width=2,height = 2,units="in",device=cairo_pdf)
ggsave(filename = "figures/growthcurve_monensin8_combined_clsa.pdf",plot=mediaall_ca_mon8.p,width=2,height = 2,units="in",device=cairo_pdf)
ggsave(filename = "figures/growthcurve_monensin8_combined_clsb.pdf",plot=mediaall_cb_mon8.p,width=2,height = 2,units="in",device=cairo_pdf)
ggsave(filename = "figures/growthcurve_monensin8_combined_clsaclsb.pdf",plot=mediaall_cab_mon8.p,width=2,height = 2,units="in",device=cairo_pdf)

# 0.1 ug/mL Nigericin
mediaall_wt_nig1.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type %in% c("bhis_0","nigericin_0.1") & data_comb.f.avg$strain %in% "wt",],aes(x=Time2,y=od,color=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_y_continuous(limits=c(0,1))+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_wt_nig1.p)

mediaall_ca_nig1.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis_0","nigericin_0.1") & data_comb.f.avg$strain_class %in% c("wt","cA"),],aes(x=Time2,y=od,color=str_med2,group=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_y_continuous(limits=c(0,1))+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_ca_nig1.p)

mediaall_cb_nig1.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis_0","nigericin_0.1") & data_comb.f.avg$strain_class %in% c("wt","cB"),],aes(x=Time2,y=od,color=str_med2,group=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_y_continuous(limits=c(0,1))+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_cb_nig1.p)

mediaall_cab_nig1.p<-ggplot(data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis_0","nigericin_0.1") & data_comb.f.avg$strain_class %in% c("wt","cAB"),],aes(x=Time2,y=od,color=str_med2,group=str_med2,fill=str_med2,linetype=media_type2))+
  stat_summary(fun=mean, geom="line", linewidth=1)+
  stat_summary(fun.data = mean_se,geom="ribbon",color="transparent",alpha=0.3)+
  scale_color_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_fill_manual(values=colors_plots,breaks=breaks_plots,labels=labels_plots)+
  scale_y_continuous(limits=c(0,1))+
  as_md_theme(my_theme)+theme(legend.position="none")+
  ylab("OD600") + xlab("Time (Hrs)") + labs(color="",fill="") + ggtitle("")
plot(mediaall_cab_nig1.p)

ggsave(filename = "figures/growthcurve_nigericin1_combined_wt.pdf",plot=mediaall_wt_nig1.p,width=2,height = 2,units="in",device=cairo_pdf)
ggsave(filename = "figures/growthcurve_nigericin1_combined_clsa.pdf",plot=mediaall_ca_nig1.p,width=2,height = 2,units="in",device=cairo_pdf)
ggsave(filename = "figures/growthcurve_nigericin1_combined_clsb.pdf",plot=mediaall_cb_nig1.p,width=2,height = 2,units="in",device=cairo_pdf)
ggsave(filename = "figures/growthcurve_nigericin1_combined_clsaclsb.pdf",plot=mediaall_cab_nig1.p,width=2,height = 2,units="in",device=cairo_pdf)

# Growth curve summary statistics (ipolygrowth)
# data_growthcurve<-data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis_0","monensin_0.8","nigericin_0.1"),]
# data_comb.f.avg.single<-aggregate(od~strain+media_type2+Time2+strain_class+str_med2,data=data_comb.f.avg,FUN=mean)

data_comb.f.avg$str_med2_rep<-paste0(data_comb.f.avg$str_med2,"_",data_comb.f.avg$replicate)
data_growthcurve<-data_comb.f.avg[data_comb.f.avg$media_type2 %in% c("bhis_0","monensin_0.8","nigericin_0.1"),]

# data_growthcurve.bhis<-data_growthcurve[data_growthcurve$media_type2 %in% "bhis" & data_growthcurve$Time2 >=0.5 & data_growthcurve$Time2 <= 14,]
# data_growthcurve.lagged<-data_growthcurve[data_growthcurve$media_type2 %!in% "bhis" & data_growthcurve$Time2 >=0.5 & data_growthcurve$Time2 <= 20,]
data_growthcurve[data_growthcurve$od <0,]$od<-0

out.multi.f.1stpass<-ipg_multisample(data=data_growthcurve,id="str_med2_rep",time.name="Time2",y.name="od",epsilon = 0.2/100)
out.multi.f.1stpass$estimates

# Trimming growth curves for better fit
growth_curve_metrics<-out.multi.f.1stpass$estimates
growth_curve_metrics$`peak growth time2`<-round(growth_curve_metrics$`peak growth time`,0)

growthcurves_nogrowth<-unique(growth_curve_metrics[growth_curve_metrics$`peak growth time` <= 2,]$str_med2_rep)

growthcurve_data.f<-as.data.table(matrix(ncol=13,dimnames = list(c(1),c(colnames(data_growthcurve)))))
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


ggsave(filename = "figures/ipolygrowth_mapped_ionophore_all.pdf",plot=ipg.plot.all,width=10, height=8,units="in",device=cairo_pdf)

write.csv(x=data_growthcurve.f,file="summarizedgrowth_ionophore.csv",row.names = FALSE)
