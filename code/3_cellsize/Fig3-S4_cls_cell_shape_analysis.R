## P207 cell shape analysis
## Author: Matthew Schnizlein
## 10/18/2024
## Fig. 3 & S4
### Contents
## Analysis of cell shape of B. fragilis P207 cls deletion strains
# Fig. S4B - Kernel probability distributions of cell lengths
# Fig. S4C - Kernel probability distributions of cell widths
# Fig. 3B - dot plot of cell widths by cls strain
# Statistical comparisons
# Fig. 3A - IQR of probability distributions found in S4B-C

library(ggplot2)
# library(reshape2)
library(growthcurver)
library(ggtext)
library(gridtext)
# library(glue)
library(mdthemes)
library(matrixStats)
library(data.table)
library(data.table)
library(magrittr)
library(readxl)
library(scales)
library(colorBlindness)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyr)
library(ggbiplot)
library(ggforce)
library(lme4)
library(spatstat)

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

breaks_1<-c("wt","ca","ca_comp","cb","cb_comp","cab","cab_acomp","cab_bcomp")
colors_1<-c("#000000","#006198","#56B4E9","#0f682c","#299764","#e38200","#7D4564","#CC79A7")
labels_1<-c("WT","\u0394*clsA*","\u0394*clsA* att::*clsA*","\u0394*clsB*","\u0394*clsB* att::*clsB*","\u0394*clsA\u0394clsB*","\u0394*clsA\u0394clsB* att::*clsA*","\u0394*clsA\u0394clsB* att::*clsB*")

setwd("~/data/3_cellsize")

####
exp0322<-read.csv(file="revised_20240322.csv",header = TRUE)
exp0329<-read.csv(file="revised_20240329.csv",header = TRUE)
exp0412<-read.csv(file="revised_20240412.csv",header = TRUE)
exp0418<-read.csv(file="revised_20240418.csv",header = TRUE)

celldata<-rbind(exp0322,exp0329,exp0412,exp0418)
celldata[is.na(celldata$association),]$image.label

celldata<-celldata[!is.na(celldata$association),]

celldata$image_num<-gsub(".*_Image","",celldata$image.label)
celldata$strain<-gsub(".*_2024","",celldata$image.label) %>% gsub("_Image.*","",.) %>% gsub("0322_","",.) %>% gsub("0329_","",.) %>% gsub("0412_","",.) %>% gsub("0418_","",.) %>% gsub("001","",.) %>%
  gsub("cab_bev","cab_b",.) %>% gsub("cab_aev","cab_a",.) %>% gsub("cbcomp","cb_comp",.) %>% gsub("cacomp","ca_comp",.) %>% gsub("ca_ev","ca",.) %>% gsub("cb_ev","cb",.)
celldata$strain2<-celldata$strain
celldata[celldata$strain %in% c("wt_a","wt_b","wt_a_2"),]$strain2<-"wt"
celldata[celldata$strain %in% c("cab_a","cab_b"),]$strain2<-"cab"

######
ggplot(celldata,aes(x=intensity,color=image))+geom_histogram(bins=200)+scale_x_log10()
ggplot(celldata,aes(x=intensity,y=shape.length,color=image))+geom_point()+scale_x_log10()
ggplot(celldata,aes(x=intensity,y=shape.width.mean,color=image))+geom_point()+scale_x_log10()

intensity.length.lm<-lm(log10(intensity)~shape.length,data=celldata)
summary(intensity.length.lm)
intensity.width.lm<-lm(log10(intensity)~shape.width.mean,data=celldata)
summary(intensity.width.lm)

kernel.intensityall.p<-ggplot(celldata, aes(x=log10(intensity),color=image))+geom_density(kernel="epanechnikov",bw="nrd",linewidth=1,position = "identity",fill=NA)+
  scale_x_continuous(name = "log10(Intensity)")+
  scale_y_continuous(name="Kernel Density")+
  as_md_theme(my_theme)
plot(kernel.intensityall.p)

kernel.widthall.p<-ggplot(celldata, aes(x=shape.width.mean,color=image))+geom_density(kernel="epanechnikov",bw="nrd",linewidth=1,position = "identity",fill=NA)+
  scale_x_continuous(limits=c(0,2),name = "Cell Width (\U00B5m)")+
  scale_y_continuous(name="Kernel Density")+
  as_md_theme(my_theme)
plot(kernel.widthall.p)
# Some experiments had higher intensity, controlling for this variability will be important

######
summary_stats_names<-c("intensity","shape.roundness","shape.length","shape.width.mean","shape.width.max","shape.width.min")
summary_stats_strains<-c("wt","ca","ca_comp","cb","cb_comp","cab","cab_acomp","cab_bcomp")
summary_stats_names_f<-vector()

for (i in 1:8){
  names_initial<-paste0(summary_stats_strains[i],"_",summary_stats_names)
  summary_stats_names_f<-c(summary_stats_names_f,names_initial)
}

summary_stats<-as.data.table(matrix(data=NA,nrow=48,ncol=8))
colnames(summary_stats)<-c("id","strain","metric","mean","median","sd","min","max")
summary_stats[,2]<-rep(summary_stats_strains,each=6)
summary_stats[,3]<-rep(summary_stats_names,times=8)

for (i in 1:48){
  summary_stats[[i,"id"]]<-summary_stats_names_f[i]
  summary_stats[[i,"mean"]]<-summary(celldata[celldata$strain2 %in% summary_stats[[i,2]],summary_stats[[i,3]]])[4]
  summary_stats[[i,"median"]]<-summary(celldata[celldata$strain2 %in% summary_stats[[i,2]],summary_stats[[i,3]]])[3]
  summary_stats[[i,"sd"]]<-sd(celldata[celldata$strain2 %in% summary_stats[[i,2]],summary_stats[[i,3]]])
  summary_stats[[i,"min"]]<-summary(celldata[celldata$strain2 %in% summary_stats[[i,2]],summary_stats[[i,3]]])[1]
  summary_stats[[i,"max"]]<-summary(celldata[celldata$strain2 %in% summary_stats[[i,2]],summary_stats[[i,3]]])[6]
}

ggplot(summary_stats[summary_stats$metric %in% "intensity",])+geom_point(aes(x=id,y=mean,color="mean"))+geom_point(aes(x=id,y=median,color="median"))+scale_color_discrete()
ggplot(summary_stats[summary_stats$metric %!in% c("intensity","shape.length","shape.roundness"),])+geom_point(aes(x=id,y=mean,color="mean"))+geom_point(aes(x=id,y=median,color="median"))+scale_color_discrete()

celldata$excluded<-FALSE
celldata[celldata$strain %in% "wt_a" & celldata$image %in% "p207_cellshape_20240322",]<-TRUE
celldata<-celldata[celldata$excluded %in% FALSE,]

cellcount_stats<-as.data.table(matrix(data=NA,nrow=8,ncol=9))
colnames(cellcount_stats)<-c("strain","cell_num","photo_num","cell_length_avg","cell_width_avg","cell_roundness_avg","cell_length_sd","cell_width_sd","cell_roundness_sd")

for (i in 1:8){
  celldata.sub.temp<-celldata[celldata$strain2 %in% summary_stats_strains[i],]
  cellcount_stats[[i,"strain"]]<-summary_stats_strains[i]
  cellcount_stats[[i,"cell_num"]]<-nrow(celldata.sub.temp)
  cellcount_stats[[i,"photo_num"]]<-length(unique(celldata.sub.temp$image_num))
  cellcount_stats[[i,"cell_length_avg"]]<-mean(celldata.sub.temp$shape.length)
  cellcount_stats[[i,"cell_width_avg"]]<-mean(celldata.sub.temp$shape.width.mean)
  cellcount_stats[[i,"cell_roundness_avg"]]<-mean(celldata.sub.temp$shape.roundness)
  cellcount_stats[[i,"cell_length_sd"]]<-sd(celldata.sub.temp$shape.length)
  cellcount_stats[[i,"cell_width_sd"]]<-sd(celldata.sub.temp$shape.width.mean)
  cellcount_stats[[i,"cell_roundness_sd"]]<-sd(celldata.sub.temp$shape.roundness)
}

summary(cellcount_stats$cell_num)


for (i in 1:8){
  subsetted.data<-celldata[celldata$strain2 %in% summary_stats_strains[i],]
  subsetted.data.sample<-subsetted.data[sample(nrow(subsetted.data),size=323,replace=FALSE),]
  assign(x = paste0(summary_stats_strains[i],"_sample"),subsetted.data.sample,envir = as.environment(.GlobalEnv))
}

data_cv_sampled<-rbind(wt_sample,ca_sample,ca_comp_sample,cb_sample,cb_comp_sample,cab_sample,cab_acomp_sample,cab_bcomp_sample)
data_cv_sampled$strain2<-factor(data_cv_sampled$strain2,levels=c("wt","ca","ca_comp","cb","cb_comp","cab","cab_acomp","cab_bcomp" ))

# For cell volume calculation:	Loferer-Krössbacher M, Klima J, Psenner R. Determination of bacterial cell dry mass by transmission electron microscopy and densitometric image analysis. Appl Environ Microbiol. 1998;64(2):688-94. Epub 1998/02/17. doi: 10.1128/aem.64.2.688-694.1998. PubMed PMID: 9464409; PMCID: PMC106103.
## V = [(w^2 · π/4) · (l − w)] + (π · w^3/6)
## where V is volume, w is width and l is length

data_cv_sampled$shape_volume<-((data_cv_sampled$shape.width.mean^2 * pi/4)*
                                 (data_cv_sampled$shape.length-
                                    data_cv_sampled$shape.width.mean))+
  (pi * (data_cv_sampled$shape.width.mean^3)/6)

ggplot(data_cv_sampled,aes(x=shape_volume))+geom_histogram()

## Kernal by length (Fig. S4B)
kernel.lengthall_ca.p<-ggplot(data_cv_sampled[data_cv_sampled$strain2 %in% c("wt","ca","ca_comp"),], aes(x=shape.length,color=strain2))+geom_density(kernel="epanechnikov",bw="nrd",linewidth=0.5,position = "identity",fill=NA)+
  scale_x_continuous(limits=c(1,5),name = "Cell Length (\U00B5m)")+
  scale_y_continuous(limits=c(0,0.8),name="Kernel Density")+
  scale_color_manual(values=colors_1,breaks=breaks_1)+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(kernel.lengthall_ca.p)

kernel.lengthall_cb.p<-ggplot(data_cv_sampled[data_cv_sampled$strain2 %in% c("wt","cb","cb_comp"),], aes(x=shape.length,color=strain2))+geom_density(kernel="epanechnikov",bw="nrd",linewidth=0.5,position = "identity",fill=NA)+
  scale_x_continuous(limits=c(1,5),name = "Cell Length (\U00B5m)")+
  scale_y_continuous(limits=c(0,0.8),name="Kernel Density")+
  scale_color_manual(values=colors_1,breaks=breaks_1)+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(kernel.lengthall_cb.p)

kernel.lengthall_cab.p<-ggplot(data_cv_sampled[data_cv_sampled$strain2 %in% c("wt","cab","cab_acomp","cab_bcomp"),], aes(x=shape.length,color=strain2))+geom_density(kernel="epanechnikov",bw="nrd",linewidth=0.5,position = "identity",fill=NA)+
  scale_x_continuous(limits=c(1,5),name = "Cell Length (\U00B5m)")+
  scale_y_continuous(limits=c(0,0.8),name="Kernel Density")+
  scale_color_manual(values=colors_1,breaks=breaks_1)+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(kernel.lengthall_cab.p)

## Kernal by Width (Fig. S4C)
kernel.widthall_ca.p<-ggplot(data_cv_sampled[data_cv_sampled$strain2 %in% c("wt","ca","ca_comp"),], aes(x=shape.width.mean,color=strain2))+geom_density(kernel="epanechnikov",bw="nrd",linewidth=0.5,position = "identity",fill=NA)+
  scale_x_continuous(limits=c(0.5,1),name = "Cell Width (\U00B5m)")+
  scale_y_continuous(limits=c(0,7),name="Kernel Density")+
  scale_color_manual(values=colors_1,breaks=breaks_1)+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(kernel.widthall_ca.p)

kernel.widthall_cb.p<-ggplot(data_cv_sampled[data_cv_sampled$strain2 %in% c("wt","cb","cb_comp"),], aes(x=shape.width.mean,color=strain2))+geom_density(kernel="epanechnikov",bw="nrd",linewidth=0.5,position = "identity",fill=NA)+
  scale_x_continuous(limits=c(0.5,1),name = "Cell Width (\U00B5m)")+
  scale_y_continuous(limits=c(0,7),name="Kernel Density")+
  scale_color_manual(values=colors_1,breaks=breaks_1)+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(kernel.widthall_cb.p)

kernel.widthall_cab.p<-ggplot(data_cv_sampled[data_cv_sampled$strain2 %in% c("wt","cab","cab_acomp","cab_bcomp"),], aes(x=shape.width.mean,color=strain2))+geom_density(kernel="epanechnikov",bw="nrd",linewidth=0.5,position = "identity",fill=NA)+
  scale_x_continuous(limits=c(0.5,1),name = "Cell Width (\U00B5m)")+
  scale_y_continuous(limits=c(0,7),name="Kernel Density")+
  scale_color_manual(values=colors_1,breaks=breaks_1)+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(kernel.widthall_cab.p)

ggsave(filename = "figures/kernel_intensity_all.pdf",plot=kernel.intensityall.p,height=6,width=8,units="in",useDingbats=FALSE)
ggsave(filename = "figures/kernel_width_all_ca.pdf",plot=kernel.widthall_ca.p,height=2,width=2,units="in",device=cairo_pdf)
ggsave(filename = "figures/kernel_width_all_cb.pdf",plot=kernel.widthall_cb.p,height=2,width=2,units="in",device=cairo_pdf)
ggsave(filename = "figures/kernel_width_all_cab.pdf",plot=kernel.widthall_cab.p,height=2,width=2,units="in",device=cairo_pdf)
ggsave(filename = "figures/kernel_length_ca.pdf",plot=kernel.lengthall_ca.p,height=2,width=2,units="in",device=cairo_pdf)
ggsave(filename = "figures/kernel_length_cb.pdf",plot=kernel.lengthall_cb.p,height=2,width=2,units="in",device=cairo_pdf)
ggsave(filename = "figures/kernel_length_cab.pdf",plot=kernel.lengthall_cab.p,height=2,width=2,units="in",device=cairo_pdf)

######


cv_length.p<-ggplot(data_cv_sampled,aes(x=strain2,y=shape.length,color=strain2))+geom_violin()+
  stat_summary(fun=mean,fun.min=mean,fun.max=mean,geom="crossbar",width=0.5,color="black") +
  scale_y_continuous(name="Cell Length (\u03bcm)")+
  scale_x_discrete(name="Strain",limits=c("wt","ca","ca_comp","cb","cb_comp","cab","cab_acomp","cab_bcomp"),labels=c("WT","\u0394*clsA*","\u0394*clsA* att::*clsA*","\u0394*clsB*","\u0394*clsB* att::*clsB*","\u0394*clsA\u0394clsB*","\u0394*clsA\u0394clsB* att::*clsA*","\u0394*clsA\u0394clsB* att::*clsB*"))+
  scale_color_manual(values=colors_1,breaks=breaks_1)+
  as_md_theme(my_theme2)+theme(legend.position="none")
plot(cv_length.p)

ggsave(filename="figures/cell_length_plot.pdf",plot=cv_length.p,width = 3,height=3,units="in",device=cairo_pdf)

## Fig. 3B
cv_width.p<-ggplot(data_cv_sampled,aes(x=strain2,y=shape.width.mean,color=strain2))+geom_jitter(height=0,width=0.3,alpha=0.2,stroke=NA)+
  stat_summary(fun=mean,fun.min=mean,fun.max=mean,geom="crossbar",width=0.5,color="black") +
  scale_y_continuous(name="Cell Width (\u03bcm)")+
  scale_x_discrete(name="",limits=rev(c("wt","ca","ca_comp","cb","cb_comp","cab","cab_acomp","cab_bcomp")),labels=rev(c("WT","\u0394*clsA*","\u0394*clsA* att::*clsA*","\u0394*clsB*","\u0394*clsB* att::*clsB*","\u0394*clsA\u0394clsB*","\u0394*clsA\u0394clsB* att::*clsA*","\u0394*clsA\u0394clsB* att::*clsB*")))+
  scale_color_manual(values=colors_1,breaks=breaks_1)+
  as_md_theme(my_theme)+theme(legend.position="none")+
  coord_flip()
plot(cv_width.p)


ggsave(filename="figures/cell_width_plot.pdf",plot=cv_width.p,width =3.5,height=2,units="in",device=cairo_pdf)

cv_volume.p<-ggplot(data_cv_sampled,aes(x=strain2,y=shape_volume,color=strain2))+geom_jitter(height=0,width=0.3,alpha=0.2)+
  stat_summary(fun=mean,fun.min=mean,fun.max=mean,geom="crossbar",width=0.5,color="black") +
  scale_y_continuous(limits=c(0.25,2.5),name="Cell Area (\u03bcm^2)")+
  scale_x_discrete(name="Strain",limits=c("wt","ca","ca_comp","cb","cb_comp","cab","cab_acomp","cab_bcomp"),breaks=c("wt","ca","ca_comp","cb","cb_comp","cab","cab_acomp","cab_bcomp"),labels=c("WT","\u0394*clsA*","\u0394*clsA* att::*clsA*","\u0394*clsB*","\u0394*clsB* att::*clsB*","\u0394*clsA\u0394clsB*","\u0394*clsA\u0394clsB* att::*clsA*","\u0394*clsA\u0394clsB* att::*clsB*"))+
  scale_color_manual(values=colors_1,breaks=breaks_1)+
  as_md_theme(my_theme2)+theme(legend.position="none")
plot(cv_volume.p)

ggsave(filename="figures/cell_vol_plot.pdf",plot=cv_volume.p,width = 3,height=3,units="in",device=cairo_pdf)

# By experiment
cv_width_byexp.p<-ggplot(data_cv_sampled,aes(x=strain2,y=shape.width.mean,color=strain2))+geom_jitter(height=0,width=0.3,alpha=0.2)+
  stat_summary(fun=mean,fun.min=mean,fun.max=mean,geom="crossbar",width=0.5,color="black") +
  scale_y_continuous(name="Cell Width (\u03bcm)")+
  scale_x_discrete(name="Strain",limits=c("wt","ca","ca_comp","cb","cb_comp","cab","cab_acomp","cab_bcomp"),labels=c("WT","\u0394*clsA*","\u0394*clsA* att::*clsA*","\u0394*clsB*","\u0394*clsB* att::*clsB*","\u0394*clsA\u0394clsB*","\u0394*clsA\u0394clsB* att::*clsA*","\u0394*clsA\u0394clsB* att::*clsB*"))+
  scale_color_manual(values=colors_1,breaks=breaks_1)+
  as_md_theme(my_theme2)+theme(legend.position="none")+facet_wrap(~image)
plot(cv_width_byexp.p)

ggsave(filename="figures/cell_width_plot_byexp.pdf",plot=cv_width_byexp.p,width =6,height=6,units="in",device=cairo_pdf)

# Volume by cell strain

summary(data_cv_sampled[data_cv_sampled$strain2 %in% "wt",]$shape_volume) # 0.999 um^3
summary(data_cv_sampled[data_cv_sampled$strain2 %in% "ca",]$shape_volume) # 0.918 um^3
summary(data_cv_sampled[data_cv_sampled$strain2 %in% "ca_comp",]$shape_volume) # 0.975 um^3
summary(data_cv_sampled[data_cv_sampled$strain2 %in% "cb",]$shape_volume) # 0.985 um^3
summary(data_cv_sampled[data_cv_sampled$strain2 %in% "cb_comp",]$shape_volume) # 0.963 um^3
summary(data_cv_sampled[data_cv_sampled$strain2 %in% "cab",]$shape_volume) # 0.871 um^3
summary(data_cv_sampled[data_cv_sampled$strain2 %in% "cab_acomp",]$shape_volume) # 0.951 um^3
summary(data_cv_sampled[data_cv_sampled$strain2 %in% "cab_bcomp",]$shape_volume) # 0.902 um^3

cellvol<-as.data.table(matrix(data="NA",nrow=8,ncol=2))
colnames(cellvol)<-c("strain","cell_vol")
cellvol$strain<-breaks_strain
cellvol$cell_vol<-c(0.999,0.918,0.975,0.985,0.963,0.871,0.951,0.902)

write.csv(file = "cell_volume_mean.csv",x = cellvol,row.names = FALSE)


### stats

width.mlm<-lmerTest::lmer(log10(shape.width.mean)~strain2+(1|image),data=data_cv_sampled)
summary(width.mlm)
plot(width.mlm)

hist(log10(data_cv_sampled$shape.width.mean))
hist(summary(width.mlm)$residuals)
shapiro.test(summary(width.mlm)$residuals)

ks.test(log10(data_cv_sampled[data_cv_sampled$strain2 %in% "wt",]$shape.width.mean),log10(data_cv_sampled[data_cv_sampled$strain2 %in% "cb",]$shape.width.mean))
ks.test(log10(data_cv_sampled[data_cv_sampled$strain2 %in% "cab_bcomp",]$shape.width.mean),log10(data_cv_sampled[data_cv_sampled$strain2 %in% "wt",]$shape.width.mean))
# t.test(log10(data_cv_sampled[data_cv_sampled$strain2 %in% "cab_bcomp",]$shape.width.mean),log10(data_cv_sampled[data_cv_sampled$strain2 %in% "wt",]$shape.width.mean))

var.test(log10(data_cv_sampled[data_cv_sampled$strain2 %in% "cb",]$shape.width.mean),log10(data_cv_sampled[data_cv_sampled$strain2 %in% "wt",]$shape.width.mean))
var.test(log10(data_cv_sampled[data_cv_sampled$strain2 %in% "ca",]$shape.width.mean),log10(data_cv_sampled[data_cv_sampled$strain2 %in% "wt",]$shape.width.mean))
var.test(log10(data_cv_sampled[data_cv_sampled$strain2 %in% "cb_comp",]$shape.width.mean),log10(data_cv_sampled[data_cv_sampled$strain2 %in% "wt",]$shape.width.mean))
var.test(log10(data_cv_sampled[data_cv_sampled$strain2 %in% "cb_comp",]$shape.width.mean),log10(data_cv_sampled[data_cv_sampled$strain2 %in% "cb",]$shape.width.mean))
var.test(log10(data_cv_sampled[data_cv_sampled$strain2 %in% "cab_acomp",]$shape.width.mean),log10(data_cv_sampled[data_cv_sampled$strain2 %in% "wt",]$shape.width.mean))
var.test(log10(data_cv_sampled[data_cv_sampled$strain2 %in% "cab",]$shape.width.mean),log10(data_cv_sampled[data_cv_sampled$strain2 %in% "wt",]$shape.width.mean))

# Residuals are non-normal but observations per group are >300 -> going to assume normality

## Random effects:
# Groups   Name        Variance  Std.Dev.
# image    (Intercept) 0.0003245 0.01802
# Residual             0.0017437 0.04176
# Number of obs: 2584, groups:  image, 4

## Fixed effects:
#   Estimate Std. Error         df t value Pr(>|t|)
#   (Intercept)      -1.207e-01  9.311e-03  3.366e+00 -12.967 0.000544 ***
#   strain2ca        -2.056e-02  3.304e-03  2.573e+03  -6.225 5.62e-10 ***
#   strain2ca_comp   -4.245e-03  3.320e-03  2.574e+03  -1.278 0.201220
#   strain2cb        -3.502e-03  3.291e-03  2.573e+03  -1.064 0.287413
#   strain2cb_comp   -6.977e-03  3.314e-03  2.574e+03  -2.105 0.035371 *
#   strain2cab       -1.925e-02  3.330e-03  2.574e+03  -5.780 8.35e-09 ***
#   strain2cab_acomp -1.114e-02  3.295e-03  2.573e+03  -3.382 0.000730 ***
#   strain2cab_bcomp -2.478e-02  3.334e-03  2.574e+03  -7.433 1.44e-13 ***

data_cv_sampled$strain3<-factor(data_cv_sampled$strain2,levels=c("cab","wt","ca","ca_comp","cb","cb_comp","cab_acomp","cab_bcomp"))
width.mlm2<-lmerTest::lmer(log10(shape.width.mean)~strain3+(1|image),data=data_cv_sampled)
summary(width.mlm2)

## Random effects:
#  Groups   Name        Variance  Std.Dev.
# image    (Intercept) 0.0003245 0.01802
# Residual             0.0017437 0.04176
# Number of obs: 2584, groups:  image, 4

## Fixed effects:
#  Estimate Std. Error         df t value Pr(>|t|)
# (Intercept)      -1.400e-01  9.304e-03  3.355e+00 -15.047 0.000339 ***
#  strain3wt         1.925e-02  3.330e-03  2.574e+03   5.780 8.35e-09 ***
#  strain3ca        -1.317e-03  3.301e-03  2.573e+03  -0.399 0.690016
#  strain3ca_comp    1.500e-02  3.305e-03  2.574e+03   4.539 5.92e-06 ***
#  strain3cb         1.574e-02  3.315e-03  2.574e+03   4.749 2.15e-06 ***
#  strain3cb_comp    1.227e-02  3.289e-03  2.573e+03   3.730 0.000195 ***
#  strain3cab_acomp  8.105e-03  3.305e-03  2.573e+03   2.453 0.014252 *
#  strain3cab_bcomp -5.535e-03  3.321e-03  2.574e+03  -1.667 0.095650 .


# Calculating IQR on Data and Kernal Density Distribution (Fig. 3A)

fivenum.width<-as.data.frame(matrix(data=NA,nrow=8,ncol=8),)
colnames(fivenum.width)<-c("strain","quartile1","quartile2","quartile3","kernal_quartile1","kernal_quartile2","kernal_quartile3","kernal_mode")
i<-1
for (i in 1:8){
  fivenum.width[i,1]<-summary_stats_strains[i]
  data_cv_sampled.sub<-data_cv_sampled[data_cv_sampled$strain2 %in% summary_stats_strains[i],]

  fivenum.width[i,2:4]<-fivenum(data_cv_sampled.sub$shape.width.mean)[2:4]

  width_density.temp<-density(data_cv_sampled.sub$shape.width.mean,kernel=c("epanechnikov"))
  fivenum.width[i,5:7]<-quantile.density(width_density.temp)[2:4]
  fivenum.width[i,8]<-width_density.temp$x[which.max(width_density.temp$y)]

}

fivenum.width$kernal_range_q1q3<-fivenum.width$kernal_quartile3-fivenum.width$kernal_quartile1
fivenum.width$kernal_range_q1q3_percwt<-(fivenum.width$kernal_range_q1q3/0.08457507
)*100
fivenum.width[,2:10]<-round(fivenum.width[,2:10],digits=3)


ggplot(fivenum.width)+geom_pointrange(aes(x=strain,ymin=quartile1,y=quartile2,ymax=quartile3,color="Data Quartiles"),position=position_nudge(x=0.1))+geom_pointrange(aes(x=strain,ymin=kernal_quartile1,y=kernal_quartile2,ymax=kernal_quartile3,color="Kernal Quartiles"),position=position_nudge(x=-0.1))+geom_point(aes(x=strain,y=kernal_mode,color="Kernal Mode"),position=position_nudge(x=-0.1))+
  scale_x_discrete(limits=rev(summary_stats_strains))+
  scale_y_continuous(name="Cell Width",limits=c(0.5,1))+
  coord_flip()+
  as_md_theme(my_theme)

fivenum.width.p<-ggplot(fivenum.width)+geom_point(aes(x=strain,y=kernal_mode,shape="Mode",color=strain),size=1)+geom_pointrange(aes(x=strain,ymin=kernal_quartile1,y=kernal_quartile2,ymax=kernal_quartile3,shape="Kernal Quartiles",color=strain),fatten=2)+
  scale_x_discrete(name="",limits=rev(summary_stats_strains),labels=rev(labels_1))+
  scale_y_continuous(name="Cell Width (\u03bcm)",limits=c(0.65,0.85))+
  scale_color_manual(breaks=breaks_1,labels=labels_1,values=colors_1)+
  coord_flip()+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(fivenum.width.p)

ggsave(filename="figures/kernal_density_iqr.pdf",plot=fivenum.width.p,width =3,height=2,units="in",device=cairo_pdf)
