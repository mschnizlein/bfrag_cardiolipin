# Lipidomics Analysis
# UT-Knoxville Metabolite Core
# 5/1/24

#####

library(ggplot2)
library(tidyverse)
# library(reshape2)
# library(growthcurver)
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
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyr)
library(ggbiplot)
library(ggforce)
library(stringr)
library(lme4)
library(ggpubr)

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

setwd("C:/Users/mksch/Dropbox/matt_data/crosson_lab_msu_schnizlein/bfrag_bile_acid_stress/lipidomics/UTknoxville/")

#####
# Initial Data Processing

breaks_1<-c("WT","clsA","clsB","clsAclsB")
labels_1<-c("WT","Δ*clsA*","Δ*clsB*","Δ*clsA*Δ*clsB*")
colors_1<-c("#000000","#006198","#8BD1CB","#e38200")

km_cluster_colors<-c("km1"="#000000","km2"="#505050","km3"="#7a7a7a")

breaks_strains<-c("WT","WT_1","WT_2","clsA","clsA_1","clsA_2","clsB","clsB_1","clsB_2","clsA_clsB","clsA_clsB_1","clsA_clsB_2")
labels_strains<-c("WT 1","WT 2","WT 3","Δ*clsA* 1","Δ*clsA* 2","Δ*clsA* 3","Δ*clsB* 1","Δ*clsB* 2","Δ*clsB* 3","Δ*clsA*Δ*clsB* 1","Δ*clsA*Δ*clsB* 2","Δ*clsA*Δ*clsB* 3")

lipid_ontology_factors_full<-c(## Glycerophosphoglycerols
                                "BMP", # bismonoacylglycerophosphate
                                ## Glycosyldiradyl- and Diradyl-glycerols
                                "DG", # diacylglycerol
                                "MGDG", # monogalactosyldiacylglycerol
                                ## Acylaminosugars
                                "SL", # saccharolipids
                                ## Fatty amides
                                "NAGly", # N-acyl glycine
                                "NAGlySer", # N-acyl glycyl serine
                                "NAOrn", # N-acyl ornithine
                                "NAE", # N-acyl ethanolamines,
                                ## Fatty ester
                                "CAR", # acylcarnitine
                                ## Sphingolipids
                                "Cer", # ceramide
                                ### Acidic glycosphingolipids
                                "HexCer", # hex ceramide
                                ### Sphingoid bases
                                "SPB", # sphingoid bases
                                ### Phosphosphingolipids
                                "SM", # sphingomyelin
                                ## Phospholipids
                                "FA", # fatty acid
                                "PA", # phosphatidic acid
                                "LPA", # L-PA
                                "PC", # p-choline
                                "LPC", # lyso-PC
                                "PE", # p-ethanolamine
                                "LPE", # lyso-PE
                                "PS", # p-serine
                                "LPS", # lyso-PS
                                "PG", # p-glycerol
                                "LPG", # L-PG
                                "CL", # cardiolipin
                                "MLCL", # ML-CL
                                "DLCL") # DL-CL
# "AA" # amino acid and "Nucl" # nucleotide excluded from POS mode dataset
lipid_ontology_factors_pos<-c("BMP","CAR","DG","MGDG","SL","NAGly","NAGlySer","NAOrn","NAE","Cer","HexCer","SPB","SM","PC","LPC","PE","LPE","PS","LPS","CL")

lipid_ontology_factors<-c("FA","LPC","PA","LPA","PE","LPE","PG","LPG","PS","LPS","CL","MLCL","DLCL")

lipid_ontology_colors_full<-c("#ff6db6",
                         "#004949","#009292",
                         "#ffb6db",
                         "#006ddb","#b66dff","#6db6ff","#b6dbff",
                         "#920000",
                         "#924900","#db6d00","#24ff24","#ffff6d",
                         "#000000","darkgray","lightgray",SteppedSequential5Steps[c(1,4,6,9,11,14,16,19,21,23,25)])

lipid_ontology_colors_pos<-c("#000000","darkgray",
                             "#004949","#009292",
                             "#b66dff",
                             "#0F6B99","#2C85B2","#7EC3E5","#B2E5FF",
                             "#ff6db6","#CC5151","#E57E7E","#FFB2B2",
                             SteppedSequential5Steps[c(1,4,6,9,11,14,21)])

lipid_ontology_colors<-c("#000000","darkgray",SteppedSequential5Steps[c(1,4,6,9,11,14,16,19,21,23,25)])

lipid_data_pos<-as.data.table(read_xlsx(path="data/20092024_KAJ_Crosson_pos_lipids.xlsx",sheet="normalized")) # manually added some of the Ontology calls after receipt of data per lipid name
lipid_data_neg<-as.data.table(read_xlsx(path="data/29042024_KAJ_Crosson_initial_neg_lipids.xlsx",sheet="normalized"))
heatmap_labelling.dt<-read.csv(file="best_kmeans_cluster_log2_best.csv",header = TRUE) # heatmap labelling .csv with Kmeans clusters was generated below in the heatmap section. Moved here to add to the rest of the analysis

colnames(lipid_data_pos)<-c(as.character(lipid_data_pos[3,1:4]),colnames(lipid_data_pos)[5:ncol(lipid_data_pos)]) %>% gsub("[0-9+.]","",.) %>% make.unique(sep="_") %>% gsub("Δ","",.) %>% gsub(" ","_",.)
lipid_data_pos<-lipid_data_pos[4:nrow(lipid_data_pos),]
lipid_data_pos[,5:ncol(lipid_data_pos)]<-lapply(lipid_data_pos[,5:ncol(lipid_data_pos)],as.numeric)
lipid_data_pos<-lipid_data_pos[,-"Formula"]
lipid_data_pos$detect_mode <- "POS"

colnames(lipid_data_neg)<-c(as.character(lipid_data_neg[3,1:3]),colnames(lipid_data_neg)[4:ncol(lipid_data_neg)]) %>% gsub("[0-9+.]","",.) %>% make.unique(sep="_") %>% gsub("Δ","",.) %>% gsub(" ","_",.)
lipid_data_neg<-lipid_data_neg[4:nrow(lipid_data_neg),]
lipid_data_neg[,4:ncol(lipid_data_neg)]<-lapply(lipid_data_neg[,4:ncol(lipid_data_neg)],as.numeric)
lipid_data_neg$detect_mode <- "NEG"

lipid_data<-rbind(lipid_data_pos,lipid_data_neg)
lipid_data[rowMeans2(as.matrix(lipid_data[,c("WT","WT_1","WT_2")])) %in% 0,] # it's a very low abundant MLCL that appears in the clsB and clsAclsB strains but is absent in the WT strains, going to ignore it for the sake of mathematical simplicity
lipid_data<-lipid_data[lipid_data$`Metabolite_name_from_MS-DIAL` %!in% "w/o MS2: MLCL 22:6_12:0_19:2",]
lipid_data<-lipid_data[lipid_data$Ontology %!in% c("AA","Nucl"),]

colnames(lipid_data)<-gsub("-","",colnames(lipid_data))
sampleIDs<-colnames(lipid_data)[4:15]

# Summing for all peaks that have the same metabolite name from MSDIAL (these came out of the HPLC with unique peaks, but had similar fragmentation patterns)
lipid_data<-reframe(across(.cols = all_of(breaks_strains),.fns=sum), .by=c("Metabolite_name_from_MSDIAL","Manually_confirmed_lipid","Ontology","detect_mode"), .data=lipid_data)

# Adding Kmeans clusters from NEG mode
lipid_data<-merge(lipid_data,heatmap_labelling.dt[,c("Metabolite_name_from_MSDIAL","km_cluster")],by.all="Metabolite_name_from_MSDIAL",all.x = TRUE)

# Elaborating on additional metadata
lipid_data$lipid_size<-gsub("w/o MS2: ","",lipid_data$Manually_confirmed_lipid) %>% gsub("[A-Z]+","",ignore.case = T,.) %>% gsub(" ","",.) %>% gsub("\\|.*","",.) %>% gsub("-","",.)
lipid_data$Ontology<-factor(lipid_data$Ontology, levels=lipid_ontology_factors_full)

# Setting up final data.frame for analysis
lipid_data.f<-cbind(lipid_data,as.data.table(matrix(data="NA",nrow=nrow(lipid_data),ncol=3,dimnames=list(NULL,c("chain_length","chain_unsaturations","chain_oxidation")))))

# Generate length, unsaturation and oxidation from lipid size i.e., 18:0:2 or 18:0;3
for (i in 1:nrow(lipid_data.f)){
  lipid_size_temp<-strsplit(lipid_data.f$lipid_size,split="\\:|\\;") %>% lapply(.,as.data.frame)
  lipid_data.f[i,"chain_length"]<-lipid_size_temp[[i]][1,] %>% as.numeric(.)
  lipid_data.f[i,"chain_unsaturations"]<-lipid_size_temp[[i]][2,] %>% as.numeric(.)
  lipid_data.f[i,"chain_oxidation"]<-lipid_size_temp[[i]][3,] %>% as.numeric(.)
}

## Evaluating WT mean area compared to std deviation
lipid_data.f$WT_area_mean<-rowMeans2(as.matrix(lipid_data.f[,c("WT","WT_1","WT_2")]))
# lipid_data.f$WT_area_sd<-rowSds(as.matrix(lipid_data.f[,c("WT","WT_1","WT_2")]))
summary(lipid_data.f$WT_area_mean)
summary(rowSds(as.matrix(lipid_data.f[,c("WT","WT_1","WT_2")])))

ggplot(lipid_data.f,aes(x=log10(WT_area_mean),y=log10(rowSds(as.matrix(lipid_data.f[,c("WT","WT_1","WT_2")])))))+geom_point()
# Mean looks correlated with standard deviation of WT area under the peaks. This suggests that even though some have much higher standard deviation than the rest of the lipid contents, this is due to high area overall compared to those lower areas

lipid_data.f.m<-melt.data.table(as.data.table(lipid_data.f),id.vars = c("Metabolite_name_from_MSDIAL","Manually_confirmed_lipid","Ontology","km_cluster","lipid_size","chain_length","chain_unsaturations","chain_oxidation","detect_mode"))
lipid_data.f.m[,c("chain_length","chain_unsaturations","chain_oxidation")]<-lapply(lipid_data.f.m[,c("chain_length","chain_unsaturations","chain_oxidation")],function(x) as.numeric(x))

lipid_data.f.m$strain<-gsub("[0-9+]","",lipid_data.f.m$variable) %>%
  gsub("_","",.) %>%
  gsub("norm","",.) %>%
  factor(.,levels=c("WT","clsA","clsB","clsAclsB"))
lipid_data.f.m$value<-as.numeric(lipid_data.f.m$value)
lipid_data.f.m$value_logical<-as.logical(lipid_data.f.m$value)

## Making fold-change relative to wt as well as relative abundance values to sample and ontology
# fold change
lipid_data.f.m$value_fc<-0
for (i in 1:length(breaks_strains)){
  for (j in c("POS","NEG")){
    lipid_data.f.m[lipid_data.f.m$variable %in% breaks_strains[i] & lipid_data.f.m$detect_mode %in% j,]$value_fc<- lipid_data.f.m[lipid_data.f.m$variable %in% breaks_strains[i] & lipid_data.f.m$detect_mode %in% j,]$value/lipid_data.f.m[lipid_data.f.m$variable %in% "WT_area_mean" & lipid_data.f.m$detect_mode %in% j,]$value
  }
}

summary(lipid_data.f.m$value_fc)

# to sample
lipid_data.f.m$value_relabund<-0
for (i in 1:length(breaks_strains)){
  for (j in c("POS","NEG")){
    lipid_data.f.m[lipid_data.f.m$variable %in% breaks_strains[i] & lipid_data.f.m$detect_mode %in% j,]$value_relabund<- (lipid_data.f.m[lipid_data.f.m$variable %in% breaks_strains[i] & lipid_data.f.m$detect_mode %in% j,]$value/sum(lipid_data.f.m[lipid_data.f.m$variable %in% breaks_strains[i] & lipid_data.f.m$detect_mode %in% j,]$value))*100
  }
}

# to ontology
lipid_data.f.m$value_relabund_ont<-0
for (i in 1:length(breaks_strains)){
  for (j in c("POS","NEG")){
    for (k in lipid_ontology_factors_full){
      lipid_data.f.m[lipid_data.f.m$variable %in% breaks_strains[i] & lipid_data.f.m$detect_mode %in% j & lipid_data.f.m$Ontology %in% k,]$value_relabund_ont<- (lipid_data.f.m[lipid_data.f.m$variable %in% breaks_strains[i] & lipid_data.f.m$detect_mode %in% j & lipid_data.f.m$Ontology %in% k,]$value/sum(lipid_data.f.m[lipid_data.f.m$variable %in% breaks_strains[i] & lipid_data.f.m$detect_mode %in% j & lipid_data.f.m$Ontology %in% k,]$value))*100
    }
  }
}

ggplot(lipid_data.f.m,aes(x=Manually_confirmed_lipid,y=log2(value)))+
  geom_point()+
  geom_smooth(method="lm")+facet_wrap(~variable)

ggplot(lipid_data.f.m,aes(x=Manually_confirmed_lipid,y=log2(value_fc)))+
  geom_point()+
  geom_smooth(method="lm")+facet_wrap(~variable)

# Making data frames with single entries for ontology or per strain (or per strain/ontology)
lipid_data.f.m.ontology<-reframe(across(.cols=c("value","value_relabund","value_relabund_ont"),.fns=sum),.by=c("variable","Ontology","strain","detect_mode"),.data=lipid_data.f.m)

# Fold change of ontologies relative to wt mean
lipid_data.f.m.ontology$value_fc<-0
for (i in 1:length(breaks_strains)){
  for (j in c("POS","NEG")){
    lipid_data.f.m.ontology[lipid_data.f.m.ontology$variable %in% breaks_strains[i] & lipid_data.f.m.ontology$detect_mode %in% j,]$value_fc<- lipid_data.f.m.ontology[lipid_data.f.m.ontology$variable %in% breaks_strains[i] & lipid_data.f.m.ontology$detect_mode %in% j,]$value/lipid_data.f.m.ontology[lipid_data.f.m.ontology$variable %in% "WT_area_mean" & lipid_data.f.m.ontology$detect_mode %in% j,]$value
  }
}
summary(lipid_data.f.m.ontology$value_fc)

lipid_data.f.m<-lipid_data.f.m[lipid_data.f.m$variable %!in% "WT_area_mean",]
lipid_data.f.m.ontology<-lipid_data.f.m.ontology[lipid_data.f.m.ontology$variable %!in% "WT_area_mean",]

lipid_data.f.m.avg<-reframe(across(.cols=c("value","value_relabund","value_relabund_ont","value_fc"),.fns=mean), .by=c("Metabolite_name_from_MSDIAL","Manually_confirmed_lipid","Ontology","km_cluster","lipid_size","chain_length","chain_unsaturations","chain_oxidation","detect_mode","strain"), .data=lipid_data.f.m)
lipid_data.f.m.avg.ontology<-reframe(across(.cols=c("value","value_relabund","value_relabund_ont"),.fns=mean),.by=c("Ontology","strain","detect_mode"),.data=lipid_data.f.m.ontology)

sum(lipid_data.f.m.ontology$value_relabund)
sum(lipid_data.f.m.avg.ontology$value_relabund)

lipid_ontology_factors_1chain<-c("FA","LPC","LPE","LPS","LPA")
lipid_ontology_factors_2chain<-c("PA","PE","PS","PG")
lipid_ontology_factors_dlcl<-c("DLCL")
lipid_ontology_factors_3.4chain<-c("CL","MLCL")

lipid_data.f.m.avg$lipid_size2<-as.character(lipid_data.f.m.avg$lipid_size)
lipid_data.f.m.avg[lipid_data.f.m.avg$value_relabund_ont < 1 & !is.na(lipid_data.f.m.avg$value_relabund_ont),]$lipid_size2<-"Other"
lipid_data.f.m.avg$lipid_size3<-gsub("\\|.*","",lipid_data.f.m.avg$lipid_size2)

lipid_data.f.m.avg.cons<-reframe(across(.cols = "value_relabund_ont",.fns=sum), .by=c("lipid_size3","Ontology","strain","chain_length","chain_unsaturations","km_cluster","detect_mode"), .data=lipid_data.f.m.avg)

#####
# POS mode relabund analysis
## Rel abund in pos mode color by ontology present
## Rel abund in each class per strain
## Peak area corr between ontology and strain

lipid_relabund_posneg.p<-ggplot(lipid_data.f.m,aes(x = variable,y=value_relabund,fill=Ontology))+geom_bar(stat="identity")+
  scale_x_discrete(limits=breaks_strains,breaks=breaks_strains,labels=labels_strains,name="Strains")+
  scale_fill_manual(breaks=lipid_ontology_factors_full,values=lipid_ontology_colors_full)+
  ylab("Lipid Relative Abundance (%)")+
  as_md_theme(my_theme2)+facet_wrap(~detect_mode)
plot(lipid_relabund_posneg.p)

ggsave(filename = "figures/lipidomics_relabund_posneg.pdf",plot=lipid_relabund_posneg.p,width=6,height=4,units="in",device=cairo_pdf)

lipid_relabund_pos.p<-ggplot(lipid_data.f.m[lipid_data.f.m$detect_mode %in% "POS",],aes(x = variable,y=value_relabund,fill=Ontology))+geom_bar(stat="identity")+
  scale_x_discrete(limits=breaks_strains,breaks=breaks_strains,labels=labels_strains,name="Strains")+
  scale_fill_manual(breaks=lipid_ontology_factors_full,values=lipid_ontology_colors_full)+
  ylab("Lipid Relative Abundance (%)")+
  as_md_theme(my_theme2)
plot(lipid_relabund_pos.p)
ggsave(filename = "figures/lipidomics_relabund_pos.pdf",plot=lipid_relabund_pos.p,width=4,height=8,units="in",device=cairo_pdf)


lipid_relabund_neg.p<-ggplot(lipid_data.f.m[lipid_data.f.m$detect_mode %in% "NEG",],aes(x = variable,y=value_relabund,fill=Ontology))+geom_bar(stat="identity")+
  scale_x_discrete(limits=breaks_strains,breaks=breaks_strains,labels=labels_strains,name="Strains")+
  scale_fill_manual(breaks=lipid_ontology_factors_full,values=lipid_ontology_colors_full)+
  ylab("Lipid Relative Abundance (%)")+
  as_md_theme(my_theme2)
plot(lipid_relabund_neg.p)
ggsave(filename = "figures/lipidomics_relabund_neg.pdf",plot=lipid_relabund_neg.p,width=4,height=3,units="in",device=cairo_pdf)


##
ontology_peakarea_pos_high.p<-ggplot(lipid_data.f.m.ontology[lipid_data.f.m.ontology$detect_mode %in% "POS" & lipid_data.f.m.ontology$Ontology %in% c("NAGly","NAE","NAGlySer"),],aes(x = strain,y=value))+geom_bar(stat="summary",fun="mean",position=position_dodge(width=0.9))+stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.5, position=position_dodge(width=0.9)) +
  scale_x_discrete(breaks=breaks_1,labels=labels_1,name="")+
  ylab("Lipid Peak Area")+
  as_md_theme(my_theme)+theme(legend.position = "none")+facet_wrap(~Ontology)
plot(ontology_peakarea_pos_high.p)

ontology_peakarea_pos_low.p<-ggplot(lipid_data.f.m.ontology[lipid_data.f.m.ontology$detect_mode %in% "POS" & lipid_data.f.m.ontology$Ontology %!in% c("PE","LPE","PS","LPS","PC","LPC","CL","NAGly","NAE","NAGlySer"),],aes(x = strain,y=value))+
  geom_bar(stat="summary",fun="mean",position=position_dodge(width=0.9))+stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.5, position=position_dodge(width=0.9)) +
  scale_x_discrete(breaks=breaks_1,labels=labels_1,name="")+
  ylab("Lipid Peak Area")+
  as_md_theme(my_theme)+theme(legend.position = "none")+facet_wrap(~Ontology)
plot(ontology_peakarea_pos_low.p)

ontology_fc_pos.p<-ggplot(lipid_data.f.m.ontology[lipid_data.f.m.ontology$detect_mode %in% "POS" & lipid_data.f.m.ontology$Ontology %!in% c("PE","LPE","PS","LPS","PC","LPC","CL"),],aes(x = strain,y=log2(value_fc)))+
  geom_bar(stat="summary",fun="mean",position=position_dodge(width=0.9),fill="black")+stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.5, position=position_dodge(width=0.9)) +
  scale_x_discrete(breaks=breaks_1,labels=labels_1,name="")+
  ylab("Lipid log2(Fold-Change)")+
  as_md_theme(my_theme)+theme(legend.position = "none")+facet_wrap(~Ontology,nrow=2)
plot(ontology_fc_pos.p)
ggsave(filename = "figures/lipidomics_posmode_fc_ontology.pdf",plot=ontology_fc_pos.p,width=7,height=4,units="in",device=cairo_pdf)

dg_strain.lm<-glm(value~strain,data=lipid_data.f.m.ontology[lipid_data.f.m.ontology$Ontology %in% "DG",])
summary(dg_strain.lm) # cB and cAB 0.025 and 0.00873
mgdg_strain.lm<-glm(value~strain,data=lipid_data.f.m.ontology[lipid_data.f.m.ontology$Ontology %in% "MGDG",])
summary(mgdg_strain.lm) # cA cB and cAB 0.000217, 0.049 and 4.16E5
sl_strain.lm<-glm(value~strain,data=lipid_data.f.m.ontology[lipid_data.f.m.ontology$Ontology %in% "SL",])
summary(sl_strain.lm) # cB and cAB 0.0201 and 0.0177
naglyser_strain.lm<-glm(value~strain,data=lipid_data.f.m.ontology[lipid_data.f.m.ontology$Ontology %in% "NAGlySer",])
summary(naglyser_strain.lm) # cB and cAB 0.00017 and 0.000687
naorn_strain.lm<-glm(value~strain,data=lipid_data.f.m.ontology[lipid_data.f.m.ontology$Ontology %in% "NAOrn",])
summary(naorn_strain.lm) # cB and cAB 0.0341 and 0.00583
nae_strain.lm<-glm(value~strain,data=lipid_data.f.m.ontology[lipid_data.f.m.ontology$Ontology %in% "NAE",])
summary(nae_strain.lm) # cA cB and cAB 0.048 0.0335 and 0.0018
car_strain.lm<-glm(value~strain,data=lipid_data.f.m.ontology[lipid_data.f.m.ontology$Ontology %in% "CAR",])
summary(car_strain.lm) # none
cer_strain.lm<-glm(value~strain,data=lipid_data.f.m.ontology[lipid_data.f.m.ontology$Ontology %in% "Cer",])
summary(cer_strain.lm) # cB 0.0141
hexcer_strain.lm<-glm(value~strain,data=lipid_data.f.m.ontology[lipid_data.f.m.ontology$Ontology %in% "HexCer",])
summary(hexcer_strain.lm) # cA and cAB 0.000727 0.000143
spb_strain.lm<-glm(value~strain,data=lipid_data.f.m.ontology[lipid_data.f.m.ontology$Ontology %in% "SPB",])
summary(spb_strain.lm) # cAB 0.0325
sm_strain.lm<-glm(value~strain,data=lipid_data.f.m.ontology[lipid_data.f.m.ontology$Ontology %in% "SM",])
summary(sm_strain.lm) # cB 0.000265



cl_relabund_pos.p<-ggplot(lipid_data.f.m.ontology[lipid_data.f.m.ontology$Ontology %in% "CL" & lipid_data.f.m.ontology$detect_mode %in% "POS",],aes(x = strain,y=value_relabund))+geom_bar(stat="summary",fun="mean",position=position_dodge(width=0.9))+stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.5, position=position_dodge(width=0.9)) +
  scale_x_discrete(breaks=breaks_1,labels=labels_1,name="")+
  ylab("Cardiolipin Peak Area POS")+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(cl_relabund_pos.p)
sum(lipid_data.f.m.avg[lipid_data.f.m.avg$detect_mode %in% "POS" & lipid_data.f.m.avg$strain %in% "WT" & lipid_data.f.m.avg$Ontology %in% "CL",]$value_relabund)
sum(lipid_data.f.m.avg[lipid_data.f.m.avg$detect_mode %in% "NEG" & lipid_data.f.m.avg$strain %in% "WT" & lipid_data.f.m.avg$Ontology %in% "CL",]$value_relabund)

#####
# Negative mode relabund analysis

# data table for lefse analysis
# lipid_data_relabund.cons<-aggregate(.~Manually_confirmed_lipid+Ontology, data=lipid_data_relabund[,2:15],FUN=sum)
# rownames(lipid_data_relabund.cons)<-lipid_data_relabund.cons$Manually_confirmed_lipid
## This didn't turn anything up

## How to analyze lipidomics data set
# Ratio of lipids present across mutants
# - how do they change?
# - what is their composition? i.e. ontology (DONE)
# - Chain length of membrane lipids (particularly classes such as cardiolipin, PE, PG, PE) - histogram by ontology (DONE)
# - Heatmap by lipid ontology class - average fold change relative to WT cells (DONE)
# - Does distribution of individual lipid classes reflect cardiolipin distribution?
#     + create sets of chain lengths to see if that is tied to cardiolipin abundances (or are different populations chosen selectively)
# Use K-means clusters to characterize contributions of each cardiolipin synthase

# Presence/absence in each of the mutants (convert to logical)
# - which absences overlap, which do not
# - what types of species are present in each of those backgrounds
# - which lipid classes were undetectable in WT conditions? Did any of those become detectable

# Biplot of lipids
# - with cardiolipin species
# - without cardiolipin species


# To research
# - mono and dilyso cardiolipin

# Total peak height per sample
# lipid_data.f.m.summed<-aggregate(value~chain_length+strain+variable+norm_status,data=lipid_data.f.m[lipid_data.f.m$Ontology %in% c("PA","PG","PS","PE") & lipid_data.f.m$strain %in% c("WT","clsAclsB"),],FUN=sum)

ggplot(lipid_data.f.m[lipid_data.f.m$Ontology %in% c("PE","PS"),], aes(x=variable,y=value,fill=as.factor(chain_length),color=as.factor(chain_length)))+geom_bar(stat="identity")+facet_wrap(~Ontology)+scale_fill_discrete()

ggplot(lipid_data.f.m[lipid_data.f.m$Ontology %!in% c("PE","LPC","LPG","PS","CL","MLCL","DLCL"),], aes(x=variable,y=value,fill=as.factor(chain_length),color=as.factor(chain_length)))+geom_bar(stat="identity")+facet_wrap(~Ontology)+scale_fill_discrete()

ggplot(lipid_data.f.m[lipid_data.f.m$Ontology %in% "PG",], aes(x=variable,y=value,fill=lipid_size,color=lipid_size))+geom_bar(stat="identity")+facet_wrap(~Ontology)


#####
## Overview of lipid classes
lipid_relabund.p<-ggplot(lipid_data.f.m[lipid_data.f.m$detect_mode %in% "NEG",],aes(x = variable,y=value_relabund,fill=Ontology))+geom_bar(stat="identity")+
  scale_x_discrete(limits=breaks_strains,breaks=breaks_strains,labels=labels_strains,name="Strains")+
  scale_fill_manual(breaks=lipid_ontology_factors,values=lipid_ontology_colors)+
  ylab("Lipid Relative Abundance (%)")+
  as_md_theme(my_theme2)
plot(lipid_relabund.p)

ggsave(filename = "figures/lipidomics_relabund.pdf",plot=lipid_relabund.p,width=4,height=3,units="in",device=cairo_pdf)

# average ontology in each strain background
lipid_data.f.m.avg.ontology[lipid_data.f.m.avg.ontology$strain %in% "WT",]
lipid_data.f.m.avg.ontology[lipid_data.f.m.avg.ontology$strain %in% "clsA",]
lipid_data.f.m.avg.ontology[lipid_data.f.m.avg.ontology$strain %in% "clsB",]
lipid_data.f.m.avg.ontology[lipid_data.f.m.avg.ontology$strain %in% "clsAclsB",]

cl_relabund.p<-ggplot(lipid_data.f.m.ontology[lipid_data.f.m.ontology$Ontology %in% c("CL","MLCL","DLCL") & lipid_data.f.m.ontology$detect_mode %in% "NEG",],aes(x = Ontology,y=value_relabund,fill=strain))+geom_boxplot(position="dodge")+
  scale_x_discrete()+
  scale_fill_manual(breaks=breaks_1,values=colors_1,labels=labels_1)+
  ylab("Lipid Relative Abundance (%)")+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(cl_relabund.p)

ggsave(filename = "figures/lipidomics_relabund_cl.pdf",plot=cl_relabund.p,width=3,height=2.7,units="in",device=cairo_pdf)

# lipid_data_test<-dcast(data=lipid_data.f.m.ontology,formula = strain+variable~Ontology)

cl_strain.lm<-lm(value_relabund~strain,data=lipid_data.f.m.ontology[lipid_data.f.m.ontology$Ontology %in% "CL"& lipid_data.f.m.ontology$detect_mode %in% "NEG",])
mlcl_strain.lm<-lm(value_relabund~strain,data=lipid_data.f.m.ontology[lipid_data.f.m.ontology$Ontology %in% "MLCL"& lipid_data.f.m.ontology$detect_mode %in% "NEG",])
dlcl_strain.lm<-lm(value_relabund~strain,data=lipid_data.f.m.ontology[lipid_data.f.m.ontology$Ontology %in% "DLCL"& lipid_data.f.m.ontology$detect_mode %in% "NEG",])

summary(cl_strain.lm)
# all significant (clsA, clsB and clsAclsB) p = 3.3E-5, 2.5E-7 and 7.6E-10 respectively
summary(mlcl_strain.lm)
# clsB and clsAclsB are significant p = 2.7E-6 and 6.15E-7 respectively
summary(dlcl_strain.lm)
# none significant

#####
## Heatmap
lipid_data.f.m$value_fc_log2<-log2(lipid_data.f.m$value_fc)
lipid_data.f.m[lipid_data.f.m$value_fc_log2 %in% "-Inf",]$value_fc_log2<- -10
lipid_data.f.m$value_fc_curt<-(lipid_data.f.m$value_fc)^(1/3)
lipid_data_heatmap<-lipid_data.f.m[lipid_data.f.m$strain %!in% "WT" & lipid_data.f.m$detect_mode %in% "NEG",]
lipid_data_heatmap_log2<-dcast(lipid_data_heatmap,formula = Metabolite_name_from_MSDIAL+Manually_confirmed_lipid+Ontology+lipid_size+chain_length+chain_unsaturations+km_cluster+lipid_size ~ variable, value.var='value_fc_log2')
lipid_data_heatmap_curt<-dcast(lipid_data_heatmap,formula = Metabolite_name_from_MSDIAL+Manually_confirmed_lipid+Ontology+lipid_size+chain_length+chain_unsaturations+km_cluster+lipid_size ~ variable, value.var='value_fc_curt')

# Based on hist of data, cubed root results in better normalization
hist(melt.data.table(as.data.table(lipid_data_heatmap_log2[,9:ncol(lipid_data_heatmap_log2)]))$value)
summary(melt.data.table(as.data.table(lipid_data_heatmap_log2[,9:ncol(lipid_data_heatmap_log2)]))$value)
hist(melt.data.table(as.data.table(lipid_data_heatmap_curt[,9:ncol(lipid_data_heatmap_curt)]))$value)

lipid_data_heatmap_log2[,9:ncol(lipid_data_heatmap_log2)]<-lapply(lipid_data_heatmap_log2[,9:ncol(lipid_data_heatmap_log2)],function(x) gsub(-Inf, -10,x)) %>% lapply(.,as.numeric)

## Adding Kmeans clusters to datasets
set.seed(2)
# Used this code to generate K means (then wrote it to a csv)
# Determine Ideal number of K-Means clusters by Elbow Method, which plots the within-cluster sum of squares based on K cluster n
# k.max<-15
#  wss <- sapply(1:k.max,
#                function(k){kmeans(lipid_data_heatmap_log2[,9:ncol(lipid_data_heatmap_log2)], k, nstart=50,iter.max = 15 )$tot.withinss})
#    plot(1:k.max, wss,
#       type="b", pch = 19, frame = FALSE,
#       xlab="Number of clusters K",
#       ylab="Total within-clusters sum of squares (Cubed Root)")

#    pdf(file = "figures/lipidomics_kmeans_elbowmethod_log2fc.pdf",width = 8,height = 6) # units in in
#    plot(1:k.max, wss,
#         type="b", pch = 19, frame = FALSE,
#         xlab="Number of clusters K",
#         ylab="Total within-clusters sum of squares (log2 FC)")
#    dev.off()

# 4 appears to be the optimal number of clusters for both square and cubed roots, so I performed analysis with 4 clusters; Log2 fold change might have ideal at 3 but 4 breaks them more evenly across strains
# km.out<-kmeans(lipid_data_heatmap_log2[,c(9:17)],centers=3)
# lipid_data_heatmap_log2$km_cluster<-paste0("km",km.out$cluster)
# write.csv(lipid_data_heatmap_log2,file="best_kmeans_cluster_log2_best.csv",row.names = FALSE)

names(lipid_ontology_colors) <- lipid_ontology_factors

col <- list(Ontology = lipid_ontology_colors,
            km_cluster=km_cluster_colors)

row_ha<-rowAnnotation(Ontology=lipid_data_heatmap_log2$Ontology, km_cluster=lipid_data_heatmap_log2$km_cluster, col=col)

abund_cols.log2<-colorRamp2(c(3, 0, -3), brewer.pal(n=3, name="RdBu"))
abund_cols.cu<-colorRamp2(c(2, 1, 0), brewer.pal(n=3, name="RdBu"))
labels_strains_heatmap<-gsub("[[:punct:]]","",labels_strains)

lipidfc.log2.heatmap.p<-Heatmap(as.matrix(lipid_data_heatmap_log2[,9:ncol(lipid_data_heatmap_log2)]),name="Log2(Fold-Change)",left_annotation = row_ha,k=3,col=abund_cols.log2,column_labels = labels_strains_heatmap[4:12])
plot(lipidfc.log2.heatmap.p)

lipidfc.curt.heatmap.p<-Heatmap(as.matrix(lipid_data_heatmap_curt[,9:ncol(lipid_data_heatmap_curt)]),name="CuRt(Fold-Change) Strain/WT",left_annotation = row_ha,k=4,col=abund_cols.cu,column_labels = labels_strains_heatmap[4:12])
plot(lipidfc.curt.heatmap.p)

# K means is not consistent between generation methods, so may need to regenerate a plot until it matches the clustering from a previous attempt

pdf(file = "figures/lipidomics_lipid_log2_fc_heatmap.pdf",
    width = 5,height = 5) # size in inches
 plot(lipidfc.log2.heatmap.p)
 dev.off()

# pdf(file = "figures/lipidomics_lipid_curt_fc_heatmap.pdf",
#    width = 6,height = 5) # size in inches
# plot(lipidfc.curt.heatmap.p)
#  dev.off()

 # Revisit below with New dataframe format
#####
## Further exploration of lipids in data set by strain
unique(lipid_data.f.m.avg.cons[lipid_data.f.m.avg.cons$Ontology %in% lipid_ontology_factors_1chain & lipid_data.f.m.avg.cons$detect_mode %in% "NEG",]$lipid_size)
unique(lipid_data.f.m.avg.cons[lipid_data.f.m.avg.cons$Ontology %in% lipid_ontology_factors_2chain & lipid_data.f.m.avg.cons$detect_mode %in% "NEG",]$lipid_size3)
unique(lipid_data.f.m.avg.cons[lipid_data.f.m.avg.cons$Ontology %in% lipid_ontology_factors_3.4chain & lipid_data.f.m.avg.cons$detect_mode %in% "NEG",]$lipid_size3)

ontology_details.p<-ggplot(lipid_data.f.m.avg.cons[lipid_data.f.m.avg.cons$strain %in% "WT",],aes(x = Ontology,y=value_relabund_ont,fill=lipid_size3))+geom_bar(stat="identity")+
  ylab("Lipid Relative Abundance (%)")+
  as_md_theme(my_theme2)
plot(ontology_details.p)

# Summarizing lipids that have less than 1% abundance in particular ontology pool
length(unique(lipid_data.f.m.avg[lipid_data.f.m.avg$Ontology %in% lipid_ontology_factors_1chain & lipid_data.f.m.avg$lipid_size2 %in% "Other" & lipid_data.f.m.avg$detect_mode %in% "NEG",]$lipid_size)) # 10 lipids
length(unique(lipid_data.f.m.avg[lipid_data.f.m.avg$Ontology %in% lipid_ontology_factors_2chain & lipid_data.f.m.avg$lipid_size2 %in% "Other" & lipid_data.f.m.avg$detect_mode %in% "NEG",]$lipid_size)) # 28 lipids
length(unique(lipid_data.f.m.avg[lipid_data.f.m.avg$Ontology %in% lipid_ontology_factors_dlcl & lipid_data.f.m.avg$lipid_size2 %in% "Other" & lipid_data.f.m.avg$detect_mode %in% "NEG",]$lipid_size)) # 14 lipids
length(unique(lipid_data.f.m.avg[lipid_data.f.m.avg$Ontology %in% "MLCL" & lipid_data.f.m.avg$lipid_size2 %in% "Other" & lipid_data.f.m.avg$detect_mode %in% "NEG",]$lipid_size)) # 27 lipids
length(unique(lipid_data.f.m.avg[lipid_data.f.m.avg$Ontology %in% "CL" & lipid_data.f.m.avg$lipid_size2 %in% "Other" & lipid_data.f.m.avg$detect_mode %in% "NEG" & !is.na(lipid_data.f.m.avg$value_relabund_ont),]$lipid_size)) # 60 lipids


# ontology_1chain_breaks<-c("11:0","12:0","12:1","13:0","14:0","14:2","15:0","15:1","15:2","15:3","16:0","16:1","16:2","17:2","17:3","17:4","18:0","18:1","18:3","19:0","19:4","20:1","20:3","20:5")
ontology_1chain_breaks<-c("11:0","12:0","13:0","14:0","14:2","15:0","15:2","15:3","16:0","16:1","16:2","17:2","17:3","17:4","18:0","18:3","19:4","20:1","Other")
ontology_1chain_labels<-c("11:0","12:0","13:0","14:0","14:2","15:0","15:2","15:3","16:0","16:1","16:2","17:2","17:3","17:4","18:0","18:3","19:4","20:1","Other (10 lipids)")
ontology_1chain_colors<-c("#000000","#1e1d1d","#353334","#4e4b4e","#676469","#6B990F","#85B22C","#A3CC51","#0F6B99","#2C85B2","#51A3CC","#260F99","#422CB2","#6551CC","#990F0F","#B22C2C","#059494","#99540F","lightgrey")

ontology_details_1chain.p<-ggplot(lipid_data.f.m.avg.cons[lipid_data.f.m.avg.cons$Ontology %in% lipid_ontology_factors_1chain & lipid_data.f.m.avg.cons$detect_mode %in% "NEG",],aes(x = strain,y=value_relabund_ont,fill=lipid_size3))+geom_bar(stat="identity")+
  ylab("Lipid Relative Abundance (%)")+
  scale_fill_manual(values=ontology_1chain_colors,breaks=ontology_1chain_breaks,labels=ontology_1chain_labels)+
  as_md_theme(my_theme2)+theme(legend.position="none")+facet_wrap(~Ontology,nrow = 1)
plot(ontology_details_1chain.p)

ontology_2chain_colors<-c("#000000","#1e1d1d","#4e4b4e","#676469","#6B990F","#A3CC51","#C3E57E","#0F6B99","#260F99","#422CB2","#990F0F","#B22C2C","#CC5151","lightgrey")
ontology_2chain_breaks<-c("27:0","28:0","29:0","29:2","30:0","30:2","30:3","31:0","32:0","32:3","33:0","33:2","33:3","Other")
ontology_2chain_labels<-c("27:0","28:0","29:0","29:2","30:0","30:2","30:3","31:0","32:0","32:3","33:0","33:2","33:3","Other (28 lipids)")

ontology_details_2chain.p<-ggplot(lipid_data.f.m.avg.cons[lipid_data.f.m.avg.cons$Ontology %in% lipid_ontology_factors_2chain & lipid_data.f.m.avg.cons$detect_mode %in% "NEG",],aes(x = strain,y=value_relabund_ont,fill=lipid_size3))+geom_bar(stat="identity")+
  ylab("Lipid Relative Abundance (%)")+
  scale_fill_manual(values=ontology_2chain_colors,breaks=ontology_2chain_breaks,labels=ontology_2chain_labels)+
  as_md_theme(my_theme2)+facet_wrap(~Ontology,nrow=1)
plot(ontology_details_2chain.p)

ontology_dlcl_breaks<-c("25:0","27:3","28:0","29:1","30:4","31:2","31:6","32:4","32:6","32:5","33:3","33:5","35:8","36:8","37:8","37:9","38:7","38:8","38:10","39:8","40:8","40:10","Other")
ontology_dlcl_labels<-c("25:0","27:3","28:0","29:1","30:4","31:2","31:6","32:4","32:6","32:5","33:3","33:5","35:8","36:8","37:8","37:9","38:7","38:8","38:10","39:8","40:8","40:10","Other (14 lipids)")
ontology_dlcl_colors<-c("#000000","#1e1d1d","#353334","#4e4b4e","#676469","#6B990F","#85B22C","#0F6B99","#2C85B2","#51A3CC","#260F99","#422CB2","#6551CC","#8F7EE5","#990F0F","#B22C2C","#CC5151","#E57E7E","#FFB2B2","#059494","#99540F","#CC8E51","lightgrey")

ontology_details_dlcl.p<-ggplot(lipid_data.f.m.avg.cons[lipid_data.f.m.avg.cons$Ontology %in% lipid_ontology_factors_dlcl & lipid_data.f.m.avg.cons$detect_mode %in% "NEG",],aes(x = strain,y=value_relabund_ont,fill=lipid_size3))+geom_bar(stat="identity")+
  ylab("Lipid Relative Abundance (%)")+
  scale_fill_manual(values=ontology_dlcl_colors,breaks=ontology_dlcl_breaks,labels=ontology_dlcl_labels)+
  as_md_theme(my_theme2)+theme(legend.position="none")+facet_wrap(~Ontology,nrow=1)
plot(ontology_details_dlcl.p)

ontology_mlcl_breaks<-c("51:8","52:6","52:8","53:6","59:14","60:0","61:0","61:1","62:0","62:1","62:5","63:0","63:1","63:2","63:3","63:4","63:5","64:4","64:6","64:7","64:8","65:6","65:8","66:4","66:5","66:7","66:8","66:10","Other")
ontology_mlcl_labels<-c("51:8","52:6","52:8","53:6","59:14","60:0","61:0","61:1","62:0","62:1","62:5","63:0","63:1","63:2","63:3","63:4","63:5","64:4","64:6","64:7","64:8","65:6","65:8","66:4","66:5","66:7","66:8","66:10","Other (27 lipids)")
ontology_mlcl_colors<-c("#000000","#1e1d1d","#353334","#4e4b4e","#676469","#9a9aa2","#6B990F","#85B22C","#0F6B99","#2C85B2","#51A3CC","#260F99","#422CB2","#6551CC","#8F7EE5","#BFB2FF","#d8d1ff","#990F0F","#B22C2C","#CC5151","#E57E7E","#059494","#54aac8","#99540F","#B26F2C","#CC8E51","#E5B17E","#FFD8B2","lightgrey")

ontology_details_mlcl.p<-ggplot(lipid_data.f.m.avg.cons[lipid_data.f.m.avg.cons$Ontology %in% "MLCL" & lipid_data.f.m.avg.cons$detect_mode %in% "NEG",],aes(x = strain,y=value_relabund_ont,fill=lipid_size3))+geom_bar(stat="identity")+
  ylab("Lipid Relative Abundance (%)")+
  scale_fill_manual(values=ontology_mlcl_colors,breaks=ontology_mlcl_breaks,labels=ontology_mlcl_labels)+
  as_md_theme(my_theme2)+theme(legend.position="none")+facet_wrap(~Ontology,nrow=1)
plot(ontology_details_mlcl.p)

ontology_cl_breaks<-c("60:0","60:3","61:0","61:3","62:0","62:3","63:1","63:2","63:3","64:1","64:3","65:1","65:3","66:2","67:11","68:11","69:9","69:11","70:9","70:10","70:11","71:9","71:11","72:4","74:14","76:14","81:19","83:19","Other")
ontology_cl_labels<-c("60:0","60:3","61:0","61:3","62:0","62:3","63:1","63:2","63:3","64:1","64:3","65:1","65:3","66:2","67:11","68:11","69:9","69:11","70:9","70:10","70:11","71:9","71:11","72:4","74:14","76:14","81:19","83:19","Other (60 lipids)")
ontology_cl_colors<-c("#000000","#1e1d1d","#353334","#4e4b4e","#676469","#9a9aa2","#6B990F","#85B22C","#A3CC51","#C3E57E","#E5FFB2","#0F6B99","#51A3CC","#B2E5FF","#260F99","#422CB2","#6551CC","#8F7EE5","#990F0F","#B22C2C","#CC5151","#E57E7E","#FFB2B2","#20A2A2","#99540F","#B26F2C","#CC8E51","#E5B17E","lightgrey")


ontology_details_cl.p<-ggplot(lipid_data.f.m.avg.cons[lipid_data.f.m.avg.cons$Ontology %in% "CL" & lipid_data.f.m.avg.cons$detect_mode %in% "NEG" & !is.na(lipid_data.f.m.avg.cons$value_relabund_ont),],aes(x = strain,y=value_relabund_ont,fill=lipid_size3))+geom_bar(stat="identity")+
  scale_
  ylab("Lipid Relative Abundance (%)")+
  scale_fill_manual(values=ontology_cl_colors,breaks=ontology_cl_breaks,labels=ontology_cl_labels)+
  as_md_theme(my_theme2)+facet_wrap(~Ontology,nrow=1)
plot(ontology_details_cl.p)

ggsave(filename = "figures/lipidomics_relabund_normtoont_1chain.pdf",plot=ontology_details_1chain.p,width=4.5,height=1.75,units="in",device=cairo_pdf)
ggsave(filename = "figures/lipidomics_relabund_normtoont_2chain.pdf",plot=ontology_details_2chain.p,width=4.5,height=1.75,units="in",device=cairo_pdf)
ggsave(filename = "figures/lipidomics_relabund_normtoont_dlcl.pdf",plot=ontology_details_dlcl.p,width=1.5,height=1.75,units="in",device=cairo_pdf)
ggsave(filename = "figures/lipidomics_relabund_normtoont_mlcl.pdf",plot=ontology_details_mlcl.p,width=1.5,height=1.75,units="in",device=cairo_pdf)
ggsave(filename = "figures/lipidomics_relabund_normtoont_cl.pdf",plot=ontology_details_cl.p,width=1.5,height=1.75,units="in",device=cairo_pdf)

lipid_data.f.m.avg.cons[lipid_data.f.m.avg.cons$Ontology %in% "LPE" & lipid_data.f.m.avg.cons$lipid_size3 %in% "15:0",]$value_relabund_ont
lipid_data.f.m.avg.cons[lipid_data.f.m.avg.cons$Ontology %in% "LPS" & lipid_data.f.m.avg.cons$lipid_size3 %in% "15:0",]$value_relabund_ont
lipid_data.f.m.avg.cons[lipid_data.f.m.avg.cons$Ontology %in% "LPA" & lipid_data.f.m.avg.cons$lipid_size3 %in% "16:0",]$value_relabund_ont

lipid_data.f.m.avg.cons[lipid_data.f.m.avg.cons$Ontology %in% "PE" & lipid_data.f.m.avg.cons$lipid_size3 %in% "30:0",]$value_relabund_ont
lipid_data.f.m.avg.cons[lipid_data.f.m.avg.cons$Ontology %in% "PA" & lipid_data.f.m.avg.cons$lipid_size3 %in% "30:0",]$value_relabund_ont
lipid_data.f.m.avg.cons[lipid_data.f.m.avg.cons$Ontology %in% "PG" & lipid_data.f.m.avg.cons$lipid_size3 %in% "32:0",]$value_relabund_ont
lipid_data.f.m.avg.cons[lipid_data.f.m.avg.cons$Ontology %in% "PS" & lipid_data.f.m.avg.cons$lipid_size3 %in% "30:0",]$value_relabund_ont

######
# Unnormalized CL, MLCL and DLCL percentage
cl_species.p<-ggplot(lipid_data.f.m[lipid_data.f.m$Ontology %in% c("CL"),],aes(x = variable,y=value_relabund,fill=lipid_size))+geom_bar(stat="identity")+
  scale_x_discrete(limits=breaks_strains,breaks=breaks_strains,labels=labels_strains,name="Strains")+
  ylab("Lipid Relative Abundance (%)")+
  as_md_theme(my_theme2)
plot(cl_species.p)

mlcl_species.p<-ggplot(lipid_data.f.m[lipid_data.f.m$Ontology %in% c("MLCL"),],aes(x = variable,y=value_relabund,fill=lipid_size))+geom_bar(stat="identity")+
  scale_x_discrete(limits=breaks_strains,breaks=breaks_strains,labels=labels_strains,name="Strains")+
  ylab("Lipid Relative Abundance (%)")+
  as_md_theme(my_theme2)
plot(mlcl_species.p)

dlcl_species.p<-ggplot(lipid_data.f.m[lipid_data.f.m$Ontology %in% c("DLCL"),],aes(x = variable,y=value_relabund,fill=lipid_size))+geom_bar(stat="identity")+
  scale_x_discrete(limits=breaks_strains,breaks=breaks_strains,labels=labels_strains,name="Strains")+
  ylab("Lipid Relative Abundance (%)")+
  as_md_theme(my_theme2)
plot(dlcl_species.p)

ggsave(filename = "figures/lipidomics_relabund_cl.pdf",plot=cl_species.p,width=8,height=6,units="in",device=cairo_pdf)

## Chain Length and Unsaturation

cl_species_chainlength.p<-ggplot(lipid_data.f.m.avg[lipid_data.f.m.avg$Ontology %in% c("CL","MLCL","DLCL") & lipid_data.f.m.avg$detect_mode %in% "NEG",],aes(x=chain_length,y=value_relabund,fill=Ontology))+geom_bar(stat="identity",position="stack")+
  scale_x_continuous(limits=c(20,100),name="Total Carbons in Lipid Chains")+
  ylab("Lipid Relative Abundance (%)")+
  as_md_theme(my_theme2)+facet_wrap(~strain)
plot(cl_species_chainlength.p)

other_species_chainlength.p<-ggplot(lipid_data.f.m.avg[lipid_data.f.m.avg$Ontology %!in% c("CL","MLCL","DLCL") & lipid_data.f.m.avg$detect_mode %in% "NEG",],aes(x=chain_length,y=value_relabund,fill=Ontology))+geom_bar(stat="identity",position="stack")+
  scale_x_continuous(name="Total Carbons in Lipid Chains")+
  ylab("Lipid Relative Abundance (%)")+
  as_md_theme(my_theme2)+facet_wrap(~strain)
plot(other_species_chainlength.p)

cl_species_unsaturations.p<-ggplot(lipid_data.f.m.avg[lipid_data.f.m.avg$Ontology %in% c("CL","MLCL","DLCL") & lipid_data.f.m.avg$detect_mode %in% "NEG",],aes(x=chain_unsaturations,y=value_relabund,fill=Ontology))+geom_bar(stat="identity",position="stack")+
  scale_x_continuous(name="Total Unsaturations in Lipid Chains")+
  ylab("Lipid Relative Abundance (%)")+
  as_md_theme(my_theme2)+facet_wrap(~strain)
plot(cl_species_unsaturations.p)

other_species_unsaturations.p<-ggplot(lipid_data.f.m.avg[lipid_data.f.m.avg$Ontology %in% c("PA","PG","PE","PC","PS","LPA","LPG","LPE","LPC") & lipid_data.f.m.avg$detect_mode %in% "NEG",],aes(x=chain_unsaturations,y=value_relabund,fill=Ontology))+geom_bar(stat="identity",position="stack")+
  scale_x_continuous(name="Total Unsaturations in Lipid Chains")+
  ylab("Lipid Relative Abundance (%)")+
  as_md_theme(my_theme2)+facet_wrap(~strain)
plot(other_species_unsaturations.p)

other_species_unsaturations_no0.2.p<-ggplot(lipid_data.f.m.avg[lipid_data.f.m.avg$Ontology %in% c("PA","PG","PE","PC","PS","LPA","LPG","LPE","LPC") & lipid_data.f.m.avg$chain_unsaturations >=3 & lipid_data.f.m.avg$detect_mode %in% "NEG",],aes(x=chain_unsaturations,y=value,fill=Ontology))+geom_bar(stat="identity",position="stack")+
  scale_x_continuous(name="Total Unsaturations in Lipid Chains")+
  ylab("Lipid Relative Abundance (%)")+
  as_md_theme(my_theme2)+facet_wrap(~strain)
plot(other_species_unsaturations_no0.2.p)

ggsave(filename = "figures/lipidomics_chainlength_cl.pdf",plot=cl_species_chainlength.p,width=8,height=6,units="in",device=cairo_pdf)
ggsave(filename = "figures/lipidomics_chainlength_other.pdf",plot=other_species_chainlength.p,width=8,height=6,units="in",device=cairo_pdf)
ggsave(filename = "figures/lipidomics_unsaturations_cl.pdf",plot=cl_species_unsaturations.p,width=8,height=6,units="in",device=cairo_pdf)
ggsave(filename = "figures/lipidomics_unsaturations_other.pdf",plot=other_species_unsaturations.p,width=8,height=6,units="in",device=cairo_pdf)
ggsave(filename = "figures/lipidomics_unsaturations_other_no1to2.pdf",plot=other_species_unsaturations_no0.2.p,width=8,height=6,units="in",device=cairo_pdf)

unique(lipid_data.f.m[lipid_data.f.m$chain_unsaturations >=5 & lipid_data.f.m$Ontology %!in% c("CL","MLCL","DLCL"),]$Ontology)

unsaturation.lm<-glm(formula= value_relabund~strain+Ontology,lipid_data.f.m[lipid_data.f.m$chain_unsaturations >=5 & lipid_data.f.m$Ontology %!in% c("CL","MLCL","DLCL") & lipid_data.f.m.avg$detect_mode %in% "NEG",], family="gaussian")
summary(unsaturation.lm)
# clsA and clsAclsB both significantly lower high-unsaturated PG/PE in the membrane (p = 0.043 and p = 0.017 respectively)

# Chain Length analysis of PG, PA, PE, and PS
lipid_data.f.m.chainsum.twoacyl<-reframe(across(.cols = "value",.fns=sum), .by=c("chain_length","strain","variable"), .data=lipid_data.f.m[lipid_data.f.m$Ontology %in% c("PA","PG","PS","PE") & lipid_data.f.m$strain %in% c("WT","clsAclsB") & lipid_data.f.m$detect_mode %in% "NEG",])
lipid_data.f.m.chainsum.twoacyl.ontology<-reframe(across(.cols = "value",.fns=sum), .by=c("chain_length","strain","variable","Ontology"), .data=lipid_data.f.m[lipid_data.f.m$Ontology %in% c("PA","PG","PS","PE") & lipid_data.f.m$strain %in% c("WT","clsAclsB") & lipid_data.f.m$detect_mode %in% "NEG",])

chainsum_twoacyl_peakarea.p<-ggplot(lipid_data.f.m.chainsum.twoacyl[lipid_data.f.m.chainsum.twoacyl$strain %in% c("WT","clsAclsB"),], aes(x=chain_length,y=value,fill=strain,color=strain))+geom_point(position=position_jitterdodge(dodge.width=0,jitter.height = 0))+stat_summary(fun=mean, geom="line", linewidth=1)+
  scale_y_log10(name="log10(Peak Area)")+
  scale_x_continuous(name="Acyl Chain Length")+
  scale_color_manual(values=colors_1,breaks=breaks_1,labels=labels_1)+
  as_md_theme(my_theme)
plot(chainsum_twoacyl_peakarea.p)

chainsum_twoacyl_peakarea_ont.p<-ggplot(lipid_data.f.m.chainsum.twoacyl.ontology[lipid_data.f.m.chainsum.twoacyl.ontology$strain %in% c("WT","clsAclsB"),], aes(x=chain_length,y=value,fill=strain,color=strain))+geom_point(position=position_jitterdodge(dodge.width=0,jitter.height = 0))+stat_summary(fun=mean, geom="line", linewidth=1)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x,n=4), labels = trans_format("log10", label_math(10^.x)),limits=c(1E-5,1E1),name="log10(Peak Area)")+
  annotation_logticks(base=10,sides = "l",scaled=TRUE,short = unit(0.05, "cm"),mid = unit(0.1, "cm"),long = unit(0.2, "cm"))+scale_x_continuous(name="Acyl Chain Length")+
  scale_color_manual(values=colors_1,breaks=breaks_1,labels=labels_1)+
  as_md_theme(my_theme)+facet_wrap(~Ontology)
plot(chainsum_twoacyl_peakarea_ont.p)

ggsave(filename="figures/chainlength_twoacyl_peakarea.png",width=3,height=3, units="in",plot=chainsum_twoacyl_peakarea.p)
ggsave(filename="figures/chainlength_twoacyl_peakarea_ontology.png",width=5,height=4, units="in",plot=chainsum_twoacyl_peakarea_ont.p)

t.test(lipid_data.f.m.chainsum.twoacyl[lipid_data.f.m.chainsum.twoacyl$strain %in% "WT" & lipid_data.f.m.chainsum.twoacyl$chain_length %in% 28,]$value,
       lipid_data.f.m.chainsum.twoacyl[lipid_data.f.m.chainsum.twoacyl$strain %in% "clsAclsB" & lipid_data.f.m.chainsum.twoacyl$chain_length %in% 28,]$value)

t.test(lipid_data.f.m.chainsum.twoacyl[lipid_data.f.m.chainsum.twoacyl$strain %in% "WT" & lipid_data.f.m.chainsum.twoacyl$chain_length %in% 38,]$value,
       lipid_data.f.m.chainsum.twoacyl[lipid_data.f.m.chainsum.twoacyl$strain %in% "clsAclsB" & lipid_data.f.m.chainsum.twoacyl$chain_length %in% 38,]$value)

lipid_data.f.m.unsatsum.twoacyl<-reframe(across(.cols = "value",.fns=sum), .by=c("chain_unsaturations","strain","variable"), .data=lipid_data.f.m[lipid_data.f.m$Ontology %in% c("PA","PG","PS","PE") & lipid_data.f.m$strain %in% c("WT","clsAclsB") & lipid_data.f.m$detect_mode %in% "NEG",])

unsat_twoacyl_peakarea.p<-ggplot(lipid_data.f.m.unsatsum.twoacyl[lipid_data.f.m.unsatsum.twoacyl$strain %in% c("WT","clsAclsB"),], aes(x=chain_unsaturations,y=value,color=strain))+geom_point(position=position_jitterdodge(dodge.width=0,jitter.height = 0))+stat_summary(fun=mean, geom="line", linewidth=1)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x,n=4), labels = trans_format("log10", label_math(10^.x)),limits=c(1E-4,1E1),name="log10(Peak Area)")+
  annotation_logticks(base=10,sides = "l",scaled=TRUE,short = unit(0.05, "cm"),mid = unit(0.1, "cm"),long = unit(0.2, "cm"))+scale_x_continuous(name="Acyl Chain Unsaturation")+
  scale_color_manual(values=colors_1,breaks=breaks_1,labels=labels_1)+
  as_md_theme(my_theme)
plot(unsat_twoacyl_peakarea.p)

unsat_twoacyl_relabund.p<-ggplot(lipid_data.f.m.unsatsum.twoacyl[lipid_data.f.m.unsatsum.twoacyl$strain %in% c("WT","clsAclsB"),], aes(x=chain_unsaturations,y=value,color=strain))+geom_point(position=position_jitterdodge(dodge.width=0,jitter.height = 0))+stat_summary(fun=mean, geom="line", linewidth=1)+
  scale_y_continuous(name="Relative Abundance (%)")+scale_x_continuous(name="Acyl Chain Unsaturation")+
  scale_color_manual(values=colors_1,breaks=breaks_1,labels=labels_1)+
  as_md_theme(my_theme)
plot(unsat_twoacyl_relabund.p)

ggsave(filename="figures/unsat_twoacyl_peakarea.png",width=3,height=3, units="in",plot=unsat_twoacyl_peakarea.p)
ggsave(filename="figures/unsat_twoacyl_relabund.png",width=3,height=3, units="in",plot=unsat_twoacyl_relabund.p)


t.test(lipid_data.f.m.unsatsum.twoacyl[lipid_data.f.m.unsatsum.twoacyl$strain %in% "WT" & lipid_data.f.m.unsatsum.twoacyl$chain_unsaturations %in% 5,]$value,
       lipid_data.f.m.unsatsum.twoacyl[lipid_data.f.m.unsatsum.twoacyl$strain %in% "clsAclsB" & lipid_data.f.m.unsatsum.twoacyl$chain_unsaturations %in% 5,]$value) # p = 0.0455

t.test(lipid_data.f.m.unsatsum.twoacyl[lipid_data.f.m.unsatsum.twoacyl$strain %in% "WT" & lipid_data.f.m.unsatsum.twoacyl$chain_unsaturations %in% 6,]$value,
       lipid_data.f.m.unsatsum.twoacyl[lipid_data.f.m.unsatsum.twoacyl$strain %in% "clsAclsB" & lipid_data.f.m.unsatsum.twoacyl$chain_unsaturations %in% 6,]$value) # p = 0.0091

chain_length.lm<-glm(cbind(value+chain_length)~strain, data=lipid_data.f.m.avg[lipid_data.f.m.avg$Ontology %in% c("PA","PG","PS","PE") & lipid_data.f.m.avg$strain %in% c("WT","clsAclsB"),])
summary(chain_length.lm)

# Difference

chainsum.datasub<-lipid_data.f.m.chainsum.twoacyl[lipid_data.f.m.chainsum.twoacyl$strain %in% c("WT","clsAclsB"),]

chainsum.datasub.dcast<-dcast(as.data.table(chainsum.datasub),formula = chain_length~variable,value.var="value")
chainsum.datasub.dcast[,c("difference_cAB1","difference_cAB2","difference_cAB3")]<-chainsum.datasub.dcast[,c("clsA_clsB","clsA_clsB_1","clsA_clsB_2")]/rowMeans2(as.matrix(chainsum.datasub.dcast[,c("WT","WT_1","WT_2")]))
chainsum.datasub.dcast.m<-melt.data.table(chainsum.datasub.dcast,measure.vars = c("difference_cAB1","difference_cAB2","difference_cAB3"))

chainsum_twoacyl.fc.p<-ggplot(chainsum.datasub.dcast.m,aes(x=chain_length,y=log2(value)))+geom_bar(stat="summary", fill="black")+geom_errorbar(stat='summary', width=0.5, linewidth=0.5)+
  as_md_theme(my_theme)+
  scale_y_continuous(name="log2(FC) (Δ*clsA*Δ*clsB*/WT)")+
  scale_x_continuous(name="Acyl Chain Length")
plot(chainsum_twoacyl.fc.p)

ggsave(filename="figures/chainsum_twoacyl_fc.pdf",plot=chainsum_twoacyl.fc.p,width=3,height = 3,units="in",device=cairo_pdf)

unsat.datasub<-lipid_data.f.m.unsatsum.twoacyl[lipid_data.f.m.unsatsum.twoacyl$strain %in% c("WT","clsAclsB"),]
unsat.datasub.dcast<-dcast(as.data.table(unsat.datasub),formula = chain_unsaturations~variable,value.var="value")
unsat.datasub.dcast[,c("difference_cAB1","difference_cAB2","difference_cAB3")]<-unsat.datasub.dcast[,c("clsA_clsB","clsA_clsB_1","clsA_clsB_2"),]/rowMeans2(as.matrix(unsat.datasub.dcast[,c("WT","WT_1","WT_2")]))
unsat.datasub.dcast.m<-melt.data.table(unsat.datasub.dcast,measure.vars = c("difference_cAB1","difference_cAB2","difference_cAB3"))

unsat_twoacyl.fc.p<-ggplot(unsat.datasub.dcast.m,aes(x=chain_unsaturations,y=log2(value)))+geom_bar(stat="summary", fill="black")+geom_errorbar(stat='summary', width=0.5, linewidth=0.5)+
  as_md_theme(my_theme)+
  scale_y_continuous(name="log2(FC) (Δ*clsA*Δ*clsB*/WT)")+
  scale_x_continuous(name="Acyl Chain Unsaturations",limits=c(-0.5,6.5))
plot(unsat_twoacyl.fc.p)

ggsave(filename="figures/unsat_twoacyl_fc.pdf",plot=unsat_twoacyl.fc.p,width=3,height = 3,units="in",device=cairo_pdf)

#####
## K-means cluster analysis
# Stats associating clusters with strains?? (look at DMM clustering)
# What are lipids/classes in each K means cluster, is there an association with a particular headgroup or size/whatever
lipid_data.f.m.avg$value_km_adj<-NA

for (i in 1:3){
  for (j in 1:4){
    lipid_data.f.m.avg[lipid_data.f.m.avg$km_cluster %in% names(km_cluster_colors)[i] & lipid_data.f.m.avg$strain %in% breaks_1[j],]$value_km_adj<-
      (lipid_data.f.m.avg[lipid_data.f.m.avg$km_cluster %in% names(km_cluster_colors)[i] & lipid_data.f.m.avg$strain %in% breaks_1[j],]$value/
      sum(lipid_data.f.m.avg[lipid_data.f.m.avg$km_cluster %in% names(km_cluster_colors)[i] & lipid_data.f.m.avg$strain %in% breaks_1[j],]$value))*100
  }
}


kmeans_species.p<-ggplot(lipid_data.f.m.avg[lipid_data.f.m.avg$detect_mode %in% "NEG",],aes(x = km_cluster,y=value_km_adj,fill=Ontology))+geom_bar(stat="identity")+
  scale_fill_manual(values=lipid_ontology_colors,breaks=lipid_ontology_factors)+
  ylab("Lipid Relative Abundance (%) within K-means Clusters")+
  as_md_theme(my_theme)+facet_wrap(~strain)
plot(kmeans_species.p)

ggsave(filename = "figures/lipidomics_kmeans_lipidabund.pdf",plot=kmeans_species.p,width=8,height=6,units="in",device=cairo_pdf)


# ld_relabund.m.avg_complete<-complete(lipid_data_relabund.m.avg,chain_unsaturations,chain_length,strain,km_cluster,Ontology,fill=list(value=0),explicit = FALSE)
lipid_data.f.m$Ontology<-lipid_data.f.m$Ontology %>% as.character(.)

ld_relabund.m_complete<-lipid_data.f.m[lipid_data.f.m$detect_mode %in% "NEG",]

ld_relabund.m_complete.clmlcl.kmavg_unsat<-reframe(across(.cols = "value_relabund",.fns=sum), .by=c("km_cluster","chain_unsaturations","variable","strain"), .data=ld_relabund.m_complete[ld_relabund.m_complete$Ontology %in% c("CL","MLCL"),])
ld_relabund.m_complete.clmlcl.kmavg_unsat<-complete(ld_relabund.m_complete.clmlcl.kmavg_unsat,chain_unsaturations,variable,km_cluster,fill=list(value_relabund=0),explicit = FALSE)

ld_relabund.m_complete.clmlcl.kmavg_chain<-reframe(across(.cols = "value_relabund",.fns=sum), .by=c("km_cluster","chain_length",,"variable","strain"), .data=ld_relabund.m_complete[ld_relabund.m_complete$Ontology %in% c("CL","MLCL"),])
ld_relabund.m_complete.clmlcl.kmavg_chain<-complete(ld_relabund.m_complete.clmlcl.kmavg_chain,chain_length,variable,km_cluster,fill=list(value_relabund=0),explicit = FALSE)

ld_relabund.m_complete.clmlcl.kmavg_chain$strain<-gsub("[0-9]+","",ld_relabund.m_complete.clmlcl.kmavg_chain$variable) %>% gsub("_","",.) %>% factor(.,levels=breaks_1)
ld_relabund.m_complete.clmlcl.kmavg_unsat$strain<-gsub("[0-9]+","",ld_relabund.m_complete.clmlcl.kmavg_unsat$variable) %>% gsub("_","",.) %>% factor(.,levels=breaks_1)

# ld_relabund.m_complete.clmlcl.kmavg_unsat<-complete(ld_relabund.m_complete.clmlcl.kmavg_unsat,km_cluster,chain_unsaturations,variable,strain,fill=list(value=0),explicit=FALSE)

cl_species_chainlength_wt_kmeans.1.2.p<-ggplot(ld_relabund.m_complete.clmlcl.kmavg_chain[ld_relabund.m_complete.clmlcl.kmavg_chain$variable %in% c("WT","WT_1","WT_2") & ld_relabund.m_complete.clmlcl.kmavg_chain$km_cluster %in% c("km1","km2"),],aes(x=chain_length,y=value_relabund,color=km_cluster,fill=km_cluster))+
  stat_summary(fun=mean, geom="line",linewidth=1) +
  geom_jitter(width=0.1,size=1)+
  scale_y_continuous()+
  scale_x_continuous(limits=c(50,90),name="Acyl Chain Carbons")+
  ylab("Lipid Relative Abundance (%)")+
  scale_fill_manual(values=c("#8BD1CB","#006198"),breaks = names(km_cluster_colors)[c(1,2)])+
  scale_color_manual(values=c("#8BD1CB","#006198"),breaks = names(km_cluster_colors)[c(1,2)])+
  as_md_theme(my_theme)+theme(legend.position = "none")+ggtitle("")
plot(cl_species_chainlength_wt_kmeans.1.2.p)

cl_species_chainlength_wt_kmeans.3.p<-ggplot(ld_relabund.m_complete.clmlcl.kmavg_chain[ld_relabund.m_complete.clmlcl.kmavg_chain$variable %in% c("WT","WT_1","WT_2") & ld_relabund.m_complete.clmlcl.kmavg_chain$km_cluster %in% c("km3"),],aes(x=chain_length,y=value_relabund,color=km_cluster,fill=km_cluster))+
  stat_summary(fun=mean, geom="line",linewidth=1) +
  geom_jitter(width=0.1,size=1)+
  scale_y_continuous(limits=c(0,0.4))+
  scale_x_continuous(limits=c(45,85),name="Total Carbons in Lipid Chains")+
  ylab("Lipid Relative Abundance (%)")+
  scale_fill_manual(values=km_cluster_colors,breaks = names(km_cluster_colors))+
  scale_color_manual(values=km_cluster_colors,breaks = names(km_cluster_colors))+
  as_md_theme(my_theme)+theme(legend.position = "none")+ggtitle("")
plot(cl_species_chainlength_wt_kmeans.3.p)

cl_species_chainlength_str_bykmean.p<-ggplot(ld_relabund.m_complete.clmlcl.kmavg_chain[ld_relabund.m_complete.clmlcl.kmavg_chain$km_cluster %in% c("km1","km2"),],aes(x=chain_length,y=value_relabund,color=strain,fill=strain))+
  stat_summary(fun=mean, geom="line",linewidth=1) +
  geom_jitter(width=0.1,size=1)+
  scale_y_continuous(limits=c(0,0.4))+
  scale_x_continuous(name="Total Carbons in Lipid Chains")+
  ylab("Lipid Relative Abundance (%)")+
  scale_fill_manual(values=colors_1,breaks = breaks_1,labels=labels_1)+
  scale_color_manual(values=colors_1,breaks = breaks_1,labels=labels_1)+
  as_md_theme(my_theme)+theme(legend.position = "none")+facet_wrap(~km_cluster,nrow=1)
plot(cl_species_chainlength_str_bykmean.p)

# ld_relabund.m.avg_complete.sum<-aggregate(value~chain_unsaturations+strain+km_cluster+chain_length,data=ld_relabund.m.avg_complete,FUN=sum)

cl_species_unsaturations_wt_bykmeans_all.p<-ggplot(ld_relabund.m_complete.clmlcl.kmavg_unsat[ld_relabund.m_complete.clmlcl.kmavg_unsat$variable %in% c("WT","WT_1","WT_2"),],aes(x=chain_unsaturations,y=value_relabund,color=km_cluster,fill=km_cluster))+
  stat_summary(fun=mean, geom="line",linewidth=1) +
  geom_jitter(width=0.1,size=1)+
  scale_y_continuous()+
  scale_x_continuous(limits=c(-1,20),name="Total Unsaturations in Lipid Chains")+
  ylab("Lipid Relative Abundance (%)")+
  scale_fill_manual(values=km_cluster_colors,breaks = names(km_cluster_colors))+
  scale_color_manual(values=km_cluster_colors,breaks = names(km_cluster_colors))+
  as_md_theme(my_theme)+theme(legend.position = "none")+ggtitle("")
plot(cl_species_unsaturations_wt_bykmeans_all.p)

cl_species_unsaturations_wt_bykmeans_1.2.p<-ggplot(ld_relabund.m_complete.clmlcl.kmavg_unsat[ld_relabund.m_complete.clmlcl.kmavg_unsat$strain %in% "WT" & ld_relabund.m_complete.clmlcl.kmavg_unsat$km_cluster %in% c("km1","km2"),],aes(x=chain_unsaturations,y=value_relabund,color=km_cluster,fill=km_cluster))+
  stat_summary(fun=mean, geom="line",linewidth=1) +
  geom_jitter(width=0.1,size=1)+
  scale_y_continuous(limits=c(0,0.4))+
  scale_x_continuous(limits=c(0,20),name="Acyl Chain Unsaturations")+
  ylab("Lipid Relative Abundance (%)")+
  scale_fill_manual(values=c("#8BD1CB","#006198"),breaks = names(km_cluster_colors)[c(1,2)])+
  scale_color_manual(values=c("#8BD1CB","#006198"),breaks = names(km_cluster_colors)[c(1,2)])+
  as_md_theme(my_theme)+theme(legend.position = "none")+ggtitle("")
plot(cl_species_unsaturations_wt_bykmeans_1.2.p)

cl_species_unsaturations_kmeans1.2_bystrain.p<-ggplot(ld_relabund.m_complete.clmlcl.kmavg_unsat[ld_relabund.m_complete.clmlcl.kmavg_unsat$km_cluster %in% c("km1","km2"),],aes(x=chain_unsaturations,y=value_relabund,color=strain,fill=strain))+
  stat_summary(fun=mean, geom="line",linewidth=1) +
  geom_jitter(width=0.1)+
  scale_y_continuous(limits=c(0,0.4))+
  scale_x_continuous(limits=c(0,20),name="Total Unsaturations in Lipid Chains")+
  ylab("Lipid Relative Abundance (%)")+
  scale_fill_manual(values=colors_1,breaks = breaks_1)+
  scale_color_manual(values=colors_1,breaks = breaks_1)+
  as_md_theme(my_theme)+theme(legend.position = "none")+facet_wrap(~km_cluster,nrow=1)+ggtitle("")
plot(cl_species_unsaturations_kmeans1.2_bystrain.p)

cl_species_unsaturations_kmeans3.p<-ggplot(ld_relabund.m_complete.clmlcl.kmavg_unsat[ld_relabund.m_complete.clmlcl.kmavg_unsat$km_cluster %in% c("km3"),],aes(x=chain_unsaturations,y=value_relabund,color=strain,fill=strain))+
  stat_summary(fun=sum, geom="line",linewidth=1) +
  stat_summary(fun=sum, geom="point") +
  scale_y_continuous(limits=c(0,0.4))+
  scale_x_continuous(limits=c(0,20),name="Total Unsaturations in Lipid Chains")+
  ylab("Lipid Relative Abundance (%)")+
  scale_fill_manual(values=colors_1,breaks = breaks_1)+
  scale_color_manual(values=colors_1,breaks = breaks_1)+
  as_md_theme(my_theme)+theme(legend.position = "none")+facet_wrap(~km_cluster,nrow=1)+ggtitle("")
plot(cl_species_unsaturations_kmeans3.p)

ggsave(filename = "figures/lipidomics_kmeans_chainlength_wt_kmeans1.2.pdf",plot=cl_species_chainlength_wt_kmeans.1.2.p,width=2.5,height=2,units="in",device=cairo_pdf)
ggsave(filename = "figures/lipidomics_kmeans_chainlength_wt_kmeans3.pdf",plot=cl_species_chainlength_wt_kmeans.1.2.p,width=2.5,height=2,units="in",device=cairo_pdf)
ggsave(filename = "figures/lipidomics_kmeans_chainlength_str_kmeans.pdf",plot=cl_species_chainlength_str_bykmean.p,width=4,height=2,units="in",device=cairo_pdf)

ggsave(filename = "figures/lipidomics_kmeans_unsaturations_wt_kmeans1.2.pdf",plot=cl_species_unsaturations_wt_bykmeans_1.2.p,width=2.5,height=2,units="in",device=cairo_pdf)
ggsave(filename = "figures/lipidomics_kmeans_unsaturations_kmeans1.2_bystrain.pdf",plot=cl_species_unsaturations_kmeans1.2_bystrain.p,width=4,height=2,units="in",device=cairo_pdf)
ggsave(filename = "figures/lipidomics_kmeans_unsaturations_kmeans3.pdf",plot=cl_species_unsaturations_kmeans3.p,width=4,height=2,units="in",device=cairo_pdf)

# Investigating the contents of Kmeans clusters 2 and 4

ggplot(lipid_data.f.m[lipid_data.f.m$km_cluster %in% c("km3") & lipid_data.f.m$strain %in% c("WT","clsAclsB"),], aes(x=chain_length,y=value,color=strain))+geom_point(position=position_jitterdodge())+facet_wrap(~km_cluster)+
  scale_y_log10()

km3.datasub<-lipid_data.f.m[lipid_data.f.m$km_cluster %in% c("km3") & lipid_data.f.m$strain %in% c("WT","clsAclsB") & lipid_data.f.m$detect_mode %in% "NEG",]


km.datasub.dcast<-dcast(as.data.table(km3.datasub),formula = Metabolite_name_from_MSDIAL+Manually_confirmed_lipid+Ontology+lipid_size+chain_length+chain_unsaturations+km_cluster~variable,value.var="value")
km.datasub.dcast$difference_cAB1<-km.datasub.dcast[,"clsA_clsB"]/rowMeans2(as.matrix(km.datasub.dcast[,c("WT","WT_1","WT_2")]))
km.datasub.dcast$difference_cAB2<-km.datasub.dcast[,"clsA_clsB_1"]/rowMeans2(as.matrix(km.datasub.dcast[,c("WT","WT_1","WT_2")]))
km.datasub.dcast$difference_cAB3<-km.datasub.dcast[,"clsA_clsB_2"]/rowMeans2(as.matrix(km.datasub.dcast[,c("WT","WT_1","WT_2")]))

km.datasub.dcast.m<-melt.data.table(km.datasub.dcast,measure.vars = c("difference_cAB1","difference_cAB2","difference_cAB3"))

km_twochain_length.fc.p<-ggplot(km.datasub.dcast.m[km.datasub.dcast.m$Ontology %in% c("PE","PG","PS"),],aes(x=chain_length,y=log2(value),fill=km_cluster))+geom_bar(stat="summary",fun=mean,position=position_dodge(width=0.9))+stat_summary(geom="errorbar",fun.data=mean_sd, width=0.5, linewidth=0.5,position=position_dodge(width=0.9))+
  scale_y_continuous(name="log2(FC) (Δ*clsA*Δ*clsB*/WT)",limits=c(-2.5,2.5))+
  scale_x_continuous(name="Acyl Chain Length")+
  scale_fill_manual(values=c("black","darkgray"))+
  as_md_theme(my_theme)+theme(legend.position = "none")+
  facet_wrap(~Ontology)
plot(km_twochain_length.fc.p)

km_twochain_unsat.fc.p<-ggplot(km.datasub.dcast.m[km.datasub.dcast.m$Ontology %in% c("PE","PG","PS"),],aes(x=chain_unsaturations,y=log2(value),fill=km_cluster))+geom_bar(stat="summary",fun=mean,position=position_dodge(width=0.9))+stat_summary(geom="errorbar",fun.data=mean_sd, width=0.5, linewidth=0.5,position=position_dodge(width=0.9))+
  scale_y_continuous(name="log2(FC) (Δ*clsA*Δ*clsB*/WT)")+
  scale_x_continuous(name="Acyl Chain Unsatuations")+
  scale_fill_manual(values=c("black","darkgray"))+
  as_md_theme(my_theme)+theme(legend.position = "none")+
  facet_wrap(~Ontology)
plot(km_twochain_unsat.fc.p)

ggsave(filename = "figures/lipidomics_kmeans_length_fc_kmeans3.pdf",plot=km_twochain_length.fc.p,width=5,height=2,units="in",device=cairo_pdf)
ggsave(filename = "figures/lipidomics_kmeans_unsat_fc_kmeans3.pdf",plot=km_twochain_unsat.fc.p,width=3,height=2,units="in",device=cairo_pdf)

sphingo_length.fc.p<-ggplot(lipid_data.f.m[lipid_data.f.m$Ontology %in% c("SPB","SM"),],aes(x=chain_length,y=log2(value_fc)))+geom_bar(stat="summary",fun=mean,position=position_dodge(width=0.9))+stat_summary(geom="errorbar",fun.data=mean_sd, width=0.5, linewidth=0.5,position=position_dodge(width=0.9))+
  scale_y_continuous(name="log2(FC) (Δ*clsA*Δ*clsB*/WT)")+
  scale_x_continuous(name="Acyl Chain Length")+
  scale_fill_manual(values=c("black","darkgray"))+
  as_md_theme(my_theme)+theme(legend.position = "none")+
  facet_wrap(~Ontology)+facet_wrap(~strain)
plot(sphingo_length.fc.p)

## Biplot
# biplot based on peak size (not rel abundance) data seems better (not much of difference though)

lipid_data.biplot<-as.data.frame(t(lipid_data.f[,c("Metabolite_name_from_MSDIAL","WT","WT_1","WT_2","clsA","clsA_1","clsA_2","clsB","clsB_1","clsB_2","clsA_clsB","clsA_clsB_1","clsA_clsB_2")]))
colnames(lipid_data.biplot)<-make.unique(unlist(lipid_data.biplot[1,]))
lipid_data.biplot<-lipid_data.biplot[-1,]
lipid_data.biplot$strain<-gsub("[0-9]+","",rownames(lipid_data.biplot)) %>% gsub("_","",.) %>% factor(.,levels=breaks_1)
lipid_data.biplot[,1:272]<-lapply(lipid_data.biplot[,1:272],as.numeric)

lipid_data_relabund.biplot<-as.data.frame(t(lipid_data_relabund[,c("Metabolite_name_from_MSDIAL","WT","WT_1","WT_2","clsA","clsA_1","clsA_2","clsB","clsB_1","clsB_2","clsA_clsB","clsA_clsB_1","clsA_clsB_2")]))
colnames(lipid_data_relabund.biplot)<-make.unique(unlist(lipid_data_relabund.biplot[1,]))
lipid_data_relabund.biplot<-lipid_data_relabund.biplot[-1,]
lipid_data_relabund.biplot$strain<-gsub("[0-9]+","",rownames(lipid_data_relabund.biplot)) %>% gsub("_","",.) %>% factor(.,levels=breaks_1)
lipid_data_relabund.biplot[,1:271]<-lapply(lipid_data_relabund.biplot[,1:271],as.numeric)


lipid_all.pca<-prcomp(lipid_data.biplot[,1:272],center=TRUE,scale.=TRUE)
summary(lipid_all.pca)
lipid_all_relabund.pca<-prcomp(lipid_data_relabund.biplot[,1:271],center=TRUE,scale.=TRUE)
summary(lipid_all_relabund.pca)

# geom_mark_ellipse uses the Khachiyan algorithm (ggforce) - calculates minimum volume to enclose objects

lipid_all.pca12.p<-ggbiplot(lipid_all.pca, choices=1:2,data=lipid_data.biplot,groups=lipid_data.biplot$strain,var.axes = FALSE)+
  geom_mark_ellipse(aes(fill = lipid_data.biplot$strain,color = lipid_data.biplot$strain)) +
  scale_x_continuous(limits=c(-3,3))+
  scale_y_continuous(limits=c(-3,3))+
  labs(color="Strain")+ylab("PC 2 (19.2%)")+xlab("PC 1 (41.1%)")+
  scale_color_manual(values=colors_1,breaks=breaks_1,labels=labels_1)+
  scale_fill_manual(values=colors_1,breaks=breaks_1,labels=labels_1)+
  as_md_theme(my_theme)+theme(legend.position="none")+ggtitle("")
plot(lipid_all.pca12.p)

# lipid_all_relabund.pca12.p<-ggbiplot(lipid_all_relabund.pca, choices=1:2,data=lipid_data_relabund.biplot,groups=lipid_data_relabund.biplot$strain,var.axes = FALSE)+
#  labs(color="Strain")+ylab("PC 2 (20.2%)")+xlab("PC 1 (41.7%)")+
#  scale_color_manual(values=colors_1,breaks=breaks_1,labels=labels_1)+
#  as_md_theme(my_theme)
# plot(lipid_all_relabund.pca12.p)

lipid_all.pca23.p<-ggbiplot(lipid_all.pca, choices=2:3,data=lipid_data.biplot,groups=lipid_data_relabund.biplot$strain,var.axes = FALSE)+
  geom_mark_ellipse(aes(fill = lipid_data.biplot$strain,color = lipid_data.biplot$strain)) +
  scale_x_continuous(limits=c(-3,3))+
  scale_y_continuous(limits=c(-3,3))+
  labs(color="Strain")+ylab("PC 3 (10.5%)")+xlab("PC 2 (19.2%)")+
  scale_color_manual(values=colors_1,breaks=breaks_1,labels=labels_1)+
  scale_fill_manual(values=colors_1,breaks=breaks_1,labels=labels_1)+
  as_md_theme(my_theme)+theme(legend.position = "none")+ggtitle("")
plot(lipid_all.pca23.p)

# Biplot excluding cardiolipin and monolysocardiolipin
lipid_data_nocl.biplot<-as.data.frame(t(lipid_data.f[lipid_data.f$Ontology %!in% c("CL","MLCL"),c("Metabolite_name_from_MSDIAL","WT","WT_1","WT_2","clsA","clsA_1","clsA_2","clsB","clsB_1","clsB_2","clsA_clsB","clsA_clsB_1","clsA_clsB_2")]))
colnames(lipid_data_nocl.biplot)<-make.unique(unlist(lipid_data_nocl.biplot[1,]))
lipid_data_nocl.biplot<-lipid_data_nocl.biplot[-1,]
lipid_data_nocl.biplot$strain<-gsub("[0-9]+","",rownames(lipid_data_nocl.biplot)) %>% gsub("_","",.) %>% factor(.,levels=breaks_1)
lipid_data_nocl.biplot[,1:161]<-lapply(lipid_data_nocl.biplot[,1:161],as.numeric)

lipid_all_nocl.pca<-prcomp(lipid_data_nocl.biplot[,1:161],center=TRUE,scale.=TRUE)
summary(lipid_all_nocl.pca)

lipid_all_nocl.pca12.p<-ggbiplot(lipid_all_nocl.pca, choices=1:2,data=lipid_data_nocl.biplot,groups=lipid_data_nocl.biplot$strain,var.axes = FALSE)+
  geom_mark_ellipse(aes(fill = lipid_data.biplot$strain,color = lipid_data.biplot$strain)) +
  scale_x_continuous(limits=c(-3,3))+
  scale_y_continuous(limits=c(-3,3))+
  labs(color="Strain")+ylab("PC 2 (16.0%)")+xlab("PC 1 (39.1%)")+
  scale_color_manual(values=colors_1,breaks=breaks_1,labels=labels_1)+
  scale_fill_manual(values=colors_1,breaks=breaks_1,labels=labels_1)+
  as_md_theme(my_theme)+theme(legend.position = "none")+ggtitle("")
plot(lipid_all_nocl.pca12.p)

lipid_all_nocl.pca14.p<-ggbiplot(lipid_all_nocl.pca, choices=c(1,4),data=lipid_data_nocl.biplot,groups=lipid_data_nocl.biplot$strain,var.axes = FALSE)+
  geom_mark_ellipse(aes(fill = lipid_data.biplot$strain,color = lipid_data.biplot$strain)) +
  scale_x_continuous(limits=c(-3,3))+
  scale_y_continuous(limits=c(-3,3))+
  labs(color="Strain")+ylab("PC 4 (9.3%)")+xlab("PC 1 (39.1%)")+
  scale_color_manual(values=colors_1,breaks=breaks_1,labels=labels_1)+
  scale_fill_manual(values=colors_1,breaks=breaks_1,labels=labels_1)+
  as_md_theme(my_theme)+theme(legend.position = "none")+ggtitle("")
plot(lipid_all_nocl.pca14.p)

ggsave(filename = "figures/lipidomics_biplot_pca12.pdf",plot=lipid_all.pca12.p,width=3,height=3,units="in",device=cairo_pdf)
ggsave(filename = "figures/lipidomics_biplot_pca23.pdf",plot=lipid_all.pca23.p,width=3,height=3,units="in",device=cairo_pdf)
ggsave(filename = "figures/lipidomics_biplot_nocl_pca12.pdf",plot=lipid_all_nocl.pca12.p,width=3,height=3,units="in",device=cairo_pdf)
ggsave(filename = "figures/lipidomics_biplot_nocl_pca14.pdf",plot=lipid_all_nocl.pca14.p,width=3,height=3,units="in",device=cairo_pdf)

lipid_data_ontology.biplot<-as.data.frame(t(lipid_data.f.ontology[,c("Ontology","WT","WT_1","WT_2","clsA","clsA_1","clsA_2","clsB","clsB_1","clsB_2","clsA_clsB","clsA_clsB_1","clsA_clsB_2")]))
colnames(lipid_data_ontology.biplot)<-make.unique(unlist(lipid_data_ontology.biplot[1,]))
lipid_data_ontology.biplot<-lipid_data_ontology.biplot[-1,]
lipid_data_ontology.biplot$strain<-gsub("[0-9]+","",rownames(lipid_data_ontology.biplot)) %>% gsub("_","",.) %>% factor(.,levels=breaks_1)
lipid_data_ontology.biplot[,1:13]<-lapply(lipid_data_ontology.biplot[,1:13],as.numeric)

lipid_ontology.pca<-prcomp(lipid_data_ontology.biplot[,1:13],center=TRUE,scale.=TRUE)
summary(lipid_ontology.pca)

lipid_ontology.pca12.p<-ggbiplot(lipid_ontology.pca, choices=1:2,data=lipid_data_ontology.biplot,groups=lipid_data_ontology.biplot$strain)+
  geom_mark_ellipse(aes(fill = lipid_data_ontology.biplot$strain,color = lipid_data_ontology.biplot$strain)) +
  scale_x_continuous(limits=c(-3,3))+
  scale_y_continuous(limits=c(-3,3))+
  labs(color="Strain")+ylab("PC 2 (26.38%)")+xlab("PC 1 (36.9%)")+
  scale_color_manual(values=colors_1,breaks=breaks_1,labels=labels_1)+
  scale_fill_manual(values=colors_1,breaks=breaks_1,labels=labels_1)+
  as_md_theme(my_theme)+ggtitle("Ontology")
plot(lipid_ontology.pca12.p)

lipid_ontology.pca14.p<-ggbiplot(lipid_ontology.pca, choices=c(1,4),data=lipid_data_ontology.biplot,groups=lipid_data_ontology.biplot$strain)+
  geom_mark_ellipse(aes(fill = lipid_data_ontology.biplot$strain,color = lipid_data_ontology.biplot$strain)) +
  scale_x_continuous(limits=c(-3,3))+
  scale_y_continuous(limits=c(-3,3))+
  labs(color="Strain")+ylab("PC 4 (7.0%)")+xlab("PC 1 (36.9%)")+
  scale_color_manual(values=colors_1,breaks=breaks_1,labels=labels_1)+
  scale_fill_manual(values=colors_1,breaks=breaks_1,labels=labels_1)+
  as_md_theme(my_theme)+ggtitle("Ontology")
plot(lipid_ontology.pca14.p)

ggsave(filename = "figures/lipidomics_biplot_ontology_pca12.pdf",plot=lipid_ontology.pca12.p,width=8,height=6,units="in",device=cairo_pdf)
ggsave(filename = "figures/lipidomics_biplot_ontology_pca14.pdf",plot=lipid_ontology.pca14.p,width=8,height=6,units="in",device=cairo_pdf)

#####
# LEfSE
# Khleborodova A (2024). lefser: R implementation of the LEfSE method for microbiome biomarker discovery. R package version 1.14.0, https://github.com/waldronlab/lefser.

str(lipid_data_relabund.cons)

wt.clsAB.design<-as.data.frame(matrix(nrow=6,ncol=2))
colnames(wt.clsAB.design)<-c("group","Class")
wt.clsAB.design$Class<-c("WT","WT","WT","clsAclsB","clsAclsB","clsAclsB")
wt.clsAB.design$group<-c("WT","WT_1","WT_2","clsA_clsB","clsA_clsB_1","clsA_clsB_2")
write.table(wt.clsAB.design,file="lefse/wt_clsAB_design.tsv")

lipid_data.f



lipid_lefse.summexp<-SummarizedExperiment(lipid_data_relabund.cons[,c("WT","WT_1","WT_2","clsA_clsB","clsA_clsB_1","clsA_clsB_2")],colData = testtest)


testtest.lefse<-lefser(lipid_lefse.summexp, groupCol = "group",kruskal.threshold = 0.1,wilcox.threshold = 0.1,lda.threshold = 0)

lefserPlot(testtest.lefse)

# No significant results from LEfSe

#####
# Presence vs. absence (MAKE NOTE OF THE LIPIDS THAT APPEAR IN OTHER STRAIN BACKGROUNDS)
# Results: These lipids that were not present in some WT replicates were present in others. I'm guessing that they are just near the limit of detection
# lipid species that are sometimes present in WT
wt_absent<-lipid_data.f.m[lipid_data.f.m$variable %in% c("WT","WT_1","WT_2") & lipid_data.f.m$value_logical %in% FALSE,]$Manually_confirmed_lipid

# fix Y axis log10 labels
wt_absent.p<-ggplot(data=lipid_data.f.m[lipid_data.f.m$Manually_confirmed_lipid %in% wt_absent,],aes(x=Metabolite_name_from_MSDIAL,y=value,color=strain))+geom_boxplot()+geom_point(position=position_jitterdodge(dodge.width=0.4,jitter.height = 0,jitter.width = 0.2))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x,n=4), labels = trans_format("log10", label_math(10^.x)),limits=c(1E-6,1E-3))+
  annotation_logticks(base=10,sides = "l",scaled=TRUE,short = unit(0.05, "cm"),mid = unit(0.1, "cm"),long = unit(0.2, "cm"),)+
  xlab("Lipid Species")+ylab("Peak Area")+
  scale_color_manual(breaks=breaks_1,labels=labels_1,values=colors_1)+
  as_md_theme(my_theme2)
plot(wt_absent.p)
ggsave(filename = "figures/lipidomics_wtabsent.pdf",width=8,height=6,units="in",device=cairo_pdf)

# Polish so that all bars are the same width (may involve some precalculations)
wt_absent_logical.p<-ggplot(data=lipid_data.f.m[lipid_data.f.m$Manually_confirmed_lipid %in% wt_absent & lipid_data.f.m$norm_status %!in% "norm" & lipid_data.f.m$value_logical %in% TRUE,],aes(x=Manually_confirmed_lipid,fill=strain))+geom_bar(position=position_dodge())+
  xlab("Lipid Species")+ylab("Positive Replicates")+
  scale_fill_manual(breaks=breaks_1,labels=labels_1,values=colors_1)+
  as_md_theme(my_theme2)
plot(wt_absent_logical.p)
ggsave(filename = "figures/lipidomics_wtabsent_logical.pdf",width=8,height=6,units="in",device=cairo_pdf)


m.dlcl.p<-ggplot(data=lipid_data.f.m[lipid_data.f.m$Ontology %in% c("DLCL") &  lipid_data.f.m$value !=0,],aes(x=reorder(Manually_confirmed_lipid,-value),y=value,fill=strain,color=NA))+geom_boxplot(position=position_dodge(width=0.7))+
  geom_boxplot(aes(color = strain),
               fatten = NULL, fill = NA, coef = 0, outlier.alpha = 0,
               show.legend = F,position=position_dodge(width=0.7))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x,n=4), labels = trans_format("log10", label_math(10^.x)),limits=c(1E-6,1E-2))+
  annotation_logticks(base=10,sides = "l",scaled=TRUE,short = unit(0.05, "cm"),mid = unit(0.1, "cm"),long = unit(0.2, "cm"),)+
  xlab("Lipid Species")+ylab("Log10(Peak Area)")+
  scale_fill_manual(breaks=breaks_1,labels=labels_1,values=colors_1)+
    scale_color_manual(breaks=breaks_1,labels=labels_1,values=colors_1)+
  as_md_theme(my_theme2)
plot(m.dlcl.p)

dlcl_regression.df<-dcast(as.data.table(lipid_data.f.m.avg),formula = Metabolite_name_from_MSDIAL+lipid_size+Ontology ~ strain, value.var='value',fun.aggregate = mean)

m.dlcl_regression.p<-ggplot(data=dlcl_regression.df[dlcl_regression.df$Ontology %in% c("DLCL"),],aes(x=WT,y=clsAclsB))+geom_point()+
  geom_abline()+
  scale_y_log10(name="DLCL Peak Area clsAclsB",breaks = trans_breaks("log10", function(x) 10^x,n=4), labels = trans_format("log10", label_math(10^.x)),limits=c(1E-6,1E-2))+
  scale_x_log10(name="DLCL Peak Area WT",breaks = trans_breaks("log10", function(x) 10^x,n=4), labels = trans_format("log10", label_math(10^.x)),limits=c(1E-6,1E-2))+
  annotation_logticks(base=10,sides = "bl",scaled=TRUE,short = unit(0.05, "cm"),mid = unit(0.1, "cm"),long = unit(0.2, "cm"),)+
  scale_fill_manual(breaks=breaks_1,labels=labels_1,values=colors_1)+
  scale_color_manual(breaks=breaks_1,labels=labels_1,values=colors_1)+
  as_md_theme(my_theme)
plot(m.dlcl_regression.p)

m.dlcl_fc.p<-ggplot(data=lipid_data.f.m[lipid_data.f.m$Ontology %in% c("DLCL") & lipid_data.f.m$strain %in% "clsAclsB",],aes(x=reorder(lipid_size,-value_fc),y=log2(value_fc)),color="black")+geom_bar(stat="summary",fun=mean,position=position_dodge(width=0.9),fill="black")+stat_summary(geom="errorbar",fun.data=mean_sd, width=0.5, linewidth=0.5,position=position_dodge(width=0.9))+
  scale_y_continuous(name="DLCL Peak Area clsAclsB",limits=c(-4,4),breaks=c(-4:4))+
  scale_fill_manual(breaks=breaks_1,labels=labels_1,values=colors_1)+
  scale_color_manual(breaks=breaks_1,labels=labels_1,values=colors_1)+
  as_md_theme(my_theme2)
plot(m.dlcl_fc.p)

ggsave(filename = "figures/dlcl_fc.pdf",plot=m.dlcl_fc.p,width=4,height=2,units="in",device=cairo_pdf)


m.dlcl_regression_fc.p<-ggplot(data=dlcl_regression.df_fc[dlcl_regression.df_fc$Ontology %in% c("DLCL"),],aes(x=WT,y=log2(clsAclsB)))+geom_point()+
  geom_abline()+

  scale_fill_manual(breaks=breaks_1,labels=labels_1,values=colors_1)+
  scale_color_manual(breaks=breaks_1,labels=labels_1,values=colors_1)+
  as_md_theme(my_theme)
plot(m.dlcl_regression_fc.p)

m.dlcl.ontology.p<-ggplot(data=lipid_data.f.m.ontology[lipid_data.f.m.ontology$Ontology %in% c("DLCL"),],aes(x=strain,y=value,color=strain))+geom_boxplot(position=position_dodge(width=0.4))+geom_point(position=position_jitterdodge(dodge.width=0.4,jitter.height = 0,jitter.width = 0.2))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x,n=4), labels = trans_format("log10", label_math(10^.x)),limits=c(1E-2,1E-1))+
  annotation_logticks(base=10,sides = "l",scaled=TRUE,short = unit(0.05, "cm"),mid = unit(0.1, "cm"),long = unit(0.2, "cm"),)+
  xlab("Strain")+ylab("Log10(DLCL Peak Area)")+
  scale_color_manual(breaks=breaks_1,labels=labels_1,values=colors_1)+
  as_md_theme(my_theme2)
plot(m.dlcl.ontology.p)

ggsave(plot=m.dlcl.p,filename = "figures/lipidomics_dlcl_bystrain.pdf",width=7,height=2,units="in",device=cairo_pdf)
ggsave(plot=m.dlcl.ontology.p,filename = "figures/lipidomics_dlcl_bystrain_ontology.pdf",width=2.5,height=2,units="in",device=cairo_pdf)



dlcl.lm<-glm(log10(value)~strain+lipid_size,data=lipid_data.f.m[lipid_data.f.m$Ontology %in% c("DLCL") & lipid_data.f.m$norm_status %!in% "norm" & lipid_data.f.m$value >0,])
summary(dlcl.lm)
drop1(dlcl.lm,test = "F")
# P = 0.00029

strain.size.lm<-glm(chain_length~strain,data=lipid_data.f.m[lipid_data.f.m$Ontology %in% c("DLCL") & lipid_data.f.m$norm_status %!in% "norm" & lipid_data.f.m$value >0,])
drop1(strain.size.lm,test="F")
# ns

######
## sub chain analysis
lipid_data_chains<-gsub(".*PS ","",lipid_data.f.m$Metabolite_name_from_MSDIAL) %>% gsub(".*PE ","",.) %>% gsub(".*PG ","",.)  %>% gsub(".*PC ","",.)   %>% gsub(".*PA ","",.) %>% gsub(".*CL ","",.) %>% str_split_fixed(.,pattern="_",n=4) %>% as.data.frame(.) %>% setnames(.,c("chain1","chain2","chain3","chain4"))
lipid_data.f.m_chain<-cbind(lipid_data.f.m,lipid_data_chains)
lipid_data.f.m_chain<-lipid_data.f.m_chain[lipid_data.f.m_chain$norm_status %in% "norm" & lipid_data.f.m_chain$Ontology %!in% "FA",c("Metabolite_name_from_MSDIAL","strain","Ontology","lipid_size","chain_length","chain_unsaturations","km_cluster","chain1","chain2","chain3","chain4","value")]
lipid_data.f.m_chain.avg<-aggregate(value ~ Metabolite_name_from_MSDIAL+Ontology+lipid_size+chain_length+chain_unsaturations+km_cluster+strain+chain1+chain2+chain3+chain4,data=lipid_data.f.m_chain,FUN=mean)

lipid_data.f.m_chain.avg[lipid_data.f.m_chain.avg$Ontology %in% "MLCL" & lipid_data.f.m_chain.avg$chain2 == "",]
# PE has three species 31:0, :31:1 and 32:0 that don't have sub species
# PG all have 2
# PS one at 31:3 without second tail
# DLCL all have 2
# all MLCL (3 chains) have three chains annotated
# - CL (4 chains) has one with only two chains "63:3|31:0_32:3"
lipid_data.f.m_chain.avg<-lipid_data.f.m_chain.avg[lipid_data.f.m_chain.avg$lipid_size %!in% c("63:3|31:0_32:3") & lipid_data.f.m_chain.avg$chain1 %!in% c("31:0","31:1","32:0","31:3"),]

lipid_chain_avg.m<-melt.data.table(as.data.table(lipid_data.f.m_chain.avg),id.vars = c("strain","Ontology","lipid_size","chain_length","chain_unsaturations","km_cluster","value"),measure.vars = c("chain1","chain2","chain3","chain4"),value.name="chain_length_ind")
lipid_chain_avg.m<-lipid_chain_avg.m[lipid_chain_avg.m$chain_length_ind !="",]
lipid_chain_avg.m$value_cur<-lipid_chain_avg.m$value^(1/3)
# rownames(lipid_chain_avg.m.curoot)<-make.unique(lipid_chain_avg.m.curoot$lipid_size)

####
subsetteddata<-lipid_chain_avg.m[lipid_chain_avg.m$strain %!in% "WT"& lipid_chain_avg.m$Ontology %in% c("CL","MLCL"),-"value"]
subsetteddata<-as.data.table(aggregate(value_cur ~ Ontology+km_cluster+strain+chain_length_ind,data=subsetteddata,FUN=mean))
subsetteddata<-dcast(subsetteddata,formula = chain_length_ind ~ strain, value.var='value_cur',fun.aggregate = mean)
subsetteddata[,2:ncol(subsetteddata)]<-lapply(subsetteddata[,2:ncol(subsetteddata)],function(x) as.numeric(gsub("NaN",-1,x)))


abund_cols.cu.ind<-colorRamp2(c(1.25, 1, 0), brewer.pal(n=3, name="RdBu"))

ind.heatmap.p<-Heatmap(as.matrix(subsetteddata[,2:ncol(subsetteddata)]),name="CuRt(Fold-Change)",row_labels = subsetteddata$chain_length_ind,col=abund_cols.cu.ind)
plot(ind.heatmap.p)

# pdf(file = "figures/lipidomics_lipid_ind_chain_heatmap.pdf",
#    width = 8,height = 6) # size in inches
#plot(ind.heatmap.p)
#dev.off()




subsetteddata<-dcast(subsetteddata,formula = strain+Ontology+chain_length+chain_unsaturations+km_cluster+lipid_size ~ chain_length_ind, value.var='value_cur',fun.aggregate = mean)
subsetteddata[,7:ncol(subsetteddata)]<-lapply(subsetteddata[,7:ncol(subsetteddata)],function(x) as.numeric(gsub("NaN",-1,x)))

names(lipid_ontology_colors) <- levels(subsetteddata$Ontology)
names(colors_1)<-breaks_1

col <- list(Ontology = lipid_ontology_colors,
            km_cluster=km_cluster_colors,strain=colors_1)

row_ha<-rowAnnotation(Ontology=subsetteddata$Ontology, km_cluster=subsetteddata$km_cluster, strain=subsetteddata$strain, col=col)

abund_cols.cu<-colorRamp2(c(2, 1, 0,-1), c(brewer.pal(n=3, name="RdBu"),"gray"))
labels_strains_heatmap<-gsub("[[:punct:]]","",labels_strains)

unsat.heatmap.p<-Heatmap(as.matrix(subsetteddata[,7:ncol(subsetteddata)]),name="CuRt(Fold-Change)",left_annotation = row_ha,col=abund_cols.cu,row_labels = subsetteddata$chain_length_ind)
plot(unsat.heatmap.p)

# pdf(file = "figures/lipidomics_lipid_indchain_bylipidsize_heatmap.pdf",
#    width = 8,height = 6) # size in inches
# plot(unsat.heatmap.p)
# dev.off()

# several MLCL and CL chains only found in the clsA mutant therefore clsB associated (20:3, 16:2, 24:0, 22:1, 18:1, 16:1)
# one chain only found in clsB mutant (clsA associated)

unique(lipid_data.f.m_chain[lipid_data.f.m_chain$lipid_size %in% "20:3",]$Ontology) # LPA
unique(lipid_data.f.m_chain[lipid_data.f.m_chain$chain1 %in% "20:3",]$Ontology) # LPA
unique(lipid_data.f.m_chain[lipid_data.f.m_chain$chain2 %in% "20:3",]$Ontology) # CL
unique(lipid_data.f.m_chain[lipid_data.f.m_chain$chain3 %in% "20:3",]$Ontology) # CL
unique(lipid_data.f.m_chain[lipid_data.f.m_chain$chain4 %in% "20:3",]$Ontology) # CL

unique(lipid_data.f.m_chain[lipid_data.f.m_chain$lipid_size %in% "16:2",]$Ontology) # LPA LPS
unique(lipid_data.f.m_chain[lipid_data.f.m_chain$chain1 %in% "16:2",]$Ontology) # DLCL LPA LPS
unique(lipid_data.f.m_chain[lipid_data.f.m_chain$chain2 %in% "16:2",]$Ontology) # CL DLCL PS
unique(lipid_data.f.m_chain[lipid_data.f.m_chain$chain3 %in% "16:2",]$Ontology) # none
unique(lipid_data.f.m_chain[lipid_data.f.m_chain$chain4 %in% "16:2",]$Ontology) # CL

unique(lipid_data.f.m_chain[lipid_data.f.m_chain$lipid_size %in% "24:0",]$Ontology) # none
unique(lipid_data.f.m_chain[lipid_data.f.m_chain$chain1 %in% "24:0",]$Ontology) # PG
unique(lipid_data.f.m_chain[lipid_data.f.m_chain$chain2 %in% "24:0",]$Ontology) # None
unique(lipid_data.f.m_chain[lipid_data.f.m_chain$chain3 %in% "24:0",]$Ontology) # CL
unique(lipid_data.f.m_chain[lipid_data.f.m_chain$chain4 %in% "24:0",]$Ontology) # none

unique(lipid_data.f.m_chain[lipid_data.f.m_chain$lipid_size %in% "22:1",]$Ontology) # none
unique(lipid_data.f.m_chain[lipid_data.f.m_chain$chain1 %in% "22:1",]$Ontology) # MLCL
unique(lipid_data.f.m_chain[lipid_data.f.m_chain$chain2 %in% "22:1",]$Ontology) # None
unique(lipid_data.f.m_chain[lipid_data.f.m_chain$chain3 %in% "22:1",]$Ontology) # MLCL
unique(lipid_data.f.m_chain[lipid_data.f.m_chain$chain4 %in% "22:1",]$Ontology) # none

unique(lipid_data.f.m_chain[lipid_data.f.m_chain$lipid_size %in% "18:1",]$Ontology) # LPE
unique(lipid_data.f.m_chain[lipid_data.f.m_chain$chain1 %in% "18:1",]$Ontology) # LPE PE
unique(lipid_data.f.m_chain[lipid_data.f.m_chain$chain2 %in% "18:1",]$Ontology) # CL PE
unique(lipid_data.f.m_chain[lipid_data.f.m_chain$chain3 %in% "18:1",]$Ontology) # CL
unique(lipid_data.f.m_chain[lipid_data.f.m_chain$chain4 %in% "18:1",]$Ontology) # CL

unique(lipid_data.f.m_chain[lipid_data.f.m_chain$lipid_size %in% "16:1",]$Ontology) # LPA LPS
unique(lipid_data.f.m_chain[lipid_data.f.m_chain$chain1 %in% "16:1",]$Ontology) # DLCL LPA LPS
unique(lipid_data.f.m_chain[lipid_data.f.m_chain$chain2 %in% "16:1",]$Ontology) # CL PE
unique(lipid_data.f.m_chain[lipid_data.f.m_chain$chain3 %in% "16:1",]$Ontology) # none
unique(lipid_data.f.m_chain[lipid_data.f.m_chain$chain4 %in% "16:1",]$Ontology) # none

# clsB mutant only (clsA associated)
unique(lipid_data.f.m_chain[lipid_data.f.m_chain$lipid_size %in% "21:0",]$Ontology) # none
unique(lipid_data.f.m_chain[lipid_data.f.m_chain$chain1 %in% "21:0",]$Ontology) # DLCL
unique(lipid_data.f.m_chain[lipid_data.f.m_chain$chain2 %in% "21:0",]$Ontology) # MLCL
unique(lipid_data.f.m_chain[lipid_data.f.m_chain$chain3 %in% "21:0",]$Ontology) # none
unique(lipid_data.f.m_chain[lipid_data.f.m_chain$chain4 %in% "21:0",]$Ontology) # none

# pdf(file = "figures/lipidomics_lipid_sqrt_fc_heatmap.pdf",
#    width = 8,height = 6) # size in inches
# plot(lipidfc.sqrt.heatmap.p)
# dev.off()

# pdf(file = "figures/lipidomics_lipid_curt_fc_heatmap.pdf",
#    width = 6,height = 5) # size in inches
# plot(lipidfc.curt.heatmap.p)
#  dev.off()




######
## Comparing NEG and POS mode PE

summary(lipid_data.f.m[lipid_data.f.m$variable %in% "WT" & lipid_data.f.m$Ontology %in% "PE" & lipid_data.f.m$detect_mode %in% "POS",]$Manually_confirmed_lipid %in% lipid_data.f.m[lipid_data.f.m$variable %in% "WT" & lipid_data.f.m$Ontology %in% "PE" & lipid_data.f.m$detect_mode %in% "NEG",]$Manually_confirmed_lipid)
# 6 POS %in% NEG, 9 not
# 6 NEG %in% POS, 16 not (NEG dataset is much larger)

lipid_pe_pos<-lipid_data.f.m[lipid_data.f.m$variable %in% "WT" & lipid_data.f.m$Ontology %in% "PE" & lipid_data.f.m$detect_mode %in% "POS",]
lipid_pe_neg<-lipid_data.f.m[lipid_data.f.m$variable %in% "WT" & lipid_data.f.m$Ontology %in% "PE" & lipid_data.f.m$detect_mode %in% "NEG",]
lipid_pe_overlap<-lipid_data.f.m[lipid_data.f.m$variable %in% "WT" & lipid_data.f.m$Ontology %in% "PE" & lipid_data.f.m$Metabolite_name_from_MSDIAL %in% lipid_pe_neg$Metabolite_name_from_MSDIAL & lipid_data.f.m$Metabolite_name_from_MSDIAL %in% lipid_pe_pos$Metabolite_name_from_MSDIAL,-c("Formula","km_cluster")]
lipid_pe_overlap<-dcast(data=lipid_pe_overlap,formula = Metabolite_name_from_MSDIAL+Manually_confirmed_lipid+lipid_size+chain_length+chain_unsaturations+norm_status+strain~detect_mode,value.var="value")
lipid_pe_overlap$NEG/lipid_pe_overlap$POS

lipid_ps_pos<-lipid_data.f.m[lipid_data.f.m$variable %in% "WT" & lipid_data.f.m$Ontology %in% "PS" & lipid_data.f.m$detect_mode %in% "POS",]
lipid_ps_neg<-lipid_data.f.m[lipid_data.f.m$variable %in% "WT" & lipid_data.f.m$Ontology %in% "PS" & lipid_data.f.m$detect_mode %in% "NEG",]
lipid_ps_overlap<-lipid_data.f.m[lipid_data.f.m$variable %in% "WT" & lipid_data.f.m$Ontology %in% "PS" & lipid_data.f.m$Metabolite_name_from_MSDIAL %in% lipid_ps_neg$Metabolite_name_from_MSDIAL & lipid_data.f.m$Metabolite_name_from_MSDIAL %in% lipid_ps_pos$Metabolite_name_from_MSDIAL,-c("Formula","km_cluster")]
lipid_ps_overlap<-dcast(data=lipid_ps_overlap,formula = Metabolite_name_from_MSDIAL+Manually_confirmed_lipid+lipid_size+chain_length+chain_unsaturations+norm_status+strain~detect_mode,value.var="value")
lipid_ps_overlap$NEG/lipid_ps_overlap$POS

lipid_cl_pos<-lipid_data.f.m[lipid_data.f.m$variable %in% "WT" & lipid_data.f.m$Ontology %in% "CL" & lipid_data.f.m$detect_mode %in% "POS",]
lipid_cl_neg<-lipid_data.f.m[lipid_data.f.m$variable %in% "WT" & lipid_data.f.m$Ontology %in% "CL" & lipid_data.f.m$detect_mode %in% "NEG",]
lipid_cl_overlap<-lipid_data.f.m[lipid_data.f.m$variable %in% "WT" & lipid_data.f.m$Ontology %in% "CL" & lipid_data.f.m$Manually_confirmed_lipid %in% lipid_cl_neg$Manually_confirmed_lipid & lipid_data.f.m$Manually_confirmed_lipid %in% lipid_cl_pos$Manually_confirmed_lipid,-c("Formula","km_cluster")]
lipid_cl_overlap<-lipid_cl_overlap[lipid_cl_overlap$Manually_confirmed_lipid %!in% "CL 63:2",]
lipid_cl_overlap<-dcast(data=lipid_cl_overlap,formula = Manually_confirmed_lipid+lipid_size+chain_length+chain_unsaturations+norm_status+strain~detect_mode,value.var="value")
lipid_cl_overlap$NEG/lipid_cl_overlap$POS
# Conclusion the difference in similarly annotated lipid between these two datasets is concerning, but definitely indicates they are not comparable
