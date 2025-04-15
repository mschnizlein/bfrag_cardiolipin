# ICP-MS Gad-DOTA
# June 2024 submission
# Two replicates no stress
# one replicate with stress

library(forcats)
library(ggplot2)
# library(reshape2)
library(growthcurver)
library(ggtext)
library(gridtext)
library(glue)
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

setwd("C:/Users/mksch/Dropbox/matt_data/crosson_lab_msu_schnizlein/bfrag_bile_acid_stress/icp-ms/")

# Data pre-processing
cellvol_means<-read.csv("../microscopy/revised_analysis/data_output/cell_volume_mean.csv")

media_atoms_1<-read_xlsx(path="icp_20240424_gad_cls/BH6-033 wt_bhis.xlsx",sheet="media",col_names = TRUE)[1:3,]
colnames(media_atoms_1)<-c("sample",colnames(media_atoms_1[2:20]))
media_atoms.m<-melt.data.table(as.data.table(media_atoms_1),id.vars = "sample")
media_atoms.m.avg<-aggregate(value~variable,data=media_atoms.m,FUN=mean)

read_icp<-function(filepath,sheetID,expID,strain){
  data_interm<-read_xlsx(path=filepath,sheet=sheetID,col_names = TRUE)[1:3,]
  colnames(data_interm)<-c("sample",colnames(data_interm[2:ncol(data_interm)]))
  data_interm$sample<-paste0(expID,data_interm$sample)
  assign(x=paste(strain,sheetID,expID,sep = "_"),value=data_interm,envir = .GlobalEnv)
}

read_icp4<-function(filepath,sheetID,expID,strain){
  data_interm<-read_xlsx(path=filepath,sheet=sheetID,col_names = TRUE)[1:4,]
  colnames(data_interm)<-c("sample",colnames(data_interm[2:ncol(data_interm)]))
  data_interm$sample<-paste0(expID,data_interm$sample)
  assign(x=paste(strain,sheetID,expID,sep = "_"),value=data_interm,envir = .GlobalEnv)
}

# Read in WT data
read_icp("icp_20240424_gad_cls/BH6-033 wt_bhis.xlsx","mp","A","wt")
read_icp("icp_20240424_gad_cls/BH6-033 wt_bhis.xlsx","atoms","A","wt")
read_icp("icp_20240424_gad_cls/BH6-033 wt_bhis.xlsx","fmol","A","wt")

read_icp("icp_20240425_gad_cls/BH6-033 wt_bhis_2nd.xlsx","mp","B","wt")
read_icp("icp_20240425_gad_cls/BH6-033 wt_bhis_2nd.xlsx","atoms","B","wt")
read_icp("icp_20240425_gad_cls/BH6-033 wt_bhis_2nd.xlsx","fmol","B","wt")

read_icp("icp_20240531_ionophore_dc/BH6-033 wt_bhis_3rd.xlsx","mp","D","wt")
read_icp("icp_20240531_ionophore_dc/BH6-033 wt_bhis_3rd.xlsx","atoms","D","wt")
read_icp("icp_20240531_ionophore_dc/BH6-033 wt_bhis_3rd.xlsx","fmol","D","wt")

read_icp4("icp_20230918_gad_cls_vtpK/BH5-157_B-Fragilis_wt_2023-10-19_Summary.xlsx","mp","C","wt")
read_icp4("icp_20230918_gad_cls_vtpK/BH5-157_B-Fragilis_wt_2023-10-19_Summary.xlsx","atoms_cell","C","wt")
read_icp4("icp_20230918_gad_cls_vtpK/BH5-157_B-Fragilis_wt_2023-10-19_Summary.xlsx","amol","C","wt")

wt_list<-list(wt_mp_A,wt_atoms_A,wt_fmol_A,wt_mp_B,wt_atoms_B,wt_fmol_B,wt_mp_D,wt_atoms_D,wt_fmol_D,wt_mp_C,wt_atoms_cell_C,wt_amol_C)

wt_list<-lapply(wt_list[], function(x) melt.data.table(data= as.data.table(x), id.vars = "sample"))
names(wt_list)<-c("mp","atcell","fmol","mp","atcell","fmol","mp","atcell","fmol","mp","atcell","fmol")
for (i in 1:3){
  setnames(wt_list[[i]], "value", paste0("wt_1",names(wt_list[i])))
}
for (i in 4:6){
  setnames(wt_list[[i]], "value", paste0("wt_2",names(wt_list[i])))
}
for (i in 7:9){
  setnames(wt_list[[i]], "value", paste0("wt_3",names(wt_list[i])))
}
for (i in 10:12){
  setnames(wt_list[[i]], "value", paste0("wt_4",names(wt_list[i])))
}

# Read in cA data
read_icp("icp_20240424_gad_cls/BH6-033 clsA_bhis.xlsx","mp","A","cA")
read_icp("icp_20240424_gad_cls/BH6-033 clsA_bhis.xlsx","atoms","A","cA")

read_icp("icp_20240425_gad_cls/BH6-033 clsA_bhis_2nd.xlsx","mp","B","cA")
read_icp("icp_20240425_gad_cls/BH6-033 clsA_bhis_2nd.xlsx","atoms","B","cA")

read_icp4("icp_20230918_gad_cls_vtpK/BH5-157_B-Fragilis_cA_2023-10-19_Summary.xlsx","mp","C","cA")
read_icp4("icp_20230918_gad_cls_vtpK/BH5-157_B-Fragilis_cA_2023-10-19_Summary.xlsx","atoms_cell","C","cA")

cA_list<-list(cA_mp_A,cA_atoms_A,cA_mp_B,cA_atoms_B,cA_mp_C,cA_atoms_cell_C)
cA_list<-lapply(cA_list[], function(x) melt.data.table(data= as.data.table(x), id.vars = "sample"))
names(cA_list)<-c("mp","atcell","mp","atcell","mp","atcell")
for (i in 1:2){
  setnames(cA_list[[i]], "value", paste0("cA_1",names(cA_list[i])))
}
for (i in 3:4){
  setnames(cA_list[[i]], "value", paste0("cA_2",names(cA_list[i])))
}
for (i in 5:6){
  setnames(cA_list[[i]], "value", paste0("cA_3",names(cA_list[i])))
}

# Read in cA_comp data
read_icp("icp_20240424_gad_cls/BH6-033 clsA-comp_bhis.xlsx","mp","A","cA_comp")
read_icp("icp_20240424_gad_cls/BH6-033 clsA-comp_bhis.xlsx","atoms","A","cA_comp")

read_icp("icp_20240425_gad_cls/BH6-033 clsA-comp_bhis_2nd.xlsx","mp","B","cA_comp")
read_icp("icp_20240425_gad_cls/BH6-033 clsA-comp_bhis_2nd.xlsx","atoms","B","cA_comp")

cA_comp_list<-list(cA_comp_mp_A,cA_comp_atoms_A,cA_comp_mp_B,cA_comp_atoms_B)
cA_comp_list<-lapply(cA_comp_list[], function(x) melt.data.table(data= as.data.table(x), id.vars = "sample"))
names(cA_comp_list)<-c("mp","atcell","mp","atcell")
for (i in 1:2){
  setnames(cA_comp_list[[i]], "value", paste0("cA_comp_1",names(cA_comp_list[i])))
}
for (i in 3:4){
  setnames(cA_comp_list[[i]], "value", paste0("cA_comp_2",names(cA_comp_list[i])))
}

# Read in cB data
read_icp("icp_20240424_gad_cls/BH6-033 clsB_bhis.xlsx","mp","A","cB")
read_icp("icp_20240424_gad_cls/BH6-033 clsB_bhis.xlsx","atoms","A","cB")

read_icp("icp_20240425_gad_cls/BH6-033 clsB_bhis_2nd.xlsx","mp","B","cB")
read_icp("icp_20240425_gad_cls/BH6-033 clsB_bhis_2nd.xlsx","atoms","B","cB")

read_icp4("icp_20230918_gad_cls_vtpK/BH5-157_B-Fragilis_cB_2023-10-19_Summary.xlsx","mp","C","cB")
read_icp4("icp_20230918_gad_cls_vtpK/BH5-157_B-Fragilis_cB_2023-10-19_Summary.xlsx","atoms_cell","C","cB")

cB_list<-list(cB_mp_A,cB_atoms_A,cB_mp_B,cB_atoms_B,cB_mp_C,cB_atoms_cell_C)
cB_list<-lapply(cB_list[], function(x) melt.data.table(data= as.data.table(x), id.vars = "sample"))
names(cB_list)<-c("mp","atcell","mp","atcell","mp","atcell")
for (i in 1:2){
  setnames(cB_list[[i]], "value", paste0("cB_1",names(cB_list[i])))
}
for (i in 3:4){
  setnames(cB_list[[i]], "value", paste0("cB_2",names(cB_list[i])))
}
for (i in 5:6){
  setnames(cB_list[[i]], "value", paste0("cB_3",names(cB_list[i])))
}

# Read in cB_comp data
read_icp("icp_20240424_gad_cls/BH6-033 clsB-comp_bhis.xlsx","mp","A","cB_comp")
read_icp("icp_20240424_gad_cls/BH6-033 clsB-comp_bhis.xlsx","atoms","A","cB_comp")

read_icp("icp_20240425_gad_cls/BH6-033 clsB-comp_bhis_2nd.xlsx","mp","B","cB_comp")
read_icp("icp_20240425_gad_cls/BH6-033 clsB-comp_bhis_2nd.xlsx","atoms","B","cB_comp")

cB_comp_list<-list(cB_comp_mp_A,cB_comp_atoms_A,cB_comp_mp_B,cB_comp_atoms_B)
cB_comp_list<-lapply(cB_comp_list[], function(x) melt.data.table(data= as.data.table(x), id.vars = "sample"))
names(cB_comp_list)<-c("mp","atcell","mp","atcell")
for (i in 1:2){
  setnames(cB_comp_list[[i]], "value", paste0("cB_comp_1",names(cB_comp_list[i])))
}
for (i in 3:4){
  setnames(cB_comp_list[[i]], "value", paste0("cB_comp_2",names(cB_comp_list[i])))
}

# Read in cAB data
read_icp("icp_20240424_gad_cls/BH6-033 clsAB_bhis.xlsx","mp","A","cAB")
read_icp("icp_20240424_gad_cls/BH6-033 clsAB_bhis.xlsx","atoms","A","cAB")

read_icp("icp_20240425_gad_cls/BH6-033 clsAB_bhis_2nd.xlsx","mp","B","cAB")
read_icp("icp_20240425_gad_cls/BH6-033 clsAB_bhis_2nd.xlsx","atoms","B","cAB")

read_icp4("icp_20230918_gad_cls_vtpK/BH5-157_B-Fragilis_cAcB_2023-10-19_Summary.xlsx","mp","C","cAB")
read_icp4("icp_20230918_gad_cls_vtpK/BH5-157_B-Fragilis_cAcB_2023-10-19_Summary.xlsx","atoms_cell","C","cAB")

cAB_list<-list(cAB_mp_A,cAB_atoms_A,cAB_mp_B,cAB_atoms_B,cAB_mp_C,cAB_atoms_cell_C)
cAB_list<-lapply(cAB_list[], function(x) melt.data.table(data= as.data.table(x), id.vars = "sample"))
names(cAB_list)<-c("mp","atcell","mp","atcell","mp","atcell")
for (i in 1:2){
  setnames(cAB_list[[i]], "value", paste0("cAB_1",names(cAB_list[i])))
}
for (i in 3:4){
  setnames(cAB_list[[i]], "value", paste0("cAB_2",names(cAB_list[i])))
}
for (i in 5:6){
  setnames(cAB_list[[i]], "value", paste0("cAB_3",names(cAB_list[i])))
}

# Read in treatment data
read_icp("icp_20240531_ionophore_dc/BH6-033 wt_dc.xlsx","mp","D","wt_dc")
read_icp("icp_20240531_ionophore_dc/BH6-033 wt_dc.xlsx","atoms","D","wt_dc")

read_icp("icp_20240531_ionophore_dc/BH6-033 wt_monensin.xlsx","mp","D","wt_mon")
read_icp("icp_20240531_ionophore_dc/BH6-033 wt_monensin.xlsx","atoms","D","wt_mon")

read_icp("icp_20240531_ionophore_dc/BH6-033 wt_nigericin.xlsx","mp","D","wt_nig")
read_icp("icp_20240531_ionophore_dc/BH6-033 wt_nigericin.xlsx","atoms","D","wt_nig")

treat_list<-list(wt_dc_mp_D,wt_dc_atoms_D,wt_mon_mp_D,wt_mon_atoms_D,wt_nig_mp_D,wt_nig_atoms_D)
treat_list<-lapply(treat_list[], function(x) melt.data.table(data= as.data.table(x), id.vars = "sample"))
names(treat_list)<-c("mp","atcell","mp","atcell","mp","atcell")
for (i in 1:2){
  setnames(treat_list[[i]], "value", paste0("wt_dc_",names(treat_list[i])))
}
for (i in 3:4){
  setnames(treat_list[[i]], "value", paste0("wt_mon_",names(treat_list[i])))
}
for (i in 5:6){
  setnames(treat_list[[i]], "value", paste0("wt_nig_",names(treat_list[i])))
}

list.f<-c(wt_list,cA_list,cA_comp_list,cB_list,cB_comp_list,cAB_list,treat_list)
# wt_list[[i]][,2]<-paste0(wt_list[[i]][[2]],"_",wt_list[[i]][[1]])
icpms.f<-Reduce(function(x,y) merge(x, y, by.all="variable", all= TRUE), list.f)
setnames(icpms.f,"variable","metal")

icpms.m<-melt.data.table(icpms.f,id.vars = c("sample","metal"))

icpms.m$metal<-factor(icpms.m$metal,levels=c("Na","Li","K","Mg","Fe","Co","Ni","Cu","Zn","P","S","Ca","V","Cr","Mn","As","Mo","Se","Cd"))

icpms.m<-icpms.m[!is.na(icpms.m$value),]
icpms.m<-icpms.m[icpms.m$value >=0,]

icpms.m$strain<-substrLeft(icpms.m$variable,8) %>% gsub("_","",.) %>% gsub("[0-9].*","",.) %>% gsub("dc.*","",.) %>% gsub("mon.*","",.) %>% gsub("nig.*","",.) %>% factor(.,levels=c("wt","cA","cAcomp","cB","cBcomp","cAB"))
icpms.m$strain2<-factor(icpms.m$strain, levels=c("cAB","wt","cA","cAcomp","cB","cBcomp"))
icpms.m$biorep<-substrRight(icpms.m$sample,1)
icpms.m$exprep<-substrLeft(icpms.m$sample,1)
icpms.m$data_type<-substrRight(as.character(icpms.m$variable),2)
icpms.m[icpms.m$data_type %in% "ll",]$data_type<-"atcell"
icpms.m$treatment<-substrLeft(icpms.m$variable,6) %>% substrRight(.,3) %>% gsub("_","",.)
icpms.m[icpms.m$treatment %!in% c("dc","mon","nig"),]$treatment<-"bhis"

ggplot(icpms.m[icpms.m$variable %in% c("wt_1atcell","wt_2atcell","wt_3atcell")],aes(x=variable,y=value))+geom_point()+
  facet_wrap(~metal)


breaks_strain<-c("wt","cA","cAcomp","cB","cBcomp","cAB","cAB_acomp","cAB_bcomp")
labels_strain<-c("WT",
            "\u0394 *clsA*","\u0394 *clsA*  att::*clsA*",
            "\u0394 *clsB*","\u0394 *clsB* att::*clsB*",
            "\u0394 *clsAclsB*","\u0394 *clsAclsB* att::*clsA*"
            ,"\u0394 *clsAclsB* att::*clsB*")
colors_strain<-c("#000000","#006198","#56B4E9","#0f682c","#299764","#e38200","#7D4564","#CC79A7")

breaks_strain_trim<-c("wt","cA","cAcomp","cB","cBcomp")
labels_strain_trim<-c("WT",
                 "\u0394 *clsA*","\u0394 *clsA*  att::*clsA*",
                 "\u0394 *clsB*","\u0394 *clsB* att::*clsB*")
colors_strain_trim<-c("#000000","#006198","#56B4E9","#0f682c","#299764")

###
icpms.m$value_normvol<-"NA"
for (i in 1:8){
  icpms.m[icpms.m$data_type %in% "atcell" & icpms.m$strain %in% breaks_strain[i],]$value_normvol<-icpms.m[icpms.m$data_type %in% "atcell" & icpms.m$strain %in% breaks_strain[i],]$value/cellvol_means[cellvol_means$strain %in% breaks_strain[i],]$cell_vol
}
icpms.m$value_normvol<-as.numeric(icpms.m$value_normvol)

# K elements (number)/cell -> mols K -> divide by volume to get molarity
# 1 mol = 6.02214076 * 10^23 number
# 1 cubic micron = 1e-12 mL = 1e-15 L
# i<-1
# elements/cellx1E6 -> mol/cell -> mol/cell * cell/cubic micron -> * cubic micron/1E-15 L = M
# i<-1
icpms.m$value_M<-"NA"
for (i in 1:8){
icpms.m[icpms.m$strain %in% breaks_strain[i],]$value_M<-(icpms.m[icpms.m$strain %in% breaks_strain[i],]$value*1E6)*(1/6.022E23)*
  (1/cellvol_means[cellvol_means$strain %in% breaks_strain[i],]$cell_vol)*(1/1E-15)
}
icpms.m$value_M<-as.numeric(icpms.m$value_M)
icpms.m$value_mM<-icpms.m$value_M*1E3
icpms.m$value_uM<-icpms.m$value_M*1E6
icpms.m$value_nM<-icpms.m$value_M*1E9

# atoms -> *1 mol/6.022 atoms -> *1/vol (mL to L)
# media_atoms.m.avg
breaks_metals<-unique(media_atoms.m.avg$variable)
icpms.m<-icpms.m[icpms.m$metal %!in% "Li",]

icpms.m$value_mM_deltamedia<-"NA"
for (i in 1:19){
  icpms.m[icpms.m$metal %in% breaks_metals[i],]$value_mM_deltamedia<-icpms.m[icpms.m$metal %in% breaks_metals[i],]$value_mM-media_atoms.m.avg[media_atoms.m.avg$variable %in% breaks_metals[i],]$value
}
icpms.m$value_mM_deltamedia<-as.numeric(icpms.m$value_mM_deltamedia)

####
# Relative abundance per sample
icpms.m$var_bio<-paste(icpms.m$variable,icpms.m$biorep,sep="_")
breaks_variables<-unique(as.character(icpms.m[icpms.m$data_type %in% "atcell",]$var_bio))

icpms.m$relabund_mM<-0
for (i in 1:length(breaks_variables)){
  icpms.m[icpms.m$data_type %in% "atcell" & icpms.m$var_bio %in% breaks_variables[i],]$relabund_mM<-(icpms.m[icpms.m$data_type %in% "atcell" & icpms.m$var_bio %in% breaks_variables[i],]$value_mM/sum(icpms.m[icpms.m$data_type %in% "atcell" & icpms.m$var_bio %in% breaks_variables[i],]$value_mM))*100
}

icpms.m$metal_1perc<-icpms.m$metal
icpms.m[icpms.m$relabund_mM <=0.75,]$metal_1perc<-"Other"
icpms.m_sub<-icpms.m[icpms.m$data_type %in% "atcell",c("strain","treatment","data_type","metal","metal_1perc","relabund_mM")]

# write.csv(file = "cls_manuscript_icpms_dataset.csv",x=icpms.m, row.names = FALSE)
# icpms.m<-read.csv(file="cls_manuscript_icpms_dataset.csv",header=TRUE)

icpms.m_avgrel<-NULL
correction_cb<-icpms.m_sub[icpms.m_sub$metal %in% "Na" & icpms.m_sub$strain %in% "cB" & icpms.m_sub$relabund_mM >=5,] # relabund condition just to subset 3 rows
correction_cb$relabund_mM<-0
icpms.m_avgrel<-aggregate(relabund_mM~strain+treatment+data_type+metal+metal_1perc,data=rbind(icpms.m_sub,correction_cb),FUN=mean)
# Aggregate has trouble with one set of values not being the same as the rest (i.e., Na); several of these had negative values due to the Gd-Dota normalization and were excluded above. If I sub in 3 extra rows with value of 0 it fixes whatever overcorrection may be happening due to those missing values

#####
summary(icpms.m[icpms.m$data_type %in% "atcell" & icpms.m$strain %in% "wt" & icpms.m$treatment %in% "bhis" & icpms.m$metal %in% "K",]$value_mM) # mean = 434.1
summary(icpms.m[icpms.m$data_type %in% "atcell" & icpms.m$strain %in% "wt" & icpms.m$treatment %in% "bhis" & icpms.m$metal %in% "Na",]$value_mM) # mean = 62.76
summary(icpms.m[icpms.m$data_type %in% "atcell" & icpms.m$strain %in% "wt" & icpms.m$treatment %in% "bhis" & icpms.m$metal %in% "P",]$value_mM) # mean = 544.5
summary(icpms.m[icpms.m$data_type %in% "atcell" & icpms.m$strain %in% "wt" & icpms.m$treatment %in% "bhis" & icpms.m$metal %in% "S",]$value_mM) # mean = 97.8
summary(icpms.m[icpms.m$data_type %in% "atcell" & icpms.m$strain %in% "wt" & icpms.m$treatment %in% "bhis" & icpms.m$metal %in% "Mg",]$value_mM) # mean = 68.97
summary(icpms.m[icpms.m$data_type %in% "atcell" & icpms.m$strain %in% "wt" & icpms.m$treatment %in% "bhis" & icpms.m$metal %in% "Fe",]$value_mM) # mean = 4.149
summary(icpms.m[icpms.m$data_type %in% "atcell" & icpms.m$strain %in% "wt" & icpms.m$treatment %in% "bhis" & icpms.m$metal %in% "Ca",]$value_mM) # mean = 4.138

summary(icpms.m[icpms.m$data_type %in% "atcell" & icpms.m$strain %in% "wt" & icpms.m$exprep %in% "D" & icpms.m$treatment %in% "bhis" & icpms.m$metal %in% "Na",]$value_mM) - summary(icpms.m[icpms.m$data_type %in% "atcell" & icpms.m$strain %in% "wt" & icpms.m$exprep %in% "D" & icpms.m$treatment %in% "mon" & icpms.m$metal %in% "Na",]$value_mM) # mean = 32.51

# Make Na, P, S, K, and Fe plots (note they are times 10^6)

allstrains_k_atcell.p<-ggplot(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("K") & icpms.m$data_type %in% "atcell",], aes(x=strain, y=value,fill=strain))+geom_boxplot(outliers = FALSE)+
  scale_y_continuous(limits=c(100,400))+
  scale_x_discrete(breaks=breaks_strain,labels=labels_strain)+
  ylab("K/cell x10^6")+
  xlab("")+
  scale_fill_manual(values=colors_strain,breaks=breaks_strain,labels=labels_strain)+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(allstrains_k_atcell.p)

allstrains_k_mM.p<-ggplot(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("K") & icpms.m$data_type %in% "atcell",], aes(x=strain, y=value_mM,fill=strain))+geom_boxplot(outliers = FALSE)+geom_point()+
  scale_y_continuous(name="K (mM)",limits=c(200,700))+
  scale_x_discrete(name="",breaks=breaks_strain,labels=labels_strain)+
  scale_fill_manual(values=colors_strain,breaks=breaks_strain,labels=labels_strain)+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(allstrains_k_mM.p)

k.atcell.lm<-lm(value_mM_deltamedia~strain,icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("K") & icpms.m$data_type %in% "atcell",])
summary(k.atcell.lm)
# cA and cB sig (0.03066 and 0.00167)
k.atcell.lm<-lm(value_mM_deltamedia~strain2,icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("K") & icpms.m$data_type %in% "atcell",])
summary(k.atcell.lm)

allstrains_k_mp.p<-ggplot(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("K") & icpms.m$data_type %in% "mp",], aes(x=strain, y=value,fill=strain))+geom_boxplot()+
  scale_y_continuous(limits=c(0.7,0.9))+
  scale_x_discrete(breaks=breaks_strain,labels=labels_strain)+
  ylab("K/Phosphorus")+
  xlab("")+
  scale_fill_manual(values=colors_strain,breaks=breaks_strain,labels=labels_strain)+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(allstrains_k_mp.p)

k.mp.lm<-lm(value~strain,icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("K") & icpms.m$data_type %in% "mp",])
summary(k.mp.lm)
k.mp2.lm<-lm(value~strain2,icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("K") & icpms.m$data_type %in% "mp",])
summary(k.mp2.lm)

ggsave(filename="figures/revised_202407/allstrains_k_atcell.pdf",plot=allstrains_k_atcell.p,width=3,height=3,units="in",device=cairo_pdf)
ggsave(filename="figures/revised_202407/allstrains_k_mM.pdf",plot=allstrains_k_mM.p,width=3,height=3,units="in",device=cairo_pdf)
ggsave(filename="figures/revised_202407/allstrains_k_mp.pdf",plot=allstrains_k_mp.p,width=3,height=3,units="in",device=cairo_pdf)



allstrains_na_atcell.p<-ggplot(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("Na") & icpms.m$data_type %in% "atcell",], aes(x=strain, y=value,fill=strain))+geom_boxplot()+
  scale_y_continuous(limits=c(0,70))+
  scale_x_discrete(breaks=breaks_strain,labels=labels_strain)+
  ylab("Na/cell x10^6")+
  xlab("")+
  scale_fill_manual(values=colors_strain,breaks=breaks_strain,labels=labels_strain)+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(allstrains_na_atcell.p)

allstrains_na_mM.p<-ggplot(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("Na") & icpms.m$data_type %in% "atcell",], aes(x=strain, y=value_mM,fill=strain))+geom_boxplot()+
  scale_y_continuous(name="Na (mM)",limits=c(0,150))+
  scale_x_discrete(name="",breaks=breaks_strain,labels=labels_strain)+
  scale_fill_manual(values=colors_strain,breaks=breaks_strain,labels=labels_strain)+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(allstrains_na_mM.p)

allstrains_na_mp.p<-ggplot(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("Na") & icpms.m$data_type %in% "mp",], aes(x=strain, y=value,fill=strain))+geom_boxplot()+
  scale_y_continuous(limits=c(0,0.25))+
  scale_x_discrete(breaks=breaks_strain,labels=labels_strain)+
  ylab("Na/Phosphorus")+
  xlab("")+
  scale_fill_manual(values=colors_strain,breaks=breaks_strain,labels=labels_strain)+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(allstrains_na_mp.p)

na.atcell.lm<-lm(value_mM_deltamedia~strain,icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("Na") & icpms.m$data_type %in% "atcell",])
summary(na.atcell.lm)
# none significant

na.mp.lm<-lm(value~strain,icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("Na") & icpms.m$data_type %in% "mp",])
summary(na.mp.lm)
# none significant

ggsave(filename="figures/revised_202407/allstrains_na_atcell.pdf",plot=allstrains_na_atcell.p,width=3,height=3,units="in",device=cairo_pdf)
ggsave(filename="figures/revised_202407/allstrains_na_mM.pdf",plot=allstrains_na_mM.p,width=3,height=3,units="in",device=cairo_pdf)
ggsave(filename="figures/revised_202407/allstrains_na_mp.pdf",plot=allstrains_na_mp.p,width=3,height=3,units="in",device=cairo_pdf)

allstrains_p_mM.p<-ggplot(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("P") & icpms.m$data_type %in% "atcell",], aes(x=strain, y=value_mM,fill=strain))+geom_boxplot()+
  scale_y_continuous(name="P (mM)",limits=c(300,900))+
  scale_x_discrete(name="",breaks=breaks_strain,labels=labels_strain)+
  scale_fill_manual(values=colors_strain,breaks=breaks_strain,labels=labels_strain)+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(allstrains_p_mM.p)

p.atcell.lm<-lm(value_mM~strain,icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("P") & icpms.m$data_type %in% "atcell",])
summary(p.atcell.lm)
# cB sig (0.00233)
p.atcell2.lm<-lm(value_mM~strain2,icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("P") & icpms.m$data_type %in% "atcell",])
summary(p.atcell2.lm)
# cA and cB sig (0.0031 and 4.7E-5)

ggsave(filename="figures/revised_202407/allstrains_p_mM.pdf",plot=allstrains_p_mM.p,width=3,height=3,units="in",device=cairo_pdf)

treat_p_mM.p<-ggplot(icpms.m[icpms.m$strain %in% "wt" & icpms.m$metal %in% "P" & icpms.m$data_type %in% "atcell"& icpms.m$exprep %in% "D",], aes(x=treatment, y=value_mM,fill=treatment))+geom_boxplot(fill="black")+
  scale_y_continuous(name="P (mM)",limits=c(400,800))+
  scale_x_discrete(breaks=c("bhis","dc","mon","nig"),labels=c("BHIS","0.01% Deoxycholate","0.8 ug/mL Monensin A","0.1 ug/mL Nigericin"))+
  ylab("P/cell x10^6")+
  xlab("")+
  as_md_theme(my_theme2)+theme(legend.position = "none")
plot(treat_p_mM.p)

treat.p.atcell.lm<-lm(value_mM~treatment,icpms.m[icpms.m$strain %in% "wt" & icpms.m$metal %in% "P" & icpms.m$data_type %in% "atcell",])
summary(treat.p.atcell.lm)
# Dc sig (2.96E-5)

treat_k_mM.p<-ggplot(icpms.m[icpms.m$strain %in% "wt" & icpms.m$metal %in% "K" & icpms.m$data_type %in% "atcell" & icpms.m$exprep %in% "D",], aes(x=treatment, y=value_mM,fill=treatment))+geom_boxplot(fill="black")+
  scale_y_continuous(name="K (mM)",limits=c(300,600))+
  scale_x_discrete(name="",breaks=c("bhis","dc","mon","nig"),labels=c("BHIS","0.01% Deoxycholate","0.8 ug/mL Monensin A","0.1 ug/mL Nigericin"))+
  as_md_theme(my_theme2)+theme(legend.position = "none")
plot(treat_k_mM.p)

treat_k_mM_delta.p<-ggplot(icpms.m[icpms.m$strain %in% "wt" & icpms.m$metal %in% "K" & icpms.m$data_type %in% "atcell" & icpms.m$exprep %in% "D",], aes(x=treatment, y=value_mM_deltamedia,fill=treatment))+geom_boxplot(fill="black")+
  scale_y_continuous(name="Intracellular - Extracellular Element (\u0394\U00B5M)",limits=c(300,600))+
  scale_x_discrete(name="",breaks=c("bhis","dc","mon","nig"),labels=c("BHIS","0.01% Deoxycholate","0.8 ug/mL Monensin A","0.1 ug/mL Nigericin"))+
  as_md_theme(my_theme2)+theme(legend.position = "none")
plot(treat_k_mM_delta.p)

treat.k.atcell.lm<-lm(value_mM_deltamedia~treatment,icpms.m[icpms.m$strain %in% "wt" & icpms.m$metal %in% "K" & icpms.m$data_type %in% "atcell" & icpms.m$exprep %in% "D",])
summary(treat.k.atcell.lm)
# Dc sig (0.00683)

treat_k_mp.p<-ggplot(icpms.m[icpms.m$strain %in% "wt" & icpms.m$metal %in% "K" & icpms.m$data_type %in% "mp",], aes(x=treatment, y=value,fill=treatment))+geom_boxplot(fill="black")+
  scale_y_continuous(limits=c(0.4,0.9))+
  scale_x_discrete(breaks=c("bhis","dc","mon","nig"),labels=c("BHIS","0.01% Deoxycholate","0.8 ug/mL Monensin A","0.1 ug/mL Nigericin"))+
  ylab("K/phosphorus")+
  xlab("")+
  as_md_theme(my_theme2)+theme(legend.position = "none")
plot(treat_k_mp.p)

treat_na_mM.p<-ggplot(icpms.m[icpms.m$strain %in% "wt" & icpms.m$metal %in% "Na" & icpms.m$data_type %in% "atcell"& icpms.m$exprep %in% "D",], aes(x=treatment, y=value_mM,fill=treatment))+geom_boxplot(fill="black")+
  scale_y_continuous(name="Na (mM)",limits=c(0,150))+
  scale_x_discrete(breaks=c("bhis","dc","mon","nig"),labels=c("BHIS","0.01% Deoxycholate","0.8 ug/mL Monensin A","0.1 ug/mL Nigericin"))+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(treat_na_mM.p)

treat_na_mM_delta.p<-ggplot(icpms.m[icpms.m$strain %in% "wt" & icpms.m$metal %in% "Na" & icpms.m$data_type %in% "atcell" & icpms.m$exprep %in% "D",], aes(x=treatment, y=value_mM_deltamedia,fill=treatment))+geom_boxplot(fill="black")+
  scale_y_continuous(name="Intracellular - Extracellular Element (\u0394\U00B5M)",limits=c(-100,0))+
  scale_x_discrete(breaks=c("bhis","dc","mon","nig"),labels=c("BHIS","0.01% Deoxycholate","0.8 ug/mL Monensin A","0.1 ug/mL Nigericin"))+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(treat_na_mM_delta.p)

treat.na.atcell.lm<-lm(value_mM_deltamedia~treatment,icpms.m[icpms.m$strain %in% "wt" & icpms.m$metal %in% "Na" & icpms.m$data_type %in% "atcell" & icpms.m$exprep %in% "D",])
summary(treat.na.atcell.lm)
# Dc and Mon sig (0.000535 and 0.0197)

treat_na_mp.p<-ggplot(icpms.m[icpms.m$strain %in% "wt" & icpms.m$metal %in% "Na" & icpms.m$data_type %in% "mp",], aes(x=treatment, y=value,fill=treatment))+geom_boxplot(fill="black")+
  scale_y_continuous(limits=c(0.05,0.2))+
  scale_x_discrete(breaks=c("bhis","dc","mon","nig"),labels=c("BHIS","0.01% Deoxycholate","0.8 ug/mL Monensin A","0.1 ug/mL Nigericin"))+
  ylab("Na/Phosphorus")+
  xlab("")+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(treat_na_mp.p)

ggsave(filename="figures/revised_202407/treat_p_mM.pdf",plot=treat_p_mM.p,width=3,height=3,units="in",device=cairo_pdf)
ggsave(filename="figures/revised_202407/treat_k_mM.pdf",plot=treat_k_mM.p,width=3,height=3,units="in",device=cairo_pdf)
ggsave(filename="figures/revised_202407/treat_k_mp.pdf",plot=treat_k_mp.p,width=3,height=3,units="in",device=cairo_pdf)
ggsave(filename="figures/revised_202407/treat_na_mM.pdf",plot=treat_na_mM.p,width=3,height=3,units="in",device=cairo_pdf)
ggsave(filename="figures/revised_202407/treat_na_mp.pdf",plot=treat_na_mp.p,width=3,height=3,units="in",device=cairo_pdf)
ggsave(filename="figures/revised_202407/treat_na_mM_delta.pdf",plot=treat_na_mM_delta.p,width=3,height=3,units="in",device=cairo_pdf)
ggsave(filename="figures/revised_202407/treat_k_mM_delta.pdf",plot=treat_k_mM_delta.p,width=3,height=3,units="in",device=cairo_pdf)

treat_allelements_high_1.p<-ggplot(icpms.m[icpms.m$strain %in% "wt" & icpms.m$data_type %in% "atcell"  & icpms.m$metal %in% c("Mg","S")& icpms.m$exprep %in% "D",], aes(x=treatment, y=value_mM,fill=treatment))+geom_boxplot(fill="black")+
  scale_y_continuous(name="Element (mM)",limits=c(60,140))+
  scale_x_discrete(name="",breaks=c("bhis","dc","mon","nig"),labels=c("BHIS","0.01% Deoxycholate","0.8 ug/mL Monensin A","0.1 ug/mL Nigericin"))+
  as_md_theme(my_theme2)+theme(legend.position = "none")+facet_wrap(~metal)
plot(treat_allelements_high_1.p)

treat_allelements_high_1_mp.p<-ggplot(icpms.m[icpms.m$strain %in% "wt" & icpms.m$data_type %in% "mp"  & icpms.m$metal %in% c("Mg","S") & icpms.m$exprep %in% "D",], aes(x=treatment, y=value,fill=treatment))+geom_boxplot(fill="black")+
  scale_y_continuous(name="Element/Phosphorus")+
  scale_x_discrete(name="",breaks=c("bhis","dc","mon","nig"),labels=c("BHIS","0.01% Deoxycholate","0.8 ug/mL Monensin A","0.1 ug/mL Nigericin"))+
  as_md_theme(my_theme2)+theme(legend.position = "none")+facet_wrap(~metal)
plot(treat_allelements_high_1_mp.p)

treat.mg.atcell.lm<-lm(value_mM_deltamedia~treatment,icpms.m[icpms.m$strain %in% "wt" & icpms.m$metal %in% "Mg" & icpms.m$data_type %in% "atcell" & icpms.m$exprep %in% "D",])
summary(treat.mg.atcell.lm)
# Dc sig (0.0314)
treat.mg.mp.lm<-lm(value~treatment,icpms.m[icpms.m$strain %in% "wt" & icpms.m$metal %in% "Na" & icpms.m$data_type %in% "mp" & icpms.m$exprep %in% "D",])
summary(treat.mg.mp.lm)
# Dc and Mon sig (0.0024 and 0.0123)

treat_allelements_high_mM_2.p<-ggplot(icpms.m[icpms.m$strain %in% "wt" & icpms.m$data_type %in% "atcell"  & icpms.m$metal %in% c("Ca","Fe")& icpms.m$exprep %in% "D",], aes(x=treatment, y=value_mM,fill=treatment))+geom_boxplot(fill="black")+
  scale_y_continuous(name="Element (mM)")+
  scale_x_discrete(name="",breaks=c("bhis","dc","mon","nig"),labels=c("BHIS","0.01% Deoxycholate","0.8 ug/mL Monensin A","0.1 ug/mL Nigericin"))+
  as_md_theme(my_theme2)+theme(legend.position = "none")+facet_wrap(~metal)
plot(treat_allelements_high_mM_2.p)

treat_allelements_high_mM_2_mp.p<-ggplot(icpms.m[icpms.m$strain %in% "wt" & icpms.m$data_type %in% "mp"  & icpms.m$metal %in% c("Ca","Fe")& icpms.m$exprep %in% "D",], aes(x=treatment, y=value,fill=treatment))+geom_boxplot(fill="black")+
  scale_y_continuous(name="Element/Phosphorus")+
  scale_x_discrete(name="",breaks=c("bhis","dc","mon","nig"),labels=c("BHIS","0.01% Deoxycholate","0.8 ug/mL Monensin A","0.1 ug/mL Nigericin"))+
  as_md_theme(my_theme2)+theme(legend.position = "none")+facet_wrap(~metal)
plot(treat_allelements_high_mM_2_mp.p)

treat_allelements_high_uM_3.p<-ggplot(icpms.m[icpms.m$strain %in% "wt" & icpms.m$data_type %in% "atcell"  & icpms.m$metal %in% c("Cu","Zn")& icpms.m$exprep %in% "D",], aes(x=treatment, y=value_uM,fill=treatment))+geom_boxplot(fill="black")+
  scale_y_continuous(name="Element (\U00B5M)")+
  scale_x_discrete(breaks=c("bhis","dc","mon","nig"),labels=c("BHIS","0.01% Deoxycholate","0.8 ug/mL Monensin A","0.1 ug/mL Nigericin"))+
  ylab("Element (mM)")+
  xlab("")+
  as_md_theme(my_theme2)+theme(legend.position = "none")+facet_wrap(~metal)
plot(treat_allelements_high_uM_3.p)

treat_allelements_high_uM_3_mp.p<-ggplot(icpms.m[icpms.m$strain %in% "wt" & icpms.m$data_type %in% "mp"  & icpms.m$metal %in% c("Cu","Zn")& icpms.m$exprep %in% "D",], aes(x=treatment, y=value_uM,fill=treatment))+geom_boxplot(fill="black")+
  scale_y_continuous(name="Element (\U00B5M)")+
  scale_x_discrete(breaks=c("bhis","dc","mon","nig"),labels=c("BHIS","0.01% Deoxycholate","0.8 ug/mL Monensin A","0.1 ug/mL Nigericin"))+
  ylab("Element/Phosphorus")+
  xlab("")+
  as_md_theme(my_theme2)+theme(legend.position = "none")+facet_wrap(~metal)
plot(treat_allelements_high_uM_3_mp.p)

treat_allelements_high_uM_4.p<-ggplot(icpms.m[icpms.m$strain %in% "wt" & icpms.m$data_type %in% "atcell"  & icpms.m$metal %in% c("Mn"),], aes(x=treatment, y=value_uM,fill=treatment))+geom_boxplot(fill="black")+
  scale_y_continuous(name="Element (\U00B5M)")+
  scale_x_discrete(name="",breaks=c("bhis","dc","mon","nig"),labels=c("BHIS","0.01% Deoxycholate","0.8 ug/mL Monensin A","0.1 ug/mL Nigericin"))+
  as_md_theme(my_theme2)+theme(legend.position = "none")+facet_wrap(~metal)
plot(treat_allelements_high_uM_4.p)

treat_allelements_high_uM_5.p<-ggplot(icpms.m[icpms.m$strain %in% "wt" & icpms.m$data_type %in% "atcell"  & icpms.m$metal %in% c("Co"),], aes(x=treatment, y=value_uM,fill=treatment))+geom_boxplot(fill="black")+
  scale_y_continuous(name="Element (\U00B5M)")+
  scale_x_discrete(name="",breaks=c("bhis","dc","mon","nig"),labels=c("BHIS","0.01% Deoxycholate","0.8 ug/mL Monensin A","0.1 ug/mL Nigericin"))+
  as_md_theme(my_theme2)+theme(legend.position = "none")+facet_wrap(~metal)
plot(treat_allelements_high_uM_5.p)

treat_allelements_high_uM_6.p<-ggplot(icpms.m[icpms.m$strain %in% "wt" & icpms.m$data_type %in% "atcell"  & icpms.m$metal %in% c("Cr","As","Mo"),], aes(x=treatment, y=value_uM,fill=treatment))+geom_boxplot(fill="black")+
  scale_y_continuous(name="Element (\U00B5M)")+
  scale_x_discrete(name="",breaks=c("bhis","dc","mon","nig"),labels=c("BHIS","0.01% Deoxycholate","0.8 ug/mL Monensin A","0.1 ug/mL Nigericin"))+
  as_md_theme(my_theme2)+theme(legend.position = "none")+facet_wrap(~metal)
plot(treat_allelements_high_uM_6.p)

treat_allelements_high_uM_7.p<-ggplot(icpms.m[icpms.m$strain %in% "wt" & icpms.m$data_type %in% "atcell"  & icpms.m$metal %in% c("V","Ni","Se","Cd"),], aes(x=treatment, y=value_uM,fill=treatment))+geom_boxplot(fill="black")+
  scale_y_continuous(name="Element (\U00B5M)")+
  scale_x_discrete(name="",breaks=c("bhis","dc","mon","nig"),labels=c("BHIS","0.01% Deoxycholate","0.8 ug/mL Monensin A","0.1 ug/mL Nigericin"))+
  as_md_theme(my_theme2)+theme(legend.position = "none")+facet_wrap(~metal,nrow=1)
plot(treat_allelements_high_uM_7.p)

ggsave(filename="figures/revised_202407/treat_allelements_1.pdf",plot=treat_allelements_high_1.p,width=3,height=3,units="in",device=cairo_pdf)
ggsave(filename="figures/revised_202407/treat_allelements_2.pdf",plot=treat_allelements_high_mM_2.p,width=3,height=3,units="in",device=cairo_pdf)
ggsave(filename="figures/revised_202407/treat_allelements_3.pdf",plot=treat_allelements_high_uM_3.p,width=3,height=3,units="in",device=cairo_pdf)
ggsave(filename="figures/revised_202407/treat_allelements_4.pdf",plot=treat_allelements_high_uM_4.p,width=3,height=3,units="in",device=cairo_pdf)
ggsave(filename="figures/revised_202407/treat_allelements_5.pdf",plot=treat_allelements_high_uM_5.p,width=3,height=3,units="in",device=cairo_pdf)
ggsave(filename="figures/revised_202407/treat_allelements_6.pdf",plot=treat_allelements_high_uM_6.p,width=3,height=3,units="in",device=cairo_pdf)
ggsave(filename="figures/revised_202407/treat_allelements_7.pdf",plot=treat_allelements_high_uM_7.p,width=3,height=3,units="in",device=cairo_pdf)


###

wt_allele_log10_mM.p<-ggplot(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$strain %in% "wt" & icpms.m$data_type %in% "atcell",], aes(x=reorder(metal,-value_mM), y=value_mM))+geom_boxplot(fill="black",color="grey")+
  scale_y_log10(name="Element (mM)",breaks = trans_breaks("log10", function(x) 10^x,n=4), labels = trans_format("log10", label_math(10^.x)),limits=c(1E-4,1E4))+
  annotation_logticks(base=10,sides = "l",scaled=TRUE,short = unit(0.05, "cm"),mid = unit(0.1, "cm"),long = unit(0.2, "cm"),)+
  scale_x_discrete(name="")+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(wt_allele_log10_mM.p)

ggsave(filename="figures/revised_202407/wt_allele_log10_mM.pdf",plot=wt_allele_log10_mM.p,width=5,height=2,units="in",device=cairo_pdf)

allstrains_ele_mM_1.p<-ggplot(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("S","Mg") & icpms.m$data_type %in% "atcell",], aes(x=strain, y=value_mM,fill=strain))+geom_boxplot()+
  scale_y_continuous(name="Element (mM)")+
  scale_x_discrete(name="",breaks=breaks_strain,labels=labels_strain)+
  scale_fill_manual(values=colors_strain,breaks=breaks_strain,labels=labels_strain)+
  as_md_theme(my_theme)+theme(legend.position = "none")+facet_wrap(~metal)
plot(allstrains_ele_mM_1.p)

allstrains_ele_mp_1.p<-ggplot(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("S","Mg") & icpms.m$data_type %in% "mp",], aes(x=strain, y=value,fill=strain))+geom_boxplot()+
  scale_y_continuous(name="Element/Phosphorus",limits=c(0.1,0.25))+
  scale_x_discrete(name="",breaks=breaks_strain,labels=labels_strain)+
  scale_fill_manual(values=colors_strain,breaks=breaks_strain,labels=labels_strain)+
  as_md_theme(my_theme)+theme(legend.position = "none")+facet_wrap(~metal)
plot(allstrains_ele_mp_1.p)

strain.s.mp.lm<-lm(value~strain,icpms.m[icpms.m$treatment %in% "bhis" & icpms.m$metal %in% "S" & icpms.m$data_type %in% "mp",])
summary(strain.s.mp.lm)
# none
strain.mg.mp.lm<-lm(value~strain,icpms.m[icpms.m$treatment %in% "bhis" & icpms.m$metal %in% "Mg" & icpms.m$data_type %in% "mp",])
summary(strain.mg.mp.lm)
# cBcomp (0.0331)

strain.s.mp2.lm<-lm(value~strain2,icpms.m[icpms.m$treatment %in% "bhis" & icpms.m$metal %in% "S" & icpms.m$data_type %in% "mp",])
summary(strain.s.mp.lm)
# none
strain.mg.mp.lm<-lm(value~strain2,icpms.m[icpms.m$treatment %in% "bhis" & icpms.m$metal %in% "Mg" & icpms.m$data_type %in% "mp",])
summary(strain.mg.mp.lm)
# cAcomp and cBcomp (0.00806 and 0.0015)

allstrains_ele_uM_2.p<-ggplot(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("Ca","Fe") & icpms.m$data_type %in% "atcell",], aes(x=strain, y=value_uM,fill=strain))+geom_boxplot()+
  scale_y_log10(name="Element (\U00B5M)")+
  scale_x_discrete(name="",breaks=breaks_strain,labels=labels_strain)+
  scale_fill_manual(values=colors_strain,breaks=breaks_strain,labels=labels_strain)+
  as_md_theme(my_theme)+theme(legend.position = "none")+facet_wrap(~metal)
plot(allstrains_ele_uM_2.p)

allstrains_ele_mp_2.p<-ggplot(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("Ca","Fe") & icpms.m$data_type %in% "mp",], aes(x=strain, y=value,fill=strain))+geom_boxplot()+
  scale_y_continuous(name="Element/Phosphorus",limits=c(0.001,0.015))+
  scale_x_discrete(name="",breaks=breaks_strain,labels=labels_strain)+
  scale_fill_manual(values=colors_strain,breaks=breaks_strain,labels=labels_strain)+
  as_md_theme(my_theme)+theme(legend.position = "none")+facet_wrap(~metal)
plot(allstrains_ele_mp_2.p)

strain.fe.mp.lm<-lm(value~strain,icpms.m[icpms.m$treatment %in% "bhis" & icpms.m$metal %in% "Fe" & icpms.m$data_type %in% "mp",])
summary(strain.fe.mp.lm)
# cB, cAB (0.00497, 0.04462)
strain.ca.mp.lm<-lm(value~strain,icpms.m[icpms.m$treatment %in% "bhis" & icpms.m$metal %in% "Ca" & icpms.m$data_type %in% "mp",])
summary(strain.ca.mp.lm)
# none)

allstrains_ele_uM_3.p<-ggplot(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("Cu","Zn") & icpms.m$data_type %in% "atcell",], aes(x=strain, y=value_uM,fill=strain))+geom_boxplot()+
  scale_y_continuous(name="Element (\U00B5M)")+
  scale_x_discrete(name="",breaks=breaks_strain,labels=labels_strain)+
  scale_fill_manual(values=colors_strain,breaks=breaks_strain,labels=labels_strain)+
  as_md_theme(my_theme)+theme(legend.position = "none")+facet_wrap(~metal)
plot(allstrains_ele_uM_3.p)

allstrains_ele_mp_3.p<-ggplot(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("Cu","Zn") & icpms.m$data_type %in% "mp",], aes(x=strain, y=value,fill=strain))+geom_boxplot()+
  scale_y_continuous(name="Element/Phosphours",limits=c(0.0001,0.0025))+
  scale_x_discrete(name="",breaks=breaks_strain,labels=labels_strain)+
  scale_fill_manual(values=colors_strain,breaks=breaks_strain,labels=labels_strain)+
  as_md_theme(my_theme)+theme(legend.position = "none")+facet_wrap(~metal)
plot(allstrains_ele_mp_3.p)

strain.cu.mp.lm<-lm(value~strain,icpms.m[icpms.m$treatment %in% "bhis" & icpms.m$metal %in% "Cu" & icpms.m$data_type %in% "mp",])
summary(strain.cu.mp.lm)
# none
strain.zn.mp.lm<-lm(value~strain,icpms.m[icpms.m$treatment %in% "bhis" & icpms.m$metal %in% "Zn" & icpms.m$data_type %in% "mp",])
summary(strain.zn.mp.lm)
# none)

allstrains_ele_uM_4.p<-ggplot(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("Mn") & icpms.m$data_type %in% "atcell",], aes(x=strain, y=value_uM,fill=strain))+geom_boxplot()+
  scale_y_continuous(name="Element (\U00B5M)")+
  scale_x_discrete(name="",breaks=breaks_strain,labels=labels_strain)+
  scale_fill_manual(values=colors_strain,breaks=breaks_strain,labels=labels_strain)+
  as_md_theme(my_theme)+theme(legend.position = "none")+facet_wrap(~metal)
plot(allstrains_ele_uM_4.p)

allstrains_ele_mp_4.p<-ggplot(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("Mn") & icpms.m$data_type %in% "mp",], aes(x=strain, y=value*1000,fill=strain))+geom_boxplot()+
  scale_y_continuous(name="Element/1000 Phosphorus",limits=c(0.05,0.2))+
  scale_x_discrete(name="",breaks=breaks_strain,labels=labels_strain)+
  scale_fill_manual(values=colors_strain,breaks=breaks_strain,labels=labels_strain)+
  as_md_theme(my_theme)+theme(legend.position = "none")+facet_wrap(~metal)
plot(allstrains_ele_mp_4.p)

allstrains_ele_uM_5.p<-ggplot(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("Co") & icpms.m$data_type %in% "atcell",], aes(x=strain, y=value_uM,fill=strain))+geom_boxplot()+
  scale_y_continuous(name="Element (\U00B5M)")+
  scale_x_discrete(name="",breaks=breaks_strain,labels=labels_strain)+
  scale_fill_manual(values=colors_strain,breaks=breaks_strain,labels=labels_strain)+
  as_md_theme(my_theme)+theme(legend.position = "none")+facet_wrap(~metal)
plot(allstrains_ele_uM_5.p)

allstrains_ele_mp_5.p<-ggplot(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("Co") & icpms.m$data_type %in% "mp",], aes(x=strain, y=value * 1000,fill=strain))+geom_boxplot()+
  scale_y_continuous(name="Element/1000 Phosphorus",limits=c(0.01,0.1))+
  scale_x_discrete(name="",breaks=breaks_strain,labels=labels_strain)+
  scale_fill_manual(values=colors_strain,breaks=breaks_strain,labels=labels_strain)+
  as_md_theme(my_theme)+theme(legend.position = "none")+facet_wrap(~metal)
plot(allstrains_ele_mp_5.p)

strain.mn.mp.lm<-lm(value~strain,icpms.m[icpms.m$treatment %in% "bhis" & icpms.m$metal %in% "Mn" & icpms.m$data_type %in% "mp",])
summary(strain.mn.mp.lm)
# none
strain.co.mp.lm<-lm(value~strain,icpms.m[icpms.m$treatment %in% "bhis" & icpms.m$metal %in% "Co" & icpms.m$data_type %in% "mp",])
summary(strain.co.mp.lm)
# cAB (0.00702)

allstrains_ele_uM_6.p<-ggplot(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("V","Cr","As","Mo") & icpms.m$data_type %in% "atcell",], aes(x=strain, y=value_uM,fill=strain))+geom_boxplot()+
  scale_y_continuous(name="Element (\U00B5M)")+
  scale_x_discrete(name="",breaks=breaks_strain,labels=labels_strain)+
  scale_fill_manual(values=colors_strain,breaks=breaks_strain,labels=labels_strain)+
  as_md_theme(my_theme)+theme(legend.position = "none")+facet_wrap(~metal)
plot(allstrains_ele_uM_6.p)

allstrains_ele_mp_6.p<-ggplot(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("V","Cr","As","Mo") & icpms.m$data_type %in% "mp",], aes(x=strain, y=value,fill=strain))+geom_boxplot()+
  scale_y_log10(name="Element/Phosphorus")+
  scale_x_discrete(name="",breaks=breaks_strain,labels=labels_strain)+
  scale_fill_manual(values=colors_strain,breaks=breaks_strain,labels=labels_strain)+
  as_md_theme(my_theme)+theme(legend.position = "none")+facet_wrap(~metal)
plot(allstrains_ele_mp_6.p)

allstrains_ele_uM_7.p<-ggplot(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("Ni","Se","Cd") & icpms.m$data_type %in% "atcell",], aes(x=strain, y=value_uM,fill=strain))+geom_boxplot()+
  scale_y_continuous(name="Element (\U00B5M)")+
  scale_x_discrete(name="",breaks=breaks_strain,labels=labels_strain)+
  scale_fill_manual(values=colors_strain,breaks=breaks_strain,labels=labels_strain)+
  as_md_theme(my_theme)+theme(legend.position = "none")+facet_wrap(~metal)
plot(allstrains_ele_uM_7.p)

wt_ele_low.p<-ggplot(icpms.m[icpms.m$strain %in% "wt" & icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("Cr","As","Mo","Ni","Se","Cd") & icpms.m$data_type %in% "atcell",], aes(x=metal, y=value_uM,fill=strain))+geom_boxplot(fill="black")+
  scale_y_continuous(name="Element (\U00B5M)")+
  scale_x_discrete(name="",limits=c("Mo","As","Cr","Se","Cd","Ni"))+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(wt_ele_low.p)

wt_ele_low_delta.p<-ggplot(icpms.m[icpms.m$strain %in% "wt" & icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("Cr","As","Mo","Ni","Se","Cd") & icpms.m$data_type %in% "atcell",], aes(x=metal, y=value_mM_deltamedia*1000,fill=strain))+geom_boxplot(fill="black")+
  scale_y_continuous(name="Intracellular - Extracellular Element (\u0394\U00B5M)")+
  scale_x_discrete(name="",limits=c("Mo","As","Cr","Se","Cd","Ni"))+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(wt_ele_low_delta.p)

ggsave(filename="figures/revised_202407/allstrain_allelements_1_mp.pdf",plot=allstrains_ele_mp_1.p,width=2.5,height=2,units="in",device=cairo_pdf)
ggsave(filename="figures/revised_202407/allstrain_allelements_2_mp.pdf",plot=allstrains_ele_mp_2.p,width=2.5,height=2,units="in",device=cairo_pdf)
ggsave(filename="figures/revised_202407/allstrain_allelements_3_mp.pdf",plot=allstrains_ele_mp_3.p,width=2.5,height=2,units="in",device=cairo_pdf)
ggsave(filename="figures/revised_202407/allstrain_allelements_4_mp.pdf",plot=allstrains_ele_mp_4.p,width=2.5,height=2,units="in",device=cairo_pdf)
ggsave(filename="figures/revised_202407/allstrain_allelements_5_mp.pdf",plot=allstrains_ele_mp_5.p,width=2.5,height=2,units="in",device=cairo_pdf)
ggsave(filename="figures/revised_202407/allstrain_allelements_6_mp.pdf",plot=allstrains_ele_mp_6.p,width=2.5,height=2,units="in",device=cairo_pdf)
ggsave(filename="figures/revised_202407/allstrain_allelements_7.pdf",plot=allstrains_ele_uM_7.p,width=3,height=3,units="in",device=cairo_pdf)
ggsave(filename="figures/revised_202407/wt_trace_delta.pdf",plot=wt_ele_low_delta.p,width=3,height=3,units="in",device=cairo_pdf)
ggsave(filename="figures/revised_202407/wt_trace_concentration.pdf",plot=wt_ele_low.p,width=3,height=3,units="in",device=cairo_pdf)

####
# Delta with exterior
allstrains_delta_mM_1.p<-ggplot(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("K") & icpms.m$data_type %in% "atcell",], aes(x=strain, y=value_mM_deltamedia,fill=strain))+geom_boxplot()+
  scale_y_continuous(name="Intracellular - Media Element (\u0394mM)")+
  scale_x_discrete(name="",breaks=breaks_strain,labels=labels_strain)+
  scale_fill_manual(values=colors_strain,breaks=breaks_strain,labels=labels_strain)+
  as_md_theme(my_theme)+theme(legend.position = "none")+facet_wrap(~metal)
plot(allstrains_delta_mM_1.p)

mean(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("K") & icpms.m$data_type %in% "atcell" & icpms.m$strain2 %in% "wt",]$value_mM_deltamedia) # 401.4 delta mM
mean(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("K") & icpms.m$data_type %in% "atcell" & icpms.m$strain2 %in% "wt",]$value_mM) # 434.1 mM

mean(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("K") & icpms.m$data_type %in% "atcell" & icpms.m$strain2 %in% "wt",]$value_mM) - mean(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("K") & icpms.m$data_type %in% "atcell" & icpms.m$strain2 %in% "cA",]$value_mM) # 65.5 mM
mean(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("K") & icpms.m$data_type %in% "atcell" & icpms.m$strain2 %in% "wt",]$value_mM) - mean(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("K") & icpms.m$data_type %in% "atcell" & icpms.m$strain2 %in% "cB",]$value_mM) # 97.9 mM

mean(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("K") & icpms.m$data_type %in% "atcell" & icpms.m$strain2 %in% "cAB",]$value_mM_deltamedia) # 471.3 mM (438.5 delta mM)
mean(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("K") & icpms.m$data_type %in% "atcell" & icpms.m$strain2 %in% "cAB",]$value_mM) - mean(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("K") & icpms.m$data_type %in% "atcell" & icpms.m$strain2 %in% "cA",]$value_mM) # 102.7 mM
mean(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("K") & icpms.m$data_type %in% "atcell" & icpms.m$strain2 %in% "cAB",]$value_mM) - mean(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("K") & icpms.m$data_type %in% "atcell" & icpms.m$strain2 %in% "cB",]$value_mM) # 135 mM

allstrains_delta_mM_2.p<-ggplot(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("Na") & icpms.m$data_type %in% "atcell",], aes(x=strain, y=value_mM_deltamedia,fill=strain))+geom_boxplot()+
  scale_y_continuous(name="Intracellular - Media Element (\u0394mM)",limits=c(-110,10),breaks=c(-100,-75,-50,-25,0))+
  scale_x_discrete(name="",breaks=breaks_strain,labels=labels_strain)+
  scale_fill_manual(values=colors_strain,breaks=breaks_strain,labels=labels_strain)+
  as_md_theme(my_theme)+theme(legend.position = "none")+facet_wrap(~metal)
plot(allstrains_delta_mM_2.p)

mean(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("Na") & icpms.m$data_type %in% "atcell" & icpms.m$strain2 %in% "wt",]$value_mM_deltamedia) # -66.7 delta mM
mean(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("Na") & icpms.m$data_type %in% "atcell" & icpms.m$strain2 %in% "wt",]$value_mM) # 62 mM

allstrains_delta_mM_3.p<-ggplot(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("Mg","S") & icpms.m$data_type %in% "atcell",], aes(x=strain, y=value_mM_deltamedia,fill=strain))+geom_boxplot()+
  scale_y_continuous(name="Intracellular - Media Element (\u0394mM)")+
  scale_x_discrete(name="",breaks=breaks_strain,labels=labels_strain)+
  scale_fill_manual(values=colors_strain,breaks=breaks_strain,labels=labels_strain)+
  as_md_theme(my_theme)+theme(legend.position = "none")+facet_wrap(~metal)
plot(allstrains_delta_mM_3.p)

allstrains_delta_mM_4.p<-ggplot(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("Fe","Ca") & icpms.m$data_type %in% "atcell",], aes(x=strain, y=value_mM_deltamedia,fill=strain))+geom_boxplot()+
  scale_y_continuous(name="Intracellular - Media Element (\u0394mM)")+
  scale_x_discrete(name="",breaks=breaks_strain,labels=labels_strain)+
  scale_fill_manual(values=colors_strain,breaks=breaks_strain,labels=labels_strain)+
  as_md_theme(my_theme)+theme(legend.position = "none")+facet_wrap(~metal)
plot(allstrains_delta_mM_4.p)

allstrains_delta_mM_5.p<-ggplot(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("Cu","Zn") & icpms.m$data_type %in% "atcell",], aes(x=strain, y=(value_mM_deltamedia)*1000,fill=strain))+geom_boxplot()+
  scale_y_continuous(name="Intracellular - Media Element (\u0394\U00B5M)")+
  scale_x_discrete(name="",breaks=breaks_strain,labels=labels_strain)+
  scale_fill_manual(values=colors_strain,breaks=breaks_strain,labels=labels_strain)+
  as_md_theme(my_theme)+theme(legend.position = "none")+facet_wrap(~metal)
plot(allstrains_delta_mM_5.p)

allstrains_delta_mM_6.p<-ggplot(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% c("Mn","Co") & icpms.m$data_type %in% "atcell",], aes(x=strain, y=(value_mM_deltamedia)*1000,fill=strain))+geom_boxplot()+
  scale_y_continuous(name="Intracellular - Media Element (\u0394\U00B5M)")+
  scale_x_discrete(name="",breaks=breaks_strain,labels=labels_strain)+
  scale_fill_manual(values=colors_strain,breaks=breaks_strain,labels=labels_strain)+
  as_md_theme(my_theme)+theme(legend.position = "none")+facet_wrap(~metal)
plot(allstrains_delta_mM_6.p)

allstrains_delta_mM_6.p<-ggplot(icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %!in% c("K","P","Na","Mg","S","Fe","Ca","Cu","Zn","Mn","Co") & icpms.m$data_type %in% "atcell",], aes(x=strain, y=(value_mM_deltamedia)*1000,fill=strain))+geom_boxplot()+
  scale_y_log10(name="Intracellular - Media Element (\u0394\U00B5M)")+
  scale_x_discrete(name="",breaks=breaks_strain,labels=labels_strain)+
  scale_fill_manual(values=colors_strain,breaks=breaks_strain,labels=labels_strain)+
  as_md_theme(my_theme)+theme(legend.position = "none")+facet_wrap(~metal)
plot(allstrains_delta_mM_6.p)

ggsave(filename="figures/revised_202407/allstrain_delta_1.pdf",plot=allstrains_delta_mM_1.p,width=3,height=3,units="in",device=cairo_pdf)
ggsave(filename="figures/revised_202407/allstrain_delta_2.pdf",plot=allstrains_delta_mM_2.p,width=3,height=3,units="in",device=cairo_pdf)
ggsave(filename="figures/revised_202407/allstrain_delta_3.pdf",plot=allstrains_delta_mM_3.p,width=3,height=3,units="in",device=cairo_pdf)
ggsave(filename="figures/revised_202407/allstrain_delta_4.pdf",plot=allstrains_delta_mM_4.p,width=3,height=3,units="in",device=cairo_pdf)
ggsave(filename="figures/revised_202407/allstrain_delta_5.pdf",plot=allstrains_delta_mM_5.p,width=3,height=3,units="in",device=cairo_pdf)
ggsave(filename="figures/revised_202407/allstrain_delta_6.pdf",plot=allstrains_delta_mM_6.p,width=3,height=3,units="in",device=cairo_pdf)

####

allstrains_s_p.p<-ggplot()+geom_point(aes(x=icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% "P" & icpms.m$data_type %in% "atcell",]$value_mM,y=icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% "S" & icpms.m$data_type %in% "atcell",]$value_mM))+
  geom_smooth(aes(x=icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% "P" & icpms.m$data_type %in% "atcell",]$value_mM,y=icpms.m[icpms.m$treatment %!in% c("dc","mon","nig") & icpms.m$metal %in% "S" & icpms.m$data_type %in% "atcell",]$value_mM),method=lm, color="black")+
  scale_y_continuous(name="S (mM)")+
  scale_x_continuous(name="P (mM)")+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(allstrains_s_p.p)

allstrains_s_p_all.p<-ggplot()+geom_point(aes(x=icpms.m[icpms.m$metal %in% "P" & icpms.m$data_type %in% "atcell",]$value_mM,y=icpms.m[icpms.m$metal %in% "S" & icpms.m$data_type %in% "atcell",]$value_mM))+
  geom_smooth(aes(x=icpms.m[icpms.m$metal %in% "P" & icpms.m$data_type %in% "atcell",]$value_mM,y=icpms.m[icpms.m$metal %in% "S" & icpms.m$data_type %in% "atcell",]$value_mM),method=lm, color="black")+
  scale_y_continuous(name="S (mM)")+
  scale_x_continuous(name="P (mM)")+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(allstrains_s_p_all.p)

s_p.lm <-lm(icpms.m[icpms.m$metal %in% "P" & icpms.m$data_type %in% "atcell",]$value_mM~icpms.m[icpms.m$metal %in% "S" & icpms.m$data_type %in% "atcell",]$value_mM)
summary(s_p.lm) # p < 2E-16

ggsave(filename="figures/revised_202407/phosphate_sulfer_comparison_allstrains_treat.pdf",plot=allstrains_s_p_all.p,width=2,height=2,units="in",device=cairo_pdf)

#####
# Relative abundance plots

ggplot(icpms.m[icpms.m$data_type %in% "atcell",],aes(x=variable,y=relabund_mM,fill=metal_1perc))+geom_bar(stat = "identity")+
  as_md_theme(my_theme2)

relabund_topelement.p<-ggplot(icpms.m_avgrel[icpms.m_avgrel$treatment %in% "bhis",],aes(x=strain,y=relabund_mM,fill=fct_reorder(metal_1perc,relabund_mM)))+geom_bar(stat = "identity")+
  scale_y_continuous(name="Relative Abundance (%)")+
  scale_x_discrete(name="",breaks=breaks_strain,labels=labels_strain)+
  scale_fill_manual(values = as.character(paletteMartin[1:6]))+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(relabund_topelement.p)

relabund_lowelement.p<-ggplot(icpms.m_avgrel[icpms.m_avgrel$treatment %in% "bhis" & icpms.m_avgrel$metal_1perc %in% "Other",],aes(x=strain,y=relabund_mM,fill=fct_reorder(metal, relabund_mM)))+geom_bar(stat = "identity")+
  scale_y_continuous(name="Relative Abundance (%)")+
  scale_x_discrete(name="",breaks=breaks_strain,labels=labels_strain)+
  scale_fill_manual(values = as.character(paletteMartin[1:13]))+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(relabund_lowelement.p)

ggsave(filename="figures/revised_202407/relabund_strain_topele.pdf",plot=relabund_topelement.p,width=2,height=2,units="in",device=cairo_pdf)
ggsave(filename="figures/revised_202407/relabund_strain_lowele.pdf",plot=relabund_lowelement.p,width=2,height=2,units="in",device=cairo_pdf)

relabund_treat_topelement.p<-ggplot(icpms.m_avgrel[icpms.m_avgrel$strain %in% "wt" & icpms.m_avgrel$treatment %in% c("bhis","dc","mon","nig"),],aes(x=treatment,y=relabund_mM,fill=fct_reorder(metal_1perc,relabund_mM)))+geom_bar(stat = "identity")+
  scale_y_continuous(name="Relative Abundance (%)")+
  scale_x_discrete(name="",breaks=c("bhis","dc","mon","nig"),labels=c("BHIS","Deoxycholate","Monensin","Nigericin"))+
  scale_fill_manual(values = as.character(paletteMartin[1:6]))+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(relabund_treat_topelement.p)

relabund_treat_lowelement.p<-ggplot(icpms.m_avgrel[icpms.m_avgrel$strain %in% "wt" & icpms.m_avgrel$treatment %in% c("bhis","dc","mon","nig") & icpms.m_avgrel$metal_1perc %in% "Other",],aes(x=treatment,y=relabund_mM,fill=fct_reorder(metal, relabund_mM)))+geom_bar(stat = "identity")+
  scale_y_continuous(name="Relative Abundance (%)")+
  scale_x_discrete(name="",breaks=c("bhis","dc","mon","nig"),labels=c("BHIS","Deoxycholate","Monensin","Nigericin"))+
  scale_fill_manual(values = as.character(paletteMartin[1:13]))+
  as_md_theme(my_theme)+theme(legend.position = "none")
plot(relabund_treat_lowelement.p)

ggsave(filename="figures/revised_202407/relabund_treat_topele.pdf",plot=relabund_treat_topelement.p,width=2,height=2,units="in",device=cairo_pdf)
ggsave(filename="figures/revised_202407/relabund_treat_lowele.pdf",plot=relabund_treat_lowelement.p,width=2,height=2,units="in",device=cairo_pdf)
