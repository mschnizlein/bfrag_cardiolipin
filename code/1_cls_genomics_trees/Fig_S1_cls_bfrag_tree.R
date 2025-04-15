## P207 intra-Bacteroides cls tree
## Author: Matthew Schnizlein
## 10/18/2024

### Contents
## Fasta and taxonomy file pre-processing
## Newick data file pre-processing (newick tree alignment generated in Geneious)
## Fig. 1E - neighbor end joining tree of cls genes from prominent gut taxa



library(tidyr)
library(ggplot2)
# library(reshape2)
# library(growthcurver)
library(ggtext)
library(gridtext)
library(mdthemes)
library(matrixStats)
library(data.table)
library(magrittr)
library(readxl)
library(scales)
library(colorBlindness)
library(RColorBrewer)
library(lme4)

library(treeio)
library(ggtree)
# library(glue)
# library(ComplexHeatmap)
# library(circlize)
# library(ggbiplot)
# library(ggforce)

setwd("C:/Users/mksch/Dropbox/matt_data/crosson_lab_msu_schnizlein/bfrag_bile_acid_stress/cls_conservation/")

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
substrLeft <- function(x, n){
  substr(x, 1, n)
}


'%!in%'<-function(x,y)!('%in%'(x,y))

my_theme<-theme_bw() + theme(axis.line=element_line(color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(), text=element_text(size = 12))
my_theme2<-theme_bw() + theme(axis.line=element_line(color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(), text=element_text(size = 12),axis.text.x = element_text(angle=30,hjust = 1))

my_theme_pres<-theme_bw() + theme(axis.line=element_line(color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(), text=element_text(size = 15))
my_theme2_pres<-theme_bw() + theme(axis.line=element_line(color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(), text=element_text(size = 15),axis.text.x = element_text(angle=30,hjust = 1))

#####
# Fasta cleanup/generation for Geneious Tree Algorithm
fasta.file<-read.table(file="protein-matching-TIGR04265.fasta",sep="\t",fill=NA)
tax.file<-read.table(file="taxonomy-matching-TIGR04265.tsv",sep = "\t",fill=NA,header=TRUE)

fasta.file2<-fasta.file[substrLeft(fasta.file$V1,1) %in% ">",]
colnames(fasta.file)<-"unedited"

fasta.file3<-as.data.frame(matrix(data=NA,nrow=4,ncol=451))

fasta.file3[]<-unlist(strsplit(fasta.file2,split="\\|"))
fasta.file3.t<-as.data.table(t(fasta.file3))
colnames(fasta.file3.t)<-c("protein_accession","reviewed_status","description","taxID")
fasta.file3.t$protein_accession<-gsub(">","",fasta.file3.t$protein_accession)
fasta.file3.t$taxID<-as.character(gsub("taxID:","",fasta.file3.t$taxID))

tax.file$Accession<-as.character(tax.file$Accession)

tax.file$Name2<-gsub("Bacteroides","",tax.file$Name) %>% gsub("\\[","",.) %>% gsub("\\]","",.) %>% gsub(" ","_",.) %>% gsub("Candidatus","",.)
tax.file$Name3<-gsub("_sp..*","_uncultured",tax.file$Name2) %>% gsub("'","",.) %>% gsub("__","_",.) %>% sub(".","",.) %>% gsub("_.*","",.)


tax.fasta.merge<-merge(x=fasta.file3.t,y=tax.file,by.x="taxID",by.y="Accession")



for (i in 1:length(tax.fasta.merge$protein_accession)){
fasta.file$unedited<-gsub(tax.fasta.merge$protein_accession[i],paste0(tax.fasta.merge$Name3[i],"_",tax.fasta.merge$protein_accession[i]),fasta.file$unedited)
}

fasta.file$unedited<-gsub("\\|unreviewed","",fasta.file$unedited) %>% gsub("\\|"," \\|",.)

# write.table(fasta.file,file="protein-matching-TIGR04265_cleaned.fasta",row.names = FALSE,col.names = FALSE,quote = FALSE)

#####
## Using Geneious Export function to create .newick formatted file (with support values) for use in R package treeio and plotting with ggtree

newick_tree_data<-read.newick(file = "../genetics_bfrag/cardiolipin_synthase_atpG/cls_tree_distance/bacteroides_cleaned_small_tree_20240724.newick")
# newick_tree_data_consensus<-read.newick(file="../genetics_bfrag/cardiolipin_synthase_atpG/cls_tree_distance/cls_protein_alignment_consensus_ecoli-outgroup.newick")
newick_tree_data_consensus<-read.newick(file="bacteroides_cls_consensus_nopectino_20240814.newick")

cls_phylo<-as.phylo(newick_tree_data)
cls_phylo_cons<-as.phylo(newick_tree_data_consensus)
# get.fields(cls_phylo)
# get.data(newick_tree_data)


svl_u <- as.matrix(cls_phylo[["tip.label"]])
svl_u_cons <- as.matrix(cls_phylo_cons[["tip.label"]])
# fit <- phytools::fastAnc(cls_phylo, svl_u, vars=TRUE, CI=TRUE)
svl_u_simple<-gsub("E coli MG1655 NC_000913 - ","Ecoli",cls_phylo["tip.label"]$tip.label) %>% gsub(" ","_",.) %>% gsub("_.*","",.) %>% gsub("'","",.)
svl_u_simple_cons<-gsub("E coli MG1655 NC_000913 - ","Ecoli",cls_phylo_cons["tip.label"]$tip.label) %>% gsub(" ","_",.) %>% gsub("_.*","",.) %>% gsub("'","",.)

# cls_phylo$tip.label

unique(svl_u)
unique(svl_u_simple)

td <- data.frame(node = nodeid(cls_phylo, cls_phylo[["tip.label"]]),
               trait = svl_u_simple)
td_cons <- data.frame(node = nodeid(cls_phylo_cons, cls_phylo_cons[["tip.label"]]),
                 trait = svl_u_simple_cons)
# nd <- data.frame(node = names(fit$ace), trait = fit$ace)

# d <- rbind(td, nd)
td$node <- as.numeric(td$node)
td_cons$node <- as.numeric(td_cons$node)
td_mod<-td
colnames(td_mod)<-c("node","strain")
td_mod_cons<-td_cons
colnames(td_mod_cons)<-c("node","strain")
tree_u <- full_join(cls_phylo, td, by = 'node')
tree_u_cons <- full_join(cls_phylo_cons, td_cons, by = 'node')


breaks_tree<-c("EcoliclsA","EcoliclsC","EcoliclsB","acidifaciens","faecium","caecimuris","ovatus","xylanisolvens","caccae","faecalis","faecichinchillae","faecis","finegoldii","thetaiotaomicron","fragilis","P207clsA","P207clsB","nordii","salyersiae","reticulotermitis","cellulosilyticus","intestinalis","oleiciplenus","stercorirosoris","clarus","eggerthii","stercoris","fluxus","helcogenes","muris","uniformis","heparinolyticus","pyogenes","graminisolvens","luti","avicola","merdipullorum","intestinipullorum","pullicola","intestinavium","merdigallinarum","merdavium","coprosuis","periocalifornicus","pectinophilus","NA")

colors_tree<-c(rep("#8BD1CB",3),rep("#000000",10),"#422CB2","#CC8E51","#006198","#299764",rep("#000000",29))

p1 <- ggtree(tree_u,aes(color=trait), size=0.75,color="black") +
  scale_color_manual(breaks=breaks_tree,values=colors_tree)+
#  geom_tiplab(geom='text',hjust = -.1) +
  geom_cladelab(data = td_mod,mapping = aes(node = node, label = strain, color = strain),
                fontsize = 3)+
  geom_treescale(x=0,y=200)+
  xlim(0, 1.2) +
  theme(text=element_text(size = 8),legend.position = "none")

plot(p1)

ggsave(filename = "selected_bacteroides_cls_tree.pdf",plot=p1,height=10,width =10,units="in",device=cairo_pdf)


p2_cons <- ggtree(tree_u_cons,aes(color=trait), size=0.75,color="black") +
  scale_color_manual(breaks=breaks_tree,values=colors_tree)+
  #  geom_tiplab(geom='text',hjust = -.1) +
  geom_cladelab(data = td_mod_cons,mapping = aes(node = node, label = strain, color = strain),
                fontsize = 3)+
  geom_treescale(x=0,y=200)+
  xlim(0, 1.2) +
  theme(text=element_text(size = 8),legend.position = "none")+
  geom_label2(aes(label=label,
                  subset = !is.na(as.numeric(label)) & as.numeric(label) > 90),label.size = 0.05,label.padding = unit(0.05,"lines"),nudge_x=-0.01,nudge_y = 1,size=2, fill=NA, color="black")

plot(p2_cons)

ggsave(filename = "selected_bacteroides_cls_tree_consensus.pdf",plot=p2_cons,height=10,width =8,units="in",device=cairo_pdf)







###
# Example Tree

library(ggtree)
library(treeio)
library(tidytree)
library(ggplot2)
library(TDbook)
## ref: http://www.phytools.org/eqg2015/asr.html
##
## load `tree_anole` and `df_svl` from 'TDbook'
svl <- as.matrix(df_svl)[,1]
fit <- phytools::fastAnc(tree_anole, svl, vars=TRUE, CI=TRUE)

td <- data.frame(node = nodeid(tree_anole, names(svl)),
                 trait = svl)
nd <- data.frame(node = names(fit$ace), trait = fit$ace)

d <- rbind(td, nd)
d$node <- as.numeric(d$node)
tree <- full_join(tree_anole, d, by = 'node')

p1 <- ggtree(tree, aes(color=trait), layout = 'circular',
             ladderize = FALSE, continuous = 'colour', size=2) +
  scale_color_gradientn(colours=c("red", 'orange', 'green', 'cyan', 'blue')) +
  geom_tiplab(hjust = -.1) +
  xlim(0, 1.2) +
  theme(legend.position.inside =  c(.05, .85))

p2 <- ggtree(tree, layout='circular', ladderize = FALSE, size=2.8) +
  geom_tree(aes(color=trait), continuous = 'colour', size=2) +
  scale_color_gradientn(colours=c("red", 'orange', 'green', 'cyan', 'blue')) +
  geom_tiplab(aes(color=trait), hjust = -.1) +
  xlim(0, 1.2) +
  theme(legend.position = c(.05, .85))

plot_list(p1, p2, ncol=2, tag_levels="A")

