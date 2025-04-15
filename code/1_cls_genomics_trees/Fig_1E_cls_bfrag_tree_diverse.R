## P207 inter-phylum cls tree
## Author: Matthew Schnizlein
## 10/18/2024

### Contents
## Newick data file pre-processing (alignment file generated in Geneious)
## Fig. 1E - neighbor end joining tree of cls genes from prominent gut taxa

library(tidyr)
library(ggplot2)
# library(reshape2)
library(growthcurver)
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

setwd("~/data/01_cls_conservation")

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

## Using Geneious Export function to create .newick formatted file (with support values) for use in R package treeio and plotting with ggtree

newick_tree_data_diverse_consensus<-read.newick(file = "diverse_cls_consensus_tree.newick")
newick_tree_data_diverse<-read.newick(file = "diverse_cls_tree.newick")

cls_phylo_diverse<-as.phylo(newick_tree_data_diverse)
cls_phylo_diverse_consensus<-as.phylo(newick_tree_data_diverse_consensus)
# get.fields(cls_phylo)
# get.data(newick_tree_data)


svl_u_simple<-gsub("NC_010729 Porphyromonas gingivalis","P. gingivalis",as.matrix(cls_phylo_diverse[["tip.label"]])) %>%
  gsub("NZ_CP117963 - F prausnitzii","F. prausnitzii",.) %>%
  gsub("NZ_CP070367 - Sphingomonas paucimobilis","S. paucimobilis",.) %>%
  gsub("NZ_CP056117","E. cloacae",.) %>%
  gsub("NZ_CP053187","T. sanguinis",.) %>%
  gsub("NZ_CP053028","P. aeruginosa",.) %>%
  gsub("NZ_CP036170","C. scindens",.) %>%
  gsub("NZ_CP028109 - Fusobacterium nucleatum","F. nucleatum",.) %>%
  gsub("NZ_CP026010","B. subtilis",.) %>%
  gsub("NZ_CP020020","S. aureus",.) %>%
  gsub("NZ_CP008926 Strep pyogenes","S. pyogenes",.) %>%
  gsub("NZ_CP008816","E. faecalis",.) %>%
  gsub("NZ_CP006659","K. pneumoniae",.) %>%
  gsub("NC_021040 Roseburia","R. intestinalis",.) %>%
  gsub("Bfragilis_NCTC9343_NC_003228","B. fragilis NCTC9343",.) %>%
  gsub("Bfragilis_p207 extraction", "B. fragilis P207",.) %>%
  gsub("Btheta_VPI-5482","B. thetaiotaomicron",.) %>%
  gsub("Dielma fastidiosa","D. fastidiosa",.) %>%
  gsub("A0A2X2JCM9\\|sphingobacterium multivorum\\|","S. multivorum ",.) %>%
  gsub("E coli MG1655 NC_000913","E. coli",.) %>%
  gsub("NC_002935 - Corynebacterium diphteriae","C. diptherieae",.) %>%
  gsub("F0EU63\\|haemophilus parainfluenza\\|","H. parainfluenza ",.) %>%
  gsub("D2ZZQ4\\|Neisseria mucosa ATCC25996\\|","N. mucosa ",.) %>%
  gsub("NC_008054 - Lactobacillus delbreuckii","L. delbreuckii",.) %>%
  gsub("'","",.) %>%
  gsub("cls CDS translation 1","cls1",.) %>%
  gsub("cls CDS translation 2","cls2",.) %>%
  gsub("cls CDS translation 3","cls3",.) %>%
  gsub("cls CDS translation 4","cls4",.) %>%
  gsub("cls CDS translation 5","cls5",.) %>%
  gsub("cls CDS translation","cls1",.) %>%
  gsub("phosphatidylserine/phosphatidylglycerophosphate/cardiolipin synthase family protein CDS translation","cls1",.) %>%
  gsub(" 3 CDS translation","3",.) %>%
  gsub(" 2 CDS translation","2",.) %>%
  gsub(" CDS translation","",.) %>%
  gsub("(cls 1) - cls1","cls1",.) %>%
  gsub("(cls 2) - cls1","cls2",.) %>%
  gsub("cardiolipin","cls1",.,ignore.case=TRUE)


svl_u_simple_consensus<-gsub("NC_010729 Porphyromonas gingivalis","P. gingivalis",as.matrix(cls_phylo_diverse_consensus[["tip.label"]])) %>%
  gsub("NZ_CP117963 - F prausnitzii","F. prausnitzii",.) %>%
  gsub("NZ_CP070367 - Sphingomonas paucimobilis","S. paucimobilis",.) %>%
  gsub("NZ_CP056117","E. cloacae",.) %>%
  gsub("NZ_CP053187","T. sanguinis",.) %>%
  gsub("NZ_CP053028","P. aeruginosa",.) %>%
  gsub("NZ_CP036170","C. scindens",.) %>%
  gsub("NZ_CP028109 - Fusobacterium nucleatum","F. nucleatum",.) %>%
  gsub("NZ_CP026010","B. subtilis",.) %>%
  gsub("NZ_CP020020","S. aureus",.) %>%
  gsub("NZ_CP008926 Strep pyogenes","S. pyogenes",.) %>%
  gsub("NZ_CP008816","E. faecalis",.) %>%
  gsub("NZ_CP006659","K. pneumoniae",.) %>%
  gsub("NC_021040 Roseburia","R. intestinalis",.) %>%
  gsub("Bfragilis_NCTC9343_NC_003228","B. fragilis NCTC9343",.) %>%
  gsub("Bfragilis_p207 extraction", "B. fragilis P207",.) %>%
  gsub("Btheta_VPI-5482","B. thetaiotaomicron",.) %>%
  gsub("Dielma fastidiosa","D. fastidiosa",.) %>%
  gsub("A0A2X2JCM9\\|sphingobacterium multivorum\\|","S. multivorum ",.) %>%
  gsub("E coli MG1655 NC_000913","E. coli",.) %>%
  gsub("NC_002935 - Corynebacterium diphteriae","C. diptherieae",.) %>%
  gsub("F0EU63\\|haemophilus parainfluenza\\|","H. parainfluenza ",.) %>%
  gsub("D2ZZQ4\\|Neisseria mucosa ATCC25996\\|","N. mucosa ",.) %>%
  gsub("NC_008054 - Lactobacillus delbreuckii","L. delbreuckii",.) %>%
  gsub("'","",.) %>%
  gsub("cls CDS translation 1","cls1",.) %>%
  gsub("cls CDS translation 2","cls2",.) %>%
  gsub("cls CDS translation 3","cls3",.) %>%
  gsub("cls CDS translation 4","cls4",.) %>%
  gsub("cls CDS translation 5","cls5",.) %>%
  gsub("cls CDS translation","cls1",.) %>%
  gsub("phosphatidylserine/phosphatidylglycerophosphate/cardiolipin synthase family protein CDS translation","cls1",.) %>%
  gsub(" 3 CDS translation","3",.) %>%
  gsub(" 2 CDS translation","2",.) %>%
  gsub(" CDS translation","",.) %>%
  gsub("(cls 1) - cls1","cls1",.) %>%
  gsub("(cls 2) - cls1","cls2",.) %>%
  gsub("cardiolipin","cls1",.,ignore.case=TRUE)
# cls_phylo$tip.label

# unique(svl_u)
unique(svl_u_simple)
unique(svl_u_simple_consensus)

td <- data.frame(node = nodeid(cls_phylo_diverse, cls_phylo_diverse[["tip.label"]]),
               trait = svl_u_simple)

td_consensus <- data.frame(node = nodeid(cls_phylo_diverse_consensus, cls_phylo_diverse_consensus[["tip.label"]]),
                 trait = svl_u_simple_consensus)
# nd <- data.frame(node = names(fit$ace), trait = fit$ace)

# d <- rbind(td, nd)
td$node <- as.numeric(td$node)
td_consensus$node <- as.numeric(td_consensus$node)
td_mod<-td
colnames(td_mod)<-c("node","strain")
td_mod_consensus<-td_consensus
colnames(td_mod_consensus)<-c("node","strain")
tree_u <- full_join(cls_phylo_diverse, td, by = 'node')
tree_u_consensus <- full_join(cls_phylo_diverse_consensus, td_consensus, by = 'node')

breaks_tree<-c(
  # Bacteroidota
  ## Bacteroidia (8)
  "P. gingivalis - cls1",
  "B. thetaiotaomicron - cls1","B. thetaiotaomicron - cls2","B. thetaiotaomicron - cls3",
  "B. fragilis NCTC9343 - cls1","B. fragilis NCTC9343 - cls2","B. fragilis P207 (cls 1) - cls1","B. fragilis P207 (cls 2) - cls1",
  ## Sphingobacteria (2)
  "S. multivorum cls1",
  "S. paucimobilis - cls1",
  # Bacillota
  ## Bacilli (9)
  "S. aureus - cls1","S. aureus - cls2",
  "B. subtilis - cls1","B. subtilis - cls2","B. subtilis - cls3",
  "S. pyogenes - cls1",
  "E. faecalis - cls1","E. faecalis - cls2",
  "L. delbreuckii cls1",
  ## Clostridia (6)
  "F. prausnitzii cls1",
  "R. intestinalis - cls1","R. intestinalis - cls2","R. intestinalis - cls3",
  "C. scindens - cls1","C. scindens - cls2",
  ## Erysipelotrichia (6)
  "D. fastidiosa cls",
  "T. sanguinis - cls1","T. sanguinis - cls2","T. sanguinis - cls3","T. sanguinis - cls4","T. sanguinis - cls5",
  # Actinomycetota
  ## Actinomycetia (1)
  "C. diptherieae cls1",
  # Fusobacteriota
  ## Fusobacteria (1)
  "F. nucleatum cls1",
  # Pseudomonadota
  ## Betaproteobacteria (1)
  "N. mucosa cls1",
  ## Gammaproteobacteria
  ### Enterbacteriales (6)
  "E. coli - clsA","E. coli - clsB",
  "E. cloacae - cls1","E. cloacae - clsB",
  "K. pneumoniae - cls1","K. pneumoniae - clsB",
  ### Pasteurellales (1)
  "H. parainfluenza cls1",
  ### Pseudomonadales (2)
  "P. aeruginosa - cls1","P. aeruginosa - clsB")

colors_tree<-c(
  # Bacteroidota (green)
  rep("#074700",8),
  rep("#479d3d",2),
  # Bacillota (blue)
  rep("#00164d",9),
  rep("#17679d",6),
  rep("#78d4f2",6),
  # Actinomycetota (orange)
  "#a66608",
  # Fusobacteriota (yellow)
  "#e28425",
  # Pseudomonadota (red)
  "#730105",
  rep("#9b2220",6),
  "#b73431",
  rep("#c53c3a",2))


# colors_tree<-c(rep("#8BD1CB",3),rep("#000000",10),"#422CB2","#CC8E51","#006198","#299764",rep("#000000",29))

p1_diverse <- ggtree(tree_u,aes(color=trait), size=0.75,color="black") +
  scale_color_manual(breaks=breaks_tree,values=colors_tree)+
#  geom_tiplab(geom='text',hjust = -.1) +
  geom_cladelab(data = td_mod,mapping = aes(node = node, label = strain, color = strain),
                fontsize = 3)+
  geom_treescale()+
  theme(text=element_text(size = 8),legend.position = "none")

plot(p1_diverse)

ggsave(filename = "diverse_cls_tree.pdf",plot=p1_diverse,height=10,width =10,units="in",device=cairo_pdf)

p1_diverse_consensus <- ggtree(tree_u_consensus,aes(color=trait), size=0.75,color="black") +
  scale_color_manual(breaks=breaks_tree,values=colors_tree)+
  #  geom_tiplab(geom='text',hjust = -.1) +
  geom_cladelab(data = td_mod_consensus,mapping = aes(node = node, label = strain, color = strain),
                fontsize = 3)+
  geom_treescale()+
  geom_label2(aes(label=label,
                  subset = !is.na(as.numeric(label)) & as.numeric(label) > 80))+
  theme(text=element_text(size = 8),legend.position = "none")

plot(p1_diverse_consensus)

ggsave(filename = "diverse_cls_tree_consensus.pdf",plot=p1_diverse_consensus,height=6,width =7,units="in",device=cairo_pdf)






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

