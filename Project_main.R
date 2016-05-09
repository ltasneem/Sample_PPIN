##libraries
library(biomaRt)
library(dplyr)
library(igraph)
library(MASS)
library(Hmisc)
library(clValid)
library(reshape2)
library(ggplot2)
library(gplots)
library(pheatmap) 
library(RColorBrewer)

source('Project_ssn.R') 
source('Project_feature_cluster.R')

k <- read.table("protein_links.txt", sep = " ", header = TRUE)
g <- create_sample_specific_network(k)
## gene * sample 2D matrix
pb <- load("prob.BRCA.N.exp.RData")
a <- prob$mRNAq2
dimnames(a) = list( 
  +       c(prob$genes),   c(prob$samples))
saveRDS(a,"gene_1080_1047.Rda")

##gene_symboltoppid list
gs <- data.frame(genes = prob$genes)
G_list <- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","ensembl_gene_id"),values=gs,mart= mart)
d <- readRDS("sample_gene2protein.Rda")
d3 = merge(G_list,d)
saveRDS(d3,"gid_pid_gsymbol.Rda")

d3 <- readRDS("gid_pid_gsymbol.Rda")
a <- readRDS("gene_1080_1047.Rda")

graphlist = graph_list(a,g)
saveRDS(graphlist,"Graph_list.Rda")
gdensity = g_prop(graphlist,'dens')
g_transitivity = g_prop(graphlist,'trans')
g_pathlength = g_prop(graphlist,'pathlen');
gmt = g_prop(graphlist,'modularity')
gbt = g_prop(graphlist,'centrality')

tmp=prob$mRNA.n
tmp1=tmp> apply(tmp, 1, quantile, probs = 1/3,  na.rm = TRUE)
tmp2=tmp> apply(tmp, 1, quantile, probs = 2/3,  na.rm = TRUE)
a_n=tmp1+tmp2

rownames(prob$mRNA) = prob$genes
colnames(prob$mRNA) = prob$samples
rownames(prob$mRNA.n) = prob$genes
colnames(prob$mRNA.n) = prob$samples.n

glist_nor = graph_list(a_n,g)
saveRDS(glist_nor,"Graph_list_normal.Rda")
gdensity_nor = g_prop(glist_nor,'dens')
gtrans_nor = g_prop(glist_nor,'trans')
gpathlen_nor = g_prop(glist_nor,'pathlen');
gmodularity_nor = g_prop(glist_nor,'modularity')
gcentrality_nor = g_prop(glist_nor,'centrality')

##getting gene by sample matrix of degree centrality
##in base of network centrality

ssn_degree_cen()
mt = g100_50_mat(d3)
g100_50_mat(mt)

##heatmap of hierarchical and k-means clustering 

mt1 = mt
rownames(mt1) = rownames(atop100)
colnames(mt1) = colnames(atop100)
mt1_100_50 = mt1[row.names(mt1)%in% gtop50list,]
mt1_100_48 = mt1_100_50[1:48,1:100]
mt1_matrix <- data.matrix(mt1_100_48)

heatmap_hclust(mt1_matrix)
heatmap_kmeans(mt1_matrix)

##validation

validation_cluster( mt1_matrix)



