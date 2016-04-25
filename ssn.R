##libraries
library(biomaRt)
library(dplyr)
library(igraph)
library(MASS)
library(Hmisc)

##creating ppi v*v adjacency matrix

k <- read.table("protein_links.txt", sep = " ", header = TRUE)
k12 <- k[,c(1,2,10)]
write.csv(k12,"edge_list_weight_8548002.csv")
saveRDS(k12,"edge_list_weight_8548002.Rda")
k12_new <- k12[ k12$combined_score >= 700, ]
k12_new <- k12[ k12$combined_score >= 950, ]
k12_new<-mutate(k12_new,protein1=as.character(protein1))
df2<-mutate(k12_new,protein1=sapply(strsplit(k12_new$protein1, split='.', fixed=TRUE),function(x) (x[2])))
df2<-mutate(df2,protein2=as.character(protein2))
df3<-mutate(df2,protein2=sapply(strsplit(df2$protein2, split='.', fixed=TRUE),function(x) (x[2])))
write.csv(k12,"edge_list_weight_700.csv")
el <- df3
el[,1]=as.character(el$protein1)
el[,2]=as.character(el$protein2)
el=as.matrix(el)
g=graph.edgelist(el[,1:2],directed=FALSE)
f <- df3$combined_score
E(g)$combined_score=as.numeric(f)
adj=get.adjacency(g,attr='combined_score',sparse=FALSE)

## gene * sample 2D matrix
pb <- load("prob.BRCA.N.exp.RData")
a <- prob$mRNAq2
dimnames(a) = list( 
  +       c(prob$genes),   c(prob$samples))

##gene_symboltoppid list
gs <- data.frame(genes = prob$genes)
G_list <- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","ensembl_gene_id"),values=gs,mart= mart)
d <- readRDS("sample_gene2protein.Rda")
d3 = merge(G_list,d)

################# ssn matrix building

d3 <- readRDS("gid_pid_gsymbol.Rda")
a <- readRDS("gene_1080_1047.Rda")
ppi <- readRDS("ppi_950.Rda")


## Alternate ssn building & feature extraction
g_density <- list()
for(i in 1:ncol(a))
{
s1 <- a[,c(i)]
s2 <- s1[which(s1==0)]
op <- t(t(s2))
p1 <- as.data.frame(as.table(op))
keeps <- c("Var1", "Freq")
p2 <- p1[keeps]
colnames(p2)[1] <- "gene_name"
colnames(p2)[2] <- "expression"
zero_ex <- merge(p2, d3, by.x='gene_name', by.y='hgnc_symbol')
g1 <- g
inter1 <- intersect(df3$protein1,zero_ex$ensembl_protein_id)
g1 <- delete_vertices(g1, inter1)
x2 <- graph.density(g1)
g_density[[i]] <- x2
}

> saveRDS(g_density,"Network_density.Rda")
> write.csv(g_density,"Network_density.csv")

########################################################

g_transitivity=list()

for(i in 1:ncol(a))
{
s1 <- a[,c(i)]
s2 <- s1[which(s1==0)]
op <- t(t(s2))
p1 <- as.data.frame(as.table(op))
keeps <- c("Var1", "Freq")
p2 <- p1[keeps]
colnames(p2)[1] <- "gene_name"
colnames(p2)[2] <- "expression"
zero_ex <- merge(p2, d3, by.x='gene_name', by.y='hgnc_symbol')
g1 <- g
inter1 <- intersect(df3$protein1,zero_ex$ensembl_protein_id)
g1 <- delete_vertices(g1, inter1)
x2 <- transitivity(g1)
g_transitivity[[i]] <- x2
}

> saveRDS(g_transitivity,"Network_transitivity.Rda")
> write.csv(g_transitivity,"Network_transitivity.csv")

##############################################################
g_density=list()
transitivity=list()
betweenness_centrality=list()

for(i in 1:ncol(a))
{
s1 <- a[,c(i)]
s2 <- s1[which(s1==0)]
op <- t(t(s2))
p1 <- as.data.frame(as.table(op))
keeps <- c("Var1", "Freq")
p2 <- p1[keeps]
colnames(p2)[1] <- "gene_name"
colnames(p2)[2] <- "expression"
zero_ex <- merge(p2, d3, by.x='gene_name', by.y='hgnc_symbol')
g1 <- g
inter1 <- intersect(df3$protein1,zero_ex$ensembl_protein_id)
g1 <- delete_vertices(g1, inter1)
x2 <- graph.density(g1)
g_density[[i]] <- x2
x2 <- transitivity(g1)
transitivity[[i]] <- x2
ls <- centr_betw(g1)
x2 <- ls$centralization
betweenness_centrality[[i]] <- x2
}

