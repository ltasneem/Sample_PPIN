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

##extracting sample column from gene * sample 2D matrix
##and list those 0 expression value genes
s1 <- a[,c(1)]
x <- names(s1[which(s1==0)])

## List of ppid from gene_symboltoppid list
##for 0 expression value genes

x2 <- d3[which(d3$hgnc_symbol == x[1]),]
x2 <- x2$ensembl_protein_id
for(i in 2:length(x))
{
  x1 <- d3[which(d3$hgnc_symbol == x[i]),]
  x1 <- x1$ensembl_protein_id
  if(length(x1) != 0)
  {
    x3 <- list(x2,x1)
    vec <- unlist(x3)
    vec <- vec[which(c(1,diff(vec)) != 0)]
    x2 <- vec
  }
}

## removing those ppid from row and column of 
##ppi v*v adjacency matrix.

ppiwithoutgr = ppi[!row.names(ppi)%in% x2,]
ppiwithoutgcr = ppiwithoutgr[,!colnames(ppiwithoutgr)%in% x2]

