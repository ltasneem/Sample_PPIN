ssn_degree_cen = function(d3,g){
  
  gbt <- readRDS("Network_btwns_centrality.Rda")
  a <- readRDS("gene_1080_1047.Rda")
  gd1 <- data.frame(matrix(unlist(gbt), nrow=1047, byrow=T))
  colnames(gd1)[1] <- "btwns_cntrl"
  gd1$num <- c(1:1047)
  gd = gd1[order(gd1$btwns_cntrl,decreasing = TRUE),]
  gd_test = gd[1:100,1:2]
  saveRDS(gd_test,"Ntwk_btwn_cntrl_100_top.Rda")
  slist = list()
  for(i in 1:length(gd_test$num))
  {
    slist[[i]] = prob$samples[gd_test$num[i]]
  }
  slist = unlist(slist)
  atop100 =  a[,colnames(a)%in% slist]
  saveRDS(atop100,"atop100.Rda")
  
  g100_vname_list = list()
  d100_list=list()
  
  for(i in 1:100)
  {
    s1 <-atop100[,c(i)]
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
    g100_vname_list[[i]] <-V(g1)$name
    ls <- centr_degree(g1,normalized = TRUE)
    d100_list[[i]] <-ls$res
  }
  
  saveRDS(g100_vname_list,"g100_vname_list.Rda")
  saveRDS(d100_list,"Ntwk_deglist_100_top.Rda")
}

g100_50_mat = function(d3)
{
  g100_vname_list = readRDS("g100_vname_list.Rda")
  pgene = data.frame(gene_name = prob$genes)
  d100_list= readRDS("Ntwk_deglist_100_top.Rda")
  
  g100_gene_degree = list()
  
  for(i in 1:100)
  {
    vname_degree <- do.call(rbind, Map(data.frame, vname=g100_vname_list[[i]], degree=d100_list[[i]]))
    vdegree_gene_protein_map<- merge(vname_degree, d3, by.x='vname', by.y='ensembl_protein_id')
    vdegree_gene_map= vdegree_gene_protein_map[,c("hgnc_symbol","degree")]
    vdegree_gene <- vdegree_gene_map[!duplicated(vdegree_gene_map$hgnc_symbol),]
    zz <- merge(vdegree_gene, pgene,by.x='hgnc_symbol', by.y='gene_name',all = TRUE)
    zz[is.na(zz)] <- 0
    g100_gene_degree[[i]]=zz
  }
  
  saveRDS(g100_gene_degree,"g100_gene_degree.Rda")
  
  mt = sapply(g100_gene_degree, function( x ) x$degree )
  
  return(mt)
  
}

graph_top50_gene = function(mt){
  var_mt = apply(mt, 1,var )
  saveRDS(var_mt,"g_variance.Rda")
  
  vmt <- data.frame(matrix(var_mt, nrow=length(var_mt), byrow=T))
  
  colnames(vmt)[1] <- "var"
  vmt$num <- c(1:length(var_mt))
  vmg =vmt[order(vmt$var,decreasing = TRUE),]
  vmg_test = vmg[1:50,1:2]
  saveRDS(vmg_test,"top_50_var.Rda")
  gtop50list = list()
  for(i in 1:length(vmg_test$num))
  {
    gtop50list[[i]] = prob$genes[vmg_test$num[i]]
  }
  gtop50list = unlist(gtop50list)
  
  saveRDS(gtop50list,"gtop50list.Rda")
  return(gtop50list)
}

heatmap_hclust = function(mt1_matrix){
 
  y<-mt1_matrix  
  hr <- hclust(as.dist(1-cor(t(y), method="pearson")), method="centroid"); 
  hc <- hclust(as.dist(1-cor(y, method="spearman")), method="centroid")  
  mycl <- cutree(hr,k=4, h=max(hr$height)/1.5); 
  mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)] ; myheatcol <- bluered(75)
  
  heatmap.2(y, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=myheatcol, scale="row", density.info="none", 
            trace="none", RowSideColors=mycolhc,margins=c(10,10))
  
}

heatmap_kmeans = function( mt1_matrix){
  set.seed(100)
  m = mt1_matrix
  km = kmeans(m, 4)
  m2 <- cbind(m,km$cluster)
  o <- order(m2[, 101])
  m2 <- m2[o, ]
  pheatmap(m2[,1:100], cluster_rows=F,cluster_cols=F, col=brewer.pal(4,"Set3"),border_color=NA)
  
}

validation_cluster = function( mt1_matrix){
  clmethods <- c("hierarchical","kmeans")
  intern <- clValid(mt1_matrix, nClust = 2:10,
                    clMethods = clmethods, validation = "internal")
  # Summary
  summary(intern)
  
  stab <- clValid(mt1_matrix, nClust = 2:10, clMethods = clmethods,
                  validation = "stability")
  # Display only optimal Scores
  optimalScores(stab)
  
}