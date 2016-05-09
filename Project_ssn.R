create_sample_specific_network = function(k){
  
  ##creating ppi v*v adjacency matrix
  
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
  
  return(g)
}
  
symbol_id_conversion = function(dataset){
  bf<- row.names(dataset)
  gn<-t(t(bf))
  colnames(gn)[1] <- "gene"
  gn
  saveRDS(gn,"sample_genename.Rda")
  ab<- readRDS("sample_genename.Rda")
  G_list <- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","ensembl_gene_id"),values=ab,mart= mart)
  saveRDS(G_list,"sample_genesymbol2idmap.Rda")
  return(G_list)
}
  


graph_list = function(g_ex,g)
{
  glist=list()
  
  for(i in 1:ncol(g_ex))
  {
    s1 <-g_ex[,c(i)]
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
    glist[[i]] <- g1
  }
  
  return(glist)
}



g_prop = function(graph1, s)
{
  g_g = list()
  if(s == 'trans')
  {
    for(i in 1:length(graph1))
    {
      g_g[[i]] <- transitivity(graph1[[i]]);
    }
  }
  
  if(s == 'dens')
  {
    for(i in 1:length(graph1))
    {
      g_g[[i]] <- graph.density(graph1[[i]]);
    }
  }
  
  if(s == 'pathlen')
  {
    for(i in 1:length(graph1))
    {
      g_g[[i]] <- average.path.length(graph1[[i]]);
    }
  }
  
  if(s == 'modularity')
  {
    for(i in 1:length(graph1))
    {
      wtc <- cluster_walktrap(graph1[[i]])
      g_g[[i]] <- modularity(wtc)
      
    }
  }
  
  if(s == 'centrality')
  {
    for(i in 1:length(graph1))
    {
      ls <- centr_betw(graph1[[i]])
      g_g[[i]] <- ls$centralization
      
    }
  }
  
  return(g_g);
}
