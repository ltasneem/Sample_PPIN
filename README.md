
Steps to run the project:

1.These 3 files are needed to run this project:

Project_main.R
Project_feature_cluster.R
Project_ssn.R

2.protein_links.txt and prob.BRCA.N.exp.RData these
two data file must be downloaded from the original source
mentioned in the paper.

k <- read.table("protein_links.txt", sep = " ", header = TRUE)
pb <- load("prob.BRCA.N.exp.RData")

3.These three files must be kept in the same workspace
where the source code is. These files hold data with some
basic calculations needed throughout the project.

Network_btwns_centrality.Rda
gene_1080_1047.Rda
gid_pid_gsymbol.Rda

Other necessary file will be generated by the source code.

