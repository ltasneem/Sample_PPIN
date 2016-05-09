
Steps to run the project:

1.These 3 files are needed to run this project:

(1) Project_main.R
(2) Project_feature_cluster.R
(3) Project_ssn.R

2.protein_links.txt and  prob.BRCA.N.exp.RData, these
two data file must be downloaded from the original source
mentioned in the paper.These two files are too big to add
here.

1) k <- read.table("protein_links.txt", sep = " ", header = TRUE)
2) pb <- load("prob.BRCA.N.exp.RData")

3.These three files must be kept in the same workspace
where the source code is. These files hold data with some
basic calculations needed throughout the project.

1) Network_btwns_centrality.Rda
2) gene_1080_1047.Rda
3) gid_pid_gsymbol.Rda

Other necessary file will be generated by the source code.

