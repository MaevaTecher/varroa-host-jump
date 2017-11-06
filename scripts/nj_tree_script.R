library("ape")
varroa <- read.table("varroa.phylip", header=TRUE, sep="\t", row.names=1)
varroa1 <- as.dist(m, "matrix")
njvarroa <- nj(varroa1)

pdf(file="varroa_nj_tree.pdf", width=10, height=10)
plot(njvarroa, type = "phylogram")
dev.off()

