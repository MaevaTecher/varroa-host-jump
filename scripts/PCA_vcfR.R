#!/usr/bin/env Rscript
.libPaths("/home/m/maeva-techer/R/x86_64-pc-linux-gnu-library/3.4")

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}

infile <- args[1]
listpop <- args[2]
outfile <- args[3]


library(adegenet)
library(vcfR)
library(ggplot2)
library(ggmap) 
library(poppr)
library(RColorBrewer)
library(gridExtra)
library(igraph)
library(StAMPP)
library(lattice)

## Import vcf file
vcf <- read.vcfR(infile)

## Convert into genlight
gl.varroa <- vcfR2genlight(vcf)
gl.varroa

## assign population
meta <- read.table(listpop, header = TRUE)
head(meta)
gl.varroa$pop <- meta$POP

##perform the pca
pca1 <- glPca(gl.varroa, nf = 5)
pca1

varroa.pca.scores <- as.data.frame(pca1$scores)
varroa.pca.scores$pop <- pop(gl.varroa)

myCol <- c("#FF1300", 	#VDAC
           "#ff8880",	#VDAM
           "darkred",	#VDMIS
           "#0040FF",	#VJAC
           "skyblue1",	#VJAM
           "green") 	#VJMIS

axe1 <- paste(round(pca1$eig[1]/sum(pca1$eig)*100, digits = 2))
axe1
axe2 <- paste(round(pca1$eig[2]/sum(pca1$eig)*100, digits = 2))
axe2

##export pdf
pdf(outfile, height=5, width= 7)
p <- ggplot(varroa.pca.scores, aes(x=PC1, y=PC2, colour=pop)) 
p <- p + geom_point(aes(fill=pop), colour="black", shape = 21, size=4, alpha = 0.9)
p <- p + scale_fill_manual(values = myCol) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw()
p <- p + theme(plot.background = element_blank()
  ,panel.grid.major = element_blank()
  ,panel.grid.minor = element_blank())
p <- p + labs(title="PCA on 44 varroa samples with LD prune SNPs")
p

barplot(pca1$eig, col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Proportion of variance explained")
title(xlab="Eigenvalue")

dev.off()

