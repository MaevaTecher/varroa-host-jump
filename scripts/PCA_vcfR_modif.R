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
#coloroption <- args[4]
#shapeoption <- args[5]


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
pca.var <- glPca(gl.varroa, nf = 4)
pca.var

varroa.pca.scores <- as.data.frame(pca.var$scores)
varroa.pca.scores$pop <- pop(gl.varroa)

## For VDAC "#FF1300"
## For VDAM "#ff8880"
## For VDMIS "darkred"
## For VJAC "#0040FF"
## For VJAM "skyblue1"
## For VJMIS "green"

myCol.all <- c("#FF1300",   #VDAC-CHN  
           "#FF1300",   #VDAC-KOR  
           "#FF1300",   #VDAC-THA  
           "#FF1300",   #VDAC-VNM  
           "#ff8880",   #VDAM-CHN  
           "#ff8880",   #VDAM-KOR  
           "#ff8880",   #VDAM-VNM 
           "darkred",   #VDMIS-VNM
           "#0040FF",   #VJAC-IDN  
           "#0040FF",   #VJAC-MYS  
           "#0040FF",   #VJAC-PNG  
           "#0040FF",   #VJAC-VNM  
           "skyblue1",  #VJAM-PNG 
           "green")     #VJMIS-PNG 

myCol.vd <- c("#FF1300",   #VDAC-CHN
           "#FF1300",   #VDAC-KOR
           "#FF1300",   #VDAC-THA
           "#FF1300",   #VDAC-VNM
           "#ff8880",   #VDAM-CHN
           "#ff8880",   #VDAM-KOR
           "#ff8880",   #VDAM-VNM
           "darkred")   #VDMIS-VNM

myCol.vj <- c("#0040FF",   #VJAC-IDN
           "#0040FF",   #VJAC-MYS
           "#0040FF",   #VJAC-PNG
           "#0040FF",   #VJAC-VNM
           "skyblue1",  #VJAM-PNG
           "green")     #VJMIS-PNG


## For China, 15 
## For Korea, 17
## For Thailand, 8
## For Viet Nam, 4
## For Indonesia, 9
## For Malaysia, 10
## For Papua New Guinea, 16
## For VJMIS "green"

myShape.all <- c(22,   #VDAC-CHN  
           24,   #VDAC-KOR  
           8,   #VDAC-THA  
           9,   #VDAC-VNM  
           22,   #VDAM-CHN  
           24,   #VDAM-KOR  
           9,   #VDAM-VNM 
           9,   #VDMIS-VNM
           10,   #VJAC-IDN  
           4,   #VJAC-MYS  
           21,   #VJAC-PNG  
           9,   #VJAC-VNM  
           21,  #VJAM-PNG 
           21)     #VJMIS-PNG

myShape.vd <- c(22,   #VDAC-CHN
           24,   #VDAC-KOR
           8,   #VDAC-THA
           9,   #VDAC-VNM
           22,   #VDAM-CHN
           24,   #VDAM-KOR
           9,   #VDAM-VNM
           9)   #VDMIS-VNM

myShape.vj <- c(10,   #VJAC-IDN
           4,   #VJAC-MYS
           21,   #VJAC-PNG
           9,   #VJAC-VNM
           21,  #VJAM-PNG
           21)     #VJMIS-PNG

## Get the variance for each axix
axe1 <- paste(round(pca.var$eig[1]/sum(pca.var$eig)*100, digits = 2))
axe1
axe2 <- paste(round(pca.var$eig[2]/sum(pca.var$eig)*100, digits = 2))
axe2

### DAPC
#dapc.varroa2 <- dapc(gl.varroa, n.pca = 4, n.da = 2)
#dapc.results2 <- as.data.frame(dapc.varroa2$posterior)

#dapc.varroa3 <- dapc(gl.varroa, n.pca = 4, n.da = 3)

#dapc.varroa4 <- dapc(gl.varroa, n.pca = 4, n.da = 4)
#dapc.varroa5 <- dapc(gl.varroa, n.pca = 4, n.da = 5)

cols <- brewer.pal(n = nPop(gl.rubi), name = "Dark2")

##export pdf
pdf(outfile, height=5, width= 7)

evoli <- barplot(100*pca.var$eig/sum(pca.var$eig), col = heat.colors(50), main="PCA Eigenvalues")
evoli <- evoli + title(ylab="Proportion of variance explained")
evoli <- evoli + title(xlab="Eigenvalue")
evoli

pikachu <- ggplot(varroa.pca.scores, aes(x=PC1, y=PC2, colour=pop)) 
pikachu <- pikachu + geom_point(aes(fill=pop, shape=pop), size=4, alpha = 0.8)
pikachu <- pikachu + scale_fill_manual(values = myCol.all) 
pikachu <- pikachu + scale_shape_manual(values = myShape.all)
pikachu <- pikachu + scale_color_manual(values = myCol.all)
pikachu <- pikachu + geom_hline(yintercept = 0) 
pikachu <- pikachu + geom_vline(xintercept = 0) 
pikachu <- pikachu + theme_bw()
pikachu <- pikachu + theme(plot.background = element_blank()
  ,panel.grid.major = element_blank()
  ,panel.grid.minor = element_blank())
pikachu <- pikachu + labs(title="PCA on 44 varroa samples with LD prune SNPs")
pikachu

#miaous2 <- scatter(dapc.varroa2, col = cols, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE, posi.pca = "topleft", cleg = 0.75)
#miaous2
#persian2 <- compoplot(dapc.varroa2, col = cols, posi = 'top')
#persian2

#miaous3 <- scatter(dapc.varroa3, col = cols, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE, posi.pca = "topleft", cleg = 0.75)
#miaous3
#persian3 <- compoplot(dapc.varroa3, col = cols, posi = 'top')
#persian3

#miaous4 <- scatter(dapc.varroa4, col = cols, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE, posi.pca = "topleft", cleg = 0.75)
#miaous4
#persian4 <- compoplot(dapc.varroa4, col = cols, posi = 'top')
#persian4

#miaous5 <- scatter(dapc.varroa5, col = cols, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE, posi.pca = "topleft", cleg = 0.75)
#miaous5
#persian5 <- compoplot(dapc.varroa5, col = cols, posi = 'top')
#persian5

dev.off()

