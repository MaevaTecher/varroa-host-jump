#!/bin/bash
#SBATCH --job-name=vdvjbowtie2
#SBATCH --partition=compute
#SBATCH --mem=15G
#SBATCH --cpus-per-task=2
#SBATCH --time=7-0
#SBATCH --input=none
#SBATCH --output=%j.out
##SBATCH --array=0-57

. $HOME/.bashrc

##Specify the path for the apps I need to use
#NGSTOOLS=/apps/unit/MikheyevU/Maeva/ngsTools
ANGSD=/apps/unit/MikheyevU/Maeva/angsd
#NGSADMIX=/apps/unit/MikheyevU/NGSadmix/32/NGSadmix

#SAMTOOLS=/apps/unit/MikheyevU/Maeva/samtools/bin/samtools
#HTSLIB=/apps/unit/MikheyevU/Maeva/htslib

MTDNAREF=/work/MikheyevU/Maeva/varroa-jump/ref/destructor/mtdnamite/VDAJ493124.fasta
BAMLIST=/work/MikheyevU/Maeva/varroa-jump/data/mtdna_bowtie2/mtbam.list
OUTDIR=/work/MikheyevU/Maeva/varroa-jump/data/mtdna_var

$ANGSD/angsd -P 12 -b $BAMLIST -ref $MTDNAREF -out $OUTDIR/mite \
        -uniqueOnly 1 \
        -remove_bads 1 \
        -only_proper_pairs 1\
        -SNP_pval 1e-8\
        -minMapQ 20 \
        -minQ 20 \
	-doGlf 2 \
	-doGeno 32\
        -trim 0 \
        -C 50 \
        -baq 1 \
	#-minInd 15 \
        #-setMinDepth 60 \
        #-setMaxDepth 400 \
        -doCounts 1 \
        -GL 1 \
        -doMajorMinor 4 \
        -doMaf 1 \
        -skipTriallelic 1\
	-doPost 1


