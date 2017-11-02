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

VDREF=/work/MikheyevU/Maeva/varroa-jump/ref/destructor/vd.fasta
BAMLIST=/work/MikheyevU/Maeva/varroa-jump/data/angsd/ngmbam.list
OUTDIR=/work/MikheyevU/Maeva/varroa-jump/data/angsd

$ANGSD/angsd -P 12 -b $BAMLIST -ref $VDREF -out $OUTDIR/mitegeno \
        -uniqueOnly 1 \
        -remove_bads 1 \
        -only_proper_pairs 1\
	-trim 0 \
	-C 50\
	-baq 1\	
	-minMapQ 20
	-minQ 20
	-doCounts 1\
	-GL 1\
	-doMajorMinor 1\
	-doMaf 1\
	-skipTriallelic 1\
	-SNP_pval 1e-8
      	-doGeno 32\
        -doPost 1
