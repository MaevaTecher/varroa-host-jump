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
NGSADMIX=/apps/unit/MikheyevU/NGSadmix/32/NGSadmix

MTDNAREF=/work/MikheyevU/Maeva/varroa-jump/ref/destructor/mtdnamite/VDAJ493124.fasta
BAMLIST=/work/MikheyevU/Maeva/varroa-jump/data/mtdna_bowtie2/mtbam.list
INPUTBEAGLE=/work/MikheyevU/Maeva/varroa-jump/data/mtdna_var/mite.beagle.gz
OUTDIR=/work/MikheyevU/Maeva/varroa-jump/data/bayesian

K=2

$NGSADMIX -P 12 -likes $INPUTBEAGLE -K $K -outfiles $OUTDIR/assign_$K -minMaf 0
