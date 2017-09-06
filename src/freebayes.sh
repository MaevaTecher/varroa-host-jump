#!/bin/bash
#SBATCH --job-name=freebayes
#SBATCH --partition=compute
#SBATCH --mem=20G
#SBATCH --cpus-per-task=2
#SBATCH --time=7-0
#SBATCH --ntasks=1
##SBATCH --mail-user=%u@oist.jp
##SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --input=none
#SBATCH --output=%j.out
##SBATCH --error=job_%j.err

. $HOME/.bashrc 

module load freebayes

freebayes --fasta-reference 454LargeContigs.fasta --bam VD149_S1.bam VD646_1_S17.bam VJ325_1_S30.bam VJ950_1_S44.bam -v freebayes4.vcf

