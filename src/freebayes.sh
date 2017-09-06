#!/bin/bash
#SBATCH --job-name=freebayes
#SBATCH --partition=compute
#SBATCH --mem=80G
#SBATCH --cpus-per-task=1
#SBATCH --time=7-0
#SBATCH --ntasks=1
##SBATCH --mail-user=%u@oist.jp
##SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --input=none
#SBATCH --output=%j.out
##SBATCH --error=job_%j.err

. $HOME/.bashrc 

module load freebayes

freebayes --fasta-reference /work/MikheyevU/Maeva/varroahost/ref454/LargeContigs.fasta --use-best-n-alleles 4 --bam /work/MikheyevU/Maeva/varroahost/freebayes/VD149_S1.bam /work/MikheyevU/Maeva/varroahost/freebayes/VD646_1_S17.bam /work/MikheyevU/Maeva/varroahost/freebayes/VJ325_1_S30.bam /work/MikheyevU/Maeva/varroahost/freebayes/VJ950_1_S44.bam -v freebayes4.vcf

