#!/bin/bash
#SBATCH --job-name=gatk
#SBATCH --partition=compute
#SBATCH --mem=30G
#SBATCH --cpus-per-task=1
#SBATCH --time=1-0:00:00
#SBATCH --ntasks=2
##SBATCH --mail-user=%u@oist.jp
##SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --input=none
#SBATCH --output=%j.out
##SBATCH --error=job_%j.err

. $HOME/.bashrc 

reference=/work/MikheyevU/Maeva/varroa-jump/ref/destructor/vd.fasta
varroavcf=/work/MikheyevU/Maeva/varroa-jump/scripts/VJ333_1.vcf
indiv="VJ333_1"

java -jar /apps/unit/MikheyevU/Maeva/GATK/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R $reference -L "BEIS01000001.1:10000-20000" -V $varroavcf -o $indiv.fasta -IUPAC $indiv
