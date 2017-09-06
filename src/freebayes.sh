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

freebayes --fasta-reference /work/MikheyevU/Maeva/varroahost/ref/454LargeContigs.fasta \
--use-best-n-alleles 4 \
--bam /work/MikheyevU/Maeva/varroahost/freebayes/VD149_S1.bam \
/work/MikheyevU/Maeva/varroahost/freebayes/VD212_S6.bam \
/work/MikheyevU/Maeva/varroahost/freebayes/VD651_1_S18.bam \
/work/MikheyevU/Maeva/varroahost/freebayes/VD658_2_S21.bam \
/work/MikheyevU/Maeva/varroahost/freebayes/VD787_2_S3.bam \
/work/MikheyevU/Maeva/varroahost/freebayes/VD153_2_S3.bam \
/work/MikheyevU/Maeva/varroahost/freebayes/VD150_2_S6.bam \
/work/MikheyevU/Maeva/varroahost/freebayes/VD622_1_S12.bam \
/work/MikheyevU/Maeva/varroahost/freebayes/VD639_1_S14.bam \
/work/MikheyevU/Maeva/varroahost/freebayes/VD646_1_S17.bam \
/work/MikheyevU/Maeva/varroahost/freebayes/VJ028_S28.bam \
/work/MikheyevU/Maeva/varroahost/freebayes/VJ325_1_S30.bam \
/work/MikheyevU/Maeva/varroahost/freebayes/VJ347_1_S33.bam \
/work/MikheyevU/Maeva/varroahost/freebayes/VJ853_4_S9.bam \
/work/MikheyevU/Maeva/varroahost/freebayes/VJ983_1_S51.bam \
/work/MikheyevU/Maeva/varroahost/freebayes/VJ847_S38.bam  \
/work/MikheyevU/Maeva/varroahost/freebayes/VJ856_1_S42.bam \
/work/MikheyevU/Maeva/varroahost/freebayes/VJ952_1_S46.bam \
/work/MikheyevU/Maeva/varroahost/freebayes/VJ954_2_S48.bam \
/work/MikheyevU/Maeva/varroahost/freebayes/VJ956_4_S10.bam \
-v freebayes10.vcf

