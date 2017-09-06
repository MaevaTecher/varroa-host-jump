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
--bam /work/MikheyevU/Maeva/varroahost/bam/alignvd454Large/VD149_S1.bam \
/work/MikheyevU/Maeva/varroahost/bam/alignvd454Large/VD212_S6.bam \
/work/MikheyevU/Maeva/varroahost/bam/alignvd454Large/VD651_1_S18.bam \
/work/MikheyevU/Maeva/varroahost/bam/alignvd454Large/VD658_2_S21.bam \
/work/MikheyevU/Maeva/varroahost/bam/alignvd454Large/VD787_2_S3.bam \
/work/MikheyevU/Maeva/varroahost/bam/alignvd454Large/VD153_2_S3.bam \
/work/MikheyevU/Maeva/varroahost/bam/alignvd454Large/VD150_2_S6.bam \
/work/MikheyevU/Maeva/varroahost/bam/alignvd454Large/VD622_1_S12.bam \
/work/MikheyevU/Maeva/varroahost/bam/alignvd454Large/VD639_1_S14.bam \
/work/MikheyevU/Maeva/varroahost/bam/alignvd454Large/VD646_1_S17.bam \
/work/MikheyevU/Maeva/varroahost/bam/alignvd454Large/VJ028_S28.bam \
/work/MikheyevU/Maeva/varroahost/bam/alignvd454Large/VJ325_1_S30.bam \
/work/MikheyevU/Maeva/varroahost/bam/alignvd454Large/VJ347_1_S33.bam \
/work/MikheyevU/Maeva/varroahost/bam/alignvd454Large/VJ853_4_S9.bam \
/work/MikheyevU/Maeva/varroahost/bam/alignvd454Large/VJ983_1_S51.bam \
/work/MikheyevU/Maeva/varroahost/bam/alignvd454Large/VJ847_S38.bam  \
/work/MikheyevU/Maeva/varroahost/bam/alignvd454Large/VJ856_1_S42.bam \
/work/MikheyevU/Maeva/varroahost/bam/alignvd454Large/VJ952_1_S46.bam \
/work/MikheyevU/Maeva/varroahost/bam/alignvd454Large/VJ954_2_S48.bam \
/work/MikheyevU/Maeva/varroahost/bam/alignvd454Large/VJ956_4_S10.bam \
-v /work/MikheyevU/Maeva/varroahost/variant/freebayes10.vcf

