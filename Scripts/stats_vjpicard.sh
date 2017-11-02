#!/bin/bash
#SBATCH --job-name=statsVJ
#SBATCH --partition=compute
#SBATCH --mem=30G
#SBATCH --cpus-per-task=1
#SBATCH --time=1-0:00:00
#SBATCH --ntasks=1
##SBATCH --mail-user=%u@oist.jp
##SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --input=none
#SBATCH --output=%j.out
##SBATCH --error=job_%j.err
#SBATCH --array=0-43

. $HOME/.bashrc 

a=(`for i in /work/MikheyevU/Maeva/varroa-jump/data/map_vj/*.bam; do basename $i; done | sort -u`) # 44

f=${a[$SLURM_ARRAY_TASK_ID]}
reference=/work/MikheyevU/Maeva/varroa-jump/ref/jacobsoni/vj_454LargeContigs.fna

java -jar /apps/unit/MikheyevU/Maeva/picard/picard.jar CollectAlignmentSummaryMetrics R= $reference I=/work/MikheyevU/Maeva/varroa-jump/data/map_vj/$f O=/work/MikheyevU/Maeva/varroa-jump/data/picardstat/`basename $f .bam`_vjreport.txt

