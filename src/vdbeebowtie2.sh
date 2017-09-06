#!/bin/bash
#SBATCH --job-name=vdvjbowtie2
#SBATCH --partition=compute
#SBATCH --mem=15G
#SBATCH --cpus-per-task=2
#SBATCH --time=7-0
##SBATCH --mail-user=%u@oist.jp
##SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --input=none
#SBATCH --output=%j.out
##SBATCH --error=job_%j.err
#SBATCH --array=0-57

. $HOME/.bashrc

a=(`ls -1 /work/MikheyevU/Maeva/varroahost/data/*-R1_001.fastq.gz`) #51
b=(`ls -1 /work/MikheyevU/Maeva/varroahost/data/*-R2_001.fastq.gz`)

base=$(basename ${a["SLURM_ARRAY_TASK_ID"]} "-R1_001.fastq.gz")
f=${a["SLURM_ARRAY_TASK_ID"]}
r=${b["SLURM_ARRAY_TASK_ID"]}
echo $f
echo $r

bowbase=/work/MikheyevU/Maeva/varroahost/ref/bee/Amel45
ref=/work/MikheyevU/Maeva/varroahost/ref/bee/Amel_4.5_genomic.fasta

outfile=/work/MikheyevU/Maeva/varroahost/bam/alignvsbee/`basename $f -R1_001.fastq.gz`BEE.bam

if [ -s "$outfile" ]
then
    echo "$outfile has some data."
    sleep 10s
else
    echo "$oufile is empty."
    bowtie2 -p 2 --very-sensitive-local --sam-rg ID:$base --sam-rg LB:Nextera --sam-rg SM:$base --sam-rg PL:ILLUMINA -x $bowbase -1 $f -2 $r | samtools view -Su -F4 - | novosort  -c 2 -m 10G -i -o /work/MikheyevU/Maeva/varroahost/bam/alignvsbee/`basename $f -R1_001.fastq.gz`BEE.bam -
fi

