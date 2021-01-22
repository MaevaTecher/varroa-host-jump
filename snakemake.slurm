#!/bin/bash
#SBATCH --job-name=snakemake
#SBATCH --partition=short
#SBATCH --time="0-02:00:00"
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
##SBATCH --mail-user=%u@oist.jp
##SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --input=none
#SBATCH --output=%j.out
##SBATCH --error=job_%j.err

#module load python/3.5.2
#conda activate popgenomics

snakemake --restart-times 1 -j 500 -p --max-jobs-per-second 1 --cluster-config cluster.json --cluster "sbatch  -p {cluster.partition} --cpus-per-task {cluster.cpus-per-task} -t {cluster.time} --mem {cluster.mem}" --rerun-incomplete --notemp --nolock 

