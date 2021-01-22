#!/bin/bash
#SBATCH --job-name=vdM1
#SBATCH --partition=compute
#SBATCH --mem=20G
#SBATCH -c 10
#SBATCH --time=4:00:00

number=$SLURM_ARRAY_TASK_ID

. $HOME/.bashrc 

/flash/MikheyevU/Maeva/varroa-jump-asia/tools/fsc26_linux64/fsc26 --tplfile div_nobot_"$number".tpl --estfile div_nobot_"$number".est --msfs --numsims 1000000 --maxlhood 0.001 --minnumloops 20 --numloops 150 -c 10 

