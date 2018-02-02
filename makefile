rule ?= all		# run rule all by default, which is the default behavior anyway
restart ?= 0	# by default do not restart failed jobs
SHELL := /bin/bash	# allows the running of modules

dry:
	snakemake -p -n $(rule)
srun:
	module load R/3.4.2
	srun -t 7-0 snakemake --restart-times $(restart) -j 500 -p --max-jobs-per-second 1  --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} --cpus-per-task {threads} -t {resources.time} --mem {resources.mem}G" --rerun-incomplete --notemp --nolock $(rule) 
rmout: 
	rm -f *.out
