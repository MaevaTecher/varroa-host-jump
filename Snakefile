# Population genetic analysis of Varroa on native and introudced hosts
 
hostBeeBowtieIndex = "ref/bees/hostbee"
SAMPLES, = glob_wildcards("data/reads/{sample}-R1_001.fastq.gz")

rule all:
	input:
		expand("data/reads/sketches/{sample}.msh", sample = SAMPLES)

rule varroaSketch:
	input:
		read1 = "data/reads/{sample}-R1_001.fastq.gz",
		read2 = "data/reads/{sample}-R2_001.fastq.gz",
	threads: 12
	output: "data/reads/sketches/{sample}.msh"
	shell: 
		"""
		mash sketch -p 4 -m 2 -o data/reads/sketches/{wildcards.sample} -r  <(bowtie2 -p {threads} -x {hostBeeBowtieIndex} -1  {input.read1} -2 {input.read2}  | samtools view -S -f12 | awk '{{print ">"$1"\\n"$10}}' )
		"""

