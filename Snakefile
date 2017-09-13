# Population genetic analysis of Varroa on native and introudced hosts
 
hostBeeBowtieIndex = "ref/bees/hostbee"
SAMPLES, = glob_wildcards("data/reads/{sample}-R1_001.fastq.gz")

rule all:
	input: "data/sketches/varroa.dnd"
		
rule removeHost:
	input:
		read1 = "data/reads/{sample}-R1_001.fastq.gz",
		read2 = "data/reads/{sample}-R2_001.fastq.gz",
	threads: 12
	output: temp("data/sketches/{sample}.fastq.gz")
	shell: 
		"""
		bowtie2 -p {threads} -x {hostBeeBowtieIndex} -1  {input.read1} -2 {input.read2}  | samtools view -S -f12 | awk '{{print "@"$1"\\n"$10"\\n+\\n"$11}}' | gzip > {output}
		"""

rule mashtree:
	input: expand("data/sketches/{sample}.fastq.gz", sample = SAMPLES)
	output: "data/sketches/varroa.dnd"
	threads: 12
	shell: "mashtree.pl --genomesize 500000000 --mindepth 2 --tmpdir /tmp/sasha --numcpus {threads} --outmatrix data/sketches/*fastq.gz > {output}"
