# Population genetic analysis of Varroa on native and introudced hosts
outDir = "/work/MikheyevU/Maeva/varroahost/data"
refDir = "/work/MikheyevU/Maeva/varroahost/ref" 
hostBeeBowtieIndex = refDir + "/bees/hostbee"
SAMPLES, = glob_wildcards(outDir + "/reads/{sample}-R1_001.fastq.gz")

rule all:
	input: outDir + "/sketches/varroa.dnd"
		
rule removeHost:
	input:
		read1 = outDir + "/reads/{sample}-R1_001.fastq.gz",
		read2 = outDir + "/reads/{sample}-R2_001.fastq.gz",
	threads: 12
	output: temp(outDir + "/sketches/{sample}.fastq.gz")
	shell: 
		"""
		bowtie2 -p {threads} -x {hostBeeBowtieIndex} -1  {input.read1} -2 {input.read2}  | samtools view -S -f12 | awk '{{print "@"$1"\\n"$10"\\n+\\n"$11}}' | gzip > {output}
		"""

rule mashtree:
	input: expand(outDir + "/sketches/{sample}.fastq.gz", sample = SAMPLES)
	output: outDir + "/sketches/varroa.dnd"
	threads: 12
	shell: "mashtree.pl --genomesize 500000000 --mindepth 2 --tmpdir /tmp/sasha --numcpus {threads} --outmatrix {input} > {output}"
