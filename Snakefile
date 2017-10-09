# Population genetic analysis of Varroa on native and introudced hosts
from scripts.split_fasta_regions import split_fasta

outDir = "data"
refDir = "ref" 
hostBeeBowtieIndex = refDir + "/bees/hostbee"
varroaBowtieIndex = refDir + "/destructor/vd"
vdRef = refDir + "/destructor/vd.fasta"

SPLITS = range(200)
REGIONS = split_fasta(vdRef, len(SPLITS))  # dictionary with regions to be called, with keys in SPLITS
for region in REGIONS:
	for idx,i in enumerate(REGIONS[region]):
		REGIONS[region][idx] = " -r " + str(i)


SAMPLES, = glob_wildcards(outDir + "/reads/{sample}-R1_001.fastq.gz")

rule all:
	input: outDir + "/sketches/varroa.dnd",
		expand(outDir + "/var/bowtie2/split/freebayes.{region}.vcf", region = REGIONS),
		expand(outDir + "/var/ngm/split/freebayes.{region}.vcf", region = REGIONS)
		
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
	output: tree = outDir + "/sketches/varroa.dnd", matrix = outDir + "/sketches/varroa.phylip"
	threads: 12
	shell: "mashtree.pl --genomesize 500000000 --mindepth 2 --tempdir /work/MikheyevU/Maeva/varroahost/scratch --numcpus {threads} --outmatrix {output.matrix} {input} > {output.tree}"

rule bowtie2:
	input:
		read1 = outDir + "/reads/{sample}-R1_001.fastq.gz",
		read2 = outDir + "/reads/{sample}-R2_001.fastq.gz",
	threads: 12
	output: 
		alignment = temp(outDir + "/alignments/bowtie2/{sample}.bam"), 
		index = temp(outDir + "/alignments/bowtie2/{sample}.bai")
	shell:
		"""
		module load bowtie2/2.2.6 samtools/1.3.1
		bowtie2 -p {threads} --very-sensitive-local --sam-rg ID:{wildcards.sample} --sam-rg LB:Nextera --sam-rg SM:{wildcards.sample} --sam-rg PL:ILLUMINA -x {varroaBowtieIndex} -1 {input.read1} -2 {input.read2} | samtools view -Su - | samtools sort - -m 55G -o - | samtools rmdup - {output.alignment}
		samtools index {output.alignment}
		"""

rule nextgenmap:
	input:
		read1 = outDir + "/reads/{sample}-R1_001.fastq.gz",
		read2 = outDir + "/reads/{sample}-R2_001.fastq.gz",
	threads: 12
	output: 
		alignment = temp(outDir + "/alignments/ngm/{sample}.bam"), 
		index = temp(outDir + "/alignments/ngm/{sample}.bai")
	shell:
		"""
		module load NextGenMap/0.5.0 samtools/1.3.1
		ngm -t {threads} -b  -1 {input.read1} -2 {input.read2} -r {vdRef} --rg-id {wildcards.sample} --rg-sm {wildcards.sample} --rg-pl ILLUMINA --rg-lb {wildcards.sample} | samtools sort - -m 55G -o - | samtools rmdup - {output.alignment}
		samtools index {output.alignment}
		"""
	
rule freeBayesNGM:
	input: 
		expand(outDir + "/alignments/ngm/{sample}.bam", sample = SAMPLES)
	output: 
		temp(outDir + "/var/ngm/split/freebayes.{region}.vcf")
	params: 
		span = lambda wildcards: REGIONS[wildcards.region],
		bams = lambda wildcards, input: os.path.dirname(input[0]) + "/*.bam"
	shell:
		"""
		module load freebayes/1.1.0 vcftools/0.1.12b
		freebayes --min-alternate-fraction 0.2 --use-best-n-alleles 3 -b {params.bams} {params.span}  -f {vdRef} | vcftools  --vcf - --minQ 40 --recode --stdout > {output}
		"""
rule freeBayesBowtie2:
	input: 
		expand(outDir + "/alignments/ngm/{sample}.bam", sample = SAMPLES)
	output: 
		temp(outDir + "/var/bowtie2/split/freebayes.{region}.vcf")
	params: 
		span = lambda wildcards: REGIONS[wildcards.region],
		bams = lambda wildcards, input: os.path.dirname(input[0]) + "/*.bam"
	shell:
		"""
		module load freebayes/1.1.0 vcftools/0.1.12b
		freebayes --min-alternate-fraction 0.2 --use-best-n-alleles 3 -b {params.bams} {params.span}  -f {vdRef} | vcftools  --vcf - --minQ 40 --recode --stdout > {output}
		"""
