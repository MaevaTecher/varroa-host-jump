# Population genetic analysis of Varroa on native and introudced hosts
outDir = "data"
refDir = "ref" 
SCRATCH  = "/work/scratch/maeva"
hostBeeBowtieIndex = refDir + "/bees/hostbee"
hostBeeMtBowtieIndex = refDir + "/bees/mtdna"
hostBowtiemellifera = refDir + "/bees/mellifera"
hostBowtiecerana = refDir + "/bees/cerana"
varroaBowtieIndex = refDir + "/destructor/vd"
vdRef = refDir + "/destructor/vd.fasta"

SAMPLES, = glob_wildcards(outDir + "/reads/{sample}-R1_001.fastq.gz")

rule all:
	input: expand(outDir + "/hostbee/mellifera/{sample}.bam", sample = SAMPLES),
		expand(outDir + "/hostbee/cerana/{sample}.bam", sample = SAMPLES)

# use all genome to verify host identity
rule checkmellifera:
	input:
		read1 = outDir + "/reads/{sample}-R1_001.fastq.gz",
		read2 = outDir + "/reads/{sample}-R2_001.fastq.gz"
	output: temp(outDir + "/hostbee/mellifera/{sample}.bam")
	shell: "bowtie2 -p {threads} -x {hostBowtiemellifera} -1 {input.read1} -2 {input.read2}  | samtools view -Su -F4 | novosort -c 2 -m 20G -i -o {output} -"
		
rule checkcerana:
	input:
		read1 = outDir + "/reads/{sample}-R1_001.fastq.gz",
		read2 = outDir + "/reads/{sample}-R2_001.fastq.gz"
	output: temp(outDir + "/hostbee/cerana/{sample}.bam")
	shell: "bowtie2 -p {threads} -x {hostBowtiecerana} -1 {input.read1} -2 {input.read2}  | samtools view -Su -F4 | novosort -c 2 -m 20G -i -o {output} -"
