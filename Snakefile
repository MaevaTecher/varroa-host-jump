# Population genetic analysis of Varroa on native and introudced hosts
outDir = "/work/MikheyevU/Maeva/varroahost/data"
refDir = "/work/MikheyevU/Maeva/varroahost/ref" 
hostBeeBowtieIndex = refDir + "/bees/hostbee"
varroaBowtieIndex = refDir + "/destructor/vdjellytrim"
VDREF = refDir + "/destructor/vd_assembly_filled.trimmed.fasta"
jacobBowtieIndex = refDir + "/jacobsoni/vjassembly"
VJREF = refDir + "/jacobsoni/vj_454LargeContigs.fna"
SAMPLES, = glob_wildcards(outDir + "/reads/{sample}-R1_001.fastq.gz")

rule all:
	outDir + "/picardstat/{sample}_picardVD.txt", outDir + "/sketches/variant_destructor.vcf"
		
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

rule map2destructor:
	input:
		read1 = outDir + "/reads/{sample}-R1_001.fastq.gz",
		read2 = outDir + "/reads/{sample}-R2_001.fastq.gz",
	threads: 12
	output: temp(outDir + "/mapbam/{sample}.bam")
	shell: "bowtie2 -p {threads} -x {varroaBowtieIndex} -1 {input.read1} -2 {input.read2}  | samtools view -Su -F4 | novosort -c 2 -m 20G -i -o {output} -"	

rule destructorstat:
	input: expand(outDir + "/mapbam/{sample}.bam", sample = SAMPLES)
	output: temp(outDir + "/picardstat/{sample}_picardVD.txt")
	shell:  "picard CollectAlignmentSummaryMetrics R= {VDREF} I= {input} O= {output}"
		
rule freebayes_destructor:
	input: expand(outDir + "/mapbam/{sample}.bam", sample = SAMPLES)
	output: outDir + "/sketches/variant_destructor.vcf"
	shell:  "freebayes --use-best-n-alleles 4 --bam {input} -v {output} -f {VDREF}"

rule map2jacobsoni:
	input:
		read1 = outDir + "/reads/{sample}-R1_001.fastq.gz",
		read2 = outDir + "/reads/{sample}-R2_001.fastq.gz",
	threads: 12
	output: temp(outDir + "/map_vj/{sample}.bam")
	shell: "bowtie2 -p {threads} -x {jacobBowtieIndex} -1 {input.read1} -2 {input.read2}  | samtools view -Su -F4 | novosort -c 2 -m 20G -i -o {output} -"	
		
rule freebayes_jacobsoni:
	input: expand(outDir + "/map_vj/{sample}.bam", sample = SAMPLES)
	output: outDir + "/sketches/variant_jacobsoni.vcf"
	shell:  "freebayes --use-best-n-alleles 4 --bam {input} -v {output} -f {VJREF}"
