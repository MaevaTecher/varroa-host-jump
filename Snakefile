# Population genetic analysis of Varroa on native and introduced hosts
from scripts.split_fasta_regions import split_fasta

outDir = "data"
refDir = "ref" 
SCRATCH  = "/work/scratch/maeva"
##Bee references
hostBeeBowtieIndex = refDir + "/bees/hostbee"
hostBeeMtBowtieIndex = refDir + "/bees/mtdna"
hostBowtiemellifera = refDir + "/bees/mellifera"
hostBowtiecerana = refDir + "/bees/cerana"

#Varroa destructor references
varroaBowtieIndex = refDir + "/destructor/vd"
vdRef = refDir + "/destructor/vd.fasta"
vdmtDNABowtieIndex = refDir + "/destructor/mtdnamite/vdnavajas"
vdmtDNA = refDir + "/destructor/mtdnamite/VDAJ493124.fasta"

SAMPLES, = glob_wildcards(outDir + "/reads/{sample}-R1_001.fastq.gz")

SPLITS = range(10)
REGIONS = split_fasta(vdmtDNA, len(SPLITS))  # dictionary with regions to be called, with keys in SPLITS
Q = (20, 40) # 99 and 99.99% mapping accuracy
for region in REGIONS:
	for idx,i in enumerate(REGIONS[region]):
		REGIONS[region][idx] = " -r " + str(i)


rule all:
	input: expand(outDir + "/mtdna_bowtie2/{sample}.bam", sample = SAMPLES),
		expand(outDir + "mtdna_bowtie2/{sample}.bam.bai", sample = SAMPLES),
		"data/mtdna_var/mtdnafiltered.vcf"

# use all genome to verify host identity
#rule checkmellifera:
#	input:
#		read1 = outDir + "/reads/{sample}-R1_001.fastq.gz",
#		read2 = outDir + "/reads/{sample}-R2_001.fastq.gz"
#	output: temp(outDir + "/hostbee/mellifera/{sample}.bam")
#	shell: "bowtie2 -p {threads} -x {hostBowtiemellifera} -1 {input.read1} -2 {input.read2}  | samtools view -Su -F4 | novosort -c 2 -m 20G -i -o {output} -"
		
#rule checkcerana:
#	input:
#		read1 = outDir + "/reads/{sample}-R1_001.fastq.gz",
#		read2 = outDir + "/reads/{sample}-R2_001.fastq.gz"
#	output: temp(outDir + "/hostbee/cerana/{sample}.bam")
#	shell: "bowtie2 -p {threads} -x {hostBowtiecerana} -1 {input.read1} -2 {input.read2}  | samtools view -Su -F4 | novosort -c 2 -m 20G -i -o {output} -"

rule bowtie2mtdna:
	input:
		read1 = outDir + "/reads/{sample}-R1_001.fastq.gz",
		read2 = outDir + "/reads/{sample}-R2_001.fastq.gz",
	threads: 12
	output: 
		alignment = temp(outDir + "/mtdna_bowtie2/{sample}.bam"), 
		index = temp(outDir + "/mtdna_bowtie2/{sample}.bam.bai")
	shell:
		"""
		module load bowtie2/2.2.6 samtools/1.3.1 VariantBam/1.4.3
		bowtie2 -p {threads} --very-sensitive-local --sam-rg ID:{wildcards.sample} --sam-rg LB:Nextera --sam-rg SM:{wildcards.sample} --sam-rg PL:ILLUMINA -x {vdmtDNABowtieIndex} -1 {input.read1} -2 {input.read2} | samtools view -Su - | samtools sort - -m 55G -T {SCRATCH}/{wildcards.sample} -o - | samtools rmdup - - | variant - -m 500 -b -o {output.alignment}
		samtools index {output.alignment}
		"""

rule freeBayesmtdna:
	input: 
		expand(outDir + outDir + "/mtdna_bowtie2/{sample}.bam", sample = SAMPLES)
	output: 
		temp(outDir + "/mtdna_var/mtdna.{region}.vcf")
	params: 
		span = lambda wildcards: REGIONS[wildcards.region],
		bams = lambda wildcards, input: os.path.dirname(input[0]) + "/*.bam",
		missing = lambda wildcards, input: len(input) * 0.9
	shell:
		"""
		module load freebayes/1.1.0 vcftools/0.1.12b vcflib/1.0.0-rc1
		for i in {params.bams}; do name=$(basename $i .bam); if [[ $name == VJ* ]] ; then echo $name VJ; else echo $name VD; fi ; done > {outDir}/mtdna_var/pops.txt
		freebayes --min-alternate-fraction 0.2 --use-best-n-alleles 4 -m 5 -q 5 --populations {outDir}/var/pops.txt -b {params.bams} {params.span}  -f {vdmtDNA} | vcffilter  -f "QUAL > 20 & NS > {params.missing}" > {output}
		"""

rule mergeVCFmtdna:
	input: 
		expand(outDir + "/mtdna_var/mtdna.{region}.vcf", region = REGIONS)
	output:
		temp(outDir + "/mtdna_var/mtdnaraw.vcf")
	shell:
		"""
		module load vcflib/1.0.0-rc1
		(grep "^#" {input[0]} ; cat {input} | grep -v "^#" ) | vcfuniq  > {output}
		"""

rule filterVCFmtdna:
	# see http://ddocent.com//filtering/
	# filter DP  + 3*sqrt(DP) https://arxiv.org/pdf/1404.0929.pdf
	# also sites with more than two variants
	input:
		rules.mergeVCFmtdna.output
	output:
		vcf = outDir + "/mtdna_var/mtdnafiltered.vcf"
	shell:
		"""
		module load vcftools/0.1.12b vcflib/1.0.0-rc1 eplot/2.08
		perl  -ne 'print "$1\\n" if /DP=(\d+)/'  {input} > {outDir}/mtdna_var/depth.txt
		sort -n {outDir}/mtdna_var/depth.txt | uniq -c | awk '{{print $2,$1}}' | eplot -d -r [200:2000] 2>/dev/null | tee 
		Nind=$(grep -m1 "^#C" {input}  | cut -f10- |wc -w)
		coverageCutoff=$(awk -v Nind=$Nind '{{sum+=$1}} END {{print "DP < "(sum / NR / Nind + sqrt(sum / NR / Nind) * 3 ) * Nind}}' {outDir}/mtdna_var/depth.txt)
		echo Using coverage cutoff $coverageCutoff
		vcffilter -s -f \"$coverageCutoff\" {input} | vcftools --vcf - --max-alleles 2 --recode --stdout > {output}
		"""
	
