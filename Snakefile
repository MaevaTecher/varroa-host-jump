# Population genetic analysis of Varroa on native and introudced hosts
from scripts.split_fasta_regions import split_fasta

outDir = "data"
refDir = "ref" 
SCRATCH  = "/work/scratch/sasha"
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
	input: 
		"data/sketches/varroa.dnd",
		"data/var/filtered.vcf.gz"
		
rule removeHost:
	input:
		read1 = outDir + "/reads/{sample}-R1_001.fastq.gz",
		read2 = outDir + "/reads/{sample}-R2_001.fastq.gz",
	threads: 12
	output: temp(outDir + "/sketches/{sample}.fastq.gz")
	shell: 
		"""
		module load bowtie2/2.2.6 samtools/1.3.1
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
		index = temp(outDir + "/alignments/bowtie2/{sample}.bam.bai")
	shell:
		"""
		module load bowtie2/2.2.6 samtools/1.3.1 VariantBam/1.4.3
		bowtie2 -p {threads} --very-sensitive-local --sam-rg ID:{wildcards.sample} --sam-rg LB:Nextera --sam-rg SM:{wildcards.sample} --sam-rg PL:ILLUMINA -x {varroaBowtieIndex} -1 {input.read1} -2 {input.read2} | samtools view -Su - | samtools sort - -m 55G -T {SCRATCH}/bowtie/{wildcards.sample} -o - | samtools rmdup - - | variant - -m 500 -b -o {output.alignment}
		samtools index {output.alignment}
		"""

rule nextgenmap:
	input:
		read1 = outDir + "/reads/{sample}-R1_001.fastq.gz",
		read2 = outDir + "/reads/{sample}-R2_001.fastq.gz",
	threads: 12
	output: 
		alignment = temp(outDir + "/alignments/ngm/{sample}.bam"), 
		index = temp(outDir + "/alignments/ngm/{sample}.bam.bai")
	shell:
		"""
		module load NextGenMap/0.5.0 samtools/1.3.1 VariantBam/1.4.3
		ngm -t {threads} -b  -1 {input.read1} -2 {input.read2} -r {vdRef} --rg-id {wildcards.sample} --rg-sm {wildcards.sample} --rg-pl ILLUMINA --rg-lb {wildcards.sample} | samtools sort - -m 55G -T {SCRATCH}/ngm/{wildcards.sample} -o - | samtools rmdup - - | variant - -m 500 -b -o {output.alignment}
		samtools index {output.alignment}
		"""
	
rule freeBayes:
	input: 
		expand(outDir + "/alignments/{{aligner}}/{sample}.bam", sample = SAMPLES)
	output: 
		temp(outDir + "/var/{aligner}/split/freebayes.{region}.vcf")
	params: 
		span = lambda wildcards: REGIONS[wildcards.region],
		bams = lambda wildcards, input: os.path.dirname(input[0]) + "/*.bam",
		missing = lambda wildcards, input: len(input) * 0.9
	shell:
		"""
		module load freebayes/1.1.0 vcftools/0.1.12b vcflib/1.0.0-rc1
		for i in {params.bams}; do name=$(basename $i .bam); if [[ $name == VJ* ]] ; then echo $name VJ; else echo $name VD; fi ; done > {outDir}/var/pops.txt
		freebayes --min-alternate-fraction 0.2 --use-best-n-alleles 4 -m 5 -q 5 --populations {outDir}/var/pops.txt -b {params.bams} {params.span}  -f {vdRef} | vcffilter  -f "QUAL > 20 & NS > {params.missing}" > {output}
		"""

rule mergeVCF:
	input: 
		expand(outDir + "/var/{{aligner}}/split/freebayes.{region}.vcf", region = REGIONS)
	output:
		temp(outDir + "/var/{aligner}/raw.vcf")
	shell:
		"""
		module load vcflib/1.0.0-rc1
		(grep "^#" {input[0]} ; cat {input} | grep -v "^#" ) | vcfuniq  > {output}
		"""

rule filterVCF:
	# see http://ddocent.com//filtering/
	# filter DP  + 3*sqrt(DP) https://arxiv.org/pdf/1404.0929.pdf
	# also removing indels and sites with more than two variants
	input:
		rules.mergeVCF.output
	output:
		vcf = outDir + "/var/{aligner}/filtered.vcf"
	shell:
		"""
		module load vcftools/0.1.12b vcflib/1.0.0-rc1 eplot/2.08
		perl  -ne 'print "$1\\n" if /DP=(\d+)/'  {input} > {outDir}/var/{wildcards.aligner}/depth.txt
		sort -n {outDir}/var/{wildcards.aligner}/depth.txt | uniq -c | awk '{{print $2,$1}}' | eplot -d -r [200:2000] 2>/dev/null | tee 
		Nind=$(grep -m1 "^#C" {input}  | cut -f10- |wc -w)
		coverageCutoff=$(awk -v Nind=$Nind '{{sum+=$1}} END {{print "DP < "(sum / NR / Nind + sqrt(sum / NR / Nind) * 3 ) * Nind}}' {outDir}/var/{wildcards.aligner}/depth.txt)
		echo Using coverage cutoff $coverageCutoff
		vcffilter -s -f \"$coverageCutoff\" {input} | vcftools --vcf - --exclude-bed ref/destructor/masking_coordinates --max-alleles 2 --remove-indels --recode --stdout > {output}
		"""

rule chooseMapper:
	input:
		ngm = outDir + "/var/ngm/filtered.vcf", 
		bowtie2 = outDir + "/var/bowtie2/filtered.vcf", 
	output:
		bgzip = outDir + "/var/filtered.vcf.gz",
		primitives = outDir + "/var/primitives.vcf.gz"
	shell:
		"""
		module load samtools/1.3.1 vcflib/1.0.0-rc1
		ngm=$(grep -vc "^#" {input.ngm})
		bowtie2=$(grep -vc "^#" {input.bowtie2})
		echo ngm has $ngm snps vs $bowtie2 for botwie2 
		if [[ $ngm -gt $bowtie2 ]]; then
			echo choosing ngm
			bgzip -c {input.ngm} > {output.bgzip}
			vcfallelicprimitives {input.ngm} | bgzip > {output.primitives}
		else
			echo choosing bowtie2
			bgzip -c {input.bowtie2} > {output.bgzip}
			vcfallelicprimitives {input.bowtie2} | bgzip > {output.primitives}
		fi
		tabix -p vcf {output.bgzip} && tabix -p vcf {output.primitives}
		"""

rule wcFst:
	input:
	shell:
	"""
	

	dc=$(awk -v host=cerana -v species=VD -v ORS=, '(NR==FNR) {a[$1]=NR; next} ($2==host) && ($3==species) {print a[$1]}'  <(zcat  data/var/filtered.vcf.gz |grep -m1 "^#C" | cut -f10- | tr "\t" "\n") Rdata/varroa.txt | sed 's/,$//')
	dm=$(awk -v host=mellifera -v species=VD -v ORS=, '(NR==FNR) {a[$1]=NR; next} ($2==host) && ($3==species) {print a[$1]}'  <(zcat  data/var/filtered.vcf.gz |grep -m1 "^#C" | cut -f10- | tr "\t" "\n") Rdata/varroa.txt | sed 's/,$//')
	jm=$(awk -v host=mellifera -v species=VJ -v ORS=, '(NR==FNR) {a[$1]=NR; next} ($2==host) && ($3==species) {print a[$1]}'  <(zcat  data/var/filtered.vcf.gz |grep -m1 "^#C" | cut -f10- | tr "\t" "\n") Rdata/varroa.txt | sed 's/,$//')
	jc=$(awk -v host=cerana -v species=VJ -v ORS=, '(NR==FNR) {a[$1]=NR; next} ($2==host) && ($3==species) {print a[$1]}'  <(zcat  data/var/filtered.vcf.gz |grep -m1 "^#C" | cut -f10- | tr "\t" "\n") Rdata/varroa.txt | sed 's/,$//')
	"""


# --background <(awk -v host=cerana -v species=VD '($2==host) && ($3==species) {a[$1]==""}  ($1 in a) {out=out","NR} END {out} ' Rdata/varroa.txt | tr "\n" ",") --target <(awk -v host=mellifera -v species=VD '($2==host) && ($3==species) {a[$1]==""}  ($1 in a) {print NR} ' Rdata/varroa.txt | tr "\n" ",")

#vcftools --gzvcf data/var/filtered.vcf.gz --weir-fst-pop <(awk '($2=="mellifera") && ($3=="VD") {print $1}' Rdata/varroa.txt) --weir-fst-pop <(awk '($2=="cerana") && ($3=="VD") {print $1}' Rdata/varroa.txt) --fst-window-size 1000 --fst-window-step 500

# At this point we go with ngm, which produces a bit more variants


# # estimate SNP effects
# rule snpEff:
# 	input: rules.consensusFilter.output
# 	output: "../data/popgen/var/snpEff.txt"
# 	shell: "java -Xmx7g -jar /apps/unit/MikheyevU/sasha/snpEff4/snpEff.jar -no-utr -no-upstream -no-intron -no-intergenic -no-downstream pmuc {input} >  {output}"
# 	""" python parse_silentReplacement.py ../ref/csd.fa temp.txt > {output} && rm temp.txt """

# rule getCDS:
# 	input: GFF, REF
# 	output: "../ref/cds.fa"
# 	shell: "gffread {input[0]} -g {input[1]} -x {output}"

# rule filterLongest:
# 	input: rules.getCDS.output
# 	output: "../ref/longest.fa"
# 	shell: "python filter_longest.py {input} > {output}"

# # determine which SNPs are fixed and which are polymorphic
# # for this we remove the outgroup and compute frequencies
# rule fixedPolymorphic:	
# 	input: rules.consensusFilter.output
# 	output: "../data/popgen/var/snps.csv"
# 	shell: """module load zlib; vcftools --vcf {input} --remove-indv Pflavoviridis --freq; \
#     awk -v OFS="," ' NR>1 {{split($5,a,":"); if((a[2]=="1") || (a[2]=="0")) state="F"; else state="P"; print $1,$2,state}}' out.frq > {output} """

# # exports silent and replacement sites from snpEff
# rule parseSilentReplacement:
# 	input: rules.filterLongest.output, rules.snpEff.output
# 	output: "../data/popgen/var/annotation.csv"
# 	shell: ". ~/python2/bin/activate ; python parse_silentReplacement.py {input} > {output}"

# # calculate how many synonymous vs_non-synonymous changes are possible
# rule silentReplacement:
# 	input: rules.filterLongest.output
# 	output: "../data/popgen/var/silentReplacement.csv"
# 	shell: ". ~/python2/bin/activate; python silent_replacement.py {input} > {output}"

# rule snipre:
# 	input: rules.silentReplacement.output, rules.fixedPolymorphic.output, rules.parseSilentReplacement.output
# 	output: "../out/bayesian_results.csv"