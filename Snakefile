##############################################################################
##### POPULATION GENOMICS AND DEMOGRAPHIC INFERENCES SNAKEMAKE PIPELINE  #####
##### HOST SWITCHES IN APIS HONEY BEES BY PARASITIC VARROA MITES 	 #####
##### Ecology & Evolution OIST, Maeva TECHER and Alexander MIKHEYEV      #####
##############################################################################

from scripts.split_fasta_regions import split_fasta
from snakemake.utils import R
import getpass

localrules: getHaps, all

### SET DIRECTORY PATHS FOR REFERENCE AND OUTPUT DATA
READDir = "/bucket/MikheyevU/Maeva_data/world-varroa/data/reads"
BAM19Dir = "/flash/MikheyevU/Maeva/varroa-jump-asia/data_2019/alignments"
BAM20Dir = "/flash/MikheyevU/Maeva/varroa-jump-asia/data_2020/alignments"
OUT19Dir = "/flash/MikheyevU/Maeva/varroa-jump-asia/data_2019"
OUT20Dir = "/flash/MikheyevU/Maeva/varroa-jump-asia/data_2020"
REFDir = "/flash/MikheyevU/Maeva/varroa-jump-asia/ref2020"
SCRATCH  = "/flash/MikheyevU/Maeva/scratch" 
#easySFS = "/flash/MikheyevU/Maeva/demography_2020/tools/easySFS/easySFS.py"

KRAKENSTD = REFDir + "/metagenomics/minikraken_8GB_20200312"
KRAKENVK = REFDir + "/metagenomics/custom_vienna"

### PATHS FOR VARROA DESTRUCTOR GENOME AND REGIONS SPLIT
VDESRef = REFDir + "/destructor/vdes_3_refseq.fasta"
VDESBowtieIndex = REFDir + "/destructor/vdes_3_refseq"
VDESmtDNA = REFDir + "/destructor/mtDNA/NC_004454.fasta"

SPLITS = range(400)
REGIONS = split_fasta(VDESRef, len(SPLITS))  # dictionary with regions to be called, with keys in SPLITS
for region in REGIONS:
        for idx,i in enumerate(REGIONS[region]):
                REGIONS[region][idx] = " -r " + str(i)

SPLITSMT = range(20)
REGIONSMT = split_fasta(VDESmtDNA, len(SPLITSMT))
for regionmt in REGIONSMT:
        for idx,i in enumerate(REGIONSMT[regionmt]):
                REGIONSMT[regionmt][idx] = " --regions " + str(i)


### PATHS FOR HONEY BEE HOST GENOMES
hostBeeMtBowtieIndex = REFDir + "/bees/mtdna/hostbeemito"
hostBeeBowtieIndex = REFDir + "/bees/hostbee"

### SAMPLES LIST AND OTHER PARAMETERS
SAMPLES, = glob_wildcards(READDir + "/{sample}_R1_001.fastq.gz")
SAMPLES19, = glob_wildcards(BAM19Dir + "/ngm/{sample}.bam")
SAMPLES20, = glob_wildcards(BAM20Dir + "/ngm/{sample}.bam")

Q = (20, 40) # 99 and 99.99% mapping accuracy
CHROMOSOMES = ["NW_019211454.1", "NW_019211455.1", "NW_019211456.1", "NW_019211457.1", "NW_019211458.1", "NW_019211459.1", "NW_019211460.1"]
KCLUSTERS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
RUNS = ["run1", "run2", "run3", "run4", "run5", "run6", "run7", "run8", "run9", "run10"]
ALIGNERS = ["ngm", "bowtie2"]


#############################################################
##### TARGET PSEUDO-RULES DEFINITION                   ######
#############################################################

rule all:
	input: 	#expand(OUT20Dir + "/alignments/ngm/{sample}.bam", sample = SAMPLES20),
		#expand(OUT20Dir + "/var/ngm/allsites63ind/allsites63_mpileupM2020.vcf"),
		#expand(OUT20Dir + "/var/ngm/snponly63/snponly63_freebayes2020.vcf")
		#expand(OUT20Dir + "/ngsadmix/vdes31_maf001/{run}/sevenchr_maf001_{kcluster}.fopt.gz", run = RUNS, kcluster = KCLUSTERS)
		#expand(OUT20Dir + "/kraken_full/raw_fastq/{sample}_GENUS.bracken", sample = SAMPLES20)
		OUT20Dir + "/meta/hosts/hosts-20.txt"


#############################################################
##### HOST IDENTITY AND DIRECT DIET CHECK	       ######
#############################################################

rule checkHost:
        input:
                read1 = READDir + "/{sample}_R1_001.fastq.gz",
                read2 = READDir + "/{sample}_R2_001.fastq.gz",
        output:
                temp(OUT20Dir + "/meta/hosts/{sample}-{q}.txt")
        threads: 12
        shell:
                """
                bowtie2 -p {threads} -x {hostBeeMtBowtieIndex} -1 {input.read1} -2 {input.read2} | samtools view -S -q {wildcards.q} -F4 - | awk -v mellifera=0 -v cerana=0 -v sample={wildcards.sample} '$3~/^L/ {{mellifera++; next}}  {{cerana++}} END {{if(mellifera>cerana) print sample"\\tmellifera\\t"cerana"\\t"mellifera ; else print sample"\\tcerana\\t"cerana"\\t"mellifera}}' > {output}
                """

rule combineHost:
        input:
                expand(OUT20Dir + "/meta/hosts/{sample}-{{q}}.txt", sample = SAMPLES)
        output:
                OUT20Dir + "/meta/hosts/hosts-{q}.txt"
        shell:
                """
                (echo -ne "id\\thost\\tcerana\\tmellifera\\n"; cat {input}) > {output}
                """

rule removeHost:
        input:
                read1 = READDir + "/{sample}_R1_001.fastq.gz",
                read2 = READDir + "/{sample}_R2_001.fastq.gz",
        threads: 12
        output: temp(OUT20Dir + "/sketches/{sample}.fastq.gz")
        shell:
                """
                bowtie2 -p {threads} -x {hostBeeBowtieIndex} -1 {input.read1} -2 {input.read2} | samtools view -S -f12 | awk '{{print "@"$1"\\n"$10"\\n+\\n"$11}}' | gzip > {output}
                """


#############################################################
##### READS MAPPING (NGM/BOWTIE2) AND QUALITY CHECK    ######
#############################################################

rule nextgenmap:
	input:
		read1 = READDir + "/{sample}_R1_001.fastq.gz",
		read2 = READDir + "/{sample}_R2_001.fastq.gz",
	threads: 6
	output: 
		alignment = temp(OUT20Dir + "/alignments/ngm/{sample}.bam"), 
		index = temp(OUT20Dir + "/alignments/ngm/{sample}.bam.bai")
	shell:
                """
		ngm -t {threads} --qry1 {input.read1} --qry2 {input.read2} --paired -r {VDESRef} --local --very-sensitive --rg-id {wildcards.sample} --rg-sm {wildcards.sample} --rg-pl ILLUMINA --rg-lb NEXTERA --rg-cn OIST | samtools view -Su - | samtools sort - -m 20G -T {SCRATCH}/ngm/{wildcards.sample} -o - | samtools rmdup - - | variant - -m 200 --bam -o {output.alignment}
		samtools index {output.alignment}	
		"""

rule bowtie2:
	input:
		read1 = READDir + "/{sample}_R1_001.fastq.gz",
		read2 = READDir + "/{sample}_R2_001.fastq.gz",
	threads: 12
	output: 
		alignment = temp(OUT20Dir + "/alignments/bowtie2/{sample}.bam"), 
		index = temp(OUT20Dir + "/alignments/bowtie2/{sample}.bam.bai"),
		read1 = OUT20Dir + "/reads_unmapped/{sample}.1",
		read2 = OUT20Dir + "/reads_unmapped/{sample}.2"

	shell:
		"""
		bowtie2 -p {threads} --very-sensitive-local --sam-rg ID:{wildcards.sample} --sam-rg LB:NEXTERA --sam-rg SM:{wildcards.sample} --sam-rg PL:ILLUMINA --un-conc-gz {OUT20Dir}/reads_unmapped/{wildcards.sample} -x {VDESBowtieIndex} -1 {input.read1} -2 {input.read2} | samtools view -Su - | samtools sort - -m 20G -T {SCRATCH}/bowtie2/{wildcards.sample} -o - | samtools rmdup - - | variant - -m 200 --bam -o {output.alignment}
		samtools index {output.alignment}
		"""


rule kraken:
        input:
                read1 = READDir + "/{sample}_R1_001.fastq.gz",
                read2 = READDir + "/{sample}_R2_001.fastq.gz"
        threads: 4
        output:
                report = temp(OUT20Dir + "/kraken_full/raw_fastq/{sample}_report"),
                taxoreport = temp(OUT20Dir + "/kraken_full/raw_fastq/{sample}.kraken"),

        shell:
                """
                kraken2 --use-names --threads 4 --db {KRAKENVK} --report {output.report} --gzip-compressed --paired {input.read1} {input.read2} > {output.taxoreport}
		"""


rule braken:
        input:
                report = OUT20Dir + "/kraken_full/raw_fastq/{sample}_report"
        output:
                taxonomy = temp(OUT20Dir + "/kraken_full/raw_fastq/{sample}_GENUS.bracken")
        shell:
                """
                bracken -d {KRAKENVK} -i {input.report} -l G -o {output.taxonomy}
                """


rule statsbam:
        input:
                alignment = OUT20Dir + "/alignments/{aligner}/{sample}.bam"
        output:
                temp(OUT20Dir + "/meta/align_{aligner}/{sample}.txt")
        shell:
                """
		echo {wildcards.sample} > {output}
		samtools depth -a {input.alignment} | awk '{{sum+=$3}} END {{ print "Mean Average Coverage on all sites = ",sum/NR}}' >> {output}
		samtools depth -a -r NW_019211454.1 -r NW_019211455.1 -r NW_019211456.1 -r NW_019211457.1 -r NW_019211458.1 -r NW_019211459.1 -r NW_019211460.1 {input.alignment} | awk '{{sum+=$3}} END {{ print "Mean Average Coverage on 7 chromosomes = ",sum/NR}}' >> {output}
		samtools depth -a -r NC_004454.2 {input.alignment} | awk '{{sum+=$3}} END {{ print "Mean Average Coverage mtDNA = ",sum/NR}}' >> {output}
		samtools flagstat {input.alignment} >> {output}
		"""


#############################################################
##### VARIANT AND INVARIANT CALL (all-sites dataset)    #####
#############################################################

## Here we will use the all-site dataset for analysis regarding:
## Genetic diversity within a population (π)
## Genetic divergence between populations (dXY and dA)
## Demographic inferences using SFS (or AFS) with fastsimcoal2


rule mpileup_call:
	input:	expand(BAM20Dir + "/ngm/{sample}.bam", sample = SAMPLES20)
	output: temp(OUT20Dir + "/var/ngm/allsites63ind_split/allsites63_mpileupM2020.{region}.vcf")
	params:	bams = lambda wildcards, input: os.path.dirname(input[0]) + "/*.bam",
		span = lambda wildcards: REGIONS[wildcards.region],
		mapqual = "--min-MQ 10",
		idlist = OUT20Dir + "/list/varroa_63.txt"
	shell:
		"""
		bcftools mpileup {params.mapqual} --samples-file {params.idlist} --fasta-ref {VDESRef} {params.span} {params.bams} | bcftools call --multiallelic-caller -Ov -o {output}
		"""


rule merge_allsites:
	input: 
		expand(OUT20Dir + "/var/ngm/allsites63ind_split/allsites63_mpileupM2020.{region}.vcf", region = REGIONS)
	output:
		temp(OUT20Dir + "/var/ngm/allsites63ind/allsites63_mpileupM2020.vcf")
	shell:
		"""
		(grep "^#" {input[0]} ; cat {input} | grep -v "^#" ) | vcfuniq  > {output}
		"""


rule filter_allsites:
	input:
		rawvcf = OUT20Dir + "/var/ngm/allsites/allsites85_mpileupM2020.vcf",
                list = OUT20Dir + "/list/varroa_85ind.txt"
	output:
		OUT20Dir + "/var/ngm/allsites/allsites85_filtered2020"
	params:	cutoff = "--max-meanDP 38.2",
		filters = "--remove-indels --minDP 8 --minQ 30 --max-missing 1 --max-alleles 2",
		span = "--chr NW_019211454.1 --chr NW_019211455.1 --chr NW_019211456.1 --chr NW_019211457.1 --chr NW_019211458.1 --chr NW_019211459.1 --chr NW_019211460.1"
	shell:
                """
                vcftools --vcf {input.rawvcf} --keep {input.list} {params.span} {params.filters} {params.cutoff} --recode --recode-INFO-all --out {output}
                """


##### SUBSETTING THE DATA FOR FASTSIMCOAL2 DEMOGRAPHIC INFERENCES


rule finalfilter_50kb:
        input:
                rawvcf = OUT20Dir + "/var/ngm/allsites/allsites85_filtered2020.recode.vcf",
                list = OUT20Dir + "/list/gene_Vdes_3_50kb.txt"
        output:
                OUT20Dir + "/var/ngm/allsites/allsites43_50kbaway_2020"
        shell:
                """
                vcftools --vcf {input.rawvcf} --exclude-bed {input.list} --recode --recode-INFO-all --out {output}
                """


rule sortsubset_VDES:
        input:
                rawvcf = OUT19Dir + "/var/ngm/allsites/allsites43_50kbaway_2019.recode.vcf",
                list = OUT19Dir + "/list/VDES19_FSC_2019.txt"
        output:
                OUT19Dir + "/var/ngm/allsites/VDES19sort_50kb_2019_FSC26.vcf"
        shell:
                """
                bcftools view --samples-file {input.list} -Ov -o {output} {input.rawvcf}
                """

rule sortsubset_VJAC:
        input:
                rawvcf = OUT19Dir + "/var/ngm/allsites/allsites43_50kbaway_2019.recode.vcf",
                list = OUT19Dir + "/list/VJAC13_FSC_2019.txt"
        output:
                OUT19Dir + "/var/ngm/allsites/VJAC13sort_50kb_2019_FSC26.vcf"
        shell:
                """
                bcftools view --samples-file {input.list} -Ov -o {output} {input.rawvcf}
                """


#############################################################
##### VARIANT SITES CALL (SNP-only dataset)    		#####
#############################################################

## Here we will use the SNP-only dataset for analysis regarding:
## Exploratory PCA for population structure
## Admixture 
## Weir and Cockerham's FST
## Introgression test using ABBA-BABA

rule freebayes:
	input: 
		expand(OUT20Dir + "/alignments/ngm/{sample}.bam", sample = SAMPLES20)
	output: 
		temp(OUT20Dir + "/var/ngm/snponly63_split/freebayes.{region}.vcf")
	params: 
		span = lambda wildcards: REGIONS[wildcards.region],
		bams = lambda wildcards, input: os.path.dirname(input[0]) + "/*.bam",
		filtering = "--min-mapping-quality 10 --min-base-quality 5 --min-alternate-fraction 0.2 --use-best-n-alleles 4",
		idlist = OUT20Dir + "/list/varroa_63.txt"
	shell:
		"""
		freebayes {params.filtering} --fasta-reference {VDESRef} --samples {params.idlist} {params.span} {params.bams} | vcfallelicprimitives > {output}
		"""

rule merge_snponly:
        input:
                expand(OUT20Dir + "/var/ngm/snponly63_split/freebayes.{region}.vcf", region = REGIONS)
        output:
                temp(OUT20Dir + "/var/ngm/snponly63/snponly63_freebayes2020.vcf")
        shell:
                """
                (grep "^#" {input[0]} ; cat {input} | grep -v "^#" ) | vcfuniq  > {output}
                """


rule filter_snponly:
        input:
                rawvcf = OUT20Dir + "/var/ngm/snponly/snponly85_freebayes2020.vcf",
                list = OUT20Dir + "/list/varroa_85ind.txt"
        output:
                temp(OUT20Dir + "/var/ngm/snponly/snponly85_filtered2020")
        params: cutoff = "--max-meanDP 26.5", minor = "--maf 0.05",
                filters = "--remove-indels --minDP 8 --minQ 30 --max-missing 0.7 --max-alleles 2",
                span = "--chr NW_019211454.1 --chr NW_019211455.1 --chr NW_019211456.1 --chr NW_019211457.1 --chr NW_019211458.1 --chr NW_019211459.1 --chr NW_019211460.1"
        shell:
                """
                vcftools --vcf {input.rawvcf} --keep {input.list} {params.span} {params.minor} {params.filters} {params.cutoff} --recode --recode-INFO-all --out {output}
                """



#rule angsd:
#        input:
#                rawvcf = temp(OUT19Dir + "/var/ngm/raw_mpileupM2019.vcf"),
#                list = temp(OUT19Dir + "/list/varroa_43ind.txt")
#        output:
#angsd -b data_2020/alignments/ngm/*bam -ref ref2020/destructor/vdes_3_refseq.fasta -out ALL -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 3 -setMinDepth 7 -setMaxDepth 30 -doCounts 1 -GL 2 -doGlf 4



##### ADMIXTURE ANALYSIS

rule vcf2plink:
        input:	OUT20Dir + "/var/ngm/snponly/snponly85_filtered2020.recode.vcf"
        output: rawvcf = temp(OUT20Dir + "/LDprune/snponly85_ready4plink2020.vcf"), plink = temp(OUT20Dir + "/LDprune/snponly85_plink")
        shell:
                """
                cp {input} {output.rawvcf}
		sed -i 's/NW_019211454.1/1/g' {output.rawvcf}
               	sed -i 's/NW_019211455.1/2/g' {output.rawvcf} 
		sed -i 's/NW_019211456.1/3/g' {output.rawvcf}
		sed -i 's/NW_019211457.1/4/g' {output.rawvcf}
		sed -i 's/NW_019211458.1/5/g' {output.rawvcf}
		sed -i 's/NW_019211459.1/6/g' {output.rawvcf}
		sed -i 's/NW_019211460.1/7/g' {output.rawvcf}
		vcftools --vcf {output.rawvcf} --out {output.plink} --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --plink
		"""

rule plink_LDlist:
        input:  OUT20Dir + "/LDprune/snponly85_plink"
        output: temp(OUT20Dir + "/LDprune/snponly85_pruned")
        shell:
                """
                plink --file {input} --indep-pairwise 10 5 0.5 --out {output}
                """


rule LDpruning:
        input:  listin = OUT20Dir + "/LDprune/snponly85_pruned.prune.in", vcfin = OUT20Dir + "/var/ngm/snponly/snponly85_filtered2020.recode.vcf",
        output: listout = temp(OUT20Dir + "/LDprune/prune2keep.txt"), vcfout = temp(OUT20Dir + "/LDprune/snponly85_LDpruned")
        shell: 	"""
		cp {input.listin} {output.listout}
               	sed -i 's/1:/NW_019211454.1\t/g' {output.listout}
               	sed -i 's/2:/NW_019211455.1\t/g' {output.listout}
               	sed -i 's/3:/NW_019211456.1\t/g' {output.listout}
               	sed -i 's/4:/NW_019211457.1\t/g' {output.listout}
               	sed -i 's/5:/NW_019211458.1\t/g' {output.listout}
               	sed -i 's/6:/NW_019211459.1\t/g' {output.listout}
               	sed -i 's/7:/NW_019211460.1\t/g' {output.listout}
               	vcftools --vcf {input.vcfin} --positions {output.listout} --recode --recode-INFO-all --out {output.vcfout}
		"""


rule vcf2GL:
	input:  OUT20Dir + "/LDprune/snponly85_sorted_LDpruned.vcf"
        output: temp(OUT20Dir + "/ngsadmix/{chromosome}.BEAGLE.GL")
        shell:
                """
                vcftools --vcf {input} --chr {wildcards.chromosome} --out {OUT20Dir}/ngsadmix/{wildcards.chromosome} --BEAGLE-GL
                """


rule mergeGL:
        input: expand(OUT20Dir + "/ngsadmix/{chromosome}.BEAGLE.GL", chromosome = CHROMOSOMES)
        output: OUT20Dir + "/ngsadmix/snponly85_sevenchr.BEAGLE.GL"
        shell:
                """
                (head -1 {input[0]}; for i in {input}; do cat $i | sed 1d; done) > {output}
                """

rule NGSadmix:
        input: OUT20Dir + "/ngsadmix/vdes31_maf001/sevenchr_maf001.BEAGLE.GL"
        threads: 12
        output: temp(OUT20Dir + "/ngsadmix/vdes31_maf001/{run}/sevenchr_maf001_{kcluster}.fopt.gz")
        shell:
                """
                NGSadmix -P {threads} -likes {input} -K {wildcards.kcluster} -outfiles {OUT20Dir}/ngsadmix/vdes31_maf001/{wildcards.run}/sevenchr_maf001_{wildcards.kcluster}
                """


#############################################################
##### MITOGENOMES VARIANT CALL FOR SPECIES DIVERGENCE   #####
#############################################################

rule freebayes_mtDNA:
        input:
                expand(OUT20Dir + "/alignments/ngm/{sample}.bam", sample = SAMPLES20)
        output:
                rawvcf = temp(OUT20Dir + "/var/ngm/mtDNA/freebayes_mtDNA2020.vcf"),
                bcfout = temp(OUT20Dir + "/var/ngm/mtDNA/freebayes_mtDNA2020.vcf.gz")
        params:
                bams = lambda wildcards, input: os.path.dirname(input[0]) + "/*.bam"
        shell:
                """
                freebayes --min-mapping-quality 20 --min-base-quality 10 --min-alternate-count 3 --ploidy 1 --region NC_004454.2 --fasta-reference {VDESRef} {params.bams} | vcfallelicprimitives > {output.rawvcf}
                bgzip -c {output.rawvcf} > {output.bcfout}
                tabix -p vcf {output.bcfout}
                """

rule filter_mtDNA:
        input:  OUT20Dir + "/var/ngm/mtDNA/freebayes_mtDNA2020.vcf"
        output: temp(OUT20Dir + "/var/ngm/mtDNA/noindels_freebayes_mtDNA2020")
        params: filtering = "--remove-indels --minDP 8 --minQ 30 --max-missing 1"
        shell:
                """
                vcftools --vcf {input} {params.filtering} --recode --recode-INFO-all --out {output}
                """


### Here we use the consensus for haplotype identification
rule mtDNA_consensus:
        input:  OUT20Dir + "/var/ngm/mtDNA/noindels_freebayes_mtDNA2020.recode.vcf.gz"
        output: temp(OUT20Dir + "/fasta/mtDNA_consensus/{sample}_freebayesnoindels.fasta")
        shell:
                """
                samtools faidx {VDESRef} NC_004454.2 | bcftools consensus --iupac-codes --sample {wildcards.sample} {input} > {output}
                """


