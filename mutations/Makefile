# This makefile powers an analysis pipeline for
# finding possible mutations in the Varroa samples.
# TL;DR there are no good mutations.

# paths to input files
INPUT=all.vcf.gz
PED=varroa.ped

# paths to binaries used in this pipeline
DNGCALL=dng call
BCFTOOLS=bcftools
BEDTOOLS=bedtools
RSCRIPT=Rscript --vanilla

DNGARGS= --theta=0.0001 \
	  --lib-error=0.0005 \
	  --ref-bias-hom=0 \
	  --ref-bias-het=0 \
	  --mu=1e-7 \
	  --lib-overdisp-hom=0.0 \
	  --lib-overdisp-het=0.0 \
	  --lib-bias=1

# housekeeping stuff

default: all

all: genome.txt \
	0_variants.vcf.gz 1_denovos.vcf.gz \
	2_denovos.libs.vcf.gz 3_deduped.vcf.gz \
	4_filtered.vcf.gz

all: random_sites.bed

all: stats_rle_dnl.csv stats.csv

.PHONY: default all

.DELETE_ON_ERROR:

# a bedtools genome file for the input contigs
genome.txt: $(INPUT)
	$(BCFTOOLS) view -h $< | awk -F '[=,<>]' '$$1 == "##contig" {print $$4 "\t" $$6}' > $@

# use only contigs that are at least 1mb in size
contigs.txt: genome.txt
	awk '$$2 > 1000000 {print $$1 "\t" $$2}' $< > $@

# create a bed file for the contigs
contigs.bed: contigs.txt
	awk '{print $$1 "\t0\t" $$2}' $< > $@

# Sample 6k random sites from our contigs, can be used to estimate denominator
random_sites.bed: contigs.txt
	$(BEDTOOLS) random -seed 6721 -l 1 -n 6000 -g $< | $(BEDTOOLS) sort -g $< > $@

# Return only the variants that are with our target contigs.
# And remove those that are within 100bp of an indel.
0_variants.vcf.gz: $(INPUT) contigs.bed
	$(BCFTOOLS) view -T contigs.bed $< | $(BCFTOOLS) filter -g 100 -O z -o $@ 

# Call denovo mutations using `dng call`
1_denovos.vcf.gz: 0_variants.vcf.gz
	$(DNGCALL) -p varroa.ped -o $@ $(DNGARGS) $<

# Identify the samples that correspond to libraries.
libs.txt: 1_denovos.vcf.gz
	$(BCFTOOLS) query -l $< | grep '^LB/' > $@

# Process calls extracting only the columns for libraries
# Add AC etc. tags, and return only biallelic sites
2_denovos.libs.vcf.gz: 1_denovos.vcf.gz libs.txt
	$(BCFTOOLS) view -S libs.txt $< | $(BCFTOOLS) +fill-tags | $(BCFTOOLS) view -m 2 -M 2 -O z -o $@

# The low-coverage samples tend to produce denovo calls in "clusters"
# Anything with more then 7 is unlikely to occur randomly and is thus suspect
# (assumes equal mutation rate and sample power)
stats_rle_dnl.csv: 2_denovos.libs.vcf.gz
	$(RSCRIPT) scripts/rle_dnl.R

# A bed file containing denovo calls
2_denovos.bed: 2_denovos.libs.vcf.gz
	$(BCFTOOLS) query -f '%CHROM\t%POS0\t%POS\n' $< > $@

# Identify denovo calls that are not a part of runs of 7 or more calls in the same sample
2_deruns.bed: 2_denovos.libs.vcf.gz
	$(RSCRIPT) scripts/rle_dnl_filter.R

# Remove denovos clustered within 100 nt of one another
# and not in long runs of calls in the same sample
2_deduped.bed: 2_denovos.libs.vcf.gz 2_deruns.bed
	$(BEDTOOLS) merge -c 1 -o count -d 100 -i $< | awk '$$4==1' \
		| $(BEDTOOLS) intersect -g contigs.txt -a - -b 2_deruns.bed > $@

# Generate a deduped file
3_deduped.vcf.gz: 2_denovos.libs.vcf.gz 2_deduped.bed
	$(BCFTOOLS) view -o $@ -O z -T 2_deduped.bed $<

# Select only samples with (1) one mutant allele, (2) P(denovo) >= 0.75, and (3) fit to model >= -3
4_filtered.vcf.gz: 3_deduped.vcf.gz
	$(BCFTOOLS) view -o $@ -O z -i 'AC==1 && INFO/DNP >= 0.75 && LLS >= -3' $<

# The result after all of this is that 4_filtered.vcf.gz contains only two sites.
# Both sites have a mutation in Son7, which has low coverage.

# Generate a statistics file for our analysis
#   Indiv: Sample
#   AvgDP: Average value of DP in the 0_variants file.
#   GTest: p-value for a test that mutations are 2_denovo calls are distributed evenly to chromosomes
#	MeanSpacing: average spacing between calls
#   MedianSpacing: median spacing between calls (the spacings of some samples show questionable skew)
#	AdTest: p-value for an AD test that the spacing follow an exponential distribution
# 	Count2: number of calls for the sample in 2_denovos
#	Count3: number of calls in 3_deduped
# 	Count4: number of calls in 4_filtered
# This will show that mutation calls come from low-coverage samples.
stats.csv: 0_variants.vcf.gz 2_denovos.libs.vcf.gz 3_deduped.vcf.gz 4_filtered.vcf.gz
	$(RSCRIPT) scripts/stats.R


############################################################################################
# Denominator / Recall rate

recall: recall_0_genotype.vcf.gz recall_1_mutated.vcf.gz \
        recall_2_denovos.vcf.gz  \
        recall_7_filtered.vcf.gz

.PHONY: recall

recall_0_genotype.vcf.gz: all.random.vcf.gz libs.txt
	$(DNGCALL) --all -m 0 -p varroa.ped -o - $(DNGARGS) $< \
	| $(BCFTOOLS) view -S libs.txt -i 'REF != "N" && REF != "."' \
	| $(BCFTOOLS) annotate -O z -x '^INFO/DNL,^FMT/GT,FMT/AD' > $@

recall_1_mutated.vcf.gz: recall_0_genotype.vcf.gz random_sites.bed
	$(RSCRIPT) scripts/make_mutations.R
	gunzip $@
	bgzip $(patsubst %.gz,%,$@)

recall_2_denovos.vcf.gz: recall_1_mutated.vcf.gz
	$(DNGCALL) -m 0 -p varroa.ped -o $@ $(DNGARGS) $<

recall_3_merged.vcf.gz: 1_denovos.vcf.gz recall_2_denovos.vcf.gz 1_denovos.vcf.gz.csi recall_2_denovos.vcf.gz.csi
	$(BCFTOOLS) concat -a 1_denovos.vcf.gz recall_2_denovos.vcf.gz | \
	$(BCFTOOLS) sort -O z -o $@

recall_4_denovos.vcf.gz: recall_3_merged.vcf.gz
	$(BCFTOOLS) view -S libs.txt $< | $(BCFTOOLS) +fill-tags | $(BCFTOOLS) view -m 2 -M 2 -O z -o $@

recall_4_deruns.bed: recall_4_denovos.vcf.gz
	$(RSCRIPT) scripts/rle_dnl_filter_recall.R

# Remove denovos clustered within 100 nt of one another
# and not in long runs of calls in the same sample
recall_4_deduped.bed: recall_4_denovos.vcf.gz recall_4_deruns.bed
	$(BEDTOOLS) merge -c 1 -o count -d 100 -i $< | awk '$$4==1' \
		| $(BEDTOOLS) intersect -g contigs.txt -a - -b recall_4_deruns.bed > $@

# Generate a deduped file
recall_5_deduped.vcf.gz: recall_4_denovos.vcf.gz recall_4_deduped.bed
	$(BCFTOOLS) view -o $@ -O z -T recall_4_deduped.bed $<

recall_6_filtered.vcf.gz: recall_5_deduped.vcf.gz
	$(BCFTOOLS) view -o $@ -O z -i 'AC==1 && INFO/DNP >= 0.75 && LLS >= -3' $<

recall_7_filtered.vcf.gz: recall_6_filtered.vcf.gz
	$(BEDTOOLS) subtract -header -a $< -b 4_filtered.vcf.gz | bgzip > $@

%.vcf.gz.csi: %.vcf.gz
	$(BCFTOOLS) index $<

#######################3
# 3464/6000 simulated mutations make it through the pipeline
# recall is 57.7%
# haploid contig length is 362,833,908
# mutation rate is estimated to be less than 1 in 1,256,856,657
