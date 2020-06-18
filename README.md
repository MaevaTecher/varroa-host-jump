# Population genomics and demographic history of successful host switches in Varroa mites

This repository aims at providing an online resource for the manuscript "Tracing the history of successful host switches in the beekeeping Varroa mite pests". 

Using whole-genome sequencing, (i) we compared sympatric populations genetic diversity, and (ii) estimated the demographic parameters of independent host switches in the _Apis_ honey bee ectoparasites mites: _Varroa destructor_ and _V. jacobsoni_. The following sections described how the bioinformatics analysis were processed and how to replicate them with available genome data reads, Snakemake pipeline and other file formats or codes necessary.


## Where can you find the data?

1. Varroa mite sequencing reads are available under the DDBJ _DNA Data Bank of Japan_ **BioProject PRJDB9195**  
2. Interactive visualization of sampling distribution and details is available [here](https://MaevaTecher.github.io/varroa-host-jump)  
3. Exploratory R codes used for the maps and other analysis such as PCAs can be found in `R_data`.
4. Snakemake pipeline is available in the file `Snakefile`, along with parameters file `cluster.json` and launcher `snakemake.slurm`.  
5. Raw, filtered and LD_pruned variant calling file `vcf` can be downloaded on **TO DO DRYAD REPO**  
6. Fastsimcoal2 input files for demographic scenarios and SFS subsets are available in `demography`.

## Customed Snakemake workflow for Varroa population genomics

### Software and dependencies necessary to run the Snakemake pipeline :  
(Sasha conda environment?)  

[`Samtools`](http://www.htslib.org/) : the version used here was samtools/1.3.1.
[`Mashtree`](https://github.com/lskatz/mashtree)

For reads mapping and filtering:  
[`Bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml): version bowtie2/2.2 used here.   
[`NextGenMap`](https://cibiv.github.io/NextGenMap/): version NextGenMap/0.5.0 used here.   
[`VariantBam`](https://github.com/broadinstitute/VariantBam): version VariantBam/1.4.3 used here.   

For variant calling:
[`Freebayes`](https://github.com/ekg/freebayes): version freebayes/1.1.0 used here.  

For .vcf files manipulation and population genomics analysis:
[`VCFtools`](https://vcftools.github.io/index.html): version vcftools/0.1.12b used here.  
[`vcflib`](https://github.com/vcflib/vcflib): version vcflib/1.0.0-rc1 used here.  
[`PLINK`](https://www.cog-genomics.org/plink/): version plink/1.90b3.36 used here.  
[`R`](https://www.r-project.org/): version R/3.4.2 used here.  
[`NGSadmix`](http://www.popgen.dk/software/index.php/NgsAdmix): version NGSadmix/32 used here.  

### How to run the pipeline?  

The first 72 lines of the `Snakemake` file are setting in which environment we ran our analysis (i.e., which folders the raw reads data can be found, and the paths of future generated resulting analysis data). We computed most analysis on Sango, OIST Super Computer Center.

You can change the targets of `rule all` depending on which particular steps of the pipeline you want to proceed following [Snakemake's documentation](https://snakemake.readthedocs.io/en/v3.9.1/).   

For example, if you want to generate the mapped read .bam files from the original .fastq files downloaded on NCBI:  
1. Place the fastq files in a `reads` folder  
2. Samples name will be recognized according to lines 45 for which fastq R1 files end with `{sample}-R1_001.fastq.gz`  
3. Apply rule all as `expand(outDir + "/alignments/ngm/{sample}.bam", sample = SAMPLES)`  

Or, if you want to generate the variant call file from NextGenMap reads apply:
`outDir + "/var/ngm/filtered.vcf"` for which `outDir` is the path of your working folder.


## Estimation of the mutation rate

TBD


## Demographic scenarios testing and parameters estimation with fastsimcoal2

We designed six evolutionary scenario using [`fastsimcoal2`](http://cmpg.unibe.ch/software/fastsimcoal2/) to choose the most likely and estimate demographic parameters such as the size of founding population at time of independent host switch events by _V. destructor_ and _V. jacobsoni_ from _A. cerana_ to _A. mellifera_.

### Useful links and tools for running your own demographic inferences

Additionally to the complete [manual available for `fastsimcoal2`](http://cmpg.unibe.ch/software/fastsimcoal2/man/fastsimcoal26.pdf) written by Laurent Excoffier, we found help on the active [google group forum](https://groups.google.com/forum/?nomobile=true#!forum/fastsimcoal). We are extremely grateful to the excellent tutorial that can be found [here](https://speciationgenomics.github.io/fastsimcoal2/), hosted on the[speciationgenomics](https://github.com/speciationgenomics) Github page by Mark Ravinet & Joana Meier. 

The SFS input for fsc26 were generated from our vcf files using the [easySFS](https://github.com/isaacovercast/easySFS) conversion scripts developped by Isaac Overcast.
For 

2D-SFS were plotted from the observed and expected SFS under each scenario model using [SFS-scripts](https://github.com/marqueda/SFS-scripts) developped by David Marques.

### Templates files in `demography` folder:  

Here I described the command line used for one demographic model (`VDNOV475_105000/scenario6_IMGwt`) with the SFS generated from the largest SNP subset for _V. destructor_ (`dataset_VD4_38108snps.obs`). The `M6_downsample_IMGwt_1_jointMAFpop1_0.obs` file should be place in the same folder as the same named .tpl and .est files.

Briefly the `M6_downsample_IMGwt_1.tpl` file draw the scenario with :  
1. A two populations model between _V. destructor_ mites on original host _A. cerana_ with a current effective population size `NVDAC1` and `NVDAM0` for the novel host _A. mellifera_.
2. We projeted the data using easySFS with 17 haploid genomes for both _A. mellifera_ (population 0) and _A. cerana_ (population 1) mites.
3. We assumed a growth rate GAM since the host switch event only for the novel host (expansion biologically known).
4. Two migration matrices are given for the population split and after.
5. We considered a single historical event `TBOT` from which a number of haploid mites `JUMP` splited from _A. cerana_ population to found the new _A. mellifera_ population.
6. The mutation rate was proposed following preliminary _de novo_ mutations estimations.  
  
In the `M6_downsample_IMGwt_1.est` file, most parameters were sampled in a uniform distribution (NB: upper limit does not constitute a maximum bound for fsc26, see manual).  
  
We copied-named `M6_downsample_IMGwt_XXX.est`, `M6_downsample_IMGwt_XXX.tpl`, `M6_downsample_IMGwt_XXX_jointMAFpop1_0.obs` 100 times for which XXX is in {1..100}. Using an array bash script we then ran the following command for each replicate.  
  
==================   
`#!/bin/bash`  
`#SBATCH --job-name=vdM6`  
`#SBATCH --partition=XXX`  
`#SBATCH --mem=5G`  
`#SBATCH -c 10`  
`#SBATCH --time=3:00:00`  
  
`number=$SLURM_ARRAY_TASK_ID`  
`[PATH_TO_FASTSIMCOAL2]/fsc26 	--tplfile M6_downsample_IMGwt_"$number".tpl \  
								--estfile M6_downsample_IMGwt_"$number".est \  
								-m \  
								--numsims 1000000 \  
								--maxlhood 0.001 \  
								--minnumloops 20 \  
								--numloops 100 \  
								-c 10`  
   
================  

Then to extract the best results from each run, we simply apply the following command lines:  
`cat M6_downsample_IMGwt_1/M6_downsample_IMGwt_1.bestlhoods >> scenario6_vd.txt`  
`for i in {2..100}; do sed -n 2p M6_downsample_IMGwt_"$i"/M6_downsample_IMGwt_"$i".bestlhoods >> scenario6_vd.txt; done`  
`for i in {1..100}; do cat M6_downsample_IMGwt_"$i"/M6_downsample_IMGwt_"$i".bestlhoods >> scenario6_vd.txt; done`  
`cat scenario6_vd.txt | wc -l` # to check that we have all replicates results  
`cat scenario6_vd.txt | sort -k 10nr` # the first line is then the scenatio with the lowest MaxEstLhood  

The bootstraps from the best replicate in the best scenario were performed following tutorial in the manual (page 56-57). 

## Contact
Questions about the data or scripts? please contact either:  

Maeva Techer, JSPS Postdoctoral Fellow, Ecology and Evolution lab @OIST  
OIST Email: maeva.techer@oist.jp  
Gmail: maeva.angelique.techer@oist.jp  

Alexander (Sasha) Mikheyev, Adjunct Professor Ecology and Evolution @OIST, Associate Professor Evolutionary Genomics@ANU
OIST Email: alexander.mikheyev@oist.jp  
ANU Email: alexander.mikheyev@anu.edu.au 

## TO do:
Dryad repo to add once the title is finalized




