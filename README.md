# Population genomics and demographic history of successful host switches in Varroa mites

This repository aims at providing an online resource for the manuscript "Tracing the history of successful host switches in the beekeeping Varroa mite pests". 

Using whole-genome sequencing, (i) we compared sympatric populations genetic diversity, and (ii) estimated the demographic parameters of independent host switches in the _Apis_ honey bee ectoparasites mites: _Varroa destructor_ and _V. jacobsoni_. The following sections described how the bioinformatics analysis were processed and how to replicate them with available genome data reads, Snakemake pipeline and other file formats or codes necessary.

<img src="/images/Varroabanner.jpg" alt="Honey bee Varroa mites"/>

## Where can you find the data?

1. Varroa mite sequencing reads are available under the DDBJ _DNA Data Bank of Japan_ **BioProject PRJDB9195**  
2. Interactive visualization of sampling distribution and details is available [here](https://MaevaTecher.github.io/varroa-host-jump)  
3. Exploratory R codes used for the maps and other analysis such as PCAs can be found in `R_data`.
4. Snakemake pipeline is available in the file `Snakefile`, along with parameters file `cluster.json` and launcher `snakemake.slurm`.  
5. Raw, filtered and LD_pruned variant calling file `vcf` can be downloaded on **DRYAD DOI to add**  
6. Fastsimcoal2 input files for demographic scenarios and SFS subsets are available in `demography`.

## Customed Snakemake workflow for Varroa population genomics

### Software and dependencies necessary to run the Snakemake pipeline :  

[`Samtools`](http://www.htslib.org/) : version samtools/1.3.1 used here.  
[`Mashtree`](https://github.com/lskatz/mashtree)   
[`Bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml): version bowtie2/2.2 used here.   
[`NextGenMap`](https://cibiv.github.io/NextGenMap/): version NextGenMap/0.5.0 used here.   
[`VariantBam`](https://github.com/broadinstitute/VariantBam): version VariantBam/1.4.3 used here.   
[`Freebayes`](https://github.com/ekg/freebayes): version freebayes/1.1.0 used here.  
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
2. Samples name will be recognized according to line 45 for which fastq R1 files end with `{sample}-R1_001.fastq.gz`  
3. Apply rule all as `expand(outDir + "/alignments/ngm/{sample}.bam", sample = SAMPLES)`  

Or, if you want to generate the variant call file from NextGenMap reads apply:
`outDir + "/var/ngm/filtered.vcf"` for which `outDir` is the path of your working folder.


## Estimation of the mutation rate

TBD


## Demographic scenarios testing and parameters estimation with fastsimcoal2

We designed six evolutionary scenario using [`fastsimcoal2`](http://cmpg.unibe.ch/software/fastsimcoal2/) to choose the most likely and estimate demographic parameters such as the size of founding population at time of independent host switch events by _V. destructor_ and _V. jacobsoni_ from _A. cerana_ to _A. mellifera_.

### Useful links and tools for running your own demographic inferences

Additionally to the complete [manual available](http://cmpg.unibe.ch/software/fastsimcoal2/man/fastsimcoal26.pdf) for `fastsimcoal2` written by Laurent Excoffier, we recommend the helpful active [google group forum](https://groups.google.com/forum/?nomobile=true#!forum/fastsimcoal). 

An excellent tutorial on fastsimcoal2 usage on genomics data can be found [here](https://speciationgenomics.github.io/fastsimcoal2/), hosted on the [speciationgenomics](https://github.com/speciationgenomics) Github page by Mark Ravinet & Joana Meier. 

[easySFS](https://github.com/isaacovercast/easySFS) was used to generate SFS input files from our vcf files (developed by Isaac Overcast).

[SFS-scripts](https://github.com/marqueda/SFS-scripts) were used to plot observed and expected 2D-SFS under each scenario model (developed by David Marques).

### Templates files in `demography` folder:  

Here, we described the command lines used for one demographic model (`VDNOV475_105000/scenario6_IMGwt`) with the SFS generated from the largest SNP subset for _V. destructor_ (`dataset_VD4_38108snps.obs`). The `M6_downsample_IMGwt_1_jointMAFpop1_0.obs` file should be place in the same folder as the same named .tpl and .est files.

Briefly the `M6_downsample_IMGwt_1.tpl` file draw the scenario with :  
1. A two populations model between _V. destructor_ mites on original host _A. cerana_ with a current effective population size `NVDAC1` and `NVDAM0` for the novel host _A. mellifera_.
2. Observed SFS Data were projected using easySFS with 17 haploid genomes for both _A. mellifera_ (population 0) and _A. cerana_ (population 1) mites.
3. A growth rate GAM since the host switch event only for the novel host (expansion biologically known).
4. Two migration matrices were given for the population split and after.
5. We considered a single historical event `TBOT` from which a number of haploid mites `JUMP` splited from _A. cerana_ population to found the new _A. mellifera_ population.
6. The mutation rate was proposed following preliminary _de novo_ mutations estimations.  
  
In the `M6_downsample_IMGwt_1.est` file, most parameters were sampled in a uniform distribution (NB: upper limit does not constitute a maximum bound for fsc26, see manual).  
  
We copied-named `M6_downsample_IMGwt_XXX.est`, `M6_downsample_IMGwt_XXX.tpl`, `M6_downsample_IMGwt_XXX_jointMAFpop1_0.obs` 100 times for which XXX is in {1..100}. Using an array bash script we then ran the following command for each replicate.  

__________________________
`#!/bin/bash`  
`#SBATCH --job-name=vdM6`  
`#SBATCH --partition=XXX`  
`#SBATCH --mem=5G`  
`#SBATCH -c 10`  
`#SBATCH --time=3:00:00`  
  
`number=$SLURM_ARRAY_TASK_ID`  
`[PATH_TO_FASTSIMCOAL2]/fsc26 --tplfile M6_downsample_IMGwt_"$number".tpl --estfile M6_downsample_IMGwt_"$number".est -m --numsims 1000000 --maxlhood 0.001 --minnumloops 20 --numloops 100 -c 10`  
__________________________


Then to extract the best results from each run, we simply apply the following command lines:  
__________________________
`cat M6_downsample_IMGwt_1/M6_downsample_IMGwt_1.bestlhoods >> scenario6_vd.txt`  
`for i in {2..100}; do sed -n 2p M6_downsample_IMGwt_"$i"/M6_downsample_IMGwt_"$i".bestlhoods >> scenario6_vd.txt; done`  
`for i in {1..100}; do cat M6_downsample_IMGwt_"$i"/M6_downsample_IMGwt_"$i".bestlhoods >> scenario6_vd.txt; done`  
`cat scenario6_vd.txt | wc -l` # to check that we have the results for the 100 replicates 
`cat scenario6_vd.txt | sort -k 10nr` # the first line is then the scenatio with the lowest MaxEstLhood  
__________________________

For the best scenario, bootstraps from the best replicate were performed following the tutorial in the manual (page 56-57). We simulated 100 SFS datasets from the best output `.par` file after modifying it to generate a DNA sequence data, using:   
`[PATH_TO_FASTSIMCOAL2]/NUMBER=YYY` # where NUMBER is the replicate number with the lowest MaxEstLhood
`fsc26 -i M6_downsample_IMGwt_${NUMBER}_boot.par -n100 -j -m -s0 -x -I -q`

We then repeat the parameters estimation 100 times for each of the 100 simulated SFS `M6_downsample_IMGwt_${NUMBER}_boot_jointMAFpop1_0.obs` using the same initial command lines. 

The best replicate for each simulated SFS run is then used to obtain the confidence interval of `NVDAM0`, `NVAC1`, `TBOT`, `JUMP`, `GAM` and migration rates.


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




