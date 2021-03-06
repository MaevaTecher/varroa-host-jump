# Online resource "The first steps toward a global pandemic: Reconstructing the demographic history of parasite host switches in its native range"

This repository aims at providing an online resource for the [preprint](https://www.biorxiv.org/content/10.1101/2020.07.30.228320v1.full)

Using whole-genome sequencing, (i) we compared sympatric populations genetic diversity, and (ii) estimated the demographic parameters of independent host switches in the _Apis_ honey bee ectoparasites mites: _Varroa destructor_ and _V. jacobsoni_. The following sections described how the bioinformatics analysis were processed and how to replicate them with available genome data reads, Snakemake pipeline and other file formats or codes necessary.

<img src="/images/Varroabanner.jpg" alt="Honey bee Varroa mites"/>

## Where can you find the data?

### Varroa mites whole-genome sequences

1. Varroa mite sequencing reads are available under the DDBJ _DNA Data Bank of Japan_ [BioProject PRJDB9195](https://www.ncbi.nlm.nih.gov/bioproject/PRJDB9195)   
2. [GCF_002443255.1 Vdes_3.0](https://www.ncbi.nlm.nih.gov/genome/?term=txid109461[orgn]) Reference genome can be directly downloaded from NCBI developed by [Techer et al. 2020. Comms. Biology](https://www.nature.com/articles/s42003-019-0606-0)  
3. Honey bees and Varroa mites fasta sequences and Bowtie2 index build are available in the `ref2020` folder on **DRYAD DOI TO ADD**  

### Codes

1. Genomics analysis are summarized into a Snakemake pipeline is available in the file `Snakefile`, along with parameters file `cluster.json` and launcher `snakemake.slurm`.
2. All scripts called in `Snakefile` are present in the `scripts` folder.  
3. Interactive visualization of sampling distribution and details is available [here](https://MaevaTecher.github.io/varroa-host-jump).   
4. R markdown including the code lines used to generate the interactive maps, PCAs analysis and bootstraping for the demographic inferences can be found in `R_data`.

### Input or output files 

1. Variant calling files (Raw, filtered and LD_pruned) as well as some lists (samples ID, SNP list, ) **DRYAD DOI to add**  
2. Fastsimcoal2 input files for demographic scenarios and SFS subsets are available in `demography`.
3. FASTA files generated for mitochondrial analysis (COX1 and COX1-COX3-ATP6-CYTB) and compared to _Varroa_ reference sequences are available in `alignment_mtDNA`

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

All files are available in the DRYAD repository.


## Demographic scenarios testing and parameters estimation with fastsimcoal2

We designed six evolutionary scenario using [`fastsimcoal2`](http://cmpg.unibe.ch/software/fastsimcoal2/) to choose the most likely and estimate demographic parameters such as the size of founding population at time of independent host switch events by _V. destructor_ and _V. jacobsoni_ from _A. cerana_ to _A. mellifera_.

### Useful links and tools for running your own demographic inferences

Additionally to the complete [manual available](http://cmpg.unibe.ch/software/fastsimcoal2/man/fastsimcoal26.pdf) for `fastsimcoal2` written by Laurent Excoffier, we recommend the helpful active [google group forum](https://groups.google.com/forum/?nomobile=true#!forum/fastsimcoal). 

An excellent tutorial on fastsimcoal2 usage on genomics data can be found [here](https://speciationgenomics.github.io/fastsimcoal2/), hosted on the [speciationgenomics](https://github.com/speciationgenomics) Github page by Mark Ravinet & Joana Meier. 

[easySFS](https://github.com/isaacovercast/easySFS) was used to generate SFS input files from our vcf files (developed by Isaac Overcast).

[SFS-scripts](https://github.com/marqueda/SFS-scripts) were used to plot observed and expected 2D-SFS under each scenario model (developed by David Marques).

### Templates files in `demography` folder:  

Here, we described the command lines used for one demographic model (`templates_vjac_MUTINB/14_mig_botgwt/mut_8e-10`) with the SFS generated from all sites subset for _V. destructor_ (`VDES32by22_folded_50kbvcf2sfs.txt`). The `mig_botgwt_TEMPLATE_jointMAFpop1_0.obs` file should be place in the same folder as the same named .tpl and .est files.

Briefly the `mig_botgwt_TEMPLATE.tpl` file draw the scenario with :  
1. A two populations model between _V. destructor_ mites on original host _A. cerana_ with a current effective population size `NVAC1` and `NVAM0` for the novel host _A. mellifera_.
2. Observed SFS Data were projected using easySFS with 17 haploid genomes for both _A. mellifera_ (population 0) and _A. cerana_ (population 1) mites.
3. A growth rate GAM since the host switch event only for the novel host (expansion biologically known).
4. Two migration matrices were given for the population split and after.
5. We considered a single historical event `TJUMP` ending at `TBOTEND` from which a number of haploid mites `NBOTAM` splited from _A. cerana_ population to found the new _A. mellifera_ population.
6. The mutation rate was proposed following preliminary _de novo_ mutations estimations.  
  
In the `mig_botgwt_TEMPLATE.est` file, most parameters were sampled in a uniform distribution (NB: upper limit does not constitute a maximum bound for fsc26, see manual).  
  
We copied-named `mig_botgwt_TEMPLATE_XXX.est`, `mig_botgwt_TEMPLATE_XXX.tpl`, `mig_botgwt_TEMPLATE_XXX_jointMAFpop1_0.obs` 100 times for which XXX is in {1..100}. Using an array bash script we then ran the following command for each replicate.  

__________________________
`#!/bin/bash`  
`#SBATCH --job-name=vd-migbotgwt`  
`#SBATCH --partition=XXX`  
`#SBATCH --mem=5G`  
`#SBATCH -c 10`  
`#SBATCH --time=3:00:00`  
  
`number=$SLURM_ARRAY_TASK_ID`  
`[PATH_TO_FASTSIMCOAL2]/fsc26 --tplfile mig_botgwt_TEMPLATE_"$number".tpl --estfile mig_botgwt_TEMPLATE_"$number".est -m --numsims 1000000 --maxlhood 0.001 --minnumloops 20 --numloops 100 -c 10`  
__________________________


Then to extract the best results from each run, we simply apply the following command lines:  
__________________________
`cat mig_botgwt_TEMPLATE_1/mig_botgwt_TEMPLATE_1.bestlhoods >> scenario_migbotgwt.txt`  
`for i in {2..100}; do sed -n 2p mig_botgwt_TEMPLATE_"$i"/mig_botgwt_TEMPLATE_"$i".bestlhoods >> scenario_migbotgwt.txt; done`  
`for i in {1..100}; do cat mig_botgwt_TEMPLATE_"$i"/mig_botgwt_TEMPLATE_"$i".bestlhoods >> scenario_migbotgwt.txt; done`  
`cat scenario_migbotgwt.txt | wc -l` # to check that we have the results for the 100 replicates 
`cat scenario_migbotgwt.txt | sort -k 10nr` # the first line is then the scenatio with the lowest MaxEstLhood  
__________________________

For the best scenario, bootstraps from the best replicate were performed following the tutorial in the manual (page 56-57). We simulated 100 SFS datasets from the best output `.par` file after modifying it to generate a DNA sequence data, using:   
`[PATH_TO_FASTSIMCOAL2]/NUMBER=YYY` # where NUMBER is the replicate number with the lowest MaxEstLhood
`fsc26 -i mig_botgwt_TEMPLATE_${NUMBER}_boot.par -n100 -j -m -s0 -x -I -q`

We then repeat the parameters estimation 100 times for each of the 100 simulated SFS `mig_botgwt_TEMPLATE_${NUMBER}_boot_jointMAFpop1_0.obs` using the same initial command lines. 

The best replicate for each simulated SFS run is then used to obtain the confidence interval of `NVAM0`, `NVAC1`, `TJUMP`, `NBOTAM`, `GAM` and migration rates.


## Contact
Questions about the data or scripts? please contact either:  

Maeva Techer, JSPS Postdoctoral Fellow, Ecology and Evolution lab @OIST  
OIST Email: maeva.techer@oist.jp  
Gmail: maeva.angelique.techer@oist.jp  

Alexander (Sasha) Mikheyev, Adjunct Professor Ecology and Evolution @OIST, Associate Professor Evolutionary Genomics@ANU
OIST Email: alexander.mikheyev@oist.jp  
ANU Email: alexander.mikheyev@anu.edu.au 


