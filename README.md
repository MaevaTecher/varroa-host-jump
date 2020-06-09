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



## Estimation of the mutation rate


## Demographic scenarios testing and parameters estimation with fastsimcoal2


## Contact
Questions about the data or scripts? please contact either:  

Maeva Techer, JSPS Postdoctoral Fellow, Ecology and Evolution lab @OIST  
OIST Email: maeva.techer@oist.jp  
Gmail: maeva.angelique.techer@oist.jp  

Alexander (Sasha) Mikheyev, Adjunct Professor Ecology and Evolution @OIST, Associate Professor Evolutionary Genomics@ANU
OIST Email: alexander.mikheyev@oist.jp  
ANU Email: alexander.mikheyev@anu.edu.au 



## TO do:
Add the vcf data to another repo or just leave the snakemake pipeline and explain using the DRA data uploaded on DDBJ
Add the tpl and est file for demographic inferences


