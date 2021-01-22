//Parameters for the coalescence simulation program : fsimcoal2.exe
2 samples to simulate :
//Population effective sizes (number of genes)
NVAM0
NVAC1
//Haploid Samples sizes and samples age 
18 0 0.8 
6 0 0.8
//Growth rates	: negative growth implies population expansion
GAM
0
//Number of migration matrices : 0 implies no migration between demes
2
//Migration matrix 0
0 MAMtoAC
MACtoAM 0
//Migration matrix 1: No migration
0 0
0 0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
2  historical event
TBOTEND 0 0 0 ResBOT 0 1
TJUMP 0 1 1 ResJUMP 0 1
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ  1   0   4.0e-10 OUTEXP
