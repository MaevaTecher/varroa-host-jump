#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}

infile <- args[1]
outfile <- args[2]

library(OutFLANK)
library(data.table)
variants <- fread(infile)
FstDataFrame <- MakeDiploidFSTMat(variants[,-c(1:3), with = FALSE], names(variants)[-c(1:3)], variants$host)
OF <- OutFLANK(FstDataFrame, NumberOfSamples = 2, qthreshold = 0.05, RightTrimFraction = 0.05)
saveRDS(list(fst = FstDataFrame, of = OF), outfile)
