### libraries ###


library(Biostrings)
library(Rsamtools)
library(foreach)
library(doParallel)
library(GenomicRanges)
library(GenomicAlignments)
library(rtracklayer)
library(ggplot2)
library(DESeq2)

### processing functions

functions = list.files(file.path(prewd,"data","Utility"))
for (i in functions){
  source(file.path(prewd,"data","Utility",i))}

source(file.path(prewd, "data", "processing.functions.R"))

functions = list.files(file.path(prewd,"data","GatheringRoutines"))
for (i in functions){
  source(file.path(prewd,"data","GatheringRoutines",i))}
