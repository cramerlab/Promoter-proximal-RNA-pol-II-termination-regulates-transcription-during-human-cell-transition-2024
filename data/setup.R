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
  source(file.path("data","Utility",i))}

source(file.path("data", "processing.functions.R"))

functions = list.files(file.path("data","GatheringRoutines"))
for (i in functions){
  source(file.path("data","GatheringRoutines",i))}
