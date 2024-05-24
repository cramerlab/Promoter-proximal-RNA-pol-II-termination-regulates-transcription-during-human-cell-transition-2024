

### rna.to.dna function ###


rna.to.dna = function(h5){
  fastq = strsplit(h5,split = "\n")
  fastq[[1]][2] = gsub("U","T",fastq[[1]][2])
  paste(c(unlist(fastq),""),collapse = "\n")
}



