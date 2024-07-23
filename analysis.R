########################################################################################
# Promoter-proximal RNA polymerase II termination regulates transcription 
# during human cell type transition


# L Kseniia, D Arjun et al., Nat. Struct. Mol. Biol, 2024
# DOI : 

# Script for analysis as performed in the article
# Authors : Arjun Devadas, MPI-NAT
########################################################################################


#####################################################
# running prerequisites
#####################################################

# loads the required libraries and custom functions used to run this script
# in case of missing libraries, please install the appropriate version for the R used
# version details of packages used in this script is available in info.txt


# setting the directory for running Setup.R as path to the location of the directory data/
scripts.dir = ""       

# running the setup file
source(file.path(scripts.dir, "setup.R"))

mc.cores = 8 # detectCores()

#####################################################
# setting directories for input and output
#####################################################
# inputs are the bam files of TT-seq, mNET-seq, various ChIP-seqs and ChIP-nexus seq
# annotation can be used as in the study or any annotation with fields "chr" for chromosome, "width" for feature length and "id" as a unique identifier (also preferably the rowname) 
# output file location can be specified

{
  # bam file locations (please provide full path for smooth running)
  tt.bam.dir = ""     # processed as in Choi et al., eLife 2021
  mnet.bam.dir = ""   # this study, please check the mapping.sh script for preprocessing and generating bam files
  cdk.bam.dir = ""    # cdk9 ChIP, please check the mapping.sh script for preprocessing
  cyc.bam.dir = ""    # cyc T1 ChIP, please check the mapping.sh script for preprocessing
  pol2.bam.dir = ""   # pol II ChIP nexus, please check the mapping.sh script for preprocessing
  
  # output folder (please provide full path for smooth running)
  out.dir = "outputs"
  dir.exists(out.dir) || dir.create(out.dir, recursive = TRUE)
  
  # other output directories
  tt.out.dir = paste0(out.dir, "/tt/")
  tt.rle.dir = paste0(tt.out.dir, "/rle/")
  tt.data.dir = paste0(tt.out.dir, "/data/")
  dir.create(tt.out.dir)
  dir.create(tt.rle.dir)
  dir.create(tt.data.dir)
  mnet.out.dir = paste0(out.dir, "/mnet/")
  mnet.rle.dir = paste0(mnet.out.dir, "/rle/")
  mnet.data.dir = paste0(mnet.out.dir, "/data/")
  dir.create(mnet.out.dir)
  dir.create(mnet.rle.dir)
  dir.create(mnet.data.dir)
  cdk.out.dir = paste0(out.dir, "/cdk9_chip/")
  cdk.rle.dir = paste0(cdk.out.dir, "/rle/")
  cdk.data.dir = paste0(cdk.out.dir, "/data/")
  dir.create(cdk.out.dir)
  dir.create(cdk.rle.dir)
  dir.create(cdk.data.dir)
  cyc.out.dir = paste0(out.dir, "/cyc_chip/")
  cyc.rle.dir = paste0(cyc.out.dir, "/rle/")
  cyc.data.dir = paste0(cyc.out.dir, "/data/")
  dir.create(cyc.out.dir)
  dir.create(cyc.rle.dir)
  dir.create(cyc.data.dir)
  pol2.out.dir = paste0(out.dir, "/pol2_nexus/")
  pol2.rle.dir = paste0(pol2.out.dir, "/rle/")
  pol2.data.dir = paste0(pol2.out.dir, "/data/")
  dir.create(pol2.out.dir)
  dir.create(pol2.rle.dir)
  dir.create(pol2.data.dir)
  
  # chr names and lengths (if you use a different version of genome, please change accordingly)
  hs.chrs.lengths = get(load("data/3_hs_chrs_names_lengths.RData"))
  hs.chrs = names(hs_chrs_lengths)
  
  # loading annotation files
  anno = get(load("data/1_GENCODE_v24_major_isoform_annotation.RData")) # annotation of gene regions
  extended.anno = get(load("data/2_GENCODE_v24_major_isoform_extended_annotation.RData")) # annotation with details of exon, intron and other genomic regions
  
  
  # creating other associated annotation data
  {
    protein = rownames(anno[which(anno$gene_type == "protein_coding"), ])
    
    # exon annotation for calculating coverages for productive initiation frequency estimation
    exon.anno = extended.anno[which(extended.anno$type == "exon"), ]
    
    single.exon.anno = exon.anno %>% 
      group_by(transcript_name) %>% 
      filter(n() == 1) %>% 
      ungroup()
    multi.exon.anno = exon.anno[-which(exon.anno$gene_name %in% single.exon.anno$gene_name), ]
    
    single.exon.anno = single.exon.anno[-which(single.exon.anno$width < 301), ]
    single.exon.anno = data.frame(resize(makeGRangesFromDataFrame(single.exon.anno, keep.extra.columns = T), fix = "end", width = (width(makeGRangesFromDataFrame(single.exon.anno, keep.extra.columns = T)) - 300)))
    
    # TSS to TSS+250 annotation for promoter proximal coverage calculation
    tss.anno = data.frame(promoters(makeGRangesFromDataFrame(anno, keep.extra.columns = T), upstream = 0, downstream = 250))
    colnames(tss.anno)[1] = "chr" # since some code snippets explicitly use this field
    rownames(tss.anno) = rownames(anno)
  }
}

#####################################################
# generating rle tracks
#####################################################
# tt-seq
{
  tt.bam.files = list.files(path = tt.bam.dir, pattern = ".bam$", full.names = F)
  
  # generating fragment mid coverage tracks
  create.fragment.mid.rle.tracks.anno(bam.files = tt.bam.files, 
                                      bam.input.folder = tt.bam.dir, 
                                      anno = anno,
                                      rle.out.folder = tt.rle.dir, 
                                      prefix = "", 
                                      human.chrs = hs.chrs,
                                      human.chrs.lengths = hs.chrs.lengths, 
                                      strand.specific = TRUE,
                                      remove.duplicates = FALSE, 
                                      size.selection = FALSE)
  
  # generating whole coverage tracks
  create.transcribed.bases.rle.tracks.anno(bam.files = tt.bam.files, 
                                           bam.input.folder = tt.bam.dir, 
                                           anno = anno,
                                           rle.out.folder = tt.rle.dir, 
                                           prefix = "", 
                                           human.chrs = hs.chrs,
                                           human.chrs.lengths = hs.chrs.lengths, 
                                           strand.specific = TRUE,
                                           remove.duplicates = FALSE, 
                                           size.selection = FALSE)
  
  # the above functions are loaded from the processing.functions.R file
}

# mnet
{
  mnet.bam.files = list.files(path = mnet.bam.dir, pattern = ".bam$", full.names = F)
  
  # generating fragment end (last Pol II position) coverage tracks
  create.fragment.end.rle.tracks.anno(bam.files = mnet.bam.files, 
                                      bam.input.folder = mnet.bam.dir, 
                                      anno = anno,
                                      rle.out.folder = mnet.rle.dir, 
                                      prefix = "", 
                                      human.chrs = hs.chrs,
                                      human.chrs.lengths = hs.chrs.lengths, 
                                      strand.specific = TRUE,
                                      remove.duplicates = FALSE, 
                                      size.selection = FALSE)
  
}

# cdk9_chip
{
  cdk.bam.files = list.files(path = cdk.bam.dir, pattern = ".bam$", full.names = F)
  
  # generating fragment mid coverage tracks
  create.fragment.mid.rle.tracks(bam.files = cdk.bam.files, 
                                 bam.input.folder = cdk.bam.dir, 
                                 anno = anno,
                                 rle.out.folder = cdk.rle.dir, 
                                 prefix = "", 
                                 human.chrs = hs.chrs,
                                 human.chrs.lengths = hs.chrs.lengths, 
                                 strand.specific = FALSE,
                                 remove.duplicates = TRUE, 
                                 size.selection = TRUE)
}

# cyc_chip
{
  cyc.bam.files = list.files(path = cyc.bam.dir, pattern = ".bam$", full.names = F) 
  
  # generating fragment mid coverage tracks
  create.fragment.mid.rle.tracks(bam.files = cyc.bam.files, 
                                 bam.input.folder = cyc.bam.dir, 
                                 anno = anno,
                                 rle.out.folder = cyc.rle.dir, 
                                 prefix = "", 
                                 human.chrs = hs.chrs,
                                 human.chrs.lengths = hs.chrs.lengths, 
                                 strand.specific = FALSE,
                                 remove.duplicates = TRUE, 
                                 size.selection = TRUE)
}

# pol2_nexus_chip
{
  pol2.bam.files = list.files(path = pol2.bam.dir, pattern = ".bam$", full.names = F) 
  
  # generating fragment mid coverage tracks
  create.nexus.rle.tracks(bam.files = pol2.bam.files, 
                          bam.input.folder = pol2.bam.dir, 
                          rle.out.folder = pol2.rle.dir, 
                          prefix = "", 
                          human.chrs  = hs.chrs, 
                          human.chrs.lengths = hs.chrs.lengths, 
                          strand.specific = FALSE, 
                          remove.duplicates = TRUE, 
                          size.selection = TRUE)
}

#####################################################
# pause positions and mnet counts for normalization
#####################################################
# transcript counts
{
  # location of antisense bias correction file (will be created in the calculate.counts step)
  mnet.antisense.bias.ratio.location = paste0(mnet.data.dir, "/antisense.bias.ratio.Rdata")
  
  # counts for the regions of interest
  calculate.counts(anno = anno, 
                   rle.location = paste0(mnet.rle.dir, "/fragment.mid.rle.tracks/"), 
                   bam.files = mnet.bam.files,
                   rle.prefix = "", 
                   chrs = hs.chrs, 
                   chrs.lengths = hs.chrs.lengths,
                   out.folder = mnet.data.dir, 
                   file.name = "unnormalized.feature.counts",
                   antisense.bias.ratio.location = tt.antisense.bias.ratio.location, 
                   antisense = T,
                   find.antisense = T)
  
  mnet.counts = get(load(paste0(mnet.data.dir, "unnormalized.counts.antisense.corrected.RData")))
  
  # size factors for normalization
  expDesign = data.frame(condition = c("0","0","12","12","24","24","72","72","96","96"), replicate = rep(c("1","2"), 5), row.names = colnames(coverage))
  dds = DESeqDataSetFromMatrix(countData = round(mnet.counts),
                               colData = expDesign,
                               design = ~ condition)
  esfObj = estimateSizeFactors(dds)
  tt.size.factors = esfObj$sizeFactor
}

# pause positions
{
  # determining max signal positions from mnet seq data
  {
    counts = list()
    
    for (bam.file in mnet.bam.files){ # loop over all samples
      print(bam.file)
      
      index.subsets = split(1:nrow(tss.anno),paste(as.character(tss.anno[,"strand"]),"_",tss.anno[,"chr"],sep = ""))
      coverage.list = list()
      
      # function to extract the position of maximum signal
      build.fragment.counts.list = function(j){
        from.transcript = strand.chr.anno[j,"start"]
        to.transcript = strand.chr.anno[j,"end"]
        
        if(strand.chr.anno[j, "strand"] == "+"){
          rle.vec = as.vector(strand.chr.fragment.counts.from.bam[from.transcript:to.transcript])
          try({rle.vec[rle.vec <= quantile(rle.vec[rle.vec != 0],0.5,na.rm = TRUE)*5] = 0},silent = TRUE) 
          max.position.transcript = which.max(rle.vec)
        }
        else{
          rle.vec = as.vector(strand.chr.fragment.counts.from.bam[to.transcript:from.transcript])
          try({rle.vec[rle.vec <= quantile(rle.vec[rle.vec != 0],0.5,na.rm = TRUE)*5] = 0},silent = TRUE) 
          max.position.transcript = which.max(rle.vec)
        }
        names(max.position.transcript) = j
        return(max.position.transcript)
      }
      
      # loop over all genes
      for (index.subset in names(index.subsets)){
        fragment.mid.track.list.chr = get(load(file.path(paste0(mnet.rle.dir, "/fragment.end.rle.tracks/",bam.file,"/",unlist(strsplit(index.subset,split = "_"))[2],".RData"))))
        
        strand.chr.fragment.counts.from.bam = fragment.mid.track.list.chr[[unlist(strsplit(index.subset,split = "_"))[1]]]
        strand.chr.anno = tss.anno[index.subsets[[index.subset]],c("start","end","width","strand")]
        
        registerDoParallel(cores = mc.cores)
        coverage.list = c(coverage.list, foreach(n = rownames(strand.chr.anno),.noexport = setdiff(ls(),c("strand.chr.fragment.counts.from.bam","strand.chr.anno","build.fragment.counts.list"))) %dopar% build.fragment.counts.list(n)) 
      }
      counts[[bam.file]] = unlist(coverage.list,recursive = TRUE,use.names = TRUE)
    }	
    
    # converting to dataframe
    pause.pos = sapply(counts,c)
    colnames(pause.pos) = bam.files
    save(pause.pos, file = paste0(mnet.data.dir, "pause.pos.RData"))
    
  }
  
  # generating promoter-proximal region annotation
  {
    # filtering
    pause.pos[pause.pos<2] = NA # if no max position is found the above function returns 1
    
    # median pause positions
    pause.pos["median"] = NA
    pause.pos["median"] = round(rowMedians(as.matrix(pause.pos[, 1:10]), na.rm = T)) # 1:10 are the columns with pause positions
    pause.pos = pause.pos[-which(is.na(pause.pos$median)), ] # removing genes with no pause positions at any time point
    
    # variation in pause positions
    pause.pos["sd"] = 0
    pause.pos["sd"] = rowSds(as.matrix(pause.pos[, 1:10]), na.rm = T)
    pause.pos[which(is.na(pause.pos$sd)), "sd"] = 0 # setting sd as 0 for genes with pause positions called at only 1 time point (sd function returns NA in this case)
    
    # filtering for less variation in positions
    pause.var.75 = rownames(pause.pos[which(pause.pos$sd < 75), ])
    pause.anno = tss.anno[pause.var.75, ]
    
    pause.anno = cbind(pause.anno[pause.var.75, ], pause.pos[pause.var.75, ]$median)
    colnames(pause.anno)[length(colnames(pause.anno))] = "pausePosition"
    
    # creating pause position +- 100 bases annotation
    pause.anno[which(pause.anno[,"strand"] == "+"),c(paste0("pause.start"),paste0("pause.end"))] = 
      cbind(pause.anno[which(pause.anno[,"strand"] == "+"),c("start")] + pause.anno[which(pause.anno[,"strand"] == "+"),"pausePosition"] - 100 + abs(pmin(pause.anno[which(pause.anno[,"strand"] == "+"),"pausePosition"] - 100,0)),
            pause.anno[which(pause.anno[,"strand"] == "+"),c("start")] + pause.anno[which(pause.anno[,"strand"] == "+"),"pausePosition"] + 100 + abs(pmin(pause.anno[which(pause.anno[,"strand"] == "+"),"pausePosition"] - 100,0)))
    pause.anno[which(pause.anno[,"strand"] == "-"),c(paste0("pause.start"),paste0("pause.end"))] = 
      cbind(pause.anno[which(pause.anno[,"strand"] == "-"),c("end")] - pause.anno[which(pause.anno[,"strand"] == "-"),"pausePosition"] - 100 - abs(pmin(pause.anno[which(pause.anno[,"strand"] == "-"),"pausePosition"] - 100,0)),
            pause.anno[which(pause.anno[,"strand"] == "-"),c("end")] - pause.anno[which(pause.anno[,"strand"] == "-"),"pausePosition"] + 100 - abs(pmin(pause.anno[which(pause.anno[,"strand"] == "-"),"pausePosition"] - 100,0)))
    pause.anno[, "start"] = pause.anno$pause.start
    pause.anno[, "end"] = pause.anno$pause.end
    pause.anno[, "width"] = pause.anno[, "end"] - pause.anno[, "start"]
  }
}

#####################################################
# counts and normalization
#####################################################
# tt-seq counts and productive initiation frequency 
{
  # transcript counts
  {
    # location of antisense bias correction file (will be created in the calculate.counts step)
    tt.antisense.bias.ratio.location = paste0(tt.data.dir, "/antisense.bias.ratio.Rdata")
    
    # counts for the regions of interest
    calculate.counts(anno = anno, 
                     rle.location = paste0(tt.rle.dir, "/fragment.mid.rle.tracks/"), 
                     bam.files = tt.bam.files,
                     rle.prefix = "", 
                     chrs = hs.chrs, 
                     chrs.lengths = hs.chrs.lengths,
                     out.folder = tt.data.dir, 
                     file.name = "unnormalized.counts",
                     antisense.bias.ratio.location = tt.antisense.bias.ratio.location, 
                     antisense = T,
                     find.antisense = T)
    
    tt.counts = get(load(paste0(tt.data.dir, "unnormalized.counts.antisense.corrected.RData")))
    
    # size factors for normalization
    expDesign = data.frame(condition = c("0","0","12","12","24","24","72","72","96","96"), replicate = rep(c("1","2"), 5), row.names = colnames(coverage))
    dds = DESeqDataSetFromMatrix(countData = round(tt.counts),
                                 colData = expDesign,
                                 design = ~ condition)
    esfObj = estimateSizeFactors(dds)
    tt.size.factors = esfObj$sizeFactor
    
    # normalizing and extracting expressed genes
    tt.counts.normalized = t(t(tt.counts)/tt.size.factors)
    tt.RPKs = tt.counts.normalized * 1e3 /anno[rownames(tt.counts.normalized), "width"]
    
    expr.genes = list()
    for (i in 1:length(colnames(tt.RPKs))) {
      expr.genes[[i]] = rownames(tt.RPKs[which(tt.RPKs[, i] > 10), ])
    }
    
    expr.genes = unique(unlist(expr.genes))
  }  
  
  # gene body coverage for productive initiation frequency estimation
  {
    calculate.counts(anno = multi.exon.anno, 
                     rle.location = paste0(tt.rle.dir, "/transcribed.bases.rle.tracks/"), 
                     bam.files = tt.bam.files,
                     rle.prefix = "", 
                     chrs = hs.chrs, 
                     chrs.lengths = hs.chrs.lengths, 
                     out.folder = tt.data.dir, 
                     file.name = "unnormalized.multi.exon.coverages",
                     antisense.bias.ratio.location = tt.antisense.bias.ratio.location, 
                     antisense = T,
                     find.antisense = F)
    calculate.counts(anno = single.exon.anno, 
                     rle.location = paste0(tt.rle.dir, "/transcribed.bases.rle.tracks/"), 
                     bam.files = tt.bam.files,
                     rle.prefix = "", 
                     chrs = hs.chrs, 
                     chrs.lengths = hs.chrs.lengths,  
                     out.folder = tt.data.dir, 
                     file.name = "unnormalized.single.exon.coverages",
                     antisense.bias.ratio.location = tt.antisense.bias.ratio.location, 
                     antisense = T,
                     find.antisense = F)
    
    # aggregating non-first exon counts
    {
      tt.coverage.extended = get(load(paste0(tt.data.dir, "unnormalized.multi.exon.coverages.antisense.corrected.RData")))
      tt.coverage.one.exon = get(load(paste0(tt.data.dir, "unnormalized.single.exon.transcript.coverages.antisense.corrected.RData")))
      
      exons.tr = exon.anno[which(!(exon.anno$exid == "")), ]$exid
      non.first.exons = exon.anno[which(!(exon.anno$exid == "") & !(exon.anno$first)), ]$exid
      tt.coverage.non.first.exons = as.data.frame(tt.coverage.extended[non.first.exons, ])
      tt.coverage.non.first.exons[, "trid"] = as.character(exon.anno[non.first.exons, ]$trid)
      
      tt.coverage.non.first.exons.aggregated = aggregate(x = tt.coverage.non.first.exons[c(mnet.bam.files)], by = tt.coverage.non.first.exons["trid"], FUN = sum)
      rownames(tt.coverage.non.first.exons.aggregated) = tt.coverage.non.first.exons.aggregated[, 1]
      tt.coverage.non.first.exons.aggregated = tt.coverage.non.first.exons.aggregated[, -1]
      lengths.non.first.exons = data.frame(row.names = non.first.exons, "length" = exon.anno[non.first.exons, ]$width, "trid" = exon.anno[non.first.exons, ]$trid)
      lengths.non.first.exons.aggregated = aggregate(x = lengths.non.first.exons["length"], by = lengths.non.first.exons["trid"], FUN = sum)
      rownames(lengths.non.first.exons.aggregated) = lengths.non.first.exons.aggregated[, 1]
      
      index = rownames(tt.coverage.non.first.exons.aggregated)
      tt.coverage.non.first.exons.aggregated.length.normalized = tt.coverage.non.first.exons.aggregated[index, ]/(lengths.non.first.exons.aggregated[index, 2])
      
      one.exon.lengths = single.exon.anno[index, ]$width
      tt.coverage.one.exon.length.normalized = tt.coverage.one.exon[index, ]/single.exon.anno[index, ]$width
      
      # genebody coverage
      tt.genebody.exon.coverage = rbind(tt.coverage.non.first.exons.aggregated, tt.coverage.one.exon)
      
      # length normalized gene body coverage
      tt.genebody.exon.coverage.length.normalized = rbind(tt.coverage.non.first.exons.aggregated.length.normalized, tt.coverage.one.exon.length.normalized)
      
      # normalizing with size factors
      tt.genebody.exon.coverage.sf.norm = t(t(tt.genebody.exon.coverage)/tt.size.factors)
      tt.genebody.exon.coverage.length.normalized.sf.norm = t(t(tt.genebody.exon.coverage.length.normalized)/tt.size.factors)
      
      # coverage per kilobases
      tt.CPKs = tt.genebody.exon.coverage.length.normalized.sf.norm * 1000
      summary(tt.CPKs)
      tt.CPKs = make.average.dataset(tt.CPKs, col.names = c(0,12,24,72,96))
      save(tt.CPKs, file = paste0(tt.data.dir, "tt.CPKs.RData"))
      
      # productive initiation frequency calculation
      initiation.rate.all = tt.genebody.exon.coverage.length.normalized.sf.norm/5 # 5 minutes of labelling time
      pif = make.average.dataset(initiation.rate.all, col.names = c(0,12,24,72,96))
      save(pif, file = paste0(tt.data.dir, "pif.RData"))
    }
  }
}

# mnet pause window coverages
{
  mnet.antisense.bias.ratio.location = paste0(mnet.data.dir, "/antisense.bias.ratio.Rdata")
  
  calculate.counts(anno = pause.anno, 
                   rle.location = paste0(mnet.rle.dir, "/fragment.end.rle.tracks/"), 
                   bam.files = mnet.bam.files,
                   rle.prefix = "", 
                   chrs = hs.chrs, 
                   chrs.lengths = hs.chrs.lengths, 
                   out.folder = mnet.data.dir, 
                   file.name = "coverage.pause.anno",
                   antisense.bias.ratio.location = mnet.antisense.bias.ratio.location, 
                   antisense = T,
                   find.antisense = T)
}

# cdk9_chip
{
  calculate.counts(anno = tss.anno, 
                   rle.location = paste0(cdk.rle.dir, "/fragment.mid.rle.tracks/"), 
                   bam.files = cdk.bam.files,
                   rle.prefix = "", 
                   chrs = hs.chrs, 
                   chrs.lengths = hs.chrs.lengths, 
                   out.folder = cdk.out.dir, 
                   file.name = "counts",
                   antisense = F)
  
  counts = paste0(cdk.out.dir, "counts.RData")
  expDesign = data.frame(condition = c("0","0","24","24","96","96"), 
                         replicate = rep(c("1","2"), 3), row.names = colnames(counts))
  dds = DESeqDataSetFromMatrix(countData = round(counts),
                               colData = expDesign,
                               design = ~ condition)
  esfObj = estimateSizeFactors(dds)
  size.factors = esfObj$sizeFactor
  save(size.factors, file = paste0(cdk.out.dir, "size.factors.RData"))
}

# cyc_chip
{
  calculate.counts(anno = tss.anno, 
                   rle.location = paste0(cyc.rle.dir, "/fragment.mid.rle.tracks/"), 
                   bam.files = cyc.bam.files,
                   rle.prefix = "", 
                   chrs = hs.chrs, 
                   chrs.lengths = hs.chrs.lengths, 
                   out.folder = cyc.out.dir, 
                   file.name = "counts",
                   antisense = F)
  counts = paste0(cyc.out.dir, "counts.RData")
  expDesign = data.frame(condition = c("0","0","24","24","96","96"), 
                         replicate = rep(c("1","2"), 3), row.names = colnames(counts))
  dds = DESeqDataSetFromMatrix(countData = round(counts),
                               colData = expDesign,
                               design = ~ condition)
  esfObj = estimateSizeFactors(dds)
  size.factors = esfObj$sizeFactor
  save(size.factors, file = paste0(cyc.out.dir, "size.factors.RData"))
}

# pol2_nexus_chip
{
  # pause coverages
  pol2.antisense.bias.ratio.location = paste0(pol2.data.dir, "/antisense.bias.ratio.Rdata") 
  antisense.bias.ratio = rep(0, length(pol2.bam.files))
  save(antisense.bias.ratio, file = pol2.antisense.bias.ratio.location)
  
  calculate.counts(anno = pause.anno, 
                   rle.location = paste0(pol2.rle.dir, "/transcribed.bases.rle.tracks/"), 
                   bam.files = pol2.bam.files,
                   rle.prefix = "", 
                   out.folder = pol2.data.dir, 
                   file.name = "coverage.pause.anno", 
                   chrs = hs.chrs, 
                   chrs.lengths = hs.chrs.lengths,
                   antisense.bias.ratio.location = pol2.antisense.bias.ratio.location, 
                   antisense = T, 
                   find.antisense = F)
  
  
  # analyzing nexus peaks and defining nexus pause regions as +- 20 bases from peaks
  {
    # peak finding function
    counts = list()
    pol2.bam.files # subset for DMSO samples only
    for (bam.file in pol2.bam.files){
      print(bam.file)
      
      coverage.list = list()
      
      build.fragment.counts.list = function(j){
        from.transcript = strand.chr.anno[j,"start"]
        to.transcript = strand.chr.anno[j,"end"]
        
        if(strand.chr.anno[j, "strand"] == "+"){
          rle.vec.plus = as.vector(fragment.mid.track.list.chr.plus[from.transcript:to.transcript])
          rle.vec.minus = as.vector(fragment.mid.track.list.chr.minus[from.transcript:to.transcript])
          max.position.transcript.plus = which.max(rle.vec.plus)
          max.position.transcript.minus = which.max(rle.vec.minus)
          max.position.transcript = round((max.position.transcript.plus + max.position.transcript.minus)/2)
        }
        else{
          rle.vec.plus = as.vector(fragment.mid.track.list.chr.plus[to.transcript:from.transcript])
          rle.vec.minus = as.vector(fragment.mid.track.list.chr.minus[to.transcript:from.transcript])
          max.position.transcript.plus = which.max(rle.vec.plus)
          max.position.transcript.minus = which.max(rle.vec.minus)
          max.position.transcript = round((max.position.transcript.plus + max.position.transcript.minus)/2)
        }
        names(max.position.transcript) = j
        return(max.position.transcript)
      }
      
      for (chr in hs.chrs){
        fragment.mid.track.list.chr = get(load(file.path(paste0(pol2.rle.dir, "/transcribed.bases.rle.tracks/",bam.file,"/",chr,".RData"))))
        fragment.mid.track.list.chr.plus = fragment.mid.track.list.chr[["+"]] 
        fragment.mid.track.list.chr.minus = fragment.mid.track.list.chr[["-"]]
        strand.chr.anno = tss.anno[which(anno$chr == chr),c("start","end","width","strand")]
        
        registerDoParallel(cores = mc.cores)
        coverage.list = c(coverage.list, foreach(n = rownames(strand.chr.anno),
                                                 .noexport = setdiff(ls(),c("fragment.mid.track.list.chr.plus","fragment.mid.track.list.chr.minus","strand.chr.anno","build.fragment.counts.list"))) %dopar% build.fragment.counts.list(n)) 
      }
      
      counts[[bam.file]] = unlist(coverage.list,recursive = TRUE,use.names = TRUE)
    }
    
    # converting peak position list into dataframe
    nexus.peaks = data.frame(sapply(counts,c))
    colnames(nexus.peaks) = pol2.bam.files
    nexus.peaks["meanPeak"] = round(rowMeans(nexus.peaks[, 1:8]))
    
    # nexus peak annotation
    {
      nexus.anno = tss.anno
      plus = rownames(tss.anno[which(tss.anno$strand == "+"),])
      minus = rownames(tss.anno[which(tss.anno$strand == "-"),])
      
      nexus.anno[plus, c("start")] = nexus.anno[plus, c("start")] + nexus.peaks[plus, "meanPeak"]
      nexus.anno[minus, c("end")] = nexus.anno[minus, c("end")] - nexus.peaks[minus, "meanPeak"]
      nexus.anno[plus, c("start")] = nexus.anno[plus, c("start")] - 20
      nexus.anno[plus, c("end")] = nexus.anno[plus, c("start")] + 40
      nexus.anno[minus, c("end")] = nexus.anno[minus, c("end")] + 20
      nexus.anno[minus, c("start")] = nexus.anno[minus, c("end")] - 40
      
      nexus.anno$width = nexus.anno$end - nexus.anno$start  
    }
    
    # nexus peak counts
    {
      calculate.counts(anno = nexus.anno, 
                       rle.location = paste0(pol2.rle.dir, "/transcribed.bases.rle.tracks/"),
                       bam.files = bam.files,
                       rle.prefix = "", 
                       out.folder = pol2.data.dir, 
                       file.name = "coverage.nexus.anno", 
                       chrs = human.chrs, 
                       chrs.lengths = human.chrs.lengths,
                       antisense.bias.ratio.location = pol2.antisense.bias.ratio.location, 
                       antisense = T, 
                       find.antisense = F)
      
    }
    
    # spike in norm for ChIP nexus (for normalizing the nexus counts before exponential fitting)
    {
      # drosophila as spike-in
      spikein.chrs = c("NT_033779.5", "NT_033778.4", "NT_037436.4", "NT_033777.3", "NC_004353.4", "NC_004354.4", "NC_024512.1")
      spikein.chrs.lengths = rep(3e7, length(spikein.chrs))
      names(spikein.chrs.lengths) = spikein.chrs
      
      create.bam.stats(bam.folder = pol2.bam.dir, 
                       bam.files = pol2.bam.files, 
                       chrs = spikein.chrs, 
                       chrs.lengths = spikein.chrs.lengths, 
                       out.file = paste0(pol2.data.dir, "/spikein.chrs.counts.RData"))
      
      bam.stats = get(load(paste0(pol2.data.dir, "/spikein.chrs.counts.RData")))
      
      pol2.spike.in.norm = colSums(bam.stats)/colSums(bam.stats)[1]
      names(pol2.spike.in.norm) = pol2.bam.files
    }
  }
  
}

#####################################################
# differential expression
#####################################################
{
  counts.table = get(load(paste0(tt.data.dir, "unnormalized.counts.antisense.corrected.RData")))
  expDesign = data.frame(condition = c("0","0","12","12","24","24", "72", "72", "96","96"), replicate = rep(c("1","2"), 5), row.names = colnames(counts.table))
  dds = DESeqDataSetFromMatrix(countData = round(counts.table),
                               colData = expDesign,
                               design = ~ condition)
  
  # dds$condition = relevel(dds$condition,ref="0")
  dds = dds[rowSums(counts(dds)) > 1,]
  dds = DESeq(dds)
  de.tt.list = list()
  de.fc.cutoff = 2
  nde.fc.cutoff = 1.5
  
  # extracting significantly differentially expressed genes
  res05 = results(dds,contrast = c("condition","0","12"),alpha=0.05)
  de.tt.list[["0_12_up"]] = rownames(subset(res05,padj < 0.05 & log2FoldChange > log2(de.fc.cutoff)))
  de.tt.list[["0_12_dn"]] = rownames(subset(res05,padj < 0.05 & log2FoldChange < -log2(de.fc.cutoff)))
  res05 = results(dds,contrast = c("condition","0","24"),alpha=0.05)
  de.tt.list[["0_24_up"]] = rownames(subset(res05,padj < 0.05 & log2FoldChange > log2(de.fc.cutoff)))
  de.tt.list[["0_24_dn"]] = rownames(subset(res05,padj < 0.05 & log2FoldChange < -log2(de.fc.cutoff)))
  res05 = results(dds,contrast = c("condition","0","72"),alpha=0.05)
  de.tt.list[["0_72_up"]] = rownames(subset(res05,padj < 0.05 & log2FoldChange > log2(de.fc.cutoff)))
  de.tt.list[["0_72_dn"]] = rownames(subset(res05,padj < 0.05 & log2FoldChange < -log2(de.fc.cutoff)))
  res05 = results(dds,contrast = c("condition","0","96"),alpha=0.05)
  de.tt.list[["0_96_up"]] = rownames(subset(res05,padj < 0.05 & log2FoldChange > log2(de.fc.cutoff)))
  de.tt.list[["0_96_dn"]] = rownames(subset(res05,padj < 0.05 & log2FoldChange < -log2(de.fc.cutoff)))
  res05 = results(dds,contrast = c("condition","12","24"),alpha=0.05)
  de.tt.list[["12_24_up"]] = rownames(subset(res05,padj < 0.05 & log2FoldChange > log2(de.fc.cutoff)))
  de.tt.list[["12_24_dn"]] = rownames(subset(res05,padj < 0.05 & log2FoldChange < -log2(de.fc.cutoff)))
  res05 = results(dds,contrast = c("condition","12","72"),alpha=0.05)
  de.tt.list[["12_72_up"]] = rownames(subset(res05,padj < 0.05 & log2FoldChange > log2(de.fc.cutoff)))
  de.tt.list[["12_72_dn"]] = rownames(subset(res05,padj < 0.05 & log2FoldChange < -log2(de.fc.cutoff)))
  res05 = results(dds,contrast = c("condition","12","96"),alpha=0.05)
  de.tt.list[["12_96_up"]] = rownames(subset(res05,padj < 0.05 & log2FoldChange > log2(de.fc.cutoff)))
  de.tt.list[["12_96_dn"]] = rownames(subset(res05,padj < 0.05 & log2FoldChange < -log2(de.fc.cutoff)))
  res05 = results(dds,contrast = c("condition","24","72"),alpha=0.05)
  de.tt.list[["24_72_up"]] = rownames(subset(res05,padj < 0.05 & log2FoldChange > log2(de.fc.cutoff)))
  de.tt.list[["24_72_dn"]] = rownames(subset(res05,padj < 0.05 & log2FoldChange < -log2(de.fc.cutoff)))
  res05 = results(dds,contrast = c("condition","24","96"),alpha=0.05)
  de.tt.list[["24_96_up"]] = rownames(subset(res05,padj < 0.05 & log2FoldChange > log2(de.fc.cutoff)))
  de.tt.list[["24_96_dn"]] = rownames(subset(res05,padj < 0.05 & log2FoldChange < -log2(de.fc.cutoff)))
  res05 = results(dds,contrast = c("condition","72","96"),alpha=0.05)
  de.tt.list[["72_96_up"]] = rownames(subset(res05,padj < 0.05 & log2FoldChange > log2(de.fc.cutoff)))
  de.tt.list[["72_96_dn"]] = rownames(subset(res05,padj < 0.05 & log2FoldChange < -log2(de.fc.cutoff)))
  save(de.tt.list, file = paste0(tt.data.dir, "de.tt.list.RData"))
}
#####################################################
# apparent pause duration
#####################################################
{
  # pif - productive initiation frequency
  # apd - apparent pause duration
  # pwc - promoter-proximal region (pause window) coverage
  
  # initiation.rate.all and pif from '# tt counts' snippet need to be loaded
  
  pwc = get(load(paste0(mnet.data.dir, "/coverage.pause.anno.antisense.corrected.Rdata")))/200
  mnet.size.factors = get(load(mnet.data.dir, "/size.factors.RData"))
  pwc = t(t(pwc)/mnet.size.factors)
  
  common = intersect(rownames(pwc), rownames(initiation.rate.all))
  
  # apparent pause duration
  pause.duration = pwc[common, mnet.bam.files]/initiation.rate.all[common, mnet.bam.files]
  apd = make.average.dataset(pause.duration, col.names = c(0,12,24,72,96))
  
  save(apd, file = paste0(mnet.out.dir, "/apd.RData"))
  
  # ifpd - dataset with robust pif and apd values for genes at all time points
  ifpd = cbind(pif[intersect(rownames(pif), rownames(apd)), ], apd[intersect(rownames(pif), rownames(apd)), ])
}
#####################################################
# classifying genes
#####################################################
{
  counts = get(load(paste0(tt.data.dir, "tt.CPKs.RData"))) # tt 0,12,24,72,96, subselect for expression cutoff
  counts.index = rownames(counts)
  
  # clustering differentially expressed genes on the basis of pif
  de.tt.list = get(load(paste0(tt.data.dir, "de.tt.list.RData"))) # differential expressed genes
  de.tt = unique(unlist(de.tt.list))
  de.tt = intersect(de.tt, expr.genes) # intersect with expression cutoff
  counts.de.tt.expressed = counts[intersect(rownames(counts), de.tt), ]
  counts = as.data.frame(counts[rownames(counts), ])
  
  # classifying into main groups based on patterns and max or min signal being at 0h or 72 and 96h (upregulated & downregulated genes)
  if(TRUE){
    clustermat = counts.de.tt.expressed
    clustermat["max"] = NA
    clustermat["min"] = NA
    
    for (i in 1:length(rownames(clustermat))) {
      clustermat[i, 6] = which.max(clustermat[i, 1:5])
      clustermat[i, 7] = which.min(clustermat[i, 1:5])
    }
    cluster.list = list()
    
    # preB
    early = rownames(clustermat[which(clustermat$max == 1 & (clustermat$min == 5 | clustermat$min == 4)), ])
    early = setdiff(early, c(de.tt.list$`72_96_dn`, de.tt.list$`0_12_dn`, de.tt.list$`12_24_dn`, de.tt.list$`24_72_dn`))
    # iMac
    late = rownames(clustermat[which(clustermat$min == 1 & (clustermat$max == 5 | clustermat$max == 4)), ])
    late = setdiff(late, c(de.tt.list$`72_96_up`, de.tt.list$`0_12_up`, de.tt.list$`12_24_up`, de.tt.list$`24_72_up`))
    
    cluster.list[["preB"]] = early
    cluster.list[["iMac"]] = late
    bend = rownames(clustermat2[which((clustermat$max == 1 | clustermat$max == 5) & (clustermat$min == 3 | clustermat$min == 2 | clustermat2$min == 4)), ])
    peak = rownames(clustermat2[which((clustermat$min == 1 | clustermat$min == 5) & (clustermat$max == 3 | clustermat$max == 2 | clustermat2$max == 4)), ])
    cluster.list[["bend"]] = setdiff(bend, unique(unlist(cluster.list)))
    cluster.list[["peak"]] = setdiff(peak, unique(unlist(cluster.list)))
    
    # subset to mrna genes
    cluster.list.protein = list()
    for (i in names(cluster.list)) {
      cluster.list.protein[[i]] = intersect(cluster.list[[i]], protein)
    }
    save(cluster.list.protein, file = paste0(tt.data.dir, "cluster.list.protein.RData"))
    
    # splitting iMac clusters on the basis of pif and apd
    data = ifpd
    data[, 1:5] = z.transform(data[,1:5]) # transforming pif values to 0-1
    data[, 6:10] = z.transform(data[,6:10]) # transforming apd values to 0-1
    clustermat = data
    nans = which(is.na(rowSums(clustermat)))
    clustermat = clustermat[-nans,]
    ifpd.tr = rownames(clustermat)
    ifpd.clustermat = ifpd[ifpd.tr, ]
    save(ifpd.clustermat, file = paste0(tt.data.dir, "ifpd.clustermat.RData"))
    
    late.ifpd = intersect(rownames(clustermat), cluster.list.protein[[2]])
    clustermat = cbind(z.transform(ifpd.clustermat[late.ifpd, 1:5]), z.transform(ifpd.clustermat[late.ifpd, 6:10]))
    
    nrcluster = 2 # number of clusters to make 
    clustercentermat = rbind()
    for (j in 1:100){
      labs = kmeans(clustermat,nrcluster)$cluster
      clustercentermat = rbind(clustercentermat,t(sapply(1:nrcluster,function(x){apply(clustermat[names(which(labs == x)),],2,median)})))
    }
    
    labs = kmeans(clustercentermat,nrcluster)$cluster
    clustercenter = t(sapply(1:nrcluster,function(x){apply(clustercentermat[which(labs == x),],2,median)}))
    labs = kmeans(clustermat,clustercenter)$cluster
    
    genecluster = list()
    for (i in 1:nrcluster){genecluster[[i]] = names(which(labs == i))}
    cluster.list.protein[["iMacI"]] = genecluster[[2]] # sometimes they cluster the other way so the indexes can change from 2 and 1 randomly
    cluster.list.protein[["iMacII"]] = genecluster[[1]]
    
    # subsetting the cluster to genes with robust values genes
    cluster.list.ifpd.protein = list()
    for (i in names(cluster.list.protein)) {
      cluster.list.ifpd.protein[[i]] = intersect((intersect(rownames(ifpd.clustermat), cluster.list.protein[[i]])), protein)
    }
    
    save(cluster.list.ifpd.protein, file = paste0(tt.data.dir, "cluster.list.protein.RData"))
  }
}

#####################################################
# exponential fitting and half-lives from ChIP nexus data
#####################################################
{
  # data
  if(TRUE){
    # load the data and normalize
    counts = get(load(paste0(pol2.data.dir, "/coverage.nexus.anno.RData")))
    counts = t(t(counts)/pol2.spike.in.norm)
    index = rownames(counts)
    for (i in 1:length(pol2.bam.files)) {
      counts[, i] = counts[, i] + counts[, i+length(pol2.bam.files)]
    }
    
    # counts dataframe now should have for the first 16 columns which are the nexus coverages but check the order of columns for the next steps
    # set up dataframes
    decay.data = data.frame(matrix(nrow = length(index), ncol = 6))
    rownames(decay.data) = index
    colnames(decay.data) = c("0.0", "0.6", "0.30", "96.0", "96.6", "96.30") # average counts 0.6 is for 0 h, 6 minutes after TRP treatment sample
    decay.data[index, "0.0"] = rowMeans(counts[index, c(1,2,5,6)]) # this assumes these column number refer to the ctrl samples of 0h, change accordingly
    decay.data[index, "0.6"] = rowMeans(counts[index, 13:14]) # 0 h, 6 min after TRP
    decay.data[index, "0.30"] = rowMeans(counts[index, 9:10]) # 0 h, 30 min after TRP
    decay.data[index, "96.0"] = rowMeans(counts[index, c(3,4,7,8)]) # 96 h, ctrl
    decay.data[index, "96.6"] = rowMeans(counts[index, 15:16]) # 96 h, 6 min after TRP
    decay.data[index, "96.30"] = rowMeans(counts[index, 11:12]) # 96 h, 30 min after TRP
    decay.data[decay.data == 0] = NA
    decay.data = decay.data[-which(is.na(rowSums(decay.data))), ]
    decay.data["0.k_LM"] = NA
    decay.data["96.k_LM"] = NA
    
    # exponential fitting
    t = c(0, 6, 30)
    for(i in 1:length(index)){
      y = unlist(decay.data[i, 1:3])
      lm.model = lm(-log(y/y[1]) ~ t)
      decay.data[i, "0.k_LM"] = coef(lm.model)[2]
      
      y = unlist(decay.data[i, 4:6])
      lm.model = lm(-log(y/y[1]) ~ t)
      decay.data[i, "96.k_LM"] = coef(lm.model)[2]
    }
    
    # half-lives
    decay.data["hl.0"] = log(2)/decay.data$`0.k_LM`
    decay.data["hl.96"] = log(2)/decay.data$`96.k_LM`
    save(decay.data, file = paste0(pol2.out.dir, "/decay.data.RData"))
  }
}
#####################################################
# termination fraction and final dataset
#####################################################
{
  # extracting 0h and 96h kinetics for analyzed transcripts
  
  tr = unlist(unique(cluster.list.ifpd.protein))
  pif.0.96 = pif[tr, ][, c(1,5)]
  colnames(pif.0.96) = paste0("pIF", colnames(pif.0.96))
  apd.0.96 = apd[tr, c(1,5)]
  colnames(apd.0.96) = paste0("aPD", colnames(apd.0.96))
  nex.pwc = decay.data[tr, ][, c(1,4)] # 0h and 96 h pause coverages for ctrl samples
  colnames(nex.pwc) = c("nPWC0", "nPWC96")
  
  data = cbind(pif.0.96[tr, ], apd.0.96[tr, ], nex.pwc[tr, ], decay.data[tr, c("hl.0", "hl.96")])
  data = data[-which(is.na(rowSums(data))), ]
  data[data < 0] = NA
  index = rownames(data)
  
  # pcr - pause clearance rate (total turnover rate at the promoter-proximal region)
  data["pcr0"] = log(2) * data[index, "0.0"]/data[index, "hl.0"]
  data["pcr96"] = log(2) * data[index, "96.0"]/data[index, "hl.96"] 
  
  # rescaling values to 0-1, note that this do not affect the distributions of the pcr or any other quantities that are estimated subsequently 
  rs = rescale(c(data$pcr0, data$pcr96))
  l = nrow(data)
  data$pcr0 = rs[1:l] # change these numbers if using any different annotation accordingly, 1865 is the nrow of the data df
  data$pcr96 = rs[(l+1):(2*l)]
  
  # rescaled pIF - productive initiation frequency
  rs = rescale(c(data$pIF0, data$pIF96)) 
  data$pIF0 = rs[1:l]
  data$pIF96 = rs[(l+1):(2*l)]
  
  # pIFcalc - elongation fraction
  data["pIF0calc"] = data[index, "pIF0"]/data[index, "pcr0"]
  data["pIF96calc"] = data[index, "pIF96"]/data[index, "pcr96"]
  rs = rescale(c(data$pIF0calc, data$pIF96calc))
  data$pIF0calc = rs[1:l]
  data$pIF96calc = rs[(l+1):(2*l)]
  
  # termination fraction
  data["drop0"] = 1 - data[index, "pIF0calc"]
  data["drop96"] = 1 - data[index, "pIF96calc"]
  rs = rescale(c(data$drop0, data$drop96))
  data$drop0 = rs[1:l]
  data$drop96 = rs[(l+1):(2*l)]
  
  # final dataframe
  save(data, file = paste0(out.dir, "/all.data.Rdata"))
}
