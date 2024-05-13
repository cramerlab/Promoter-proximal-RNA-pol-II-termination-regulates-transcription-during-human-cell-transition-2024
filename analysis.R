#####################################################
# running prerequisites
#####################################################

# loads the required libraries and custom functions used to run this script
# in case of missing libraries, please install the appropriate version for the R used
# version details of packages used in this script is available in info.txt


# setting the directory for running Setup.R as path to the location of the directory containing the CalculationUnit folder
prewd = "/usr/users/adevada/Projects/Human/"       

# Running the setup file
source(file.path(prewd,"CalculationUnit","CodeUnits","Setup","Setup.R"))
setwd(file.path(prewd))
# number of cores
mc.cores = 8 # detectCores()

# basic objects
basic.objects = c(ls(),"basic.objects")

#####################################################
# setting directories for input and output
#####################################################
# inputs are the bam files of TT-seq, mNET-seq, various ChIP-seqs and ChIP-nexus seq
# annotation can be used as in the study or any annotation with fields specified in README.md can be given
# output file location can be specified
{
  # bam file locations
  tt.bam.dir = ""
  mnet.bam.dir = ""
  cdk.bam.dir = ""    # cdk9 ChIP
  cyc.bam.dir = ""    # cyc T1 ChIP
  pol2.bam.dir = ""   # pol II chip nexus
  
  # annotation location ending with annotation file name
  anno.dir = ""
  
  # output folder
  out.dir = ""
  
  # other directories
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
  
  # associated annotation data
  human.chrs = human.chrs[1:24]
  human.chrs.lengths = human.chrs.lengths[1:24]
  anno = get(load("anno/anno.basic.complete.RData"))
  exon.anno = get(load("ProcessingScripts/anno/anno.extended.complete.RData"))
  single.exon.anno = get(load("Annotations/gencode.v24.mi.one.exon.anno.300.TTS.RData"))
  tss.anno = get(load("Annotations/gencode.v24.mi.basic.anno.TSS.250.RData"))
  trs = rownames(anno[which(anno$type == "transcript"), ]) # transcripts
  trs.exon = rownames(exon.anno[which(exon.anno$type == "transcript"), ]) # transcripts
  protein = as.character(anno[which(anno$gene_type == "protein_coding"), ]$trid)
}

#####################################################
# generating rle tracks
#####################################################
# tt
{
  tt.bam.files = paste0(tt.data.dir, "/bam.files.RData") 
  create.bam.file.names(bam.folder = tt.bam.dir, out.file = tt.bam.files)
  create.fragment.mid.rle.tracks.anno(bam.files = get(load(tt.bam.files)), 
                                      bam.input.folder = tt.bam.dir, anno = anno,
                                      rle.out.folder = tt.rle.dir, 
                                      prefix = "", 
                                      human.chrs = human.chrs,
                                      human.chrs.lengths = human.chrs.lengths, 
                                      strand.specific = TRUE,
                                      remove.duplicates = FALSE, size.selection = FALSE)
  
  create.transcribed.bases.rle.tracks.anno(bam.files = bam.files, 
                                           bam.input.folder = bam.folder, 
                                           anno = anno,
                                           rle.out.folder = rle.out.folder, 
                                           prefix = prefix, 
                                           human.chrs = human.chrs,
                                           human.chrs.lengths = human.chrs.lengths, 
                                           strand.specific = TRUE,
                                           remove.duplicates = FALSE, size.selection = FALSE)
  
}

# mnet
{
  mnet.bam.files = paste0(mnet.data.dir, "/bam.files.RData") 
  create.bam.file.names(bam.folder = mnet.bam.dir, out.file = mnet.bam.files)
  create.fragment.end.rle.tracks.anno(bam.files = get(load(mnet.bam.files)), 
                                      bam.input.folder = mnet.bam.dir, anno = anno,
                                      rle.out.folder = mnet.rle.dir, 
                                      prefix = "", 
                                      human.chrs = human.chrs,
                                      human.chrs.lengths = human.chrs.lengths, 
                                      strand.specific = TRUE,
                                      remove.duplicates = FALSE, size.selection = FALSE)
  
}

# cdk9_chip
{
  cdk.bam.files = paste0(cdk.data.dir, "/bam.files.RData") 
  create.bam.file.names(bam.folder = cdk.bam.dir, out.file = cdk.bam.files)
  create.fragment.mid.rle.tracks(bam.files = get(load(cdk.bam.files)), 
                                 bam.input.folder = cdk.bam.dir, anno = anno,
                                 rle.out.folder = cdk.rle.dir, 
                                 prefix = "", 
                                 human.chrs = human.chrs,
                                 human.chrs.lengths = human.chrs.lengths, 
                                 strand.specific = FALSE,
                                 remove.duplicates = TRUE, size.selection = TRUE)
}

# cyc_chip
{
  cyc.bam.files = paste0(cyc.data.dir, "/bam.files.RData") 
  create.bam.file.names(bam.folder = cyc.bam.dir, out.file = cyc.bam.files)
  create.fragment.mid.rle.tracks(bam.files = get(load(cyc.bam.files)), 
                                 bam.input.folder = cdk.bam.dir, anno = anno,
                                 rle.out.folder = cdk.rle.dir, 
                                 prefix = "", 
                                 human.chrs = human.chrs,
                                 human.chrs.lengths = human.chrs.lengths, 
                                 strand.specific = FALSE,
                                 remove.duplicates = TRUE, size.selection = TRUE)
}

# pol2_nexus_chip
{
  pol2.bam.files = paste0(pol2.data.dir, "/bam.files.RData") 
  create.bam.file.names(bam.folder = pol2.bam.dir, out.file = pol2.bam.files)
  create.nexus.rle.tracks(bam.files = pol2.bam.files, 
                          bam.input.folder = pol2.bam.dir, 
                          rle.out.folder = pol2.rle.dir, 
                          prefix = "", 
                          human.chrs  = human.chrs, 
                          human.chrs.lengths = human.chrs.lengths, 
                          strand.specific = FALSE, 
                          remove.duplicates = TRUE, size.selection = TRUE)
}

#####################################################
# pause positions
#####################################################
{
  # pause pos
  {
    counts = list()
    for (bam.file in mnet.bam.files){
      print(bam.file)
      index.subsets = split(1:nrow(tss.anno),paste(as.character(tss.anno[,"strand"]),"_",tss.anno[,"chr"],sep = ""))
      coverage.list = list()
      
      build.fragment.counts.list = function(j){
        from.transcript = strand.chr.anno[j,"start"]
        to.transcript = strand.chr.anno[j,"end"]
        
        if(strand.chr.anno[j, "strand"] == "+"){
          rle.vec = as.vector(strand.chr.fragment.counts.from.bam[from.transcript:to.transcript])
          try({rle.vec[rle.vec <= quantile(rle.vec[rle.vec != 0],0.5,na.rm = TRUE)*3] = 0},silent = TRUE) 
          max.position.transcript = which.max(rle.vec)
        }
        else{
          rle.vec = as.vector(strand.chr.fragment.counts.from.bam[to.transcript:from.transcript])
          try({rle.vec[rle.vec <= quantile(rle.vec[rle.vec != 0],0.5,na.rm = TRUE)*3] = 0},silent = TRUE) 
          max.position.transcript = which.max(rle.vec)
        }
        names(max.position.transcript) = j
        return(max.position.transcript)
      }
      for (index.subset in names(index.subsets)){
        fragment.mid.track.list.chr = get(load(file.path(paste0(rle.folder,bam.file,"/",unlist(strsplit(index.subset,split = "_"))[2],".RData"))))
        
        strand.chr.fragment.counts.from.bam = fragment.mid.track.list.chr[[unlist(strsplit(index.subset,split = "_"))[1]]]
        strand.chr.anno = tss.anno[index.subsets[[index.subset]],c("start","end","width","strand")]
        
        registerDoParallel(cores = mc.cores)
        coverage.list = c(coverage.list, foreach(n = rownames(strand.chr.anno),.noexport = setdiff(ls(),c("strand.chr.fragment.counts.from.bam","strand.chr.anno","build.fragment.counts.list"))) %dopar% build.fragment.counts.list(n)) 
      }
      counts[[bam.file]] = unlist(coverage.list,recursive = TRUE,use.names = TRUE)
    }	
    pause.pos = sapply(counts,c)
    colnames(pause.pos) = bam.files
    save(pause.pos, file = paste0(mnet.data.dir, "pause.pos.RData"))
    
  }
  
  # pause anno
  {
    pause.pos[pause.pos<2] = NA
    pause.pos["median"] = NA
    pause.pos["median"] = round(rowMedians(as.matrix(pause.pos[1:10]), na.rm = T)) # median pause positions for further analysis
    pause.pos = pause.pos[-which(is.na(pause.pos$median)), ]
    pause.pos["sd"] = 0
    pause.pos["sd"] = rowSds(as.matrix(pause.pos[, 1:10]), na.rm = T)
    pause.pos[which(is.na(pause.pos$sd)), "sd"] = 0
    pause.var.75 = rownames(pause.pos[which(pause.pos$sd < 75), ])
    pause.anno = tss.anno[pause.var.75, ]
    pause.pos["medianPausePosittion"] = pause.pos$median
    save(pause.pos, file = paste0(mnet.data.dir, "pause.pos.RData"))
    pause.anno = cbind(pause.anno[pause.var.75, ], pause.pos[pause.var.75, ]$medianPausePosittion)
    colnames(pause.anno)[14] = "pausePosition"
    pause.anno[which(pause.anno[,"strand"] == "+"),c(paste0("pause.start"),paste0("pause.end"))] = 
      cbind(pause.anno[which(pause.anno[,"strand"] == "+"),c("start")] + pause.anno[which(pause.anno[,"strand"] == "+"),"pausePosition"] - 100 + abs(pmin(pause.anno[which(pause.anno[,"strand"] == "+"),"pausePosition"] - 100,0)),
            pause.anno[which(pause.anno[,"strand"] == "+"),c("start")] + pause.anno[which(pause.anno[,"strand"] == "+"),"pausePosition"] + 100 + abs(pmin(pause.anno[which(pause.anno[,"strand"] == "+"),"pausePosition"] - 100,0)))
    pause.anno[which(pause.anno[,"strand"] == "-"),c(paste0("pause.start"),paste0("pause.end"))] = 
      cbind(pause.anno[which(pause.anno[,"strand"] == "-"),c("end")] - pause.anno[which(pause.anno[,"strand"] == "-"),"pausePosition"] - 100 - abs(pmin(pause.anno[which(pause.anno[,"strand"] == "-"),"pausePosition"] - 100,0)),
            pause.anno[which(pause.anno[,"strand"] == "-"),c("end")] - pause.anno[which(pause.anno[,"strand"] == "-"),"pausePosition"] + 100 - abs(pmin(pause.anno[which(pause.anno[,"strand"] == "-"),"pausePosition"] - 100,0)))
    pause.anno[, "start"] = pause.anno$pause.start
    pause.anno[, "end"] = pause.anno$pause.end
    pause.anno[, "width"] = pause.anno[, "end"] - pause.anno[, "start"]
    save(pause.anno, file = paste0(mnet.data.dir, "/pause.anno.RData"))
  }
}

#####################################################
# counts
#####################################################
# tt pif 
{
  # transcript counts
  {
    tt.antisense.bias.ratio.location = paste0(tt.data.dir, "/antisense.bias.ratio.Rdata")
    calculate.counts(anno = anno, 
                     rle.location = paste0(tt.rle.dir, "/fragment.mid.rle.tracks/"), 
                     bam.files = tt.bam.files,
                     rle.prefix = "", 
                     out.folder = tt.data.dir, 
                     file.name = "unnormalized.counts",
                     antisense.bias.ratio.location = tt.antisense.bias.ratio.location, 
                     antisense = T,
                     find.antisense = T)
    tt.counts = get(load(paste0(tt.data.dir, "unnormalized.counts.RData")))
    
    expDesign = data.frame(condition = c("0","0","12","12","24","24","30","30","36","36","72","72","96","96"), replicate = rep(c("1","2"), 7), row.names = colnames(coverage))
    dds = DESeqDataSetFromMatrix(countData = round(tt.counts),
                                 colData = expDesign,
                                 design = ~ condition)
    esfObj = estimateSizeFactors(dds)
    tt.size.factors = esfObj$sizeFactor
    save(tt.size.factors, file = paste0(tt.data.dir, "size.factors.RData"))
  }  
  
  # gene body CPKs
  {
    calculate.counts(anno = exon.anno, 
                     rle.location = paste0(tt.rle.dir, "/transcribed.bases.rle.tracks/"), 
                     bam.files = tt.bam.files,
                     rle.prefix = "", 
                     out.folder = tt.data.dir, 
                     file.name = "unnormalized.extended.coverages",
                     antisense.bias.ratio.location = tt.antisense.bias.ratio.location, 
                     antisense = T)
    calculate.counts(anno = single.exon.anno, 
                     rle.location = paste0(tt.rle.dir, "/transcribed.bases.rle.tracks/"), 
                     bam.files = tt.bam.files,
                     rle.prefix = "", 
                     out.folder = tt.data.dir, 
                     file.name = "unnormalized.single.exon.transcript.coverages",
                     antisense.bias.ratio.location = tt.antisense.bias.ratio.location, 
                     antisense = T)
    tt.coverage.extended = get(load(paste0(tt.data.dir, "unnormalized.extended.coverages.RData")))
    tt.coverage.one.exon = get(load(paste0(tt.data.dir, "unnormalized.single.exon.transcript.coverages.RData")))
    
    # aggregating non-first exon counts
    {
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
      
      tt.genebody.exon.coverage.sf.norm = t(t(tt.genebody.exon.coverage)/tt.size.factors)
      tt.genebody.exon.coverage.length.normalized.sf.norm = t(t(tt.genebody.exon.coverage.length.normalized)/tt.size.factors)
      
      tt.CPKs = tt.genebody.exon.coverage.length.normalized.sf.norm * 1000
      summary(tt.CPKs)
      tt.CPKs = make.average.dataset(tt.CPKs, col.names = c(0,12,24,30,36,72,96))
      save(tt.CPKs, file = paste0(tt.data.dir, "tt.CPKs.RData"))
      initiation.rate.all = tt.genebody.exon.coverage.length.normalized.sf.norm/5
      pif = make.average.dataset(initiation.rate.all, col.names = c(0,12,24,30,36,72,96))
      save(pif, file = paste0(tt.data.dir, "pif.RData"))
    }
  }
}

# mnet
{
  mnet.antisense.bias.ratio.location = paste0(mnet.data.dir, "/antisense.bias.ratio.Rdata")
  calculate.counts(anno = pause.anno, 
                   rle.location = paste0(mnet.rle.dir, "/fragment.end.rle.tracks/"), 
                   bam.files = mnet.bam.files,
                   rle.prefix = "", 
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
  pol2.antisense.bias.ratio.location = paste0(pol2.data.dir, "/antisense.bias.ratio.Rdata")
  calculate.counts(anno = pause.anno, 
                   rle.location = paste0(pol2.rle.dir, "/transcribed.bases.rle.tracks/"), 
                   bam.files = pol2.bam.files,
                   rle.prefix = "", 
                   out.folder = pol2.data.dir, 
                   file.name = "coverage.pause.anno", 
                   chrs = human.chrs, 
                   chrs.lengths = human.chrs.lengths,
                   antisense = T, 
                   antisense.bias.ratio.location = mnet.antisense.bias.ratio.location, 
                   find.antisense = F)
  counts = get(load(paste0(pol2.data.dir, "coverage.pause.anno.RData"))) 
  expDesign = data.frame(condition = c("0","0","96","96"), # change
                         replicate = rep(c("1","2"), 2), row.names = colnames(counts))
  dds = DESeqDataSetFromMatrix(countData = round(counts),
                               colData = expDesign,
                               design = ~ condition)
  esfObj = estimateSizeFactors(dds)
  size.factors = esfObj$sizeFactor
  save(size.factors, file = paste0(pol2.data.dir, "size.factors.RData"))
}

#####################################################
# differential expression
#####################################################
{
  counts.table = get(load(paste0(tt.data.dir, "unnormalized.counts.RData")))[, -c(7:10)]
  expDesign = data.frame(condition = c("0","0","12","12","24","24", "72", "72", "96","96"), replicate = rep(c("1","2"), 5), row.names = colnames(counts.table))
  dds = DESeqDataSetFromMatrix(countData = round(counts.table),
                               colData = expDesign,
                               design = ~ condition)
  
  #dds$condition = relevel(dds$condition,ref="0")
  dds = dds[rowSums(counts(dds)) > 1,]
  dds = DESeq(dds)
  de.tt.list = list()
  de.fc.cutoff = 2
  nde.fc.cutoff = 1.5
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
  # apd apparent pause duration
  pif = get(load(paste0(tt.out.dir, "/pif.RData")))[, -c(7,8,9,10)]
  pwc = get(load(paste0(mnet.data.dir, "/coverage.pause.anno.RData")))/200
  size.factors = get(load(mnet.data.dir, "/size.factors.RData"))
  pwc = t(t(pwc)/size.factors)
  common = intersect(rownames(pwc), rownames(init.freq))
  pause.duration = pwc[common, mnet.bam.files]/init.freq[common, mnet.bam.files]
  pause.duration = make.average.dataset(pause.duration, col.names = c(0,12,24,72,96))
  save(pause.duration.avg, file = paste0(mnet.out.dir, "apd.RData"))
  pif = get(load(paste0(tt.out.dir, "/pif.RData")))[, -c(7,8,9,10)]
  apd = get(load(paste0(mnet.out.dir, "apd.RData")))[, -c(7,8,9,10)]
  
  # ifpd - dataset with robust pif and apd values for genes at all time points
  ifpd = cbind(pif[intersect(rownames(pif), rownames(apd)), ], apd[intersect(rownames(pif), rownames(apd)), ])
}
#####################################################
# clustering
#####################################################
{
  counts = get(load(paste0(tt.data.dir, "tt.CPKs.RData")))[, -(4:5)] # tt 0,12,24,72,96, subselect for expression cutoff
  counts.index = rownames(counts)
  
  # clustering differentially expressed genes on the basis of pif
  de.tt.list = get(load(paste0(tt.data.dir, "tt.des.RData"))) # differential expressed genes
  de.tt = unique(unlist(de.tt.list))
  # de.tt = intersect(de.tt, trs) # intersect with expression cutoff
  counts.de.tt = counts[intersect(rownames(counts), de.tt), ]
  counts = as.data.frame(counts[rownames(counts), ])
  # based on patterns and significance in the changes at 72 and 96h
  if(TRUE){
    clustermat = counts.de.tt
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
    
    cluster.list.protein = list()
    for (i in names(cluster.list)) {
      cluster.list.protein[[i]] = intersect(cluster.list[[i]], protein)
    }
    save(cluster.list.protein, file = paste0(tt.data.dir, "cluster.list.protein.RData"))
    
    # Splitting late clusters on the basis of apd
    data = ifpd
    data[, 1:5] = z.transform(data[,1:5])
    data[, 6:10] = z.transform(data[,6:10])
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
    cluster.list[["iMacI"]] = genecluster[[2]]
    cluster.list[["iMacII"]] = genecluster[[1]]
    cluster.list.protein[["iMacI"]] = genecluster[[2]]
    cluster.list.protein[["iMacII"]] = genecluster[[1]]
    
    # subsetting the cluster to mRNA genes
    cluster.list.ifpd.protein = list()
    for (i in names(cluster.list.protein)) {
      cluster.list.ifpd.protein[[i]] = intersect((intersect(rownames(ifpd.clustermat), cluster.list.protein[[i]])), protein)
    }
    save(cluster.list.ifpd.protein, file = paste0(tt.data.dir, "cluster.list.protein.RData"))
    cluster.list.ifpd.protein = get(load(paste0(tt.data.dir, "cluster.list.protein.RData"))) 
  }
}

#####################################################
# exponential fitting and half-lives
#####################################################
{
  # data
  if(TRUE){
    counts = get(load(paste0(pol2.data.dir, "/coverage.pause.anno.RData")))
    norm = get(load(paste0(pol2.data.dir, "/size.factors.RData")))
    counts = t(t(counts)/norm)
    index = rownames(counts)
    decay.data = data.frame(matrix(nrow = length(index), ncol = 6))
    rownames(decay.data) = index
    colnames(decay.data) = c("0.0", "0.6", "0.30", "96.0", "96.6", "96.30") # average counts
    decay.data[index, "0.0"] = rowMeans(counts[index, c(1,2,5,6)])
    decay.data[index, "0.6"] = rowMeans(counts[index, 13:14])
    decay.data[index, "0.30"] = rowMeans(counts[index, 9:10])
    decay.data[index, "96.0"] = rowMeans(counts[index, c(3,4,7,8)])
    decay.data[index, "96.6"] = rowMeans(counts[index, 15:16])
    decay.data[index, "96.30"] = rowMeans(counts[index, 11:12])
    decay.data[decay.data == 0] = NA
    decay.data = decay.data[-which(is.na(rowSums(decay.data))), ]
    decay.data["0.k_LM"] = NA
    decay.data["96.k_LM"] = NA
    t = c(0, 6, 30)
    for(i in 1:length(index)){
      y = unlist(decay.data[i, 1:3])
      lm.model = lm(-log(y/y[1]) ~ t)
      decay.data[i, "0.k_LM"] = coef(lm.model)[2]
      
      y = unlist(decay.data[i, 4:6])
      lm.model = lm(-log(y/y[1]) ~ t)
      decay.data[i, "96.k_LM"] = coef(lm.model)[2]
    }
    decay.data["hl.0"] = log(2)/decay.data$`0.k_LM`
    decay.data["hl.96"] = log(2)/decay.data$`96.k_LM`
    save(decay.data, file = paste0(pol2.out.dir, "/decay.data.RData"))
  }
}
#####################################################
# termination fraction
#####################################################
{
  tr = unlist(unique(cluster.list.ifpd.protein))
  data = cbind(decay.data[tr, ], pif[tr, ])
  data = data[-which(is.na(rowSums(data))), ]
  data[data < 0] = NA
  index = rownames(data)
  
  # pcr - pause clearance rate(total turnover rate at the promoter-proximal region)
  data["pcr0"] = data[index, "0.0"]/data[index, "hl.0"]
  data["pcr96"] = data[index, "96.0"]/data[index, "hl.96"] 
  rs = rescale(c(data$pcr0, data$pcr96))
  data$pcr0 = rs[1:1865]
  data$pcr96 = rs[1866:3730]
  rs = rescale(c(data$pIF0, data$pIF96)) # check
  data$pIF0 = rs[1:1865]
  data$pIF96 = rs[1866:3730]
  # pIFcalc - elongation fraction
  data["pIF0calc"] = data[index, "pIF0"]/data[index, "pcr0"]
  data["pIF96calc"] = data[index, "pIF96"]/data[index, "pcr96"]
  rs = rescale(c(data$pIF0calc, data$pIF96calc))
  data$pIF0calc = rs[1:1865]
  data$pIF96calc = rs[1866:3730]
  # drop - termination fraction
  data["drop0"] = 1 - data[index, "pIF0calc"]
  data["drop96"] = 1 - data[index, "pIF96calc"]
  rs = rescale(c(data$drop0, data$drop96))
  data$drop0 = rs[1:1865]
  data$drop96 = rs[1866:3730]
  save(data, file = paste0(pol2.out.dir, "/hl.term.RData"))
}

