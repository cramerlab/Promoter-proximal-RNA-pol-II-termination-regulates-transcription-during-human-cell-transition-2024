create.bam.stats = function(bam.folder, bam.files, chrs, chrs.lengths, out.file){
  bam.lengths = c()
  bam.stats = data.frame(matrix(nrow = length(chrs)+2, ncol = length(bam.files)))
  rownames(bam.stats) = c(chrs, "totalReadsMapped", "duplicateReads")
  colnames(bam.stats) = bam.files
  
  for(bam.file in bam.files){
    get.bam.chr.lengths = function(which.chr){
      param = ScanBamParam(which=GRanges(seqnames = which.chr,ranges = IRanges(0,chrs.lengths[which.chr])))
      bam = readGAlignmentPairs(file = file.path(paste0(bam.folder, bam.file)),param = param)
      bam.chr = bam[which(seqnames(bam) == which.chr)] # Remove non-concordant mappings
      bam.chr.length = length(bam.chr)
      bam.dups.length = length(bam[duplicated(paste(start(left(bam)),end(right(bam)))),])
      bam.length.vector = c(bam.chr.length, bam.dups.length)
      names(bam.length.vector) = c(which.chr, "dup")
      return(bam.length.vector)
    }
    bam.chr.length.list = foreach(n = chrs,.noexport = setdiff(ls(), c("chrs.lengths"))) %dopar% get.bam.chr.lengths(n)
    bam.stats[chrs, bam.file] = unlist(bam.chr.length.list)[chrs]
    bam.stats["totalReadsMapped", bam.file] = sum(unlist(bam.chr.length.list)[chrs])
    bam.stats["duplicateReads", bam.file] = sum(unlist(bam.chr.length.list)) - sum(unlist(bam.chr.length.list)[chrs])
    
  }
  save(bam.stats, file = outfile)
}

create.bam.file.names = function(bam.folder, out.file){
  bam.files = list.files(bam.folder, pattern = "\\.bam$", full.names = FALSE)
  if(all(grepl(".bam", bam.files))){
    bam.files = gsub('.{4}$', '', bam.files)
  }
  save(bam.files, file = out.file)
}

create.seq.depth.from.bam = function(bam.folder, bam.files, human.chrs, human.chrs.lengths, out.file){
  bam.lengths = c()
  for(bam.file in bam.files){
    get.bam.chr.lengths = function(which.chr){
      param = ScanBamParam(which=GRanges(seqnames = which.chr,ranges = IRanges(0,human.chrs.lengths[which.chr])))
      bam = readGAlignmentPairs(file = file.path(paste0(bam.folder, bam.file,".bam")),param = param)
      bam = bam[start(left(bam)) <= end(right(bam))]
      bam = bam[which(seqnames(bam) == which.chr)] # Remove non-concordant mappings
      bam.length = length(bam)
      return(bam.length)
    }
    bam.chr.length.list = foreach(n = human.chrs,.noexport = setdiff(ls(), c("human.chrs.lengths"))) %dopar% get.bam.chr.lengths(n)
    bam.lengths = c(bam.lengths, sum(unlist(bam.chr.length.list)))
  }
  names(bam.lengths) = bam.files
  save(bam.lengths, file = out.file)
}

create.fragment.mid.rle.tracks = function(bam.files, bam.input.folder, rle.out.folder, prefix, human.chrs, human.chrs.lengths, strand.specific = TRUE, remove.duplicates = FALSE, size.selection = FALSE, strand.mode = 1) {
  dir.create(file.path(rle.out.folder, "fragment.mid.rle.tracks"))
  
  for (bam.file in bam.files){
    print(bam.file)
    dir.create(file.path(rle.out.folder, "fragment.mid.rle.tracks", bam.file))
    
    registerDoParallel(cores = mc.cores)
    build.fragment.mid.track.list = function(which.chr){
      print("Starting")
      print(which.chr)
      fragment.mid.track.list.chr = list()
      param = ScanBamParam(which=GRanges(seqnames = which.chr,ranges = IRanges(0,human.chrs.lengths[which.chr])))
      bam = readGAlignmentPairs(file = file.path(paste0(bam.input.folder, bam.file,".bam")),param = param, strandMode = strand.mode)
      bam = bam[start(left(bam)) <= end(right(bam))]
      bam = bam[which(seqnames(bam) == which.chr)] # Remove non-concordant mappings
      if(length(bam) != 0){
        if(size.selection) {
          bam = bam[(end(right(bam)) - start(left(bam))) <= 500] # Size selection
        }
        if(remove.duplicates) {
          bam = bam[!duplicated(paste(start(left(bam)),end(right(bam)))),] # Remove duplicates
        }
        
        starts = end(left(bam[strand(bam) == "+"]))
        fragment.mid.track.list.chr[["+"]] = 0
        if(length(starts) != 0){
          inner.mate.distance = start(right(bam[strand(bam) == "+"])) - starts
          mid.points = round(starts + (inner.mate.distance)/2)
          fragment.count = rep(1,length(inner.mate.distance))
          fragment.counts.aggregate = aggregate(fragment.count,list(id = mid.points),function(x){sum(x,na.rm = TRUE)})
          fragment.counts.vec = Rle(0,human.chrs.lengths[which.chr])
          fragment.counts.vec[fragment.counts.aggregate[,"id"]] = fragment.counts.aggregate[,"x"]
          fragment.mid.track.list.chr[["+"]] = fragment.counts.vec
        }
        
        starts = end(left(bam[strand(bam) == "-"]))
        fragment.mid.track.list.chr[["-"]] = 0
        if(length(starts) != 0){
          inner.mate.distance = start(right(bam[strand(bam) == "-"])) - starts
          mid.points = round(starts + (inner.mate.distance)/2)
          fragment.count = rep(1,length(inner.mate.distance))
          fragment.counts.aggregate = aggregate(fragment.count,list(id = mid.points),function(x){sum(x,na.rm = TRUE)})
          fragment.counts.vec = Rle(0,human.chrs.lengths[which.chr])
          fragment.counts.vec[fragment.counts.aggregate[,"id"]] = fragment.counts.aggregate[,"x"]
          fragment.mid.track.list.chr[["-"]] = fragment.counts.vec
        }
        
        if(!strand.specific){
          fragment.mid.track.list.chr[["+"]] = fragment.mid.track.list.chr[["+"]] + fragment.mid.track.list.chr[["-"]]
          fragment.mid.track.list.chr[["-"]] = fragment.mid.track.list.chr[["+"]]
        }
        
        save(fragment.mid.track.list.chr,file = file.path(rle.out.folder, "fragment.mid.rle.tracks", bam.file, paste0(prefix, which.chr, ".RData")))
        return()
      }else{
        fragment.mid.track.list.chr[["+"]] = 0
        if(!strand.specific){
          fragment.mid.track.list.chr[["+"]] = 0
          fragment.mid.track.list.chr[["-"]] = 0
        }
        save(fragment.mid.track.list.chr,file = file.path(rle.out.folder, "fragment.mid.rle.tracks", bam.file, paste0(prefix, which.chr, ".RData")))
        return()
      }
      print("Ending")
      print(which.chr)
      }
    fragment.mid.track.list = foreach(n = human.chrs,.noexport = setdiff(ls(), c("human.chrs.lengths"))) %dopar% build.fragment.mid.track.list(n)
  }
  
  gc()
}

create.fragment.end.rle.tracks = function(bam.files, bam.input.folder, rle.out.folder, prefix, human.chrs, human.chrs.lengths, strand.specific = TRUE, remove.duplicates = FALSE, size.selection = FALSE, strand.mode = 1) {
  dir.create(file.path(rle.out.folder, "fragment.end.rle.tracks"))
  
  for (bam.file in bam.files){
    dir.create(file.path(rle.out.folder, "fragment.end.rle.tracks", bam.file))
    
    registerDoParallel(cores = mc.cores)
    build.fragment.end.track.list = function(which.chr){
      fragment.end.track.list.chr = list()
      param = ScanBamParam(which=GRanges(seqnames = which.chr,ranges = IRanges(0,human.chrs.lengths[which.chr])))
      bam = readGAlignmentPairs(file = file.path(paste0(bam.input.folder, bam.file,".bam")),param = param)
      bam = bam[start(left(bam)) <= end(right(bam))]
      bam = bam[which(seqnames(bam) == which.chr)] # Remove non-concordant mappings
      if(size.selection) {
        bam = bam[(end(right(bam)) - start(left(bam))) <= 500] # Size selection
      }
      if(remove.duplicates) {
        bam = bam[!duplicated(paste(start(left(bam)),end(right(bam)))),] # Remove duplicates
      }
      
      end.points = end(right(bam[strand(bam) == "+"]))
      fragment.count = rep(1,length(end.points))
      fragment.counts.aggregate = aggregate(fragment.count,list(id = end.points),function(x){sum(x,na.rm = TRUE)})
      fragment.counts.vec = Rle(0,human.chrs.lengths[which.chr])
      fragment.counts.vec[fragment.counts.aggregate[,"id"]] = fragment.counts.aggregate[,"x"]
      fragment.end.track.list.chr[["+"]] = fragment.counts.vec
      
      end.points = start(left(bam[strand(bam) == "-"]))
      fragment.count = rep(1,length(end.points))
      fragment.counts.aggregate = aggregate(fragment.count,list(id = end.points),function(x){sum(x,na.rm = TRUE)})
      fragment.counts.vec = Rle(0,human.chrs.lengths[which.chr])
      fragment.counts.vec[fragment.counts.aggregate[,"id"]] = fragment.counts.aggregate[,"x"]
      fragment.end.track.list.chr[["-"]] = fragment.counts.vec
      if(!strand.specific){
        fragment.end.track.list.chr[["+"]] = fragment.end.track.list.chr[["+"]] + fragment.end.track.list.chr[["-"]]
        fragment.end.track.list.chr[["-"]] = fragment.end.track.list.chr[["+"]]
      }
      
      save(fragment.end.track.list.chr,file = file.path(rle.out.folder, "fragment.end.rle.tracks", bam.file ,paste0(prefix, which.chr, ".RData")))
      return()
    }
    
    fragment.end.track.list = foreach(n = human.chrs,.noexport = setdiff(ls(),c("human.chrs.lengths"))) %dopar% build.fragment.end.track.list(n)
  }
  
  gc()
}

create.fragment.size.rle.tracks = function(bam.files, bam.input.folder, rle.out.folder, prefix, human.chrs, human.chrs.lengths, strand.specific = TRUE, remove.duplicates = FALSE, size.selection = FALSE, strand.mode = 1) {
  dir.create(file.path(rle.out.folder, "fragment.size.rle.tracks"))
  
  for (bam.file in bam.files){
    dir.create(file.path(rle.out.folder, "fragment.size.rle.tracks", bam.file))
    
    registerDoParallel(cores = mc.cores)
    build.fragment.size.track.list = function(which.chr){
      fragment.size.track.list.chr = list()
      param = ScanBamParam(which=GRanges(seqnames = which.chr,ranges = IRanges(0,human.chrs.lengths[which.chr])))
      bam = readGAlignmentPairs(file = file.path(paste0(bam.input.folder, bam.file,".bam")),param = param)
      bam = bam[start(left(bam)) <= end(right(bam))]
      bam = bam[which(seqnames(bam) == which.chr)] # Remove non-concordant mappings
      if(size.selection) {
        bam = bam[(end(right(bam)) - start(left(bam))) <= 500] # Size selection
      }
      if(remove.duplicates) {
        bam = bam[!duplicated(paste(start(left(bam)),end(right(bam)))),] # Remove duplicates
      }
      
      left.cigars = cigar(left(bam))
      left.cigar.ops = explodeCigarOps(left.cigars)
      right.cigars = cigar(right(bam))
      right.cigar.ops = explodeCigarOps(right.cigars)
      bam = bam[(!sapply(left.cigar.ops,function(x){any(x == "N")}) & !sapply(right.cigar.ops,function(x){any(x == "N")})),]
      
      strand.bam = bam[strand(bam) == "+"]
      fragment.size.vec = Rle(0,human.chrs.lengths[which.chr])
      if (length(strand.bam) > 0){
        left.cigars = cigar(left(strand.bam))
        left.cigar.ops = explodeCigarOps(left.cigars)
        left.cigar.lengths = sapply(explodeCigarOpLengths(left.cigars,ops=c("M","I","=","X")),sum) # sapply(explodeCigarOpLengths(left.cigars,ops=c("M","I","S","=","X")),sum)
        right.cigars = cigar(right(strand.bam))
        right.cigar.ops = explodeCigarOps(right.cigars)
        right.cigar.lengths = sapply(explodeCigarOpLengths(right.cigars,ops=c("M","I","=","X")),sum) # sapply(explodeCigarOpLengths(right.cigars,ops=c("M","I","S","=","X")),sum)
        starts = end(left(strand.bam))
        inner.mate.distance = start(right(strand.bam)) - starts
        mid.points = round(starts + (inner.mate.distance)/2)			
        fragment.size = inner.mate.distance + left.cigar.lengths + right.cigar.lengths
        if (length(fragment.size[fragment.size >= 0 & fragment.size <= 500]) > 0){
          mid.points = mid.points[fragment.size >= 0 & fragment.size <= 500]
          fragment.size = fragment.size[fragment.size >= 0 & fragment.size <= 500]
          try({fragment.size.aggregate = aggregate(fragment.size,list(id = mid.points),function(x){mean(x,na.rm = TRUE)});fragment.size.vec[fragment.size.aggregate[,"id"]] = fragment.size.aggregate[,"x"]},silent = TRUE)
        }
      }
      fragment.size.track.list.chr[["+"]] = fragment.size.vec
      
      strand.bam = bam[strand(bam) == "-"]
      fragment.size.vec = Rle(0,human.chrs.lengths[which.chr])
      if (length(strand.bam) > 0){
        left.cigars = cigar(left(strand.bam))
        left.cigar.ops = explodeCigarOps(left.cigars)
        left.cigar.lengths = sapply(explodeCigarOpLengths(left.cigars,ops=c("M","I","=","X")),sum) # sapply(explodeCigarOpLengths(left.cigars,ops=c("M","I","S","=","X")),sum)
        right.cigars = cigar(right(strand.bam))
        right.cigar.ops = explodeCigarOps(right.cigars)
        right.cigar.lengths = sapply(explodeCigarOpLengths(right.cigars,ops=c("M","I","=","X")),sum) # sapply(explodeCigarOpLengths(right.cigars,ops=c("M","I","S","=","X")),sum)
        starts = end(left(strand.bam))
        inner.mate.distance = start(right(strand.bam)) - starts
        mid.points = round(starts + (inner.mate.distance)/2)			
        fragment.size = inner.mate.distance + left.cigar.lengths + right.cigar.lengths
        if (length(fragment.size[fragment.size >= 0 & fragment.size <= 500]) > 0){
          mid.points = mid.points[fragment.size >= 0 & fragment.size <= 500]
          fragment.size = fragment.size[fragment.size >= 0 & fragment.size <= 500]
          try({fragment.size.aggregate = aggregate(fragment.size,list(id = mid.points),function(x){mean(x,na.rm = TRUE)});fragment.size.vec[fragment.size.aggregate[,"id"]] = fragment.size.aggregate[,"x"]},silent = TRUE)
        }
      }
      fragment.size.track.list.chr[["-"]] = fragment.size.vec
      if(!strand.specific){
        fragment.size.track.list.chr[["+"]] = fragment.size.track.list.chr[["+"]] + fragment.size.track.list.chr[["-"]]
        fragment.size.track.list.chr[["-"]] = fragment.size.track.list.chr[["+"]]
      }
      
      save(fragment.size.track.list.chr,file = file.path(rle.out.folder, "fragment.size.rle.tracks", bam.file, paste0(prefix, which.chr, ".RData")))
      return()
    }
    
    fragment.size.track.list = foreach(n = human.chrs,.noexport = setdiff(ls(),c("human.chrs.lengths"))) %dopar% build.fragment.size.track.list(n)
  }
  
  gc()
}

create.transcribed.bases.rle.tracks = function(bam.files, bam.input.folder, rle.out.folder, prefix, human.chrs, human.chrs.lengths, strand.specific = TRUE, remove.duplicates = FALSE, size.selection = FALSE, strand.mode = 1) {
  dir.create(file.path(rle.out.folder, "transcribed.bases.rle.tracks"))
  for (bam.file in bam.files){
    dir.create(file.path(file.path(rle.out.folder, "transcribed.bases.rle.tracks", bam.file)))
    
    registerDoParallel(cores = mc.cores)
    build.transcribed.bases.coverage.track.list = function(which.chr){
      transcribed.bases.coverage.track.list.chr = list()
      param = ScanBamParam(which=GRanges(seqnames = which.chr,ranges = IRanges(0,human.chrs.lengths[which.chr])))
      bam = readGAlignmentPairs(file = file.path(paste0(bam.input.folder, bam.file,".bam")),param = param, strandMode = strand.mode)
      bam = bam[start(left(bam)) <= end(right(bam))]
      bam = bam[which(seqnames(bam) == which.chr)] # Remove non-concordant mappings
      if(length(bam) != 0){
        if(size.selection) {
          bam = bam[(end(right(bam)) - start(left(bam))) <= 500] # Size selection
        }
        if(remove.duplicates) {
          bam = bam[!duplicated(paste(start(left(bam)),end(right(bam)))),] # Remove duplicates
        }
        starts = end(left(bam[strand(bam) == "+"]))
        transcribed.bases.coverage.track.list.chr[["+"]] = 0
        if(length(starts) != 0){
          rle.vec = Rle(0,human.chrs.lengths[which.chr])
          coverage.vec = coverage(GRanges(seqnames = which.chr,ranges = IRanges(start = start(left(bam[strand(bam) == "+"])),end = end(right(bam[strand(bam) == "+"])))))[[which.chr]]
          rle.vec[1:length(coverage.vec)] = coverage.vec
          transcribed.bases.coverage.track.list.chr[["+"]] = rle.vec
        }
        
        starts = end(left(bam[strand(bam) == "-"]))
        transcribed.bases.coverage.track.list.chr[["-"]] = 0
        if(length(starts) != 0){
          rle.vec = Rle(0,human.chrs.lengths[which.chr])
          coverage.vec = coverage(GRanges(seqnames = which.chr,ranges = IRanges(start = start(left(bam[strand(bam) == "-"])),end = end(right(bam[strand(bam) == "-"])))))[[which.chr]]
          rle.vec[1:length(coverage.vec)] = coverage.vec
          transcribed.bases.coverage.track.list.chr[["-"]] = rle.vec
        }
        if(!strand.specific) {
          transcribed.bases.coverage.track.list.chr[["+"]] = transcribed.bases.coverage.track.list.chr[["+"]] + transcribed.bases.coverage.track.list.chr[["-"]]
          transcribed.bases.coverage.track.list.chr[["-"]] = transcribed.bases.coverage.track.list.chr[["+"]]
        }
        
        
        save(transcribed.bases.coverage.track.list.chr,file = file.path(rle.out.folder, "transcribed.bases.rle.tracks", bam.file, paste0(prefix, which.chr, ".RData")))
        return()
      }else{
        transcribed.bases.coverage.track.list.chr[["+"]] = 0
        if(!strand.specific){
          transcribed.bases.coverage.track.list.chr[["+"]] = 0
          transcribed.bases.coverage.track.list.chr[["-"]] = 0
        }
        save(transcribed.bases.coverage.track.list.chr,file = file.path(rle.out.folder, "transcribed.bases.rle.tracks", bam.file, paste0(prefix, which.chr, ".RData")))
        return()
      }
      
    }
    
    transcribed.bases.coverage.track.list = foreach(n = human.chrs,.noexport = setdiff(ls(),c("human.chrs.lengths"))) %dopar% build.transcribed.bases.coverage.track.list(n)
  }
}

create.fragment.mid.rle.tracks.anno = function(bam.files, bam.input.folder, rle.out.folder, prefix, anno, human.chrs, human.chrs.lengths, strand.specific = TRUE, remove.duplicates = FALSE, size.selection = FALSE, strand.mode = 1) {
  dir.create(file.path(rle.out.folder, "fragment.mid.rle.tracks"))
  anno[, c("start")] = anno[, c("start")] - 500
  anno[, c("end")] = anno[, c("end")] + 500
  anno[, "width"] = anno[, "end"] - anno[, "start"]
  for (bam.file in bam.files){
    dir.create(file.path(rle.out.folder, "fragment.mid.rle.tracks", bam.file))
    
    registerDoParallel(cores = mc.cores)
    build.fragment.mid.track.list = function(which.chr){
      fragment.mid.track.list.chr = list()
      param = ScanBamParam(which=GRanges(seqnames = which.chr,ranges = IRanges(0,human.chrs.lengths[which.chr])))
      bam = readGAlignmentPairs(file = file.path(paste0(bam.input.folder, bam.file,".bam")),param = param, strandMode = strand.mode)
      bam = bam[start(left(bam)) <= end(right(bam))]
      bam = bam[which(seqnames(bam) == which.chr)] # Remove non-concordant mappings
      if(length(bam) != 0){
        if(size.selection) {
          bam = bam[(end(right(bam)) - start(left(bam))) <= 500] # Size selection
        }
        if(remove.duplicates) {
          bam = bam[!duplicated(paste(start(left(bam)),end(right(bam)))),] # Remove duplicates
        }
        bam.chr.ranges = GRanges(seqnames = which.chr,strand = strand(bam),ranges = IRanges(start = start(left(bam)), end = end(right(bam))))
        
        strand.chr.anno = makeGRangesFromDataFrame(anno[which(anno$chr == which.chr),c("chr","start","end","width", "strand")])
        bam.anno.overlaps = findOverlaps(bam.chr.ranges, strand.chr.anno, maxgap = 0L, minoverlap = 1L, type = "within", select = "all", ignore.strand = T)
        bam = bam[unique(queryHits(bam.anno.overlaps))]
        
        starts = end(left(bam[strand(bam) == "+"]))
        fragment.mid.track.list.chr[["+"]] = 0
        if(length(starts) != 0){
          inner.mate.distance = start(right(bam[strand(bam) == "+"])) - starts
          mid.points = round(starts + (inner.mate.distance)/2)
          fragment.count = rep(1,length(inner.mate.distance))
          fragment.counts.aggregate = aggregate(fragment.count,list(id = mid.points),function(x){sum(x,na.rm = TRUE)})
          fragment.counts.vec = Rle(0,human.chrs.lengths[which.chr])
          fragment.counts.vec[fragment.counts.aggregate[,"id"]] = fragment.counts.aggregate[,"x"]
          fragment.mid.track.list.chr[["+"]] = fragment.counts.vec
        }
        
        starts = end(left(bam[strand(bam) == "-"]))
        fragment.mid.track.list.chr[["-"]] = 0
        if(length(starts) != 0){
          inner.mate.distance = start(right(bam[strand(bam) == "-"])) - starts
          mid.points = round(starts + (inner.mate.distance)/2)
          fragment.count = rep(1,length(inner.mate.distance))
          fragment.counts.aggregate = aggregate(fragment.count,list(id = mid.points),function(x){sum(x,na.rm = TRUE)})
          fragment.counts.vec = Rle(0,human.chrs.lengths[which.chr])
          fragment.counts.vec[fragment.counts.aggregate[,"id"]] = fragment.counts.aggregate[,"x"]
          fragment.mid.track.list.chr[["-"]] = fragment.counts.vec
        }
        
        if(!strand.specific){
          fragment.mid.track.list.chr[["+"]] = fragment.mid.track.list.chr[["+"]] + fragment.mid.track.list.chr[["-"]]
          fragment.mid.track.list.chr[["-"]] = fragment.mid.track.list.chr[["+"]]
        }
        
        save(fragment.mid.track.list.chr,file = file.path(rle.out.folder, "fragment.mid.rle.tracks", bam.file, paste0(prefix, which.chr, ".RData")))
        return()
      }else{
        fragment.mid.track.list.chr[["+"]] = 0
        if(!strand.specific){
          fragment.mid.track.list.chr[["+"]] = 0
          fragment.mid.track.list.chr[["-"]] = 0
        }
        save(fragment.mid.track.list.chr,file = file.path(rle.out.folder, "fragment.mid.rle.tracks", bam.file, paste0(prefix, which.chr, ".RData")))
        return()
      }
      
    }
    
    fragment.mid.track.list = foreach(n = human.chrs,.noexport = setdiff(ls(), c("human.chrs.lengths"))) %dopar% build.fragment.mid.track.list(n)
  }
  gc()
}

create.fragment.end.rle.tracks.anno = function(bam.files, bam.input.folder, rle.out.folder, prefix, anno, human.chrs, human.chrs.lengths, strand.specific = TRUE, remove.duplicates = FALSE, size.selection = FALSE, strand.mode = 1) {
  dir.create(file.path(rle.out.folder, "fragment.end.rle.tracks"))
  anno[, c("start")] = anno[, c("start")] - 500
  anno[, c("end")] = anno[, c("end")] + 500
  anno[, "width"] = anno[, "end"] - anno[, "start"]
  
  for (bam.file in bam.files){
    dir.create(file.path(rle.out.folder, "fragment.end.rle.tracks", bam.file))
    
    registerDoParallel(cores = mc.cores)
    build.fragment.end.track.list = function(which.chr){
      fragment.end.track.list.chr = list()
      param = ScanBamParam(which=GRanges(seqnames = which.chr,ranges = IRanges(0,human.chrs.lengths[which.chr])))
      bam = readGAlignmentPairs(file = file.path(paste0(bam.input.folder, bam.file,".bam")),param = param)
      bam = bam[start(left(bam)) <= end(right(bam))]
      bam = bam[which(seqnames(bam) == which.chr)] # Remove non-concordant mappings
      if(size.selection) {
        bam = bam[(end(right(bam)) - start(left(bam))) <= 500] # Size selection
      }
      if(remove.duplicates) {
        bam = bam[!duplicated(paste(start(left(bam)),end(right(bam)))),] # Remove duplicates
      }
      
      bam.chr.ranges = GRanges(seqnames = which.chr,strand = strand(bam),ranges = IRanges(start = start(left(bam)), end = end(right(bam))))
      
      strand.chr.anno = makeGRangesFromDataFrame(anno[which(anno$chr == which.chr),c("chr","start","end","width", "strand")])
      bam.anno.overlaps = findOverlaps(bam.chr.ranges, strand.chr.anno, maxgap = 0L, minoverlap = 1L, type = "within", select = "all", ignore.strand = T)
      bam = bam[unique(queryHits(bam.anno.overlaps))]
      
      end.points = end(right(bam[strand(bam) == "+"]))
      fragment.count = rep(1,length(end.points))
      fragment.counts.aggregate = aggregate(fragment.count,list(id = end.points),function(x){sum(x,na.rm = TRUE)})
      fragment.counts.vec = Rle(0,human.chrs.lengths[which.chr])
      fragment.counts.vec[fragment.counts.aggregate[,"id"]] = fragment.counts.aggregate[,"x"]
      fragment.end.track.list.chr[["+"]] = fragment.counts.vec
      
      end.points = start(left(bam[strand(bam) == "-"]))
      fragment.count = rep(1,length(end.points))
      fragment.counts.aggregate = aggregate(fragment.count,list(id = end.points),function(x){sum(x,na.rm = TRUE)})
      fragment.counts.vec = Rle(0,human.chrs.lengths[which.chr])
      fragment.counts.vec[fragment.counts.aggregate[,"id"]] = fragment.counts.aggregate[,"x"]
      fragment.end.track.list.chr[["-"]] = fragment.counts.vec
      if(!strand.specific) {
        fragment.end.track.list.chr[["+"]] = fragment.end.track.list.chr[["+"]] + fragment.end.track.list.chr[["-"]]
        fragment.end.track.list.chr[["-"]] = fragment.end.track.list.chr[["+"]]
      }
      
      save(fragment.end.track.list.chr,file = file.path(rle.out.folder, "fragment.end.rle.tracks", bam.file ,paste0(prefix, which.chr, ".RData")))
      return()
    }
    
    fragment.end.track.list = foreach(n = human.chrs,.noexport = setdiff(ls(),c("human.chrs.lengths"))) %dopar% build.fragment.end.track.list(n)
  }
  
  gc()
}

create.fragment.size.rle.tracks.anno = function(bam.files, bam.input.folder, rle.out.folder, prefix, anno, human.chrs, human.chrs.lengths, strand.specific = TRUE, remove.duplicates = FALSE, size.selection = FALSE, strand.mode = 1) {
  dir.create(file.path(rle.out.folder, "fragment.size.rle.tracks"))
  anno[, c("start")] = anno[, c("start")] - 500
  anno[, c("end")] = anno[, c("end")] + 500
  anno[, "width"] = anno[, "end"] - anno[, "start"]
  
  for (bam.file in bam.files){
    dir.create(file.path(rle.out.folder, "fragment.size.rle.tracks", bam.file))
    
    registerDoParallel(cores = mc.cores)
    build.fragment.size.track.list = function(which.chr){
      fragment.size.track.list.chr = list()
      param = ScanBamParam(which=GRanges(seqnames = which.chr,ranges = IRanges(0,human.chrs.lengths[which.chr])))
      bam = readGAlignmentPairs(file = file.path(paste0(bam.input.folder, bam.file,".bam")),param = param)
      bam = bam[start(left(bam)) <= end(right(bam))]
      bam = bam[which(seqnames(bam) == which.chr)] # Remove non-concordant mappings
      if(size.selection) {
        bam = bam[(end(right(bam)) - start(left(bam))) <= 500] # Size selection
      }
      if(remove.duplicates) {
        bam = bam[!duplicated(paste(start(left(bam)),end(right(bam)))),] # Remove duplicates
      }
      
      bam.chr.ranges = GRanges(seqnames = which.chr,strand = strand(bam),ranges = IRanges(start = start(left(bam)), end = end(right(bam))))
      
      strand.chr.anno = makeGRangesFromDataFrame(anno[which(anno$chr == which.chr),c("chr","start","end","width", "strand")])
      bam.anno.overlaps = findOverlaps(bam.chr.ranges, strand.chr.anno, maxgap = 0L, minoverlap = 1L, type = "within", select = "all", ignore.strand = T)
      bam = bam[unique(queryHits(bam.anno.overlaps))]
      
      left.cigars = cigar(left(bam))
      left.cigar.ops = explodeCigarOps(left.cigars)
      right.cigars = cigar(right(bam))
      right.cigar.ops = explodeCigarOps(right.cigars)
      bam = bam[(!sapply(left.cigar.ops,function(x){any(x == "N")}) & !sapply(right.cigar.ops,function(x){any(x == "N")})),]
      
      strand.bam = bam[strand(bam) == "+"]
      fragment.size.vec = Rle(0,human.chrs.lengths[which.chr])
      if (length(strand.bam) > 0){
        left.cigars = cigar(left(strand.bam))
        left.cigar.ops = explodeCigarOps(left.cigars)
        left.cigar.lengths = sapply(explodeCigarOpLengths(left.cigars,ops=c("M","I","=","X")),sum) # sapply(explodeCigarOpLengths(left.cigars,ops=c("M","I","S","=","X")),sum)
        right.cigars = cigar(right(strand.bam))
        right.cigar.ops = explodeCigarOps(right.cigars)
        right.cigar.lengths = sapply(explodeCigarOpLengths(right.cigars,ops=c("M","I","=","X")),sum) # sapply(explodeCigarOpLengths(right.cigars,ops=c("M","I","S","=","X")),sum)
        starts = end(left(strand.bam))
        inner.mate.distance = start(right(strand.bam)) - starts
        mid.points = round(starts + (inner.mate.distance)/2)			
        fragment.size = inner.mate.distance + left.cigar.lengths + right.cigar.lengths
        if (length(fragment.size[fragment.size >= 0 & fragment.size <= 500]) > 0){
          mid.points = mid.points[fragment.size >= 0 & fragment.size <= 500]
          fragment.size = fragment.size[fragment.size >= 0 & fragment.size <= 500]
          try({fragment.size.aggregate = aggregate(fragment.size,list(id = mid.points),function(x){mean(x,na.rm = TRUE)});fragment.size.vec[fragment.size.aggregate[,"id"]] = fragment.size.aggregate[,"x"]},silent = TRUE)
        }
      }
      fragment.size.track.list.chr[["+"]] = fragment.size.vec
      
      strand.bam = bam[strand(bam) == "-"]
      fragment.size.vec = Rle(0,human.chrs.lengths[which.chr])
      if (length(strand.bam) > 0){
        left.cigars = cigar(left(strand.bam))
        left.cigar.ops = explodeCigarOps(left.cigars)
        left.cigar.lengths = sapply(explodeCigarOpLengths(left.cigars,ops=c("M","I","=","X")),sum) # sapply(explodeCigarOpLengths(left.cigars,ops=c("M","I","S","=","X")),sum)
        right.cigars = cigar(right(strand.bam))
        right.cigar.ops = explodeCigarOps(right.cigars)
        right.cigar.lengths = sapply(explodeCigarOpLengths(right.cigars,ops=c("M","I","=","X")),sum) # sapply(explodeCigarOpLengths(right.cigars,ops=c("M","I","S","=","X")),sum)
        starts = end(left(strand.bam))
        inner.mate.distance = start(right(strand.bam)) - starts
        mid.points = round(starts + (inner.mate.distance)/2)			
        fragment.size = inner.mate.distance + left.cigar.lengths + right.cigar.lengths
        if (length(fragment.size[fragment.size >= 0 & fragment.size <= 500]) > 0){
          mid.points = mid.points[fragment.size >= 0 & fragment.size <= 500]
          fragment.size = fragment.size[fragment.size >= 0 & fragment.size <= 500]
          try({fragment.size.aggregate = aggregate(fragment.size,list(id = mid.points),function(x){mean(x,na.rm = TRUE)});fragment.size.vec[fragment.size.aggregate[,"id"]] = fragment.size.aggregate[,"x"]},silent = TRUE)
        }
      }
      fragment.size.track.list.chr[["-"]] = fragment.size.vec
      if(!strand.specific) {
        fragment.size.track.list.chr[["+"]] = fragment.size.track.list.chr[["+"]] + fragment.size.track.list.chr[["-"]]
        fragment.size.track.list.chr[["-"]] = fragment.size.track.list.chr[["+"]]
      }
      
      save(fragment.size.track.list.chr,file = file.path(rle.out.folder, "fragment.size.rle.tracks", bam.file, paste0(prefix, which.chr, ".RData")))
      return()
    }
    
    fragment.size.track.list = foreach(n = human.chrs,.noexport = setdiff(ls(),c("human.chrs.lengths"))) %dopar% build.fragment.size.track.list(n)
  }
  
  gc()
}

create.transcribed.bases.rle.tracks.anno = function(bam.files, bam.input.folder, rle.out.folder, prefix, anno, human.chrs, human.chrs.lengths, strand.specific = TRUE, remove.duplicates = FALSE, size.selection = FALSE, strand.mode = 1) {
  dir.create(file.path(rle.out.folder, "transcribed.bases.rle.tracks"))
  anno[, c("start")] = anno[, c("start")] - 500
  anno[, c("end")] = anno[, c("end")] + 500
  anno[, "width"] = anno[, "end"] - anno[, "start"]
  for (bam.file in bam.files){
    print(bam.file)
    dir.create(file.path(file.path(rle.out.folder, "transcribed.bases.rle.tracks", bam.file)))
    
    registerDoParallel(cores = mc.cores)
    build.transcribed.bases.coverage.track.list = function(which.chr){
      print("Starting")
      print(which.chr)
      transcribed.bases.coverage.track.list.chr = list()
      param = ScanBamParam(which=GRanges(seqnames = which.chr,ranges = IRanges(0,human.chrs.lengths[which.chr])))
      bam = readGAlignmentPairs(file = file.path(paste0(bam.input.folder, bam.file,".bam")),param = param)
      bam = bam[start(left(bam)) <= end(right(bam))]
      bam = bam[which(seqnames(bam) == which.chr)] # Remove non-concordant mappings
      if(size.selection) {
        bam = bam[(end(right(bam)) - start(left(bam))) <= 500] # Size selection
      }
      if(remove.duplicates) {
        bam = bam[!duplicated(paste(start(left(bam)),end(right(bam)))),] # Remove duplicates
      }
      
      bam.chr.ranges = GRanges(seqnames = which.chr,strand = strand(bam),ranges = IRanges(start = start(left(bam)), end = end(right(bam))))
      
      strand.chr.anno = makeGRangesFromDataFrame(anno[which(anno$chr == which.chr),c("chr","start","end","width", "strand")])
      bam.anno.overlaps = findOverlaps(bam.chr.ranges, strand.chr.anno, maxgap = 0L, minoverlap = 1L, type = "within", select = "all", ignore.strand = T)
      bam = bam[unique(queryHits(bam.anno.overlaps))]
      
      rle.vec = Rle(0,human.chrs.lengths[which.chr])
      coverage.vec = coverage(GRanges(seqnames = which.chr,ranges = IRanges(start = start(left(bam[strand(bam) == "+"])),end = end(right(bam[strand(bam) == "+"])))))[[which.chr]]
      rle.vec[1:length(coverage.vec)] = coverage.vec
      transcribed.bases.coverage.track.list.chr[["+"]] = rle.vec
      
      rle.vec = Rle(0,human.chrs.lengths[which.chr])
      coverage.vec = coverage(GRanges(seqnames = which.chr,ranges = IRanges(start = start(left(bam[strand(bam) == "-"])),end = end(right(bam[strand(bam) == "-"])))))[[which.chr]]
      rle.vec[1:length(coverage.vec)] = coverage.vec
      transcribed.bases.coverage.track.list.chr[["-"]] = rle.vec
      if(!strand.specific) {
        transcribed.bases.coverage.track.list.chr[["+"]] = transcribed.bases.coverage.track.list.chr[["+"]] + transcribed.bases.coverage.track.list.chr[["-"]]
        transcribed.bases.coverage.track.list.chr[["-"]] = transcribed.bases.coverage.track.list.chr[["+"]]
      }
      
      save(transcribed.bases.coverage.track.list.chr,file = file.path(rle.out.folder, "transcribed.bases.rle.tracks", bam.file, paste0(prefix, which.chr, ".RData")))
      return()
    }
    
    transcribed.bases.coverage.track.list = foreach(n = human.chrs,.noexport = setdiff(ls(),c("human.chrs.lengths"))) %do% build.transcribed.bases.coverage.track.list(n)
  }
}

create.nexus.rle.tracks = function(bam.files, bam.input.folder, rle.out.folder, prefix, human.chrs, human.chrs.lengths, 
                                               strand.specific = TRUE, remove.duplicates = FALSE, size.selection = FALSE, 
                                               strand.mode = 1) {
  dir.create(file.path(rle.out.folder, "transcribed.bases.rle.tracks"))
  for (bam.file in bam.files){
    print(bam.file)
    dir.create(file.path(file.path(rle.out.folder, "transcribed.bases.rle.tracks", bam.file)))
    
    registerDoParallel(cores = mc.cores)
    bam.all = readRDS(paste0(bam.input.folder, bam.file,".bam"))
    build.transcribed.bases.coverage.track.list = function(which.chr){
      print(which.chr)
      transcribed.bases.coverage.track.list.chr = list()
      bam = bam.all[which(seqnames(bam.all) == which.chr)]
      
      bam.p <- bam[strand(bam) == "+"]
      bam.n <- bam[strand(bam) == "-"]
      
      bam.p <- resize(bam.p, 1)
      bam.n <- resize(bam.n, 1)
      
      
      if(length(bam) != 0){
        if(size.selection) {
          bam = bam[(end(right(bam)) - start(left(bam))) <= 500] # Size selection
        }
        if(remove.duplicates) {
          bam = bam[!duplicated(paste(start(left(bam)),end(right(bam)))),] # Remove duplicates
        }
        starts = start(ranges(bam.p))
        transcribed.bases.coverage.track.list.chr[["+"]] = 0
        if(length(starts) != 0){
          rle.vec = Rle(0,human.chrs.lengths[which.chr])
          coverage.vec = coverage(GRanges(seqnames = which.chr,ranges = IRanges(start = starts,end = starts)))[[which.chr]]
          rle.vec[1:length(coverage.vec)] = coverage.vec
          transcribed.bases.coverage.track.list.chr[["+"]] = rle.vec
        }
        
        starts = start(ranges(bam.n))
        transcribed.bases.coverage.track.list.chr[["-"]] = 0
        if(length(starts) != 0){
          rle.vec = Rle(0,human.chrs.lengths[which.chr])
          coverage.vec = coverage(GRanges(seqnames = which.chr,ranges = IRanges(start = starts,end = starts)))[[which.chr]]
          rle.vec[1:length(coverage.vec)] = coverage.vec
          transcribed.bases.coverage.track.list.chr[["-"]] = rle.vec
        }
        if(!strand.specific) {
          transcribed.bases.coverage.track.list.chr[["+"]] = transcribed.bases.coverage.track.list.chr[["+"]] + transcribed.bases.coverage.track.list.chr[["-"]]
          transcribed.bases.coverage.track.list.chr[["-"]] = transcribed.bases.coverage.track.list.chr[["+"]]
        }
        
        save(transcribed.bases.coverage.track.list.chr,file = file.path(rle.out.folder, "transcribed.bases.rle.tracks", bam.file, paste0(prefix, which.chr, ".RData")))
        print(paste0("Non zero RLE created for ", which.chr))
      }else{
        transcribed.bases.coverage.track.list.chr[["+"]] = 0
        if(!strand.specific){
          transcribed.bases.coverage.track.list.chr[["+"]] = 0
          transcribed.bases.coverage.track.list.chr[["-"]] = 0
        }
        save(transcribed.bases.coverage.track.list.chr,file = file.path(rle.out.folder, "transcribed.bases.rle.tracks", bam.file, paste0(prefix, which.chr, ".RData")))
        print(paste0("Zero RLE created for ", which.chr))
      }
      
    }
    
    transcribed.bases.coverage.track.list = foreach(n = human.chrs,.noexport = setdiff(ls(),c("human.chrs.lengths"))) %dopar% build.transcribed.bases.coverage.track.list(n)
  }
}

calculate.counts = function(anno, rle.location, bam.files, rle.prefix, out.folder, file.name, antisense.bias.ratio.location = "", chrs, chrs.lengths, antisense = F, find.antisense = T){
  if (TRUE){
    counts = list()
    for (bam.file in bam.files){
      print(bam.file)
      index.subsets = split(1:nrow(anno),paste(as.character(anno[,"strand"]),"/",anno[,"chr"],sep = ""))
      coverage.list = list()
      
      build.fragment.counts.list = function(j){
        from.transcript = strand.chr.anno[j,"start"]
        to.transcript = strand.chr.anno[j,"end"]
        sum.transcript = sum(as.vector(strand.chr.fragment.counts.from.bam[from.transcript:to.transcript]))
        names(sum.transcript) = j
        return(sum.transcript)
      }
      
      for (index.subset in names(index.subsets)){
        print(index.subset)
        fragment.mid.track.list.chr = get(load(file.path(paste0(rle.location,bam.file,"/",rle.prefix,unlist(strsplit(index.subset,split = "/"))[2],".RData"))))
        
        if(unlist(strsplit(index.subset,split = '/'))[1] == "*"){
          strand.chr.fragment.counts.from.bam = (fragment.mid.track.list.chr[["+"]] + fragment.mid.track.list.chr[["-"]])/2
        }
        else{
          strand.chr.fragment.counts.from.bam = fragment.mid.track.list.chr[[unlist(strsplit(index.subset,split = '/'))[1]]]
        }
        
        strand.chr.anno = anno[index.subsets[[index.subset]],c("start","end","width")]
        
        registerDoParallel(cores = mc.cores)
        coverage.list = c(coverage.list, foreach(n = rownames(strand.chr.anno),.noexport = setdiff(ls(),c("strand.chr.fragment.counts.from.bam","strand.chr.anno","build.fragment.counts.list"))) %dopar% build.fragment.counts.list(n)) 
      }
      counts[[bam.file]] = unlist(coverage.list,recursive = TRUE,use.names = TRUE)
    }	
    counts = sapply(counts,c)
    print("counts created")
    #rownames(counts) = as.character(anno[as.numeric(unlist(index.subsets,recursive = TRUE,use.names = TRUE)),"trid"])
    colnames(counts) = bam.files
    if(antisense) {
      # Antisense bias
      if(find.antisense){
        antisense.bias.ratio.mat = cbind()
      for (bam.file in bam.files){
        registerDoParallel(cores = mc.cores)
        build.antisense.bias.ratios = function(which.chr,distance = 0){
          chr.coverage.list = get(load(file.path(paste0(rle.location,bam.file,"/",rle.prefix,which.chr,".RData"))))
          anno.positions = Rle(0,chrs.lengths[which.chr])
          check.ids = rownames(anno[which(anno[,"strand"] == "+" & anno[,"chr"] == which.chr),])
          for (check.id in check.ids){anno.positions[(anno[check.id,"start"] - distance):(anno[check.id,"end"] + distance)] = 1}
          diff.positions = Rle(1,chrs.lengths[which.chr])
          check.ids = rownames(anno[which(anno[,"strand"] == "-" & anno[,"chr"] == which.chr),])
          for (check.id in check.ids){diff.positions[(anno[check.id,"start"] - distance):(anno[check.id,"end"] + distance)] = 0}
          positions = anno.positions*diff.positions
          runValue(positions) = as.logical(runValue(positions))
          sense = chr.coverage.list[["+"]][positions]
          antisense = chr.coverage.list[["-"]][positions]
          anno.positions = Rle(0,chrs.lengths[which.chr])
          check.ids = rownames(anno[which(anno[,"strand"] == "-" & anno[,"chr"] == which.chr),])
          for (check.id in check.ids){anno.positions[(anno[check.id,"start"] - distance):(anno[check.id,"end"] + distance)] = 1}
          diff.positions = Rle(1,chrs.lengths[which.chr])
          check.ids = rownames(anno[which(anno[,"strand"] == "+" & anno[,"chr"] == which.chr),])
          for (check.id in check.ids){diff.positions[(anno[check.id,"start"] - distance):(anno[check.id,"end"] + distance)] = 0}
          positions = anno.positions*diff.positions
          runValue(positions) = as.logical(runValue(positions))
          sense = c(sense,chr.coverage.list[["-"]][positions])
          antisense = c(antisense,chr.coverage.list[["+"]][positions])
          sense[which(sense < 100)] = NA
          antisense[which(antisense == 0)] = NA
          return(median(antisense/sense,na.rm = TRUE))
        }
        antisense.bias.ratios = as.vector(foreach(n = chrs) %dopar% build.antisense.bias.ratios(n))
        
        antisense.bias.ratio.mat = cbind(antisense.bias.ratio.mat,antisense.bias.ratios)
      }
      
      antisense.bias.ratio.mat = apply(antisense.bias.ratio.mat,c(1,2),as.numeric)
      rownames(antisense.bias.ratio.mat) = chrs
      colnames(antisense.bias.ratio.mat) = bam.files
      
      antisense.bias.ratio = apply(antisense.bias.ratio.mat,2,function(x){median(x,na.rm = TRUE)})
      antisense.bias.ratio[is.na(antisense.bias.ratio)] = 0
      save(antisense.bias.ratio, file = antisense.bias.ratio.location)
      }
      # Antisense counts
      antisense.counts = list()
      for (bam.file in bam.files) {
        index.subsets = split(1:nrow(anno),paste(Vectorize(strand.switch)(as.character(anno[,"strand"])),"/",anno[,"chr"],sep = ""))
        coverage.list = list()
        
        build.fragment.counts.list = function(j){
          from.transcript = strand.chr.anno[j,"start"]
          to.transcript = strand.chr.anno[j,"end"]
          sum.transcript = sum(as.vector(strand.chr.fragment.counts.from.bam[from.transcript:to.transcript]))
          names(sum.transcript) = j
          return(sum.transcript)
        }
        
        for (index.subset in names(index.subsets)){
          fragment.mid.track.list.chr = get(load(file.path(paste0(rle.location,bam.file,"/",rle.prefix,unlist(strsplit(index.subset,split = "/"))[2],".RData"))))
          
          strand.chr.fragment.counts.from.bam = fragment.mid.track.list.chr[[unlist(strsplit(index.subset,split = "/"))[1]]]
          strand.chr.anno = anno[index.subsets[[index.subset]],c("start","end","width")]
          
          registerDoParallel(cores = mc.cores)
          coverage.list = c(coverage.list, foreach(n = rownames(strand.chr.anno),.noexport = setdiff(ls(),c("strand.chr.fragment.counts.from.bam","strand.chr.anno","build.fragment.counts.list"))) %dopar% build.fragment.counts.list(n)) 
        }
        antisense.counts[[bam.file]] = unlist(coverage.list,recursive = TRUE,use.names = TRUE)
      }	
      antisense.counts = sapply(antisense.counts,c)
      print("antisense counts created")
      #rownames(antisense.counts) = as.character(anno[as.numeric(unlist(index.subsets,recursive = TRUE,use.names = TRUE)),"trid"])
      colnames(antisense.counts) = bam.files
      
      antisense.bias.ratio = get(load(antisense.bias.ratio.location))
      antisense.bias.ratio = antisense.bias.ratio[bam.files]
      
      if(all(colnames(counts) == names(antisense.bias.ratio))) {
        print("saving counts")
        index = rownames(anno)
        counts.antisense.corrected = t(t(counts[index, ] - t(t(antisense.counts[index, ])*antisense.bias.ratio))/(1 - antisense.bias.ratio^2))
        counts.antisense.corrected[counts.antisense.corrected < 0] = 0
        counts = cbind(counts[index, ], antisense.counts[index, ])
        save(counts, file = paste0(out.folder, file.name,".RData"))
        save(counts.antisense.corrected, file = paste0(out.folder,file.name, ".antisense.corrected.RData"))
      }
      else {
        print("Incorrect antisense bias ratio file!!")
      }
    }
    else {
      print("saving counts")
      index = rownames(anno)
      counts = cbind(counts[index, ])
      save(counts, file = paste0(out.folder, file.name,".RData"))
    }
  }
}

create.normalized.rle.tracks = function(bam.files, seq.depth, rle.folder, out.folder, prefix, human.chrs){
  dir.create(file.path(out.folder))
  for (bam.file in bam.files){
    dir.create(file.path(paste0(out.folder,bam.file)))
    for (which.chr in human.chrs){
      transcribed.bases.coverage.track.list.chr = get(load(file.path(paste0(rle.folder,bam.file,"/",prefix,which.chr,".RData"))))
      transcribed.bases.coverage.track.list.chr[["+"]] = transcribed.bases.coverage.track.list.chr[["+"]]/seq.depth[bam.file]
      transcribed.bases.coverage.track.list.chr[["-"]] = transcribed.bases.coverage.track.list.chr[["-"]]/seq.depth[bam.file]
      save(transcribed.bases.coverage.track.list.chr,file = file.path(paste0(out.folder,bam.file,"/",prefix,which.chr,".RData")))
    }
  }
}

create.bigwig.from.rle = function(rle.folder, rle.prefix, out.folder, bam.files, strand.specific = TRUE, human.chrs){
  dir.create(file.path(out.folder))
  for (i in 1:length(bam.files)) {
      # +
      rle.list.chr1 =  get(load(file.path(paste0(rle.folder,bam.files[i],"/",rle.prefix,human.chrs[1],".RData"))))
      g.ranges = RleList(rle.list.chr1[["+"]])
      names(g.ranges) = "chr1"
      g.ranges = as(RleList(g.ranges),"GRanges")
      strand(g.ranges) = "+"
      score(g.ranges)[is.na(score(g.ranges))] = 0
      bw.g.ranges = g.ranges
      for (human.chr in setdiff(human.chrs,"chr1")){
        rle.list.chr1 =  get(load(file.path(paste0(rle.folder,bam.files[i],"/",rle.prefix,human.chr,".RData"))))
        g.ranges = RleList(rle.list.chr1[["+"]])
        names(g.ranges) = human.chr
        g.ranges = as(RleList(g.ranges),"GRanges")
        strand(g.ranges) = "+"
        score(g.ranges)[is.na(score(g.ranges))] = 0
        bw.g.ranges = c(bw.g.ranges,g.ranges)
      }
      export(bw.g.ranges, con = file.path(paste0(out.folder,bam.files[i],"plus.bw")))
      
      # -
      if(strand.specific){
        rle.list.chr1 =  get(load(file.path(paste0(rle.folder,bam.files[i],"/",rle.prefix,human.chrs[1],".RData"))))
        g.ranges = RleList(rle.list.chr1[["-"]])
        names(g.ranges) = "chr1"
        g.ranges = as(RleList(g.ranges),"GRanges")
        strand(g.ranges) = "-"
        score(g.ranges)[is.na(score(g.ranges))] = 0
        bw.g.ranges = g.ranges
        for (human.chr in setdiff(human.chrs,"chr1")){
          rle.list.chr1 =  get(load(file.path(paste0(rle.folder,bam.files[i],"/",rle.prefix,human.chr,".RData"))))
          g.ranges = RleList(rle.list.chr1[["-"]])
          names(g.ranges) = human.chr
          g.ranges = as(RleList(g.ranges),"GRanges")
          strand(g.ranges) = "-"
          score(g.ranges)[is.na(score(g.ranges))] = 0
          bw.g.ranges = c(bw.g.ranges,g.ranges)
        }
        export(bw.g.ranges, con = file.path(paste0(out.folder,bam.files[i],"minus.bw")))
      }
    }
  }

z.transform = function(mat){mat = (mat - apply(mat,1,mean))/apply(mat,1,sd);return(mat)}

make.average.dataset = function(df, col.names){
  index = rownames(df)
  n.row = nrow(df)
  n.col = ncol(df)
  avg.df = as.data.frame(matrix(nrow = n.row, ncol = n.col/2))
  rownames(avg.df) = index
  colnames(avg.df) = col.names
  df[is.na(df)] = NA
  #df[is.infinite(df)] = NA
  #df[is.nan(df)] = NA

  for(i in 1:(n.col/2)){
    avg.df[index, i] = rowMeans(df[index, ((2*i)-1):(2*i)], na.rm = T)
  }
  return(avg.df)
}

create.coverages = function(bam.files, replicate.per.sample, gene.anno, from, to, window, 
                            rle.prefix, rle.folder, out.folder, out.file, on, align){
  # bam files list and rle locations
  bam.files.list = list()
  for (i in 1:(length(bam.files)/replicate.per.sample)) {
    bam.files.list[[i]] = bam.files[1:replicate.per.sample]
    bam.files = bam.files[(replicate.per.sample+1):length(bam.files)]
  }
  bam.files.list.rle =lapply(bam.files.list, FUN = function(x) {paste0(rle.folder, x, "/")})
  
  # outfile
  outfile = paste0(out.folder, out.file, "/cm.", from, "_", to, "_", window, "_", on, "_", align, ".RData")
  
  # calculate and save coverage matrices
  coverage.matrix.list = list()
  for(i in 1:length(bam.files.list)){
    coverage.matrix.list[[i]] = build.coverage.profiles(gene.anno = gene.anno, bam.file.subsets = bam.files.list.rle[i], on = on, align = align,
                                                        from = from, to = to, window = window, filter.filter = NULL, mc.cores = detectCores())
  }
  
  # save output
  save(coverage.matrix.list, file = outfile)
}

create.scaled.coverages = function(bam.files, replicate.per.sample, gene.anno, from, to, median.length, window, 
                                   rle.prefix, rle.folder, out.folder, out.file, on){
  # bam files list and rle locations
  bam.files.list = list()
  for (i in 1:(length(bam.files)/replicate.per.sample)) {
    bam.files.list[[i]] = bam.files[1:replicate.per.sample]
    bam.files = bam.files[(replicate.per.sample+1):length(bam.files)]
  }
  bam.files.list.rle =lapply(bam.files.list, FUN = function(x) {paste0(rle.folder, x, "/")})
  
  # outfile
  outfile = paste0(out.folder, out.file, "/scm.", from, "_", to, "_", window, "_", on, "_", align, ".RData")
  
  # calculate and save coverage matrices
  coverage.matrix.list = list()
  for(i in 1:length(bam.files.list)){
    coverage.matrix.list[[i]] = build.scaled.coverage.profiles(gene.anno = gene.anno, bam.file.subsets = bam.files.list.rle, on = on,
                                                               from = from, to = to, window = window, median.length = median.length, 
                                                               filter.filter = NULL, mc.cores = detectCores())
  }
  
  # save output
  save(coverage.matrix.list, file = outfile)
}
