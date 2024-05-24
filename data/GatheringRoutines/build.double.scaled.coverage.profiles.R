

### build.double.scaled.coverage.profiles ###


build.double.scaled.coverage.profiles = function(
		gene.anno,
		bam.file.subsets,
		on,
		from,
		to,
		window,
		median.first.region.length,
		median.second.region.length,
		first.region.length.column.id,
		filter.filter,
		mc.cores
)
{	
  if (!is.null(filter.filter)){window = length(filter.filter);filter.filter = filter.filter/sum(filter.filter)}
  
  adjusted.from = from - ceiling(window/2)
  adjusted.to = to + ceiling(window/2)
  
  ### gather profiles ###
  
  coverage.matrix.list = list()
  strand.switch = function(which.strand){switch(which.strand,"+" = "-","-" = "+")}
  for (bam.file.subset in 1:length(bam.file.subsets)){
    coverage.list = list()
    for (cov.path in bam.file.subsets[[bam.file.subset]]){
      index.subsets = split(1:nrow(gene.anno),as.character(gene.anno[,"chr"]))
      
      build.coverage.list = function(j){
        from.transcript = chr.gene.anno[j,"start"]
        to.transcript = chr.gene.anno[j,"end"]
        strand.transcript = as.character(chr.gene.anno[j,"strand"])
        
        if(on == "antisense"){
          if(strand.transcript == "+"){
            strand.length = length(chr.coverage.list.from.bam[[strand.switch(strand.transcript)]])
            upstream.positions = (from.transcript+adjusted.from):(from.transcript-1)
            valid.upstream.positions = upstream.positions[upstream.positions > 0 & upstream.positions <= strand.length]
            positions = from.transcript:(from.transcript+anno.lengths.first[j])
            valid.positions = positions[positions > 0 & positions <= strand.length]
            positions.second = (from.transcript+anno.lengths.first[j]+1):(from.transcript+anno.lengths.first[j]+anno.lengths.second[j])
            valid.positions.second = positions.second[positions.second > 0 & positions.second <= strand.length]
            downstream.positions = (from.transcript+anno.lengths.first[j]+anno.lengths.second[j]+1):(from.transcript+anno.lengths.first[j]+anno.lengths.second[j]-median.first.region.length-median.second.region.length+adjusted.to)
            valid.downstream.positions = downstream.positions[downstream.positions > 0 & downstream.positions <= strand.length]
            vec.template = Rle(rep(0,length(adjusted.from:adjusted.to)))
            if (length(valid.upstream.positions) > 0){
              cov.vec = chr.coverage.list.from.bam[[strand.switch(strand.transcript)]][valid.upstream.positions]
              vec.template[c(upstream.positions > 0 & upstream.positions <= strand.length,rep(FALSE,(median.first.region.length+1)),rep(FALSE,(median.second.region.length)),rep(FALSE,length(downstream.positions)))] = cov.vec
            }
            if (length(valid.positions) > 0){
              cov.vec = Rle(resize.profile(as.vector(chr.coverage.list.from.bam[[strand.switch(strand.transcript)]][valid.positions]),(median.first.region.length+1)))
              vec.template[c(rep(FALSE,length(upstream.positions)),rep(TRUE,(median.first.region.length+1)),rep(FALSE,(median.second.region.length)),rep(FALSE,length(downstream.positions)))] = cov.vec
            }
            if (length(valid.positions.second) > 0){
              cov.vec = Rle(resize.profile(as.vector(chr.coverage.list.from.bam[[strand.switch(strand.transcript)]][valid.positions.second]),(median.second.region.length)))
              vec.template[c(rep(FALSE,length(upstream.positions)),rep(FALSE,(median.first.region.length+1)),rep(TRUE,(median.second.region.length)),rep(FALSE,length(downstream.positions)))] = cov.vec
            }
            if (length(valid.downstream.positions) > 0){
              cov.vec = chr.coverage.list.from.bam[[strand.switch(strand.transcript)]][valid.downstream.positions]
              vec.template[c(rep(FALSE,length(upstream.positions)),rep(FALSE,(median.first.region.length+1)),rep(FALSE,(median.second.region.length)),downstream.positions > 0 & downstream.positions <= strand.length)] = cov.vec
            }
            return(vec.template)
          } else {
            strand.length = length(chr.coverage.list.from.bam[[strand.switch(strand.transcript)]])
            upstream.positions = (to.transcript-adjusted.from):(to.transcript+1)
            valid.upstream.positions = upstream.positions[upstream.positions > 0 & upstream.positions <= strand.length]
            positions = to.transcript:(to.transcript-anno.lengths.first[j])
            valid.positions = positions[positions > 0 & positions <= strand.length]
            positions.second = (to.transcript-anno.lengths.first[j]-1):(to.transcript-anno.lengths.first[j]-anno.lengths.second[j])
            valid.positions.second = positions.second[positions.second > 0 & positions.second <= strand.length]
            downstream.positions = (to.transcript-anno.lengths.first[j]-anno.lengths.second[j]-1):(to.transcript-anno.lengths.first[j]-anno.lengths.second[j]+median.first.region.length+median.second.region.length-adjusted.to)
            valid.downstream.positions = downstream.positions[downstream.positions > 0 & downstream.positions <= strand.length]
            vec.template = Rle(rep(0,length(adjusted.from:adjusted.to)))
            if (length(valid.upstream.positions) > 0){
              cov.vec = chr.coverage.list.from.bam[[strand.switch(strand.transcript)]][valid.upstream.positions]
              vec.template[c(upstream.positions > 0 & upstream.positions <= strand.length,rep(FALSE,(median.first.region.length+1)),rep(FALSE,(median.second.region.length)),rep(FALSE,length(downstream.positions)))] = cov.vec
            }
            if (length(valid.positions) > 0){
              cov.vec = Rle(resize.profile(as.vector(chr.coverage.list.from.bam[[strand.switch(strand.transcript)]][valid.positions]),(median.first.region.length+1)))
              vec.template[c(rep(FALSE,length(upstream.positions)),rep(TRUE,(median.first.region.length+1)),rep(FALSE,(median.second.region.length)),rep(FALSE,length(downstream.positions)))] = cov.vec
            }
            if (length(valid.positions.second) > 0){
              cov.vec = Rle(resize.profile(as.vector(chr.coverage.list.from.bam[[strand.switch(strand.transcript)]][valid.positions.second]),(median.second.region.length)))
              vec.template[c(rep(FALSE,length(upstream.positions)),rep(FALSE,(median.first.region.length+1)),rep(TRUE,(median.second.region.length)),rep(FALSE,length(downstream.positions)))] = cov.vec
            }
            if (length(valid.downstream.positions) > 0){
              cov.vec = chr.coverage.list.from.bam[[strand.switch(strand.transcript)]][valid.downstream.positions]
              vec.template[c(rep(FALSE,length(upstream.positions)),rep(FALSE,(median.first.region.length+1)),rep(FALSE,(median.second.region.length)),downstream.positions > 0 & downstream.positions <= strand.length)] = cov.vec
            }
            return(vec.template)
          }
        } else if(on == "sense"){
          if(strand.transcript == "+"){
            strand.length = length(chr.coverage.list.from.bam[[strand.transcript]])
            upstream.positions = (from.transcript+adjusted.from):(from.transcript-1)
            valid.upstream.positions = upstream.positions[upstream.positions > 0 & upstream.positions <= strand.length]
            positions = from.transcript:(from.transcript+anno.lengths.first[j])
            valid.positions = positions[positions > 0 & positions <= strand.length]
            positions.second = (from.transcript+anno.lengths.first[j]+1):(from.transcript+anno.lengths.first[j]+anno.lengths.second[j])
            valid.positions.second = positions.second[positions.second > 0 & positions.second <= strand.length]
            downstream.positions = (from.transcript+anno.lengths.first[j]+anno.lengths.second[j]+1):(from.transcript+anno.lengths.first[j]+anno.lengths.second[j]-median.first.region.length-median.second.region.length+adjusted.to)
            valid.downstream.positions = downstream.positions[downstream.positions > 0 & downstream.positions <= strand.length]
            vec.template = Rle(rep(0,length(adjusted.from:adjusted.to)))
            if (length(valid.upstream.positions) > 0){
              cov.vec = chr.coverage.list.from.bam[[strand.transcript]][valid.upstream.positions]
              vec.template[c(upstream.positions > 0 & upstream.positions <= strand.length,rep(FALSE,(median.first.region.length+1)),rep(FALSE,(median.second.region.length)),rep(FALSE,length(downstream.positions)))] = cov.vec
            }
            if (length(valid.positions) > 0){
              cov.vec = Rle(resize.profile(as.vector(chr.coverage.list.from.bam[[strand.transcript]][valid.positions]),(median.first.region.length+1)))
              vec.template[c(rep(FALSE,length(upstream.positions)),rep(TRUE,(median.first.region.length+1)),rep(FALSE,(median.second.region.length)),rep(FALSE,length(downstream.positions)))] = cov.vec
            }
            if (length(valid.positions.second) > 0){
              cov.vec = Rle(resize.profile(as.vector(chr.coverage.list.from.bam[[strand.transcript]][valid.positions.second]),(median.second.region.length)))
              vec.template[c(rep(FALSE,length(upstream.positions)),rep(FALSE,(median.first.region.length+1)),rep(TRUE,(median.second.region.length)),rep(FALSE,length(downstream.positions)))] = cov.vec
            }
            if (length(valid.downstream.positions) > 0){
              cov.vec = chr.coverage.list.from.bam[[strand.transcript]][valid.downstream.positions]
              vec.template[c(rep(FALSE,length(upstream.positions)),rep(FALSE,(median.first.region.length+1)),rep(FALSE,(median.second.region.length)),downstream.positions > 0 & downstream.positions <= strand.length)] = cov.vec
            }
            return(vec.template)
          } else {
            strand.length = length(chr.coverage.list.from.bam[[strand.transcript]])
            upstream.positions = (to.transcript-adjusted.from):(to.transcript+1)
            valid.upstream.positions = upstream.positions[upstream.positions > 0 & upstream.positions <= strand.length]
            positions = to.transcript:(to.transcript-anno.lengths.first[j])
            valid.positions = positions[positions > 0 & positions <= strand.length]
            positions.second = (to.transcript-anno.lengths.first[j]-1):(to.transcript-anno.lengths.first[j]-anno.lengths.second[j])
            valid.positions.second = positions.second[positions.second > 0 & positions.second <= strand.length]
            downstream.positions = (to.transcript-anno.lengths.first[j]-anno.lengths.second[j]-1):(to.transcript-anno.lengths.first[j]-anno.lengths.second[j]+median.first.region.length+median.second.region.length-adjusted.to)
            valid.downstream.positions = downstream.positions[downstream.positions > 0 & downstream.positions <= strand.length]
            vec.template = Rle(rep(0,length(adjusted.from:adjusted.to)))
            if (length(valid.upstream.positions) > 0){
              cov.vec = chr.coverage.list.from.bam[[strand.transcript]][valid.upstream.positions]
              vec.template[c(upstream.positions > 0 & upstream.positions <= strand.length,rep(FALSE,(median.first.region.length+1)),rep(FALSE,(median.second.region.length)),rep(FALSE,length(downstream.positions)))] = cov.vec
            }
            if (length(valid.positions) > 0){
              cov.vec = Rle(resize.profile(as.vector(chr.coverage.list.from.bam[[strand.transcript]][valid.positions]),(median.first.region.length+1)))
              vec.template[c(rep(FALSE,length(upstream.positions)),rep(TRUE,(median.first.region.length+1)),rep(FALSE,(median.second.region.length)),rep(FALSE,length(downstream.positions)))] = cov.vec
            }
            if (length(valid.positions.second) > 0){
              cov.vec = Rle(resize.profile(as.vector(chr.coverage.list.from.bam[[strand.transcript]][valid.positions.second]),(median.second.region.length)))
              vec.template[c(rep(FALSE,length(upstream.positions)),rep(FALSE,(median.first.region.length+1)),rep(TRUE,(median.second.region.length)),rep(FALSE,length(downstream.positions)))] = cov.vec
            }
            if (length(valid.downstream.positions) > 0){
              cov.vec = chr.coverage.list.from.bam[[strand.transcript]][valid.downstream.positions]
              vec.template[c(rep(FALSE,length(upstream.positions)),rep(FALSE,(median.first.region.length+1)),rep(FALSE,(median.second.region.length)),downstream.positions > 0 & downstream.positions <= strand.length)] = cov.vec
            }
            return(vec.template)
          }
        }
      }
      
      sub.coverage.list = list()
      
      for (index.subset in names(index.subsets)){
        chr.coverage.list.from.bam = get(load(paste0(cov.path,index.subset,".RData")))
        chr.gene.anno = gene.anno[which(gene.anno[,"chr"] == index.subset),c("start","end","strand")]
        anno.lengths.first = gene.anno[which(gene.anno[,"chr"] == index.subset),first.region.length.column.id]
        anno.lengths.second = gene.anno[which(gene.anno[,"chr"] == index.subset),"length"] - gene.anno[which(gene.anno[,"chr"] == index.subset),first.region.length.column.id]
        
        registerDoParallel(cores = mc.cores)
        sub.coverage.list[[index.subset]] = foreach(n = 1:nrow(chr.gene.anno),.noexport = setdiff(ls(),c("chr.coverage.list.from.bam","chr.gene.anno","anno.lengths.first","anno.lengths.second","build.coverage.list"))) %dopar% build.coverage.list(n)
      }
      coverage.list[[cov.path]] = RleList(Reduce(c,sub.coverage.list))
    }
    coverage.list = Reduce('+',coverage.list)
    if (is.null(filter.filter)){filter.filter = rep(1,window)/window}
    coverage.list = mclapply(coverage.list,function(x){filter(as.vector(x),filter.filter)},mc.cores = mc.cores)
    coverage.mat = t(sapply(coverage.list,c))
    rownames(coverage.mat) = as.character(gene.anno[as.numeric(unlist(index.subsets)),"id"])
    colnames(coverage.mat) = adjusted.from:adjusted.to
    coverage.mat = coverage.mat[,as.character(from:to)]
    coverage.matrix = coverage.mat[as.character(gene.anno[,"id"]),]		
    coverage.matrix.list[[bam.file.subset]] = coverage.matrix
  }
  return(coverage.matrix.list)
}



