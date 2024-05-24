

### build.coverage.profiles ###


build.coverage.profiles = function(
		gene.anno,
		bam.file.subsets,
		align,
		on,
		from,
		to,
		window,
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
            if (align == "end"){
              positions = (to.transcript+adjusted.from):(to.transcript+adjusted.to)
              valid.positions = positions[positions > 0 & positions <= strand.length]
              vec.template = Rle(rep(NA,length(adjusted.from:adjusted.to)))
              if (length(valid.positions) > 0){
                cov.vec = chr.coverage.list.from.bam[[strand.switch(strand.transcript)]][valid.positions]
                vec.template[positions > 0 & positions <= strand.length] = cov.vec
              }
              return(vec.template)
            } else if (align == "start"){
              positions = (from.transcript+adjusted.from):(from.transcript+adjusted.to)
              valid.positions = positions[positions > 0 & positions <= strand.length]
              vec.template = Rle(rep(NA,length(adjusted.from:adjusted.to)))
              if (length(valid.positions) > 0){
                cov.vec = chr.coverage.list.from.bam[[strand.switch(strand.transcript)]][valid.positions]
                vec.template[positions > 0 & positions <= strand.length] = cov.vec
              }
              return(vec.template)
            }
          } else {
            strand.length = length(chr.coverage.list.from.bam[[strand.switch(strand.transcript)]])
            if (align == "end"){
              positions = (from.transcript-adjusted.from):(from.transcript-adjusted.to)
              valid.positions = positions[positions > 0 & positions <= strand.length]
              vec.template = Rle(rep(NA,length(adjusted.from:adjusted.to)))
              if (length(valid.positions) > 0){
                cov.vec = chr.coverage.list.from.bam[[strand.switch(strand.transcript)]][valid.positions]
                vec.template[positions > 0 & positions <= strand.length] = cov.vec
              }
              return(vec.template)
            } else if (align == "start"){
              positions = (to.transcript-adjusted.from):(to.transcript-adjusted.to)
              valid.positions = positions[positions > 0 & positions <= strand.length]
              vec.template = Rle(rep(NA,length(adjusted.from:adjusted.to)))
              if (length(valid.positions) > 0){
                cov.vec = chr.coverage.list.from.bam[[strand.switch(strand.transcript)]][valid.positions]
                vec.template[positions > 0 & positions <= strand.length] = cov.vec
              }
              return(vec.template)
            }
          }
        } else if(on == "sense"){
          if(strand.transcript == "+"){
            strand.length = length(chr.coverage.list.from.bam[[strand.transcript]])
            if (align == "end"){
              positions = (to.transcript+adjusted.from):(to.transcript+adjusted.to)
              valid.positions = positions[positions > 0 & positions <= strand.length]
              vec.template = Rle(rep(NA,length(adjusted.from:adjusted.to)))
              if (length(valid.positions) > 0){
                cov.vec = chr.coverage.list.from.bam[[strand.transcript]][valid.positions]
                vec.template[positions > 0 & positions <= strand.length] = cov.vec
              }
              return(vec.template)
            } else if (align == "start"){
              positions = (from.transcript+adjusted.from):(from.transcript+adjusted.to)
              valid.positions = positions[positions > 0 & positions <= strand.length]
              vec.template = Rle(rep(NA,length(adjusted.from:adjusted.to)))
              if (length(valid.positions) > 0){
                cov.vec = chr.coverage.list.from.bam[[strand.transcript]][valid.positions]
                vec.template[positions > 0 & positions <= strand.length] = cov.vec
              }
              return(vec.template)
            }
          } else {
            strand.length = length(chr.coverage.list.from.bam[[strand.transcript]])
            if (align == "end"){
              positions = (from.transcript-adjusted.from):(from.transcript-adjusted.to)
              valid.positions = positions[positions > 0 & positions <= strand.length]
              vec.template = Rle(rep(NA,length(adjusted.from:adjusted.to)))
              if (length(valid.positions) > 0){
                cov.vec = chr.coverage.list.from.bam[[strand.transcript]][valid.positions]
                vec.template[positions > 0 & positions <= strand.length] = cov.vec
              }
              return(vec.template)
            } else if (align == "start"){
              positions = (to.transcript-adjusted.from):(to.transcript-adjusted.to)
              valid.positions = positions[positions > 0 & positions <= strand.length]
              vec.template = Rle(rep(NA,length(adjusted.from:adjusted.to)))
              if (length(valid.positions) > 0){
                cov.vec = chr.coverage.list.from.bam[[strand.transcript]][valid.positions]
                vec.template[positions > 0 & positions <= strand.length] = cov.vec
              }
              return(vec.template)
            }
          }
        }
      }
      
      sub.coverage.list = list()
      
      for (index.subset in names(index.subsets)){
        chr.coverage.list.from.bam = get(load(paste0(cov.path,index.subset,".RData")))
        chr.gene.anno = gene.anno[which(gene.anno[,"chr"] == index.subset),c("start","end","strand")]
        
        registerDoParallel(cores = mc.cores)
        sub.coverage.list[[index.subset]] = foreach(n = 1:nrow(chr.gene.anno),.noexport = setdiff(ls(),c("chr.coverage.list.from.bam","chr.gene.anno","build.coverage.list"))) %dopar% build.coverage.list(n)
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



