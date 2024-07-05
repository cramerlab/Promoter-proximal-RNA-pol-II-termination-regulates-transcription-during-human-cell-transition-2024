

### position.weight.matrix ###


position.weight.matrix = function(
		gene.anno,
		non = 1,
		align = "TTS",
		from = -10,
		to = 10,
		pseudo = TRUE,
		frequencies = TRUE,
		mc.cores = 4
)
{	
	### gather multinucleotide frequency data ###
	
	seq.length = length((from - non):(to + non))
	
	index.subsets = split(1:nrow(gene.anno),as.character(gene.anno[,"chr"]))
	
	nucleotide.list = list()
	
	build.nucleotide.list = function(j){
		from.transcript = chr.gene.anno[j,"start"]
		to.transcript = chr.gene.anno[j,"end"]
		strand.transcript = as.character(chr.gene.anno[j,"strand"])
		
		if(strand.transcript == "+"){
			if (align == "TTS"){
				nucleotides = unlist(strsplit(as.character(subseq(chr.ref.genome,start=(to.transcript+(from - non)),width=seq.length)[[1]]),split=""))
			} else if (align == "TSS"){
				nucleotides = unlist(strsplit(as.character(subseq(chr.ref.genome,start=(from.transcript+(from - non)),width=seq.length)[[1]]),split=""))
			}
		} else {
			if (align == "TTS"){
				nucleotides = unlist(strsplit(as.character(reverseComplement(subseq(chr.ref.genome,start=(from.transcript-(to + non)),width=seq.length))[[1]]),split=""))
			} else if (align == "TSS"){
				nucleotides = unlist(strsplit(as.character(reverseComplement(subseq(chr.ref.genome,start=(to.transcript-(to + non)),width=seq.length))[[1]]),split=""))
			}
		}
		return(as.numeric(switchvector(multi.mer.list[[1]])[nucleotides]))
	}
	
	for (index.subset in names(index.subsets)){
	  chr.ref.genome = readDNAStringSet(filepath = file.path(prewd,"RawData","HumanGenomeReferenceConsortium","fasta_hg20_GRCh38_April_2014",paste0("GRCh38.",index.subset,".fa")))
		chr.gene.anno = gene.anno[which(gene.anno[,"chr"] == index.subset),c("start","end","strand")]
		
		registerDoParallel(cores = mc.cores)
		nucleotide.list[[index.subset]] = foreach(n = 1:nrow(chr.gene.anno),.noexport = setdiff(ls(),c("chr.ref.genome","chr.gene.anno","multi.mer.list"))) %dopar% build.nucleotide.list(n)
	}
	nucleotide.list = Reduce(c,nucleotide.list)
	
	if (non > 1){nucleotide.list = mclapply(nucleotide.list,function(x){as.numeric(filter(x,5^(0:8)[1:non],sides = 1)-(sum(5^(0:8)[1:non])-1))},mc.cores = mc.cores)}
	nucleotide.mat = t(sapply(nucleotide.list,c))
	colnames(nucleotide.mat) = ((from - non):(to + non))-(non-1)
	nucleotide.mat = nucleotide.mat[,as.character(from:to)]
	dim.nucleotide.mat = dim(nucleotide.mat)
	nucleotide.mat = multi.mer.list[[non]][nucleotide.mat]
	dim(nucleotide.mat) = dim.nucleotide.mat
	
	rownames(nucleotide.mat) = as.character(gene.anno[as.numeric(unlist(index.subsets)),"id"])
	colnames(nucleotide.mat) = from:to
	if (pseudo){
	  if (frequencies){
	    pwm = apply(nucleotide.mat,2,function(x){freqs = summary(factor(x,levels = multi.mer.list[[as.character(non)]]),maxsum = length(multi.mer.list[[as.character(non)]])+2)[multi.mer.list[[as.character(non)]]] + oligonucleotide.list[[as.character(non)]][multi.mer.list[[as.character(non)]]];freqs = freqs/sum(freqs)})
	  } else {
	    pwm = apply(nucleotide.mat,2,function(x){freqs = summary(factor(x,levels = multi.mer.list[[as.character(non)]]),maxsum = length(multi.mer.list[[as.character(non)]])+2)[multi.mer.list[[as.character(non)]]] + oligonucleotide.list[[as.character(non)]][multi.mer.list[[as.character(non)]]]})
	  }
	} else {
	  if (frequencies){
	    pwm = apply(nucleotide.mat,2,function(x){freqs = summary(factor(x,levels = multi.mer.list[[as.character(non)]]),maxsum = length(multi.mer.list[[as.character(non)]])+2)[multi.mer.list[[as.character(non)]]];freqs = freqs/sum(freqs)})
	  } else {
	    pwm = apply(nucleotide.mat,2,function(x){freqs = summary(factor(x,levels = multi.mer.list[[as.character(non)]]),maxsum = length(multi.mer.list[[as.character(non)]])+2)[multi.mer.list[[as.character(non)]]]})
	  }
	}
	colnames(pwm) = from:to

	return(pwm)
}



