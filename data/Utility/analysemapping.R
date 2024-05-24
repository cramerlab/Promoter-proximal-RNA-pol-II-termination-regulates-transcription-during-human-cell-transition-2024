

### Analyse the mapping ###


analysemapping = function(vect,from,to,analysis.path,check = FALSE)
{	
	feat = switchvector(vect)
	unifeat = unique(feat)
	write.table(paste0(to,"_for_",from),file = file.path(analysis.path,paste0(to,"_for_",from,".txt")),append = FALSE, sep = "\t",col.names=F,row.names=F)
	for(i in 1:length(unifeat)) {vec = feat[names(feat) %in% names(which(feat == unifeat[i]))]
		if(length(vec) > 1){vecnames = vec
			vec = names(vec)
			names(vec) = vecnames
			vec = unique(vec)
			vec = as.matrix(vec)
			colnames(vec) = unifeat[i]
			rownames(vec) = NULL
			if (!all(vec[1] == vec)){write.table(t(as.matrix(vec)),file = file.path(analysis.path,paste0(to,"_for_",from,".txt")),append = TRUE, sep = "\t",col.names=F,row.names=T);if(check){cat("+ \n")}} else if(check){cat("* \n")}
		} else if(check){cat(". \n")}
	}
	
	feat = vect
	unifeat = unique(feat)
	write.table(paste0(from,"_for_",to),file = file.path(analysis.path,paste0(from,"_for_",to,".txt")),append = FALSE, sep = "\t",col.names=F,row.names=F)
	for(i in 1:length(unifeat)) {vec = feat[names(feat) %in% names(which(feat == unifeat[i]))]
		if(length(vec) > 1){vecnames = vec
			vec = names(vec)
			names(vec) = vecnames
			vec = unique(vec)
			vec = as.matrix(vec)
			colnames(vec) = unifeat[i]
			rownames(vec) = NULL
			if (!all(vec[1] == vec)){write.table(t(as.matrix(vec)),file = file.path(analysis.path,paste0(from,"_for_",to,".txt")),append = TRUE, sep = "\t",col.names=F,row.names=T);if(check){cat("+ \n")}} else if(check){cat("* \n")}
		} else if(check){cat(". \n")}
	}
}



