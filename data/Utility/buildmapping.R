

### Build a mapping ###


buildmapping = function(martmat,from,to)
{
	martmat = martmat[,c(from,to)]
	testing = apply(martmat,1,function(x){paste(x[1],x[2])})
	rownames(martmat) = testing
	martmat = martmat[unique(testing),]
	mapping = martmat[,to]
	names(mapping) = martmat[,from]
	print(paste(from,"NAs:",sum(is.na(names(mapping)) | names(mapping) == "")))
	if (length(names(mapping)) > length(unique(names(mapping)))){print("Names of mapping are NOT unique")} else {print("Names of mapping are unique")}
	print(paste(to,"NAs:",sum(is.na(mapping) | mapping == "")))
	if (length(mapping) > length(unique(mapping))){print("Mapping is NOT unique")} else {print("Mapping is unique")}
	print(paste("Length:",length(mapping)))
	return(mapping)
}



