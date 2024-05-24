

### resize.profile function ###


resize.profile = function(profile,scale.to.length){
	length.profile = length(profile)
	ratio = max(round(length.profile/scale.to.length),1)
	if (ratio > 1){profile = filter(c(rep(profile[1],ratio),profile,rep(profile[length.profile],ratio)),rep(1,ratio)/ratio)[(ratio+1):(length.profile+ratio)]}
	resized.profile = profile[round(seq(1,length.profile,length.out=scale.to.length))]
	return(resized.profile)
}



