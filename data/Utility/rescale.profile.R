

### rescale.profile function ###


rescale.profile = function(profile,scale.to.length,summarize=function(x){median(x,na.rm=TRUE)}){
	breaks = round(seq(1,length(profile),length.out=scale.to.length+1))
	scaled.profile = c(summarize(profile[c(breaks[1]:breaks[2])]),as.numeric(running(breaks[-1],fun=function(y){summarize(profile[c((y[1]+1):y[2])])},w=2,align="left")))
	return(scaled.profile)
}



