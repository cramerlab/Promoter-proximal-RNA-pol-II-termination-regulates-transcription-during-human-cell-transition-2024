

### is.finite.min function ###


is.finite.min = function(vec,na.rm = TRUE)
{
	return(min(vec[is.finite(vec)],na.rm = na.rm))
}


### is.finite.max function ###


is.finite.max = function(vec,na.rm = TRUE)
{
  return(max(vec[is.finite(vec)],na.rm = na.rm))
}



