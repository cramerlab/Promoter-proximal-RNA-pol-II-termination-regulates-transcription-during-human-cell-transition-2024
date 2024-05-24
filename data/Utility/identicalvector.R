

### identicalvector function ###


identicalvector = function(x,y)
{
	all(all(x %in% y),all(y %in% x))
}



