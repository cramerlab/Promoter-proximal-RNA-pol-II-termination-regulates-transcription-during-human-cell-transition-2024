

### checkvector function ###


checkvector = function(vec)
{
  check = c("REAL" = sum(!is.na(vec) & is.finite(vec)),
            "NA" = sum(is.na(vec)),
            "INF" = sum(is.infinite(vec)),
            "SUM" = sum(!is.na(vec) & is.finite(vec)) + sum(is.na(vec)) + sum(is.infinite(vec)),
            "length" = length(vec))
	return(check)
}



