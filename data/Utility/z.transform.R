

### z.transform function: transform the columns of a matrix to z-scores ###

z.transform = function(mat){
  mat = t(mat)
  mat = (mat - apply(mat,1,mean,na.rm = TRUE))/apply(mat,1,sd,na.rm = TRUE)
  return(t(mat))
}



