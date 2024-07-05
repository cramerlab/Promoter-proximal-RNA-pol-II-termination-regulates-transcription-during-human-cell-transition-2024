

### auc function ###


auc = function (x,y,method = c("trapezoid","step","spline")) 
{
  idx <- order(x)
  x <- x[idx]
  y <- y[idx]
  switch(match.arg(arg = method,choices = c("trapezoid","step","spline")),
         trapezoid = {a <- sum((apply(cbind(y[-length(y)],y[-1]),1,mean))*(x[-1] - x[-length(x)]))},
         step = {a <- sum(y[-length(y)] * (x[-1] - x[-length(x)]))},
         spline = {a <- integrate(splinefun(x,y,method = "natural"),lower = min(x),upper = max(x))$value})
  return(a)
}



