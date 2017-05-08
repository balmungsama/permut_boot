f <- function(d, i){
  d2 <- d[i,]
  d2.x1 <- d2$x1
  d2.x2 <- d2$x2
  
  value <- fun(d2.x1, d2.x2)
  value <- as.numeric(value)

  return(value)
}

boot.CI <- function(X.in, fun, R, Use, ...){
  require(boot)
  boot.stat <- boot(X.in, statistic = f, R = R, strata = X.in[ , c(1,2)] )
  return(boot.stat)
}