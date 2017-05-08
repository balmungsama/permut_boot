f <- function(d, i){
  d2 <- d[i,]
  
  # if(boot.count>1) (browser())

  d2.x1 <- d2$values[which(d2$labels=='x1')]
  d2.x2 <- d2$values[which(d2$labels=='x2')]
  
  value <- fun(d2.x1, d2.x2)
  value <- as.numeric(value)
  
  # boot.count <<- boot.count + 1
  
  return(value)
}

boot.CI.mean <- function(X.in, fun, R, Use, ...){
  require(boot)
  boot.stat <- boot(data = X.in, statistic = f, R = R, strata = X.in$labels)
  return(boot.stat)
}