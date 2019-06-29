# Onesample FUNctions ------------------------------------------------------------


permu.onesample <- function(x, y = NULL, 
                            FUN = mean,
                            fetch = NULL,
                            nperm = 500, 
                            nboot = 0, 
                            ci.percent = .95, 
                            alternative = "two.sided", ...) {
  #' This function runs a one-sample permutation test to determine if the provided data is 
  #' significantly different from zero.
  #' 
  #' If just one vector is supplied (x), this function will determine if the test statistic of 
  #' interestfor that one vector is sifferent from zero. If two are provided (x & y), the test will 
  #' determine if the *paired* difference between those two vectors (x - y) differs from zero. 
  #' Therefore, if both x and y are provided, they must be of the same length.
  #' 
  #' @param x A numeric vector 
  #' @param y A second paired numeric vector (optional) 
  #' @param FUN The function to be performed on the vectors to yeild the test statistic of interest
  #' @param fetch The field containing the test statistic in the list returned from FUN (optional)
  #' @param nperm The number of permutations resamples to run (resample without replacement, default = 500)
  #' @param nboot The number of bootstrap resamples to run (resample with replacement, defautl = 0)
  #' @param ci.percent The desired confidence interval to determine both the alpha-threshold of 
  #'   the permutation test and the confidence interval of the bootstrap resampling (default = .95)
  #' @param alternative How significance of the permutation test is to be determined ("more", "less", 
  #'   default = "two.sided")
  #' @param ... Additional parameters to be passed into the user-specified FUN (not currently working)
  #' 
  #' @examples 
  #' permu.onesample(x, y, FUN = t.test, fetch = 'statistic', nperm = 1000, nboot = 500, alternative = 'more')
  
  # sanity checks
  
  alternative.methods <- c("two.sided", "less", "more")
  if (!alternative %in% alternative.methods) {
    stop("Must choose either 'two.dided' (default), 'more', or 'less' for the alternative.")
  }
  
  if (length(x) <= 1) {
    stop("x must have a length greater than 1.")
  }
  
  if (!is.numeric(x)) {
    stop("x must be a numeric vector.")
  }

  if (is.null(y)) {
    z <- x
  } else if (length(x) == length(y)) {
    z <- x - y
  } else {
    stop("Numeric vectors x and y must have the same length.")
  }
  
  # get the observed values
  observed.val <- FUN(z, paired = T)
  if (is.list(observed.val)) {
    observed.val <- observed.val[[fetch]]
  }
  
  # run the permutations
  perm.distro <- rep(NaN, nperm)
  for (pp in 1:nperm) {
    tmp.pp <- FUN(z * (1 - 2 * rbinom(size = 1, prob = .5, n = length(z))))
    
    if (is.list(tmp.pp)) {
      tmp.pp <- tmp.pp[[fetch]]
    }
    
    perm.distro[pp] <- tmp.pp
  }
  
  perm.pval <- perm.get_pval(observed = observed.val, distro = perm.distro, alternative = alternative)
  
  # running the bootstrapped tests
  
  # check.install.load.package(package = "boot", install = T, load = T)
  if (nboot > 0) {
    
    if (ci.percent <= 0 | ci.percent >= 1) {
      stop("Bootstrap conficence intervals, ci.percent, must be between 0 and 1 (non-inclusive).")
    }
    
    boot.results <- boot.onesample(x = z, nboot = nboot, FUN = FUN, ci.percent = ci.percent, fetch = fetch)
  }
  
  # compile the output
  output                          <- list()
  output$observed                 <- observed.val
  output$permutation$distribution <- perm.distro
  output$permutation$p.value      <- perm.pval
  
  if (nboot > 0) {
    output$bootstrap$distribution         <- boot.results$distro
    output$bootstrap$confidence.intervals <- boot.results$ci
  }
  
  return(output)
}

boot.onesample <- function(x, y = NULL, nboot = 1000, ci.percent = .95, FUN = mean, fetch = NULL) {
  
  percentiles <- c( (1 - ci.percent)/2,  ci.percent + (1 - ci.percent)/2)
  
  if (nboot < 1 | nboot%%1 != 0) {
    stop("nboot must be an integer greater than 1.")
  }
  
  if (is.null(y)) {
    z <- x
  } else if (length(x) == length(y)) {
    z <- x - y
  } else {
    stop("Numeric vectors x and y must have the same length.")
  }
  
  
  boot.distro <- rep(NaN, nboot)
  
  for(bb in 1:nboot) {
    
    boot.rand.z <- sample(1:length(z), size = length(z), replace = T)
    
    z.shuffle <- z[boot.rand.z]
    
    tmp.val <- FUN(z.shuffle, paired = T)
    if (is.list(tmp.val)) {
      tmp.val <- tmp.val[[fetch]]
    }
    
    boot.distro[bb] <- tmp.val
    
  }
  
  boot.ci <- quantile(boot.distro, percentiles)
  
  boot.output             <- list()
  boot.output[['distro']] <- boot.distro
  boot.output[['ci'    ]] <- boot.ci
  
  return(boot.output)
}

# Twosample functions ------------------------------------------------------------------------------

permu.twosample <- function(x, y, 
                            FUN = mean.diff,
                            fetch = NULL,
                            nperm = 500, 
                            nboot = 0, 
                            ci.percent = .95, 
                            alternative = "two.sided", ...) {
  #' This function runs a two-sample permutation test to determine if the provided data is 
  #' significantly different from zero.
  #' 
  #' If just one vector is supplied (x), this function will determine if the test statistic of 
  #' interestfor that one vector is sifferent from zero. If two are provided (x & y), the test will 
  #' determine if the difference between those two vectors (x - y) differs from zero. 
  #' 
  #' @param x A numeric vector 
  #' @param y A numeric vector
  #' @param FUN The function to be performed on the vectors to yeild the test statistic of interest
  #' @param fetch The field containing the test statistic in the list returned from FUN (optional)
  #' @param nperm The number of permutations resamples to run (resample without replacement, default = 500)
  #' @param nboot The number of bootstrap resamples to run (resample with replacement, defautl = 0)
  #' @param ci.percent The desired confidence interval to determine both the alpha-threshold of 
  #'   the permutation test and the confidence interval of the bootstrap resampling (default = .95)
  #' @param alternative How significance of the permutation test is to be determined ("more", "less", 
  #'   default = "two.sided")
  #' @param ... Additional parameters to be passed into the user-specified FUN (not currently working)
  #' 
  #' @examples 
  #' permu.onesample(x, y, FUN = t.test, fetch = 'statistic', nperm = 1000, nboot = 500, alternative = 'more')
  
  # sanity checks
  
  alternative.methods <- c("two.sided", "less", "more")
  if (!alternative %in% alternative.methods) {
    stop("Must choose either 'two.dided' (default), 'more', or 'less' for the alternative.")
  }
  
  if (length(x) <= 1) {
    stop("x must have a length greater than 1.")
  }
  
  if (length(y) <= 1) {
    stop("y must have a length greater than 1.")
  }
  
  if (!is.numeric(x)) {
    stop("x must be a numeric vector.")
  }
  
  if (!is.numeric(y)) {
    stop("y must be a numeric vector.")
  }
  
  
  # get the observed values
  observed.val <- FUN(x, y, paired = F)
  if (is.list(observed.val)) {
    observed.val <- observed.val[[fetch]]
  }
  
  # run the permutations
  perm.distro <- rep(NaN, nperm)
  for (pp in 1:nperm) {
    data <-       data.frame(label = 'x', val = x)
    data <- rbind(data.frame(label = 'y', val = y), data)
    
    data$label <- sample(data$label, size = dim(data)[1], replace = F)
    
    x.sample <- data$val[data$label=='x']
    y.sample <- data$val[data$label=='y']
    
    tmp.pp <- FUN(x = x.sample, y = y.sample, paired = F)
    if (is.list(tmp.pp)) {
      tmp.pp <- tmp.pp[[fetch]]
    }
    
    perm.distro[pp] <- tmp.pp
  }
  
  perm.pval <- perm.get_pval(observed = observed.val, distro = perm.distro, alternative = alternative)
  
  # running the bootstrapped tests
  
  # check.install.load.package(package = "boot", install = T, load = T)
  if (nboot > 0) {
    
    if (ci.percent <= 0 | ci.percent >= 1) {
      stop("Bootstrap conficence intervals, ci.percent, must be between 0 and 1 (non-inclusive).")
    }
    
    boot.results <- boot.twosample(x, y, nboot = nboot, FUN = FUN, ci.percent = ci.percent, fetch = fetch)
  }
  
  # compile the output
  output                          <- list()
  output$observed                 <- observed.val
  output$permutation$distribution <- perm.distro
  output$permutation$p.value      <- perm.pval
  
  if (nboot > 0) {
    output$bootstrap$distribution         <- boot.results$distro
    output$bootstrap$confidence.intervals <- boot.results$ci
  }
  
  return(output)
}

boot.twosample <- function(x, y, nboot = 1000, ci.percent = .95, FUN = mean.diff, fetch = NULL) {
  
  percentiles <- c( (1 - ci.percent)/2,  ci.percent + (1 - ci.percent)/2)
  
  if (nboot < 1 | nboot%%1 != 0) {
    stop("nboot must be an integer greater than 1.")
  }
  
  boot.distro <- rep(NaN, nboot)
  
  for(bb in 1:nboot) {
    
    x.sample <- sample(x, size = length(x), replace = T)
    y.sample <- sample(y, size = length(y), replace = T)
    
    tmp.val <- FUN(x = x.sample, y = y.sample, paired = F)
    if (is.list(tmp.val)) {
      tmp.val <- tmp.val[[fetch]]
    }
    
    boot.distro[bb] <- tmp.val
    
  }
  
  boot.ci <- quantile(boot.distro, percentiles)
  
  boot.output             <- list()
  boot.output[['distro']] <- boot.distro
  boot.output[['ci'    ]] <- boot.ci
  
  return(boot.output)
}

# support ------------------------------------------------------------------------------------------

mean.diff <- function(x, y, ...) {
  if (!is.numeric(x)) {
    stop('x must be a numeric vector')
  } else if (!is.numeric(y)) {
    stop('y must be a numeric vector.')
  }
  
  z <- mean(x) - mean(y)
  
  return(z)
}

perm.get_pval <- function(observed, distro, alternative = "two.sided") {
  
  p <- mean(distro >= observed)
  
  if (alternative == "more") {
    p <- p
    
  } else if (alternative == "less") {
    p <- 1 - p
    
  } else if (alternative == "two.sided") {
    p <- 2 * min(p, 1 - p)
    
  } else {
    stop("Must choose either 'two.sided' (default), 'more', or 'less' for the alternative.")
  }
  
  return(p)
  
}

check.install.load.package <- function(package, install = F, load = F) {
  
  pkg.status <- package %in% installed.packages()[, "Package"]
  
  cat("Package '", package, "' installed: ", pkg.status, "\n", sep = "")
  
  if(install & !pkg.status) {
    install.packages(package)
  } 
  
  if (load) {
    require(package, character.only = T)
  }
}
