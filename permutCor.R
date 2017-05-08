fcorrel <- function(d, i){
  d2 <- d[i,]
  d2.x1 <- d2$x1
  d2.x2 <- d2$x2
  
  value <- fun(d2.x1, d2.x2)
  value <- as.numeric(value)

  return(value)
}

boot.CI.correl <- function(X.in, fun, R, Use, ...){
  require(boot)
  boot.stat <- boot(X.in, statistic = fcorrel, R = R)
  return(boot.stat)
}

permuCor <- function(x1, x2, 
											fun = cor, 
											nsims = 1000, p = 0.05, var.names = NULL, plot = TRUE, boot.CI.correl = TRUE, CI.thresh = 0.95, plot.title = NULL, ...) {
	require(ggplot2)
	
	if(is.null(var.names)){
		var.names <- c('x1', 'x2')
	}
	
	# browser()
	na_ind <- NULL
	if(sum(is.na(x1)) > 0) {
		na_ind <- c(na_ind, which(is.na(x1)))
	}
	if(sum(is.na(x2)) > 0) {
		na_ind <- c(na_ind, which(is.na(x2)))
	}
	
	if( length(na_ind) > 0 ) {
		x1 <- x1[-na_ind]
		x2 <- x2[-na_ind]
	}
	

	fun <<- fun
	# Use <<- Use
	obs.val <- fun(x1, x2)
	
	pool.vals <- c(x1, x2)
	pool.labs <- c(rep('x1', length(x1)), rep('x2', length(x2)))
	
	out.sim  <- NULL
	for(sim in 1:nsims){
		
		# labs.sim <- sample(pool.labs, replace = FALSE)
		# tab.sim  <- data.frame(names = labs.sim, values = pool.vals)
		pool.vals[which(pool.labs == 'x2')] <- sample(pool.vals[which(pool.labs == 'x2')], replace = FALSE)
		tab.sim  <- data.frame(names = pool.labs, values = pool.vals)
		
		x1.sim   <- tab.sim$values[which(tab.sim$names=='x1')]
		x2.sim   <- tab.sim$values[which(tab.sim$names=='x2')]
		
		out.sim[sim] <- fun(x1.sim, x2.sim)
		
		# cat(x2.sim, '\n')
	}
	
	# compute the bootstrap confidence intervals for the function of interest
	if(boot.CI.correl == T){
		# pool.tab  <- data.frame(names = pool.labs, values = pool.vals)
		boot.tab  <- data.frame(x1 = x1, x2 = x2)
		
		boot.stat <- boot.CI.correl(X.in = boot.tab, fun = fun, R = nsims)
		boot.CI.correl   <- boot.CI.correl(boot.stat, type='perc', conf = CI.thresh)
	}
	rm(fun, envir = .GlobalEnv)
	
	if(plot==TRUE){
		boot.xmin <- boot.CI.correl$percent
		boot.xmax <- boot.CI.correl$percent
		boot.xmin <- as.numeric(boot.xmin[length(boot.xmin)-1])
		boot.xmax <- as.numeric(boot.xmax[length(boot.xmax)  ])
		
		xmax <- boot.xmax
		xmin <- boot.xmin
		fig1 <- ggplot() +
			geom_histogram(mapping = aes(x = out.sim), bins=30, alpha = 1) +
			geom_rect(aes(xmin = boot.xmin,
										xmax = boot.xmax,
										ymin = 0,
										ymax = Inf), fill='pink', alpha = 0.3) +
			geom_vline(xintercept = c(boot.xmin, 
																boot.xmax), 
								 size = 1, colour='pink', alpha=0.5) +
			xlab('test statistic') +
			geom_vline(xintercept = obs.val, colour = 'red', alpha = 1, size = 1) +

			geom_vline(xintercept = obs.val, colour = 'red', alpha = 1) + 
			annotate('text', x = obs.val, y = nsims/100, label = signif(obs.val, digits = 3), fontface=2, size=5, colour='gray90') +
			theme_minimal() +
			theme(axis.line.x = element_line(size = 0.5, colour = "black"),
						axis.line.y = element_line(size = 0.5, colour = "black"),
						panel.grid.major = element_line(size = 0, colour = 'white'),
						panel.grid.minor = element_line(size = 0, colour = 'white')) +
			scale_y_continuous(expand = c(0,0)) +
			ggtitle(plot.title)
		
		print(fig1)
	}
	
	# out.sim <<- out.sim
	p.val <- mean(abs(out.sim) > abs(obs.val))

	if(p.val < p){
		sig.lab <- '  *  '
	} else {
		sig.lab <- '     '
	}
	
	theList <- list(title = plot.title,
									var.names = var.names, 
									test.stat = obs.val,
									boot.conf.int = data.frame(Level = paste0(CI.thresh*100, '%   '), 
																						 Percentile = paste0('(', signif(boot.xmin, 4), ', ', signif(boot.xmax, 4), ')' ) ),
									significance = data.frame(p.val = p.val, sig = sig.lab))
	
	return(theList)
}