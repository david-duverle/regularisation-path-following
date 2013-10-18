step.length.cox <- function(corrector, depth, lambda2, exprs.all, d, rslist, wlist,
                            min.lambda, max.arclength, frac.arclength,
                            add.newvars, backshoot, approx.Gram, h0 = NULL,
                            eps = .Machine$double.eps)
{
    active <- corrector$active

    force.active <- NULL

    lambda <- corrector$lambda - min.lambda

    k <- length(active)
    Lambda2 <- rep(lambda2, length(active))
    b <- corrector$b[active]
    
    eta <- corrector$eta
    wsum <- corrector$wsum
    n <- length(d)
    A <- rep(0, n)
    if (!approx.Gram) 
        AA <- matrix(0, n, n)
    rset <- w <- NULL

    for (i in 1:sum(d == 1)) {
      if (!is.null(rslist[[i]])) {
        rset0 <- rset
        rset <- c(rset0, rslist[[i]])
      }
      w <- c(rep(1, length(rset0)), wlist[[i]]) * eta[rset]
      w1 <- w / wsum[i]
      A[rset] <- A[rset] + w1 - w1^2
      if (!approx.Gram) {
        k <- length(rset)
        AA[1:k, 1:k] <- AA[1:k, 1:k] - outer(w1[rset], w1[rset])
      }
    }
   
    if (approx.Gram)
		dhdb.x <- t(get.active.values(exprs.all, active) * A)
    else {
		diag(AA) <- A
		dhdb.x <- t(get.active.values(exprs.all, active)) %*% AA
    }

   C.active = apply(get.active.values(exprs.all, active) * corrector$a, 2, sum)
   # C = apply(exprs.all * corrector$a, 2, sum)

    if (length(active) == 1)
        db <- sign(C.active) / (dhdb.x %*% get.active.values(exprs.all, active) + lambda2)
    else {
        db <- (solve(dhdb.x %*% get.active.values(exprs.all, active) + diag(Lambda2)) %*% sign(C.active))
    }

    newa <- NULL
    if (!backshoot) {

      w.1 = - corrector$a
      w.2 = - t(db) %*% t(get.active.values(exprs.all, active)) %*% AA

      hd <- -b / db
      
      if(length(hd[hd > eps]) > 0)
         curmin = min(hd[hd > eps])
      else
         curmin = lambda
         		
      active.array = gene.set.name.to.num.array(active, colnames(exprs.all), depth)
      idxs = as.integer(rep(0, depth))
		
      res = .C("get_min_itemset", as.integer(depth), as.integer(exprs.all), as.integer(nrow(exprs.all)), as.integer(ncol(exprs.all)), as.double(w.1), as.double(w.2), as.integer(t(active.array)), as.integer(nrow(active.array)), as.double(corrector$lambda), as.double(curmin), idxs, as.integer(global.trace))
      
      h = as.double(res[10])
      min.idxs = as.integer(res[[11]])
      
      if(min.idxs[1] > 0)
         names(h) = make.gene.set.name(colnames(exprs.all)[min.idxs[min.idxs > 0]])
      
		if(global.trace) {
			cat("\rStep length (min): ", h, sep="");
			if(length(names(h)))
			 	cat(" (", names(h), ")", sep="");
			cat("\n")
		}
      if(min.idxs[1] > 0 && h < lambda) {
   		# if(names(h) %in% active) {
   		#    print("### Argmin is in active set already")
   		#    browser()
   		# }

			newa = names(h)
      }
      h <- min(h * frac.arclength, max.arclength / sum(abs(db)))
    } 
    else {
      hd <- b / db
      ii <- hd > eps & hd < -h0
      if (any(ii)) 
        h <- -max(hd[ii])
      else 
        h = 0
    }
    list(h = -h, db = db, newa = newa, arclength = h * sum(abs(db)))
}


predictor.cox <- function(b, step)
{
	b - step$db * step$h
}

get.active.values <- function(exprs.all, subsets) {
	if(class(subsets) != 'character') {
		cat("Subsets provided is not a string:", subsets, "\n")
		browser()
	}
   vals = sapply(strsplit(subsets, "*", fixed = TRUE), 
		function(feature) { 
			if(all(feature %in% colnames(exprs.all))) { 
				if(length(feature) > 1) 
						apply(exprs.all[,feature], 1, all) 
					else 
						exprs.all[,feature]; 
			} 
			else
				rep(FALSE, nrow(exprs.all)) 
		})
   colnames(vals) = subsets
   
   return(vals)
}

make.gene.set.name <- function(obj1, obj2 = NULL) {
   if(!is.null(obj2) & length(obj1) > 1)
      return(sapply(obj1, make.gene.set.name, USE.NAMES=FALSE, obj2))
   if(length(obj2) > 1)
      return(sapply(obj2, make.gene.set.name, USE.NAMES=FALSE, obj1))
      
   if(is.null(obj2))
      arr = obj1
   else
      arr = c(obj1, obj2)
   arr = arr[order(arr)]
   return(paste(arr, sep="", collapse="*"))
}

gene.set.name.to.num.array <- function(gene.sets, gene.names, depth = 3) {
	l = lapply(gene.sets, function(gene.set) {
		ids = sapply(strsplit(gene.set, "*", fixed=T)[[1]], function(gene) {
			which(gene.names == gene)
		})
		c(as.numeric(ids[order(ids)]), rep(0, depth-length(ids)))
	})
	matrix(unlist(l), ncol=depth, byrow=T)
}

corrector.cox <- function(exprs.all, depth, d, rslist, wlist, rept, method, active, tmpa,
                          force.active, lambda, lambda2, b0.tmpa,
                          bshoot.threshold, relax.lambda,
                          eps = .Machine$double.eps)
  {
    nobs <- nrow(exprs.all)
    p <- length(tmpa)
 	
    if (p > 0) {
		b2 <- c(pmax(b0.tmpa, 0), -pmin(b0.tmpa, 0))

		# xa <- x[, tmpa, drop = FALSE]
		xa = get.active.values(exprs.all, tmpa)

		penalty <- rep(1, p)
		z <- c(as.vector(xa), d, rept, penalty, rep(0, nobs), rep(0, nobs))
		mz <- c(nobs, method, lambda, lambda2, 0)
		
		sol <- .C('solve_coxpath',
		          as.integer(2 * p),
		          as.double(b2),
		          as.double(rep(0, 2 * p)),
		          as.double(rep(1e300, 2 * p)),
		          as.integer(0),
		          as.double(z),
		          as.double(mz))
		if (sol[[5]] != 0) {
		  cat('Convergence warning\n')
		}

		b0.tmpa <- sol[[2]][1:p] - sol[[2]][(p + 1):(2 * p)]
		names(b0.tmpa) = tmpa
		i <- (p + 2) * nobs + p
		eta <- sol[[6]][i + c(1:nobs)]
		i <- i + nobs
		wsum <- sol[[6]][i + c(1:nobs)][d == 1]
		lp <- sol[[7]][5]
    }
    else {
      eta <- rep(1, nobs)
      wsum <- rep(0, sum(d))
    }
    
    rset <- NULL
    a <- d == 1
    for (i in 1:sum(a)) {
      if (!is.null(rslist[[i]])) {
        rset0 <- rset
        rset <- c(rset0, rslist[[i]])
      }
      w <- c(rep(1, length(rset0)), wlist[[i]]) * eta[rset]
      if (wsum[i] == 0) 
        wsum[i] <- sum(w)
      a[rset] <- a[rset] - w / wsum[i]
    }
    
    if (p == 0) 
        lp <- -sum(log(wsum))
        	 
	 filename.prefix = paste('.temp.pid', Sys.getpid(), '.', sep='')
	 
	 write(a, paste(filename.prefix, "weights.dat", sep=''), sep="\n")
	 write(-a, paste(filename.prefix, "neg_weights.dat", sep=''), sep="\n")
    
    if (p > 0) {
		cmd = paste("./lcm C -l 1 -u ", depth, " -w ", filename.prefix, "weights.dat ", filename.prefix, "input.dat ", lambda * (1 - relax.lambda), " ", filename.prefix, "output.dat", sep='')
		output = system(cmd, intern = TRUE, ignore.stdout = TRUE, ignore.stderr = TRUE)
				
		fc <- file(paste(filename.prefix, "output.dat", sep=''))
		vals = strsplit(readLines(fc), " ")
		close(fc)
		
		vals = sapply(vals, function(item) { make.gene.set.name(colnames(exprs.all)[as.numeric(item)]) })
		
		newa <- as.character(vals[! vals %in% tmpa])
		newactive <- as.character(vals[! vals %in% active])
		
		cmd = paste("./lcm CK -l 1 -u ", depth, " -w ", filename.prefix, "neg_weights.dat ", filename.prefix, "input.dat ", lambda * (1 - relax.lambda), " ", filename.prefix, "output.dat", sep='')
	   
		output = system(cmd, intern = TRUE, ignore.stdout = TRUE, ignore.stderr = TRUE)
		
		fc <- file(paste(filename.prefix, "output.dat", sep=''))
		vals = strsplit(readLines(fc), " ")
		close(fc)
		
		vals = sapply(vals, function(item) { make.gene.set.name(colnames(exprs.all)[as.numeric(item)]) })
		if(length(vals) > 0) {
			if(class(newactive) == 'list') #DEBUG
				browser()
			
			newa <- unique(append(newa, as.character(vals[! vals %in% tmpa])))
			newactive <- unique(append(newactive, as.character(vals[! vals %in% active])))
			
			if(class(newactive) == 'list') #DEBUG
				browser()
			
   	}
   	
		if(sum(newa %in% tmpa) > 0) {
		    cat("ERROR: non-unique elements")
		    cat(newa)
		    cat(tmpa)
		    browser()
		}
		
		# i <- which(abs(corr) >= lambda * (1 - relax.lambda))
		
		i <- which(abs(b0.tmpa[active]) < eps) 
		inactive <- active[i]
        if(length(i) > 0) {
            # b0.tmpa = b0.tmpa[-i]
     		active = active[!active %in% inactive]
            b0.tmpa[i] = 0
        }
		if(length(newactive) > 0) {
		    new.b = rep(0, length(newactive))
		    names(new.b) = newactive
		    b0.tmpa = append(b0.tmpa, new.b)
		}
    	    
		active <- append(active, newactive)
    }
	else {
		cmd = paste("./lcm C -K 1 -l 1 -u 3 -w ", filename.prefix, "weights.dat ", filename.prefix, "input.dat 0", sep='') 
		
		max.val = as.numeric(system(cmd, intern = TRUE, ignore.stdout = F, ignore.stderr = TRUE))
				
		lambda = 0
		how.many.results = 1
		while(max.val > lambda + 0.001) { #HACK!
			how.many.results = how.many.results*100
			cmd = paste("./lcm Ff -# ", how.many.results, " -l 1 -u 3 -w ", filename.prefix, "weights.dat ", filename.prefix, "input.dat 0 ", filename.prefix, "output.dat", sep="") # fine-tune maximum num of results parameter (curr. 1000)
			output = system(cmd, intern = TRUE, ignore.stdout = TRUE, ignore.stderr = TRUE)
		
			fc <- file(paste(filename.prefix, "output.dat", sep=""))
			vals = strsplit(readLines(fc), " ")
			close(fc)
				
			lambdas = lapply(vals, function(val) { as.double(substr(val[[length(val)]], 2, nchar(val[[length(val)]])-1)) })
			i = which.max(lambdas)
			lambda = lambdas[[i]]
			c1 = as.integer(vals[[i]][1:length(vals[[i]])-1])
		}
		
		cmd = paste("./lcm C -K 1 -l 1 -u 3 -w ", filename.prefix, "neg_weights.dat ", filename.prefix, "input.dat 0", sep='') 
		
		max.val = as.numeric(system(cmd, intern = TRUE, ignore.stdout = F, ignore.stderr = TRUE))
		
		how.many.results = 1
		while(max.val > lambda + 0.001) { #HACK!
			how.many.results = how.many.results*100
		
			cmd = paste("./lcm Ff -# ", how.many.results, " -l 1 -u 1 -w ", filename.prefix, "neg_weights.dat ", filename.prefix, "input.dat ", lambda , " ", filename.prefix, "output.dat", sep='')
			output = system(cmd, intern = TRUE, ignore.stdout = TRUE, ignore.stderr = TRUE)

			fc <- file(paste(filename.prefix, "output.dat", sep=''))
			vals = strsplit(readLines(fc), " ")
			close(fc)
			lambdas = lapply(vals, function(val) { as.double(substr(val[[length(val)]], 2, nchar(val[[length(val)]])-1)) })
			i = which.max(lambdas)

			if(length(i) == 1 && lambdas[[i]] > lambda) {
				lambda = lambdas[[i]]
				c1 = as.integer(vals[[i]][1:length(vals[[i]])-1])
			}
		}
		
		newa <- newactive <- colnames(exprs.all)[c1]
		inactive <- NULL
		active <- append(active, newactive)
		
		if(length(c1) > 0) {
	    	new.b = rep(0, length(newactive))
    		names(new.b) = newactive
    		b0.tmpa = c(b0.tmpa, new.b)
    	}
    }
    
    df <- length(active) - length(newactive)
	
 	if(length(active) == 0) #DEBUG
		browser()
	if(class(active) == 'list') #DEBUG
		browser()
	
	if(length(newactive) == 0)
		newactive = NULL
	else if(class(newactive) == 'list')
		browser()
   backshoot <- (any(abs(b0.tmpa[newactive]) > bshoot.threshold))
    
    list(eta = eta, wsum = wsum, b = b0.tmpa, lp = lp, active = active,
         force.active = force.active, newactive = newactive, newa = newa,
         inactive = inactive, corr = NULL, lambda = lambda, df = df,
         backshoot = backshoot, a = a)
}

extract.features <- function(exprs, extra.features = NULL, deviation.threshold = NULL, quantile = NULL) {
   if(is.null(deviation.threshold))
      deviation.threshold = 1.5
   
   # extracting features sets:
	m <- ncol(exprs)

	if(is.null(quantile)) { 
		exprs.sd = apply(exprs, 2, sd, na.rm=T)
		exprs.sd = do.call("rbind", replicate(nrow(exprs), exprs.sd, simplify = FALSE))
		exprs.mean = apply(exprs, 2, mean, na.rm=T)
		exprs.mean = do.call("rbind", replicate(nrow(exprs), exprs.mean, simplify = FALSE))
	
		exprs.std = (exprs - exprs.mean)/exprs.sd
		exprs.expressed = (exprs.std >= deviation.threshold)
		exprs.repressed = (exprs.std <= -deviation.threshold)
	}
	else {
		n = nrow(exprs)
		exprs.order = apply(exprs, 2, order)
		exprs.expressed = matrix(FALSE, n, m)
		exprs.repressed = matrix(FALSE, n, m)
		for(i in 1:m) {
			exprs.repressed[exprs.order[1:(n/quantile), i], i] = TRUE;
			exprs.expressed[exprs.order[(n * (quantile-1)/quantile):n, i], i] = TRUE;
		}
	}
	
	exprs.all = cbind(exprs.repressed, exprs.expressed)
	colnames(exprs.all) = c(paste('dn.', seq(m), sep=''), paste('up.', seq(m), sep=''))
	
   s = apply(exprs.all, 2, any)
	exprs.all = exprs.all[, s > 0 & !is.na(s)]
   exprs.all[is.na(exprs.all)] = F
	
	m = ncol(exprs.all)
   
	if(!is.null(extra.features)) {
   	if(nrow(extra.features) != nrow(exprs)) {
   	   stop("Wrong input: extra features need to have same number of row as gene expression matrix")
   	}
   	exprs.all = cbind(exprs.all, extra.features)
   	colnames(exprs.all)[seq(m+1, m+ncol(extra.features))] = colnames(extra.features)
	}
   
   return(exprs.all)
}


itemset.coxpath <- function(data, gene.symbols = NULL, depth = 3, nopenalty.subset = NULL,
                    method = c('breslow', 'efron'), lambda2 = 1e-5,
                    min.steps = 10, max.steps = 10 * min(n, m), max.norm = 100 * m,
                    min.lambda = (if (m >= n) 1e-3 else 0), max.vars = Inf,
                    max.arclength = Inf, frac.arclength = 1, add.newvars = 1,
                    bshoot.threshold = 0.1, relax.lambda = 1e-7,
                    approx.Gram = FALSE, standardize = TRUE,
                    eps = .Machine$double.eps, trace = FALSE, deviation.threshold = NULL, quantile = NULL,
						  max.time.per.step = 30 * depth)
{
	dyn.load("solve_coxpath.so")
	dyn.load("explore_itemset.so")
	
   assign("global.trace", trace, envir = .GlobalEnv)
	filename.prefix = paste('.temp.pid', Sys.getpid(), '.', sep='')
   
	if(is.null(gene.symbols) && length(data) == 2 && class(data) == 'list' && ncol(data[[1]]$x) == nrow(data[[2]])) {
		gene.symbols = data[[2]]
		data = data[[1]]
	}
	call <- match.call()
	method <- match.arg(method)
	mthd <- switch(method, breslow = 1, efron = 2)
	exprs <- data$x

	time <- data$time
	status <- data$status

	n <- length(time)
		
	if(nrow(exprs) > ncol(exprs)) {
		print("Sanity check failed: input data has wrong dimensions: more rows than columns.")
		return(NULL);
	}
	
	exprs.all = extract.features(exprs, data$extra.features, deviation.threshold, quantile)
   if(sum(is.na(exprs.all)) > 0) {
      print("ERROR: exprs.all contains NA values")
      browser()
   }
   
   m = dim(exprs.all)[2]
   
	o <- order(status)
	oo <- o[order(time[o], decreasing = TRUE)]
	exprs.all <- exprs.all[oo, ]
	time <- time[oo]
	status <- status[oo]
	complete <- which(status == 1)
	nnc <- length(complete)
	rept <- rep(0, n)
   
	exprs.sets = apply(exprs.all, 1, which)
	
	unlink(paste(filename.prefix, "input.dat", sep=''))
	lapply(exprs.sets, function(line) { write(line, paste(filename.prefix, "input.dat", sep=''), n=dim(exprs.all)[2], append=TRUE)} );


    for (i in complete)
        rept[i] <- sum(time[i:n] == time[i] & status[i:n] == 1)
    rslist <- wlist <- vector('list', length = nnc)
    
    for (i in 1:nnc) {
      if (i == 1) {
        ii <- time >= time[complete[1]]
        rslist[[1]] <- which(ii)
      } else if (rept[complete[i]] >= rept[complete[i] - 1]) {
        ii <- (time >= time[complete[i]]) & (time < time[complete[i - 1]])
        rslist[[i]] <- which(ii)
      }
      
      wlist[[i]] <- rep(1, sum(ii))
      if (mthd == 2) {
        if (rept[complete[i]] > 0) {
          tie <- time[ii] == time[complete[i]] & status[ii] == 1
          di <- max(rept[ii][tie])
          wlist[[i]][tie] <- wlist[[i]][tie] - (di - rept[complete[i]]) / di
        }
      }
    }
    
    if (frac.arclength > 1 || frac.arclength <= 0) {
      frac.arclength <- 1
      cat('frac.arclength should be in (0,1]. frac.arclength is set to 1.\n')
    }
    if (max.arclength < Inf && frac.arclength < 1) {
      frac.arclength <- 1
      cat(paste('frac.arclength<1 can be used only if max.arclength=Inf.',
                'frac.arclength is set to 1.\n'))
    }
    
    n.repeat <- n.repeat1 <- ceiling(1 / frac.arclength)    

    lam.vec <- step.len <- rep(0, max.steps)
    # bmat.pred <- bmat.corr <- cmat <- matrix(0, nrow = max.steps, ncol = ncol(exprs.all))
    bmat.pred <- bmat.corr <- cmat <- list();
    lp <- df <- new.df <- rep(0, max.steps)
    new.A <- rep(FALSE, max.steps)
    actions <- vector('list', length = max.steps)
    backshoot <- FALSE
    b.tmpa <- NULL
    
    force.active <- NULL
    
	 
    corrector <- corrector.cox(exprs.all, depth, status, rslist, wlist, rept, mthd,
                               nopenalty.subset, nopenalty.subset,
                               force.active, 0, lambda2, b.tmpa,
                               bshoot.threshold, relax.lambda)

    k <- 1
    b.tmpa = bmat.pred[[k]] = bmat.corr[[k]] = corrector$b
    
    lam.vec[k] <- lambda <- corrector$lambda
    new.df[k] <- df[k] <- corrector$df
    lp[k] <- corrector$lp
    new.A[k] <- TRUE
    
    actions[[k]] <- active <- corrector$active

    cox.print.debug('Lambda = ', lambda, '. Let the first factor in.\n#', k)
    cox.print.debug('\t', active, ' added')

    if (max.steps <= 1)
      stop('Increase max.steps.')
    if (max.norm <= sum(abs(b.tmpa)))
      stop('Increase max.norm.')
    if (lambda <= min.lambda)
      stop('Decrease min.lambda')
    if (max.vars <= 1)
      stop('Increase max.vars.')

    while(TRUE) {
      if (!backshoot) {
          
			arclength <- max.arclength
			if (new.A[k])
			 frac.arclength1 <- frac.arclength
			else {
			 frac.arclength1 <- 1
			 if (n.repeat1 > 1 && max.arclength == Inf)
			   arclength <- step$arclength
			}
			k <- k + 1
                
			cox.print.debug('\n#', k, ' (lambda: ', lambda, ')') 
		   gc();
			
			elapsed = system.time(step <- step.length.cox(corrector, depth, lambda2, exprs.all, status, rslist, wlist, min.lambda, arclength, frac.arclength1, add.newvars, backshoot, approx.Gram))
        
			if(elapsed[1] > max.time.per.step & k > min.steps) {
			   cat("Step took ", elapsed[1], "s (max.time.per.step: ", max.time.per.step, ").\nEnding algorithm.", sep='')
			   max.steps = k
			}

			b.tmpa[active] <- predictor.cox(b.tmpa[active], step)
			
			bmat.pred[[k]] <- b.tmpa[active]

			step.len[k - 1] <- h <- step$h
			lam.vec[k] <- lambda <- lambda + h
        
			tmpa <- append(active, step$newa)
			new.b = rep(0, length(step$newa))
			names(new.b) = step$newa
			b.tmpa = append(b.tmpa, new.b)
			b.tmpa = b.tmpa[tmpa]
        		  
			a <- abs(b.tmpa)
      }
      else {
        cox.print.debug('\n#', k, ' (backshoot):')
        
        step <- step.length.cox(corrector, depth, lambda2, exprs.all, status, rslist, wlist, min.lambda, Inf, 1, add.newvars, backshoot, approx.Gram, h)
                                
        step.len[k - 1] <- h + step$h
        h <- step$h
        lam.vec[k] <- lambda <- lambda + h
        a <- abs(b.tmpa)
      }
      		
      corrector <- corrector.cox(exprs.all, depth, status, rslist, wlist, rept, mthd, active, tmpa, force.active, lambda, lambda2, b.tmpa, bshoot.threshold, relax.lambda)
      newa <- corrector$newa		
		
      while(length(newa) > 0) {
        cox.print.debug('\nRepeating step ', k, ':')
        tmpa <- append(tmpa, newa)
        
	    new.b = rep(0, length(newa))
	    names(new.b) = newa
	    b.tmpa = append(b.tmpa, new.b)
		
        a <- abs(b.tmpa)
        
        corrector <- corrector.cox(exprs.all, depth, status, rslist, wlist, rept, mthd,
                                   active, tmpa, force.active, lambda, lambda2,
                                   b.tmpa, bshoot.threshold, relax.lambda)
        newa <- corrector$newa
        b.tmpa = corrector$b
      }
      newaction <- corrector$newactive
      if(length(corrector$inactive) > 0)
        newaction = append(newaction, paste("-", corrector$inactive, sep=""))
      
		if(class(active) == 'list') #DEBUG
			browser()
		
      if (length(corrector$active) <= n) {
          if (length(newaction) > 0) {
            if (corrector$backshoot && !backshoot) {
              cox.print.debug('\nOvershooting occurred: increasing lambda again')
              backshoot <- TRUE
              n.repeat1 <- 1
            } 
            else {
              active <- corrector$active
              b.tmpa <- corrector$b
              actions[[k]] <- newaction
              new.df[k] <- corrector$df
              new.A[k] <- TRUE
          
              if (trace) {
                 cat("added/dropped:\n")
                cat(paste("\t", newaction))
              }
          
              backshoot <- FALSE
              n.repeat1 <- n.repeat
            }
          }
          else {
            active <- corrector$active
            b.tmpa <- corrector$b
            backshoot <- FALSE
            n.repeat1 <- max(n.repeat1 - 1, 1)
          }
      }
      
      if (!backshoot) {
        bmat.corr[[k]] <- b.tmpa[active]
        cmat[k] <- corrector$corr
        lp[k] <- corrector$lp
        df[k] <- corrector$df
        
        if (lambda <= min.lambda 
            || k == max.steps 
            || length(corrector$active) > min(n, max.vars) 
            || sum(corrector$a) >= max.norm) {
            if (length(corrector$active) > min(n, max.vars))
                k <- k - 1

            if (lambda <= min.lambda)
                cox.print.debug('\nLambda = ', min.lambda, '\n')
            else if (k == max.steps) {
                cox.print.debug('\nMaximum steps (', max.steps, ') taken.\n')
				 }
            else if (length(corrector$active) > min(n, max.vars))
                cox.print.debug('\nNumber of active variables has reached its maximum.\n')
            else
                cox.print.debug('\n|beta| >= ', max.norm, '\n')

            break
        }
      }
    }
    
    unlink(paste(filename.prefix, "input.dat", sep=''))
    unlink(paste(filename.prefix, "output.dat", sep=''))
    unlink(paste(filename.prefix, "neg_weights.dat", sep=''))
    unlink(paste(filename.prefix, "weights.dat", sep=''))
        
    bmat.pred <- bmat.pred[1:k]
    bmat.corr <- bmat.corr[1:k]
    cmat <- cmat[1:k]

    df <- df[1:k]
    lp <- lp[1:k]
    aic <- -2 * lp + 2 * df
    bic <- -2 * lp + log(n) * df

    object <- list(call = call, lambda = lam.vec[1:k], lambda2 = lambda2,
                   step.length = abs(step.len[1:(k-1)]), corr = cmat,
                   new.df = new.df[1:k], df = df, loglik = lp, aic = aic,
                   bic = bic, b.predictor = bmat.pred, b.corrector = bmat.corr,
                   new.A = new.A[1:k], actions = actions[1:k],
                   sdx = rep(1, m), xnames = colnames(exprs.all), method = method,
                   nopenalty.subset = nopenalty.subset,
                   standardize = standardize, gene.symbols = gene.symbols, deviation.threshold = deviation.threshold, quantile = quantile)
    class(object) <- 'itemset.coxpath'
    object
  }

plot.itemset.coxpath <- function(x, xvar = c('norm', 'lambda', 'step', 'stepcoeffs'),
                         type = c('coefficients', 'aic', 'bic'),
                         plot.all.steps = FALSE, xlimit = NULL, plot.until.step = NULL,
                         predictor = FALSE, omit.zero = TRUE, breaks = TRUE,
                         mar = NULL, main = NULL, eps = .Machine$double.eps,
                         ...)
  {
    object <- x
    ii <- object$new.A
    if (plot.all.steps) {
      ii[!ii] <- TRUE
    } else {
      ii[length(ii)] <- TRUE
    }
    if(!is.null(plot.until.step))
      ii[plot.until.step:length(ii)] = FALSE
    lam <- object$lambda[ii]
    xvar <- match.arg(xvar)
    type <- match.arg(type)

    coef.pred <- object$b.predictor[ii]
    coef.corr <- object$b.corrector[ii]

    
    m <- ncol(coef.pred)
    k <- length(coef.corr)
    s <- switch(xvar, 
      norm = sapply(coef.corr, function(item) { if(length(item) == 0) 0 else sum(abs(item)) }),
      lambda = lam,
      step = seq(k),
      stepcoeffs = seq(k))
    if (xvar != 'lambda') {
      if (is.null(xlimit)) xlimit <- max(s)
      else if (xlimit <= min(s)) stop('Increase xlimit.')
      xi <- s <= xlimit
    } else {
      if (is.null(xlimit)) xlimit <- min(s)
      else if (xlimit >= max(s)) stop('Decrease xlimit.')
      xi <- s >= xlimit
    }
    

    coef.names = unique(names(unlist(coef.corr)))
    
    coef.corr = t(sapply(coef.corr, function(item) {
        new.item = rep(0, length(coef.names))
        names(new.item) = coef.names
        new.item[names(item)] = item 
        return(new.item)
    }));
    
    coef.pred = t(sapply(coef.pred, function(item) {
        new.item = rep(0, length(coef.names))
        names(new.item) = coef.names
        new.item[names(item)] = item 
        return(new.item)
    }))
    
    colnames(coef.corr) = colnames(coef.pred) = coef.names = line.to.gene.symbol(coef.names, object$gene.symbols)
    
    k <- max(which(xi))
    xname <- switch(xvar, norm = '|beta|', lambda = 'lambda', step = 'step', stepcoeffs = '')
    if (!is.null(mar)) 
        par(mar = mar)
    else {
        mar = c(2.8, 2.5, 2, 4)
       if(type == "coefficients")
          mar[4] = 10.5
       if(xvar == "stepcoeffs")
          mar[1] = 9
       par(mar = mar)
     }

     if(xvar == 'stepcoeffs')
        xaxt = 'n'
     else
        xaxt = NULL
     
    if (type == 'aic') {
      aic <- object$aic[ii][xi]
      plot(s[xi], aic, xlab = xname, ylab = 'AIC', type = 'o', pch = 16,
           cex = 0.3, xaxt = xaxt, xaxs = "i", mgp=c(1.5, 0.4, 0), ...)

    } else if (type == 'bic') {
      bic <- object$bic[ii][xi]
      plot(s[xi], bic, xlab = xname, ylab = 'BIC', type = 'o', pch = 16, cex = 0.3, xaxt = xaxt, xaxs = "i", mgp=c(1.5, 0.4, 0), ...)
      # if (is.null(main)) title('BIC', line = 2.5)
      # else title(main, line = 2.5)
    } else {
      ylab <- ifelse(object$standardize, 'coefficients', 'Coefficients')
      yvals = coef.corr[xi, ]
      # yvals[c(2:dim(yvals)[1], dim(yvals)[1]),] == 0 & 
      yvals[yvals == 0 & yvals[c(1, 1:dim(yvals)[1]-1),] == 0] = NA
      for(i in 1:dim(yvals)[2]) {
         if(i > length(object$actions[ii]))
            next #DEBUG: browser()
       
         yvals[i, names(object$actions[ii][[i]])][is.na(yvals[i, names(object$actions[ii][[i]])])] = 0
      }
      matplot(s[xi], yvals, xlab = xname, type = 'o', pch = '*',
              ylab = ylab, lty = 1, xaxt = xaxt, xaxs = "i", mgp=c(1.5, 0.4, 0), ...)

      abline(h = 0, lty = 3)
      
      at = coef.corr[k,]
      at = at[!is.na(yvals[k,])]
		
      coef.names = coef.names[!is.na(yvals[k,])]         
      
      axis(4, at = at, labels = FALSE, tick = TRUE)

		if(require(TeachingDemos)) {
			mindiff = 0.015 #1/(max(at) - min(at))
			# browser()
			at = spread.labs(at, mindiff=mindiff)
		}
      axis(4, at = at, labels = coef.names, cex = 0.5, adj = 0, las = 1, cex.axis=0.55, tick = FALSE)

      if (predictor) {
        for (i in 1:m) {
          segments(s[xi][-k], coef.corr[xi, ][-k,i], s[xi][-1],
                   coef.pred[xi, ][-1,i], lty = 2, col = i)
        }
      }
    }
    
    if (breaks) {
      new <- object$new.A[ii] & xi
      axis(3, at = s[new], labels = object$new.df[ii][new], cex = 0.6)
      abline(v = s[new], col = "grey", lwd=0.5)
    }
    
    if(xvar == 'stepcoeffs') {
       added.labs = NULL
       removed.labs = NULL
		 added.at = NULL
		 removed.at = NULL
		 
       for(i in seq(1, length(s[xi]))) {
          if(i > length(object$actions[ii]) || is.null(object$actions[ii][[i]]))
             next #DEBUG: browser()
            			
			to.add = line.to.gene.symbol(object$actions[ii][[i]][substr(object$actions[ii][[i]], 1, 1) != '-'], object$gene.symbols)			
			
			if((length(to.add) > 1) || to.add != "") {
	         added.labs = append(added.labs, to.add)
				added.at = append(added.at, rep(i, length(to.add)))
			}
			
			to.remove = line.to.gene.symbol(object$actions[ii][[i]][substr(object$actions[ii][[i]], 1, 1) == '-'], object$gene.symbols)
			if((length(to.remove) > 0) || to.remove != "") {
	         removed.labs = append(removed.labs, to.remove)
				removed.at = append(removed.at, rep(i, length(to.remove)))
			}         
      }
		
		all.at = spread.labs(c(added.at, removed.at), mindiff=0.3)
      at = all.at[1:length(added.at)]
      axis(1, at = at, labels = added.labs, cex = 0.5, adj = 1, las = 3, cex.axis=0.55, tick = FALSE, col.axis=1, hadj = 0.9)
      
      at = all.at[(length(added.at)+1):(length(added.at)+length(removed.at))]
      axis(1, at = at, labels = removed.labs, cex = 0.5, adj = 1, las = 3, cex.axis=0.55, tick = FALSE, col.axis=2, hadj = 0.9)
      
		
    }
    
  }

predict.itemset.coxpath <- function(object, data, s,
                            type = c('coefficients', 'loglik', 'lp', 'risk', 'coxph'),
                            mode = c('step', 'norm.fraction', 'norm', 'lambda.fraction', 'lambda'),
                            eps = .Machine$double.eps, exprs.all = NULL, deviation.threshold = NULL, quantile = NULL, ...) {
	 mode <- match.arg(mode)
	 type <- match.arg(type)
	 if (missing(data) && type != 'coefficients') {
	     warning('No data argument; type switched to coefficients')
	     type <- 'coefficients'
	 }
	 if (!missing(s)) {
	   if (length(s) > 1 && type == 'coxph') {
	     warning('Length(s) > 1. Only the first element is used.')
	     s <- s[1]
	   }
	 }
    
	 b <- object$b.corrector
	 if(is.null(b) || length(b) == 0) {
		 cat("No corrector parameters provided. Cannot fit.\n")
		 return(NULL)
	 }
	 
	 coef.names = unique(names(unlist(b)))

	 if(is.null(exprs.all))
		 exprs.all = extract.features(data$x, data$extra.features, deviation.threshold = deviation.threshold, quantile = quantile)
	
	 x.used = get.active.values(exprs.all, coef.names)
	 one <- rep(1, nrow(x.used))
	 meanx.used <- drop(one %*% x.used) / nrow(x.used)

	 b = t(sapply(b, function(item) {
	     new.item = rep(0, length(coef.names))
	     names(new.item) = coef.names
	     new.item[names(item)] = item 
	     return(new.item)
	 }));
    
	 k <- nrow(b)
	 steps <- seq(k)
	 if (missing(s)) {
	   s <- steps[object$new.A]
	   if (mode != 'step') {
	     warning('no s argument; mode switched to step')
	     mode <- 'step'
	   }
	 }
	 sb <- switch(mode, step = {
	   if (any(s < 1) || any(s > k)) 
	     stop('Argument s out of range')
	   steps
	 }, norm.fraction = {
	   if (any(s > 1) || any(s < 0)) 
	     stop('Argument s out of range')
	   bnorm <- apply(abs(b), 1, sum)
	   bnorm / bnorm[k]
	 }, norm = {
	   bnorm <- apply(abs(b), 1, sum)
	   if (any(s > bnorm[k]) || any(s < bnorm[1])) 
	     stop('Argument s out of range')
	   bnorm
	 }, lambda.fraction = {
	   if (any(s > 1) || any(s < 0))
	     step('Argument s out of range')
	   lam <- object$lambda
	   lam[lam < eps] <- eps
	   lam <- log(lam)
	   (lam - min(lam)) / (max(lam) - min(lam))
	 }, lambda = {
	   lam <- object$lambda
	   if (any(s > lam[1]) || any(s < lam[k]))
	     stop('Argument s out of range')
	   lam
	 })
	 sfrac <- (s - sb[1]) / (sb[k] - sb[1])
	 sb <- (sb - sb[1]) / (sb[k] - sb[1])
	 usb <- unique(sb)
	 useq <- match(usb, sb)
	 sb <- sb[useq]
	 b <- b[useq, ]
	 coord <- approx(sb, seq(sb), sfrac)$y
	 left <- floor(coord)
	 right <- ceiling(coord)
	 newb <- (((sb[right] - sfrac) * b[left, , drop = FALSE] + 
	           (sfrac - sb[left]) * b[right, , drop = FALSE]) /
	          (sb[right] - sb[left]))
	 newb[left == right, ] <- b[left[left == right], ]
	 coef <- newb
    
	 if (type == 'coefficients') {
	   fit <- coef
	   dimnames(fit) <- list(s, object$xnames)
	 } else if (type == 'loglik') {
	   fit <- logplik(x.used, data$time, data$status, t(coef), object$method)
	   names(fit) <- s
	 } else if (type == 'lp' || type == 'risk') {
	   b0 <- coef %*% meanx.used 
	   fit <- scale(x.used %*% t(coef), b0, FALSE)
	   if (type == 'risk') fit <- exp(fit)
	   dimnames(fit) <- list(seq(nrow(x.used)), s)      
	 } else {
	    print("ERROR: not implemented yet")
	    browser()
	   coef <- drop(coef)
	   active <- abs(coef) > eps
	   coef <- coef[active]
	   x <- x.used[, active, drop = FALSE]
	   time <- data$time
	   status <- data$status
	   fit <- coxph(Surv(time, status) ~ x, method = object$method)
	   junk <- logplik(x, time, status, coef, object$method, TRUE)
	   w <- junk$w
	   dmat <- junk$dmat
	   oo <- junk$oo
	   a <- sum(active)
	   info <- matrix(0, a, a)
	   for (i in 1:sum(status == 1)) {
	     ind <- dmat[, i] > 0
	     xr <- x[oo[ind], , drop = FALSE]
	     wr <- w[ind, i]
	     v1 <- xr * wr
	     v2 <- apply(v1, 2, sum)
	     info <- info + t(xr) %*% v1 - outer(v2, v2)
	   }
	   fit$coefficients <- coef
	   fit$var <- solve(info)
	   fit$loglik <- c(fit$loglik[1], junk$loglik)
	   fit$iter <- fit$residuals <- NULL
	   fit$linear.predictors <- junk$eta - sum(coef * meanx.used[active])
	   fit$method <- object$method
	   fit$assign <- seq(a)
	   fit$wald.test <- sum(coef*(info %*% coef))
	 }    
	 attr(fit, 's') <- s
	 attr(fit, 'fraction') <- sfrac
	 attr(fit, 'mode') <- mode
	 return(fit)
	}

	logplik <- function(x, time, status, b, method = c('breslow', 'efron'),
	                 return.all = FALSE)
	{
	 method <- match.arg(method)
	 n <- length(time)
	 o <- order(status, decreasing=T)
	 oo <- o[order(time[o])]
	 time <- time[oo]
	 status <- status[oo]
	 rept <- rep(0, n)
	 for (i in 1:n)
	     rept[i] <- sum(time[i:n] == time[i] & status[i:n] == 1)
    
	 complete <- which(status == 1)
	 nnc <- length(complete)
	 if (nnc == 0) {
	   stop('No complete observation. Failed to compute partial likelihood.')
	   browser()
	 }
	 dmat <- matrix(0, n, nnc)
	 for (i in 1:nnc) {
	   dmat[time >= time[complete[i]], i] <- 1
	   if (method == 'efron') {
	     if (rept[complete[i]] > 0) {
	       tie <- time == time[complete[i]] & status == 1
	       di <- max(rept[tie])
	       dmat[tie, i] <- dmat[tie, i] - (di - rept[complete[i]]) / di
	     }
	   }
	 }

	 eta <- x %*% b
	 eeta <- exp(eta)
	 k <- ncol(eta)
	 loglik <- rep(0, k)
	 for (i in 1:k) {
	   w <- dmat * eeta[oo, i]
	   wsum <- apply(w, 2, sum)
	   loglik[i] <- sum(eta[oo, i][status == 1]) - sum(log(wsum))
	 }
	 if (return.all) {
	   return(list(loglik = loglik, w = scale(w, F, wsum), eta = eta,
	               dmat = dmat, oo = oo))
	 } else {
	   return(loglik)
	 }
	}


print.itemset.coxpath <- function(x, ...)
  {
    cat('Call:\n')
    dput(x$call)
    actions <- line.to.gene.symbol(x$actions, x$gene.symbols)
    
    k <- length(actions)
    for (i in 1:k) {
      if (length(actions[[i]]) > 0) {
        cat('Step', i, ':')
        # for (ii in actions[[i]]) {
            cat(paste("\t", actions[[i]]))
            # cat("\n")
            # cat(paste("\t", x$actions[[i]]))
        cat('\n')
      }
    }
  }

summary.itemset.coxpath <- function(object, ...)
  {
    cat('Call:\n')
    dput(object$call)
    ii <- object$new.A
    ii[length(ii)] <- TRUE
    M <- data.frame(Df = object$df[ii], Log.p.lik = object$loglik[ii],
                    AIC = object$aic[ii], BIC = object$bic[ii])
    dimnames(M)[[1]] <- paste('Step', which(ii), sep=' ')
    M
  }

line.to.gene.symbol <- function(object, gene.symbols, reorder=F) {
   if(class(object) == "list")
      return(lapply(object, line.to.gene.symbol, gene.symbols))
   if(length(object) > 1)
      return(sapply(object, USE.NAMES=FALSE, line.to.gene.symbol, gene.symbols))
   if(is.null(object))
      return(NULL)
   if(length(object) == 0)
           return('')
   
   object = gsub("down.", "dn.", object, fixed=T)
	
	group = sapply(strsplit(as.character(object), "[\\+\\*]", fixed = F)[[1]], function(item) {
      items = strsplit(item, ".", fixed=TRUE)[[1]]
      
      if(length(items) > 1) {
         name = gsub(" +", "", gene.symbols$Gene.Symbol[as.numeric(items[2])])
         name = gsub("\\/\\/\\/.*", "#", name)

         if (is.na(name) || is.null(name) || name == '') {
            name = paste('pb', gene.symbols$ID[as.numeric(items[2])], sep='')
         }
         return(paste(items[1], sep=".", name))
      }
      else
         return(item)
   })
	
	if(reorder)
		group = group[order(grepl(".", group, fixed=T))]
	
	paste(group, collapse="*")   
}



cox.print.debug <- function(...) {
   if(global.trace)
    cat(paste(..., "\n", sep=""))
}

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

run.cv <- function(data, gene.symbols = NULL, prefix = "", load.saved.file = T, training.set.ratio = 0.6, num.repeats = 20, trace = 3, min.lambda = 0.5, depth = 3, ...) {
	if(is.null(gene.symbols) && length(data) == 2 && class(data) == 'list' && ncol(data[[1]]$x) == nrow(data[[2]])) {
		gene.symbols = data[[2]]
		data = data[[1]]
	}
	
	n = length(data$time)
	fits = list()
	run.times = list()
	start.i = 1
	depth.str = paste('.d', depth, sep = '')
	
   file.name = paste(prefix, "cv.", as.character(sys.call()[2]), ".subset-", training.set.ratio, depth.str, '.data', sep = '')
	
	print(file.name)
	
	if(load.saved.file && file.exists(file.name)) {
		load(file.name)
		start.i = length(fits)+1
		cat("###\nSkipping: ", length(fits), " repeats\n\n")
	}
	
	if(start.i <= num.repeats) {	
		for(i in start.i:num.repeats) {
			cat('CV repeat: ', i, '/', num.repeats, "\n", sep = '')
		
			shuffled = sample(seq(n))
				   training = (shuffled <= n*training.set.ratio)
			
			fit = itemset.coxpath(data[training,], gene.symbols, depth = depth, trace = trace, min.lambda = min.lambda, ...)

			fits[[i]] = list(fit = fit, training.set = training)
			save(fits, file=file.name)
		   gc();
		}
	}
	
	
	return(fits)
}


plot.all.kaplan.meier <- function(vars, data, fit, max.fit, group.ratios, title = "KM") {
   par(mfrow = c(4, 3), oma = c(0, 0, 3, 0))
   
   for(i in 1:prod(par('mfrow'))) {
      
      if(i > length(vars)) {
			plot.new()
         break
      }
      
      plot.kaplan.meier(names(vars)[i], data, gene.symbols, group.ratios, fit$b.corrector[[max.fit]][names(vars)[i]])
   }

   mtext(paste(title, " - Step ", max.fit , "/", length(fit$aic), sep=""), outer = TRUE, cex=1.2)
   
}

gene.symbol.to.line <- function(object, gene.symbols) {
   if(class(object) == "list")
      return(lapply(object, gene.symbol.to.line, gene.symbols))
   if(length(object) > 1)
      return(sapply(object, USE.NAMES=FALSE, gene.symbol.to.line, gene.symbols))
   if(is.null(object))
      return(NULL)
   if(length(object) == 0)
           return('')
   
   filtered.symbols = gsub("\\/\\/\\/.*", "#", gsub(" +", "", gene.symbols$Gene.Symbol))
   object = gsub("\\s*", "", object) # remove whitespace
   object = gsub("dn.", "dn.", object, fixed=T)
	
   paste(sapply(strsplit(as.character(object), "[\\+\\*]", fixed = F)[[1]], function(item) {
      items = strsplit(item, ".", fixed=TRUE)[[1]]
      
      if(length(items) > 1) {
         lines = which(filtered.symbols == paste(items[-1], collapse='.'))
         if(length(lines) == 0) {
            lines = which(gene.symbols$ID == paste(items[-1], collapse='.'))
            if(length(lines) == 0)
               line = items[2]
            else
               line = lines[1]
         }
         else 
            line = lines[1]
            
         return(paste(items[1], sep=".", line))
      }
      else
         return(gsub('−', '-', item)) # replacing unicode by ascii
   }), collapse="*")   
   
}


plot.interaction.kaplan.meier.from.gene.symbols <- function(combi.genes, data.1, gene.symbols, data.2 = NULL, data.3 = NULL, data.1.title = NULL, data.2.title = NULL, data.3.title = NULL, ...) {
   return(plot.interaction.kaplan.meier(gene.symbol.to.line(combi.genes, gene.symbols), data.1 = data.1, gene.symbols, data.1.title = ifelse(is.null(data.1.title), as.character(sys.call()[3]), data.1.title),
    data.2 = data.2, data.2.title = ifelse(is.null(data.2.title), as.character(sys.call()[4]), data.2.title),
    data.3 = data.3, data.3.title = ifelse(is.null(data.3.title), as.character(sys.call()[5]), data.3.title), ...))
}

plot.interaction.kaplan.meier <- function(combi.var, data.1, gene.symbols, data.1.title = "data.1", data.2 = NULL, data.2.title = "data.2", data.3 = NULL, data.3.title = "data.3", var.coef = NULL, group.ratios = c(0.5, 0.5), main.title = NULL, pval.threshold = 0.05) {
	if(length(group.ratios) == 0 || group.ratios[1] == 0)
		stop("Error: group.ratios[1] == 0");
	
	if(length(combi.var) > 1)
		combi.vars = combi.var
	else
		combi.vars = c(combi.var)
	
	tot.per.line = max(sapply(combi.vars, function(x) { length(strsplit(x, '*', fixed=T)[[1]]) })) + 1
	# tot.per.line = max(sapply(combi.var, function(x) {length(strsplit(x, '*', fixed=T)[[1]])})) + 1
	tot.per.col = max(3, min(length(combi.var), 6))
	
   par(mfrow = c(tot.per.col, tot.per.line), oma = c(0, 0, 3, 0))
   
	for(i in 1:length(combi.vars)) {
		combi.var = combi.vars[i]
		
		if(substr(combi.var, 1, 1) == '¬') {
			combi.var.clean = substr(combi.var, 3, nchar(combi.var)-2)
		}
		else {
			combi.var.clean = combi.var
		}
		
	   pval.combi = plot.kaplan.meier(combi.var, data.1, gene.symbols, group.ratios, var.coef)
   
	   if(!is.null(data.2)) {
	      plot.kaplan.meier(combi.var, data.2, gene.symbols, group.ratios, var.coef)
	      mtext(data.2.title, adj = 3/(2*tot.per.line), outer = TRUE, line = -.8, cex = 0.8)
	   }
	   if(!is.null(data.3)) {
	      plot.kaplan.meier(combi.var, data.3, gene.symbols, group.ratios, var.coef)
	      mtext(data.3.title, adj = 5/(2*tot.per.line), outer = TRUE, line = -.8, cex = 0.8)
	   }
   
		split.genes = strsplit(combi.var.clean, '*', fixed=T)[[1]]
		
	   pvals = sapply(split.genes, function(gene) {
	      pval = plot.kaplan.meier(gene, data.1, gene.symbols, group.ratios)
      
	      if(!is.null(data.2))
	         plot.kaplan.meier(gene, data.2, gene.symbols, group.ratios)

	      if(!is.null(data.3))
	         plot.kaplan.meier(gene, data.3, gene.symbols, group.ratios)
      
	      return(pval)
	   })
		
		blanks = tot.per.line - length(split.genes) - 1
		if(blanks > 0)
			for(i in 1:blanks)
				plot.new()
			
	}
	
	if(length(combi.vars) == 1) {
	   title = line.to.gene.symbol(combi.var, gene.symbols)
	   if(!is.null(main.title))
	      title = paste(title, " (", main.title, ")", sep='')
	   mtext(title, outer = TRUE, cex=1.2, line=1)
	   mtext(data.1.title, adj = 1/(2*tot.per.line), outer = TRUE, line = -.8, cex = 0.8)
	}
	else {
	   mtext(data.1.title, outer = TRUE, cex=1.2, line=1)
	}
	
   return((pval.combi < pval.threshold) & any(pvals > pval.threshold))
}

plot.synthetic.kaplan.meier.from.gene.symbols <- function(combi.genes, data.1, gene.symbols, ...) {
   return(plot.synthetic.kaplan.meier(gene.symbol.to.line(combi.genes, gene.symbols), data.1 = data.1, gene.symbols, ...))
}

plot.synthetic.kaplan.meier <- function(combi.var, data.1, gene.symbols, data.1.title = "data.1") {

	group.ratios = c(0.5, 0.5)
	
	if(length(combi.var) == 1)
		combi.vars = c(combi.var)
	else
		combi.vars = combi.var
	
	tot.per.line = 4
	tot.per.col = max(3, min(length(combi.vars), 6))
	
   par(mfrow = c(tot.per.col, tot.per.line), oma = c(0, 0, 3, 0))
   
	for(i in 1:length(combi.vars)) {
		combi.var = combi.vars[i]

	   pval.combi = plot.kaplan.meier(combi.var, data.1, gene.symbols, group.ratios, plot.only.colours = c("red", "blue"), title="")
	   pval.combi = plot.kaplan.meier(combi.var, data.1, gene.symbols, group.ratios, plot.only.colours = c("green", "purple"), title="")

	   pval.combi = plot.kaplan.meier(combi.var, data.1, gene.symbols, group.ratios, plot.only.colours = c("red", "purple"), title="")
	   pval.combi = plot.kaplan.meier(combi.var, data.1, gene.symbols, group.ratios, plot.only.colours = c("green", "blue"), title="")
		
	}
	
   mtext(data.1.title, outer = TRUE, cex=1.2, line=1)
}

plot.kaplan.meier.from.gene.symbols <- function(var, data, gene.symbols, ...) {
	return(plot.kaplan.meier(gene.symbol.to.line(var, gene.symbols), data, gene.symbols, ...))
}

plot.kaplan.meier <- function(var, data, gene.symbols, group.ratios = c(0.5, 0.5), var.coef = NULL, plot.third.group = F, plot.only.colours = NULL, title = NULL, legend.size = 0.8) {
	if(length(group.ratios) == 0 || group.ratios[1] == 0)
		stop("Error: group.ratios[1] == 0");
	if(!is.null(plot.only.colours) && (group.ratios[[1]] != 0.5 || group.ratios[[2]] != 0.5)) {
		print("For four-group plots, group threshold needs to be 0.5")
		group.ratios = c(0.5, 0.5)
	}
	
	if(substr(var, 1, 1) == '¬') {
		neg.correl = T
		var = substr(var, 3, nchar(var)-2)
	}
	else
		neg.correl = F
	
   groups = sapply(strsplit(var, '*', fixed=T)[[1]], function(gene) {
      gene = strsplit(gene, ".", fixed=TRUE)[[1]]
      			
		if(length(gene) == 1)
         col = data$extra.features[, gene] * 2
      else {
         if(gene[1] == 'up')
            col = data$x[, as.integer(gene[2])]
         else
            col = -data$x[, as.integer(gene[2])]
         
         col = (col > col[order(col)][floor(length(data$time)*group.ratios[1])]) + (col > col[order(col)][ceiling(length(data$time)*group.ratios[2])])
      }
      return(col)
   })
		
	if(!is.null(plot.only.colours)) {	
   	groups = groups[,1] + groups[,2]/2
      colours = c("red", "blue", "purple", "green")
   	
		select.colours = colours %in% plot.only.colours
		
		colours = colours[select.colours]
		select.groups = groups %in% (0:3)[select.colours]
		
		if(length(table(groups[select.groups])) < 2) {
			plot.new()
			return(0)
		}
		
	   surv <- survfit(Surv(data$time[select.groups], data$status[select.groups] == 1)~groups[select.groups]);
	   logrank <- survdiff(Surv(data$time[select.groups], data$status[select.groups] == 1)~groups[select.groups]);
		
		genes = strsplit(var, '*', fixed=T)[[1]]
		legend = c(
			paste(opposite.expression(line.to.gene.symbol(genes[1], gene.symbols)), "AND", opposite.expression(line.to.gene.symbol(genes[2], gene.symbols))), 
			paste(opposite.expression(line.to.gene.symbol(genes[1], gene.symbols)), "AND", line.to.gene.symbol(genes[2], gene.symbols)),
			paste(line.to.gene.symbol(genes[1], gene.symbols), "AND", opposite.expression(line.to.gene.symbol(genes[2], gene.symbols))),
			paste(line.to.gene.symbol(genes[1], gene.symbols), "AND", line.to.gene.symbol(genes[2], gene.symbols))) 
		legend = legend[select.colours]
	   p.val = pchisq(logrank$chisq, 1, lower.tail = FALSE);
	}
	else {
		   groups = apply(groups, 1, sum)/ncol(groups)
			groups = (groups >= 1) + (groups >= 2)
		# }
				
	   surv <- survfit(Surv(data$time[groups != 1], data$status[groups != 1] == 1)~groups[groups != 1]);
	   logrank <- survdiff(Surv(data$time[groups != 1], data$status[groups != 1] == 1)~groups[groups != 1]);
	   p.val = pchisq(logrank$chisq, 1, lower.tail = FALSE);
		
		 if(plot.third.group) {
	      colours = c("blue", "grey", "red")
		   surv <- survfit(Surv(data$time, data$status == 1)~groups);
		   logrank <- survdiff(Surv(data$time, data$status == 1)~groups);
			legend = c(gsub("*", " AND ", opposite.expression(line.to.gene.symbol(var, gene.symbols)), "Everything else", line.to.gene.symbol(var, gene.symbols), fixed=T))
		}
	   else {
	      colours = c("blue", "red")

			legend = gsub("*", " AND ", c(opposite.expression(line.to.gene.symbol(var, gene.symbols)), line.to.gene.symbol(var, gene.symbols)), fixed=T)
		}
	}
	
	par(mar=c(2.5, 2.5, 2.8, 1))
	if(neg.correl)
		prefix = '¬ '
	else
		prefix = ''
	
	if(is.null(title)) {
		if(! is.null(plot.only.colours))
			title = paste(prefix, gsub('up.', '', gsub('dn.', '', line.to.gene.symbol(var, gene.symbols), fixed=T), fixed=T), ifelse(is.null(var.coef), "", paste(" (", round(var.coef, digits=3), ")", sep='')), sep="")
		else
			title = paste(prefix, line.to.gene.symbol(var, gene.symbols), ifelse(is.null(var.coef), "", paste(" (", round(var.coef, digits=3), ")", sep='')), sep="")
	}
   plot(surv, conf.int = F, col = colours, main= title, cex.main=0.9);

	if(p.val < 1e-5)
		mtext("p-val < 1e-5", cex = 0.65)
	else
		mtext(paste("p-val:", round(p.val, digits = 6)), cex = 0.65)
         	
	legend(x = "bottomleft", legend = legend, fill = colours, cex=legend.size)
	
	
   return(p.val)
}

get.kaplan.meier.p.value <- function(var, data, group.ratios = c(0.5, 0.5)) {
	if(length(group.ratios) == 0 || group.ratios[1] == 0)
		stop("Error: group.ratios[1] == 0");
	
	if(substr(var, 1, 1) == '¬') {
		neg.correl = T
		var = substr(var, 3, nchar(var)-2)
	}
	else
		neg.correl = F
	
   groups = sapply(strsplit(var, '*', fixed=T)[[1]], function(gene) {
      gene = strsplit(gene, ".", fixed=TRUE)[[1]]
      			
		if(length(gene) == 1)
         col = data$extra.features[, gene] * 2
      else {
         if(gene[1] == 'up')
            col = data$x[, as.integer(gene[2])]
         else
            col = -data$x[, as.integer(gene[2])]
         
         col = (col > col[order(col)][floor(length(data$time)*group.ratios[1])]) + (col > col[order(col)][ceiling(length(data$time)*group.ratios[2])])
      }
      return(col)
   })
	
   groups = apply(groups, 1, sum)/ncol(groups)
	groups = (groups >= 1) + (groups >= 2)

   logrank <- survdiff(Surv(data$time[groups != 1], data$status[groups != 1] == 1)~groups[groups != 1]);
   p.val = pchisq(logrank$chisq, 1, lower.tail = FALSE);
			
   return(p.val)
}


opposite.expression <- function(var) {
   opposites = sapply(strsplit(var, '*', fixed=T)[[1]], function(gene) {
      gene = strsplit(gene, ".", fixed=TRUE)[[1]]
      			
		if(length(gene) == 1)
         val = paste("¬", gene)
      else {
         if(gene[1] == 'up')
            val = paste('dn.', gene[2], sep='')
         else
            val = paste('up.', gene[2], sep='')
    	}
		 return(val)		
	})
	return(paste(opposites, collapse=" AND "))
}

