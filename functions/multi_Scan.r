## load the previous SaRa program code
#dyn.load("diagnosticValue64.dll")
#source("llce.r")

## Multi-SaRa starts here
multi.llce <- function(x, h, sd.est)
{
    if(is.vector(x)) x <- matrix(x, nrow = 1)
    N <- nrow(x)
    T <- ncol(x)
    diag.stat <- matrix(NA, nrow = N, ncol = T)
    for(i in 1:N) {
        diag.stat[i, ] <- llce(x[i, ], h)
    }
    return(diag.stat/sd.est * sqrt(h/2))
}

## Test multi.llce
##xx <- matrix(rnorm(10000 * 2000), 2000, 10000)
##system.time(stat.xx <- multi.llce(xx, 10))
##   user  system elapsed
##  4.062   0.398   4.459

zmat.pval <- function(x, sd.est = 1)
{
    if(is.vector(x)) x <- matrix(x, nrow = 1)
    N <- nrow(x)
    T <- ncol(x)
    p <- 2 * pnorm(abs(x), mean = 0, sd = sd.est, lower.tail = FALSE)
}

## Test llce.pval
##pval.xx <- llce.pval(stat.xx, 1)

##estimateSigma <-
##function(Y,h=10){                   #constant case
##  n     = length(Y)                                 # can we make it faster?
##  YBar  = rep(0,n)
##  for (i in 1:n) {
##     a       = min(n,i+h)
##     b       = max(1,i-h)
##     YBar[i] = mean(Y[b:a])
##  }
##  return(sd(Y-YBar))
##}

estimateSigma <- function(x, h = 10) {
  T <- length(x)
  sums <- cumsum(c(rep(0, h + 1), x, rep(0, h)))
  sum.loc <- sums[(1:T) + 2 * h + 1] - sums[1:T]
  ns <- c(h + (1:h), rep(2 * h + 1, T - (2 * h)), h + (h:1))
  xbar <- sum.loc/ns
  res <- x - xbar
  var0 <- var(res) * (2 * h + 1)/(2 * h)
  return(sqrt(var0))
}



##estimate.sd <- function(x, h = 20) {
##  sds <- apply(x, 1, estimateSigma, h = h)
##  return(sds)
##}

estimate.sd <- function(x) {
  sds <- apply(x, 1, sd)
  return(sds)
}


## Multiple U kernel starts here
multi.scanU <- function(x, l, sd.est) {
  T <- ncol(x)
  N <- nrow(x)
  res <- matrix(NA, N, T - l + 1)
  for(i in 1:N) {
    sums <- cumsum(c(0, x[i, ]))
    res[i, ] <- sums[(1 + l):(T + 1)] - sums[1:(T - l + 1)]
  }
  means <- rowMeans(x)
  res <- (res - means * l)/sd.est/sqrt(l * (1 - l/T))
  return(res)
}

multi.cumsum <- function(x){
  T <- ncol(x)
  N <- nrow(x)
  sums1 <- matrix(0, nrow = N, ncol = T)
  for(i in 1:N) {
    sums1[i, ] <- cumsum(x[i, ])
  }
  sums <- cbind(0, sums1)
  return(sums)
}

multi.scanU.2 <- function(sums, l, sd.est, means = NULL) {
  T <- ncol(sums) - 1
  N <- nrow(sums)
  local.sums <- sums[, (l + 1):(T + 1)] - sums[, 1:(T - l + 1)]
  if (is.null(means)) means <- sums[, T + 1]/T
  res <- (local.sums - means * l)/sd.est/sqrt(l * (1 - l/T))
  return(res)
}


## Higher Criticism ##
higher.criticism <- function(p,alpha0)
{
  N <- nrow(p)
  T <- ncol(p)
  HC <- rep(0, T)
  M <- rep(0, T)
  for (i in 1:T) {
    pp <- sort(p[, i])
    W <- ((1:N)/N - pp)/sqrt(pp * (1 - pp)) * sqrt(N)
    HC[i] <- max(W[alpha0:floor(N/2)])
    M[i] <- which.max(W[alpha0:floor(N/2)]) + alpha0 - 1
  }
  attr(HC, "max") <- M
  return(HC)
}

## Adaptive Fisher
af.weights <- function(k, rr, ss) {
    r <- k - ss + 1
    s <- k - rr + 1
    if(r >= 1) {
        w1 <- (s - r + 1)/(k - 1:r + 1)
    } else w1 <- NULL
    if (r + 1 < s) {
        w2 <- (s - (r + 1):s + 1)/(k - (r + 1):s + 1)
    } else w2 <- NULL
    return(c(w1, w2))
}

adaptive.fisher <- function(p, alpha0) {
    N <- nrow(p)
    T <- ncol(p)
    means <- rep(NA, N)
    sds <- rep(NA, N)
    for(i in 1:N) {
        ww <- af.weights(N, 1, i)
        means[i] <- 2 * sum(ww)
        sds[i] <- 2 * sqrt(sum(ww^2))
    }
    #pp <- apply(p, 2, sort)
    #tt <- -2 * apply(log(pp), 2, cumsum)
    #zz <- (tt - means)/sds
    zz <- (-2 * apply(log(apply(p, 2, sort)), 2, cumsum) - means)/sds
    AF <- rep(NA, T)
    attr(AF, "max") <- rep(NA, T)
    for (j in 1:T){
        w <- which.max(zz[alpha0:floor(N/2), j]) + alpha0 - 1
        attr(AF, "max")[j] <- w
        AF[j] <- zz[w, j]
    }
    #AF <- apply(zz[alpha0:floor(N/2), ], 2, max)
    #attr(AF, "max") <- apply(zz[alpha0:floor(N/2), ], 2, which.max) + alpha0 - 1
    return(AF)
    #return(apply(zz[alpha0:floor(N/2), ], 2, max))
}

## Fisher
fisher <- function(p) colSums(-log(p))

## Stouffer
stouffer <- function(p) colSums(-qnorm(p))

## Weighted sum of squared statistics ##
weighted.sum <- function (U, p0 = 1){
  U.dim1 <- dim(U)[1]
  U.dim2 <- dim(U)[2]
  U.squared <- U^2
  weights <- exp(U.squared/2)/((1 - p0)/p0 + exp(U.squared/2))
  res <- colSums(weights * U.squared/2)
  return(res)
}

## Max of squared statistics
max.square <- function(U) apply(U^2, 2, max)

## My function to estimate the trend in data ##
mytrend <- function(x, h = 100) {
  T <- length(x)
  sums <- cumsum(c(rep(0, h + 1), x, rep(0, h)))
  sum.loc <- sums[(1:T) + 2 * h + 1] - sums[1:T]
  ns <- c(h + (1:h), rep(2 * h + 1, T - (2 * h)), h + (h:1))
  xbar <- sum.loc/ns
  return(xbar)
}

detrend <- function(x, ...) {
    N <- nrow(x)
    T <- ncol(x)
    x.detrend <- matrix(0, N, T)
    for(i in 1:N) {
        x.detrend[i, ] <- x[i, ] - mytrend(x[i, ], ...)
    }
    return(x.detrend)
}


## function for thresholding the multi-sample U scan statistic
thres.U <- function(x, th) {
    wmax <- length(x)
    for(i in 1:wmax) {
        names(x[[i]]) <- 1:length(x[[i]])
        x[[i]] <- na.omit(x[[i]])
        x[[i]] <- x[[i]][x[[i]] > th]
    }
    len <- sapply(x, length)
    pool <- cbind(rep(1:wmax, len), as.numeric(unlist(sapply(x, names))), unlist(x))
    rownames(pool) <- NULL
    colnames(pool) <- c("win", "pos", "score")
    pool <- pool[order(pool[, 3], decreasing = TRUE), ]
    end <- pool[, "pos"] + pool[, "win"] - 1
    pool <- cbind(pool[, 1:2], end, pool[, 3])
    colnames(pool) <- c("win", "pos", "end", "score")
    ch.regions <- list()
    while(nrow(pool) > 0) {
        ch.regions[[length(ch.regions) + 1]] <- pool[1, ]
        start <- pool[1, 2]
        end <- pool[1, 3]
        cond1 <- pool[, 2] <= end
        cond2 <- start <= pool[, 3]
        cond <- cond1 & cond2
        pool <- pool[-which(cond), , drop = FALSE]
        cat(nrow(pool), "regions remaining in the pool!\n")
    }
    res <- do.call(rbind, ch.regions)
    colnames(res) <- c("win", "pos", "end", "score")
    return(res)
}

## function for finding the local maximum of the D statistics given the window sizes, than threshold by a given value
localmax.D <- function(x, h, th) {
    n.win <- length(h)
    T <- ncol(x)
    if(n.win != nrow(x)) stop("Error: ncol(x) != length(h)\n")
    if(n.win != length(th)) stop("Error: length(th) != length(h)\n")
    cp <- list()
    for(i in 1:n.win) {
        hh <- h[i]
        is.max <- rep(FALSE, T - 2 * hh)
        for(j in (hh + 1):(T - hh)) {
            if((hh + 1) == which.max(x[i, (j - hh):(j + hh)]) && x[i, j] > th[i]) is.max[j - hh] <- TRUE
        }
        cp[[i]] <- which(is.max) + hh
    }
    return(cp)
}

## calculate the D matrix for the given change point
calc.D.chp <- function(x, chp, h) {
    n.win <- length(h)
    if(length(chp) != n.win) stop("Error: length(chp) != n.win\n")
    N <- nrow(x)
    sds <- estimate.sd(x)
    D <- list()
    for (i in 1:n.win){
      hh <- h[i]
      n.chp <- length(chp[[i]])
      D.chp <- matrix(NA, N, n.chp)
      xx <- cbind(matrix(0, ncol = hh - 1, nrow = N), x, matrix(0, ncol = hh, nrow = N))
      for (j in 1:n.chp) {
        D.chp[, j] <- rowMeans(xx[, chp[[i]][j] + c(0:(hh - 1))]) - rowMeans(xx[, chp[[i]][j] + c(hh:(2 * hh - 1))])
      }
      D[[i]] <- D.chp/sds * sqrt(hh/2)
    }
    return(D)
}

## calculate the U matrix for the given regions
calc.U.region <- function(x, region) {
    N <- nrow(x)
    T <- ncol(x)
    R <- nrow(region)
    Umat <- matrix(NA, N, R)
    means <- rowMeans(x)
    sds <- apply(x, 1, sd)
    for(i in 1:R) {
        start <- region[i, "pos"]
        end <- region[i, "end"]
        l <- region[i, "win"]
        if(start == end) {
            U <- x[, start]
        } else {
            U <- rowSums(x[, start:end])
        }
        Umat[, i] <- (U - l * means)/(l * (1 - l/T))/sds
    }
    return(Umat)
}

## calculate the BIC
get.one.BIC <- function(x, cp, mod = FALSE) {
  T <- length(x)
  cp <- setdiff(cp, T)
  J <- length(cp)
  
  start <- c(1, cp + 1)
  end <- c(cp, T)
  len <- end - start + 1
  means <- rep(NA, J + 1)
  for (i in 1:(J + 1)) {
    means[i] <- mean(x[start[i]:end[i]])
  }
  SST <- sum(x^2)
  sigma2 <- (SST - sum(len * means^2))/T

  if (mod) {
    bic <- T/2 * log(sigma2) + 3/2 * J * log(T) + 1/2 * sum(log(len/T))
  } else {
    bic <- T/2 * log(sigma2) + J * log(T)
  }
  return(list(bic = bic, means = means, len = len, SST = SST, cp = cp, mod = mod))
}

get.all.BICs <- function(x, cp, mod = FALSE) {
    N <- nrow(x)
    if (is.list(cp) && length(cp) != N) {
        stop("length(cp) != nrow(x)\n")
    }
    bic <- list()
    if(is.list(cp)) {
        for(i in 1:N) {
            bic[[i]] <- get.one.BIC(x[i, ], cp[[i]], mod)
        }
    } else {
        for(i in 1:N) {
            bic[[i]] <- get.one.BIC(x[i, ], cp, mod)
        }
    }
    return(bic)
}

remove.one.Change.Point <- function(BIC) {
  T <- sum(BIC$len)
  J <- length(BIC$len) - 1
  if(J < 1) stop("No change point left!\n")
  if (BIC$mod) {
    delta.term3 <- rep(NA, J)
    for (i in 1:J) {
      delta.term3[i] <- 1/2 * (log((BIC$len[i] + BIC$len[i + 1])/T) - sum(log(BIC$len[i:(i + 1)]/T)))
    }
    delta0 <- -3/2 * log(T) + delta.term3[i]
  } else {
    delta0 <- -log(T)
  }
  s2 <- (BIC$SST - sum(BIC$len * BIC$means^2))/T
  delta1 <- rep(NA, J)
  for(i in 1:J) {
    means.temp <- BIC$means
    means.temp[i] <- means.temp[i + 1] <- (BIC$means[i] * BIC$len[i] + BIC$means[i + 1] * BIC$len[i + 1])/sum(BIC$len[i:(i + 1)])
    s2.new <- (BIC$SST - sum(BIC$len * means.temp^2))/T
    delta1[i] <- T/2 * (log(s2.new) - log(s2))
  }
  delta <- delta0 + delta1
  bic <- BIC$bic + delta
  bic.new <- min(bic)
  cp.remove <- which.min(bic)
  len.new <- BIC$len[-cp.remove]
  len.new[cp.remove] <- BIC$len[cp.remove] + BIC$len[cp.remove + 1]
  means.new <- BIC$means[-cp.remove]
  means.new[cp.remove] <- (BIC$means[cp.remove] * BIC$len[cp.remove] + BIC$means[cp.remove + 1] * BIC$len[cp.remove + 1])/sum(BIC$len[cp.remove:(cp.remove + 1)])
  BIC.new <- list(bic = bic.new, means = means.new, len = len.new, SST = BIC$SST, mod = BIC$mod, cp = BIC$cp[-cp.remove])
  return(list(BIC.new = BIC.new, cp.remove = cp.remove, delta = delta))
}

iterative.remove.cp <- function(BIC) {
  try <- TRUE
  while (try && length(BIC$cp) > 0) {
    BIC.next <- remove.one.Change.Point(BIC)
    if (min(BIC.next$delta) < 0) {
      BIC <- BIC.next$BIC
    } else {
      try <- FALSE
    }
  }
  return(BIC)
}

soft.thresholding.mu <- function(x, cp, delta) {
    J <- length(cp)
    T <- ncol(x)
    N <- nrow(x)
    start <- c(1, cp + 1)
    end <- c(cp, T)
    means <- matrix(NA, N, J + 1)
    for (i in 1:(J + 1)) {
        means[, i] <- rowMeans(x[, start[i]:end[i], drop = FALSE])
    }
    mu <- sign(means) * pmax(abs(means) - delta, 0)
    keep <- (mu[, 1:J] != 0) | (mu[, 2:(J + 1)] != 0)
    cp.keep <- list()
    for (i in 1:N) {
        cp.keep[[i]] <- cp[keep[i, ]]
    }
    return(cp.keep)
}

lasso.delta <- function(x, cp, lambda, prec = 1e-7) {
    T <- length(x)
    start <- c(1, cp + 1)
    end <- c(cp, T)
    len <- end - start + 1
    J <- length(cp)
    means <- rep(NA, J + 1)
    for(i in 1:(J + 1)) {
        means[i] <- mean(x[start[i]:end[i]])
    }
    mu <- means
    del <- mu - c(0, mu[1:J])
    TSS <- sum(x^2)
    target <- function(TSS, del) {
        mu <- cumsum(del)
        ##print(TSS)
        ##print(del)
        sum((x - rep(mu, len))^2)/T/2 + lambda * sum(abs(del))
    }
    tar <- target(TSS, del)
    print(tar)
    converge <- FALSE
    while (!converge) {
        for (i in 1:(J + 1)) {
            temp <- del
            temp[i] <- 0
            mu <- cumsum(temp)
            dd <- sum(len[i:(J + 1)] * (means - mu)[i:(J + 1)])/sum(len[i:(J + 1)])
            if (i == 1) {
                del[i] <- dd
            } else {
                del[i] <- sign(dd) * max(abs(dd) - lambda/sum(len[i:(J + 1)]) * T, 0) 
            }
        }
        tar.old <- tar
        tar <- target(TSS, del)
        print(tar)
        if ((tar.old - tar)/tar < prec) converge <- TRUE
    }
    return(list(del = del, cp = cp[del[2:(J + 1)] != 0]))
}

thresholding.by.h <- function(x, cps, lambda, s) {
    T <- ncol(x)
    N <- nrow(x)
    keep <- list()
    delta <- list()
    for (k in 1:length(cps)) {
        cp <- cps[[k]]
        J <- length(cp)
        start <- c(1, cp + 1)
        end <- c(cp, T)
        means <- matrix(NA, N, J + 1)
        for (i in 1:(J + 1)) {
            means[, i] <- rowMeans(x[, start[i]:end[i], drop = FALSE])
        }
        mu <- sign(means) * pmax(abs(means) - lambda, 0)
        keep[[k]] <- (mu[, 1:J, drop = FALSE] != 0) | (mu[, 2:(J + 1), drop = FALSE] != 0)
        delta[[k]] <- mu[, 2:(J + 1), drop = FALSE] - mu[, 1:J, drop = FALSE]
    }
    res <- list()
    for (i in 1:N) {
        cp <- list()
        del <- list()
        for(k in 1:length(cps)) {
            cp[[k]] <- cps[[k]][keep[[k]][i, ]]
            del[[k]] <- delta[[k]][i, keep[[k]][i, ]]
        }
        temp <- cp[[1]]
        temp.d <- del[[1]]
        if(length(cps) > 1) for(k in 2:length(cps)) {
            if(length(cp[[k]]) > 0) for(j in 1:length(cp[[k]])) {
                cc <- cp[[k]][j]
                dd <- del[[k]][j]
                ii <- which((temp %in% (cc + ((-s):s))) & (temp.d * dd > 0))
                if(length(ii) > 0) {
                    tog <- c(cc, temp[ii])
                    tog.d <- c(dd, temp.d[ii])
                    temp <- temp[-ii]
                    temp.d <- temp.d[-ii]
                    jj <- which.max(abs(tog.d))
                    temp <- c(temp, tog[jj])
                    or <- order(temp)
                    temp <- temp[or]
                    temp.d <- c(temp.d, tog.d[[jj]])[or]
                } else {
                    temp <- c(temp, cc)
                    temp.d <- c(temp.d, dd)
                    or <- order(temp)
                    temp <- temp[or]
                    temp.d <- temp.d[or]
                }
            }
        }
        res[[i]] <- temp
    }
    return(res)
}



thresholding.by.h.2 <- function(x, cps, lambda, s) {
    T <- ncol(x)
    N <- nrow(x)

    cp.all <- sort(unique(unlist(cps)))    
    keep <- list()
    delta <- list()
    for (k in 1:length(cps)) {
        cp <- cps[[k]]
        J <- length(cp)
        start <- c(1, cp + 1)
        end <- c(cp, T)
        means <- matrix(NA, N, J + 1)
        for (i in 1:(J + 1)) {
            means[, i] <- rowMeans(x[, start[i]:end[i], drop = FALSE])
        }

        kp <- matrix(FALSE, N, J)
        de <- matrix(0, N, J)
        for(i in 1:N) {
            kk <- 1:J
            mm <- means[i, ]
            ll <- end - start + 1
            dd <- mm[2:(J + 1)] - mm[1:J]
            while(length(kk) > 0 && any(abs(dd) < lambda[i, k])) {
                j <- which.min(abs(dd))
                mm[j] <- (mm[j] * ll[j] + mm[j + 1] * ll[j + 1])/(ll[j] + ll[j + 1])
                ll[j] <- ll[j] + ll[j + 1]
                ll <- ll[-(j + 1)]
                mm <- mm[-(j + 1)]
                kk <- kk[-j]
                dd <- mm[2:length(mm)] - mm[1:(length(mm) - 1)]
            }
            kp[i, kk] <- TRUE
            de[i, kk] <- dd
        }

        keep[[k]] <- kp
        delta[[k]] <- de
        #mu <- sign(means) * pmax(abs(means) - lambda, 0)
        #keep[[k]] <- (mu[, 1:J, drop = FALSE] != 0) | (mu[, 2:(J + 1), drop = FALSE] != 0)
        #delta[[k]] <- mu[, 2:(J + 1), drop = FALSE] - mu[, 1:J, drop = FALSE]
    }
    
    res <- list()
    for (i in 1:N) {
        cp <- list()
        del <- list()
        for(k in 1:length(cps)) {
            cp[[k]] <- cps[[k]][keep[[k]][i, ]]
            del[[k]] <- delta[[k]][i, keep[[k]][i, ]]
        }
        temp <- cp[[1]]
        temp.d <- del[[1]]
        if(length(cps) > 1) for(k in 2:length(cps)) {
            if(length(cp[[k]]) > 0) for(j in 1:length(cp[[k]])) {
                cc <- cp[[k]][j]
                dd <- del[[k]][j]
                ii <- which((temp %in% (cc + ((-s):s))) & (temp.d * dd > 0))
                if(length(ii) > 0) {
                    tog <- c(cc, temp[ii])
                    tog.d <- c(dd, temp.d[ii])
                    temp <- temp[-ii]
                    temp.d <- temp.d[-ii]
                    jj <- which.max(abs(tog.d))
                    temp <- c(temp, tog[jj])
                    or <- order(temp)
                    temp <- temp[or]
                    temp.d <- c(temp.d, tog.d[[jj]])[or]
                } else {
                    temp <- c(temp, cc)
                    temp.d <- c(temp.d, dd)
                    or <- order(temp)
                    temp <- temp[or]
                    temp.d <- temp.d[or]
                }
            }
        }
        res[[i]] <- temp
    }
    return(res)
}

median_polish<-function(x, flank_size = 10, multiplier = 4.5){
  N <- nrow(x)
  T <- ncol(x)
  x.mod <- array(NA, c(N, T))
  for (i in 1:N){
    y <- .Call("median_2", x[i, ], T , flank_size, as.numeric(multiplier))
    x.mod[i, ] <- y
  }
return(x.mod)
}

thresholding.by.h.3 <- function(x, H, cps, delta){
	N <- nrow(x)
	T <- ncol(x)
	h.order <- order(H, decreasing = TRUE)
	cp.all <- cps[[h.order[1]]]
	if (length(H)>1){
	for (i in 2:length(h.order)){
		cp.temp <- cps[[h.order[i]]]
		for (j in 1:length(cp.temp)){
			if (!any(cp.all %in% (cp.temp[j] + c(-H[h.order[i]]:H[h.order[i]])))){
				cp.all <- c(cp.all, cp.temp[j])
			}
		}
	}
	}
	cp.all <- sort(unique(cp.all))

        J <- length(cp.all)
        start <- c(1, cp.all + 1)
        end <- c(cp.all, T)
        means <- matrix(NA, N, J + 1)
        for (j in 1:(J + 1)) {
            means[, j] <- rowMeans(x[, start[j]:end[j], drop = FALSE])
        }
#        kp <- matrix(FALSE, N, J)
#        de <- matrix(0, N, J)
	keep<-list()
        for(i in 1:N) {
            kk <- 1:J
            mm <- means[i, ]
            ll <- end - start + 1
            dd <- mm[2:(J + 1)] - mm[1:J]
            while(length(kk) > 0 && any(abs(dd) < delta[i])) {
                j <- which.min(abs(dd))
                mm[j] <- (mm[j] * ll[j] + mm[j + 1] * ll[j + 1])/(ll[j] + ll[j + 1])
                ll[j] <- ll[j] + ll[j + 1]
                ll <- ll[-(j + 1)]
                mm <- mm[-(j + 1)]
                kk <- kk[-j]
                dd <- mm[2:length(mm)] - mm[1:(length(mm) - 1)]
            }
#            kp[i, kk] <- TRUE
#            de[i, kk] <- dd
	    keep[[i]] <- cp.all[kk]
        }
#        keep[[k]] <- kp
#        delta[[k]] <- de
return(keep)
}
	
backward_del <- function(means, candidates, seq_length, cutoff){
	num_seg <- length(means)
	if (length(candidates) != num_seg - 1){
		stop("Lengths of candidates and segments do not match.")	
        }
	seg_len <- c(candidates, seq_length) - c(0, candidates)
	if (any(seg_len <= 0)) stop("Wrong input of candidates.")
	output <- .Call("backward_del_cutoff", as.numeric(means), as.numeric(seg_len), as.integer(c(0, candidates)), as.integer(num_seg), as.numeric(cutoff))
        if (length(output[[1]])==1 & output[[1]][1]==0){
		output_alt <- list()
		for (i in 1:3){
			output_alt[[i]] <- vector(mode = "numeric")
		}
		return(output_alt)
	}else
	{	
		for (i in 1:3){
			output[[i]] <- output[[i]][-1]
		}
		return(output)
	}
}

thresholding.by.h.3.new <- function(x, H, cps, delta){
	N <- nrow(x)
	T <- ncol(x)
	h.order <- order(H, decreasing = TRUE)
	cp.all <- cps[[h.order[1]]]
	if (length(H)>1){
	for (i in 2:length(h.order)){
		cp.temp <- cps[[h.order[i]]]
		for (j in 1:length(cp.temp)){
			if (!any(cp.all %in% (cp.temp[j] + c(-H[h.order[i]]:H[h.order[i]])))){
				cp.all <- c(cp.all, cp.temp[j])
			}
		}
	}
	}
	cp.all <- sort(unique(cp.all))

        J <- length(cp.all)
        start <- c(1, cp.all + 1)
        end <- c(cp.all, T)
        means <- matrix(NA, N, J + 1)
        for (j in 1:(J + 1)) {
            means[, j] <- rowMeans(x[, start[j]:end[j], drop = FALSE])
        }
	
	keep<-list()
	for (i in 1:N){
		keep[[i]] <- backward_del(means[i, ], cp.all, T, delta[i])[[1]]
	}
#        for(i in 1:N) {
#            kk <- 1:J
#            mm <- means[i, ]
#            ll <- end - start + 1
#            dd <- mm[2:(J + 1)] - mm[1:J]
#            while(length(kk) > 0 && any(abs(dd) < delta[i])) {
#                j <- which.min(abs(dd))
#                ll[j] <- ll[j] + ll[j + 1]
#                mm[j] <- (mm[j] * ll[j] + mm[j + 1] * ll[j + 1])/(ll[j] + ll[j + 1])
#                ll <- ll[-(j + 1)]
#                mm <- mm[-(j + 1)]
#                kk <- kk[-j]
#                dd <- mm[2:length(mm)] - mm[1:(length(mm) - 1)]
#            }
#	    keep[[i]] <- cp.all[kk]
#        }
return(keep)
}


thresholding.by.h.2.new <- function(x, cps, lambda, s) {
    T <- ncol(x)
    N <- nrow(x)

    cp.all <- sort(unique(unlist(cps)))    
    keep <- list()
    delta <- list()
    for (k in 1:length(cps)) {
        cp <- cps[[k]]
        J <- length(cp)
        start <- c(1, cp + 1)
        end <- c(cp, T)
        means <- matrix(NA, N, J + 1)
        for (i in 1:(J + 1)) {
            means[, i] <- rowMeans(x[, start[i]:end[i], drop = FALSE])
        }

        kp <- matrix(FALSE, N, J)
        de <- matrix(0, N, J)
        for(i in 1:N) {
		output <- backward_del(means[i, ], cp, T, lambda[i, k])
		kk <- output[[2]]
		dd <- output[[3]]
#            kk <- 1:J
#            mm <- means[i, ]
#            ll <- end - start + 1
#            dd <- mm[2:(J + 1)] - mm[1:J]
#            while(length(kk) > 0 && any(abs(dd) < lambda[i, k])) {
#                j <- which.min(abs(dd))
#                mm[j] <- (mm[j] * ll[j] + mm[j + 1] * ll[j + 1])/(ll[j] + ll[j + 1])
#                ll[j] <- ll[j] + ll[j + 1]
#                ll <- ll[-(j + 1)]
#                mm <- mm[-(j + 1)]
#                kk <- kk[-j]
#                dd <- mm[2:length(mm)] - mm[1:(length(mm) - 1)]
#            }
            kp[i, kk] <- TRUE
            de[i, kk] <- dd
        }

        keep[[k]] <- kp
        delta[[k]] <- de
        #mu <- sign(means) * pmax(abs(means) - lambda, 0)
        #keep[[k]] <- (mu[, 1:J, drop = FALSE] != 0) | (mu[, 2:(J + 1), drop = FALSE] != 0)
        #delta[[k]] <- mu[, 2:(J + 1), drop = FALSE] - mu[, 1:J, drop = FALSE]
    }
    
    res <- list()
    for (i in 1:N) {
        cp <- list()
        del <- list()
        for(k in 1:length(cps)) {
            cp[[k]] <- cps[[k]][keep[[k]][i, ]]
            del[[k]] <- delta[[k]][i, keep[[k]][i, ]]
        }
        temp <- cp[[1]]
        temp.d <- del[[1]]
        if(length(cps) > 1) for(k in 2:length(cps)) {
            if(length(cp[[k]]) > 0) for(j in 1:length(cp[[k]])) {
                cc <- cp[[k]][j]
                dd <- del[[k]][j]
                ii <- which((temp %in% (cc + ((-s):s))) & (temp.d * dd > 0))
                if(length(ii) > 0) {
                    tog <- c(cc, temp[ii])
                    tog.d <- c(dd, temp.d[ii])
                    temp <- temp[-ii]
                    temp.d <- temp.d[-ii]
                    jj <- which.max(abs(tog.d))
                    temp <- c(temp, tog[jj])
                    or <- order(temp)
                    temp <- temp[or]
                    temp.d <- c(temp.d, tog.d[[jj]])[or]
                } else {
                    temp <- c(temp, cc)
                    temp.d <- c(temp.d, dd)
                    or <- order(temp)
                    temp <- temp[or]
                    temp.d <- temp.d[or]
                }
            }
        }
        res[[i]] <- temp
    }
    return(res)
}
