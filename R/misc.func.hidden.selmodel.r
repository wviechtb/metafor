############################################################################

.selmodel.pval <- function(yi, vi, alternative) {
   zi <- yi / sqrt(vi)
   if (alternative == "two.sided") {
      pval <- 2 * pnorm(abs(zi), lower.tail=FALSE)
   } else {
      pval <- pnorm(zi, lower.tail=alternative=="less")
   }
   return(pval)
}

.selmodel.verbose <- function(ll, beta, tau2, delta, mstyle, digits) {
   cat(mstyle$verbose(paste0("ll = ",         fmtx(ll,    digits[["fit"]], flag=" "),                "  ")))
   cat(mstyle$verbose(paste0("beta =",  paste(fmtx(beta,  digits[["est"]], flag=" "), collapse=" "), "  ")))
   cat(mstyle$verbose(paste0("tau2 =",        fmtx(tau2,  digits[["var"]], flag=" "),                "  ")))
   cat(mstyle$verbose(paste0("delta =", paste(fmtx(delta, digits[["est"]], flag=" "), collapse=" "))))
   cat("\n")
}

.mapfun <- function(x, lb, ub, fun=NA) {
   if (is.na(fun)) {
      if (lb==0 && ub==1) {
         plogis(x)
      } else {
         lb + (ub-lb) / (1 + exp(-x)) # map (-inf,inf) to (lb,ub)
      }
   } else {
      x <- sapply(x, fun)
      pmin(pmax(x, lb), ub)
   }
}

.mapinvfun <- function(x, lb, ub, fun=NA) {
   if (is.na(fun)) {
      if (lb==0 && ub==1) {
         qlogis(x)
      } else {
         log((x-lb)/(ub-x)) # map (lb,ub) to (-inf,inf)
      }
   } else {
      sapply(x, fun)
   }
}

.ptable <- function(pvals, steps, subset) {

   pvals[!subset] <- NA

   pgrp     <- sapply(pvals, function(p) which(p <= steps)[1])
   psteps.l <- as.character(c(0,steps[-length(steps)]))
   psteps.r <- as.character(steps)
   len.l    <- nchar(psteps.l)
   pad.l    <- sapply(max(len.l) - len.l, function(x) paste0(rep(" ", x), collapse=""))
   psteps.l <- paste0(psteps.l, pad.l)
   psteps   <- paste0(psteps.l, " < p <= ", psteps.r)
   ptable   <- table(factor(pgrp, levels=seq_along(steps), labels=psteps))
   ptable   <- data.frame(k=as.vector(ptable), row.names=names(ptable))

   return(list(pgrp=pgrp, ptable=ptable))

}

############################################################################

.selmodel.int <- function(yvals, yi, vi, preci, yhat, wi.fun, delta, tau2, alternative, pval.min, steps) {
   pval <- .selmodel.pval(yvals, vi, alternative)
   pval[pval < pval.min]     <- pval.min
   pval[pval > (1-pval.min)] <- 1-pval.min
   wi.fun(pval, delta, yi, vi, preci, alternative, steps) * dnorm(yvals, yhat, sqrt(vi+tau2))
}

.selmodel.ll.cont <- function(par, yi, vi, X, preci, subset, k, pX, pvals,
                              deltas, delta.arg, delta.transf, mapfun, delta.min, delta.max, decreasing,
                              tau2.arg, tau2.transf, tau2.max, beta.arg,
                              wi.fun, steps, pgrp,
                              alternative, pval.min, intCtrl, verbose, digits, dofit=FALSE) {

   mstyle <- .get.mstyle()

   beta  <- par[1:pX]
   tau2  <- par[pX+1]
   delta <- par[(pX+2):(pX+1+deltas)]

   beta <- ifelse(is.na(beta.arg), beta, beta.arg)

   if (tau2.transf)
      tau2 <- exp(tau2)

   tau2[!is.na(tau2.arg)] <- tau2.arg

   tau2[tau2 < .Machine$double.eps*10] <- 0
   tau2[tau2 > tau2.max] <- tau2.max

   if (delta.transf)
      delta <- mapply(.mapfun, delta, delta.min, delta.max, mapfun)

   delta <- ifelse(is.na(delta.arg), delta, delta.arg)

   yhat <- c(X %*% beta)

   Ai <- rep(NA_real_, k)

   for (i in seq_len(k)[subset]) {
      tmp <- try(integrate(.selmodel.int, lower=intCtrl$lower, upper=intCtrl$upper,
                           yi=yi[i], vi=vi[i], preci=preci[i], yhat=yhat[i], wi.fun=wi.fun,
                           delta=delta, tau2=tau2, alternative=alternative, pval.min=pval.min, steps=steps,
                           subdivisions=intCtrl$subdivisions, rel.tol=intCtrl$rel.tol)$value, silent=TRUE)
      #tmp <- try(cubintegrate(.selmodel.int, lower=intCtrl$lower, upper=intCtrl$upper,
      #                        yi=yi[i], vi=vi[i], preci=preci[i], yhat=yhat[i], wi.fun=wi.fun,
      #                        delta=delta, tau2=tau2, alternative=alternative, pval.min=pval.min, steps=steps)$integral, silent=TRUE)
      if (inherits(tmp, "try-error"))
         stop(mstyle$stop(paste0("Could not integrate over density in study ", i, ".")), call.=FALSE)
      Ai[i] <- tmp
   }

   #llval <- sum(log(wi.fun(pvals, delta, yi, vi, preci, alternative, steps)) + dnorm(yi, yhat, sqrt(vi+tau2), log=TRUE) - log(Ai))
   llval0 <- sum(                                                                                                   dnorm(yi[!subset], yhat[!subset], sqrt(vi[!subset]+tau2), log=TRUE))
   llval1 <- sum(log(wi.fun(pvals[ subset], delta, yi[ subset], vi[ subset], preci[ subset], alternative, steps)) + dnorm(yi[ subset], yhat[ subset], sqrt(vi[ subset]+tau2), log=TRUE) - log(Ai[subset]))
   llval  <- llval0 + llval1

   if (dofit) {
      res <- list(ll=llval, beta=beta, tau2=tau2, delta=delta)
      return(res)
   }

   if (verbose)
      .selmodel.verbose(ll=llval, beta=beta, tau2=tau2, delta=delta, mstyle=mstyle, digits=digits)

   if (verbose > 2) {
      xs <- seq(pval.min, 1-pval.min, length.out=1001)
      ys <- wi.fun(xs, delta, yi, vi, preci=1, alternative, steps)
      plot(xs, ys, type="l", lwd=2, xlab="p-value", ylab="Relative Likelihood of Selection")
   }

   return(-1*llval)

}

############################################################################

.selmodel.ll.stepfun <- function(par, yi, vi, X, preci, subset, k, pX, pvals,
                                 deltas, delta.arg, delta.transf, mapfun, delta.min, delta.max, decreasing,
                                 tau2.arg, tau2.transf, tau2.max, beta.arg,
                                 wi.fun, steps, pgrp,
                                 alternative, pval.min, intCtrl, verbose, digits, dofit=FALSE) {

   mstyle <- .get.mstyle()

   beta  <- par[1:pX]
   tau2  <- par[pX+1]
   delta <- par[(pX+2):(pX+1+deltas)]

   beta <- ifelse(is.na(beta.arg), beta, beta.arg)

   if (tau2.transf)
      tau2 <- exp(tau2)

   tau2[!is.na(tau2.arg)] <- tau2.arg

   tau2[tau2 < .Machine$double.eps*10] <- 0
   tau2[tau2 > tau2.max] <- tau2.max

   if (decreasing) {

      if (delta.transf) {
         delta <- exp(delta)
         delta <- cumsum(c(0, -delta[-1]))
         delta <- exp(delta)
      }

   } else {

      if (delta.transf)
         delta <- mapply(.mapfun, delta, delta.min, delta.max, mapfun)

   }

   delta <- ifelse(is.na(delta.arg), delta, delta.arg)

   if (decreasing && any(!is.na(delta.arg[-1])))
      delta <- rev(cummax(rev(delta)))

   yhat <- c(X %*% beta)

   sei <- sqrt(vi + tau2)

   N <- length(steps)

   Ai <- rep(NA_real_, k)

   if (alternative == "greater") {
      for (i in seq_len(k)[subset]) {
         Ai[i] <- pnorm(qnorm(steps[1], 0, sqrt(vi[i]), lower.tail=FALSE), yhat[i], sei[i], lower.tail=FALSE)
         for (j in 2:N) {
            if (j < N) {
               Ai[i] <- Ai[i] + delta[j] / preci[i] * (pnorm(qnorm(steps[j], 0, sqrt(vi[i]), lower.tail=FALSE), yhat[i], sei[i], lower.tail=FALSE) - pnorm(qnorm(steps[j-1], 0, sqrt(vi[i]), lower.tail=FALSE), yhat[i], sei[i], lower.tail=FALSE))
            } else {
               Ai[i] <- Ai[i] + delta[j] / preci[i] * pnorm(qnorm(steps[j-1], 0, sqrt(vi[i]), lower.tail=FALSE), yhat[i], sei[i], lower.tail=TRUE)
            }
         }
      }
   }

   if (alternative == "less") {
      for (i in seq_len(k)[subset]) {
         Ai[i] <- pnorm(qnorm(steps[1], 0, sqrt(vi[i]), lower.tail=TRUE), yhat[i], sei[i], lower.tail=TRUE)
         for (j in 2:N) {
            if (j < N) {
               Ai[i] <- Ai[i] + delta[j] / preci[i] * (pnorm(qnorm(steps[j], 0, sqrt(vi[i]), lower.tail=TRUE), yhat[i], sei[i], lower.tail=TRUE) - pnorm(qnorm(steps[j-1], 0, sqrt(vi[i]), lower.tail=TRUE), yhat[i], sei[i], lower.tail=TRUE))
            } else {
               Ai[i] <- Ai[i] + delta[j] / preci[i] * pnorm(qnorm(steps[j-1], 0, sqrt(vi[i]), lower.tail=TRUE), yhat[i], sei[i], lower.tail=FALSE)
            }
         }
      }
   }

   if (alternative == "two.sided") {
      for (i in seq_len(k)[subset]) {
         Ai[i] <- pnorm(qnorm(steps[1]/2, 0, sqrt(vi[i]), lower.tail=FALSE), yhat[i], sei[i], lower.tail=FALSE) + pnorm(qnorm(steps[1]/2, 0, sqrt(vi[i]), lower.tail=TRUE), yhat[i], sei[i], lower.tail=TRUE)
         for (j in 2:N) {
            if (j < N) {
               Ai[i] <- Ai[i] + delta[j] / preci[i] * ((pnorm(qnorm(steps[j]/2, 0, sqrt(vi[i]), lower.tail=FALSE), yhat[i], sei[i], lower.tail=FALSE) - pnorm(qnorm(steps[j-1]/2, 0, sqrt(vi[i]), lower.tail=FALSE), yhat[i], sei[i], lower.tail=FALSE)) + (pnorm(qnorm(steps[j]/2, 0, sqrt(vi[i]), lower.tail=TRUE), yhat[i], sei[i], lower.tail=TRUE) - pnorm(qnorm(steps[j-1]/2, 0, sqrt(vi[i]), lower.tail=TRUE), yhat[i], sei[i], lower.tail=TRUE)))
            } else {
               Ai[i] <- Ai[i] + delta[j] / preci[i] * (pnorm(qnorm(steps[j-1]/2, 0, sqrt(vi[i]), lower.tail=FALSE), yhat[i], sei[i], lower.tail=TRUE) - pnorm(qnorm(steps[j-1]/2, 0, sqrt(vi[i]), lower.tail=TRUE), yhat[i], sei[i], lower.tail=TRUE))
            }
         }
      }
   }

   #llval <- sum(log(delta[pgrp] / preci) + dnorm(yi, yhat, sei, log=TRUE) - log(Ai))
   llval0 <- sum(                                             dnorm(yi[!subset], yhat[!subset], sei[!subset], log=TRUE))
   llval1 <- sum(log(delta[pgrp[ subset]] / preci[ subset]) + dnorm(yi[ subset], yhat[ subset], sei[ subset], log=TRUE) - log(Ai[subset]))
   llval  <- llval0 + llval1

   if (dofit) {

      res <- list(ll=llval, beta=beta, tau2=tau2, delta=delta)
      return(res)

   }

   if (verbose)
      .selmodel.verbose(ll=llval, beta=beta, tau2=tau2, delta=delta, mstyle=mstyle, digits=digits)

   if (verbose > 2) {
      xs <- seq(0, 1, length.out=1001)
      ys <- wi.fun(xs, delta, yi, vi, preci=1, alternative, steps)
      plot(xs, ys, type="l", lwd=2, xlab="p-value", ylab="Relative Likelihood of Selection")
   }

   return(-1*llval)

}

############################################################################

.selmodel.ll.trunc <- function(par, yi, vi, X, preci, subset, k, pX, pvals,
                               deltas, delta.arg, delta.transf, mapfun, delta.min, delta.max, decreasing,
                               tau2.arg, tau2.transf, tau2.max, beta.arg,
                               wi.fun, steps, pgrp,
                               alternative, pval.min, intCtrl, verbose, digits, dofit=FALSE) {

   mstyle <- .get.mstyle()

   beta  <- par[1:pX]
   tau2  <- par[pX+1]
   delta <- par[(pX+2):(pX+1+deltas)]

   beta <- ifelse(is.na(beta.arg), beta, beta.arg)

   if (tau2.transf)
      tau2 <- exp(tau2)

   tau2[!is.na(tau2.arg)] <- tau2.arg

   tau2[tau2 < .Machine$double.eps*10] <- 0
   tau2[tau2 > tau2.max] <- tau2.max

   if (delta.transf)
      delta <- mapply(.mapfun, delta, delta.min, delta.max, mapfun)

   delta <- ifelse(is.na(delta.arg), delta, delta.arg)

   yhat <- c(X %*% beta)

   sei <- sqrt(vi + tau2)

   if (is.na(steps))
      steps <- delta[2]

   if (alternative == "greater") {
      ll0i <- dnorm(yi[!subset], yhat[!subset], sei[!subset], log=TRUE)
      ll1i <- ifelse(yi[subset] > steps, 0, log(delta[1])) + dnorm(yi[subset], yhat[subset], sei[subset], log=TRUE) - log(1 - (1-delta[1]) * pnorm(steps, yhat[subset], sei[subset], lower.tail=TRUE))
   }

   if (alternative == "less") {
      ll0i <- dnorm(yi[!subset], yhat[!subset], sei[!subset], log=TRUE)
      ll1i <- ifelse(yi[subset] < steps, 0, log(delta[1])) + dnorm(yi[subset], yhat[subset], sei[subset], log=TRUE) - log(1 - (1-delta[1]) * pnorm(steps, yhat[subset], sei[subset], lower.tail=FALSE))
   }

   llval <- sum(ll0i) + sum(ll1i)

   if (dofit) {

      res <- list(ll=llval, beta=beta, tau2=tau2, delta=delta)
      return(res)

   }

   if (verbose)
      .selmodel.verbose(ll=llval, beta=beta, tau2=tau2, delta=delta, mstyle=mstyle, digits=digits)

   return(-1*llval)

}

############################################################################

.rma.selmodel.ineqfun.pos <- function(par, yi, vi, X, preci, k, pX, pvals,
                                      deltas, delta.arg, delta.transf, mapfun, delta.min, delta.max, decreasing,
                                      tau2.arg, tau2.transf, tau2.max, beta.arg,
                                      wi.fun, steps, pgrp, alternative, pval.min, intCtrl, verbose, digits, dofit=FALSE) {

   delta <- par[-seq_len(pX+1)]

   if (delta.transf)
      delta <- mapply(.mapfun, delta, delta.min, delta.max, mapfun)

   delta <- ifelse(is.na(delta.arg), delta, delta.arg)

   diffs <- -diff(delta) # -1 * differences (delta1-delta2, delta2-delta3, ...) must be positive

   return(diffs)

}

.rma.selmodel.ineqfun.neg <- function(par, yi, vi, X, preci, k, pX, pvals,
                                      deltas, delta.arg, delta.transf, mapfun, delta.min, delta.max, decreasing,
                                      tau2.arg, tau2.transf, tau2.max, beta.arg,
                                      wi.fun, steps, pgrp, alternative, pval.min, intCtrl, verbose, digits, dofit=FALSE) {

   delta <- par[-seq_len(pX+1)]

   if (delta.transf)
      delta <- mapply(.mapfun, delta, delta.min, delta.max, mapfun)

   delta <- ifelse(is.na(delta.arg), delta, delta.arg)

   diffs <- diff(delta) # differences (delta1-delta2, delta2-delta3, ...) must be negative

   return(diffs)

}

############################################################################
