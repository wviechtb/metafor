############################################################################

.selmodel.pval <- function(yi, vi, alternative) {
   zi <- yi / sqrt(vi)
   if (alternative == "two.sided") {
      pval <- 2 * pnorm(abs(zi), lower.tail=FALSE)
   } else {
      pval <- pnorm(zi, lower.tail = alternative == "less")
   }
   return(pval)
}

.selmodel.verbose <- function(ll, beta, tau2, delta, mstyle, digits) {
   cat(mstyle$verbose(paste0("ll = ",         ifelse(is.na(ll),    NA, formatC(ll,    digits=digits[["fit"]], format="f", flag=" ")),                "  ")))
   cat(mstyle$verbose(paste0("beta =",  paste(ifelse(is.na(beta),  NA, formatC(beta,  digits=digits[["est"]], format="f", flag=" ")), collapse=" "), "  ")))
   cat(mstyle$verbose(paste0("tau2 =",        ifelse(is.na(tau2),  NA, formatC(tau2,  digits=digits[["var"]], format="f", flag=" ")),                "  ")))
   cat(mstyle$verbose(paste0("delta =", paste(ifelse(is.na(delta), NA, formatC(delta, digits=digits[["est"]], format="f", flag=" ")), collapse=" "))))
   cat("\n")
}

.mapfun <- function(x, lb, ub, fun=NA) {
   if (is.na(fun)) {
      lb + (ub-lb) / (1 + exp(-x))
   } else {
      x <- sapply(x, fun)
      pmin(pmax(x, lb), ub)
   }
}

.mapinvfun <- function(x, lb, ub, fun=NA) {
   if (is.na(fun)) {
      log((x-lb)/(ub-x))
   } else {
      sapply(x, fun)
   }
}

############################################################################

.selmodel.int <- function(yvals, yi, vi, preci, yhat, wi.fun, delta, tau2, alternative, pval.min, steps) {
   pval <- .selmodel.pval(yvals, vi, alternative)
   pval[pval < pval.min]     <- pval.min
   pval[pval > (1-pval.min)] <- 1-pval.min
   wi.fun(pval, delta, yi, vi, preci, alternative, steps) * dnorm(yvals, yhat, sqrt(vi+tau2))
}

.selmodel.ll.cont <- function(par, yi, vi, X.fit, preci, k, pX, pvals, deltas, delta.val, delta.transf, mapfun, delta.min, delta.max, tau2.val, tau2.transf, tau2.max, beta.val, wi.fun, steps, pgrp, alternative, pval.min, intCtrl, verbose, digits, dofit=FALSE) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   beta  <- par[1:pX]
   tau2  <- par[pX+1]
   delta <- par[(pX+2):(pX+1+deltas)]

   beta[!is.na(beta.val)] <- beta.val

   if (tau2.transf)
      tau2 <- exp(tau2)

   tau2[!is.na(tau2.val)] <- tau2.val

   tau2[tau2 < .Machine$double.eps*10] <- 0
   tau2[tau2 > tau2.max] <- tau2.max

   if (delta.transf)
      delta <- mapply(.mapfun, delta, delta.min, delta.max, mapfun)

   delta <- ifelse(is.na(delta.val), delta, delta.val)

   yhat <- c(X.fit %*% beta)

   Ai <- rep(NA_real_, k)
   for (i in 1:k) {
      tmp <- try(integrate(.selmodel.int, lower=intCtrl$lower, upper=intCtrl$upper, yi=yi[i], vi=vi[i], preci=preci[i], yhat=yhat[i], wi.fun=wi.fun, delta=delta, tau2=tau2, alternative=alternative, pval.min=pval.min, steps=steps, subdivisions=intCtrl$subdivisions, rel.tol=intCtrl$rel.tol)$value, silent=TRUE)
      if (inherits(tmp, "try-error"))
         stop(mstyle$stop(paste0("Could not integrate over density in study ", i, ".")), call.=FALSE)
      Ai[i] <- tmp
   }

   ll <- sum(log(wi.fun(pvals, delta, yi, vi, preci, alternative, steps)) + dnorm(yi, yhat, sqrt(vi+tau2), log=TRUE) - log(Ai))

   if (dofit) {

      res <- list(ll=ll, beta=beta, tau2=tau2, delta=delta)
      return(res)

   }

   if (verbose)
      .selmodel.verbose(ll=ll, beta=beta, tau2=tau2, delta=delta, mstyle=mstyle, digits=digits)

   if (verbose > 2) {
      xs <- seq(pval.min, 1-pval.min, length=101)
      ys <- wi.fun(xs, delta, yi, vi, preci=1, alternative, steps)
      plot(xs, ys, type="l", lwd=2, xlab="p-value", ylab="Relative Likelihood of Selection")
   }

   return(-1*ll)

}

############################################################################

.selmodel.ll.stepfun <- function(par, yi, vi, X.fit, preci, k, pX, pvals, deltas, delta.val, delta.transf, mapfun, delta.min, delta.max, tau2.val, tau2.transf, tau2.max, beta.val, wi.fun, steps, pgrp, alternative, pval.min, intCtrl, verbose, digits, dofit=FALSE) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   beta  <- par[1:pX]
   tau2  <- par[pX+1]
   delta <- par[(pX+2):(pX+1+deltas)]

   beta[!is.na(beta.val)] <- beta.val

   if (tau2.transf)
      tau2 <- exp(tau2)

   tau2[!is.na(tau2.val)] <- tau2.val

   tau2[tau2 < .Machine$double.eps*10] <- 0
   tau2[tau2 > tau2.max] <- tau2.max

   if (delta.transf)
      delta <- mapply(.mapfun, delta, delta.min, delta.max, mapfun)

   delta <- ifelse(is.na(delta.val), delta, delta.val)

   yhat <- c(X.fit %*% beta)

   N <- length(steps)

   Ai <- rep(NA_real_, k)
   if (alternative == "greater") {
      for (i in 1:k) {
         sei <- sqrt(vi[i]+tau2)
         Ai[i] <- pnorm(qnorm(steps[1], 0, sqrt(vi[i]), lower.tail=FALSE), yhat[i], sei, lower.tail=FALSE)
         for (j in 2:N) {
            if (j < N) {
               Ai[i] <- Ai[i] + delta[j] / preci[i] * (pnorm(qnorm(steps[j], 0, sqrt(vi[i]), lower.tail=FALSE), yhat[i], sei, lower.tail=FALSE) - pnorm(qnorm(steps[j-1], 0, sqrt(vi[i]), lower.tail=FALSE), yhat[i], sei, lower.tail=FALSE))
            } else {
               Ai[i] <- Ai[i] + delta[j] / preci[i] * pnorm(qnorm(steps[j-1], 0, sqrt(vi[i]), lower.tail=FALSE), yhat[i], sei, lower.tail=TRUE)
            }
         }
      }
   }
   if (alternative == "less") {
      for (i in 1:k) {
         sei <- sqrt(vi[i]+tau2)
         Ai[i] <- pnorm(qnorm(steps[1], 0, sqrt(vi[i]), lower.tail=TRUE), yhat[i], sei, lower.tail=TRUE)
         for (j in 2:N) {
            if (j < N) {
               Ai[i] <- Ai[i] + delta[j] / preci[i] * (pnorm(qnorm(steps[j], 0, sqrt(vi[i]), lower.tail=TRUE), yhat[i], sei, lower.tail=TRUE) - pnorm(qnorm(steps[j-1], 0, sqrt(vi[i]), lower.tail=TRUE), yhat[i], sei, lower.tail=TRUE))
            } else {
               Ai[i] <- Ai[i] + delta[j] / preci[i] * pnorm(qnorm(steps[j-1], 0, sqrt(vi[i]), lower.tail=TRUE), yhat[i], sei, lower.tail=FALSE)
            }
         }
      }
   }
   if (alternative == "two.sided") {
      for (i in 1:k) {
         sei <- sqrt(vi[i]+tau2)
         Ai[i] <- pnorm(qnorm(steps[1]/2, 0, sqrt(vi[i]), lower.tail=FALSE), yhat[i], sei, lower.tail=FALSE) + pnorm(qnorm(steps[1]/2, 0, sqrt(vi[i]), lower.tail=TRUE), yhat[i], sei, lower.tail=TRUE)
         for (j in 2:N) {
            if (j < N) {
               Ai[i] <- Ai[i] + delta[j] / preci[i] * ((pnorm(qnorm(steps[j]/2, 0, sqrt(vi[i]), lower.tail=FALSE), yhat[i], sei, lower.tail=FALSE) - pnorm(qnorm(steps[j-1]/2, 0, sqrt(vi[i]), lower.tail=FALSE), yhat[i], sei, lower.tail=FALSE)) + (pnorm(qnorm(steps[j]/2, 0, sqrt(vi[i]), lower.tail=TRUE), yhat[i], sei, lower.tail=TRUE) - pnorm(qnorm(steps[j-1]/2, 0, sqrt(vi[i]), lower.tail=TRUE), yhat[i], sei, lower.tail=TRUE)))
            } else {
               Ai[i] <- Ai[i] + delta[j] / preci[i] * (pnorm(qnorm(steps[j-1]/2, 0, sqrt(vi[i]), lower.tail=FALSE), yhat[i], sei, lower.tail=TRUE) - pnorm(qnorm(steps[j-1]/2, 0, sqrt(vi[i]), lower.tail=TRUE), yhat[i], sei, lower.tail=TRUE))
            }
         }
      }
   }

   ll <- sum(log(delta[pgrp] / preci) + dnorm(yi, yhat, sqrt(vi+tau2), log=TRUE) - log(Ai))

   if (dofit) {

      res <- list(ll=ll, beta=beta, tau2=tau2, delta=delta)
      return(res)

   }

   if (verbose)
      .selmodel.verbose(ll=ll, beta=beta, tau2=tau2, delta=delta, mstyle=mstyle, digits=digits)

   if (verbose > 2) {
      xs <- seq(0, 1, length=101)
      ys <- wi.fun(xs, delta, yi, vi, preci=1, alternative, steps)
      plot(xs, ys, type="l", lwd=2, xlab="p-value", ylab="Relative Likelihood of Selection")
   }

   return(-1*ll)

}

############################################################################
