.fsn.fisher <- function(fsnum, pi, alpha) {
   k <- length(pi)
   X2 <- -2*sum(log(c(pi, rep(0.5, fsnum))))
   return(pchisq(X2, df=2*(k+fsnum), lower.tail=FALSE) - alpha)
}

.fsn.fitre <- function(yi, vi) {

   k     <- length(yi)
   wi    <- 1/vi
   sumwi <- sum(wi)
   est   <- sum(wi*yi)/sumwi
   Q     <- sum(wi * (yi - est)^2)
   tau2  <- max(0, (Q - (k-1)) / (sumwi - sum(wi^2)/sumwi))
   wi    <- 1 / (vi + tau2)
   sumwi <- sum(wi)
   est   <- sum(wi*yi)/sumwi
   se    <- sqrt(1 / sumwi)
   zval  <- est / se
   pval  <- 2*pnorm(abs(zval), lower.tail=FALSE)

   return(list(est=est, se=se, zval=zval, pval=pval, tau2=tau2))

}

.fsn.fitnew <- function(new, yi, vi, vnew, tau2, alpha, iters) {

   new <- ceiling(new)

   mus   <- rep(NA_real_, iters)
   pvals <- rep(NA_real_, iters)

   for (j in seq_len(iters)) {
      yinew <- c(yi, rnorm(new, 0, sqrt(vnew+tau2)))
      vinew <- c(vi, rep(vnew, new))
      tmp <- .fsn.fitre(yinew, vinew)
      mus[j] <- tmp$est
      pvals[j] <- tmp$pval
   }

   return(list(mean = mean(mus), rejrate = mean(pvals <= alpha)))

}

.fsn.re <- function(fsnum, yi, vi, vnew, tau2, target, alpha, iters, verbose=FALSE) {

   fsnum <- ceiling(fsnum)
   tmp <- .fsn.fitnew(fsnum, yi, vi, vnew, tau2, alpha, iters)
   est <- tmp$mean
   diff <- est - target
   if (verbose)
      cat("fsnum =", formatC(fsnum, width=4, format="d"), "  est =", fmtx(est), "  target =", fmtx(target), "  diff =", fmtx(diff, flag=" "), "\n")
   return(diff)

}
