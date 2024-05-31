############################################################################

.fsn.fisher <- function(fsnum, pi, alpha) {
   k <- length(pi)
   X2 <- -2*sum(log(c(pi, rep(0.5, fsnum))))
   return(pchisq(X2, df=2*(k+fsnum), lower.tail=FALSE) - alpha)
}

############################################################################

.fsn.scale <- function(x, k) {
   if (k == 0) return(x)
   if (k == 1) return(0)
   if (k >= 2) return((x-mean(x))/sd(x))
}

.fsn.gen <- function(fsnum, yi, vi, vt, est, tau2, tau2fix, test, weighted, target, alpha,
                     exact, method, mumiss, upperint, maxint, verbose=FALSE, newest=FALSE) {

   fsnum <- floor(fsnum)

   if (fsnum > maxint)
      fsnum <- maxint

   yinew <- c(yi, .fsn.scale(rnorm(fsnum), fsnum)*sqrt(vt+tau2) + mumiss)
   vinew <- c(vi, rep(vt,fsnum))

   if (is.null(target)) {

      if (exact && fsnum <= 5000) {

         tmp <- suppressWarnings(try(rma(yinew, vinew, method=method, tau2=tau2fix, test=test, weighted=weighted), silent=TRUE))

         if (inherits(tmp, "try-error"))
            stop()

         est.fsn  <- tmp$beta[1]
         tau2.fsn <- tmp$tau2
         pval.fsn <- tmp$pval

         if (mumiss != 0 && sign(est.fsn) == sign(mumiss))
            pval.fsn <- 1

      } else {

         k <- length(yi)

         if (is.element(method, c("FE","EE","CE"))) {

            tau2.fsn <- 0

         } else {

            est.fsn <- (k*est + fsnum*mumiss) / (k + fsnum)

            if (is.null(tau2fix)) {
               tau2.fsn <- max(0, ((k-1)*tau2 + max(0,(fsnum-1))*tau2 + k*(est-est.fsn)^2 + fsnum*(mumiss-est.fsn)^2) / (k + fsnum - 1))
            } else {
               tau2.fsn <- tau2
            }

         }

         if (isTRUE(weighted)) {
            est.fsn  <- weighted.mean(yinew, 1 / (vinew + tau2.fsn))
            zval.new <- est.fsn / sqrt(1 / (sum(1 / (vi + tau2.fsn)) + fsnum / (vt + tau2.fsn)))
         } else {
            est.fsn  <- mean(yinew)
            zval.new <- (k + fsnum) * est.fsn / sqrt(sum(vi + tau2.fsn) + fsnum * (vt + tau2.fsn))
         }

         pval.fsn <- 2*pnorm(abs(zval.new), lower.tail=FALSE)

         if (mumiss != 0 && sign(est.fsn * mumiss) == 1)
            pval.fsn <- 1

      }

      if (newest) {
         return(list(est.fsn=est.fsn, tau2.fsn=tau2.fsn, pval.fsn=pval.fsn))
      } else {
         if (fsnum == maxint) {
            diff <- 0
         } else {
            diff <- pval.fsn - alpha
         }
      }

      if (verbose)
         cat("fsnum =", formatC(fsnum, width=nchar(upperint)+1, format="d"), " est =", fmtx(est.fsn, flag=" "), " tau2 =", fmtx(tau2.fsn), " pval =", fmtx(pval.fsn), " alpha =", fmtx(alpha), " diff =", fmtx(diff, flag=" "), "\n")

   } else {

      if (exact && fsnum <= 5000) {

         tmp <- suppressWarnings(try(rma(yinew, vinew, method=method, tau2=tau2fix, test=test, weighted=weighted), silent=TRUE))

         if (inherits(tmp, "try-error"))
            stop()

         est.fsn  <- tmp$beta[1]
         tau2.fsn <- tmp$tau2
         pval.fsn <- tmp$pval

      } else {

         k <- length(yi)

         if (is.element(method, c("FE","EE","CE"))) {

            tau2.fsn <- 0

         } else {

            est.fsn <- (k*est + fsnum*mumiss) / (k + fsnum)

            if (is.null(tau2fix)) {
               tau2.fsn <- ((k-1)*tau2 + max(0,(fsnum-1))*tau2 + k*(est-est.fsn)^2 + fsnum*(mumiss-est.fsn)^2) / (k + fsnum - 1)
            } else {
               tau2.fsn <- tau2
            }

         }

         if (isTRUE(weighted)) {
            est.fsn  <- weighted.mean(yinew, 1 / (vinew + tau2.fsn))
            zval.new <- est.fsn / sqrt(1 / (sum(1 / (vi + tau2.fsn)) + fsnum / (vt + tau2.fsn)))
         } else {
            est.fsn  <- mean(yinew)
            zval.new <- (k + fsnum) * est.fsn / sqrt(sum(vi + tau2.fsn) + fsnum * (vt + tau2.fsn))
         }
         pval.fsn <- 2*pnorm(abs(zval.new), lower.tail=FALSE)

      }

      if (newest) {
         return(list(est.fsn=est.fsn, tau2.fsn=tau2.fsn, pval.fsn=pval.fsn))
      } else {
         if (fsnum == maxint) {
            diff <- 0
         } else {
            diff <- est.fsn - target
         }
      }

      if (verbose)
         cat("fsnum =", formatC(fsnum, width=nchar(upperint)+1, format="d"), " est =", fmtx(est.fsn, flag=" "), " tau2 =", fmtx(tau2.fsn), " target =", fmtx(target), " diff =", fmtx(diff, flag=" "), "\n")

   }

   return(diff)

}

############################################################################

.rnd.fsn <- function(fsnum) {

   if (is.finite(fsnum) && abs(fsnum - round(fsnum)) >= .Machine$double.eps^0.5) {
      fsnum <- ceiling(fsnum)
   } else {
      fsnum <- round(fsnum)
   }

   return(fsnum)

}

############################################################################
