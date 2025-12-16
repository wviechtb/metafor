############################################################################

.mapfun.alpha <- function(x, lb, ub) {
   if (is.infinite(lb) || is.infinite(ub)) {
      x
   } else {
      lb + (ub-lb) / (1 + exp(-x)) # map (-inf,inf) to (lb,ub)
   }
}

.mapinvfun.alpha <- function(x, lb, ub) {
   if (is.infinite(lb) || is.infinite(ub)) {
      x
   } else {
      log((x-lb)/(ub-x))
   }
}

############################################################################

# -1 times the log-likelihood (regular or restricted) for location-scale models

.ll.rma.ls <- function(par, yi, vi, X, Z, reml, k, pX,
                       alpha.arg, beta.arg, verbose, digits,
                       REMLf, link, mZ, alpha.min, alpha.max, alpha.transf,
                       tau2.min, tau2.max, optbeta) {

   mstyle <- .get.mstyle()

   if (optbeta) {
      beta  <- par[seq_len(pX)]
      beta  <- ifelse(is.na(beta.arg), beta, beta.arg)
      alpha <- par[-seq_len(pX)]
   } else {
      alpha <- par
   }

   if (alpha.transf)
      alpha <- mapply(.mapfun.alpha, alpha, alpha.min, alpha.max)

   alpha <- ifelse(is.na(alpha.arg), alpha, alpha.arg)

   # compute predicted tau2 values

   if (link == "log") {
      tau2 <- exp(c(Z %*% alpha))
   } else {
      tau2 <- c(Z %*% alpha)
   }

   if (any(is.na(tau2)) || any(tau2 < tau2.min) || any(tau2 > tau2.max)) {

      llval <- -Inf
      llcomp <- FALSE

   } else {

      llcomp <- TRUE

      if (any(tau2 < 0)) {

         llval  <- -Inf
         llcomp <- FALSE

      } else {

         # compute weights / weight matrix

         wi <- 1/(vi + tau2)
         W <- diag(wi, nrow=k, ncol=k)

         if (!optbeta) {

            stXWX <- try(.invcalc(X=X, W=W, k=k), silent=TRUE)

            if (inherits(stXWX, "try-error")) {

               llval  <- -Inf
               llcomp <- FALSE

            } else {

               beta <- stXWX %*% crossprod(X,W) %*% as.matrix(yi)

            }

         }

      }

   }

   if (llcomp) {

      # compute residual sum of squares

      RSS <- sum(wi*c(yi - X %*% beta)^2)

      # compute log-likelihood

      if (!reml) {
         llval <- -1/2 * (k) * log(2*base::pi) - 1/2 * sum(log(vi + tau2)) - 1/2 * RSS
      } else {
         llval <- -1/2 * (k-pX) * log(2*base::pi) + ifelse(REMLf, 1/2 * determinant(crossprod(X), logarithm=TRUE)$modulus, 0) +
                  -1/2 * sum(log(vi + tau2)) - 1/2 * determinant(crossprod(X,W) %*% X, logarithm=TRUE)$modulus - 1/2 * RSS
      }

   }

   if (!is.null(mZ))
      alpha <- mZ %*% alpha

   if (verbose) {
      cat(mstyle$verbose(paste0("ll = ",            fmtx(llval, digits[["fit"]], flag=" "), "  ")))
      if (optbeta)
         cat(mstyle$verbose(paste0("beta = ", paste(fmtx(beta,  digits[["est"]], flag=" "), collapse=" "), "  ")))
      cat(mstyle$verbose(paste0("alpha = ",   paste(fmtx(alpha, digits[["est"]], flag=" "), collapse=" "))))
      cat("\n")
   }

   return(-1 * llval)

}

.rma.ls.ineqfun.pos <- function(par, yi, vi, X, Z, reml, k, pX, alpha.arg, beta.arg, verbose, digits, REMLf, link, mZ, alpha.min, alpha.max, alpha.transf, tau2.min, tau2.max, optbeta) {

   if (optbeta) {
      alpha <- par[-seq_len(pX)]
   } else {
      alpha <- par
   }

   if (alpha.transf)
      alpha <- mapply(.mapfun.alpha, alpha, alpha.min, alpha.max)

   alpha <- ifelse(is.na(alpha.arg), alpha, alpha.arg)

   tau2 <- c(Z %*% alpha)

   return(tau2)

}

.rma.ls.ineqfun.neg <- function(par, yi, vi, X, Z, reml, k, pX, alpha.arg, beta.arg, verbose, digits, REMLf, link, mZ, alpha.min, alpha.max, alpha.transf, tau2.min, tau2.max, optbeta) {

   if (optbeta) {
      alpha <- par[-seq_len(pX)]
   } else {
      alpha <- par
   }

   if (alpha.transf)
      alpha <- mapply(.mapfun.alpha, alpha, alpha.min, alpha.max)

   alpha <- ifelse(is.na(alpha.arg), alpha, alpha.arg)

   tau2 <- -c(Z %*% alpha)

   return(tau2)

}

############################################################################
