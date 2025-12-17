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

.ll.rma.ls <- function(par, yi, vi, X, Z, reml,
                       alpha.arg, beta.arg, omega2.arg, verbose, digits,
                       REMLf, link, mZ, alpha.min, alpha.max, alpha.transf, omega2.transf,
                       tau2.min, tau2.max, optbeta, reshet, mfmaxit) {

   mstyle <- .get.mstyle()

   k <- length(yi)
   p <- ncol(X)

   if (reshet) {

      omega2 <- par[length(par)]

      if (omega2.transf)
         omega2 <- exp(omega2)

      omega2[!is.na(omega2.arg)] <- omega2.arg
      omega2[omega2 < .Machine$double.eps*10] <- 0

      par <- par[-length(par)]

   } else {

      omega2 <- 0

   }

   if (optbeta) {
      beta  <- par[seq_len(p)]
      beta  <- ifelse(is.na(beta.arg), beta, beta.arg)
      alpha <- par[-seq_len(p)]
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

   if (any(is.na(tau2)) || any(tau2 < tau2.min) || any(tau2 > tau2.max) || is.na(omega2)) {

      llval <- -Inf
      llcomp <- FALSE

   } else {

      llcomp <- TRUE

      if (any(tau2 < 0)) {

         llval  <- -Inf
         llcomp <- FALSE

      } else {

         if (omega2 <= sqrt(.Machine$double.eps))
            reshet <- FALSE

         if (reshet) {

            llcomp <- FALSE

            lli <- rep(NA_real_, k)

            intfun <- function(hi, yi, vi, Xi, Zi, beta, alpha, omega2)
               dnorm(yi, mean = c(Xi %*% beta), sd = sqrt(vi + exp(c(Zi %*% alpha) + hi))) * dnorm(hi, mean = 0, sd = sqrt(omega2))

            for (i in 1:k) {
               tmp <- try(integrate(intfun, lower=-Inf, upper=Inf, stop.on.error=FALSE, yi=yi[i], vi=vi[i], Xi=X[i,,drop=FALSE], Zi=Z[i,,drop=FALSE], beta=beta, alpha=alpha, omega2=omega2), silent=TRUE)
               if (inherits(tmp, "try-error")) {
                  lli[i] <- Inf
                  break
               } else {
                  lli[i] <- tmp$value
               }
            }

            llval <- sum(log(lli), na.rm=TRUE)

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

   }

   if (llcomp) {

      # compute residual sum of squares

      RSS <- sum(wi*c(yi - X %*% beta)^2)

      # compute log-likelihood

      if (!reml) {
         llval <- -1/2 * (k) * log(2*base::pi) - 1/2 * sum(log(vi + tau2)) - 1/2 * RSS
      } else {
         llval <- -1/2 * (k-p) * log(2*base::pi) + ifelse(REMLf, 1/2 * determinant(crossprod(X), logarithm=TRUE)$modulus, 0) +
                  -1/2 * sum(log(vi + tau2)) - 1/2 * determinant(crossprod(X,W) %*% X, logarithm=TRUE)$modulus - 1/2 * RSS
      }

   }

   if (!is.null(mZ))
      alpha <- mZ %*% alpha

   iteration <- .getfromenv("iteration", default=NULL)

   if (isTRUE(iteration > mfmaxit))
      stop(mstyle$stop(paste0("Maximum number of iterations (mfmaxit=", mfmaxit, ") reached.")), call.=FALSE)

   if (verbose) {
      if (!is.null(iteration))
         cat(mstyle$verbose(paste0("Iteration ", formatC(iteration, width=5, flag="-", format="f", digits=0), " ")))
      cat(mstyle$verbose(paste0("ll = ", fmtx(llval, digits[["fit"]], flag=" "))))
      if (optbeta)
         cat(mstyle$verbose(paste0("  beta = ", paste(fmtx(beta,  digits[["est"]], flag=" "), collapse=" "))))
      cat(mstyle$verbose(paste0("  alpha = ",   paste(fmtx(alpha, digits[["est"]], flag=" "), collapse=" "))))
      if (reshet)
         cat(mstyle$verbose(paste0("  omega2 = ", paste(fmtx(omega2, digits[["var"]]), collapse=" "))))
      cat("\n")
   }

   try(assign("iteration", iteration+1, envir=.metafor), silent=TRUE)

   return(-1 * llval)

}

.rma.ls.ineqfun.pos <- function(par, yi, vi, X, Z, reml, alpha.arg, beta.arg, omega2.arg, verbose, digits, REMLf, link, mZ, alpha.min, alpha.max, alpha.transf, omega2.transf, tau2.min, tau2.max, optbeta, reshet, mfmaxit) {

   if (optbeta) {
      alpha <- par[-seq_len(ncol(X))]
   } else {
      alpha <- par
   }

   if (alpha.transf)
      alpha <- mapply(.mapfun.alpha, alpha, alpha.min, alpha.max)

   alpha <- ifelse(is.na(alpha.arg), alpha, alpha.arg)

   tau2 <- c(Z %*% alpha)

   return(tau2)

}

.rma.ls.ineqfun.neg <- function(par, yi, vi, X, Z, reml, alpha.arg, beta.arg, omega2.arg, verbose, digits, REMLf, link, mZ, alpha.min, alpha.max, alpha.transf, omega2.transf, tau2.min, tau2.max, optbeta, reshet, mfmaxit) {

   if (optbeta) {
      alpha <- par[-seq_len(ncol(X))]
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
