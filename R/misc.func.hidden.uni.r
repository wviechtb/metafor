############################################################################

# function to calculate
# solve(t(X) %*% W %*% X) = .invcalc(X=X, W=W, k=k)
# solve(t(X) %*% X)       = .invcalc(X=X, W=diag(k), k=k)
# via the QR decomposition

.invcalc <- function(X, W, k) {

   sWX <- sqrt(W) %*% X
   res.qrs <- qr.solve(sWX, diag(k))
   #res.qrs <- try(qr.solve(sWX, diag(k)), silent=TRUE)
   #if (inherits(res.qrs, "try-error"))
   #   stop("Cannot compute QR decomposition.")
   return(tcrossprod(res.qrs))

}

############################################################################

# function for confint.rma.uni() with the Q-profile method and for the PM estimator

.QE.func <- function(tau2val, Y, vi, X, k, objective, verbose=FALSE, digits=4) {

   mstyle <- .get.mstyle()

   if (any(tau2val + vi < 0))
      stop(mstyle$stop("Some marginal variances are negative."), call.=FALSE)

   W     <- diag(1/(vi + tau2val), nrow=k, ncol=k)
   stXWX <- .invcalc(X=X, W=W, k=k)
   P     <- W - W %*% X %*% stXWX %*% crossprod(X,W)
   RSS   <- crossprod(Y,P) %*% Y

   if (verbose)
      cat(mstyle$verbose(paste("tau2 =", fmtx(tau2val, digits[["var"]], addwidth=4), "  RSS - objective =", fmtx(RSS - objective, digits[["var"]], flag=" "), "\n")))

   return(RSS - objective)

}

############################################################################

# function for confint.rma.uni() with method="GENQ"

.GENQ.func <- function(tau2val, P, vi, Q, level, k, p, getlower, verbose=FALSE, digits=4) {

   mstyle <- .get.mstyle()

   S <- diag(sqrt(vi + tau2val), nrow=k, ncol=k)
   lambda <- Re(eigen(S %*% P %*% S, symmetric=TRUE, only.values=TRUE)$values)
   tmp <- CompQuadForm::farebrother(Q, lambda[seq_len(k-p)])

   # starting with version 1.4.2 of CompQuadForm, the element is called 'Qq' (before it was called 'res')
   # this way, things should work regardless of the version of CompQuadForm that is installed

   if (exists("res", tmp))
      tmp$Qq <- tmp$res

   if (getlower) {
      res <- tmp$Qq - level
   } else {
      res <- (1 - tmp$Qq) - level
   }

   if (verbose)
      cat(mstyle$verbose(paste("tau2 =", fmtx(tau2val, digits[["var"]], addwidth=4), "  objective =", fmtx(res, digits[["var"]], flag=" "), "\n")))

   return(res)

}

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
