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

.re.fit.quick <- function(yi, vi, method, threshold=10^-5, maxiter=100, addei=TRUE, addhi=TRUE) {

   k <- length(yi)
   wi <- 1 / vi
   sumwi <- sum(wi)
   theta <- sum(wi * yi) / sumwi
   Q <- sum(wi * (yi - theta)^2)

   if (method == "EE")
      tau2 <- 0

   if (method == "HS")
      tau2 <- max(0, (Q - k) / sumwi)

   if (method == "HSk")
      tau2 <- max(0, (k/(k-1)*Q - k) / sumwi)

   if (method == "HE")
      tau2 <- max(0, var(yi) - mean(vi))

   if (method %in% c("DL","ML","REML","EB"))
      tau2 <- max(0, (Q - (k-1)) / (sumwi - sum(wi^2) / sumwi))

   if (method == "SJ") {
      tau2.0 <- var(yi) * (k-1) / k
      wi <- 1 / (vi / tau2.0 + 1)
      mu <- sum(wi*yi) / sum(wi)
      tau2 <- sum(wi * (yi - mu)^2) / (k-1)
   }

   if (method %in% c("ML","REML","EB")) {

      tau2 <- max(0.01, tau2)

      diff <- 1
      conv <- 1
      iter <- 1

      while (diff > threshold) {
         if (iter > maxiter) {
            conv <- 0
            break
         }
         tau2.old <- tau2
         wi <- 1 / (tau2 + vi)
         sumwi <- sum(wi)
         mu <- sum(wi*yi) / sumwi
         if (method == "ML")
            tau2 <- sum(wi^2 * ((yi - mu)^2 - vi)) / sum(wi^2)
         if (method == "REML")
            tau2 <- sum(wi^2 * ((yi - mu)^2 - vi)) / sum(wi^2) + 1 / sumwi
         if (method == "EB")
            tau2 <- sum(wi * (k/(k-1) * ((yi - mu)^2) - vi)) / sumwi
         tau2 <- max(0, tau2)
         diff <- abs(tau2 - tau2.old)
         iter <- iter + 1
      }

      if (conv == 0)
         stop()

   }

   wi <- 1 / (tau2 + vi)
   sumwi <- sum(wi)
   mu <- sum(wi*yi) / sumwi

   if (method == "REML") {
      ll <- sum(dnorm(yi, mean=mu, sd=sqrt(vi + tau2), log=TRUE)) + 1/2 * log(2*k*pi) - 1/2 * log(sumwi)
   } else {
      method <- "ML"
      ll <- sum(dnorm(yi, mean=mu, sd=sqrt(vi + tau2), log=TRUE))
   }

   fit.stats <- matrix(NA_real_, nrow=5, ncol=2)
   dimnames(fit.stats) <- list(c("ll","dev","AIC","BIC","AICc"), c("ML","REML"))
   fit.stats["ll",method] <- ll

   out <- list(beta=mu, vb=1/sumwi, tau2=tau2, Q=Q, k=k, fit.stats=fit.stats)

   if (addei)
      out$ei <- yi - mu

   if (addhi)
      out$hi <- wi / sumwi

   return(out)

}

############################################################################
