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
