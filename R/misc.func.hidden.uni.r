############################################################################

### function to calculate:
### solve(t(X) %*% W %*% X) = .invcalc(X=X, W=W, k=k)
### solve(t(X) %*% X)       = .invcalc(X=X, W=diag(k), k=k)
### without taking the actual inverse

.invcalc <- function(X, W, k) {

   sWX <- sqrt(W) %*% X
   res.qrs <- qr.solve(sWX, diag(k))
   #res.qrs <- try(qr.solve(sWX, diag(k)), silent=TRUE)
   #if (inherits(res.qrs, "try-error"))
   #   stop("Cannot compute QR decomposition.")
   return(tcrossprod(res.qrs))

}

############################################################################

### function for confint.rma.uni() with Q-profile method and for the PM estimator

.QE.func <- function(tau2val, Y, vi, X, k, objective, verbose=FALSE, digits=4) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (any(tau2val + vi < 0))
      stop(mstyle$stop("Some marginal variances are negative."))

   W     <- diag(1/(vi + tau2val), nrow=k, ncol=k)
   stXWX <- .invcalc(X=X, W=W, k=k)
   P     <- W - W %*% X %*% stXWX %*% crossprod(X,W)
   RSS   <- crossprod(Y,P) %*% Y

   if (verbose)
      cat(mstyle$verbose(paste("tau2 =", formatC(tau2val, digits=digits[["var"]], width=digits[["var"]]+4, format="f"), "  RSS - objective =", formatC(RSS - objective, format="f", digits=digits[["var"]], flag=" "), "\n")))

   return(RSS - objective)

}

############################################################################

### function for confint.rma.uni() with method="GENQ"

.GENQ.func <- function(tau2val, P, vi, Q, level, k, p, getlower, verbose=FALSE, digits=4) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   S <- diag(sqrt(vi + tau2val), nrow=k, ncol=k)
   lambda <- Re(eigen(S %*% P %*% S, symmetric=TRUE, only.values=TRUE)$values)
   tmp <- CompQuadForm::farebrother(Q, lambda[seq_len(k-p)])

   ### starting with version 1.4.2 of CompQuadForm, the element is called 'Qq' (before it was called 'res')
   ### this way, things should work regardless of the version of CompQuadForm that is installed

   if (exists("res", tmp))
      tmp$Qq <- tmp$res

   if (getlower) {
      res <- tmp$Qq - level
   } else {
      res <- (1 - tmp$Qq) - level
   }

   if (verbose)
      cat(mstyle$verbose(paste("tau2 =", formatC(tau2val, digits=digits[["var"]], width=digits[["var"]]+4, format="f"), "  objective =", formatC(res, format="f", digits=digits[["var"]], flag=" "), "\n")))

   return(res)

}

############################################################################

### generate all possible permutations

# .genperms <- function(k) {
#
#    v <- seq_len(k)
#
#    sub <- function(k, v) {
#       if (k==1L) {
#          matrix(v,1,k)
#       } else {
#          X  <-  NULL
#          for(i in seq_len(k)) {
#             X <- rbind(X, cbind(v[i], Recall(k-1, v[-i])))
#          }
#       X
#       }
#    }
#
#    return(sub(k, v[seq_len(k)]))
#
# }

### generate all possible unique permutations

.genuperms <- function(x) {

   z <- NULL

   sub <- function(x, y) {
      len.x <- length(x)
      if (len.x == 0L) {
         return(y)
      } else {
         prev.num <- 0
         for (i in seq_len(len.x)) {
            num <- x[i]
            if (num > prev.num) {
               prev.num <- num
               z <- rbind(z, Recall(x[-i], c(y,num)))
            }
         }
         return(z)
      }
   }

   return(sub(x, y=NULL))

}

.permci <- function(val, obj, j, exact, iter, progbar, level, digits, control) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   ### fit model with shifted outcome
   res <- try(suppressWarnings(rma.uni(obj$yi - c(val*obj$X[,j]), obj$vi, weights=obj$weights, mods=obj$X, intercept=FALSE, method=obj$method, weighted=obj$weighted, test=obj$test, tau2=ifelse(obj$tau2.fix, obj$tau2, NA), control=obj$control, skipr2=TRUE)), silent=TRUE)

   if (inherits(res, "try-error"))
      stop()

   ### p-value based on permutation test
   pval <- permutest(res, exact=exact, iter=iter, progbar=FALSE, control=control)$pval[j]

   ### get difference between p-value and level
   diff <- pval - level / ifelse(control$alternative == "two.sided", 1, 2)

   ### show progress
   if (progbar)
      cat(mstyle$verbose(paste("pval =", formatC(pval, format="f", digits=digits[["pval"]]), " diff =", formatC(diff, format="f", digits=digits[["pval"]], flag=" "), " val =", formatC(val, format="f", digits=digits[["est"]], flag=" "), "\n")))

   ### penalize negative differences, which should force the CI bound to correspond to a p-value of *at least* level
   diff <- ifelse(diff < 0, diff*10, diff)

   return(diff)

}

############################################################################

### -1 times the log likelihood (regular or restricted) for location-scale model

.ll.rma.ls <- function(par, yi, vi, X, Z, reml, k, pX, verbose, digits, REMLf, link) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   #beta  <- par[seq_len(pX)]
   #alpha <- par[-seq_len(pX)]

   alpha <- par

   ### compute predicted tau2 values

   if (link == "log")
      tau2 <- exp(c(Z %*% alpha))
   if (link == "identity")
      tau2 <- c(Z %*% alpha)

   if (any(tau2 < 0)) {

      llval <- -Inf

   } else {

      ### compute weights
      wi <- 1/(vi + tau2)

      ### when using this, the optimization only pertains to the parameter(s) in 'alpha', as 'beta' is then fully
      ### determined by the current value(s) of 'alpha'; this is actually also how the standard RE/ME model is fitted;
      ### but is this really the best way of doing this? one could also optimize over beta and alpha jointly
      W <- diag(wi, nrow=k, ncol=k)
      stXWX <- try(.invcalc(X=X, W=W, k=k), silent=TRUE)

      if (inherits(stXWX, "try-error")) {

         llval <- -Inf

      } else {

         beta <- stXWX %*% crossprod(X,W) %*% as.matrix(yi)

         ### compute residual sum of squares
         RSS <- sum(wi*(yi - X %*% beta)^2)

         ### log-likelihood (could leave out additive constants)
         if (!reml) {
            llval <- -1/2 * (k) * log(2*base::pi) - 1/2 * sum(log(vi + tau2)) - 1/2 * RSS
         } else {
            llval <- -1/2 * (k-pX) * log(2*base::pi) + ifelse(REMLf, 1/2 * determinant(crossprod(X), logarithm=TRUE)$modulus, 0) +
                     -1/2 * sum(log(vi + tau2)) - 1/2 * determinant(crossprod(X,W) %*% X, logarithm=TRUE)$modulus - 1/2 * RSS
         }

      }

   }

   if (verbose) {
      cat(mstyle$verbose(paste0("ll = ",          ifelse(is.na(llval), NA, formatC(llval, digits=digits[["fit"]], format="f", flag=" ")), "  ")))
      cat(mstyle$verbose(paste0("alpha = ", paste(ifelse(is.na(alpha), NA, formatC(alpha, digits=digits[["est"]], format="f", flag=" ")), collapse=" "))))
      cat("\n")
   }

   return(-1 * llval)

}

############################################################################
