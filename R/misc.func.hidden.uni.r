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
      cat(mstyle$verbose(paste("tau2 =", fmtx(tau2val, digits[["var"]], addwidth=4), "  objective =", fmtx(res, digits[["var"]], flag=" "), "\n")))

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
   args <- list(yi=obj$yi - c(val*obj$X[,j]), vi=obj$vi, weights=obj$weights, mods=obj$X, intercept=FALSE, method=obj$method, weighted=obj$weighted, test=obj$test, tau2=ifelse(obj$tau2.fix, obj$tau2, NA), control=obj$control, skipr2=TRUE)
   res <- try(suppressWarnings(.do.call(rma.uni, args)), silent=TRUE)

   if (inherits(res, "try-error"))
      stop()

   ### p-value based on permutation test
   pval <- permutest(res, exact=exact, iter=iter, progbar=FALSE, control=control)$pval[j]

   ### get difference between p-value and level
   diff <- pval - level / ifelse(control$alternative == "two.sided", 1, 2)

   ### show progress
   if (progbar)
      cat(mstyle$verbose(paste("pval =", fmtx(pval, digits[["pval"]]), " diff =", fmtx(diff, digits[["pval"]], flag=" "), " val =", fmtx(val, digits[["est"]], flag=" "), "\n")))

   ### penalize negative differences, which should force the CI bound to correspond to a p-value of *at least* level
   diff <- ifelse(diff < 0, diff*10, diff)

   return(diff)

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

### -1 times the log likelihood (regular or restricted) for location-scale model

.ll.rma.ls <- function(par, yi, vi, X, Z, reml, k, pX,
                       alpha.val, beta.val, verbose, digits,
                       REMLf, link, mZ, alpha.min, alpha.max, alpha.transf,
                       tau2.min, tau2.max, optbeta) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (optbeta) {
      beta  <- par[seq_len(pX)]
      beta  <- ifelse(is.na(beta.val), beta, beta.val)
      alpha <- par[-seq_len(pX)]
   } else {
      alpha <- par
   }

   if (alpha.transf)
      alpha <- mapply(.mapfun.alpha, alpha, alpha.min, alpha.max)

   alpha <- ifelse(is.na(alpha.val), alpha, alpha.val)

   ### compute predicted tau2 values

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

         ### compute weights / weights matrix
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

      ### compute residual sum of squares
      RSS <- sum(wi*c(yi - X %*% beta)^2)

      ### compute log-likelihood
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

.rma.ls.ineqfun.pos <- function(par, yi, vi, X, Z, reml, k, pX, alpha.val, beta.val, verbose, digits, REMLf, link, mZ, alpha.min, alpha.max, alpha.transf, tau2.min, tau2.max, optbeta) {

   if (optbeta) {
      alpha <- par[-seq_len(pX)]
   } else {
      alpha <- par
   }

   if (alpha.transf)
      alpha <- mapply(.mapfun.alpha, alpha, alpha.min, alpha.max)

   alpha <- ifelse(is.na(alpha.val), alpha, alpha.val)

   tau2 <- c(Z %*% alpha)

   return(tau2)

}

.rma.ls.ineqfun.neg <- function(par, yi, vi, X, Z, reml, k, pX, alpha.val, beta.val, verbose, digits, REMLf, link, mZ, alpha.min, alpha.max, alpha.transf, tau2.min, tau2.max, optbeta) {

   if (optbeta) {
      alpha <- par[-seq_len(pX)]
   } else {
      alpha <- par
   }

   if (alpha.transf)
      alpha <- mapply(.mapfun.alpha, alpha, alpha.min, alpha.max)

   alpha <- ifelse(is.na(alpha.val), alpha, alpha.val)

   tau2 <- -c(Z %*% alpha)

   return(tau2)

}

############################################################################
