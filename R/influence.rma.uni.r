### Note: The definitions used for dffits, dfbetas, and cook.d below give the same results as
### influence.measures(lm(...)) when all vi=0 (except for cook.d which is not scaled by 1/p).

influence.rma.uni <- function(model, digits, ...) {

   if (!is.element("rma.uni", class(model)))
      stop("Argument 'model' must be an object of class \"rma.uni\".")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop("Unknown 'na.action' specified under options().")

   x <- model

   if (missing(digits))
      digits <- x$digits

   if (x$k == 1)
      stop("Stopped because k = 1.")

   #########################################################################

   tau2.del <- rep(NA_real_, x$k.f)
   delpred  <- rep(NA_real_, x$k.f)
   vdelpred <- rep(NA_real_, x$k.f)
   QE.del   <- rep(NA_real_, x$k.f)
   dffits   <- rep(NA_real_, x$k.f)
   dfbs     <- matrix(NA_real_, nrow=x$k.f, ncol=length(x$b))
   cook.d   <- rep(NA_real_, x$k.f)
   cov.r    <- rep(NA_real_, x$k.f)
   weight   <- rep(NA_real_, x$k.f)

   ### predicted values under the full model

   pred.full <- x$X.f %*% x$b

   ### calculate inverse of variance-covariance matrix under the full model (needed for the Cook's distances)
   ### also need H matrix for dffits calculation (when not using the standard weights)

   if (x$weighted) {
      if (is.null(x$weights)) {
         W   <- diag(1/(x$vi + x$tau2), nrow=x$k, ncol=x$k)
         svb <- crossprod(x$X,W) %*% x$X / x$s2w
      } else {
         svb   <- chol2inv(chol(x$vb))
         A     <- diag(x$weights, nrow=x$k, ncol=x$k)
         stXAX <- .invcalc(X=x$X, W=A, k=x$k)
         H     <- x$X %*% stXAX %*% t(x$X) %*% A
      }
   } else {
      svb  <- chol2inv(chol(x$vb))
      stXX <- .invcalc(X=x$X, W=diag(x$k), k=x$k)
      H    <- x$X %*% stXX %*% t(x$X)
   }

   ### hat values

   options(na.action = "na.pass")
   hii <- hatvalues(x)
   options(na.action = na.act)

   ### note: skipping NA cases
   ### also: it is possible that model fitting fails, so that generates more NAs (these NAs will always be shown in output)

   for (i in seq_len(x$k.f)[x$not.na]) {

      res <- try(suppressWarnings(rma(x$yi.f[-i], x$vi.f[-i], weights=x$weights.f[-i], mods=cbind(x$X.f[-i,]), method=x$method, weighted=x$weighted, intercept=FALSE, knha=x$knha, control=x$control)), silent=TRUE)

      if (inherits(res, "try-error"))
         next

      ### removing an observation could lead to a model coefficient becoming inestimable

      if (any(res$coef.na))
         next

      ### save tau2.del and QE.del values

      tau2.del[i] <- res$tau2
      QE.del[i]   <- res$QE

      ### 'deleted' predicted value for the ith observation based on the model without the ith observation included

      Xi          <- matrix(x$X.f[i,], nrow=1)
      delpred[i]  <- Xi %*% res$b
      vdelpred[i] <- Xi %*% tcrossprod(res$vb,Xi)

      ### compute dffits

      #dffits[i]  <- (pred.full[i] - delpred[i]) / sqrt(vdelpred[i])
      if (x$weighted) {
         if (is.null(x$weights)) {
            dffits[i] <- (pred.full[i] - delpred[i]) / sqrt(res$s2w * hii[i] * (tau2.del[i] + x$vi.f[i]))
         } else {
            dffits[i] <- (pred.full[i] - delpred[i]) / sqrt(res$s2w * diag(H %*% diag(tau2.del[i] + x$vi, nrow=x$k, ncol=x$k) %*% t(H)))[i-sum(!x$not.na[1:i])]
         }
      } else {
         dffits[i] <- (pred.full[i] - delpred[i]) / sqrt(res$s2w * diag(H %*% diag(tau2.del[i] + x$vi, nrow=x$k, ncol=x$k) %*% t(H)))[i-sum(!x$not.na[1:i])]
      }

      ### compute var-cov matrix of the fixed effects for the full model, but with tau2.del[i] plugged in

      if (x$weighted) {
         if (is.null(x$weights)) {
            vb.del <- .invcalc(X=x$X, W=diag(1/(x$vi+tau2.del[i]), nrow=x$k, ncol=x$k), k=x$k)
         } else {
            vb.del <- tcrossprod(stXAX,x$X) %*% A %*% diag(x$vi+tau2.del[i], nrow=x$k, ncol=x$k) %*% A %*% x$X %*% stXAX
         }
      } else {
         vb.del <- tcrossprod(stXX,x$X) %*% diag(x$vi+tau2.del[i], nrow=x$k, ncol=x$k) %*% x$X %*% stXX
      }

      ### compute dbeta and dfbetas value(s)

      dfb <- x$b - res$b
      dfbs[i,] <- dfb / sqrt(res$s2w * diag(vb.del))
      #dfbs[i,] <- dfb / sqrt(diag(res$vb))

      ### compute Cook's distance

      cook.d[i]  <- crossprod(dfb,svb) %*% dfb # / x$p
      #cook.d[i] <- sum(1/(x$vi.f+tau2.del[i]) * (pred.full - x$X.f %*% res$b)^2) / x$p

      ### compute covariance ratio

      cov.r[i]   <- det(res$vb) / det(x$vb)

   }

   ### calculate studentized residual

   delresid   <- x$yi.f - delpred
   sedelresid <- sqrt(x$vi.f + vdelpred + tau2.del)
   standelres <- delresid / sedelresid

   ### extract weights

   options(na.action="na.omit")
   weight[x$not.na] <- weights(x)
   options(na.action = na.act)

   #########################################################################

   inf <- cbind(standelres, dffits, cook.d, cov.r, tau2.del, QE.del, hii, weight)
   dfbs <- cbind(dfbs)

   inf <- data.frame(inf)
   dfbs <- data.frame(dfbs)

   #########################################################################

   ### determine "influential" cases

   is.infl <-
      #abs(inf$standelres) > qnorm(.975) |
      abs(inf$dffits) > 3*sqrt(x$p/(x$k-x$p)) |
      pchisq(inf$cook.d, df=x$p) > .50 |
      #inf$cov.r > 1 + 3*x$p/(x$k-x$p) |
      #inf$cov.r < 1 - 3*x$p/(x$k-x$p) |
      inf$hii > 3*x$p/x$k |
      apply(abs(dfbs) > 1, 1, any) ### consider using rowAnys() from matrixStats package

   #print(ids.infl)

   #########################################################################

   out <- list(inf=inf, dfbs=dfbs, tau2=x$tau2, QE=x$QE, ids=x$ids, not.na=x$not.na, is.infl=is.infl, k=x$k, p=x$p, digits=digits)
   rownames(out$inf) <- x$slab
   rownames(out$dfbs) <- x$slab

   colnames(out$dfbs) <- rownames(x$b)
   colnames(out$inf) <- c("rstudent", "dffits", "cook.d", "cov.r", "tau2.del", "QE.del", "hat", "weight")

   class(out) <- "infl.rma.uni"
   return(out)

}
