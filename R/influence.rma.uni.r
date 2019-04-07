influence.rma.uni <- function(model, digits, progbar=FALSE, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(model, "rma.uni"))
      stop(mstyle$stop("Argument 'model' must be an object of class \"rma.uni\"."))

   if (inherits(model, "rma.ls"))
      stop(mstyle$stop("Method not available for objects of class \"rma.ls\"."))

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   x <- model

   if (x$k == 1)
      stop(mstyle$stop("Stopped because k = 1."))

   ddd <- list(...)

   btt <- .set.btt(ddd$btt, x$p, int.incl=FALSE)
   m <- length(btt)

   if (is.null(ddd$measure)) {
      measure <- "all"
   } else {
      measure <- ddd$measure
   }

   if (missing(digits)) {
      digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
   }

   if (!measure == "cooks.distance" && inherits(model, "robust.rma"))
      stop(mstyle$stop("Method not available for objects of class \"robust.rma\"."))

   #########################################################################

   tau2.del <- rep(NA_real_, x$k)
   delpred  <- rep(NA_real_, x$k)
   vdelpred <- rep(NA_real_, x$k)
   s2w.del  <- rep(NA_real_, x$k)
   QE.del   <- rep(NA_real_, x$k)
   dffits   <- rep(NA_real_, x$k)
   dfbs     <- matrix(NA_real_, nrow=x$k, ncol=x$p)
   cook.d   <- rep(NA_real_, x$k)
   cov.r    <- rep(NA_real_, x$k)
   weight   <- rep(NA_real_, x$k)

   ### predicted values under the full model

   pred.full <- x$X %*% x$beta

   ### calculate inverse of variance-covariance matrix under the full model (needed for the Cook's distances)

   svb <- chol2inv(chol(x$vb[btt,btt,drop=FALSE]))

   ### also need stXAX/stXX and H matrix for DFFITS calculation when not using the standard weights

   if (x$weighted) {
      if (!is.null(x$weights)) {
         A <- diag(x$weights, nrow=x$k, ncol=x$k)
         stXAX <- .invcalc(X=x$X, W=A, k=x$k)
         H <- x$X %*% stXAX %*% t(x$X) %*% A
      }
   } else {
      stXX <- .invcalc(X=x$X, W=diag(x$k), k=x$k)
      H <- x$X %*% stXX %*% t(x$X)
   }

   ### hat values

   options(na.action = "na.omit")
   hat <- hatvalues(x)
   options(na.action = na.act)

   ### note: skipping NA cases
   ### also: it is possible that model fitting fails, so that generates more NAs (these NAs will always be shown in output)

   if (progbar)
      pbar <- txtProgressBar(min=0, max=x$k, style=3)

   for (i in seq_len(x$k)) {

      if (progbar)
         setTxtProgressBar(pbar, i)

      res <- try(suppressWarnings(rma.uni(x$yi, x$vi, weights=x$weights, mods=x$X, intercept=FALSE, method=x$method, weighted=x$weighted, test=x$test, tau2=ifelse(x$tau2.fix, x$tau2, NA), control=x$control, subset=-i)), silent=TRUE)

      if (inherits(res, "try-error"))
         next

      ### removing an observation could lead to a model coefficient becoming inestimable

      if (any(res$coef.na))
         next

      ### save tau2.del and QE.del values

      tau2.del[i] <- res$tau2
      QE.del[i]   <- res$QE

      ### 'deleted' predicted value for the ith observation based on the model without the ith observation included

      Xi          <- matrix(x$X[i,], nrow=1)
      delpred[i]  <- Xi %*% res$beta
      vdelpred[i] <- Xi %*% tcrossprod(res$vb,Xi)
      s2w.del[i]  <- res$s2w

      ### compute DFFITS

      if (x$weighted) {
         if (is.null(x$weights)) {
            dffits[i] <- (pred.full[i] - delpred[i]) / sqrt(res$s2w * hat[i] * (tau2.del[i] + x$vi[i]))
         } else {
            dffits[i] <- (pred.full[i] - delpred[i]) / sqrt(res$s2w * diag(H %*% diag(tau2.del[i] + x$vi, nrow=x$k, ncol=x$k) %*% t(H)))[i]
         }
      } else {
         dffits[i] <- (pred.full[i] - delpred[i]) / sqrt(res$s2w * diag(H %*% diag(tau2.del[i] + x$vi, nrow=x$k, ncol=x$k) %*% t(H)))[i]
      }

      #dffits[i]  <- (pred.full[i] - delpred[i]) / sqrt(vdelpred[i])

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

      ### compute DFBETA and DFBETAS

      dfb <- x$beta - res$beta
      dfbs[i,] <- dfb / sqrt(res$s2w * diag(vb.del))

      #dfbs[i,] <- dfb / sqrt(diag(res$vb))

      ### compute DFBETA (including coefficients as specified via btt)

      dfb <- x$beta[btt] - res$beta[btt]

      ### compute Cook's distance

      cook.d[i] <- crossprod(dfb,svb) %*% dfb # / x$p

      #cook.d[i] <- sum(1/(x$vi+tau2.del[i]) * (pred.full - x$X %*% res$beta)^2) / x$p

      ### compute covariance ratio

      cov.r[i] <- det(res$vb[btt,btt,drop=FALSE]) / det(x$vb[btt,btt,drop=FALSE])

   }

   if (progbar)
      close(pbar)

   ### calculate studentized residual

   resid <- x$yi - delpred
   resid[abs(resid) < 100 * .Machine$double.eps] <- 0
   #resid[abs(resid) < 100 * .Machine$double.eps * median(abs(resid), na.rm=TRUE)] <- 0 ### see lm.influence
   #seresid <- sqrt(x$vi + vdelpred + tau2.del)
   seresid <- sqrt(x$vi * s2w.del + vdelpred + tau2.del * s2w.del) ### this yields the same results as a mean shift outlier model when using test="knha"
   stresid <- resid / seresid

   ### extract weights

   options(na.action="na.omit")
   weight <- weights(x)
   options(na.action = na.act)

   #########################################################################

   inf <- matrix(NA_real_, nrow=x$k.f, ncol=8)
   inf[x$not.na,] <- cbind(stresid, dffits, cook.d, cov.r, tau2.del, QE.del, hat, weight)
   colnames(inf) <- c("rstudent", "dffits", "cook.d", "cov.r", "tau2.del", "QE.del", "hat", "weight")
   inf <- data.frame(inf)

   tmp <- dfbs
   dfbs <- matrix(NA_real_, nrow=x$k.f, ncol=x$p)
   dfbs[x$not.na,] <- tmp
   colnames(dfbs) <- rownames(x$beta)
   dfbs <- data.frame(dfbs)

   #########################################################################

   ### determine "influential" cases

   is.infl <-
      #abs(inf$rstudent) > qnorm(.975) |
      abs(inf$dffits) > 3*sqrt(x$p/(x$k-x$p)) |
      pchisq(inf$cook.d, df=m) > .50 |
      #inf$cov.r > 1 + 3*m/(x$k-m) |
      #inf$cov.r < 1 - 3*m/(x$k-m) |
      inf$hat > 3*x$p/x$k |
      apply(abs(dfbs) > 1, 1, any) ### consider using rowAnys() from matrixStats package

   #print(is.infl)

   #########################################################################

   if (na.act == "na.omit") {

      out <- list(rstudent=inf$rstudent[x$not.na], dffits=inf$dffits[x$not.na], cook.d=inf$cook.d[x$not.na], cov.r=inf$cov.r[x$not.na], tau2.del=inf$tau2.del[x$not.na], QE.del=inf$QE.del[x$not.na], hat=inf$hat[x$not.na], weight=inf$weight[x$not.na], inf=ifelse(is.infl & !is.na(is.infl), "*", "")[x$not.na], slab=x$slab[x$not.na], digits=digits)

      out <- list(inf=out)
      out$dfbs <- lapply(dfbs, function(z) z[x$not.na])
      out$dfbs <- c(out$dfbs, list(slab=x$slab[x$not.na], digits=digits))

      out <- c(out, list(ids=x$ids[x$not.na], not.na=x$not.na[x$not.na], is.infl=is.infl[x$not.na]))

   }

   if (na.act == "na.exclude" || na.act == "na.pass") {

      out <- list(rstudent=inf$rstudent, dffits=inf$dffits, cook.d=inf$cook.d, cov.r=inf$cov.r, tau2.del=inf$tau2.del, QE.del=inf$QE.del, hat=inf$hat, weight=inf$weight, inf=ifelse(is.infl & !is.na(is.infl), "*", ""), slab=x$slab, digits=digits)

      out <- list(inf=out)
      out$dfbs <- lapply(dfbs, function(z) z)
      out$dfbs <- c(out$dfbs, list(slab=x$slab, digits=digits))

      out <- c(out, list(ids=x$ids, not.na=x$not.na, is.infl=is.infl))

   }

   out <- c(out, list(tau2=x$tau2, QE=x$QE, k=x$k, p=x$p, m=m, digits=digits))

   if (na.act == "na.fail" && any(!x$not.na))
      stop(mstyle$stop("Missing values in results."))

   class(out$inf)  <- c("list.rma")
   class(out$dfbs) <- c("list.rma")
   class(out) <- c("infl.rma.uni")

   if (measure == "cooks.distance") {
      names(out$inf$cook.d) <- out$inf$slab
      out <- out$inf$cook.d
   }

   if (measure == "dfbetas")
      out <- out$dfbs

   if (measure == "rstudent") {

      if (na.act == "na.omit") {
         resid.f   <- resid
         seresid.f <- seresid
         stresid.f <- stresid
      }

      if (na.act == "na.exclude" || na.act == "na.pass") {
         resid.f   <- rep(NA_real_, x$k.f)
         seresid.f <- rep(NA_real_, x$k.f)
         stresid.f <- rep(NA_real_, x$k.f)
         resid.f[x$not.na]   <- resid
         seresid.f[x$not.na] <- seresid
         stresid.f[x$not.na] <- stresid
      }

      out <- list(resid=resid.f, se=seresid.f, z=stresid.f, slab=out$inf$slab, digits=digits)

      class(out) <- c("list.rma")

   }

   return(out)

}
