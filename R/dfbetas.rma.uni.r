dfbetas.rma.uni <- function(model, ...) {

   if (!inherits(model, "rma.uni"))
      stop("Argument 'model' must be an object of class \"rma.uni\".")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop("Unknown 'na.action' specified under options().")

   x <- model

   if (x$k == 1)
      stop("Stopped because k = 1.")

   #########################################################################

   tau2.del <- rep(NA_real_, x$k.f)
   dfbs <- matrix(NA_real_, nrow=x$k.f, ncol=x$p)

   if (x$weighted && !is.null(x$weights)) {
      A     <- diag(x$weights, nrow=x$k, ncol=x$k)
      stXAX <- .invcalc(X=x$X, W=A, k=x$k)
   } else {
      stXX <- .invcalc(X=x$X, W=diag(x$k), k=x$k)
   }

   ### note: skipping NA cases
   ### also: it is possible that model fitting fails, so that generates more NAs (these NAs will always be shown in output)

   for (i in seq_len(x$k.f)[x$not.na]) {

      res <- try(suppressWarnings(rma.uni(x$yi.f, x$vi.f, weights=x$weights.f, mods=x$X.f, intercept=FALSE, method=x$method, weighted=x$weighted, test=x$test, tau2=ifelse(x$tau2.fix, x$tau2, NA), control=x$control, subset=-i)), silent=TRUE)

      if (inherits(res, "try-error"))
         next

      ### removing an observation could lead to a model coefficient becoming inestimable

      if (any(res$coef.na))
         next

      ### save tau2.del

      tau2.del[i] <- res$tau2

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

   }

   #########################################################################

   if (na.act == "na.omit") {
      out <- dfbs[x$not.na,,drop=FALSE]
      rownames(out) <- x$slab[x$not.na]
   }

   if (na.act == "na.exclude" || na.act == "na.pass") {
      out <- dfbs
      rownames(out) <- x$slab
   }

   if (na.act == "na.fail" && any(!x$not.na))
      stop("Missing values in results.")

   colnames(out) <- rownames(x$b)

   out <- data.frame(out)

   return(out)

}
