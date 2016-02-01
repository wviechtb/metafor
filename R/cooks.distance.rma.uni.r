cooks.distance.rma.uni <- function(model, ...) {

   if (!inherits(model, "rma.uni"))
      stop("Argument 'model' must be an object of class \"rma.uni\".")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop("Unknown 'na.action' specified under options().")

   x <- model

   if (x$k == 1)
      stop("Stopped because k = 1.")

   #########################################################################

   cook.d <- rep(NA_real_, x$k.f)

   ### calculate inverse of variance-covariance matrix under the full model (needed for the Cook's distances)

   svb <- chol2inv(chol(x$vb))

   ### note: skipping NA cases
   ### also: it is possible that model fitting fails, so that generates more NAs (these NAs will always be shown in output)

   for (i in seq_len(x$k.f)[x$not.na]) {

      res <- try(suppressWarnings(rma.uni(x$yi.f, x$vi.f, weights=x$weights.f, mods=x$X.f, intercept=FALSE, method=x$method, weighted=x$weighted, knha=x$knha, tau2=ifelse(x$tau2.fix, x$tau2, NA), control=x$control, subset=-i)), silent=TRUE)

      if (inherits(res, "try-error"))
         next

      ### removing an observation could lead to a model coefficient becoming inestimable

      if (any(res$coef.na))
         next

      ### compute dfbeta value(s)

      dfb <- x$b - res$b

      ### compute Cook's distance

      cook.d[i]  <- crossprod(dfb,svb) %*% dfb

   }

   #########################################################################

   if (na.act == "na.omit") {
      out <- cook.d[x$not.na]
      names(out) <- x$slab[x$not.na]
   }

   if (na.act == "na.exclude" || na.act == "na.pass") {
      out <- cook.d
      names(out) <- x$slab
   }

   if (na.act == "na.fail" && any(!x$not.na))
      stop("Missing values in results.")

   return(out)

}
