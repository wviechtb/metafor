rstudent.rma.uni <- function(model, digits, ...) {

   if (!inherits(model, "rma.uni"))
      stop("Argument 'model' must be an object of class \"rma.uni\".")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop("Unknown 'na.action' specified under options().")

   x <- model

   if (missing(digits))
      digits <- x$digits

   #########################################################################

   tau2.del <- rep(NA_real_, x$k.f)
   delpred  <- rep(NA_real_, x$k.f)
   vdelpred <- rep(NA_real_, x$k.f)

   ### note: skipping NA cases
   ### also: it is possible that model fitting fails, so that generates more NAs (these NAs will always be shown in output)

   for (i in seq_len(x$k.f)[x$not.na]) {

      res <- try(suppressWarnings(rma.uni(x$yi.f, x$vi.f, weights=x$weights.f, mods=x$X.f, intercept=FALSE, method=x$method, weighted=x$weighted, test=x$test, tau2=ifelse(x$tau2.fix, x$tau2, NA), control=x$control, subset=-i)), silent=TRUE)

      if (inherits(res, "try-error"))
         next

      ### removing an observation could lead to a model coefficient becoming inestimable

      if (any(res$coef.na))
         next

      tau2.del[i] <- res$tau2
      Xi          <- matrix(x$X.f[i,], nrow=1)
      delpred[i]  <- Xi %*% res$b
      vdelpred[i] <- Xi %*% tcrossprod(res$vb,Xi)

   }

   delresid <- x$yi.f - delpred
   delresid[abs(delresid) < 100 * .Machine$double.eps] <- 0
   #delresid[abs(delresid) < 100 * .Machine$double.eps * median(abs(delresid), na.rm=TRUE)] <- 0 ### see lm.influence
   sedelresid <- sqrt(x$vi.f + vdelpred + tau2.del)
   standelres <- delresid / sedelresid

   #########################################################################

   if (na.act == "na.omit") {
      out <- list(resid=delresid[x$not.na], se=sedelresid[x$not.na], z=standelres[x$not.na])
      out$slab <- x$slab[x$not.na]
   }

   if (na.act == "na.exclude" || na.act == "na.pass") {
      out <- list(resid=delresid, se=sedelresid, z=standelres)
      out$slab <- x$slab
   }

   if (na.act == "na.fail" && any(!x$not.na))
      stop("Missing values in results.")

   out$digits <- digits

   class(out) <- "list.rma"
   return(out)

}
