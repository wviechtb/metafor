plot.profile.rma <- function(x, ylim, pch=19, ylab, cline=FALSE, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(x, "profile.rma"))
      stop(mstyle$stop("Argument 'x' must be an object of class \"profile.rma\"."))

   if (dev.cur() == 1) {
      par(mfrow=c(x$comps, 1))
      #on.exit(par(mfrow=c(1,1)))
   }

   if (missing(ylim)) {
      missing.ylim <- TRUE
   } else {
      missing.ylim <- FALSE
   }

   if (missing(ylab)) {
      missing.ylab <- TRUE
   } else {
      missing.ylab <- FALSE
   }

   #########################################################################

   if (x$comps == 1) {

      if (missing.ylim)
         ylim <- x$ylim

      if (missing.ylab)
         ylab <- paste(ifelse(x$method=="REML", "Restricted ", ""), "Log-Likelihood", sep="")

      plot(x[[1]], x[[2]], type="o", xlab=x$xlab, ylab=ylab, main=x$title, bty="l", pch=pch, ylim=ylim, ...)
      abline(v=x$vc, lty="dotted")
      abline(h=x$maxll, lty="dotted")

      if (cline)
         abline(h=x$maxll - qchisq(0.95, df=1)/2, lty="dotted")

   } else {

      for (j in seq_len(x$comps)) {

         if (missing.ylim)
            ylim <- x[[j]]$ylim

      if (missing.ylab)
         ylab <- paste(ifelse(x[[j]]$method=="REML", "Restricted ", ""), "Log-Likelihood", sep="")

         plot(x[[j]], ylim=ylim, pch=pch, ylab=ylab, cline=cline, ...)

      }

   }

}
