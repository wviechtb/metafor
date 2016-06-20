plot.profile.rma <- function(x, ylim, pch=19, ...) {

   #########################################################################

   if (!inherits(x, "profile.rma"))
      stop("Argument 'x' must be an object of class \"profile.rma\".")

   if (dev.cur() == 1) {
      par(mfrow=c(x$comps, 1))
      on.exit(par(mfrow=c(1,1)))
   }

   if (missing(ylim)) {
      missing.ylim <- TRUE
   } else {
      missing.ylim <- FALSE
   }

   if (x$comps == 1) {

      if (missing.ylim)
         ylim <- x$ylim

      plot(x[[1]], x[[2]], type="o", xlab=x$xlab, ylab=paste(ifelse(x$method=="REML", "Restricted", ""), " Log-Likelihood", sep=""), main=x$title, bty="l", pch=pch, ylim=ylim, ...)
      abline(v=x$vc, lty="dotted")
      abline(h=x$maxll, lty="dotted")
      #abline(h=max(lls, na.rm=TRUE), lty="dotted")

   } else {

      for (j in 1:x$comps) {

         if (missing.ylim)
            ylim <- x[[j]]$ylim

         plot(x[[j]], ylim=ylim, pch=pch, ...)

      }

   }

}
