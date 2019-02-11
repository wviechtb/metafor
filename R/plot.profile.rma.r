plot.profile.rma <- function(x, xlim, ylim, pch=19, xlab, ylab, main, cline=FALSE, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(x, "profile.rma"))
      stop(mstyle$stop("Argument 'x' must be an object of class \"profile.rma\"."))

   if (dev.cur() == 1) {
      par(mfrow=c(x$comps, 1))
      #on.exit(par(mfrow=c(1,1)))
   }

   missing.xlim <- missing(xlim)
   missing.ylim <- missing(ylim)
   missing.xlab <- missing(xlab)
   missing.ylab <- missing(ylab)
   missing.main <- missing(main)

   #########################################################################

   if (x$comps == 1) {

      if (missing.xlim)
         xlim <- x$xlim

      if (missing.ylim)
         ylim <- x$ylim

      if (missing.xlab)
         xlab <- x$xlab

      if (missing.ylab)
         ylab <- paste(ifelse(x$method=="REML", "Restricted ", ""), "Log-Likelihood", sep="")

      if (missing.main)
         main <- x$title

      plot(x[[1]], x[[2]], type="o", xlab=xlab, ylab=ylab, main=main, bty="l", pch=pch, xlim=xlim, ylim=ylim, ...)
      abline(v=x$vc, lty="dotted")
      abline(h=x$maxll, lty="dotted")

      if (cline)
         abline(h=x$maxll - qchisq(0.95, df=1)/2, lty="dotted")

   } else {

      for (j in seq_len(x$comps)) {

         if (missing.xlim)
            xlim <- x[[j]]$xlim

         if (missing.ylim)
            ylim <- x[[j]]$ylim

         if (missing.xlab) {
            xlab <- x[[j]]$xlab
         } else {
            if (length(xlab) == 1) {
               xlab <- rep(xlab, x$comps)
            }
         }

         if (missing.ylab) {
            ylab <- paste(ifelse(x[[j]]$method=="REML", "Restricted ", ""), "Log-Likelihood", sep="")
         } else {
            if (length(ylab) == 1) {
               ylab <- rep(ylab, x$comps)
            }
         }

         if (missing.main) {
            main <- x[[j]]$title
         } else {
            if (length(main) == 1) {
               main <- rep(main, x$comps)
            }
         }

         plot(x[[j]], xlim=xlim, ylim=ylim, main=if (missing.main) main else main[j], pch=pch, xlab=if (missing.xlab) xlab else xlab[j], ylab=if (missing.ylab) ylab else ylab[j], cline=cline, ...)

      }

   }

}
