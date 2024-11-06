plot.profile.rma <- function(x, xlim, ylim, pch=19, xlab, ylab, main, refline=TRUE, cline=FALSE, ...) {

   #########################################################################

   mstyle <- .get.mstyle()

   .chkclass(class(x), must="profile.rma")

   .start.plot()

   if (x$comps > 1L) {
      # if no plotting device is open or mfrow is too small, set mfrow appropriately
      if (dev.cur() == 1L || prod(par("mfrow")) < x$comps)
         par(mfrow=n2mfrow(x$comps))
      on.exit(par(mfrow=c(1L,1L)), add=TRUE)
   }

   missing.xlim <- missing(xlim)
   missing.ylim <- missing(ylim)
   missing.xlab <- missing(xlab)
   missing.ylab <- missing(ylab)
   missing.main <- missing(main)

   lplot   <- function(..., time, LB, startmethod, sub1,      sqrt, exp, pred, blup, callback) plot(...)
   lpoints <- function(..., time, LB, startmethod, sub1, log, sqrt, exp, pred, blup, callback) points(...) # need 'log' here so profile(res, log="x") doesn't throw a warning

   #########################################################################

   if (x$comps == 1) {

      if (missing.xlim)
         xlim <- x$xlim

      if (missing.ylim)
         ylim <- x$ylim

      if (missing.xlab)
         xlab <- x$xlab

      if (missing.ylab) {
         if (isTRUE(x$exp)) {
            ylab <- paste0(ifelse(x$method=="REML", "Restricted ", ""), "Likelihood")
         } else {
            ylab <- paste0(ifelse(x$method=="REML", "Restricted ", ""), "Log-Likelihood")
         }
      }

      if (missing.main)
         main <- x$title

      ### add the actual vc value to the profile

      if (min(x[[1]]) <= x$vc && max(x[[1]]) >= x$vc) {
         pos <- which(x[[1]] >= x$vc)[1]
         x[[1]] <- c(x[[1]][seq_len(pos-1)], x$vc,    x[[1]][pos:length(x[[1]])])
         x[[2]] <- c(x[[2]][seq_len(pos-1)], x$maxll, x[[2]][pos:length(x[[2]])])
      }

      lplot(x[[1]], x[[2]], type="n", xlab=xlab, ylab=ylab, main=main, bty="l", xlim=xlim, ylim=ylim, ...)

      if (refline) {
         abline(v=x$vc, lty="dotted")
         abline(h=x$maxll, lty="dotted")
      }

      if (isTRUE(cline))
         cline <- 0.05

      if (is.numeric(cline)) {
         cline <- .level(cline, argname="cline")
         if (isTRUE(x$exp)) {
            hval <- exp(log(x$maxll) - qchisq(1-cline, df=1)/2)
         } else {
            hval <- x$maxll - qchisq(1-cline, df=1)/2
         }
         abline(h=hval, lty="dotted")
      }

      lpoints(x[[1]], x[[2]], type="o", pch=pch, ...)

   } else {

      for (j in seq_len(x$comps)) {

         if (missing.xlim)
            xlim <- x[[j]]$xlim

         if (missing.ylim)
            ylim <- x[[j]]$ylim

         if (missing.xlab) {
            xlab <- x[[j]]$xlab
         } else {
            xlab <- .expand1(xlab, x$comps)
         }

         if (missing.ylab) {
            if (isTRUE(x$exp)) {
               ylab <- paste0(ifelse(x$method=="REML", "Restricted ", ""), "Likelihood")
            } else {
               ylab <- paste0(ifelse(x$method=="REML", "Restricted ", ""), "Log-Likelihood")
            }
         } else {
            ylab <- .expand1(ylab, x$comps)
         }

         if (missing.main) {
            main <- x[[j]]$title
         } else {
            main <- .expand1(main, x$comps)
         }

         lplot(x[[j]], xlim=xlim, ylim=ylim, pch=pch,
               xlab=if (missing.xlab) xlab else xlab[j],
               ylab=if (missing.ylab) ylab else ylab[j],
               main=if (missing.main) main else main[j],
               cline=cline, ...)

      }

   }

}
