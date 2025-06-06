plot.rma.peto <- function(x, qqplot=FALSE, ...) {

   #########################################################################

   mstyle <- .get.mstyle()

   .chkclass(class(x), must="rma.peto")

   na.act <- getOption("na.action")
   on.exit(options(na.action=na.act), add=TRUE)

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   .start.plot()

   # if no plotting device is open or mfrow is too small, set mfrow appropriately
   if (dev.cur() == 1L || prod(par("mfrow")) < 4L)
      par(mfrow=n2mfrow(4))
   on.exit(par(mfrow=c(1L,1L)), add=TRUE)

   bg <- .coladj(par("bg","fg"), dark=0.35, light=-0.35)
   col.na <- .coladj(par("bg","fg"), dark=0.2, light=-0.2)

   #########################################################################

   forest(x, ...)
   title("Forest Plot", ...)

   #########################################################################

   funnel(x, ...)
   title("Funnel Plot", ...)

   #########################################################################

   radial(x, ...)
   title("Radial Plot", ...)

   #########################################################################

   if (qqplot) {

      qqnorm(x, ...)

   } else {

      options(na.action = "na.pass")
      z <- rstandard(x)$z
      options(na.action = na.act)

      not.na <- !is.na(z)

      if (na.act == "na.omit") {
         z   <- z[not.na]
         ids <- x$ids[not.na]
         not.na <- not.na[not.na]
      }

      if (na.act == "na.exclude" || na.act == "na.pass")
         ids <- x$ids

      k <- length(z)

      plot(NA, NA, xlim=c(1,k), ylim=c(min(z, -2, na.rm=TRUE), max(z, 2, na.rm=TRUE)), xaxt="n", xlab="Study", ylab="", bty="l", ...)
      lines(seq_len(k)[not.na], z[not.na], col=col.na, ...)
      lines(seq_len(k), z, ...)
      points(x=seq_len(k), y=z, pch=21, bg=bg, ...)
      axis(side=1, at=seq_len(k), labels=ids, ...)
      abline(h=0, lty="dashed", ...)
      abline(h=c(qnorm(0.025),qnorm(0.975)), lty="dotted", ...)

      title("Standardized Residuals", ...)

   }

   #########################################################################

   invisible()

}
