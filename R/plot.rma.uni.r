plot.rma.uni <- function(x, qqplot=FALSE, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="rma.uni", notav=c("robust.rma", "rma.ls", "rma.nn", "rma.uni.selmodel"))

   na.act <- getOption("na.action")
   on.exit(options(na.action=na.act), add=TRUE)

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   par.mfrow <- par("mfrow")
   par(mfrow=c(2,2))
   on.exit(par(mfrow = par.mfrow), add=TRUE)

   #########################################################################

   if (x$int.only) {

      ######################################################################

      forest(x, ...)
      title("Forest Plot", ...)

      ######################################################################

      funnel(x, ...)
      title("Funnel Plot", ...)

      ######################################################################

      radial(x, ...)
      title("Radial Plot", ...)

      ######################################################################

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
         lines(seq_len(k)[not.na], z[not.na], col="lightgray", ...)
         lines(seq_len(k), z, ...)
         points(x=seq_len(k), y=z, pch=21, bg="black", ...)
         axis(side=1, at=seq_len(k), labels=ids, ...)
         abline(h=0, lty="dashed", ...)
         abline(h=c(qnorm(.025),qnorm(.975)), lty="dotted", ...)

         title("Standardized Residuals", ...)

      }

   } else {

      ######################################################################

      forest(x, ...)
      title("Forest Plot", ...)

      ######################################################################

      funnel(x, ...)
      title("Residual Funnel Plot", ...)

      ######################################################################

      options(na.action = "na.pass")
      z    <- rstandard(x)$z
      pred <- fitted(x)
      options(na.action = na.act)

      plot(pred, z, ylim=c(min(z, -2, na.rm=TRUE), max(z, 2, na.rm=TRUE)), pch=19, bty="l", xlab="Fitted Value", ylab="Standardized Residual", ...)
      abline(h=0, lty="dashed", ...)
      abline(h=c(qnorm(.025),qnorm(.975)), lty="dotted", ...)
      title("Fitted vs. Standardized Residuals", ...)

      ######################################################################

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
         lines(seq_len(k)[not.na], z[not.na], col="lightgray", ...)
         lines(seq_len(k), z, ...)
         points(x=seq_len(k), y=z, pch=21, bg="black", ...)
         axis(side=1, at=seq_len(k), labels=ids, ...)
         abline(h=0, lty="dashed", ...)
         abline(h=c(qnorm(.025),qnorm(.975)), lty="dotted", ...)

         title("Standardized Residuals", ...)

      }

      ######################################################################

   }

   invisible()

}
