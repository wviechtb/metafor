tes.rma <- function(x,
   H0=0, alternative="two.sided", alpha=.05,
   test, tes.alternative="greater", progbar=TRUE, tes.alpha=.10,
   digits, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="rma", notav=c("rma.glmm", "rma.mv", "robust.rma", "rma.ls", "rma.uni.selmodel"))

   na.act <- getOption("na.action")
   on.exit(options(na.action=na.act))

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   ### set defaults for digits

   if (missing(digits)) {
      digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
   }

   if (missing(test))
      test <- NULL

   #########################################################################

   if (x$int.only) {
      theta <- c(x$beta)
   } else {
      options(na.action="na.omit")
      theta <- fitted(x)
      options(na.action = na.act)
   }

   tes.default(c(x$yi), vi=x$vi, H0=H0, alternative=alternative, alpha=alpha, theta=theta, tau2=x$tau2,
   test=test, tes.alternative=tes.alternative, progbar=progbar, tes.alpha=tes.alpha,
   digits=digits, ...)

}
