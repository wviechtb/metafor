addpoly.rma <- function(x, row=-2, level=x$level, annotate=TRUE,
addpred=FALSE, digits=2, width, mlab, transf, atransf, targs,
efac=1, col, border, lty, fonts, cex, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="rma")

   if (!x$int.only)
      stop(mstyle$stop("Fitted model should not contain moderators."))

   if (missing(width))
      width <- NULL

   if (missing(mlab))
      mlab <- NULL

   if (missing(transf))
      transf <- FALSE

   if (missing(atransf))
      atransf <- FALSE

   if (missing(targs))
      targs <- NULL

   if (missing(col))
      col <- "black"

   if (missing(border))
      border <- "black"

   if (missing(lty))
      lty <- "dotted"

   if (missing(fonts))
      fonts <- NULL

   if (missing(cex))
      cex <- NULL

   ddd <- list(...)

   if (!is.null(ddd$addcred))
      addpred <- ddd$addcred

   if (is.null(ddd$pi.type)) {
      pi.type <- "default"
   } else {
      pi.type <- ddd$pi.type
   }

   pred <- predict(x, level=level, pi.type=pi.type)

   ci.lb <- pred$ci.lb
   ci.ub <- pred$ci.ub

   if (addpred) {
      pi.lb <- pred$pi.lb
      pi.ub <- pred$pi.ub
   } else {
      pi.lb <- NA
      pi.ub <- NA
   }

   #########################################################################

   ### label for model estimate (if not specified)

   if (is.null(mlab))
      mlab <- sapply(x$method, switch, "FE"="FE Model", "EE"="EE Model", "CE"="CE Model", "RE Model", USE.NAMES=FALSE)

   ### passing ci.lb and ci.ub, so that the bounds are correct when the model was fitted with test="knha"

   addpoly(x$beta, ci.lb=ci.lb, ci.ub=ci.ub, pi.lb=pi.lb, pi.ub=pi.ub,
           rows=row, level=level, annotate=annotate, digits=digits, width=width,
           mlab=mlab, transf=transf, atransf=atransf, targs=targs,
           efac=efac, col=col, border=border, lty=lty, fonts=fonts, cex=cex, ...)

}
