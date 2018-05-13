addpoly.rma <- function(x, row=-2, level=x$level, annotate=TRUE,
addcred=FALSE, digits=2, width, mlab, transf, atransf, targs,
efac=1, col, border, fonts, cex, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(x, "rma"))
      stop(mstyle$stop("Argument 'x' must be an object of class \"rma\"."))

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

   if (missing(fonts))
      fonts <- NULL

   if (missing(cex))
      cex <- NULL

   if (addcred) {
      temp <- predict(x, level=level)
      cr.lb <- temp$cr.lb
      cr.ub <- temp$cr.ub
   } else {
      cr.lb <- NA
      cr.ub <- NA
   }

   #########################################################################

   ### label for model estimate (if not specified)

   if (is.null(mlab))
      mlab <- ifelse((x$method=="FE"), "FE Model", "RE Model")

   ### passing ci.lb and ci.ub, so that the bounds are correct when the model was fitted with test="knha"

   addpoly(x$beta, ci.lb=x$ci.lb, ci.ub=x$ci.ub, cr.lb=cr.lb, cr.ub=cr.ub,
           rows=row, level=level, annotate=annotate, digits=digits, width=width,
           mlab=mlab, transf=transf, atransf=atransf, targs=targs,
           efac=efac, col=col, border=border, fonts=fonts, cex=cex, ...)

}
