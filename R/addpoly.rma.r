addpoly.rma <- function(x, row=-2, level=x$level, annotate=TRUE,
addpred=FALSE, digits=2, width, mlab, transf, atransf, targs,
efac=1, col, border, fonts, cex, ...) {

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

   if (addpred) {
      temp <- predict(x, level=level, pi.type=pi.type)
      pi.lb <- temp$pi.lb
      pi.ub <- temp$pi.ub
   } else {
      pi.lb <- NA
      pi.ub <- NA
   }

   #########################################################################

   ### label for model estimate (if not specified)

   if (is.null(mlab))
      mlab <- ifelse((x$method=="FE"), "FE Model", "RE Model")

   ### passing ci.lb and ci.ub, so that the bounds are correct when the model was fitted with test="knha"

   addpoly(x$beta, ci.lb=x$ci.lb, ci.ub=x$ci.ub, pi.lb=pi.lb, pi.ub=pi.ub,
           rows=row, level=level, annotate=annotate, digits=digits, width=width,
           mlab=mlab, transf=transf, atransf=atransf, targs=targs,
           efac=efac, col=col, border=border, fonts=fonts, cex=cex, ...)

}
