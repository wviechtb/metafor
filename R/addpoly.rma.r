# Note: Works with "robust.rma" objects.

addpoly.rma <- function(x, row=-2, level=x$level,
annotate=TRUE, digits=2, width, mlab, transf, atransf, targs,
efac=1, col, border, cex, ...) {

   #########################################################################

   if (!inherits(x, "rma"))
      stop("Argument 'x' must be an object of class \"rma\".")

   if (!x$int.only)
      stop("Fitted model should not contain moderators.")

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

   if (missing(cex))
      cex <- NULL

   if (missing(col))
      col <- "black"

   if (missing(border))
      border <- "black"

   #########################################################################

   ### label for model estimate (if not specified)

   if (is.null(mlab))
      mlab <- ifelse((x$method=="FE"), "FE Model", "RE Model")

   ### passing ci.lb and ci.ub, so that the bounds are correct when the model was fitted with knha=TRUE

   addpoly(x$b, ci.lb=x$ci.lb, ci.ub=x$ci.ub, rows=row, level=level,
           annotate=annotate, digits=digits, width=width, mlab=mlab,
           transf=transf, atransf=atransf, targs=targs,
           efac=efac, col=col, border=border, cex=cex, ...)

}
