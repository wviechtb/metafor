addpoly.predict.rma <- function(x,
rows=-2,                annotate, addpred=FALSE, predstyle, predlim, digits, width, mlab,
transf, atransf, targs, efac, col, border, lty, fonts, cex, constarea=FALSE, ...) {

   #########################################################################

   mstyle <- .get.mstyle()

   .chkclass(class(x), must="predict.rma")

   if (x$pred.type == "scale")
      stop(mstyle$stop("Cannot add polygons based on predicted scale values."))

   if (missing(annotate))
      annotate <- .getfromenv("forest", "annotate", default=TRUE)

   if (missing(predstyle)) {
      predstyle <- "line"
   } else {
      predstyle <- match.arg(predstyle, c("line", "polygon", "bar", "shade", "dist"))
      addpred <- TRUE
   }

   if (missing(predlim))
      predlim <- NULL

   if (missing(digits))
      digits <- .getfromenv("forest", "digits", default=2)

   if (missing(width))
      width <- .getfromenv("forest", "width")

   if (missing(mlab))
      mlab <- NULL

   if (missing(transf))
      transf <- .getfromenv("forest", "transf", default=FALSE)

   if (missing(atransf))
      atransf <- .getfromenv("forest", "atransf", default=FALSE)

   if (missing(targs))
      targs <- .getfromenv("forest", "targs")

   if (missing(efac))
      efac <- .getfromenv("forest", "efac")

   if (missing(col))
      col <- par("fg")

   if (missing(border))
      border <- par("fg")

   if (missing(lty))
      lty <- "dotted"

   if (missing(fonts))
      fonts <- .getfromenv("forest", "fonts")

   if (missing(cex))
      cex <- .getfromenv("forest", "cex")

   if (addpred) {
      pi.lb <- x$pi.lb
      pi.ub <- x$pi.ub
      if (is.null(pi.lb) || is.null(pi.ub))
         warning(mstyle$warning("Could not extract prediction interval bounds."), call.=FALSE)
   } else {
      pi.lb <- rep(NA_real_, length(x$pred))
      pi.ub <- rep(NA_real_, length(x$pred))
   }

   #########################################################################

   addpoly(x$pred, ci.lb=x$ci.lb, ci.ub=x$ci.ub, pi.lb=pi.lb, pi.ub=pi.ub,
           rows=rows,             annotate=annotate, predstyle=predstyle, predlim=predlim,
           digits=digits, width=width, mlab=mlab,
           transf=transf, atransf=atransf, targs=targs,
           efac=efac, col=col, border=border, lty=lty, fonts=fonts, cex=cex,
           constarea=constarea, ...)

}
