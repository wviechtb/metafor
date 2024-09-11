addpoly.rma         <- function(x,
row=-2,  level=x$level, annotate, addpred=FALSE, predstyle, predlim, digits, width, mlab,
transf, atransf, targs, efac, col, border, lty, fonts, cex, ...) {

   #########################################################################

   mstyle <- .get.mstyle()

   .chkclass(class(x), must="rma")

   if (!x$int.only)
      stop(mstyle$stop("Fitted model should not contain moderators."))

   if (missing(annotate))
      annotate <- .getfromenv("forest", "annotate", default=TRUE)

   if (missing(predstyle)) {
      predstyle <- "line"
   } else {
      predstyle <- match.arg(predstyle, c("line", "bar", "shade", "dist"))
      addpred <- TRUE
   }

   if (missing(predlim))
      predlim <- NULL

   if (missing(digits))
      digits <- .getfromenv("forest", "digits", default=2)

   if (missing(width))
      width <- .getfromenv("forest", "width", default=NULL)

   if (missing(mlab))
      mlab <- NULL

   if (missing(transf))
      transf <- .getfromenv("forest", "transf", default=FALSE)

   if (missing(atransf))
      atransf <- .getfromenv("forest", "atransf", default=FALSE)

   if (missing(targs))
      targs <- .getfromenv("forest", "targs", default=NULL)

   if (missing(efac))
      efac <- .getfromenv("forest", "efac", default=1)

   if (missing(col))
      col <- par("fg")

   if (missing(border))
      border <- par("fg")

   if (missing(lty))
      lty <- "dotted"

   if (missing(fonts))
      fonts <- .getfromenv("forest", "fonts", default=NULL)

   if (missing(cex))
      cex <- .getfromenv("forest", "cex", default=NULL)

   ddd <- list(...)

   if (!is.null(ddd$addcred))
      addpred <- ddd$addcred

   pi.type <- .chkddd(ddd$pi.type, "default")

   predres <- predict(x, level=level, pi.type=pi.type)

   ci.lb <- predres$ci.lb
   ci.ub <- predres$ci.ub

   if (addpred) {
      pi.lb <- predres$pi.lb
      pi.ub <- predres$pi.ub
      if (is.null(pi.lb) || is.null(pi.ub))
         warning(mstyle$warning("Could not extract prediction interval bounds."), call.=FALSE)
   } else {
      pi.lb <- NA_real_
      pi.ub <- NA_real_
   }

   #########################################################################

   ### label for model estimate (if not specified)

   if (is.null(mlab))
      mlab <- sapply(x$method, switch, "FE"="Fixed-Effect Model", "EE"="Equal-Effects Model", "CE"="Common-Effect Model", "Random-Effects Model", USE.NAMES=FALSE)
      #mlab <- sapply(x$method, switch, "FE"="FE Model", "EE"="EE Model", "CE"="CE Model", "RE Model", USE.NAMES=FALSE)

   ### passing ci.lb and ci.ub, so that the bounds are correct when the model was fitted with test="knha"

   addpoly(x$beta, ci.lb=ci.lb, ci.ub=ci.ub, pi.lb=pi.lb, pi.ub=pi.ub,
           rows=row, level=level, annotate=annotate, predstyle=predstyle, predlim=predlim,
           digits=digits, width=width, mlab=mlab,
           transf=transf, atransf=atransf, targs=targs,
           efac=efac, col=col, border=border, lty=lty, fonts=fonts, cex=cex, ...)

}
