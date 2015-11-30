baujat.rma <- function(x, xlim, ylim, xlab, ylab, cex, symbol, grid=TRUE, ...) {

   if (!is.element("rma", class(x)))
      stop("Argument 'x' must be an object of class \"rma\".")

   if (is.element("rma.glmm", class(x)))
      stop("Method not yet implemented for objects of class \"rma.glmm\". Sorry!")

   if (is.element("rma.mv", class(x)))
      stop("Method not yet implemented for objects of class \"rma.mv\". Sorry!")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop("Unknown 'na.action' specified under options().")

   if (x$k == 1)
      stop("Stopped because k = 1.")

   #########################################################################

   ### set up vectors to store results in

   delpred  <- rep(NA_real_, x$k.f)
   vdelpred <- rep(NA_real_, x$k.f)

   ### predicted values under the full model

   pred.full <- x$X.f %*% x$b

   ### note: skipping NA cases
   ### also: it is possible that model fitting fails, so that generates more NAs (these NAs will always be shown in output)

   for (i in seq_len(x$k.f)[x$not.na]) {

      if (is.element("rma.uni", class(x)))
         res <- try(suppressWarnings(rma(x$yi.f[-i], x$vi.f[-i], weights=x$weights.f[-i], mods=cbind(x$X.f[-i,]), method=x$method, weighted=x$weighted, intercept=FALSE, knha=x$knha, control=x$control)), silent=TRUE)

      if (is.element("rma.mh", class(x))) {
         if (is.element(x$measure, c("RR","OR","RD"))) {
            res <- try(suppressWarnings(rma.mh(ai=x$ai.f[-i], bi=x$bi.f[-i], ci=x$ci.f[-i], di=x$di.f[-i], measure=x$measure, add=x$add, to=x$to, drop00=x$drop00, correct=x$correct)), silent=TRUE)
         } else {
            res <- try(suppressWarnings(rma.mh(x1i=x$x1i.f[-i], x2i=x$x2i.f[-i], t1i=x$t1i.f[-i], t2i=x$t2i.f[-i], measure=x$measure, add=x$add, to=x$to, drop00=x$drop00, correct=x$correct)), silent=TRUE)
         }
      }

      if (is.element("rma.peto", class(x)))
         res <- try(suppressWarnings(rma.peto(ai=x$ai.f[-i], bi=x$bi.f[-i], ci=x$ci.f[-i], di=x$di.f[-i], add=x$add, to=x$to, drop00=x$drop00)), silent=TRUE)

      if (inherits(res, "try-error"))
         next

      ### removing an observation could lead to a model coefficient becoming inestimable (for 'rma.uni' objects)

      if (any(res$coef.na))
         next

      Xi          <- matrix(x$X.f[i,], nrow=1)
      delpred[i]  <- Xi %*% res$b
      vdelpred[i] <- Xi %*% tcrossprod(res$vb,Xi)

   }

   yhati <- (delpred - pred.full)^2 / vdelpred

   #########################################################################

   ### x-axis values (use 'na.pass' to make sure we get a vector of length k.f)

   options(na.action = "na.pass")
   xhati <- 1/(x$tau2 + x$vi.f) * resid(x)^2
   options(na.action = na.act)

   #########################################################################

   ### set some defaults (if not specified)

   if (missing(cex))
      cex <- 0.8

   if (missing(xlab)) {
      if (x$method == "FE") {
         xlab <- ifelse(x$int.only, "Contribution to Overall Heterogeneity", "Contribution to Residual Heterogeneity")
      } else {
         xlab <- "Squared Pearson Residual"
      }
   }

   if (missing(ylab))
      ylab <- ifelse(x$int.only, "Influence on Overall Result", "Influence on Fitted Value")

   if (missing(xlim))
      xlim <- range(xhati, na.rm=TRUE)

   if (missing(ylim))
      ylim <- range(yhati, na.rm=TRUE)

   #########################################################################

   ### draw empty plot

   plot(NA, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...)

   ### add grid (if requested)

   if (grid)
      grid()

   ### add points/symbols

   if (missing(symbol))
      symbol <- "ids"

   if (is.numeric(symbol))
      points(xhati, yhati, cex=cex, pch=symbol, ...)

   if (is.character(symbol) && symbol == "ids")
      text(xhati, yhati, x$ids, cex=cex, ...)

   if (is.character(symbol) && symbol == "slab")
      text(xhati, yhati, x$slab, cex=cex, ...)

   #########################################################################

   invisible(data.frame(x=xhati[x$not.na], y=yhati[x$not.na], ids=x$ids[x$not.na], slab=x$slab[x$not.na]))

}
