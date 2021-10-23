baujat.rma <- function(x, xlim, ylim, xlab, ylab, cex, symbol="ids", grid=TRUE, progbar=FALSE, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="rma", notav=c("rma.glmm", "rma.mv", "robust.rma", "rma.ls", "rma.uni.selmodel"))

   na.act <- getOption("na.action")
   on.exit(options(na.action=na.act))

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   if (x$k == 1)
      stop(mstyle$stop("Stopped because k = 1."))

   ### grid argument can either be a logical or a color

   if (is.logical(grid))
      gridcol <- "lightgray"
   if (is.character(grid)) {
      gridcol <- grid
      grid <- TRUE
   }

   #########################################################################

   ### set up vectors to store results in

   delpred  <- rep(NA_real_, x$k.f)
   vdelpred <- rep(NA_real_, x$k.f)

   ### predicted values under the full model

   pred.full <- x$X.f %*% x$beta

   ### note: skipping NA cases
   ### also: it is possible that model fitting fails, so that generates more NAs (these NAs will always be shown in output)

   if (progbar)
      pbar <- pbapply::startpb(min=0, max=x$k.f)

   for (i in seq_len(x$k.f)) {

      if (progbar)
         pbapply::setpb(pbar, i)

      if (!x$not.na[i])
         next

      if (inherits(x, "rma.uni"))
         res <- try(suppressWarnings(.do.call(rma.uni, yi=x$yi.f, vi=x$vi.f, weights=x$weights.f, mods=x$X.f, intercept=FALSE, method=x$method, weighted=x$weighted, test=x$test, level=x$level, tau2=ifelse(x$tau2.fix, x$tau2, NA), control=x$control, subset=-i, skipr2=TRUE)), silent=TRUE)
      if (inherits(x, "rma.mh")) {
         if (is.element(x$measure, c("RR","OR","RD"))) {
            res <- try(suppressWarnings(.do.call(rma.mh, ai=x$ai.f, bi=x$bi.f, ci=x$ci.f, di=x$di.f, measure=x$measure, add=x$add, to=x$to, drop00=x$drop00, correct=x$correct, level=x$level, subset=-i)), silent=TRUE)
         } else {
            res <- try(suppressWarnings(.do.call(rma.mh, x1i=x$x1i.f, x2i=x$x2i.f, t1i=x$t1i.f, t2i=x$t2i.f, measure=x$measure, add=x$add, to=x$to, drop00=x$drop00, correct=x$correct, level=x$level, subset=-i)), silent=TRUE)
         }
      }

      if (inherits(x, "rma.peto"))
         res <- try(suppressWarnings(.do.call(rma.peto, ai=x$ai.f, bi=x$bi.f, ci=x$ci.f, di=x$di.f, add=x$add, to=x$to, drop00=x$drop00, level=x$level, subset=-i)), silent=TRUE)

      if (inherits(res, "try-error"))
         next

      ### removing an observation could lead to a model coefficient becoming inestimable (for 'rma.uni' objects)

      if (any(res$coef.na))
         next

      Xi          <- matrix(x$X.f[i,], nrow=1)
      delpred[i]  <- Xi %*% res$beta
      vdelpred[i] <- Xi %*% tcrossprod(res$vb,Xi)

   }

   if (progbar)
      pbapply::closepb(pbar)

   yhati <- (delpred - pred.full)^2 / vdelpred

   #########################################################################

   ### x-axis values (use 'na.pass' to make sure we get a vector of length k.f)

   options(na.action = "na.pass")
   xhati <- resid(x)^2 / (x$tau2.f + x$vi.f)
   options(na.action = na.act)

   #########################################################################

   ### set some defaults (if not specified)

   if (missing(cex))
      cex <- 0.8

   if (missing(xlab)) {
      if (is.element(x$method, c("FE","EE","CE"))) {
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

   ### add grid (and redraw box)

   if (.isTRUE(grid)) {
      grid(col=gridcol)
      box(...)
   }

   if (is.numeric(symbol)) {

      if (length(symbol) == 1L)
         symbol <- rep(symbol, x$k.all)

      if (length(symbol) != x$k.all)
         stop(mstyle$stop(paste0("Length of the 'symbol' argument (", length(symbol), ") does not correspond to the size of the original dataset (", x$k.all, ").")))

      if (!is.null(x$subset))
         symbol <- symbol[x$subset]

      points(x=xhati, y=yhati, cex=cex, pch=symbol, ...)

   }

   if (is.character(symbol) && symbol=="ids")
      text(xhati, yhati, x$ids, cex=cex, ...)

   if (is.character(symbol) && symbol=="slab")
      text(xhati, yhati, x$slab, cex=cex, ...)

   #########################################################################

   sav <- data.frame(x=xhati[x$not.na], y=yhati[x$not.na], ids=x$ids[x$not.na], slab=x$slab[x$not.na], stringsAsFactors=FALSE)

   invisible(sav)

}
