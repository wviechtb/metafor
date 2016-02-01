# Note: Works with "robust.rma" objects.

forest.rma <- function(x, annotate=TRUE, addfit=TRUE, addcred=FALSE, showweights=FALSE,
xlim, alim, clim, ylim, at, steps=5, level=x$level, refline=0, digits=2, width,
xlab, slab, mlab, ilab, ilab.xpos, ilab.pos, order,
transf, atransf, targs, rows,
efac=1, pch=15, psize, col, border, lty,
cex, cex.lab, cex.axis, ...) {

   #########################################################################

   if (!inherits(x, "rma"))
      stop("Argument 'x' must be an object of class \"rma\".")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop("Unknown 'na.action' specified under options().")

   #if (!is.null(order))
   #   order <- match.arg(order, c("obs", "fit", "prec", "resid", "rstandard", "abs.resid", "abs.rstandard"))

   if (missing(transf))
      transf <- FALSE

   if (missing(atransf))
      atransf <- FALSE

   transf.char  <- deparse(substitute(transf))
   atransf.char <- deparse(substitute(atransf))

   if (is.function(transf) && is.function(atransf))
      stop("Use either 'transf' or 'atransf' to specify a transformation (not both).")

   if (missing(targs))
      targs <- NULL

   if (missing(at))
      at <- NULL

   if (missing(ilab))
      ilab <- NULL

   if (missing(ilab.xpos))
      ilab.xpos <- NULL

   if (missing(ilab.pos))
      ilab.pos <- NULL

   if (missing(order))
      order <- NULL

   if (missing(psize))
      psize <- NULL

   if (missing(cex))
      cex <- NULL

   if (missing(cex.lab))
      cex.lab <- NULL

   if (missing(cex.axis))
      cex.axis <- NULL

   ### set default colors if user has not specified 'col' and 'border' arguments

   if (x$int.only) {

      if (missing(col)) {
         col <- c("black", "gray50") ### 1st color for summary polygon, 2nd color for credibility interval
      } else {
         if (length(col) == 1L)      ### if user only specified one value, assume it is for summary polygon
            col <- c(col, "gray50")
      }

      if (missing(border))
         border <- "black"           ### border color of summary polygon

   } else {

      if (missing(col))
         col <- "gray"               ### color of fitted values

      if (missing(border))
         border <- "gray"            ### border color of fitted values

   }

   if (missing(lty)) {
      lty <- c("solid", "dotted", "solid") ### 1st value = CIs, 2nd value = credibility interval, 3rd = horizontal line(s)
   } else {
      if (length(lty) == 1L)
         lty <- c(lty, "dotted", "solid")
      if (length(lty) == 2L)
         lty <- c(lty, "solid")
   }

   if (length(efac) == 1L)
      efac <- rep(efac, 2) ### 1st = segments/arrows, 2nd = summary polygon or fitted polygons

   measure <- x$measure

   ### TODO: remove this when there is a weights() function for 'rma.glmm' objects
   if (inherits(x, "rma.glmm") && showweights)
      stop("Option 'showweights=TRUE' currently not possible for 'rma.glmm' objects. Sorry!")

   #########################################################################

   ### digits[1] for annotations, digits[2] for x-axis labels

   if (length(digits) == 1L)
      digits <- c(digits,digits)

   alpha <- ifelse(level > 1, (100-level)/100, 1-level)

   ### extract data and study labels
   ### note: yi.f/vi.f and pred may contain NAs

   yi <- x$yi.f
   vi <- x$vi.f
   X  <- x$X.f

   k <- length(yi)                              ### length of yi.f

   if (missing(slab)) {
      if (x$slab.null) {
         slab <- paste("Study", x$slab)        ### x$slab is always of length yi.f (i.e., NAs also have an slab)
      } else {
         slab <- x$slab                         ### note: slab must have same length as yi.f in rma object
      }                                         ### even when fewer studies used for model fitting (due to NAs)
   }

   if (length(yi) != length(slab))
      stop("Number of outcomes does not correspond to the length of the 'slab' argument.")

   if (is.vector(ilab) || NCOL(ilab) == 1)      ### note: ilab must have same length as yi.f in rma object
      ilab <- cbind(ilab)                       ### even when fewer studies used for model fitting

   ### note: is.vector() is FALSE if ilab is a factor, so use NCOL() to test for this case

   if (length(pch) == 1L)                       ### note: pch must have same length as yi.f in rma object
      pch <- rep(pch, k)                        ### or be equal to a single value (which is then repeated)

   if (length(pch) != length(yi))
      stop("Number of outcomes does not correspond to the length of the 'pch' argument.")

   ### extract fitted values

   options(na.action = "na.pass")               ### using na.exclude to get the entire vector (length of yi.f)

      if (x$int.only) {
         pred <- fitted(x)
         pred.ci.lb <- rep(NA_real_, k)
         pred.ci.ub <- rep(NA_real_, k)
      } else {
         temp <- predict(x, level=level)
         pred <- temp$pred
         if (addcred) {
            pred.ci.lb <- temp$cr.lb
            pred.ci.ub <- temp$cr.ub
         } else {
            pred.ci.lb <- temp$ci.lb
            pred.ci.ub <- temp$ci.ub
         }
      }

      if (inherits(x, "rma.glmm")) {            ### TODO: change this when there is a weights() function for 'rma.glmm' objects
         weights <- NULL
      } else {
         weights <- weights(x)                  ### these are the weights used for the actual model fitting
      }

   options(na.action = na.act)

   ### if user has set the point sizes

   if (!is.null(psize)) {                       ### note: psize must have same length as yi.f (including NAs)
      if (length(psize) == 1L)                  ### or be equal to a single value (which is then repeated)
         psize <- rep(psize, k)
      if (length(psize) != length(yi))
         stop("Number of outcomes does not correspond to the length of the 'psize' argument.")
   }

   ### sort the data if requested

   if (!is.null(order)) {

      if (is.character(order)) {

         if (length(order) != 1)
            stop("Incorrect length of 'order' argument.")

         if (order == "obs")
            sort.vec <- order(yi)
         if (order == "fit")
            sort.vec <- order(pred)
         if (order == "prec")
            sort.vec <- order(vi, yi)
         if (order == "resid")
            sort.vec <- order(yi-pred, yi)
         if (order == "rstandard")
            sort.vec <- order(rstandard(x)$z, yi)
         if (order == "abs.resid")
            sort.vec <- order(abs(yi-pred), yi)
         if (order == "abs.rstandard")
            sort.vec <- order(abs(rstandard(x)$z), yi)

      } else {
         sort.vec <- order                      ### in principle, can also subset with the order argument
      }

      yi         <- yi[sort.vec]
      vi         <- vi[sort.vec]
      X          <- X[sort.vec,,drop=FALSE]
      slab       <- slab[sort.vec]
      ilab       <- ilab[sort.vec,,drop=FALSE]  ### if ilab is still NULL, then this remains NULL
      pred       <- pred[sort.vec]
      pred.ci.lb <- pred.ci.lb[sort.vec]
      pred.ci.ub <- pred.ci.ub[sort.vec]
      weights    <- weights[sort.vec]
      pch        <- pch[sort.vec]
      psize      <- psize[sort.vec]             ### if psize is still NULL, then this remains NULL

   }

   k <- length(yi)                              ### in case length of k has changed

   ### set rows value

   if (missing(rows)) {
      rows <- k:1
   } else {
      if (length(rows) == 1L) {                 ### note: rows must be a single value or the same
         rows <- rows:(rows-k+1)                ### length of yi.f (including NAs) *after ordering/subsetting*
      }
   }

   if (length(rows) != length(yi))
      stop("Number of outcomes does not correspond to the length of the 'rows' argument.")

   ### reverse order

   yi         <- yi[k:1]
   vi         <- vi[k:1]
   X          <- X[k:1,,drop=FALSE]
   slab       <- slab[k:1]
   ilab       <- ilab[k:1,,drop=FALSE]          ### if ilab is still NULL, then this remains NULL
   pred       <- pred[k:1]
   pred.ci.lb <- pred.ci.lb[k:1]
   pred.ci.ub <- pred.ci.ub[k:1]
   weights    <- weights[k:1]
   pch        <- pch[k:1]
   psize      <- psize[k:1]                     ### if psize is still NULL, then this remains NULL
   rows       <- rows[k:1]

   ### check for NAs in yi/vi and act accordingly

   yiviX.na <- is.na(yi) | is.na(vi) | apply(is.na(X), 1, any)

   if (any(yiviX.na)) {

      not.na <- !yiviX.na

      if (na.act == "na.omit") {
         yi         <- yi[not.na]
         vi         <- vi[not.na]
         X          <- X[not.na,,drop=FALSE]
         slab       <- slab[not.na]
         ilab       <- ilab[not.na,,drop=FALSE] ### if ilab is still NULL, then this remains NULL
         pred       <- pred[not.na]
         pred.ci.lb <- pred.ci.lb[not.na]
         pred.ci.ub <- pred.ci.ub[not.na]
         weights    <- weights[not.na]
         pch        <- pch[not.na]
         psize      <- psize[not.na]            ### if psize is still NULL, then this remains NULL

         rows.new <- rows                       ### rearrange rows due to NAs being omitted from plot
         rows.na  <- rows[!not.na]              ### shift higher rows down according to number of NAs omitted
         for (j in seq_len(length(rows.na))) {
            rows.new[rows >= rows.na[j]] <- rows.new[rows >= rows.na[j]] - 1
         }
         rows <- rows.new[not.na]

      }

      if (na.act == "na.fail")
         stop("Missing values in results.")

   }                                            ### note: yi/vi may be NA if na.act == "na.exclude" or "na.pass"

   k <- length(yi)                              ### in case length of k has changed

   ### calculate individual CI bounds

   ci.lb <- yi - qnorm(alpha/2, lower.tail=FALSE) * sqrt(vi)
   ci.ub <- yi + qnorm(alpha/2, lower.tail=FALSE) * sqrt(vi)

   ### if requested, apply transformation to yi's and CI bounds

   if (is.function(transf)) {
      if (is.null(targs)) {
         yi         <- sapply(yi, transf)
         ci.lb      <- sapply(ci.lb, transf)
         ci.ub      <- sapply(ci.ub, transf)
         pred       <- sapply(pred, transf)
         pred.ci.lb <- sapply(pred.ci.lb, transf)
         pred.ci.ub <- sapply(pred.ci.ub, transf)
      } else {
         yi         <- sapply(yi, transf, targs)
         ci.lb      <- sapply(ci.lb, transf, targs)
         ci.ub      <- sapply(ci.ub, transf, targs)
         pred       <- sapply(pred, transf, targs)
         pred.ci.lb <- sapply(pred.ci.lb, transf, targs)
         pred.ci.ub <- sapply(pred.ci.ub, transf, targs)
      }
   }

   ### make sure order of intervals is always increasing

   tmp <- .psort(ci.lb, ci.ub)
   ci.lb <- tmp[,1]
   ci.ub <- tmp[,2]

   tmp <- .psort(pred.ci.lb, pred.ci.ub)
   pred.ci.lb <- tmp[,1]
   pred.ci.ub <- tmp[,2]

   ### apply ci limits if specified

   if (!missing(clim)) {
      clim <- sort(clim)
      if (length(clim) != 2L)
         stop("Argument 'clim' must be of length 2.")
      ci.lb[ci.lb < clim[1]] <- clim[1]
      ci.ub[ci.ub > clim[2]] <- clim[2]
      pred.ci.lb[pred.ci.lb < clim[1]] <- clim[1]
      pred.ci.ub[pred.ci.ub > clim[2]] <- clim[2]
   }

   ### set default point sizes (if not specified by user)

   if (is.null(psize)) {
      if (is.null(weights)) {
         if (any(vi <= 0, na.rm=TRUE)) {           ### in case any vi value is zero
            psize <- rep(1, k)
         } else {                                  ### default psize is proportional to inverse standard error
            wi    <- 1/sqrt(vi)                    ### note: vi's that are NA are ignored (but vi's whose yi is
            psize <- wi/sum(wi, na.rm=TRUE)        ### NA are NOT ignored; an unlikely case in practice)
            psize <- (psize - min(psize, na.rm=TRUE)) / (max(psize, na.rm=TRUE) - min(psize, na.rm=TRUE))
            psize <- (psize * 1.0) + 0.5           ### note: only vi's that are still in the subset are used for determining the default point sizes
            if (all(is.na(psize)))                 ### if k=1, then psize is NA, so catch this (and maybe some other problems)
               psize <- rep(1, k)
         }
      } else {
         wi <- weights
         psize <- wi/sum(wi, na.rm=TRUE)
         psize <- (psize - min(psize, na.rm=TRUE)) / (max(psize, na.rm=TRUE) - min(psize, na.rm=TRUE))
         psize <- (psize * 1.0) + 0.5
         if (all(is.na(psize)))
            psize <- rep(1, k)
      }
   }

   #########################################################################

   ### total range of CI bounds

   rng <- max(ci.ub, na.rm=TRUE) - min(ci.lb, na.rm=TRUE)

   if (annotate) {
      if (showweights) {
         plot.multp.l <- 2.00
         plot.multp.r <- 2.00
      } else {
         plot.multp.l <- 1.20
         plot.multp.r <- 1.20
      }
   } else {
      plot.multp.l <- 1.20
      plot.multp.r <- 0.40
   }

   ### set plot limits

   if (missing(xlim)) {
      xlim <- c(min(ci.lb, na.rm=TRUE) - rng * plot.multp.l, max(ci.ub, na.rm=TRUE) + rng * plot.multp.r)
      xlim <- round(xlim, digits[2])
      #xlim[1] <- xlim[1]*max(1, digits[2]/2)
      #xlim[2] <- xlim[2]*max(1, digits[2]/2)
   }

   ### set x axis limits (at argument overrides alim argument)

   alim.spec <- TRUE

   if (missing(alim)) {
      if (is.null(at)) {
         alim <- range(pretty(x=c(min(ci.lb, na.rm=TRUE), max(ci.ub, na.rm=TRUE)), n=steps-1))
         alim.spec <- FALSE
      } else {
         alim <- range(at)
      }
   }

   ### make sure the plot and x axis limits are sorted

   alim <- sort(alim)
   xlim <- sort(xlim)

   ### plot limits must always encompass the yi values

   if (xlim[1] > min(yi, na.rm=TRUE)) { xlim[1] <- min(yi, na.rm=TRUE) }
   if (xlim[2] < max(yi, na.rm=TRUE)) { xlim[2] <- max(yi, na.rm=TRUE) }

   ### x axis limits must always encompass the yi values (no longer required)

   #if (alim[1] > min(yi, na.rm=TRUE)) { alim[1] <- min(yi, na.rm=TRUE) }
   #if (alim[2] < max(yi, na.rm=TRUE)) { alim[2] <- max(yi, na.rm=TRUE) }

   ### plot limits must always encompass the x axis limits

   if (alim[1] < xlim[1]) { xlim[1] <- alim[1] }
   if (alim[2] > xlim[2]) { xlim[2] <- alim[2] }

   ### set y axis limits

   if (missing(ylim)) {
      if (x$int.only && addfit) {
         ylim <- c(-1.5, k+3)
      } else {
         ylim <- c(0.5, k+3)
      }
   } else {
      ylim <- sort(ylim)
   }

   ### generate x axis positions if none are specified

   if (is.null(at)) {
      if (alim.spec) {
         at <- seq(from=alim[1], to=alim[2], length.out=steps)
      } else {
         at <- pretty(x=c(min(ci.lb, na.rm=TRUE), max(ci.ub, na.rm=TRUE)), n=steps-1)
      }
   } else {
      at[at < alim[1]] <- alim[1] ### remove at values that are below or above the axis limits
      at[at > alim[2]] <- alim[2]
      at <- unique(at)
   }

   ### x axis labels (apply transformation to axis labels if requested)

   at.lab <- at

   if (is.function(atransf)) {
      if (is.null(targs)) {
         at.lab <- formatC(sapply(at.lab, atransf), digits=digits[2], format="f", drop0trailing=TRUE)
      } else {
         at.lab <- formatC(sapply(at.lab, atransf, targs), digits=digits[2], format="f", drop0trailing=TRUE)
      }
   } else {
      at.lab <- formatC(at.lab, digits=digits[2], format="f", drop0trailing=TRUE)
   }

   #########################################################################

   ### adjust margins

   par.mar <- par("mar")
   par.mar.adj <- par.mar - c(0,3,1,1)
   par.mar.adj[par.mar.adj < 0] <- 0
   par(mar = par.mar.adj)
   on.exit(par(mar = par.mar))

   ### start plot

   plot(NA, NA, xlim=xlim, ylim=ylim, xlab="", ylab="", yaxt="n", xaxt="n", xaxs="i", bty="n", ...)

   ### horizontal title line

   abline(h=ylim[2]-2, lty=lty[3], ...)

   ### add reference line

   if (is.numeric(refline))
      segments(refline, ylim[1]-5, refline, ylim[2]-2, lty="dotted", ...)

   ### set cex, cex.lab, and cex.axis sizes as a function of the height of the figure

   par.usr <- par("usr")
   height  <- par.usr[4] - par.usr[3]
   lheight <- strheight("O")
   cex.adj <- ifelse(k * lheight > height * 0.8, height/(1.25 * k * lheight), 1)

   if (is.null(cex)) {
      cex <- par("cex") * cex.adj
   } else {
      if (is.null(cex.lab))
         cex.lab <- cex
      if (is.null(cex.axis))
         cex.axis <- cex
   }
   if (is.null(cex.lab))
      cex.lab <- par("cex") * cex.adj
   if (is.null(cex.axis))
      cex.axis <- par("cex") * cex.adj

   #########################################################################

   ### if addfit and not an intercept-only model, add fitted polygons

   if (addfit && !x$int.only) {

      for (i in seq_len(k)) {

         if (is.na(pred[i]))
            next

         polygon(x=c(max(pred.ci.lb[i], alim[1]), pred[i], min(pred.ci.ub[i], alim[2]), pred[i]), y=c(rows[i], rows[i]+(height/100)*cex*efac[2], rows[i], rows[i]-(height/100)*cex*efac[2]), col=col, border=border, ...)

         ### this would only draw intervals if bounds fall within alim range
         #if ((pred.ci.lb[i] > alim[1]) && (pred.ci.ub[i] < alim[2]))
         #   polygon(x=c(pred.ci.lb[i], pred[i], pred.ci.ub[i], pred[i]), y=c(rows[i], rows[i]+(height/100)*cex*efac[2], rows[i], rows[i]-(height/100)*cex*efac[2]), col=col, border=border, ...)

      }

   }

   #########################################################################

   ### if addfit and intercept-only model, add fixed/random-effects model polygon

   if (addfit && x$int.only) {

      if (inherits(x, "rma.mv") && x$withG && x$tau2s > 1) {

         if (!is.logical(addcred)) {
            ### for multiple tau^2 (and gamma^2) values, need to specify level(s) of the inner factor(s) to compute the credibility interval
            ### this can be done via the addcred argument (i.e., instead of using a logical, one specifies the level(s))
            if (length(addcred) == 1)
               addcred <- c(addcred, addcred)
            temp <- predict(x, level=level, tau2.levels=addcred[1], gamma2.levels=addcred[2])
            addcred <- TRUE ### set addcred to TRUE, so if (x$method != "FE" && addcred) further below works
         } else {
            if (addcred) {
               ### here addcred=TRUE, but user has not specified the level, so throw an error
               stop("Need to specify the level of the inner factor(s) via the 'addcred' argument.")
            } else {
               ### here addcred=FALSE, so just use the first tau^2 and gamma^2 arbitrarily (so predict() works)
               temp <- predict(x, level=level, tau2.levels=1, gamma2.levels=1)
            }
         }

      } else {

         temp <- predict(x, level=level)

      }

      b       <- temp$pred
      b.ci.lb <- temp$ci.lb
      b.ci.ub <- temp$ci.ub
      b.cr.lb <- temp$cr.lb
      b.cr.ub <- temp$cr.ub

      if (is.function(transf)) {
         if (is.null(targs)) {
            b       <- sapply(b, transf)
            b.ci.lb <- sapply(b.ci.lb, transf)
            b.ci.ub <- sapply(b.ci.ub, transf)
            b.cr.lb <- sapply(b.cr.lb, transf)
            b.cr.ub <- sapply(b.cr.ub, transf)
         } else {
            b       <- sapply(b, transf, targs)
            b.ci.lb <- sapply(b.ci.lb, transf, targs)
            b.ci.ub <- sapply(b.ci.ub, transf, targs)
            b.cr.lb <- sapply(b.cr.lb, transf, targs)
            b.cr.ub <- sapply(b.cr.ub, transf, targs)
         }
      }

      ### make sure order of intervals is always increasing

      tmp <- .psort(b.ci.lb, b.ci.ub)
      b.ci.lb <- tmp[,1]
      b.ci.ub <- tmp[,2]

      tmp <- .psort(b.cr.lb, b.cr.ub)
      b.cr.lb <- tmp[,1]
      b.cr.ub <- tmp[,2]

      ### apply ci limits if specified

      if (!missing(clim)) {
         b.ci.lb[b.ci.lb < clim[1]] <- clim[1]
         b.ci.ub[b.ci.ub > clim[2]] <- clim[2]
         b.cr.lb[b.cr.lb < clim[1]] <- clim[1]
         b.cr.ub[b.cr.ub > clim[2]] <- clim[2]
      }

      ### add credibility interval

      if (x$method != "FE" && addcred) {

         segments(max(b.cr.lb, alim[1]), -1, min(b.cr.ub, alim[2]), -1, lty=lty[2], col=col[2], ...)

         if (b.cr.lb >= alim[1]) {
            segments(b.cr.lb, -1-(height/150)*cex*efac[1], b.cr.lb, -1+(height/150)*cex*efac[1], col=col[2], ...)
         } else {
            polygon(x=c(alim[1], alim[1]+(1.4/100)*cex*(xlim[2]-xlim[1]), alim[1]+(1.4/100)*cex*(xlim[2]-xlim[1]), alim[1]), y=c(-1, -1+(height/150)*cex*efac[1], -1-(height/150)*cex*efac[1], -1), col=col[2], border=col[2], ...)
         }

         if (b.cr.ub <= alim[2]) {
            segments(b.cr.ub, -1-(height/150)*cex*efac[1], b.cr.ub, -1+(height/150)*cex*efac[1], col=col[2], ...)
         } else {
            polygon(x=c(alim[2], alim[2]-(1.4/100)*cex*(xlim[2]-xlim[1]), alim[2]-(1.4/100)*cex*(xlim[2]-xlim[1]), alim[2]), y=c(-1, -1+(height/150)*cex*efac[1], -1-(height/150)*cex*efac[1], -1), col=col[2], border=col[2], ...)
         }

      }

      ### polygon for the summary estimate

      polygon(x=c(b.ci.lb, b, b.ci.ub, b), y=c(-1, -1+(height/100)*cex*efac[2], -1, -1-(height/100)*cex*efac[2]), col=col[1], border=border, ...)

      ### add label for model estimate

      if (missing(mlab))
         mlab <- ifelse((x$method=="FE"), "FE Model", "RE Model")

      text(xlim[1], -1, mlab, pos=4, cex=cex, ...)

   }

   #########################################################################

   ### add x axis

   axis(side=1, at=at, labels=at.lab, cex.axis=cex.axis, ...)

   ### add x axis label

   if (missing(xlab))
      xlab <- .setlab(measure, transf.char, atransf.char, gentype=1)

   mtext(xlab, side=1, at=min(at) + (max(at)-min(at))/2, line=par("mgp")[1]-0.5, cex=cex.lab, ...)

   ### add CI ends (either | or <> if outside of axis limits)

   for (i in seq_len(k)) {

      ### need to skip missings, as if() check below will otherwise throw an error
      if (is.na(yi[i]) || is.na(vi[i]))
         next

      ### if the lower bound is actually larger than upper x-axis limit, then everything is to the right and just draw a polygon pointing in that direction
      if (ci.lb[i] >= alim[2]) {
         polygon(x=c(alim[2], alim[2]-(1.4/100)*cex*(xlim[2]-xlim[1]), alim[2]-(1.4/100)*cex*(xlim[2]-xlim[1]), alim[2]), y=c(rows[i], rows[i]+(height/150)*cex*efac[1], rows[i]-(height/150)*cex*efac[1], rows[i]), col="black", ...)
         next
      }

      ### if the upper bound is actually lower than lower x-axis limit, then everything is to the left and just draw a polygon pointing in that direction
      if (ci.ub[i] <= alim[1]) {
         polygon(x=c(alim[1], alim[1]+(1.4/100)*cex*(xlim[2]-xlim[1]), alim[1]+(1.4/100)*cex*(xlim[2]-xlim[1]), alim[1]), y=c(rows[i], rows[i]+(height/150)*cex*efac[1], rows[i]-(height/150)*cex*efac[1], rows[i]), col="black", ...)
         next
      }

      segments(max(ci.lb[i], alim[1]), rows[i], min(ci.ub[i], alim[2]), rows[i], lty=lty[1], ...)

      if (ci.lb[i] >= alim[1]) {
         segments(ci.lb[i], rows[i]-(height/150)*cex*efac[1], ci.lb[i], rows[i]+(height/150)*cex*efac[1], ...)
      } else {
         polygon(x=c(alim[1], alim[1]+(1.4/100)*cex*(xlim[2]-xlim[1]), alim[1]+(1.4/100)*cex*(xlim[2]-xlim[1]), alim[1]), y=c(rows[i], rows[i]+(height/150)*cex*efac[1], rows[i]-(height/150)*cex*efac[1], rows[i]), col="black", ...)
      }

      if (ci.ub[i] <= alim[2]) {
         segments(ci.ub[i], rows[i]-(height/150)*cex*efac[1], ci.ub[i], rows[i]+(height/150)*cex*efac[1], ...)
      } else {
         polygon(x=c(alim[2], alim[2]-(1.4/100)*cex*(xlim[2]-xlim[1]), alim[2]-(1.4/100)*cex*(xlim[2]-xlim[1]), alim[2]), y=c(rows[i], rows[i]+(height/150)*cex*efac[1], rows[i]-(height/150)*cex*efac[1], rows[i]), col="black", ...)
      }

   }

   ### add study labels on the left

   text(xlim[1], rows, slab, pos=4, cex=cex, ...)

   ### add info labels

   if (!is.null(ilab)) {
      if (is.null(ilab.xpos))
         stop("Must specify 'ilab.xpos' argument when adding information with 'ilab'.")
      if (length(ilab.xpos) != NCOL(ilab))
         stop("Number of 'ilab' columns does not match length of 'ilab.xpos' argument.")
      for (l in seq_len(NCOL(ilab))) {
         text(ilab.xpos[l], rows, ilab[,l], pos=ilab.pos[l], cex=cex, ...)
      }
   }

   ### add study annotations on the right: yi [LB, UB]
   ### and add model fit annotations if requested: b [LB, UB]
   ### (have to add this here, so that alignment is correct)

   if (annotate) {

      if (is.function(atransf)) {

         if (is.null(targs)) {
            if (addfit && x$int.only) {
               annotext <- cbind(sapply(c(yi, b), atransf), sapply(c(ci.lb, b.ci.lb), atransf), sapply(c(ci.ub, b.ci.ub), atransf))
            } else {
               annotext <- cbind(sapply(yi, atransf), sapply(ci.lb, atransf), sapply(ci.ub, atransf))
            }
         } else {
            if (addfit && x$int.only) {
               annotext <- cbind(sapply(c(yi, b), atransf, targs), sapply(c(ci.lb, b.ci.lb), atransf, targs), sapply(c(ci.ub, b.ci.ub), atransf, targs))
            } else {
               annotext <- cbind(sapply(yi, atransf, targs), sapply(ci.lb, atransf, targs), sapply(ci.ub, atransf, targs))
            }
         }

         ### make sure order of intervals is always increasing

         tmp <- .psort(annotext[,2:3])
         annotext[,2:3] <- tmp

      } else {

         if (addfit && x$int.only) {
            annotext <- cbind(c(yi, b), c(ci.lb, b.ci.lb), c(ci.ub, b.ci.ub))
         } else {
            annotext <- cbind(yi, ci.lb, ci.ub)
         }

      }

      if (showweights) {
         if (addfit && x$int.only) {
            annotext <- cbind(c(weights,100), annotext)
         } else {
            annotext <- cbind(weights, annotext)
         }
      }

      annotext <- formatC(annotext, format="f", digits=digits[1])

      if (missing(width)) {
         width <- apply(annotext, 2, function(x) max(nchar(x)))
      } else {
         if (length(width) == 1L)
            width <- rep(width, ncol(annotext))
      }

      for (j in 1:ncol(annotext)) {
         annotext[,j] <- formatC(annotext[,j], width=width[j])
      }

      if (showweights) {
         annotext <- cbind(annotext[,1], "%   ", annotext[,2], " [", annotext[,3], ", ", annotext[,4], "]")
      } else {
         annotext <- cbind(annotext[,1], " [", annotext[,2], ", ", annotext[,3], "]")
      }

      annotext <- apply(annotext, 1, paste, collapse="")

      if (addfit && x$int.only) {
         text(x=xlim[2], c(rows,-1), labels=annotext, pos=2, cex=cex, ...)
      } else {
         text(x=xlim[2], rows, labels=annotext, pos=2, cex=cex, ...)
      }

   }

   ### add yi points

   for (i in seq_len(k)) {

      ### need to skip missings, as if() check below will otherwise throw an error
      if (is.na(yi[i]))
         next

      if (yi[i] >= alim[1] && yi[i] <= alim[2])
         points(yi[i], rows[i], pch=pch[i], cex=cex*psize[i], ...)

   }

   #points(yi, rows, pch=pch, cex=cex*psize, ...)

   ### add horizontal line at 0 for the standard FE/RE model display

   if (x$int.only && addfit)
      abline(h=0, lty=lty[3], ...)

   #########################################################################

   ### return some information about plot invisibly

   res <- list('xlim'=par("usr")[1:2], 'alim'=alim, 'at'=at, 'ylim'=ylim, 'rows'=rows, 'cex'=cex, 'cex.lab'=cex.lab, 'cex.axis'=cex.axis)

   invisible(res)

}
