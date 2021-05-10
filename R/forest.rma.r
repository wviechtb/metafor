forest.rma <- function(x,
annotate=TRUE, addfit=TRUE, addpred=FALSE, showweights=FALSE, header=FALSE,
xlim, alim, olim, ylim, top=3, at, steps=5, level=x$level, refline=0, digits=2L, width,
xlab, slab, mlab, ilab, ilab.xpos, ilab.pos, order,
transf, atransf, targs, rows,
efac=1, pch=15, psize, plim=c(0.5,1.5), colout, col, border,
lty, fonts, cex, cex.lab, cex.axis, annosym, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="rma", notav="rma.ls")

   na.act <- getOption("na.action")
   on.exit(options(na.action=na.act))

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   if (missing(transf))
      transf <- FALSE

   if (missing(atransf))
      atransf <- FALSE

   transf.char  <- deparse(substitute(transf))
   atransf.char <- deparse(substitute(atransf))

   if (is.function(transf) && is.function(atransf))
      stop(mstyle$stop("Use either 'transf' or 'atransf' to specify a transformation (not both)."))

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

   if (missing(colout))
      colout <- "black"

   if (missing(psize))
      psize <- NULL

   if (missing(cex))
      cex <- NULL

   if (missing(cex.lab))
      cex.lab <- NULL

   if (missing(cex.axis))
      cex.axis <- NULL

   level <- ifelse(level == 0, 1, ifelse(level >= 1, (100-level)/100, ifelse(level > .5, 1-level, level)))

   ### digits[1] for annotations, digits[2] for x-axis labels
   ### note: digits can also be a list (e.g., digits=list(2L,3)); trailing 0's are dropped for intergers

   if (length(digits) == 1L)
      digits <- c(digits,digits)

   ############################################################################

   ### set default colors if user has not specified 'col' and 'border' arguments

   if (x$int.only) {

      if (missing(col)) {
         col <- c("black", "gray50") # 1st color for summary polygon, 2nd color for prediction interval line
      } else {
         if (length(col) == 1L)      # if user only specified one value, assume it is for the summary polygon
            col <- c(col, "gray50")
      }

      if (missing(border))
         border <- "black"           # border color of summary polygon

   } else {

      if (missing(col))
         col <- "gray"               # color of fitted values

      if (missing(border))
         border <- "gray"            # border color of fitted values

   }

   ### set default line types if user has not specified 'lty' argument

   if (missing(lty)) {
      lty <- c("solid", "dotted", "solid") # 1st value = CIs, 2nd value = prediction interval, 3rd = horizontal line(s)
   } else {
      if (length(lty) == 1L)
         lty <- c(lty, "dotted", "solid")
      if (length(lty) == 2L)
         lty <- c(lty, "solid")
   }

   ### vertical expansion factor: 1st = CI end lines, 2nd = arrows, 3rd = summary polygon or fitted polygons

   if (length(efac) == 1L)
      efac <- rep(efac, 3)

   if (length(efac) == 2L)
      efac <- c(efac[1], efac[1], efac[2]) # if 2 values specified: 1st = CI end lines and arrows, 2nd = summary polygon or fitted polygons

   ### annotation symbols vector

   if (missing(annosym))
      annosym <- c(" [", ", ", "]", "-") # 4th element for minus sign symbol
   if (length(annosym) == 3L)
      annosym <- c(annosym, "-")
   if (length(annosym) != 4L)
      stop(mstyle$stop("Argument 'annosym' must be a vector of length 3."))

   ### get measure from object

   measure <- x$measure

   ### column header

   estlab <- .setlab(measure, transf.char, atransf.char, gentype=3, short=TRUE)
   if (is.expression(estlab)) {
      header.right <- parse(text=paste0("bold(", estlab, " * '", annosym[1], "' * '", 100*(1-level), "% CI'", " * '", annosym[3], "')"))
   } else {
      header.right <- paste0(estlab, annosym[1], 100*(1-level), "% CI", annosym[3])
   }

   if (is.logical(header)) {
      if (header) {
         header.left <- "Study"
      } else {
         header.left <- NULL
         header.right <- NULL
      }
   } else {
      if (!is.character(header))
         stop(mstyle$stop("Argument 'header' must either be a logical or character vector."))
      if (length(header) == 1L) {
         header.left <- header
      } else {
         header.left <- header[1]
         header.right <- header[2]
      }
   }

   if (!annotate)
      header.right <- NULL

   ddd <- list(...)

   if (!is.null(ddd$addcred))
      addpred <- ddd$addcred

   if (is.null(ddd$pi.type)) {
      pi.type <- "default"
   } else {
      pi.type <- ddd$pi.type
   }

   if (is.null(ddd$decreasing)) {
      decreasing <- FALSE
   } else {
      decreasing <- ddd$decreasing
   }

   if (!is.null(ddd$clim))
      olim <- ddd$clim

   lplot     <- function(..., textpos, addcred, pi.type, decreasing, clim) plot(...)
   labline   <- function(..., textpos, addcred, pi.type, decreasing, clim) abline(...)
   lsegments <- function(..., textpos, addcred, pi.type, decreasing, clim) segments(...)
   laxis     <- function(..., textpos, addcred, pi.type, decreasing, clim) axis(...)
   lmtext    <- function(..., textpos, addcred, pi.type, decreasing, clim) mtext(...)
   lpolygon  <- function(..., textpos, addcred, pi.type, decreasing, clim) polygon(...)
   ltext     <- function(..., textpos, addcred, pi.type, decreasing, clim) text(...)
   lpoints   <- function(..., textpos, addcred, pi.type, decreasing, clim) points(...)

   ### TODO: remove this when there is a weights() function for 'rma.glmm' objects
   if (inherits(x, "rma.glmm") && showweights)
      stop(mstyle$stop("Option 'showweights=TRUE' not possible for 'rma.glmm' objects."))

   ### TODO: remove this when there is a weights() function for 'rma.uni.selmodel' objects
   if (inherits(x, "rma.uni.selmodel") && showweights)
      stop(mstyle$stop("Option 'showweights=TRUE' not possible for 'rma.uni.selmodel' objects."))

   if (!is.null(ddd$subset))
      stop(mstyle$stop("Function does not have a 'subset' argument."))

   #########################################################################

   ### extract data and study labels
   ### note: yi.f/vi.f and pred may contain NAs

   yi <- x$yi.f
   vi <- x$vi.f
   X  <- x$X.f

   k <- length(yi)                              # length of yi.f

   ### note: slab (if specified), ilab (if specified), pch (if vector), psize (if
   ###       vector), colout (if vector), order (if vector) must have the same
   ###       length as the original dataset

   if (missing(slab)) {

      if (x$slab.null) {
         slab <- paste("Study", x$ids)          # x$ids is always of length yi.f (i.e., NAs also have an id)
      } else {
         slab <- x$slab                         # x$slab is always of length yi.f (i.e., NAs also have a study label)
      }

   } else {

      if (is.null(slab) || (length(slab) == 1L && is.na(slab))) # slab=NULL or slab=NA can be used to suppress study labels
         slab <- rep("", x$k.all)

      if (length(slab) != x$k.all)
         stop(mstyle$stop(paste0("Length of the 'slab' argument (", length(slab), ") does not correspond to the size of the original dataset (", x$k.all, ").")))

      if (!is.null(x$subset))
         slab <- slab[x$subset]

   }

   if (!is.null(ilab)) {

      if (is.null(dim(ilab)))
         ilab <- cbind(ilab)

      if (nrow(ilab) != x$k.all)
         stop(mstyle$stop(paste0("Length of the 'ilab' argument (", nrow(ilab), ") does not correspond to the size of the original dataset (", x$k.all, ").")))

      if (!is.null(x$subset))
         ilab <- ilab[x$subset,,drop=FALSE]

   }

   if (length(pch) == 1L)
      pch <- rep(pch, x$k.all)

   if (length(pch) != x$k.all)
      stop(mstyle$stop(paste0("Length of the 'pch' argument (", length(pch), ") does not correspond to the size of the original dataset (", x$k.all, ").")))

   if (!is.null(x$subset))
      pch <- pch[x$subset]

   if (!is.null(psize)) {

      if (length(psize) == 1L)
         psize <- rep(psize, x$k.all)

      if (length(psize) != x$k.all)
         stop(mstyle$stop(paste0("Length of the 'psize' argument (", length(psize), ") does not correspond to the size of the original dataset (", x$k.all, ").")))

      if (!is.null(x$subset))
         psize <- psize[x$subset]

   }

   if (length(colout) == 1L)
      colout <- rep(colout, x$k.all)

   if (length(colout) != x$k.all)
      stop(mstyle$stop(paste0("Length of the 'colout' argument (", length(colout), ") does not correspond to the size of the original dataset (", x$k.all, ").")))

   if (!is.null(x$subset))
      colout <- colout[x$subset]

   ### extract fitted values

   options(na.action = "na.pass") # using na.pass to get the entire vector (length of yi.f)

      if (x$int.only) {
         pred <- fitted(x)
         pred.ci.lb <- rep(NA_real_, k)
         pred.ci.ub <- rep(NA_real_, k)
      } else {
         temp <- predict(x, level=level, pi.type=pi.type)
         pred <- temp$pred
         if (addpred) {
            pred.ci.lb <- temp$pi.lb
            pred.ci.ub <- temp$pi.ub
         } else {
            pred.ci.lb <- temp$ci.lb
            pred.ci.ub <- temp$ci.ub
         }
      }

      weights <- try(weights(x), silent=TRUE) # does not work for rma.glmm and rma.uni.selmodel objects

      if (inherits(weights, "try-error"))
         weights <- rep(1, k)

   ### sort the data if requested

   if (!is.null(order)) {

      if (length(order) == 1L) {

         order <- match.arg(order, c("obs", "yi", "fit", "prec", "vi", "resid", "rstandard", "abs.resid", "abs.rstandard"))

         if (order == "obs" || order == "yi")
            sort.vec <- order(yi)
         if (order == "fit")
            sort.vec <- order(pred)
         if (order == "prec" || order == "vi")
            sort.vec <- order(vi, yi)
         if (order == "resid")
            sort.vec <- order(yi-pred, yi)
         if (order == "rstandard")
            sort.vec <- order(rstandard(x)$z, yi)      # need options(na.action = "na.pass") here as well
         if (order == "abs.resid")
            sort.vec <- order(abs(yi-pred), yi)
         if (order == "abs.rstandard")
            sort.vec <- order(abs(rstandard(x)$z), yi) # need options(na.action = "na.pass") here as well

      } else {

         if (length(order) != x$k.all)
            stop(mstyle$stop(paste0("Length of the 'order' argument (", length(order), ") does not correspond to the size of the original dataset (", x$k.all, ").")))

         if (grepl("^order\\(", deparse(substitute(order)))) {
            sort.vec <- order
         } else {
            sort.vec <- order(order, decreasing=decreasing)
         }

         if (!is.null(x$subset))
            sort.vec <- sort.vec[x$subset] - sum(!x$subset)

      }

      yi         <- yi[sort.vec]
      vi         <- vi[sort.vec]
      X          <- X[sort.vec,,drop=FALSE]
      slab       <- slab[sort.vec]
      ilab       <- ilab[sort.vec,,drop=FALSE]  # if NULL, remains NULL
      pred       <- pred[sort.vec]
      pred.ci.lb <- pred.ci.lb[sort.vec]
      pred.ci.ub <- pred.ci.ub[sort.vec]
      weights    <- weights[sort.vec]
      pch        <- pch[sort.vec]
      psize      <- psize[sort.vec]             # if NULL, remains NULL
      colout     <- colout[sort.vec]

   }

   options(na.action = na.act)

   k <- length(yi)                              # in case length of k has changed

   ### set rows value

   if (missing(rows)) {
      rows <- k:1
   } else {
      if (length(rows) == 1L) {                 # note: rows must be a single value or the same
         rows <- rows:(rows-k+1)                # length of yi.f (including NAs) *after ordering*
      }
   }

   if (length(rows) != k)
      stop(mstyle$stop(paste0("Length of the 'rows' argument (", length(rows), ") does not correspond to the number of outcomes (", k, ").")))

   ### reverse order

   yi         <- yi[k:1]
   vi         <- vi[k:1]
   X          <- X[k:1,,drop=FALSE]
   slab       <- slab[k:1]
   ilab       <- ilab[k:1,,drop=FALSE]          # if NULL, remains NULL
   pred       <- pred[k:1]
   pred.ci.lb <- pred.ci.lb[k:1]
   pred.ci.ub <- pred.ci.ub[k:1]
   weights    <- weights[k:1]
   pch        <- pch[k:1]
   psize      <- psize[k:1]                     # if NULL, remains NULL
   colout     <- colout[k:1]
   rows       <- rows[k:1]

   ### check for NAs in yi/vi/X and act accordingly

   yiviX.na <- is.na(yi) | is.na(vi) | apply(is.na(X), 1, any)

   if (any(yiviX.na)) {

      not.na <- !yiviX.na

      if (na.act == "na.omit") {
         yi         <- yi[not.na]
         vi         <- vi[not.na]
         X          <- X[not.na,,drop=FALSE]
         slab       <- slab[not.na]
         ilab       <- ilab[not.na,,drop=FALSE] # if NULL, remains NULL
         pred       <- pred[not.na]
         pred.ci.lb <- pred.ci.lb[not.na]
         pred.ci.ub <- pred.ci.ub[not.na]
         weights    <- weights[not.na]
         pch        <- pch[not.na]
         psize      <- psize[not.na]            # if NULL, remains NULL
         colout     <- colout[not.na]

         rows.new <- rows                       # rearrange rows due to NAs being omitted from plot
         rows.na  <- rows[!not.na]              # shift higher rows down according to number of NAs omitted
         for (j in seq_len(length(rows.na))) {
            rows.new[rows >= rows.na[j]] <- rows.new[rows >= rows.na[j]] - 1
         }
         rows <- rows.new[not.na]

      }

      if (na.act == "na.fail")
         stop(mstyle$stop("Missing values in results."))

   }                                            # note: yi/vi may be NA if na.act == "na.exclude" or "na.pass"

   k <- length(yi)                              # in case length of k has changed

   ### calculate individual CI bounds

   ci.lb <- yi - qnorm(level/2, lower.tail=FALSE) * sqrt(vi)
   ci.ub <- yi + qnorm(level/2, lower.tail=FALSE) * sqrt(vi)

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

   ### apply observation/outcome limits if specified

   if (!missing(olim)) {
      if (length(olim) != 2L)
         stop(mstyle$stop("Argument 'olim' must be of length 2."))
      olim <- sort(olim)
      yi[yi < olim[1]] <- olim[1]
      yi[yi > olim[2]] <- olim[2]
      ci.lb[ci.lb < olim[1]] <- olim[1]
      ci.ub[ci.ub > olim[2]] <- olim[2]
      pred.ci.lb[pred.ci.lb < olim[1]] <- olim[1]
      pred.ci.ub[pred.ci.ub > olim[2]] <- olim[2]
   }

   ### set default point sizes (if not specified by user)

   if (is.null(psize)) {
      if (length(plim) < 2L)
         stop(mstyle$stop("Argument 'plim' must be of length 2 or 3."))
      wi <- sqrt(weights)
      if (!is.na(plim[1]) && !is.na(plim[2])) {
         rng <- max(wi, na.rm=TRUE) - min(wi, na.rm=TRUE)
         if (rng <= .Machine$double.eps^0.5) {
            psize <- rep(1, k)
         } else {
            psize <- (wi - min(wi, na.rm=TRUE)) / rng
            psize <- (psize * (plim[2] - plim[1])) + plim[1]
         }
      }
      if (is.na(plim[1]) && !is.na(plim[2])) {
         psize <- wi / max(wi, na.rm=TRUE) * plim[2]
         if (length(plim) == 3L)
            psize[psize <= plim[3]] <- plim[3]
      }
      if (!is.na(plim[1]) && is.na(plim[2])) {
         psize <- wi / min(wi, na.rm=TRUE) * plim[1]
         if (length(plim) == 3L)
            psize[psize >= plim[3]] <- plim[3]
      }
      if (all(is.na(psize)))
         psize <- rep(1, k)
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
      xlim <- round(xlim, digits[[2]])
      #xlim[1] <- xlim[1]*max(1, digits[[2]]/2)
      #xlim[2] <- xlim[2]*max(1, digits[[2]]/2)
   }

   ### set x-axis limits (at argument overrides alim argument)

   alim.spec <- TRUE

   if (missing(alim)) {
      if (is.null(at)) {
         alim <- range(pretty(x=c(min(ci.lb, na.rm=TRUE), max(ci.ub, na.rm=TRUE)), n=steps-1))
         alim.spec <- FALSE
      } else {
         alim <- range(at)
      }
   }

   ### make sure the plot and x-axis limits are sorted

   alim <- sort(alim)
   xlim <- sort(xlim)

   ### plot limits must always encompass the yi values

   if (xlim[1] > min(yi, na.rm=TRUE)) { xlim[1] <- min(yi, na.rm=TRUE) }
   if (xlim[2] < max(yi, na.rm=TRUE)) { xlim[2] <- max(yi, na.rm=TRUE) }

   ### x-axis limits must always encompass the yi values (no longer required)

   #if (alim[1] > min(yi, na.rm=TRUE)) { alim[1] <- min(yi, na.rm=TRUE) }
   #if (alim[2] < max(yi, na.rm=TRUE)) { alim[2] <- max(yi, na.rm=TRUE) }

   ### plot limits must always encompass the x-axis limits

   if (alim[1] < xlim[1]) { xlim[1] <- alim[1] }
   if (alim[2] > xlim[2]) { xlim[2] <- alim[2] }

   ### allow adjustment of position of study labels and annotations via textpos argument

   if (is.null(ddd$textpos))
      ddd$textpos <- c(xlim[1], xlim[2])

   if (length(ddd$textpos) != 2L)
      stop(mstyle$stop("Argument 'textpos' must be of length 2."))

   if (is.na(ddd$textpos[1]))
      ddd$textpos[1] <- xlim[1]

   if (is.na(ddd$textpos[2]))
      ddd$textpos[2] <- xlim[2]

   ### set y-axis limits

   if (missing(ylim)) {
      if (x$int.only && addfit) {
         ylim <- c(-1.5, max(rows, na.rm=TRUE)+top)
      } else {
         ylim <- c(0.5, max(rows, na.rm=TRUE)+top)
      }
   } else {
      ylim <- sort(ylim)
   }

   ### generate x-axis positions if none are specified

   if (is.null(at)) {
      if (alim.spec) {
         at <- seq(from=alim[1], to=alim[2], length.out=steps)
      } else {
         at <- pretty(x=c(min(ci.lb, na.rm=TRUE), max(ci.ub, na.rm=TRUE)), n=steps-1)
      }
   } else {
      at[at < alim[1]] <- alim[1] # remove at values that are below or above the axis limits
      at[at > alim[2]] <- alim[2]
      at <- unique(at)
   }

   ### x-axis labels (apply transformation to axis labels if requested)

   at.lab <- at

   if (is.function(atransf)) {
      if (is.null(targs)) {
         at.lab <- formatC(sapply(at.lab, atransf), digits=digits[[2]], format="f", drop0trailing=is.integer(digits[[2]]))
      } else {
         at.lab <- formatC(sapply(at.lab, atransf, targs), digits=digits[[2]], format="f", drop0trailing=is.integer(digits[[2]]))
      }
   } else {
      at.lab <- formatC(at.lab, digits=digits[[2]], format="f", drop0trailing=is.integer(digits[[2]]))
   }

   #########################################################################

   ### set/get fonts (1st for study labels, 2nd for annotations, 3rd for ilab)
   ### when passing a named vector, the names are for 'family' and the values are for 'font'

   if (missing(fonts)) {
      fonts <- rep(par("family"), 3)
   } else {
      if (length(fonts) == 1L)
         fonts <- rep(fonts, 3)
      if (length(fonts) == 2L)
         fonts <- c(fonts, fonts[1])
   }

   if (is.null(names(fonts)))
      fonts <- structure(c(1L,1L,1L), names=fonts)

   par(family=names(fonts)[1], font=fonts[1])

   ### adjust margins

   par.mar <- par("mar")
   par.mar.adj <- par.mar - c(0,3,1,1)
   par.mar.adj[par.mar.adj < 0] <- 0
   par(mar = par.mar.adj)
   on.exit(par(mar = par.mar), add=TRUE)

   ### start plot

   lplot(NA, NA, xlim=xlim, ylim=ylim, xlab="", ylab="", yaxt="n", xaxt="n", xaxs="i", bty="n", ...)

   ### horizontal title line

   labline(h=ylim[2]-(top-1), lty=lty[3], ...)

   ### get coordinates of the plotting region

   par.usr <- par("usr")

   ### add reference line

   if (is.numeric(refline))
      lsegments(refline, par.usr[3], refline, ylim[2]-(top-1), lty="dotted", ...)

   ### set cex, cex.lab, and cex.axis sizes as a function of the height of the figure

   height  <- par.usr[4] - par.usr[3]

   if (is.null(cex)) {
      lheight <- strheight("O")
      cex.adj <- ifelse(k * lheight > height * 0.8, height/(1.25 * k * lheight), 1)
   }

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

         lpolygon(x=c(max(pred.ci.lb[i], alim[1]), pred[i], min(pred.ci.ub[i], alim[2]), pred[i]), y=c(rows[i], rows[i]+(height/100)*cex*efac[3], rows[i], rows[i]-(height/100)*cex*efac[3]), col=col, border=border, ...)

         ### this would only draw intervals if bounds fall within alim range
         #if ((pred.ci.lb[i] > alim[1]) && (pred.ci.ub[i] < alim[2]))
         #   lpolygon(x=c(pred.ci.lb[i], pred[i], pred.ci.ub[i], pred[i]), y=c(rows[i], rows[i]+(height/100)*cex*efac[3], rows[i], rows[i]-(height/100)*cex*efac[3]), col=col, border=border, ...)

      }

   }

   #########################################################################

   ### if addfit and intercept-only model, add fixed/random-effects model polygon

   if (addfit && x$int.only) {

      if (inherits(x, "rma.mv") && x$withG && x$tau2s > 1) {

         if (!is.logical(addpred)) {
            ### for multiple tau^2 (and gamma^2) values, need to specify level(s) of the inner factor(s) to compute the prediction interval
            ### this can be done via the addpred argument (i.e., instead of using a logical, one specifies the level(s))
            if (length(addpred) == 1L)
               addpred <- c(addpred, addpred)
            temp <- predict(x, level=level, tau2.levels=addpred[1], gamma2.levels=addpred[2], pi.type=pi.type)
            addpred <- TRUE ### set addpred to TRUE, so if (!is.element(x$method, c("FE","EE")) && addpred) further below works
         } else {
            if (addpred) {
               ### here addpred=TRUE, but user has not specified the level, so throw an error
               stop(mstyle$stop("Must specify the level of the inner factor(s) via the 'addpred' argument."))
            } else {
               ### here addpred=FALSE, so just use the first tau^2 and gamma^2 arbitrarily (so predict() works)
               temp <- predict(x, level=level, tau2.levels=1, gamma2.levels=1, pi.type=pi.type)
            }
         }

      } else {

         temp <- predict(x, level=level, pi.type=pi.type)

      }

      beta       <- temp$pred
      beta.ci.lb <- temp$ci.lb
      beta.ci.ub <- temp$ci.ub
      beta.pi.lb <- temp$pi.lb
      beta.pi.ub <- temp$pi.ub

      if (is.function(transf)) {
         if (is.null(targs)) {
            beta       <- sapply(beta, transf)
            beta.ci.lb <- sapply(beta.ci.lb, transf)
            beta.ci.ub <- sapply(beta.ci.ub, transf)
            beta.pi.lb <- sapply(beta.pi.lb, transf)
            beta.pi.ub <- sapply(beta.pi.ub, transf)
         } else {
            beta       <- sapply(beta, transf, targs)
            beta.ci.lb <- sapply(beta.ci.lb, transf, targs)
            beta.ci.ub <- sapply(beta.ci.ub, transf, targs)
            beta.pi.lb <- sapply(beta.pi.lb, transf, targs)
            beta.pi.ub <- sapply(beta.pi.ub, transf, targs)
         }
      }

      ### make sure order of intervals is always increasing

      tmp <- .psort(beta.ci.lb, beta.ci.ub)
      beta.ci.lb <- tmp[,1]
      beta.ci.ub <- tmp[,2]

      tmp <- .psort(beta.pi.lb, beta.pi.ub)
      beta.pi.lb <- tmp[,1]
      beta.pi.ub <- tmp[,2]

      ### apply observation/outcome limits if specified

      if (!missing(olim)) {
         pred[pred < olim[1]] <- olim[1]
         pred[pred > olim[2]] <- olim[2]
         beta.ci.lb[beta.ci.lb < olim[1]] <- olim[1]
         beta.ci.ub[beta.ci.ub > olim[2]] <- olim[2]
         beta.pi.lb[beta.pi.lb < olim[1]] <- olim[1]
         beta.pi.ub[beta.pi.ub > olim[2]] <- olim[2]
      }

      ### add prediction interval

      if (!is.element(x$method, c("FE","EE")) && addpred) {

         lsegments(max(beta.pi.lb, alim[1]), -1, min(beta.pi.ub, alim[2]), -1, lty=lty[2], col=col[2], ...)

         if (beta.pi.lb >= alim[1]) {
            lsegments(beta.pi.lb, -1-(height/150)*cex*efac[1], beta.pi.lb, -1+(height/150)*cex*efac[1], col=col[2], ...)
         } else {
            lpolygon(x=c(alim[1], alim[1]+(1.4/100)*cex*(xlim[2]-xlim[1]), alim[1]+(1.4/100)*cex*(xlim[2]-xlim[1]), alim[1]), y=c(-1, -1+(height/150)*cex*efac[2], -1-(height/150)*cex*efac[2], -1), col=col[2], border=col[2], ...)
         }

         if (beta.pi.ub <= alim[2]) {
            lsegments(beta.pi.ub, -1-(height/150)*cex*efac[1], beta.pi.ub, -1+(height/150)*cex*efac[1], col=col[2], ...)
         } else {
            lpolygon(x=c(alim[2], alim[2]-(1.4/100)*cex*(xlim[2]-xlim[1]), alim[2]-(1.4/100)*cex*(xlim[2]-xlim[1]), alim[2]), y=c(-1, -1+(height/150)*cex*efac[2], -1-(height/150)*cex*efac[2], -1), col=col[2], border=col[2], ...)
         }

      }

      ### polygon for the summary estimate

      lpolygon(x=c(beta.ci.lb, beta, beta.ci.ub, beta), y=c(-1, -1+(height/100)*cex*efac[3], -1, -1-(height/100)*cex*efac[3]), col=col[1], border=border, ...)

      ### add label for model estimate

      if (missing(mlab))
         mlab <- sapply(x$method, switch, "FE"="FE Model", "EE"="EE Model", "RE Model", USE.NAMES=FALSE)

      ltext(ddd$textpos[1], -1, mlab, pos=4, cex=cex, ...)

   }

   #########################################################################

   ### add x-axis

   laxis(side=1, at=at, labels=at.lab, cex.axis=cex.axis, ...)

   ### add x-axis label

   if (missing(xlab))
      xlab <- .setlab(measure, transf.char, atransf.char, gentype=1)

   lmtext(xlab, side=1, at=min(at) + (max(at)-min(at))/2, line=par("mgp")[1]-0.5, cex=cex.lab, ...)

   ### add CI ends (either | or <> if outside of axis limits)

   for (i in seq_len(k)) {

      ### need to skip missings (if check below will otherwise throw an error)
      if (is.na(yi[i]) || is.na(vi[i]))
         next

      ### if the lower bound is actually larger than upper x-axis limit, then everything is to the right and just draw a polygon pointing in that direction
      if (ci.lb[i] >= alim[2]) {
         lpolygon(x=c(alim[2], alim[2]-(1.4/100)*cex*(xlim[2]-xlim[1]), alim[2]-(1.4/100)*cex*(xlim[2]-xlim[1]), alim[2]), y=c(rows[i], rows[i]+(height/150)*cex*efac[2], rows[i]-(height/150)*cex*efac[2], rows[i]), col=colout[i], ...)
         next
      }

      ### if the upper bound is actually lower than lower x-axis limit, then everything is to the left and just draw a polygon pointing in that direction
      if (ci.ub[i] <= alim[1]) {
         lpolygon(x=c(alim[1], alim[1]+(1.4/100)*cex*(xlim[2]-xlim[1]), alim[1]+(1.4/100)*cex*(xlim[2]-xlim[1]), alim[1]), y=c(rows[i], rows[i]+(height/150)*cex*efac[2], rows[i]-(height/150)*cex*efac[2], rows[i]), col=colout[i], ...)
         next
      }

      lsegments(max(ci.lb[i], alim[1]), rows[i], min(ci.ub[i], alim[2]), rows[i], lty=lty[1], col=colout[i], ...)

      if (ci.lb[i] >= alim[1]) {
         lsegments(ci.lb[i], rows[i]-(height/150)*cex*efac[1], ci.lb[i], rows[i]+(height/150)*cex*efac[1], col=colout[i], ...)
      } else {
         lpolygon(x=c(alim[1], alim[1]+(1.4/100)*cex*(xlim[2]-xlim[1]), alim[1]+(1.4/100)*cex*(xlim[2]-xlim[1]), alim[1]), y=c(rows[i], rows[i]+(height/150)*cex*efac[2], rows[i]-(height/150)*cex*efac[2], rows[i]), col=colout[i], ...)
      }

      if (ci.ub[i] <= alim[2]) {
         lsegments(ci.ub[i], rows[i]-(height/150)*cex*efac[1], ci.ub[i], rows[i]+(height/150)*cex*efac[1], col=colout[i], ...)
      } else {
         lpolygon(x=c(alim[2], alim[2]-(1.4/100)*cex*(xlim[2]-xlim[1]), alim[2]-(1.4/100)*cex*(xlim[2]-xlim[1]), alim[2]), y=c(rows[i], rows[i]+(height/150)*cex*efac[2], rows[i]-(height/150)*cex*efac[2], rows[i]), col=colout[i], ...)
      }

   }

   ### add study labels on the left

   ltext(ddd$textpos[1], rows, slab, pos=4, cex=cex, ...)

   ### add info labels

   if (!is.null(ilab)) {
      if (is.null(ilab.xpos))
         stop(mstyle$stop("Must specify 'ilab.xpos' argument when adding information with 'ilab'."))
      if (length(ilab.xpos) != ncol(ilab))
         stop(mstyle$stop(paste0("Number of 'ilab' columns (", ncol(ilab), ") does not match length of 'ilab.xpos' argument (", length(ilab.xpos), ").")))
      if (!is.null(ilab.pos) && length(ilab.pos) == 1L)
         ilab.pos <- rep(ilab.pos, ncol(ilab))
      par(family=names(fonts)[3], font=fonts[3])
      for (l in seq_len(ncol(ilab))) {
         ltext(ilab.xpos[l], rows, ilab[,l], pos=ilab.pos[l], cex=cex, ...)
      }
      par(family=names(fonts)[1], font=fonts[1])
   }

   ### add study annotations on the right: yi [LB, UB]
   ### and add model fit annotations if requested: b [LB, UB]
   ### (have to add this here, so that alignment is correct)

   if (annotate) {

      if (is.function(atransf)) {

         if (is.null(targs)) {
            if (addfit && x$int.only) {
               annotext <- cbind(sapply(c(yi, beta), atransf), sapply(c(ci.lb, beta.ci.lb), atransf), sapply(c(ci.ub, beta.ci.ub), atransf))
            } else {
               annotext <- cbind(sapply(yi, atransf), sapply(ci.lb, atransf), sapply(ci.ub, atransf))
            }
         } else {
            if (addfit && x$int.only) {
               annotext <- cbind(sapply(c(yi, beta), atransf, targs), sapply(c(ci.lb, beta.ci.lb), atransf, targs), sapply(c(ci.ub, beta.ci.ub), atransf, targs))
            } else {
               annotext <- cbind(sapply(yi, atransf, targs), sapply(ci.lb, atransf, targs), sapply(ci.ub, atransf, targs))
            }
         }

         ### make sure order of intervals is always increasing

         tmp <- .psort(annotext[,2:3])
         annotext[,2:3] <- tmp

      } else {

         if (addfit && x$int.only) {
            annotext <- cbind(c(yi, beta), c(ci.lb, beta.ci.lb), c(ci.ub, beta.ci.ub))
         } else {
            annotext <- cbind(yi, ci.lb, ci.ub)
         }

      }

      if (showweights) {
         if (addfit && x$int.only) {
            annotext <- cbind(c(unname(weights),100), annotext)
         } else {
            annotext <- cbind(unname(weights), annotext)
         }
      }

      annotext <- .fcf(annotext, digits[[1]])
      annotext <- sub("-", annosym[4], annotext, fixed=TRUE)

      if (missing(width)) {
         width <- apply(annotext, 2, function(x) max(nchar(x)))
      } else {
         if (length(width) == 1L)
            width <- rep(width, ncol(annotext))
      }

      for (j in seq_len(ncol(annotext))) {
         annotext[,j] <- formatC(annotext[,j], width=width[j])
      }

      if (showweights) {
         annotext <- cbind(annotext[,1], "%   ", annotext[,2], annosym[1], annotext[,3], annosym[2], annotext[,4], annosym[3])
      } else {
         annotext <- cbind(annotext[,1], annosym[1], annotext[,2], annosym[2], annotext[,3], annosym[3])
      }

      annotext <- apply(annotext, 1, paste, collapse="")
      annotext[grepl("NA", annotext, fixed=TRUE)] <- ""

      par(family=names(fonts)[2], font=fonts[2])
      if (addfit && x$int.only) {
         ltext(ddd$textpos[2], c(rows,-1), labels=annotext, pos=2, cex=cex, ...)
      } else {
         ltext(ddd$textpos[2], rows, labels=annotext, pos=2, cex=cex, ...)
      }
      par(family=names(fonts)[1], font=fonts[1])

   }

   ### add yi points

   for (i in seq_len(k)) {

      ### need to skip missings, as if () check below will otherwise throw an error
      if (is.na(yi[i]))
         next

      if (yi[i] >= alim[1] && yi[i] <= alim[2])
         lpoints(yi[i], rows[i], pch=pch[i], col=colout[i], cex=cex*psize[i], ...)

   }

   #lpoints(yi, rows, pch=pch, cex=cex*psize, ...)

   ### add horizontal line at 0 for the standard FE/RE model display

   if (x$int.only && addfit)
      labline(h=0, lty=lty[3], ...)

   ### add header

   ltext(ddd$textpos[1], ylim[2]-(top-1)+1, header.left,  pos=4, font=2, cex=cex, ...)
   ltext(ddd$textpos[2], ylim[2]-(top-1)+1, header.right, pos=2, font=2, cex=cex, ...)

   #########################################################################

   ### return some information about plot invisibly

   res <- list(xlim=par("usr")[1:2], alim=alim, at=at, ylim=ylim, rows=rows, cex=cex, cex.lab=cex.lab, cex.axis=cex.axis)

   invisible(res)

}
