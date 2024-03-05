forest.rma       <- function(x, y,
annotate=TRUE, addfit=TRUE, addpred=FALSE, showweights=FALSE, header=FALSE,
xlim, alim, olim, ylim, at, steps=5, level=x$level, refline=0, digits=2L, width,
xlab, slab, mlab, ilab, ilab.xpos, ilab.pos, order,
transf, atransf, targs, rows,
efac=1, pch, psize, plim=c(0.5,1.5), colout, col, border, shade, colshade,
lty, fonts, cex, cex.lab, cex.axis, ...) {

   #########################################################################

   mstyle <- .get.mstyle()
   
   .chkclass(class(x), must="rma", notav=c("rma.ls", "rma.gen"))
   if (!missing(y))
     .chkclass(class(y), must="rma", notav=c("rma.ls", "rma.gen"))

   na.act <- getOption("na.action")
   on.exit(options(na.action=na.act), add=TRUE)

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   if (missing(transf))
      transf <- FALSE

   if (missing(atransf))
      atransf <- FALSE

   transf.char  <- deparse(transf)
   atransf.char <- deparse(atransf)

   if (is.function(transf) && is.function(atransf))
      stop(mstyle$stop("Use either 'transf' or 'atransf' to specify a transformation (not both)."))

   .start.plot()

   if (missing(targs))
      targs <- NULL

   if (missing(at))
      at <- NULL

   mf <- match.call()

   if (missing(ilab)) {
      ilab <- NULL
   } else {
      ilab <- .getx("ilab", mf=mf, data=x$data)
   }

   if (missing(ilab.xpos))
      ilab.xpos <- NULL

   if (missing(ilab.pos))
      ilab.pos <- NULL

   if (missing(order)) {
      order <- NULL
   } else {
      order <- .getx("order", mf=mf, data=x$data)
   }

   if (missing(colout)) {
      colout <- par("fg")
   } else {
      colout <- .getx("colout", mf=mf, data=x$data)
   }

   if (missing(shade)) {
      shade <- NULL
   } else {
      shade <- .getx("shade", mf=mf, data=x$data)
   }

   if (missing(colshade))
      colshade <- .coladj(par("bg","fg"), dark=0.1, light=-0.1)

   if (missing(pch)) {
      pch <- 15
   } else {
      pch <- .getx("pch", mf=mf, data=x$data)
   }

   if (missing(psize)) {
      psize <- NULL
   } else {
      psize <- .getx("psize", mf=mf, data=x$data)
   }

   if (missing(cex))
      cex <- NULL

   if (missing(cex.lab))
      cex.lab <- NULL

   if (missing(cex.axis))
      cex.axis <- NULL

   level <- .level(level)

   ### digits[1] for annotations, digits[2] for x-axis labels, digits[3] (if specified) for weights
   ### note: digits can also be a list (e.g., digits=list(2,3L)); trailing 0's on the x-axis labels
   ### are dropped if the value is an integer

   if (length(digits) == 1L)
      digits <- c(digits,digits,digits)

   if (length(digits) == 2L)
      digits <- c(digits,digits[[1]])

   ddd <- list(...)

   ############################################################################

   ### set default colors if user has not specified 'col' and 'border' arguments

   if (x$int.only) {

      if (missing(col)) {
         col <- c(par("fg"), .coladj(par("fg"), dark=-0.3, light=0.3)) # 1st = summary polygon, 2nd = PI
      } else {
         if (length(col) == 1L) # if user only specified one value, assume it is for the summary polygon
            col <- c(col, .coladj(col, dark=-0.3, light=0.3))
      }

      if (missing(border))
         border <- par("fg") # border color of summary polygon

   } else {

      if (missing(col))
         col <- .coladj(par("bg","fg"), dark=0.2, light=-0.2) # color for fitted value polygons

      if (missing(border))
         border <- .coladj(par("bg","fg"), dark=0.3, light=-0.3)

   }

   ### set default line types if user has not specified 'lty' argument

   if (missing(lty)) {
      lty <- c("solid", "dotted", "solid") # 1st = CIs, 2nd = PI, 3rd = horizontal line(s)
   } else {
      if (length(lty) == 1L)
         lty <- c(lty, "dotted", "solid")
      if (length(lty) == 2L)
         lty <- c(lty, "solid")
   }

   ### vertical expansion factor: 1st = CI/PI end lines, 2nd = arrows, 3rd = summary polygon or fitted polygons

   if (length(efac) == 1L)
      efac <- rep(efac, 3L)

   if (length(efac) == 2L)
      efac <- c(efac[1], efac[1], efac[2]) # if 2 values specified: 1st = CI/PI end lines and arrows, 2nd = summary polygon or fitted polygons

   efac[efac == 0] <- NA

   ### annotation symbols vector

   if (is.null(ddd$annosym)) {
      annosym <- c(" [", ", ", "]", "-", " ") # 4th element for minus sign symbol; 5th for space (in place of numbers and +); see [a]
   } else {
      annosym <- ddd$annosym
      if (length(annosym) == 3L)
         annosym <- c(annosym, "-", " ")
      if (length(annosym) == 4L)
         annosym <- c(annosym, " ")
      if (length(annosym) != 5L)
         stop(mstyle$stop("Argument 'annosym' must be a vector of length 3 (or 4 or 5)."))
   }

   ### adjust annosym for tabular figures

   if (isTRUE(ddd$tabfig == 1))
      annosym <- c("\u2009[", ",\u2009", "]", "\u2212", "\u2002") # \u2009 thin space; \u2212 minus, \u2002 en space
   if (isTRUE(ddd$tabfig == 2))
      annosym <- c("\u2009[", ",\u2009", "]", "\u2013", "\u2002") # \u2009 thin space; \u2013 en dash, \u2002 en space
   if (isTRUE(ddd$tabfig == 3))
      annosym <- c("\u2009[", ",\u2009", "]", "\u2212", "\u2007") # \u2009 thin space; \u2212 minus, \u2007 figure space

   ### get measure from object

   measure <- x$measure

   ### column header

   estlab <- .setlab(measure, transf.char, atransf.char, gentype=3, short=TRUE)
   if (is.expression(estlab)) {
      header.right <- str2lang(paste0("bold(", estlab, " * '", annosym[1], "' * '", 100*(1-level), "% CI'", " * '", annosym[3], "')"))
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

   if (!is.null(ddd$addcred))
      addpred <- ddd$addcred

   pi.type <- .chkddd(ddd$pi.type, "default")

   decreasing <- .chkddd(ddd$decreasing, FALSE)

   if (!is.null(ddd$clim))
      olim <- ddd$clim

   ### row adjustments for 1) study labels, 2) annotations, and 3) ilab elements

   if (is.null(ddd$rowadj)) {
      rowadj <- rep(0,3)
   } else {
      rowadj <- ddd$rowadj
      if (length(rowadj) == 1L)
         rowadj <- c(rowadj,rowadj,0) # if one value is specified, use it for both 1&2
      if (length(rowadj) == 2L)
         rowadj <- c(rowadj,0) # if two values are specified, use them for 1&2
   }

   top <- .chkddd(ddd$top, 3)

   if (is.null(ddd$xlabadj)) {
      xlabadj <- c(NA,NA)
   } else {
      xlabadj <- ddd$xlabadj
      if (length(xlabadj) == 1L)
         xlabadj <- c(xlabadj, 1-xlabadj)
   }

   xlabfont <- .chkddd(ddd$xlabfont, 1)

   lplot     <- function(..., textpos, addcred, pi.type, decreasing, clim, rowadj, annosym, tabfig, top, xlabadj, xlabfont, at.lab) plot(...)
   labline   <- function(..., textpos, addcred, pi.type, decreasing, clim, rowadj, annosym, tabfig, top, xlabadj, xlabfont, at.lab) abline(...)
   lsegments <- function(..., textpos, addcred, pi.type, decreasing, clim, rowadj, annosym, tabfig, top, xlabadj, xlabfont, at.lab) segments(...)
   laxis     <- function(..., textpos, addcred, pi.type, decreasing, clim, rowadj, annosym, tabfig, top, xlabadj, xlabfont, at.lab) axis(...)
   lmtext    <- function(..., textpos, addcred, pi.type, decreasing, clim, rowadj, annosym, tabfig, top, xlabadj, xlabfont, at.lab) mtext(...)
   lpolygon  <- function(..., textpos, addcred, pi.type, decreasing, clim, rowadj, annosym, tabfig, top, xlabadj, xlabfont, at.lab) polygon(...)
   ltext     <- function(..., textpos, addcred, pi.type, decreasing, clim, rowadj, annosym, tabfig, top, xlabadj, xlabfont, at.lab) text(...)
   lpoints   <- function(..., textpos, addcred, pi.type, decreasing, clim, rowadj, annosym, tabfig, top, xlabadj, xlabfont, at.lab) points(...)

   if (is.character(showweights)) {
      weighttype  <- match.arg(showweights, c("diagonal", "rowsum"))
      if (weighttype == "rowsum" && !inherits(x, "rma.mv"))
         weighttype <- "diagonal"
      if (weighttype == "rowsum" && !x$int.only)
         stop(mstyle$stop("Row-sum weights are only meaningful for intercept-only models."))
      showweights <- TRUE
   } else {
      weighttype <- "diagonal"
   }

   if (!is.logical(showweights))
      stop(mstyle$stop("Argument 'showweights' must be a logical."))

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

   slab.null <- FALSE

   if (missing(slab)) {

      if (x$slab.null) {
         slab <- paste("Study", x$ids)          # x$ids is always of length yi.f (i.e., NAs also have an id)
         slab.null <- TRUE
      } else {
         slab <- x$slab                         # x$slab is always of length yi.f (i.e., NAs also have a study label)
      }

   } else {

      slab <- .getx("slab", mf=mf, data=x$data)

      if (length(slab) == 1L && is.na(slab)) {  # slab=NA can be used to suppress study labels
         slab <- rep("", x$k.all)
         slab.null <- TRUE
      }

      if (length(slab) != x$k.all)
         stop(mstyle$stop(paste0("Length of the 'slab' argument (", length(slab), ") does not correspond to the size of the original dataset (", x$k.all, ").")))

      slab <- .getsubset(slab, x$subset)

   }

   if (!is.null(ilab)) {

      if (is.null(dim(ilab)))
         ilab <- cbind(ilab)

      if (nrow(ilab) != x$k.all)
         stop(mstyle$stop(paste0("Length of the 'ilab' argument (", nrow(ilab), ") does not correspond to the size of the original dataset (", x$k.all, ").")))

      ilab <- .getsubset(ilab, x$subset)

   }

   if (length(pch) == 1L)
      pch <- rep(pch, x$k.all)

   if (length(pch) != x$k.all)
      stop(mstyle$stop(paste0("Length of the 'pch' argument (", length(pch), ") does not correspond to the size of the original dataset (", x$k.all, ").")))

   pch <- .getsubset(pch, x$subset)

   if (!is.null(psize)) {

      if (length(psize) == 1L)
         psize <- rep(psize, x$k.all)

      if (length(psize) != x$k.all)
         stop(mstyle$stop(paste0("Length of the 'psize' argument (", length(psize), ") does not correspond to the size of the original dataset (", x$k.all, ").")))

      psize <- .getsubset(psize, x$subset)

   }

   if (length(colout) == 1L)
      colout <- rep(colout, x$k.all)

   if (length(colout) != x$k.all)
      stop(mstyle$stop(paste0("Length of the 'colout' argument (", length(colout), ") does not correspond to the size of the original dataset (", x$k.all, ").")))

   colout <- .getsubset(colout, x$subset)

   shade.type <- "none"

   if (is.character(shade)) {

      shade.type <- "character"
      shade <- shade[1]

      if (!is.element(shade, c("zebra", "zebra1", "zebra2", "all")))
         stop(mstyle$stop("Unknown option specified for 'shade' argument."))

   }

   if (is.logical(shade)) {

      if (length(shade) == 1L) {
         shade <- "zebra"
         shade.type <- "character"
      } else {
         shade.type <- "logical"
         shade <- .chksubset(shade, x$k.all, stoponk0=FALSE)
         shade <- .getsubset(shade, x$subset)
      }

   }

   if (is.numeric(shade))
      shade.type <- "numeric"

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

      weights <- try(weights(x, type=weighttype), silent=TRUE) # does not work for rma.glmm and rma.uni.selmodel objects

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

         if (grepl("^order\\(", deparse1(substitute(order)))) {
            sort.vec <- order
         } else {
            sort.vec <- order(order, decreasing=decreasing)
         }

         if (!is.null(x$subset))
            sort.vec <- .getsubset(sort.vec, x$subset) - sum(!x$subset)

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
      if (shade.type == "logical")
         shade <- shade[sort.vec]

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
   if (shade.type == "logical")
      shade <- shade[k:1]

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
         for (j in seq_along(rows.na)) {
            rows.new[rows >= rows.na[j]] <- rows.new[rows >= rows.na[j]] - 1
         }
         rows <- rows.new[not.na]

         if (shade.type == "logical")
            shade <- shade[not.na]

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

   if (!is.null(at)) {
      if (anyNA(at))
         stop(mstyle$stop("Argument 'at' cannot contain NAs."))
      if (any(is.infinite(at)))
         stop(mstyle$stop("Argument 'at' cannot contain +-Inf values."))
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

   alim <- sort(alim)[1:2]

   if (anyNA(alim))
      stop(mstyle$stop("Argument 'alim' cannot contain NAs."))

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

   if (is.null(ddd$at.lab)) {

      at.lab <- at

      if (is.function(atransf)) {
         if (is.null(targs)) {
            at.lab <- fmtx(sapply(at.lab, atransf), digits[[2]], drop0ifint=TRUE)
         } else {
            at.lab <- fmtx(sapply(at.lab, atransf, targs), digits[[2]], drop0ifint=TRUE)
         }
      } else {
         at.lab <- fmtx(at.lab, digits[[2]], drop0ifint=TRUE)
      }

   } else {

      at.lab <- ddd$at.lab

   }

   ### set plot limits (xlim)

   ncol.ilab <- ifelse(is.null(ilab), 0, ncol(ilab))

   if (slab.null) {
      area.slab <- 25
   } else {
      area.slab <- 40
   }

   if (annotate) {
      if (showweights) {
         area.anno <- 30
      } else {
         area.anno <- 25
      }
   } else {
      area.anno <- 10
   }

   iadd <- 5

   area.slab   <- area.slab + iadd*ncol.ilab
   #area.anno  <- area.anno
   area.forest <- 100 + iadd*ncol.ilab - area.slab - area.anno

   area.slab   <- area.slab   / (100 + iadd*ncol.ilab)
   area.anno   <- area.anno   / (100 + iadd*ncol.ilab)
   area.forest <- area.forest / (100 + iadd*ncol.ilab)

   plot.multp.l <- area.slab / area.forest
   plot.multp.r <- area.anno / area.forest

   if (missing(xlim)) {
      if (min(ci.lb, na.rm=TRUE) < alim[1]) {
         f.1 <- alim[1]
      } else {
         f.1 <- min(ci.lb, na.rm=TRUE)
      }
      if (max(ci.ub, na.rm=TRUE) > alim[2]) {
         f.2 <- alim[2]
      } else {
         f.2 <- max(ci.ub, na.rm=TRUE)
      }
      rng <- f.2 - f.1
      xlim <- c(f.1 - rng * plot.multp.l, f.2 + rng * plot.multp.r)
      xlim <- round(xlim, digits[[2]])
      #xlim[1] <- xlim[1]*max(1, digits[[2]]/2)
      #xlim[2] <- xlim[2]*max(1, digits[[2]]/2)
   }

   xlim <- sort(xlim)

   ### plot limits must always encompass the yi values (no longer done)

   #if (xlim[1] > min(yi, na.rm=TRUE)) { xlim[1] <- min(yi, na.rm=TRUE) }
   #if (xlim[2] < max(yi, na.rm=TRUE)) { xlim[2] <- max(yi, na.rm=TRUE) }

   ### x-axis limits must always encompass the yi values (no longer done)

   #if (alim[1] > min(yi, na.rm=TRUE)) { alim[1] <- min(yi, na.rm=TRUE) }
   #if (alim[2] < max(yi, na.rm=TRUE)) { alim[2] <- max(yi, na.rm=TRUE) }

   ### plot limits must always encompass the x-axis limits (no longer done)

   #if (alim[1] < xlim[1]) { xlim[1] <- alim[1] }
   #if (alim[2] > xlim[2]) { xlim[2] <- alim[2] }

   ### allow adjustment of position of study labels and annotations via textpos argument

   textpos <- .chkddd(ddd$textpos, xlim)

   if (length(textpos) != 2L)
      stop(mstyle$stop("Argument 'textpos' must be of length 2."))

   if (is.na(textpos[1]))
      textpos[1] <- xlim[1]

   if (is.na(textpos[2]))
      textpos[2] <- xlim[2]

   ### set y-axis limits

   if (missing(ylim)) {
      if (x$int.only && addfit) {
        if (missing(y)) {
          ylim <- c(-1.5, max(rows, na.rm=TRUE)+top)
        } else {
          ylim <- c(-2.5, max(rows, na.rm=TRUE)+top)
        }
      } else {
         ylim <- c(0.5, max(rows, na.rm=TRUE)+top)
      }
   } else {
      ylim <- sort(ylim)
   }

   #########################################################################

   ### set/get fonts (1st for study labels, 2nd for annotations, 3rd for ilab)
   ### when passing a named vector, the names are for 'family' and the values are for 'font'

   if (missing(fonts)) {
      fonts <- rep(par("family"), 3L)
   } else {
      if (length(fonts) == 1L)
         fonts <- rep(fonts, 3L)
      if (length(fonts) == 2L)
         fonts <- c(fonts, fonts[1])
   }

   if (is.null(names(fonts)))
      fonts <- setNames(c(1L,1L,1L), nm=fonts)

   par(family=names(fonts)[1], font=fonts[1])

   ### adjust margins

   par.mar <- par("mar")
   par.mar.adj <- par.mar - c(0,3,1,1)
   par.mar.adj[par.mar.adj < 0] <- 0
   par(mar = par.mar.adj)
   on.exit(par(mar = par.mar), add=TRUE)

   #if (identical(par("mar"), c(5.1,4.1,4.1,2.1)))
   #   par(mar = c(5.1,1.1,3.1,1.1))

   ### start plot

   lplot(NA, NA, xlim=xlim, ylim=ylim, xlab="", ylab="", yaxt="n", xaxt="n", xaxs="i", bty="n", ...)

   ### add shading

   if (shade.type == "character") {
      if (shade == "zebra" || shade == "zebra1")
         tmp <- rep_len(c(TRUE,FALSE), k)
      if (shade == "zebra2")
         tmp <- rep_len(c(FALSE,TRUE), k)
      if (shade == "all")
         tmp <- rep_len(TRUE, k)
      shade <- tmp
   }

   if (shade.type %in% c("character","logical")) {
      for (i in seq_len(k)) {
         if (shade[i])
            rect(xlim[1], rows[i]-0.5, xlim[2], rows[i]+0.5, border=colshade, col=colshade)
      }
   }

   if (shade.type == "numeric") {
      for (i in seq_along(shade)) {
         rect(xlim[1], shade[i]-0.5, xlim[2], shade[i]+0.5, border=colshade, col=colshade)
      }
   }

   ### horizontal title line

   labline(h=ylim[2]-(top-1), lty=lty[3], ...)

   ### get coordinates of the plotting region

   par.usr <- par("usr")

   ### add reference line

   if (is.numeric(refline))
      lsegments(refline, par.usr[3], refline, ylim[2]-(top-1), lty="dotted", ...)

   ### set cex, cex.lab, and cex.axis sizes as a function of the height of the figure

   height <- par.usr[4] - par.usr[3]

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
            ### for multiple tau^2 (and gamma^2) values, need to specify level(s) of the inner factor(s) to compute the PI
            ### this can be done via the addpred argument (i.e., instead of using a logical, one specifies the level(s))
            if (length(addpred) == 1L)
               addpred <- c(addpred, addpred)
            temp <- predict(x, level=level, tau2.levels=addpred[1], gamma2.levels=addpred[2], pi.type=pi.type)
            addpred <- TRUE ### set addpred to TRUE, so if (!is.element(x$method, c("FE","EE","CE")) && addpred) further below works
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

      if (!is.element(x$method, c("FE","EE","CE")) && addpred) {

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
         mlab <- sapply(x$method, switch, "FE"="FE Model", "EE"="EE Model", "CE"="CE Model", "RE Model", USE.NAMES=FALSE)

      if (is.list(mlab)) {
         ltext(textpos[1], -1+rowadj[1], mlab[[1]], pos=4, cex=cex, ...)
      } else {
         ltext(textpos[1], -1+rowadj[1], mlab, pos=4, cex=cex, ...)
      }
      
      ### add polygon and label for model estimate for other fixed/random effects model 
      
      if (!missing(y)) {
        
        # obtain estimates for other model
        temp <- predict(y, level=level, pi.type=pi.type)
        beta.y       <- temp$pred
        beta.ci.lb.y <- temp$ci.lb
        beta.ci.ub.y <- temp$ci.ub
        if (is.function(transf)) {
          if (is.null(targs)) {
            beta.y       <- sapply(beta.y, transf)
            beta.ci.lb.y <- sapply(beta.ci.lb.y, transf)
            beta.ci.ub.y <- sapply(beta.ci.ub.y, transf)
          } else {
            beta.y       <- sapply(beta.y, transf, targs)
            beta.ci.lb.y <- sapply(beta.ci.lb.y, transf, targs)
            beta.ci.ub.y <- sapply(beta.ci.ub.y, transf, targs)
          }
        }
        # make sure order of intervals is always increasing
        tmp <- .psort(beta.ci.lb.y, beta.ci.ub.y)
        beta.ci.lb.y <- tmp[,1]
        beta.ci.ub.y <- tmp[,2]
        # apply observation/outcome limits if specified
        if (!missing(olim)) {
          beta.ci.lb.y[beta.ci.lb.y < olim[1]] <- olim[1]
          beta.ci.ub.y[beta.ci.ub.y > olim[2]] <- olim[2]
        }
        
        # add polygon #
        lpolygon(x=c(beta.ci.lb.y, beta.y, beta.ci.ub.y, beta.y), y=c(-2, -2+(height/100)*cex*efac[3], -2, -2-(height/100)*cex*efac[3]), col=col[1], border=border, ...)
        
        # add labels
        if (missing(mlab.y))
          mlab.y <- sapply(y$method, switch, "FE"="FE Model", "EE"="EE Model", "CE"="CE Model", "RE Model", USE.NAMES=FALSE)
        if (is.list(mlab.y)) {
          ltext(ddd$textpos[1], -2, mlab.y[[1]], pos=4, cex=cex, ...)
        } else {
          ltext(ddd$textpos[1], -2, mlab.y, pos=4, cex=cex, ...)
        }
      }

   }

   #########################################################################

   ### add x-axis

   laxis(side=1, at=at, labels=at.lab, cex.axis=cex.axis, ...)

   ### add x-axis label

   if (missing(xlab))
      xlab <- .setlab(measure, transf.char, atransf.char, gentype=1)

   if (!is.element(length(xlab), 1:3))
      stop(mstyle$stop("Argument 'xlab' argument must be of length 1, 2, or 3."))

   if (length(xlab) == 1L)
      lmtext(xlab, side=1, at=min(at) + (max(at)-min(at))/2, line=par("mgp")[1]-0.5, cex=cex.lab, font=xlabfont[1], ...)
   if (length(xlab) == 2L) {
      lmtext(xlab[1], side=1, at=min(at), line=par("mgp")[1]-0.5, cex=cex.lab, adj=xlabadj[1], font=xlabfont[1], ...)
      lmtext(xlab[2], side=1, at=max(at), line=par("mgp")[1]-0.5, cex=cex.lab, adj=xlabadj[2], font=xlabfont[1], ...)
   }
   if (length(xlab) == 3L) {
      lmtext(xlab[1], side=1, at=min(at), line=par("mgp")[1]-0.5, cex=cex.lab, adj=xlabadj[1], font=xlabfont[1], ...)
      lmtext(xlab[2], side=1, at=min(at) + (max(at)-min(at))/2, line=par("mgp")[1]-0.5, cex=cex.lab, font=xlabfont[2], ...)
      lmtext(xlab[3], side=1, at=max(at), line=par("mgp")[1]-0.5, cex=cex.lab, adj=xlabadj[2], font=xlabfont[1], ...)
   }

   ### add CI ends (either | or <> if outside of axis limits)

   for (i in seq_len(k)) {

      ### need to skip missings (if check below will otherwise throw an error)
      if (is.na(yi[i]) || is.na(vi[i]))
         next

      ### if the lower bound is actually larger than upper x-axis limit, then everything is to the right and just draw a polygon pointing in that direction
      if (ci.lb[i] >= alim[2]) {
         lpolygon(x=c(alim[2], alim[2]-(1.4/100)*cex*(xlim[2]-xlim[1]), alim[2]-(1.4/100)*cex*(xlim[2]-xlim[1]), alim[2]), y=c(rows[i], rows[i]+(height/150)*cex*efac[2], rows[i]-(height/150)*cex*efac[2], rows[i]), col=colout[i], border=colout[i], ...)
         next
      }

      ### if the upper bound is actually lower than lower x-axis limit, then everything is to the left and just draw a polygon pointing in that direction
      if (ci.ub[i] <= alim[1]) {
         lpolygon(x=c(alim[1], alim[1]+(1.4/100)*cex*(xlim[2]-xlim[1]), alim[1]+(1.4/100)*cex*(xlim[2]-xlim[1]), alim[1]), y=c(rows[i], rows[i]+(height/150)*cex*efac[2], rows[i]-(height/150)*cex*efac[2], rows[i]), col=colout[i], border=colout[i], ...)
         next
      }

      lsegments(max(ci.lb[i], alim[1]), rows[i], min(ci.ub[i], alim[2]), rows[i], lty=lty[1], col=colout[i], ...)

      if (ci.lb[i] >= alim[1]) {
         lsegments(ci.lb[i], rows[i]-(height/150)*cex*efac[1], ci.lb[i], rows[i]+(height/150)*cex*efac[1], col=colout[i], ...)
      } else {
         lpolygon(x=c(alim[1], alim[1]+(1.4/100)*cex*(xlim[2]-xlim[1]), alim[1]+(1.4/100)*cex*(xlim[2]-xlim[1]), alim[1]), y=c(rows[i], rows[i]+(height/150)*cex*efac[2], rows[i]-(height/150)*cex*efac[2], rows[i]), col=colout[i], border=colout[i], ...)
      }

      if (ci.ub[i] <= alim[2]) {
         lsegments(ci.ub[i], rows[i]-(height/150)*cex*efac[1], ci.ub[i], rows[i]+(height/150)*cex*efac[1], col=colout[i], ...)
      } else {
         lpolygon(x=c(alim[2], alim[2]-(1.4/100)*cex*(xlim[2]-xlim[1]), alim[2]-(1.4/100)*cex*(xlim[2]-xlim[1]), alim[2]), y=c(rows[i], rows[i]+(height/150)*cex*efac[2], rows[i]-(height/150)*cex*efac[2], rows[i]), col=colout[i], border=colout[i], ...)
      }

   }

   ### add study labels on the left

   ltext(textpos[1], rows+rowadj[1], slab, pos=4, cex=cex, ...)

   ### add info labels

   if (!is.null(ilab)) {
      if (is.null(ilab.xpos)) {
         #stop(mstyle$stop("Must specify 'ilab.xpos' argument when adding information with 'ilab'."))
         dist <- min(ci.lb, na.rm=TRUE) - xlim[1]
         if (ncol.ilab == 1L)
            ilab.xpos <- xlim[1] + dist*0.75
         if (ncol.ilab == 2L)
            ilab.xpos <- xlim[1] + dist*c(0.65, 0.85)
         if (ncol.ilab == 3L)
            ilab.xpos <- xlim[1] + dist*c(0.60, 0.75, 0.90)
         if (ncol.ilab >= 4L)
            ilab.xpos <- seq(xlim[1] + dist*0.5, xlim[1] + dist*0.9, length.out=ncol.ilab)
      }
      if (length(ilab.xpos) != ncol.ilab)
         stop(mstyle$stop(paste0("Number of 'ilab' columns (", ncol.ilab, ") does not match length of 'ilab.xpos' argument (", length(ilab.xpos), ").")))
      if (!is.null(ilab.pos) && length(ilab.pos) == 1L)
         ilab.pos <- rep(ilab.pos, ncol.ilab)
      par(family=names(fonts)[3], font=fonts[3])
      for (l in seq_len(ncol.ilab)) {
         ltext(ilab.xpos[l], rows+rowadj[3], ilab[,l], pos=ilab.pos[l], cex=cex, ...)
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
              if (missing(y)) {
                annotext <- cbind(sapply(c(yi, beta), atransf), sapply(c(ci.lb, beta.ci.lb), atransf), sapply(c(ci.ub, beta.ci.ub), atransf))
              } else {
                annotext <- cbind(sapply(c(yi, beta, beta.y), atransf), sapply(c(ci.lb, beta.ci.lb, beta.ci.lb.y), atransf), sapply(c(ci.ub, beta.ci.ub, beta.ci.ub.y), atransf))
              }
            } else {
               annotext <- cbind(sapply(yi, atransf), sapply(ci.lb, atransf), sapply(ci.ub, atransf))
            }
         } else {
            if (addfit && x$int.only) {
              if (missing(y)) {
                annotext <- cbind(sapply(c(yi, beta), atransf, targs), sapply(c(ci.lb, beta.ci.lb), atransf, targs), sapply(c(ci.ub, beta.ci.ub), atransf, targs))
              } else {
                annotext <- cbind(sapply(c(yi, beta, beta.y), atransf, targs), sapply(c(ci.lb, beta.ci.lb, beta.ci.lb.y), atransf, targs), sapply(c(ci.ub, beta.ci.ub, beta.ci.ub.y), atransf, targs))
              }
            } else {
               annotext <- cbind(sapply(yi, atransf, targs), sapply(ci.lb, atransf, targs), sapply(ci.ub, atransf, targs))
            }
         }

         ### make sure order of intervals is always increasing

         tmp <- .psort(annotext[,2:3])
         annotext[,2:3] <- tmp

      } else {

         if (addfit && x$int.only) {
           if (missing(y)) {
             annotext <- cbind(c(yi, beta), c(ci.lb, beta.ci.lb), c(ci.ub, beta.ci.ub))
           } else {
             annotext <- cbind(c(yi, beta, beta.y), c(ci.lb, beta.ci.lb, beta.ci.lb.y), c(ci.ub, beta.ci.ub, beta.ci.ub.y))
           }
         } else {
            annotext <- cbind(yi, ci.lb, ci.ub)
         }

      }

      if (showweights) {
         if (addfit && x$int.only) {
            annotext <- cbind(c(unname(weights),100), annotext)
            annotext <- fmtx(annotext, c(digits[[3]], digits[[1]], digits[[1]], digits[[1]]))
            annotext[nrow(annotext),1] <- "100"
         } else {
            annotext <- cbind(unname(weights), annotext)
            annotext <- fmtx(annotext, c(digits[[3]], digits[[1]], digits[[1]], digits[[1]]))
         }
      } else {
         annotext <- fmtx(annotext, digits[[1]])
      }

      if (missing(width)) {
         width <- apply(annotext, 2, function(x) max(nchar(x)))
      } else {
         if (length(width) == 1L)
            width <- rep(width, ncol(annotext))
         if (length(width) != ncol(annotext))
            stop(mstyle$stop(paste0("Length of 'width' argument (", length(width), ") does not match the number of annotation columns (", ncol(annotext), ").")))
      }

      for (j in seq_len(ncol(annotext))) {
         annotext[,j] <- formatC(annotext[,j], width=width[j])
      }

      if (showweights)
         width <- width[-1] # remove the first entry for the weights (so this can be used by addpoly() via .metafor)

      if (showweights) {
         annotext <- cbind(annotext[,1], paste0("%", paste0(rep(substr(annosym[1],1,1),3), collapse="")), annotext[,2], annosym[1], annotext[,3], annosym[2], annotext[,4], annosym[3])
      } else {
         annotext <- cbind(annotext[,1], annosym[1], annotext[,2], annosym[2], annotext[,3], annosym[3])
      }

      annotext <- apply(annotext, 1, paste, collapse="")
      annotext[grepl("NA", annotext, fixed=TRUE)] <- ""
      annotext <- gsub("-", annosym[4], annotext, fixed=TRUE) # [a]
      annotext <- gsub(" ", annosym[5], annotext, fixed=TRUE)

      par(family=names(fonts)[2], font=fonts[2])

      if (addfit && x$int.only) {
        if (missing(y)) {
          ltext(ddd$textpos[2], c(rows,-1), labels=annotext, pos=2, cex=cex, ...)
        } else {
          ltext(ddd$textpos[2], c(rows,-1,-2), labels=annotext, pos=2, cex=cex, ...)
        }
      } else {
         ltext(textpos[2], rows+rowadj[2], labels=annotext, pos=2, cex=cex, ...)
      }

      par(family=names(fonts)[1], font=fonts[1])

   } else {
      width <- NULL
   }

   ### add yi points

   for (i in seq_len(k)) {

      ### need to skip missings, as if () check below will otherwise throw an error
      if (is.na(yi[i]))
         next

      if (yi[i] >= alim[1] && yi[i] <= alim[2])
         lpoints(x=yi[i], y=rows[i], pch=pch[i], col=colout[i], cex=cex*psize[i], ...)

   }

   #lpoints(x=yi, y=rows, pch=pch, cex=cex*psize, ...)

   ### add horizontal line at 0 for the standard FE/RE model display

   if (x$int.only && addfit)
      labline(h=0, lty=lty[3], ...)

   ### add header

   ltext(textpos[1], ylim[2]-(top-1)+1, header.left,  pos=4, font=2, cex=cex, ...)
   ltext(textpos[2], ylim[2]-(top-1)+1, header.right, pos=2, font=2, cex=cex, ...)

   #########################################################################

   ### return some information about plot invisibly

   res <- list(xlim=par("usr")[1:2], alim=alim, at=at, ylim=ylim, rows=rows, cex=cex, cex.lab=cex.lab, cex.axis=cex.axis, ilab.xpos=ilab.xpos, ilab.pos=ilab.pos, textpos=textpos, areas=c(area.slab, area.forest, area.anno))

   ### put stuff into the .metafor environment, so that it can be used by addpoly()

   sav <- c(res, list(level=level, annotate=annotate, digits=digits[[1]], width=width, transf=transf, atransf=atransf, targs=targs, efac=efac, fonts=fonts[1:2], annosym=annosym))
   try(assign("forest", sav, envir=.metafor), silent=TRUE)

   invisible(res)

}
