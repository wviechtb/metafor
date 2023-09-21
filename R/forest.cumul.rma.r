forest.cumul.rma <- function(x,
annotate=TRUE,                                                header=FALSE,
xlim, alim, olim, ylim, at, steps=5, level=x$level, refline=0, digits=2L, width,
xlab,             ilab, ilab.xpos, ilab.pos,
transf, atransf, targs, rows,
efac=1, pch, psize,                           col,         shade, colshade,
lty, fonts, cex, cex.lab, cex.axis, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="cumul.rma")

   na.act <- getOption("na.action")

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

   yi <- x$estimate

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

   if (missing(col)) {
      col <- par("fg")
   } else {
      col <- .getx("col", mf=mf, data=x$data)
   }

   if (missing(pch)) {
      pch <- 15
   } else {
      pch <- .getx("pch", mf=mf, data=x$data)
   }

   if (missing(psize)) {
      psize <- 1
   } else {
      psize <- .getx("psize", mf=mf, data=x$data)
   }

   if (missing(shade)) {
      shade <- NULL
   } else {
      shade <- .getx("shade", mf=mf, data=x$data)
   }

   if (missing(colshade))
      colshade <- .coladj(par("bg","fg"), dark=0.1, light=-0.1)

   if (missing(cex))
      cex <- NULL

   if (missing(cex.lab))
      cex.lab <- NULL

   if (missing(cex.axis))
      cex.axis <- NULL

   level <- .level(level)

   ### digits[1] for annotations, digits[2] for x-axis labels
   ### note: digits can also be a list (e.g., digits=list(2,3L)); trailing 0's on the x-axis labels
   ### are dropped if the value is an integer

   if (length(digits) == 1L)
      digits <- c(digits,digits)

   ddd <- list(...)

   ############################################################################

   ### set default line types if user has not specified 'lty' argument

   if (missing(lty)) {
      lty <- c("solid", "solid") # 1st = CIs, 2nd = horizontal line(s)
   } else {
      if (length(lty) == 1L)
         lty <- c(lty, "solid")
   }

   ### vertical expansion factor: 1st = CI end lines, 2nd = arrows

   if (length(efac) == 1L)
      efac <- rep(efac, 2L)

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

   if (is.null(ddd$top)) {
      top <- 3
   } else {
      top <- ddd$top
   }

   if (is.null(ddd$xlabadj)) {
      xlabadj <- c(NA,NA)
   } else {
      xlabadj <- ddd$xlabadj
      if (length(xlabadj) == 1L)
         xlabadj <- c(xlabadj, 1-xlabadj)
   }

   if (is.null(ddd$xlabfont)) {
      xlabfont <- 1
   } else {
      xlabfont <- ddd$xlabfont
   }

   lplot     <- function(..., textpos, clim, rowadj, annosym, tabfig, top, xlabadj, xlabfont) plot(...)
   labline   <- function(..., textpos, clim, rowadj, annosym, tabfig, top, xlabadj, xlabfont) abline(...)
   lsegments <- function(..., textpos, clim, rowadj, annosym, tabfig, top, xlabadj, xlabfont) segments(...)
   laxis     <- function(..., textpos, clim, rowadj, annosym, tabfig, top, xlabadj, xlabfont) axis(...)
   lmtext    <- function(..., textpos, clim, rowadj, annosym, tabfig, top, xlabadj, xlabfont) mtext(...)
   lpolygon  <- function(..., textpos, clim, rowadj, annosym, tabfig, top, xlabadj, xlabfont) polygon(...)
   ltext     <- function(..., textpos, clim, rowadj, annosym, tabfig, top, xlabadj, xlabfont) text(...)
   lpoints   <- function(..., textpos, clim, rowadj, annosym, tabfig, top, xlabadj, xlabfont) points(...)

   #########################################################################

   ### extract data / results and other arguments

   vi <- x$se^2
   ci.lb <- x$ci.lb
   ci.ub <- x$ci.ub

   ### check length of yi and vi

   k <- length(yi) # either of length k when na.action="na.omit" or k.f otherwise

   if (length(vi) != k)
      stop(mstyle$stop("Length of 'yi' and 'vi' (or 'sei') is not the same."))

   ### note: ilab, pch, psize, col must be of the same length as yi (which may
   ###       or may not contain NAs; this is different than the other forest()
   ###       functions but it would be tricky to make this fully consistent now

   if (x$slab.null) {
      slab    <- paste("+ Study", x$ids)  # cumul() removes the studies with NAs when na.action="na.omit"
      slab[1] <- paste("Study", x$ids[1])
      slab.null <- TRUE
   } else {
      slab    <- paste("+", x$slab)       # cumul() removes the studies with NAs when na.action="na.omit"
      slab[1] <- paste(x$slab[1])
      slab.null <- FALSE
   }

   if (!is.null(ilab)) {

      if (is.null(dim(ilab)))
         ilab <- cbind(ilab)

      if (nrow(ilab) != k)
         stop(mstyle$stop(paste0("Length of the 'ilab' argument (", nrow(ilab), ") does not correspond to the number of outcomes (", k, ").")))

   }

   if (length(pch) == 1L)
      pch <- rep(pch, k)

   if (length(pch) != k)
      stop(mstyle$stop(paste0("Length of the 'pch' argument (", length(pch), ") does not correspond to the number of outcomes (", k, ").")))

   if (length(psize) == 1L)
      psize <- rep(psize, k)

   if (length(psize) != k)
      stop(mstyle$stop(paste0("Length of the 'psize' argument (", length(psize), ") does not correspond to the number of outcomes (", k, ").")))

   if (length(col) == 1L)
      col <- rep(col, k)

   if (length(col) != k)
      stop(mstyle$stop(paste0("Length of the 'col' argument (", length(col), ") does not correspond to the number of outcomes (", k, ").")))

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
         shade <- .chksubset(shade, k, stoponk0=FALSE)
      }

   }

   if (is.numeric(shade))
      shade.type <- "numeric"

   ### set rows value

   if (missing(rows)) {
      rows <- k:1
   } else {
      if (length(rows) == 1L)
         rows <- rows:(rows-k+1)
   }

   if (length(rows) != k)
      stop(mstyle$stop(paste0("Length of the 'rows' argument (", length(rows), ") does not correspond to the number of outcomes (", k, ").")))

   ### reverse order

   yi      <- yi[k:1]
   vi      <- vi[k:1]
   ci.lb   <- ci.lb[k:1]
   ci.ub   <- ci.ub[k:1]
   slab    <- slab[k:1]
   ilab    <- ilab[k:1,,drop=FALSE]               # if NULL, remains NULL
   pch     <- pch[k:1]
   psize   <- psize[k:1]                          # if NULL, remains NULL
   col     <- col[k:1]
   rows    <- rows[k:1]
   if (shade.type == "logical")
      shade <- shade[k:1]

   ### check for NAs in yi/vi and act accordingly

   yivi.na <- is.na(yi) | is.na(vi)

   if (any(yivi.na)) {

      not.na <- !yivi.na

      if (na.act == "na.omit") {
         yi      <- yi[not.na]
         vi      <- vi[not.na]
         ci.lb   <- ci.lb[not.na]
         ci.ub   <- ci.ub[not.na]
         slab    <- slab[not.na]
         ilab    <- ilab[not.na,,drop=FALSE]    # if NULL, remains NULL
         pch     <- pch[not.na]
         psize   <- psize[not.na]               # if NULL, remains NULL
         col     <- col[not.na]
         if (shade.type == "logical")
            shade <- shade[not.na]

         rows.new <- rows                       # rearrange rows due to NAs being omitted from plot
         rows.na  <- rows[!not.na]              # shift higher rows down according to number of NAs omitted
         for (j in seq_along(rows.na)) {
            rows.new[rows >= rows.na[j]] <- rows.new[rows >= rows.na[j]] - 1
         }
         rows <- rows.new[not.na]

      }

      if (na.act == "na.fail")
         stop(mstyle$stop("Missing values in results."))

   }                                            # note: yi/vi may be NA if na.act == "na.exclude" or "na.pass"

   k <- length(yi)                              # in case length of k has changed

   ### if requested, apply transformation to yi's and CI bounds

   if (is.function(transf)) {
      if (is.null(targs)) {
         yi    <- sapply(yi, transf)
         ci.lb <- sapply(ci.lb, transf)
         ci.ub <- sapply(ci.ub, transf)
      } else {
         yi    <- sapply(yi, transf, targs)
         ci.lb <- sapply(ci.lb, transf, targs)
         ci.ub <- sapply(ci.ub, transf, targs)
      }
   }

   ### make sure order of intervals is always increasing

   tmp <- .psort(ci.lb, ci.ub)
   ci.lb <- tmp[,1]
   ci.ub <- tmp[,2]

   ### apply observation/outcome limits if specified

   if (!missing(olim)) {
      if (length(olim) != 2L)
         stop(mstyle$stop("Argument 'olim' must be of length 2."))
      olim <- sort(olim)
      yi[yi < olim[1]] <- olim[1]
      yi[yi > olim[2]] <- olim[2]
      ci.lb[ci.lb < olim[1]] <- olim[1]
      ci.ub[ci.ub > olim[2]] <- olim[2]
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

   ### set plot limits (xlim)

   ncol.ilab <- ifelse(is.null(ilab), 0, ncol(ilab))

   if (slab.null) {
      area.slab <- 25
   } else {
      area.slab <- 40
   }

   if (annotate) {
      area.anno <- 25
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

   if (is.null(ddd$textpos)) {
      textpos <- xlim
   } else {
      textpos <- ddd$textpos
   }

   if (length(textpos) != 2L)
      stop(mstyle$stop("Argument 'textpos' must be of length 2."))

   if (is.na(textpos[1]))
      textpos[1] <- xlim[1]

   if (is.na(textpos[2]))
      textpos[2] <- xlim[2]

   ### set y-axis limits

   if (missing(ylim)) {
      ylim <- c(0.5, max(rows, na.rm=TRUE)+top)
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

   ### start plot

   lplot(NA, NA, xlim=xlim, ylim=ylim, xlab="", ylab="", yaxt="n", xaxt="n", xaxs="i", bty="n", ...)

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

   labline(h=ylim[2]-(top-1), lty=lty[2], ...)

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

   ### add x-axis

   laxis(side=1, at=at, labels=at.lab, cex.axis=cex.axis, ...)

   ### add x-axis label

   if (missing(xlab))
      xlab <- .setlab(measure, transf.char, atransf.char, gentype=2)

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
         lpolygon(x=c(alim[2], alim[2]-(1.4/100)*cex*(xlim[2]-xlim[1]), alim[2]-(1.4/100)*cex*(xlim[2]-xlim[1]), alim[2]), y=c(rows[i], rows[i]+(height/150)*cex*efac[2], rows[i]-(height/150)*cex*efac[2], rows[i]), col=col[i], border=col[i], ...)
         next
      }

      ### if the upper bound is actually lower than lower x-axis limit, then everything is to the left and just draw a polygon pointing in that direction
      if (ci.ub[i] <= alim[1]) {
         lpolygon(x=c(alim[1], alim[1]+(1.4/100)*cex*(xlim[2]-xlim[1]), alim[1]+(1.4/100)*cex*(xlim[2]-xlim[1]), alim[1]), y=c(rows[i], rows[i]+(height/150)*cex*efac[2], rows[i]-(height/150)*cex*efac[2], rows[i]), col=col[i], border=col[i], ...)
         next
      }

      lsegments(max(ci.lb[i], alim[1]), rows[i], min(ci.ub[i], alim[2]), rows[i], lty=lty[1], col=col[i], ...)

      if (ci.lb[i] >= alim[1]) {
         lsegments(ci.lb[i], rows[i]-(height/150)*cex*efac[1], ci.lb[i], rows[i]+(height/150)*cex*efac[1], col=col[i], ...)
      } else {
         lpolygon(x=c(alim[1], alim[1]+(1.4/100)*cex*(xlim[2]-xlim[1]), alim[1]+(1.4/100)*cex*(xlim[2]-xlim[1]), alim[1]), y=c(rows[i], rows[i]+(height/150)*cex*efac[2], rows[i]-(height/150)*cex*efac[2], rows[i]), col=col[i], border=col[i], ...)
      }

      if (ci.ub[i] <= alim[2]) {
         lsegments(ci.ub[i], rows[i]-(height/150)*cex*efac[1], ci.ub[i], rows[i]+(height/150)*cex*efac[1], col=col[i], ...)
      } else {
         lpolygon(x=c(alim[2], alim[2]-(1.4/100)*cex*(xlim[2]-xlim[1]), alim[2]-(1.4/100)*cex*(xlim[2]-xlim[1]), alim[2]), y=c(rows[i], rows[i]+(height/150)*cex*efac[2], rows[i]-(height/150)*cex*efac[2], rows[i]), col=col[i], border=col[i], ...)
      }

   }

   ### add study labels on the left

   ltext(textpos[1], rows+rowadj[1], slab, pos=4, cex=cex, col=col, ...)

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
      if (length(ilab.xpos) != ncol(ilab))
         stop(mstyle$stop(paste0("Number of 'ilab' columns (", ncol(ilab), ") does not match length of 'ilab.xpos' argument (", length(ilab.xpos), ").")))
      if (!is.null(ilab.pos) && length(ilab.pos) == 1L)
         ilab.pos <- rep(ilab.pos, ncol(ilab))
      par(family=names(fonts)[3], font=fonts[3])
      for (l in seq_len(ncol(ilab))) {
         ltext(ilab.xpos[l], rows+rowadj[3], ilab[,l], pos=ilab.pos[l], cex=cex, ...)
      }
      par(family=names(fonts)[1], font=fonts[1])
   }

   ### add study annotations on the right: yi [LB, UB]

   if (annotate) {

      if (is.function(atransf)) {
         if (is.null(targs)) {
            annotext <- cbind(sapply(yi, atransf), sapply(ci.lb, atransf), sapply(ci.ub, atransf))
         } else {
            annotext <- cbind(sapply(yi, atransf, targs), sapply(ci.lb, atransf, targs), sapply(ci.ub, atransf, targs))
         }

         ### make sure order of intervals is always increasing

         tmp <- .psort(annotext[,2:3])
         annotext[,2:3] <- tmp

      } else {

         annotext <- cbind(yi, ci.lb, ci.ub)

      }

      annotext <- fmtx(annotext, digits[[1]])

      if (missing(width)) {
         width <- apply(annotext, 2, function(x) max(nchar(x)))
      } else {
         if (length(width) == 1L)
            width <- rep(width, ncol(annotext))
         if (length(width) != ncol(annotext))
            stop(mstyle$stop(paste0("Length of 'width' argument (", length(width), ") does not the match number of annotation columns (", ncol(annotext), ").")))
      }

      for (j in seq_len(ncol(annotext))) {
         annotext[,j] <- formatC(annotext[,j], width=width[j])
      }

      annotext <- cbind(annotext[,1], annosym[1], annotext[,2], annosym[2], annotext[,3], annosym[3])

      annotext <- apply(annotext, 1, paste, collapse="")
      annotext[grepl("NA", annotext, fixed=TRUE)] <- ""
      annotext <- gsub("-", annosym[4], annotext, fixed=TRUE) # [a]
      annotext <- gsub(" ", annosym[5], annotext, fixed=TRUE)

      par(family=names(fonts)[2], font=fonts[2])
      ltext(textpos[2], rows+rowadj[2], labels=annotext, pos=2, cex=cex, col=col, ...)
      par(family=names(fonts)[1], font=fonts[1])

   } else {
      width <- NULL
   }

   ### add yi points

   for (i in seq_len(k)) {

      ### need to skip missings (if check below will otherwise throw an error)
      if (is.na(yi[i]))
         next

      if (yi[i] >= alim[1] && yi[i] <= alim[2])
         lpoints(x=yi[i], y=rows[i], pch=pch[i], cex=cex*psize[i], col=col[i], ...)

   }

   #lpoints(x=yi, y=rows, pch=pch, cex=cex*psize, ...)

   ### add header

   ltext(textpos[1], ylim[2]-(top-1)+1, header.left,  pos=4, font=2, cex=cex, ...)
   ltext(textpos[2], ylim[2]-(top-1)+1, header.right, pos=2, font=2, cex=cex, ...)

   #########################################################################

   ### return some information about plot invisibly

   res <- list(xlim=par("usr")[1:2], alim=alim, at=at, ylim=ylim, rows=rows, cex=cex, cex.lab=cex.lab, cex.axis=cex.axis, ilab.xpos=ilab.xpos, ilab.pos=ilab.pos, textpos=textpos)

   ### add some additional stuff to be put into .metafor environment, so that it can be used by addpoly()

   sav <- c(res, list(level=level, annotate=annotate, digits=digits[[1]], width=width, transf=transf, atransf=atransf, targs=targs, fonts=fonts[1:2], annosym=annosym))
   try(assign("forest", sav, envir=.metafor), silent=TRUE)

   invisible(res)

}
