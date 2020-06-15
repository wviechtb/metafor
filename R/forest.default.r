forest.default <- function(x, vi, sei, ci.lb, ci.ub,
annotate=TRUE,                             showweights=FALSE, header=FALSE,
xlim, alim, clim, ylim, top=3, at, steps=5, level=95,      refline=0, digits=2L, width,
xlab, slab,       ilab, ilab.xpos, ilab.pos, subset,
transf, atransf, targs, rows,
efac=1, pch=15, psize, col,                 lty, fonts,
cex, cex.lab, cex.axis, annosym, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   na.act <- getOption("na.action")

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

   yi <- x

   ### note: transf and atransf must be function names and cannot, for example, be arguments
   ### passed down from other functions (i.e., deparse(substitute(...)) will grab exactly what
   ### is specified for the argument), so the following function would not work:
   ###
   ### misc <- function(x, vi, tfunction=FALSE)
   ###   forest.default(x, vi, atransf=tfunction)

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

   if (missing(subset))
      subset <- NULL

   if (missing(psize))
      psize <- NULL

   if (missing(col))
      col <- NULL

   if (missing(cex))
      cex <- NULL

   if (missing(cex.lab))
      cex.lab <- NULL

   if (missing(cex.axis))
      cex.axis <- NULL

   if (missing(lty)) {
      lty <- c("solid", "solid") ### 1st value = CIs, 2nd value = horizontal line(s)
   } else {
      if (length(lty) == 1L)
         lty <- c(lty, "solid")
   }

   ### vertical expansion factor: 1st = CI end lines, 2nd = arrows

   if (length(efac) == 1L)
      efac <- rep(efac, 2)

   ### annotation symbols vector

   if (missing(annosym))
      annosym <- c(" [", ", ", "]")
   if (length(annosym) != 3L)
      stop(mstyle$stop("Argument 'annosym' must be a vector of length 3."))

   level <- ifelse(level == 0, 1, ifelse(level >= 1, (100-level)/100, ifelse(level > .5, 1-level, level)))

   ### set measure based on the measure attribute of yi

   if (is.null(attr(yi, "measure"))) {
      measure <- "GEN"
   } else {
      measure <- attr(yi, "measure")
   }

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

   lplot     <- function(..., textpos) plot(...)
   labline   <- function(..., textpos) abline(...)
   lsegments <- function(..., textpos) segments(...)
   laxis     <- function(..., textpos) axis(...)
   lmtext    <- function(..., textpos) mtext(...)
   lpolygon  <- function(..., textpos) polygon(...)
   ltext     <- function(..., textpos) text(...)
   lpoints   <- function(..., textpos) points(...)

   if (!is.null(ddd$order))
      stop(mstyle$stop("Function does not have an 'order' argument (could use 'subset' argument instead)."))

   #########################################################################

   ### digits[1] for annotations, digits[2] for x-axis labels
   ### note: digits can also be a list (e.g., digits=list(2L,3))

   if (length(digits) == 1L)
      digits <- c(digits,digits)

   if (hasArg(ci.lb) && hasArg(ci.ub)) {     ### CI bounds are specified by user
      if (length(ci.lb) != length(ci.ub))
         stop(mstyle$stop("Length of 'ci.lb' and 'ci.ub' is not the same."))
      if (missing(vi) && missing(sei)) {     ### vi/sei not specified, so calculate vi based on CI
         vi <- ((ci.ub - ci.lb) / (2*qnorm(level/2, lower.tail=FALSE)))^2
      } else {
         if (missing(vi))                    ### vi not specified, but sei is, so set vi = sei^2
            vi <- sei^2
      }
      if (length(ci.lb) != length(vi))
         stop(mstyle$stop("Length of 'vi' (or 'sei') does not match length of ('ci.lb', 'ci.ub') pairs."))
   } else {                                  ### CI bounds are not specified by user
      if (missing(vi)) {
         if (missing(sei)) {
            stop(mstyle$stop("Must specify either 'vi', 'sei', or ('ci.lb', 'ci.ub') pairs."))
         } else {
            vi <- sei^2
            ci.lb <- yi - qnorm(level/2, lower.tail=FALSE) * sei
            ci.ub <- yi + qnorm(level/2, lower.tail=FALSE) * sei
         }
      } else {
         ci.lb <- yi - qnorm(level/2, lower.tail=FALSE) * sqrt(vi)
         ci.ub <- yi + qnorm(level/2, lower.tail=FALSE) * sqrt(vi)
      }
   }

   ### check length of yi and vi

   if (length(yi) != length(vi))
      stop(mstyle$stop("Length of 'yi' does not match the length of 'vi', 'sei', or the ('ci.lb', 'ci.ub') pairs."))

   k <- length(yi)

   if (missing(slab)) {                         ### note: slab must have same length as yi (including NAs)
      if (!is.null(attr(yi, "slab"))) {         ### even when subsetting eventually
         slab <- attr(yi, "slab")               ### check if slab info can be found in slab attribute of yi
      } else {
         slab <- paste("Study", seq_len(k))
      }
   } else {
      if (length(slab) == 1L && is.na(slab))
         slab <- rep("", k)
   }

   if (length(yi) != length(slab))
      stop(mstyle$stop("Number of outcomes does not correspond to the length of the 'slab' argument."))

   if (is.null(dim(ilab)))                      ### note: ilab must have same length as yi (including NAs)
      ilab <- cbind(ilab)                       ### even when subsetting eventually

   if (length(pch) == 1L)                       ### note: pch must have same length as yi (including NAs)
      pch <- rep(pch, k)                        ### or be equal to a single value (which is then repeated)

   if (length(pch) != length(yi))
      stop(mstyle$stop("Number of outcomes does not correspond to the length of the 'pch' argument."))

   ### if user has set the point sizes

   if (!is.null(psize)) {                       ### note: psize must have same length as yi (including NAs)
      if (length(psize) == 1L)                  ### or be equal to a single value (which is then repeated)
         psize <- rep(psize, k)
      if (length(psize) != length(yi))
         stop(mstyle$stop("Number of outcomes does not correspond to the length of the 'psize' argument."))
   }

   ### if user has set the col argument

   if (!is.null(col)) {                         ### note: col must have same length as yi (including NAs)
      if (length(col) == 1L)                    ### or be equal to a single value (which is then repeated)
         col <- rep(col, k)
      if (length(col) != length(yi))
         stop(mstyle$stop("Number of outcomes does not correspond to the length of the 'col' argument."))
   } else {
      col <- rep("black", k)
   }

   ### if a subset of studies is specified (can also be used to reorder the data)

   if (!is.null(subset)) {
      yi    <- yi[subset]
      vi    <- vi[subset]
      ci.lb <- ci.lb[subset]
      ci.ub <- ci.ub[subset]
      slab  <- slab[subset]
      ilab  <- ilab[subset,,drop=FALSE]         ### if ilab is still NULL, then this remains NULL
      pch   <- pch[subset]
      psize <- psize[subset]                    ### if psize is still NULL, then this remains NULL
      col   <- col[subset]
   }

   k <- length(yi)                              ### in case length of k has changed

   ### set rows value

   if (missing(rows)) {
      rows <- k:1
   } else {
      if (length(rows) == 1L)                   ### note: rows must be a single value or the same
         rows <- rows:(rows-k+1)                ### length of yi (including NAs) *after ordering/subsetting*
   }

   if (length(rows) != length(yi))
      stop(mstyle$stop("Number of outcomes does not correspond to the length of the 'rows' argument."))

   ### reverse order

   yi    <- yi[k:1]
   vi    <- vi[k:1]
   ci.lb <- ci.lb[k:1]
   ci.ub <- ci.ub[k:1]
   slab  <- slab[k:1]
   ilab  <- ilab[k:1,,drop=FALSE]               ### if ilab is still NULL, then this remains NULL
   pch   <- pch[k:1]
   psize <- psize[k:1]                          ### if psize is still NULL, then this remains NULL
   col   <- col[k:1]
   rows  <- rows[k:1]

   ### check for NAs in yi/vi and act accordingly

   yivi.na <- is.na(yi) | is.na(vi)

   if (any(yivi.na)) {

      not.na <- !yivi.na

      if (na.act == "na.omit") {
         yi    <- yi[not.na]
         vi    <- vi[not.na]
         ci.lb <- ci.lb[not.na]
         ci.ub <- ci.ub[not.na]
         slab  <- slab[not.na]
         ilab  <- ilab[not.na,,drop=FALSE]      ### if ilab is still NULL, then this remains NULL
         pch   <- pch[not.na]
         psize <- psize[not.na]                 ### if psize is still NULL, then this remains NULL
         col   <- col[not.na]

         rows.new <- rows                       ### rearrange rows due to NAs being omitted from plot
         rows.na  <- rows[!not.na]              ### shift higher rows down according to number of NAs omitted
         for (j in seq_len(length(rows.na))) {
            rows.new[rows >= rows.na[j]] <- rows.new[rows >= rows.na[j]] - 1
         }
         rows <- rows.new[not.na]

      }

      if (na.act == "na.fail")
         stop(mstyle$stop("Missing values in results."))

   }                                            ### note: yi/vi may be NA if na.act == "na.exclude" or "na.pass"

   k <- length(yi)                              ### in case length of k has changed

   ### if requested, apply transformation to yi's and CI bounds

   if (is.function(transf)) {
      if (is.null(targs)) {
         yi         <- sapply(yi, transf)
         ci.lb      <- sapply(ci.lb, transf)
         ci.ub      <- sapply(ci.ub, transf)
      } else {
         yi         <- sapply(yi, transf, targs)
         ci.lb      <- sapply(ci.lb, transf, targs)
         ci.ub      <- sapply(ci.ub, transf, targs)
      }
   }

   ### make sure order of intervals is always increasing

   tmp <- .psort(ci.lb, ci.ub)
   ci.lb <- tmp[,1]
   ci.ub <- tmp[,2]

   ### apply ci limits if specified

   if (!missing(clim)) {
      clim <- sort(clim)
      if (length(clim) != 2L)
         stop(mstyle$stop("Argument 'clim' must be of length 2."))
      ci.lb[ci.lb < clim[1]] <- clim[1]
      ci.ub[ci.ub > clim[2]] <- clim[2]
   }

   if (showweights) {                           ### inverse variance weights after ordering/subsetting and
      weights <- 1/vi                           ### omitting NAs (so these weights always add up to 100%)
      weights <- 100 * weights/sum(weights, na.rm=TRUE)
   }

   ### set default point sizes (if not specified by user)

   if (is.null(psize)) {
      if (any(vi <= 0, na.rm=TRUE)) {           ### in case any vi value is zero
         psize <- rep(1, k)
      } else {                                  ### default psize is proportional to inverse standard error
         wi    <- 1/sqrt(vi)                    ### note: vi's that are NA are ignored (but vi's whose yi is
         psize <- wi/sum(wi, na.rm=TRUE)        ### NA are NOT ignored; an unlikely case in practice)
         rng   <- max(psize, na.rm=TRUE) - min(psize, na.rm=TRUE)
         if (rng <= .Machine$double.eps^0.5) {
            psize <- rep(1, k)
         } else {
            psize <- (psize - min(psize, na.rm=TRUE)) / rng
            psize <- (psize * 1.0) + 0.5        ### note: only vi's that are still in the subset are used for determining the default point sizes
         }
         if (all(is.na(psize)))                 ### if k=1, then psize is NA, so catch this (and maybe some other problems)
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
      xlim <- round(xlim, digits[[2]])
      #xlim[1] <- xlim[1]*max(1, digits[[2]]/2)
      #xlim[2] <- xlim[2]*max(1, digits[[2]]/2)
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

   ### allow adjustment of position of study labels and annotations via textpos argument

   if (is.null(ddd$textpos))
      ddd$textpos <- c(xlim[1], xlim[2])

   if (length(ddd$textpos) != 2L)
      stop(mstyle$stop("Argument 'textpos' must be of length 2."))

   if (is.na(ddd$textpos[1]))
      ddd$textpos[1] <- xlim[1]

   if (is.na(ddd$textpos[2]))
      ddd$textpos[2] <- xlim[2]

   ### set y axis limits

   if (missing(ylim)) {
      ylim <- c(0.5, k+top)
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
         at.lab <- formatC(sapply(at.lab, atransf), digits=digits[[2]], format="f", drop0trailing=ifelse(class(digits[[2]]) == "integer", TRUE, FALSE))
      } else {
         at.lab <- formatC(sapply(at.lab, atransf, targs), digits=digits[[2]], format="f", drop0trailing=ifelse(class(digits[[2]]) == "integer", TRUE, FALSE))
      }
   } else {
      at.lab <- formatC(at.lab, digits=digits[[2]], format="f", drop0trailing=ifelse(class(digits[[2]]) == "integer", TRUE, FALSE))
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
   on.exit(par(mar = par.mar))

   ### start plot

   lplot(NA, NA, xlim=xlim, ylim=ylim, xlab="", ylab="", yaxt="n", xaxt="n", xaxs="i", bty="n", col="black", ...)

   ### horizontal title line

   labline(h=ylim[2]-(top-1), lty=lty[2], col="black", ...)

   ### add reference line

   if (is.numeric(refline))
      lsegments(refline, ylim[1]-5, refline, ylim[2]-(top-1), lty="dotted", col="black", ...)

   ### set cex, cex.lab, and cex.axis sizes as a function of the height of the figure

   par.usr <- par("usr")
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

   ### add x axis

   laxis(side=1, at=at, labels=at.lab, cex.axis=cex.axis, col="black", ...)

   ### add x axis label

   if (missing(xlab))
      xlab <- .setlab(measure, transf.char, atransf.char, gentype=1)

   lmtext(xlab, side=1, at=min(at) + (max(at)-min(at))/2, line=par("mgp")[1]-0.5, cex=cex.lab, col="black", ...)

   ### add CI ends (either | or <> if outside of axis limits)

   for (i in seq_len(k)) {

      ### need to skip missings (if check below will otherwise throw an error)
      if (is.na(yi[i]) || is.na(ci.lb[i]) || is.na(ci.ub[i]))
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

   ltext(ddd$textpos[1], rows, slab, pos=4, cex=cex, col=col, ...)

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

      if (showweights)
         annotext <- cbind(weights, annotext)

      annotext <- formatC(annotext, format="f", digits=digits[[1]])

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
      ltext(ddd$textpos[2], rows, labels=annotext, pos=2, cex=cex, col=col, ...)
      par(family=names(fonts)[1], font=fonts[1])

   }

   ### add yi points

   for (i in seq_len(k)) {

      ### need to skip missings (if check below will otherwise throw an error)
      if (is.na(yi[i]))
         next

      if (yi[i] >= alim[1] && yi[i] <= alim[2])
         lpoints(yi[i], rows[i], pch=pch[i], cex=cex*psize[i], col=col[i], ...)

   }

   #lpoints(yi, rows, pch=pch, cex=cex*psize, ...)

   ### add header

   ltext(ddd$textpos[1], ylim[2]-(top-1)+1, header.left, pos=4, font=2, cex=cex, ...)
   ltext(ddd$textpos[2], ylim[2]-(top-1)+1, header.right, pos=2, font=2, cex=cex, ...)

   #########################################################################

   ### return some information about plot invisibly

   res <- list('xlim'=par("usr")[1:2], 'alim'=alim, 'at'=at, 'ylim'=ylim, 'rows'=rows, 'cex'=cex, 'cex.lab'=cex.lab, 'cex.axis'=cex.axis)

   invisible(res)

}
