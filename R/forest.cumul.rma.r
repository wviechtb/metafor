forest.cumul.rma <- function(x,          annotate=TRUE,
xlim, alim, clim, ylim, top=3, at, steps=5, level=x$level, refline=0, digits=2L, width,
xlab,             ilab, ilab.xpos, ilab.pos,
transf, atransf, targs, rows,
efac=1, pch=15, psize=1, col,       lty, fonts,
cex, cex.lab, cex.axis, annosym, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   if (!inherits(x, "cumul.rma"))
      stop(mstyle$stop("Argument 'x' must be an object of class \"cumul.rma\"."))

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

   if (missing(annosym))
      annosym <- c(" [", ", ", "]")
   if (length(annosym) != 3)
      stop(mstyle$stop("Argument 'annosym' must be a vector of length 3."))

   ddd <- list(...)

   lplot     <- function(..., textpos) plot(...)
   labline   <- function(..., textpos) abline(...)
   lsegments <- function(..., textpos) segments(...)
   laxis     <- function(..., textpos) axis(...)
   lmtext    <- function(..., textpos) mtext(...)
   lpolygon  <- function(..., textpos) polygon(...)
   ltext     <- function(..., textpos) text(...)
   lpoints   <- function(..., textpos) points(...)

   #########################################################################

   ### digits[1] for annotations, digits[2] for x-axis labels
   ### note: digits can also be a list (e.g., digits=list(2L,3))

   if (length(digits) == 1L)
      digits <- c(digits,digits)

   level <- ifelse(level == 0, 1, ifelse(level >= 1, (100-level)/100, ifelse(level > .5, 1-level, level)))

   yi <- x$estimate

   ### set measure based on the measure attribute of yi

   if (is.null(attr(yi, "measure"))) {
      measure <- "GEN"
   } else {
      measure <- attr(yi, "measure")
   }

   vi <- x$se^2
   ci.lb <- x$ci.lb
   ci.ub <- x$ci.ub

   ### check length of yi and vi

   if (length(yi) != length(vi))
      stop(mstyle$stop("Length of 'yi' and 'vi' (or 'sei') is not the same."))

   k <- length(yi)

   if (x$slab.null) {
      slab    <- paste("+ Study", x$slab)
      slab[1] <- paste("Study", x$slab[1])
   } else {
      slab    <- paste("+", x$slab)
      slab[1] <- paste(x$slab[1])
   }

   if (is.null(dim(ilab)))
      ilab <- cbind(ilab)

   if (length(pch) == 1L)
      pch <- rep(pch, k)
   if (length(pch) != length(yi))
      stop(mstyle$stop("Number of outcomes does not correspond to the length of the 'pch' argument."))

   if (length(psize) == 1L)
      psize <- rep(psize, k)
   if (length(psize) != length(yi))
      stop(mstyle$stop("Number of outcomes does not correspond to the length of the 'psize' argument."))

   ### if user has set the col argument

   if (!is.null(col)) {
      if (length(col) == 1L)
         col <- rep(col, k)
      if (length(col) != length(yi))
         stop(mstyle$stop("Number of outcomes does not correspond to the length of the 'col' argument."))
   } else {
      col <- rep("black", k)
   }

   ### set rows value

   if (missing(rows)) {
      rows <- k:1
   } else {
      if (length(rows) == 1L)
         rows <- rows:(rows-k+1)
   }

   if (length(rows) != length(yi))
      stop(mstyle$stop("Number of outcomes does not correspond to the length of the 'rows' argument."))

   ### reverse order

   yi      <- yi[k:1]
   vi      <- vi[k:1]
   ci.lb   <- ci.lb[k:1]
   ci.ub   <- ci.ub[k:1]
   slab    <- slab[k:1]
   ilab    <- ilab[k:1,,drop=FALSE]               ### if ilab is still NULL, then this remains NULL
   pch     <- pch[k:1]
   psize   <- psize[k:1]                          ### if psize is still NULL, then this remains NULL
   col     <- col[k:1]
   rows    <- rows[k:1]

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
         ilab    <- ilab[not.na,,drop=FALSE]    ### if ilab is still NULL, then this remains NULL
         pch     <- pch[not.na]
         psize   <- psize[not.na]               ### if psize is still NULL, then this remains NULL
         col     <- col[not.na]

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

   ### calculate individual CI bounds (skipped: CI bounds are already extracted above)

   #ci.lb <- yi - qnorm(level/2, lower.tail=FALSE) * sqrt(vi)
   #ci.ub <- yi + qnorm(level/2, lower.tail=FALSE) * sqrt(vi)

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

   #########################################################################

   ### total range of CI bounds

   rng <- max(ci.ub, na.rm=TRUE) - min(ci.lb, na.rm=TRUE)

   if (annotate) {
      plot.multp.l <- 1.20
      plot.multp.r <- 1.20
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

   if (length(ddd$textpos) != 2)
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

   if (missing(fonts)) {
      fonts <- rep(par("family"), 3)
   } else {
      if (length(fonts) == 1L)
         fonts <- rep(fonts, 3)
      if (length(fonts) == 2L)
         fonts <- c(fonts, fonts[1])
   }

   par(family=fonts[1])

   ### adjust margins

   par.mar <- par("mar")
   par.mar.adj <- par.mar - c(0,3,1,1)
   par.mar.adj[par.mar.adj < 0] <- 0
   par(mar = par.mar.adj)
   on.exit(par(mar = par.mar))

   ### start plot

   lplot(NA, NA, xlim=xlim, ylim=ylim, xlab="", ylab="", yaxt="n", xaxt="n", xaxs="i", bty="n", ...)

   ### horizontal title line

   labline(h=ylim[2]-(top-1), lty=lty[2], ...)

   ### add reference line

   if (is.numeric(refline))
      lsegments(refline, ylim[1]-5, refline, ylim[2]-(top-1), lty="dotted", ...)

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

   laxis(side=1, at=at, labels=at.lab, cex.axis=cex.axis, ...)

   ### add x axis label

   if (missing(xlab))
      xlab <- .setlab(measure, transf.char, atransf.char, gentype=2)

   lmtext(xlab, side=1, at=min(at) + (max(at)-min(at))/2, line=par("mgp")[1]-0.5, cex=cex.lab, ...)

   ### add CI ends (either | or <> if outside of axis limits)

   for (i in seq_len(k)) {

      ### need to skip missings, as if() check below will otherwise throw an error
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

   ltext(ddd$textpos[1], rows, slab, pos=4, cex=cex, col=col, ...)

   ### add info labels

   if (!is.null(ilab)) {
      if (is.null(ilab.xpos))
         stop(mstyle$stop("Must specify 'ilab.xpos' argument when adding information with 'ilab'."))
      if (length(ilab.xpos) != ncol(ilab))
         stop(mstyle$stop(paste0("Number of 'ilab' columns (", ncol(ilab), ") does not match length of 'ilab.xpos' argument (", length(ilab.xpos), ").")))
      if (!is.null(ilab.pos) && length(ilab.pos) == 1)
         ilab.pos <- rep(ilab.pos, ncol(ilab))
      par(family=fonts[3])
      for (l in seq_len(ncol(ilab))) {
         ltext(ilab.xpos[l], rows, ilab[,l], pos=ilab.pos[l], cex=cex, ...)
      }
      par(family=fonts[1])
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

      annotext <- cbind(annotext[,1], annosym[1], annotext[,2], annosym[2], annotext[,3], annosym[3])
      annotext <- apply(annotext, 1, paste, collapse="")
      annotext[grepl("NA", annotext, fixed=TRUE)] <- ""
      par(family=fonts[2])
      ltext(ddd$textpos[2], rows, labels=annotext, pos=2, cex=cex, col=col, ...)
      par(family=fonts[1])

   }

   ### add yi points

   for (i in seq_len(k)) {

      ### need to skip missings, as if() check below will otherwise throw an error
      if (is.na(yi[i]))
         next

      if (yi[i] >= alim[1] && yi[i] <= alim[2])
         lpoints(yi[i], rows[i], pch=pch[i], cex=cex*psize[i], col=col[i], ...)

   }

   #lpoints(yi, rows, pch=pch, cex=cex*psize, ...)

   #########################################################################

   ### return some information about plot invisibly

   res <- list('xlim'=par("usr")[1:2], 'alim'=alim, 'at'=at, 'ylim'=ylim, 'rows'=rows, 'cex'=cex, 'cex.lab'=cex.lab, 'cex.axis'=cex.axis)

   invisible(res)

}
