labbe.rma <- function(x, xlim, ylim, lim, xlab, ylab, flip=FALSE,
ci=FALSE, pi=FALSE, grid=FALSE, legend=FALSE,
add=x$add, to=x$to, transf, targs, pch=21, psize, plim=c(0.5,3.5), col, bg, lty, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(x), must="rma", notav=c("rma.mv", "rma.ls", "rma.gen", "rma.uni.selmodel"))

   if (!x$int.only)
      stop(mstyle$stop("L'Abbe plots can only be drawn for models without moderators."))

   if (!is.element(x$measure, c("RR","OR","RD","AS","IRR","IRD","IRSD")))
      stop(mstyle$stop("Argument 'measure' must have been set to one of the following: 'RR','OR','RD','AS','IRR','IRD','IRSD'."))

   if (is.null(x$outdat.f))
      stop(mstyle$stop("Information needed to construct the plot is not available in the model object."))

   na.act <- getOption("na.action")
   on.exit(options(na.action=na.act), add=TRUE)

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   if (length(add) == 2L) # for rma.mh and rma.peto objects (1st 'add' value applies to the individual outcomes)
      add <- add[1]

   if (length(to) == 2L)  # for rma.mh and rma.peto objects (1st 'to' value applies to the individual outcomes)
      to <- to[1]

   if (!is.element(to, c("all","only0","if0all","none")))
      stop(mstyle$stop("Unknown 'to' argument specified."))

   .start.plot()

   if (missing(transf))
      transf <- FALSE

   transf.char <- deparse(transf)

   if (missing(targs))
      targs <- NULL

   if (missing(psize))
      psize <- NULL

   if (missing(lty)) {
      lty <- c("solid", "dashed") # 1 = diagonal line, 2 = estimated effect line
   } else {
      if (length(lty) == 1L)
         lty <- c(lty, lty)
   }

   if (is.logical(ci))
      cicol <- .coladj(par("bg","fg"), dark=0.15, light=-0.15)

   if (is.character(ci)) {
      cicol <- ci
      ci <- TRUE
   }

   if (is.logical(pi))
      picol <- .coladj(par("bg","fg"), dark=0.05, light=-0.05)

   if (is.character(pi)) {
      picol <- pi
      pi <- TRUE
   }

   ### get ... argument

   ddd <- list(...)

   ### set defaults or get addyi and addvi arguments

   addyi <- .chkddd(ddd$addyi, TRUE)
   addvi <- .chkddd(ddd$addvi, TRUE)

   ### grid argument can either be a logical or a color

   if (is.logical(grid))
      gridcol <- .coladj(par("bg","fg"), dark=c(0.2,-0.6), light=c(-0.2,0.6))

   if (is.character(grid)) {
      gridcol <- grid
      grid <- TRUE
   }

   llim <- ddd$llim

   lplot     <- function(..., addyi, addvi, llim) plot(...)
   lbox      <- function(..., addyi, addvi, llim) box(...)
   lsegments <- function(..., addyi, addvi, llim) segments(...)
   llines    <- function(..., addyi, addvi, llim) lines(...)
   lpoints   <- function(..., addyi, addvi, llim) points(...)
   lpolygon  <- function(..., addyi, addvi, llim) polygon(...)

   #########################################################################

   ### note: pch, psize, col, and bg (if vectors) must be of the same length as the original dataset
   ###       so we have to apply the same subsetting (if necessary) and removing of NAs as was
   ###       done during the model fitting (note: NAs are removed further below)

   pch <- .expand1(pch, x$k.all)

   if (length(pch) != x$k.all)
      stop(mstyle$stop(paste0("Length of the 'pch' argument (", length(pch), ") does not correspond to the size of the original dataset (", x$k.all, ").")))

   pch <- .getsubset(pch, x$subset)

   ### if user has set the point sizes

   if (!is.null(psize)) {
      psize <- .expand1(psize, x$k.all)
      if (length(psize) != x$k.all)
         stop(mstyle$stop(paste0("Length of the 'psize' argument (", length(psize), ") does not correspond to the size of the original dataset (", x$k.all, ").")))
      psize <- .getsubset(psize, x$subset)
   }

   if (missing(col))
      col <- par("fg")

   col <- .expand1(col, x$k.all)

   if (length(col) != x$k.all)
      stop(mstyle$stop(paste0("Length of the 'col' argument (", length(col), ") does not correspond to the size of the original dataset (", x$k.all, ").")))

   col <- .getsubset(col, x$subset)

   if (missing(bg))
      bg <- .coladj(par("bg","fg"), dark=0.35, light=-0.35)

   bg <- .expand1(bg, x$k.all)

   if (length(bg) != x$k.all)
      stop(mstyle$stop(paste0("Length of the 'bg' argument (", length(bg), ") does not correspond to the size of the original dataset (", x$k.all, ").")))

   bg <- .getsubset(bg, x$subset)

   #########################################################################

   ### these vectors may contain NAs

   dat.ai  <- x$outdat.f$ai
   dat.bi  <- x$outdat.f$bi
   dat.ci  <- x$outdat.f$ci
   dat.di  <- x$outdat.f$di
   dat.x1i <- x$outdat.f$x1i
   dat.x2i <- x$outdat.f$x2i
   dat.y1i <- x$outdat.f$t1i
   dat.y2i <- x$outdat.f$t2i

   ### drop00=TRUE may induce that the contrast-based yi value is NA; so
   ### make sure that the corresponding arm-based yi values are also NA

   yi.is.na <- is.na(x$yi.f)
   dat.ai[yi.is.na]  <- NA_real_
   dat.bi[yi.is.na]  <- NA_real_
   dat.ci[yi.is.na]  <- NA_real_
   dat.di[yi.is.na]  <- NA_real_
   dat.x1i[yi.is.na] <- NA_real_
   dat.x2i[yi.is.na] <- NA_real_
   dat.y1i[yi.is.na] <- NA_real_
   dat.y2i[yi.is.na] <- NA_real_

   options(na.action = "na.pass") # to make sure dat.x and dat.y are of the same length

   measure <- switch(x$measure, "RR"="PLN", "OR"="PLO", "RD"="PR", "AS"="PAS", "IRR"="IRLN", "IRD"="IR", "IRSD"="IRS")

   if (is.element(x$measure, c("RR","OR","RD","AS"))) {
      if (flip) {
         args.x <- list(measure=measure, xi=dat.ai, mi=dat.bi, add=add, to=to, addyi=addyi, addvi=addvi)
         args.y <- list(measure=measure, xi=dat.ci, mi=dat.di, add=add, to=to, addyi=addyi, addvi=addvi)
      } else {
         args.x <- list(measure=measure, xi=dat.ci, mi=dat.di, add=add, to=to, addyi=addyi, addvi=addvi)
         args.y <- list(measure=measure, xi=dat.ai, mi=dat.bi, add=add, to=to, addyi=addyi, addvi=addvi)
      }
   }

   if (is.element(x$measure, c("IRR","IRD","IRSD"))) {
      if (flip) {
         args.x <- list(measure=measure, xi=dat.x1i, ti=dat.y1i, add=add, to=to, addyi=addyi, addvi=addvi)
         args.y <- list(measure=measure, xi=dat.x2i, ti=dat.y2i, add=add, to=to, addyi=addyi, addvi=addvi)
      } else {
         args.x <- list(measure=measure, xi=dat.x2i, ti=dat.y2i, add=add, to=to, addyi=addyi, addvi=addvi)
         args.y <- list(measure=measure, xi=dat.x1i, ti=dat.y1i, add=add, to=to, addyi=addyi, addvi=addvi)
      }
   }

   dat.x <- .do.call(escalc, args.x)
   dat.y <- .do.call(escalc, args.y)

   options(na.action = na.act)

   ### check for NAs in yi/vi pairs and filter out

   has.na <- apply(is.na(dat.x), 1, any) | apply(is.na(dat.y), 1, any)
   not.na <- !has.na

   if (any(has.na)) {

      dat.x <- dat.x[not.na,]
      dat.y <- dat.y[not.na,]
      pch   <- pch[not.na]
      col   <- col[not.na]
      bg    <- bg[not.na]

      if (is.null(psize))
         psize <- psize[not.na]

   }

   if (length(dat.x$yi)==0L || length(dat.y$yi)==0L)
      stop(mstyle$stop("No information in object to compute the arm-level outcomes."))

   #########################################################################

   ### determine point sizes

   vi <- dat.x$vi + dat.y$vi

   k <- length(vi)

   if (is.null(psize)) {
      if (length(plim) < 2L)
         stop(mstyle$stop("Argument 'plim' must be of length 2 or 3."))
      wi <- sqrt(1/vi)
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

   ### determine x/y values for line that indicates the estimated effect

   min.yi <- min(c(dat.x$yi, dat.y$yi))
   max.yi <- max(c(dat.x$yi, dat.y$yi))
   rng.yi <- max.yi - min.yi

   len <- 10000

   intrcpt <- x$beta[1]

   if (is.null(llim)) {

      if (x$measure == "RD")
         x.vals <- seq(ifelse(intrcpt>0, 0, -intrcpt), ifelse(intrcpt>0, 1-intrcpt, 1), length.out=len)
      if (x$measure == "RR")
         x.vals <- seq(min.yi-rng.yi, ifelse(intrcpt>0, -intrcpt, 0), length.out=len)
      if (x$measure == "OR")
         x.vals <- seq(min.yi-rng.yi, max.yi+rng.yi, length.out=len)
      if (x$measure == "AS")
         x.vals <- seq(ifelse(intrcpt>0, 0, -intrcpt), ifelse(intrcpt>0, asin(sqrt(1))-intrcpt, asin(sqrt(1))), length.out=len)
      if (x$measure == "IRR")
         x.vals <- seq(min.yi-rng.yi, ifelse(intrcpt>0, -intrcpt, 0), length.out=len)
      if (x$measure == "IRD")
         x.vals <- seq(ifelse(intrcpt>0, 0, -intrcpt), ifelse(intrcpt>0, 1-intrcpt, 1), length.out=len)
      if (x$measure == "IRSD")
         x.vals <- seq(ifelse(intrcpt>0, 0, -intrcpt), ifelse(intrcpt>0, 1-intrcpt, 1), length.out=len)

   } else {

      if (length(llim) != 2L)
         stop(mstyle$stop("Argument 'llim' must be of length 2."))

      x.vals <- seq(llim[1], llim[2], length.out=len)

   }

   y.vals <- intrcpt + 1*x.vals

   if (ci || pi) {

      predres <- predict(x)

      y.vals.ci.lb <- predres$ci.lb + 1*x.vals
      y.vals.ci.ub <- predres$ci.ub + 1*x.vals
      y.vals.pi.lb <- predres$pi.lb + 1*x.vals
      y.vals.pi.ub <- predres$pi.ub + 1*x.vals

   } else {

      y.vals.ci.lb <- y.vals.ci.ub <- y.vals.pi.lb <- y.vals.pi.ub <- NULL

   }

   if (is.function(transf)) {
      if (is.null(targs)) {
         dat.x$yi     <- sapply(dat.x$yi, transf)
         dat.y$yi     <- sapply(dat.y$yi, transf)
         x.vals       <- sapply(x.vals, transf)
         y.vals       <- sapply(y.vals, transf)
         y.vals.ci.lb <- sapply(y.vals.ci.lb, transf)
         y.vals.ci.ub <- sapply(y.vals.ci.ub, transf)
         y.vals.pi.lb <- sapply(y.vals.pi.lb, transf)
         y.vals.pi.ub <- sapply(y.vals.pi.ub, transf)
      } else {
         if (!is.primitive(transf) && !is.null(targs) && length(formals(transf)) == 1L)
            stop(mstyle$stop("Function specified via 'transf' does not appear to have an argument for 'targs'."))
         dat.x$yi     <- sapply(dat.x$yi, transf, targs)
         dat.y$yi     <- sapply(dat.y$yi, transf, targs)
         x.vals       <- sapply(x.vals, transf, targs)
         y.vals       <- sapply(y.vals, transf, targs)
         y.vals.ci.lb <- sapply(y.vals.ci.lb, transf, targs)
         y.vals.ci.ub <- sapply(y.vals.ci.ub, transf, targs)
         y.vals.pi.lb <- sapply(y.vals.pi.lb, transf, targs)
         y.vals.pi.ub <- sapply(y.vals.pi.ub, transf, targs)
      }
   }

   min.yi <- min(c(dat.x$yi, dat.y$yi))
   max.yi <- max(c(dat.x$yi, dat.y$yi))

   if (missing(lim)) {

      if (missing(xlim))
         xlim <- c(min.yi, max.yi)

      if (missing(ylim))
         ylim <- c(min.yi, max.yi)

   } else {

      xlim <- lim
      ylim <- lim

   }

   ### order points by psize

   order.vec <- order(psize, decreasing=TRUE)

   dat.x$yi.o  <- dat.x$yi[order.vec]
   dat.y$yi.o  <- dat.y$yi[order.vec]
   pch.o       <- pch[order.vec]
   col.o       <- col[order.vec]
   bg.o        <- bg[order.vec]
   psize.o     <- psize[order.vec]

   ### add x-axis label

   if (missing(xlab)) {
      xlab <- .setlab(measure, transf.char, atransf.char="FALSE", gentype=1)
      xlab <- paste(xlab, "(Group 1)")
   }

   ### add y-axis label

   if (missing(ylab)) {
      ylab <- .setlab(measure, transf.char, atransf.char="FALSE", gentype=1)
      ylab <- paste(ylab, "(Group 2)")
   }

   lplot(NA, NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)

   ### add PI bounds

   if (pi)
      lpolygon(c(x.vals,rev(x.vals)), c(y.vals.pi.lb,rev(y.vals.pi.ub)), col=picol, border=NA, ...)

   ### add CI bounds

   if (ci)
      lpolygon(c(x.vals,rev(x.vals)), c(y.vals.ci.lb,rev(y.vals.ci.ub)), col=cicol, border=NA, ...)

   ### add grid (and redraw box)

   if (.isTRUE(grid)) {
      grid(col=gridcol)
      lbox(...)
   }

   ### add diagonal reference line

   #abline(a=0, b=1, lty=lty[1], ...)
   lsegments(min(x.vals), min(x.vals), max(x.vals), max(x.vals), lty=lty[1], ...)

   ### add estimated effects line

   llines(x.vals, y.vals, lty=lty[2], ...)

   ### add points

   lpoints(x=dat.x$yi.o, y=dat.y$yi.o, cex=psize.o, pch=pch.o, col=col.o, bg=bg.o, ...)

   ### add legend

   lopts <- list(x      = ifelse(intrcpt > 0, "bottomright", "topleft"),
                 y      = NULL,
                 inset  = 0.01,
                 cex    = 1,
                 pt.cex = 2.5)

   if (is.list(legend)) {

      # replace defaults with any user-defined values
      lopts.pos <- pmatch(names(legend), names(lopts))
      lopts[c(na.omit(lopts.pos))] <- legend[!is.na(lopts.pos)]

      # rescale pt.cex based on cex (if pt.cex was not specified)
      if (!is.null(legend$cex) && is.null(legend$pt.cex))
         lopts$pt.cex <- lopts$pt.cex * lopts$cex

      legend <- TRUE

   } else {

      if (is.character(legend)) {
         lopts$x <- legend
         legend <- TRUE
      } else {
         if (!is.logical(legend))
            stop(mstyle$stop("Argument 'legend' must either be logical, a string, or a list."), call.=FALSE)
      }

   }

   if (legend) {

      lvl <- round(100*(1-x$level), x$digits[["ci"]])
      ltxt <- c("Reference Line of No Effect", "Line for the Estimated Effect",
                paste0(lvl, "% Confidence Interval"), paste0(lvl, "% Prediction Interval"))
      lpch <- c(NA,NA,22,22)
      if (is.numeric(lty)) {
         llty <- c(lty[1],lty[2],0,0)
      } else {
         llty <- c(lty[1],lty[2],"blank","blank")
      }
      lpt.bg <- c(NA,NA,cicol,picol)

      sel <- c(lty != "blank" & lty != 0, ci, pi)

      if (any(sel)) {
         legend(x=lopts$x, y=lopts$y, inset=lopts$inset, bg=.coladj(par("bg"), dark=0, light=0),
                pch=lpch[sel], pt.cex=lopts$pt.cex, pt.lwd=0, pt.bg=lpt.bg[sel],
                lty=llty[sel], legend=ltxt[sel], cex=lopts$cex)
      }

   }

   #########################################################################

   ### prepare data frame to return

   sav <- data.frame(x=dat.x$yi, y=dat.y$yi, cex=psize, pch=pch, col=col, bg=bg, ids=x$ids[not.na], slab=x$slab[not.na])

   invisible(sav)

}
