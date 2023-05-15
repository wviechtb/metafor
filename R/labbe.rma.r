labbe.rma <- function(x, xlim, ylim, xlab, ylab,
add=x$add, to=x$to, transf, targs, pch=21, psize, plim=c(0.5,3.5), col, bg, grid=FALSE, lty, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="rma", notav=c("rma.ls", "rma.gen", "rma.uni.selmodel"))

   if (!x$int.only)
      stop(mstyle$stop("L'Abbe plots can only be drawn for models without moderators."))

   if (!is.element(x$measure, c("RR","OR","RD","AS","IRR","IRD","IRSD")))
      stop(mstyle$stop("Argument 'measure' must be set to one of the following: 'RR','OR','RD','AS','IRR','IRD','IRSD'."))

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

   if (exists(".darkplots") && .isTRUE(.darkplots))
      par(fg="gray95", bg="gray10", col="gray95", col.axis="gray95", col.lab="gray95", col.main="gray95", col.sub="gray95")

   if (missing(transf))
      transf <- FALSE

   transf.char <- deparse(transf)

   if (missing(targs))
      targs <- NULL

   if (missing(psize))
      psize <- NULL

   if (missing(lty)) {
      lty <- c("solid", "dashed") # 1st value = diagonal line, 2nd value = estimated effect line
   } else {
      if (length(lty) == 1L)
         lty <- c(lty, lty)
   }

   ### get ... argument

   ddd <- list(...)

   ### set defaults or get addyi and addvi arguments

   addyi <- ifelse(is.null(ddd$addyi), TRUE, ddd$addyi)
   addvi <- ifelse(is.null(ddd$addvi), TRUE, ddd$addvi)

   ### grid argument can either be a logical or a color

   if (is.logical(grid)) {
      if (.is.dark(par("bg"))) {
         gridcol <- "gray30"
      } else {
         gridcol <- "gray70"
      }
   }
   if (is.character(grid)) {
      gridcol <- grid
      grid <- TRUE
   }

   #########################################################################

   ### note: 'pch', 'psize', 'col', and 'bg' must be of the same length as the original data passed to rma()
   ###       so we have to apply the same subsetting (if necessary) and removing of NAs as done during the
   ###       model fitting (note: NAs are removed further below)

   if (length(pch) == 1L)
      pch <- rep(pch, x$k.all)

   if (length(pch) != x$k.all)
      stop(mstyle$stop(paste0("Length of the 'pch' argument (", length(pch), ") does not correspond to the size of the original dataset (", x$k.all, ").")))

   pch <- .getsubset(pch, x$subset)

   ### if user has set the point sizes

   if (!is.null(psize)) {
      if (length(psize) == 1L)
         psize <- rep(psize, x$k.all)
      if (length(psize) != x$k.all)
         stop(mstyle$stop(paste0("Length of the 'psize' argument (", length(psize), ") does not correspond to the size of the original dataset (", x$k.all, ").")))
      psize <- .getsubset(psize, x$subset)
   }

   if (missing(col))
      col <- par("fg")

   if (length(col) == 1L)
      col <- rep(col, x$k.all)

   if (length(col) != x$k.all)
      stop(mstyle$stop(paste0("Length of the 'col' argument (", length(col), ") does not correspond to the size of the original dataset (", x$k.all, ").")))

   col <- .getsubset(col, x$subset)

   if (missing(bg)) {
      if (.is.dark(par("bg"))) {
         bg <- "gray40"
      } else {
         bg <- "gray"
      }
   }

   if (length(bg) == 1L)
      bg <- rep(bg, x$k.all)

   if (length(bg) != x$k.all)
      stop(mstyle$stop(paste0("Length of the 'bg' argument (", length(bg), ") does not correspond to the size of the original dataset (", x$k.all, ").")))

   bg <- .getsubset(bg, x$subset)

   #########################################################################

   ### these vectors may contain NAs

   ai  <- x$outdat.f$ai
   bi  <- x$outdat.f$bi
   ci  <- x$outdat.f$ci
   di  <- x$outdat.f$di
   x1i <- x$outdat.f$x1i
   x2i <- x$outdat.f$x2i
   t1i <- x$outdat.f$t1i
   t2i <- x$outdat.f$t2i

   ### drop00=TRUE may induce that the contrast-based yi value is NA; so
   ### make sure that the corresponding arm-based yi values are also NA

   yi.is.na <- is.na(x$yi.f)
   ai[yi.is.na]  <- NA_real_
   bi[yi.is.na]  <- NA_real_
   ci[yi.is.na]  <- NA_real_
   di[yi.is.na]  <- NA_real_
   x1i[yi.is.na] <- NA_real_
   x2i[yi.is.na] <- NA_real_
   t1i[yi.is.na] <- NA_real_
   t2i[yi.is.na] <- NA_real_

   options(na.action = "na.pass") ### to make sure dat.t and dat.c are of the same length

   measure <- switch(x$measure, "RR"="PLN", "OR"="PLO", "RD"="PR", "AS"="PAS", "IRR"="IRLN", "IRD"="IR", "IRSD"="IRS")

   if (is.element(x$measure, c("RR","OR","RD","AS"))) {
      args.t <- list(measure=measure, xi=ai, mi=bi, add=add, to=to, addyi=addyi, addvi=addvi)
      args.c <- list(measure=measure, xi=ci, mi=di, add=add, to=to, addyi=addyi, addvi=addvi)
   }

   if (is.element(x$measure, c("IRR","IRD","IRSD"))) {
      args.t <- list(measure=measure, xi=x1i, ti=t1i, add=add, to=to, addyi=addyi, addvi=addvi)
      args.c <- list(measure=measure, xi=x2i, ti=t2i, add=add, to=to, addyi=addyi, addvi=addvi)
   }

   dat.t <- .do.call(escalc, args.t)
   dat.c <- .do.call(escalc, args.c)

   options(na.action = na.act)

   ### check for NAs in yi/vi pairs and filter out

   has.na <- apply(is.na(dat.t), 1, any) | apply(is.na(dat.c), 1, any)
   not.na <- !has.na

   if (any(has.na)) {

      dat.t <- dat.t[not.na,]
      dat.c <- dat.c[not.na,]
      pch   <- pch[not.na]
      col   <- col[not.na]
      bg    <- bg[not.na]

      if (is.null(psize))
         psize <- psize[not.na]

   }

   if (length(dat.t$yi)==0L || length(dat.c$yi)==0L)
      stop(mstyle$stop("No information in object to compute arm-level outcomes."))

   #########################################################################

   ### determine point sizes

   vi <- dat.t$vi + dat.c$vi

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

   min.yi <- min(c(dat.t$yi, dat.c$yi))
   max.yi <- max(c(dat.t$yi, dat.c$yi))
   rng.yi <- max.yi - min.yi

   len <- 1000

   intrcpt <- x$beta[1]

   if (x$measure == "RD")
      c.vals <- seq(ifelse(intrcpt>0, 0, -intrcpt), ifelse(intrcpt>0, 1-intrcpt, 1), length.out=len)
   if (x$measure == "RR")
      c.vals <- seq(min.yi-rng.yi, ifelse(intrcpt>0, 0-intrcpt, 0), length.out=len)
   if (x$measure == "OR")
      c.vals <- seq(min.yi-rng.yi, max.yi+rng.yi, length.out=len)
   if (x$measure == "AS")
      c.vals <- seq(ifelse(intrcpt>0, 0, -intrcpt), ifelse(intrcpt>0, asin(sqrt(1))-intrcpt, asin(sqrt(1))), length.out=len)
   if (x$measure == "IRR")
      c.vals <- seq(min.yi-rng.yi, ifelse(intrcpt>0, 0-intrcpt, 0), length.out=len)
   if (x$measure == "IRD")
      c.vals <- seq(ifelse(intrcpt>0, 0, -intrcpt), ifelse(intrcpt>0, 1-intrcpt, 1), length.out=len)
   if (x$measure == "IRSD")
      c.vals <- seq(ifelse(intrcpt>0, 0, -intrcpt), ifelse(intrcpt>0, 1-intrcpt, 1), length.out=len)

   t.vals <- intrcpt + 1*c.vals

   if (is.function(transf)) {
      if (is.null(targs)) {
         dat.t$yi <- sapply(dat.t$yi, transf)
         dat.c$yi <- sapply(dat.c$yi, transf)
         c.vals   <- sapply(c.vals, transf)
         t.vals   <- sapply(t.vals, transf)
      } else {
         dat.t$yi <- sapply(dat.t$yi, transf, targs)
         dat.c$yi <- sapply(dat.c$yi, transf, targs)
         c.vals   <- sapply(c.vals, transf, targs)
         t.vals   <- sapply(t.vals, transf, targs)
      }
   }

   min.yi <- min(c(dat.t$yi, dat.c$yi))
   max.yi <- max(c(dat.t$yi, dat.c$yi))

   if (missing(xlim))
      xlim <- c(min.yi, max.yi)

   if (missing(ylim))
      ylim <- c(min.yi, max.yi)

   ### order points by psize

   order.vec <- order(psize, decreasing=TRUE)

   dat.t$yi.o  <- dat.t$yi[order.vec]
   dat.c$yi.o  <- dat.c$yi[order.vec]
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

   plot(NA, NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)

   ### add grid (and redraw box)

   if (.isTRUE(grid)) {
      grid(col=gridcol)
      box(...)
   }

   ### add diagonal and estimated effects lines

   abline(a=0, b=1, lty=lty[1], ...)
   lines(c.vals, t.vals, lty=lty[2], ...)

   ### add points

   points(x=dat.c$yi.o, y=dat.t$yi.o, cex=psize.o, pch=pch.o, col=col.o, bg=bg.o, ...)

   #########################################################################

   ### prepare data frame to return

   sav <- data.frame(x=dat.c$yi, y=dat.t$yi, cex=psize, pch=pch, col=col, bg=bg, ids=x$ids[not.na], slab=x$slab[not.na], stringsAsFactors=FALSE)

   invisible(sav)

}
