labbe.rma <- function(x, xlim, ylim, xlab, ylab,
add=x$add, to=x$to, transf, targs, pch=21, psize, bg="gray", grid=FALSE, lty, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(x, "rma"))
      stop(mstyle$stop("Argument 'x' must be an object of class \"rma\"."))

   if (inherits(x, "rma.ls"))
      stop(mstyle$stop("Method not available for objects of class \"rma.ls\"."))

   if (inherits(x, "rma.uni.selmodel"))
      stop(mstyle$stop("Method not available for objects of class \"rma.uni.selmodel\"."))

   if (!x$int.only)
      stop(mstyle$stop("L'Abbe plot only applicable for models without moderators."))

   if (!is.element(x$measure, c("RR","OR","RD","AS","IRR","IRD","IRSD")))
      stop(mstyle$stop("Argument 'measure' must be one of the following: 'RR','OR','RD','AS','IRR','IRD','IRSD'."))

   na.act <- getOption("na.action")
   on.exit(options(na.action=na.act))

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   if (length(add) == 2L) ### for rma.mh and rma.peto objects (1st 'add' value applies to the individual outcomes)
      add <- add[1]

   if (length(to) == 2L)  ### for rma.mh and rma.peto objects (1st 'to' value applies to the individual outcomes)
      to <- to[1]

   if (!is.element(to, c("all","only0","if0all","none")))
      stop(mstyle$stop("Unknown 'to' argument specified."))

   if (missing(transf))
      transf <- FALSE

   transf.char <- deparse(substitute(transf))

   if (missing(targs))
      targs <- NULL

   if (missing(psize))
      psize <- NULL

   if (missing(lty)) {
      lty <- c("solid", "dashed") ### 1st value = diagonal line, 2nd value = estimated effect line
   } else {
      if (length(lty) == 1L)
         lty <- c(lty, lty)
   }

   ### get ... argument

   ddd <- list(...)

   ### set defaults or get addyi and addvi arguments

   addyi <- ifelse(is.null(ddd$addyi), TRUE, ddd$addyi)
   addvi <- ifelse(is.null(ddd$addvi), TRUE, ddd$addvi)

   #########################################################################

   k <- x$k.f

   if (length(pch) == 1L)                       ### note: pch must have same length as number of tables (including NAs)
      pch <- rep(pch, k)                        ### or be equal to a single value (which is then repeated)

   if (length(pch) != k)
      stop(mstyle$stop(paste0("Number of tables (", k, ") does not correspond to the length of the 'pch' argument (", length(pch), ").")))

   ### if user has set the point sizes

   if (!is.null(psize)) {                       ### note: psize must have same length as number of tables (including NAs)
      if (length(psize) == 1L)                  ### or be equal to a single value (which is then repeated)
         psize <- rep(psize, k)
      if (length(psize) != k)
         stop(mstyle$stop(paste0("Number of tables (", k, ") does not correspond to the length of the 'psize' argument (", length(psize), ").")))
   }

   if (length(bg) == 1L)                        ### note: bg must have same length as number of tables (including NAs)
      bg <- rep(bg, k)                          ### or be equal to a single value (which is then repeated)

   if (length(bg) != k)
      stop(mstyle$stop(paste0("Number of tables (", k, ") does not correspond to the length of the 'bg' argument (", length(bg), ").")))

   #########################################################################

   ### these vectors may contain NAs

   ai  <- x$ai.f
   bi  <- x$bi.f
   ci  <- x$ci.f
   di  <- x$di.f
   x1i <- x$x1i.f
   x2i <- x$x2i.f
   t1i <- x$t1i.f
   t2i <- x$t2i.f

   ### drop00=TRUE may induce that the contrast-based yi value is NA; so
   ### make sure that the corresponding arm-based yi values are also NA

   yi.is.na <- is.na(x$yi.f)
   ai[yi.is.na]  <- NA
   bi[yi.is.na]  <- NA
   ci[yi.is.na]  <- NA
   di[yi.is.na]  <- NA
   x1i[yi.is.na] <- NA
   x2i[yi.is.na] <- NA
   t1i[yi.is.na] <- NA
   t2i[yi.is.na] <- NA

   options(na.action = "na.pass") ### to make sure dat.t and dat.c are of the same length

   if (x$measure == "RR") {
      measure <- "PLN"
      dat.t <- escalc(measure=measure, xi=ai, mi=bi, add=add, to=to, addyi=addyi, addvi=addvi)
      dat.c <- escalc(measure=measure, xi=ci, mi=di, add=add, to=to, addyi=addyi, addvi=addvi)
   }

   if (x$measure == "OR") {
      measure <- "PLO"
      dat.t <- escalc(measure=measure, xi=ai, mi=bi, add=add, to=to, addyi=addyi, addvi=addvi)
      dat.c <- escalc(measure=measure, xi=ci, mi=di, add=add, to=to, addyi=addyi, addvi=addvi)
   }

   if (x$measure == "RD") {
      measure <- "PR"
      dat.t <- escalc(measure=measure, xi=ai, mi=bi, add=add, to=to, addyi=addyi, addvi=addvi)
      dat.c <- escalc(measure=measure, xi=ci, mi=di, add=add, to=to, addyi=addyi, addvi=addvi)
   }

   if (x$measure == "AS") {
      measure <- "PAS"
      dat.t <- escalc(measure=measure, xi=ai, mi=bi, add=add, to=to, addyi=addyi, addvi=addvi)
      dat.c <- escalc(measure=measure, xi=ci, mi=di, add=add, to=to, addyi=addyi, addvi=addvi)
   }

   if (x$measure == "IRR") {
      measure <- "IRLN"
      dat.t <- escalc(measure=measure, xi=x1i, ti=t1i, add=add, to=to, addyi=addyi, addvi=addvi)
      dat.c <- escalc(measure=measure, xi=x2i, ti=t2i, add=add, to=to, addyi=addyi, addvi=addvi)
   }

   if (x$measure == "IRD") {
      measure <- "IR"
      dat.t <- escalc(measure=measure, xi=x1i, ti=t1i, add=add, to=to, addyi=addyi, addvi=addvi)
      dat.c <- escalc(measure=measure, xi=x2i, ti=t2i, add=add, to=to, addyi=addyi, addvi=addvi)
   }

   if (x$measure == "IRSD") {
      measure <-
      dat.t <- escalc(measure=measure, xi=x1i, ti=t1i, add=add, to=to, addyi=addyi, addvi=addvi)
      dat.c <- escalc(measure=measure, xi=x2i, ti=t2i, add=add, to=to, addyi=addyi, addvi=addvi)
   }

   options(na.action = na.act)

   ### check for NAs in yi/vi pairs and filter out

   has.na <- apply(is.na(dat.t), 1, any) | apply(is.na(dat.c), 1, any)
   not.na <- !has.na

   if (any(has.na)) {

      dat.t <- dat.t[not.na,]
      dat.c <- dat.c[not.na,]
      pch   <- pch[not.na]
      bg    <- bg[not.na]

   }

   if (length(dat.t$yi)==0L || length(dat.c$yi)==0L)
      stop(mstyle$stop("No information in object to compute arm-level outcomes."))

   #########################################################################

   ### determine point sizes

   if (is.null(psize)) {
      vi <- dat.t$vi + dat.c$vi
      wi <- 1/sqrt(vi)
      rng <- max(wi) - min(wi)
      if (rng <= .Machine$double.eps^0.5) {
         psize <- rep(1, length(wi))
      } else {
         psize <- 0.5 + 3 * (wi - min(wi))/rng
      }
   } else {
      psize <- psize[not.na]
   }

   ### determine x/y values for line that indicates the estimated effect

   min.yi <- min(c(dat.t$yi, dat.c$yi))
   max.yi <- max(c(dat.t$yi, dat.c$yi))
   rng.yi <- max.yi - min.yi

   len <- 1000

   intrcpt <- c(x$beta)

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
   psize.o     <- psize[order.vec]
   pch.o       <- pch[order.vec]
   bg.o        <- bg[order.vec]

   ### add x axis label

   if (missing(xlab)) {
      xlab <- .setlab(measure, transf.char, atransf.char="FALSE", gentype=1)
      xlab <- paste(xlab, "(Group 1)")
   }

   ### add y axis label

   if (missing(ylab)) {
      ylab <- .setlab(measure, transf.char, atransf.char="FALSE", gentype=1)
      ylab <- paste(ylab, "(Group 2)")
   }

   plot(NA, NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)

   ### add grid (and redraw box)

   if (grid) {
      grid()
      box(...)
   }

   abline(a=0, b=1, lty=lty[1], ...)
   lines(c.vals, t.vals, lty=lty[2], ...)
   points(dat.c$yi.o, dat.t$yi.o, cex=psize.o, pch=pch.o, bg=bg.o, ...)

   #########################################################################

   ### prepare data frame to return
   sav <- data.frame(x=dat.c$yi, y=dat.t$yi, cex=psize, pch=pch, bg=bg, ids=x$ids[not.na], slab=x$slab[not.na], stringsAsFactors=FALSE)

   invisible(sav)

}
