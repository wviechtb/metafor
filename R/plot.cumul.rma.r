plot.cumul.rma <- function(x, yaxis, xlim, ylim, xlab, ylab, at, transf, atransf, targs,
digits, cols=c("gray80","gray10"), grid=TRUE, pch=19, cex=1, lwd=2, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="cumul.rma")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   if (missing(yaxis)) {
      if (is.null(x$tau2)) {
         yaxis <- "I2"
      } else {
         yaxis <- "tau2"
      }
   } else {
      yaxis <- match.arg(yaxis, c("tau2","I2","H2"))
      if (is.null(x$tau2))
         stop(mstyle$stop("Cannot use yaxis=\"tau2\" for fixed-effects models."))
   }

   if (missing(transf))
      transf <- FALSE

   if (missing(atransf))
      atransf <- FALSE

   transf.char  <- deparse(substitute(transf))
   atransf.char <- deparse(substitute(atransf))

   if (is.function(transf) && is.function(atransf))
      stop(mstyle$stop("Use either 'transf' or 'atransf' to specify a transformation (not both)."))

   if (missing(xlab))
      xlab <- .setlab(x$measure, transf.char, atransf.char, gentype=2)

   if (missing(ylab)) {
      if (yaxis == "tau2")
         #ylab <- "Amount of Heterogeneity (tau^2)"
         ylab <- expression(paste("Amount of Heterogeneity ", (tau^2)))
      if (yaxis == "I2")
         #ylab <- "Percentage of Variability due to Heterogeneity (I^2)"
         ylab <- expression(paste("Percentage of Variability due to Heterogeneity ", (I^2)))
      if (yaxis == "H2")
         #ylab <- "Ratio of Total Variability to Sampling Variability (H^2)"
         ylab <- expression(paste("Ratio of Total Variability to Sampling Variability ", (H^2)))
   }

   par.mar <- par("mar")
   par.mar.adj <- par.mar + c(0,0.5,0,0) # need a bit more space on the right for the y-axis label
   par(mar = par.mar.adj)
   on.exit(par(mar = par.mar))

   if (missing(at))
      at <- NULL

   if (missing(targs))
      targs <- NULL

   if (missing(digits)) {
      if (yaxis == "tau2")
         digits <- c(2L,3L)
      if (yaxis == "I2")
         digits <- c(2L,1L)
      if (yaxis == "H2")
         digits <- c(2L,1L)
   } else {
      if (length(digits) == 1L)     ### digits[1] for x-axis labels
         digits <- c(digits,digits) ### digits[2] for y-axis labels
   }

   ### note: digits can also be a list (e.g., digits=list(2L,3)); trailing 0's are dropped for intergers

   ddd <- list(...)

   if (!is.null(ddd$addgrid))
      grid <- ddd$addgrid

   ### grid argument can either be a logical or a color

   if (is.logical(grid))
      gridcol <- "lightgray"
   if (is.character(grid)) {
      gridcol <- grid
      grid <- TRUE
   }

   lplot <- function(..., addgrid) plot(...)
   laxis <- function(..., addgrid) axis(...)

   #########################################################################

   ### set up data frame with the values to be plotted

   dat <- data.frame(estim=x$estimate)

   if (yaxis == "tau2")
      dat$yval <- x$tau2
   if (yaxis == "I2")
      dat$yval <- x$I2
   if (yaxis == "H2")
      dat$yval <- x$H2

   ### apply chosen na.action

   if (na.act == "na.fail" && anyNA(dat))
      stop(mstyle$stop("Missing values in results."))

   if (na.act == "na.omit")
      dat <- na.omit(dat)

   ### number of remaining rows/points

   k <- nrow(dat)

   ### if requested, apply transformation to estimates

   if (is.function(transf)) {
      if (is.null(targs)) {
         dat$estim <- sapply(dat$estim, transf)
      } else {
         dat$estim <- sapply(dat$estim, transf, targs)
      }
   }

   ### set xlim and ylim values

   if (missing(xlim)) {
      xlim <- range(dat$estim, na.rm=TRUE)
   } else {
      xlim <- sort(xlim) ### just in case the user supplies the limits in the wrong order
   }

   if (missing(ylim)) {
      ylim <- range(dat$yval, na.rm=TRUE)
   } else {
      ylim <- sort(ylim) ### just in case the user supplies the limits in the wrong order
   }

   ### if user has specified 'at' argument, make sure xlim actually contains the min and max 'at' values

   if (!is.null(at)) {
      xlim[1] <- min(c(xlim[1], at), na.rm=TRUE)
      xlim[2] <- max(c(xlim[2], at), na.rm=TRUE)
   }

   ### set up plot

   lplot(NA, NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, xaxt="n", yaxt="n", ...)

   ### generate x-axis positions if none are specified

   if (is.null(at)) {
      at <- axTicks(side=1)
   } else {
      at <- at[at > par("usr")[1]]
      at <- at[at < par("usr")[2]]
   }

   at.lab <- at

   if (is.function(atransf)) {
      if (is.null(targs)) {
         at.lab <- formatC(sapply(at.lab, atransf), digits=digits[[1]], format="f", drop0trailing=is.integer(digits[[1]]))
      } else {
         at.lab <- formatC(sapply(at.lab, atransf, targs), digits=digits[[1]], format="f", drop0trailing=is.integer(digits[[1]]))
      }
   } else {
      at.lab <- formatC(at.lab, digits=digits[[1]], format="f", drop0trailing=is.integer(digits[[1]]))
   }

   ### add x-axis

   laxis(side=1, at=at, labels=at.lab, ...)

   ### add y-axis

   aty <- axTicks(side=2)
   laxis(side=2, at=aty, labels=formatC(aty, digits=digits[[2]], format="f", drop0trailing=is.integer(digits[[2]])), ...)

   ### add grid

   if (.isTRUE(grid)) {
      abline(v=at,  lty="dotted", col=gridcol)
      abline(h=aty, lty="dotted", col=gridcol)
   }

   ### vector with color gradient for points

   cols.points <- colorRampPalette(cols)(k)

   #gray.vals.points <- seq(from=.9, to=.1, length.out=k)
   #cols.points <- gray(gray.vals.points)
   #cols <- colorRampPalette(c("yellow","red"))(k)
   #cols <- colorRampPalette(c("blue","red"))(k)
   #cols <- rev(heat.colors(k+4))[-c(1:2,(k+1):(k+2)]
   #cols <- rev(topo.colors(k))
   #cols <- rainbow(k, start=.2, end=.4)

   ### add lines that have a gradient (by interpolating values)
   ### looks better this way, especially when k is low

   for (i in seq_len(k-1)) {
      if (is.na(dat$estim[i]) || is.na(dat$estim[i+1]) || is.na(dat$yval[i]) || is.na(dat$yval[i+1]))
         next
      estims <- approx(c(dat$estim[i], dat$estim[i+1]), n=50)$y
      yvals  <- approx(c(dat$yval[i], dat$yval[i+1]), n=50)$y
      cols.lines <- colorRampPalette(c(cols.points[i], cols.points[i+1]))(50)
      #gray.vals.lines <- approx(c(gray.vals.points[i], gray.vals.points[i+1]), n=50)$y
      #cols.lines <- gray(gray.vals.lines)
      segments(estims[-50], yvals[-50], estims[-1], yvals[-1], col=cols.lines, lwd=lwd)
   }

   ### add lines (this does no interpolation)
   #segments(dat$estim[-k], dat$yval[-k], dat$estim[-1], dat$yval[-1], col=cols.points, lwd=lwd)

   ### add points

   points(x=dat$estim, y=dat$yval, pch=pch, col=cols.points, cex=cex)

   ### redraw box around plot

   box()

   ### return data frame invisibly

   dat$col <- cols.points

   invisible(dat)

}
