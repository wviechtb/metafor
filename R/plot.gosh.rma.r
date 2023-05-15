plot.gosh.rma <- function(x, het="I2", pch=16, cex, out, col, alpha, border,
xlim, ylim, xhist=TRUE, yhist=TRUE, hh=0.3, breaks, adjust, lwd, labels, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="gosh.rma")

   het <- match.arg(het, c("QE", "I2", "H2", "tau2"))

   if (het == "tau2" && is.element(x$method, c("FE","EE","CE")))
      stop(mstyle$stop("Cannot plot 'tau2' for equal/fixed-effects models."))

   if (missing(cex)) {
      cex <- par("cex") * 0.5
   } else {
      cex <- par("cex") * cex
   }

   ddd <- list(...)

   if (!is.null(ddd$trim)) {

      trim <- ddd$trim

      if (!is.list(trim)) {
         if (length(trim) == 1L)
            trim <- rep(trim, ncol(x$res)-4L)
         trim <- as.list(trim)
      }

      X <- cbind(x$res[,het], x$res[,6:ncol(x$res)])
      del <- rep(FALSE, nrow(X))
      for (i in seq_len(ncol(X))) {
         del[X[,i] < quantile(X[,i], trim[[i]][1], na.rm=TRUE) | X[,i] > quantile(X[,i], 1-trim[[i]][length(trim[[i]])], na.rm=TRUE)] <- TRUE
      }

      del[is.na(del)] <- TRUE
      x$res <- x$res[!del,]
      x$incl <- x$incl[!del,]

   }

   if (exists(".darkplots") && .isTRUE(.darkplots))
      par(fg="gray95", bg="gray10", col="gray95", col.axis="gray95", col.lab="gray95", col.main="gray95", col.sub="gray95")

   lplot  <- function(..., trim) plot(...)
   lpairs <- function(..., trim) pairs(...)

   if (missing(alpha))
      alpha <- nrow(x$res)^(-0.2)

   if (length(alpha) == 1L)
      alpha <- c(alpha, 0.5, 0.9) ### 1st for points, 2nd for histograms, 3rd for density lines

   if (length(alpha) == 2L)
      alpha <- c(alpha[1], alpha[2], 0.9)

   missout <- ifelse(missing(out), TRUE, FALSE) ### need this for panel.hist()

   if (missout) {

      if (missing(col))
         col <- par("fg")

      col <- col2rgb(col) / 255
      col.pnts <- rgb(col[1], col[2], col[3], alpha[1])
      col.hist <- rgb(col[1], col[2], col[3], alpha[2])
      col.line <- rgb(col[1], col[2], col[3], alpha[3])

   } else {

      if (length(out) != 1L)
         stop(mstyle$stop("Argument 'out' should only specify a single study."))

      out <- round(out)

      if (out > x$k || out < 1)
         stop(mstyle$stop("Non-existing study chosen as potential outlier."))

      if (missing(col)) {
         if (.is.dark(par("bg"))) {
            col <- c("firebrick", "dodgerblue")
         } else {
            col <- c("red", "blue")
         }
      }

      if (length(col) != 2L)
         stop(mstyle$stop("Argument 'col' should specify two colors when argument 'out' is used."))

      col.o <- col2rgb(col[1]) / 255
      col.i <- col2rgb(col[2]) / 255
      col.pnts.o <- rgb(col.o[1], col.o[2], col.o[3], alpha[1])
      col.pnts.i <- rgb(col.i[1], col.i[2], col.i[3], alpha[1])
      col.pnts   <- ifelse(x$incl[,out], col.pnts.o, col.pnts.i)
      col.hist.o <- rgb(col.o[1], col.o[2], col.o[3], alpha[2])
      col.hist.i <- rgb(col.i[1], col.i[2], col.i[3], alpha[2])
      col.line.o <- rgb(col.o[1], col.o[2], col.o[3], alpha[3])
      col.line.i <- rgb(col.i[1], col.i[2], col.i[3], alpha[3])

   }

   if (missing(border)) {
      if (.is.dark(par("bg"))) {
         border <- par("bg")
      } else {
         border <- "white"
      }
   }

   if (length(border) == 1L)
      border <- c(border, border)

   if (length(hh) == 1L)
      hh <- c(hh, hh)

   if (x$int.only && (any(hh < 0) | any(hh > 1)))
      stop(mstyle$stop("Invalid value(s) specified for 'hh' argument."))

   if (missing(breaks))
      breaks <- "Sturges"

   if (length(breaks) == 1L)
      breaks <- list(breaks, breaks) # use list so can also specify two vectors (or two functions)

   if (missing(adjust))
      adjust <- 1

   if (length(adjust) == 1L)
      adjust <- c(adjust, adjust)

   if (missing(lwd))
      lwd <- 2

   if (length(lwd) == 1L)
      lwd <- c(lwd, lwd)

   if (missing(labels)) {

      if (het == "QE" && x$int.only)
         labels <- expression(Q)
      if (het == "QE" && !x$int.only)
         labels <- expression(Q[E])
      if (het == "I2")
         labels <- expression(I^2)
      if (het == "H2")
         labels <- expression(H^2)
      if (het == "tau2")
         labels <- expression(tau^2)

      if (x$int.only) {
         labels <- c(.setlab(x$measure, transf.char="FALSE", atransf.char="FALSE", gentype=2), labels)
      } else {
         labels <- c(labels, colnames(x$res)[-seq_len(5)])
      }

   }

   #########################################################################

   if (x$int.only) {

      par.mar <- par("mar")
      par.mar.adj <- par.mar - c(0,-1,3.1,1.1)
      par.mar.adj[par.mar.adj < 0] <- 0
      on.exit(par(mar = par.mar), add=TRUE)

      if (xhist & yhist)
         layout(mat=matrix(c(1,2,3,4), nrow=2, byrow=TRUE), widths=c(1-hh[2],hh[2]), heights=c(hh[1],1-hh[1]))
      if (xhist & !yhist)
         layout(mat=matrix(c(1,2), nrow=2, byrow=TRUE), heights=c(hh[1],1-hh[1]))
      if (!xhist & yhist)
         layout(mat=matrix(c(1,2), nrow=1, byrow=TRUE), widths=c(1-hh[2],hh[2]))

      hx <- hist(x$res[,6],   breaks=breaks[[1]], plot=FALSE)
      hy <- hist(x$res[,het], breaks=breaks[[2]], plot=FALSE)

      if (missout) {

         if (missing(xlim))
            xlim <- range(hx$breaks)
         if (missing(ylim))
            ylim <- range(hy$breaks)

         if (xhist) {
            d <- density(x$res[,6], adjust=adjust[1], na.rm=TRUE)
            brks <- hx$breaks
            nB <- length(brks)
            y <- hx$density
            par(mar=c(0,par.mar.adj[2:4]))
            plot(NULL, xlim=xlim, ylim=c(0,max(hx$density,d$y)), xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
            rect(brks[-nB], 0, brks[-1], y, col=col.hist, border=border[1])
            if (lwd[1] > 0)
               lines(d$x, d$y, lwd=lwd[1], col=col.line)
         }

      } else {

         isout <- x$incl[,out]

         hx.o <- hist(x$res[isout,6],  breaks=hx$breaks, plot=FALSE)
         hx.i <- hist(x$res[!isout,6], breaks=hx$breaks, plot=FALSE)

         hy.o <- hist(x$res[isout,het],  breaks=hy$breaks, plot=FALSE)
         hy.i <- hist(x$res[!isout,het], breaks=hy$breaks, plot=FALSE)

         if (missing(xlim))
            xlim <- c(min(hx.o$breaks, hx.i$breaks), max(hx.o$breaks, hx.i$breaks))
         if (missing(ylim))
            ylim <- c(min(hy.o$breaks, hy.i$breaks), max(hy.o$breaks, hy.i$breaks))

         if (xhist) {
            d.o <- density(x$res[isout,6],  adjust=adjust[1], na.rm=TRUE)
            d.i <- density(x$res[!isout,6], adjust=adjust[1], na.rm=TRUE)
            brks.o <- hx.o$breaks
            brks.i <- hx.i$breaks
            nB.o <- length(brks.o)
            nB.i <- length(brks.i)
            y.o <- hx.o$density
            y.i <- hx.i$density
            par(mar=c(0,par.mar.adj[2:4]))
            plot(NULL, xlim=xlim, ylim=c(0,max(hx.o$density,hx.i$density,d.o$y,d.i$y)), xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
            rect(brks.i[-nB.i], 0, brks.i[-1], y.i, col=col.hist.i, border=border[1])
            rect(brks.o[-nB.o], 0, brks.o[-1], y.o, col=col.hist.o, border=border[1])
            if (lwd[1] > 0) {
               lines(d.i$x, d.i$y, lwd=lwd[1], col=col.line.i)
               lines(d.o$x, d.o$y, lwd=lwd[1], col=col.line.o)
            }
         }

      }

      if (xhist & yhist)
         plot.new()

      par(mar = par.mar.adj)
      lplot(x$res[,6], x$res[,het], xlim=xlim, ylim=ylim, pch=pch, cex=cex, col=col.pnts, bty="l", xlab=labels[1], ylab=labels[2], ...)

      if (missout) {

         if (yhist) {
            d <- density(x$res[,het], adjust=adjust[2], na.rm=TRUE)
            brks <- hy$breaks
            nB <- length(brks)
            y <- hy$density
            par(mar=c(par.mar.adj[1],0,par.mar.adj[3:4]))
            plot(NULL, xlim=c(0,max(hy$density,d$y)), ylim=ylim, xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
            rect(0, brks[-nB], y, brks[-1], col=col.hist, border=border[2])
            if (lwd[2] > 0)
               lines(d$y, d$x, lwd=lwd[2], col=col.line)
         }

      } else {

         if (yhist) {
            d.o <- density(x$res[isout,het],  adjust=adjust[2], na.rm=TRUE)
            d.i <- density(x$res[!isout,het], adjust=adjust[2], na.rm=TRUE)
            brks.o <- hy.o$breaks
            brks.i <- hy.i$breaks
            nB.o <- length(brks.o)
            nB.i <- length(brks.i)
            y.o <- hy.o$density
            y.i <- hy.i$density
            par(mar=c(par.mar.adj[1],0,par.mar.adj[3:4]))
            plot(NULL, xlim=c(0,max(hy.o$density,hy.i$density,d.o$y,d.i$y)), ylim=ylim, xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
            rect(0, brks.i[-nB.i], y.i, brks.i[-1], col=col.hist.i, border=border[2])
            rect(0, brks.o[-nB.o], y.o, brks.o[-1], col=col.hist.o, border=border[2])
            if (lwd[2] > 0) {
               lines(d.i$y, d.i$x, lwd=lwd[2], col=col.line.i)
               lines(d.o$y, d.o$x, lwd=lwd[2], col=col.line.o)
            }
         }

      }

      ### reset to a single figure
      if (xhist | yhist)
         layout(matrix(1))

   } else {

      isout <- x$incl[,out]

      ### function for histograms with kernel density estimates on the diagonal

      panel.hist <- function(x, ...) {
         usr <- par("usr")
         on.exit(par(usr=usr))
         par(usr = c(usr[1:2], 0, 1.2 + hh[1]))
         h <- hist(x, plot=FALSE, breaks=breaks[[1]])
         if (missout) {
            brks <- h$breaks
            nB <- length(brks)
            y <- h$density
            z <- y / max(y)
            rect(brks[-nB], 0, brks[-1], z, col=col.hist, border=border[1])
            res <- density(x, adjust=adjust[1], na.rm=TRUE)
            res$y <- res$y / max(y)
            if (lwd[1] > 0)
               lines(res, lwd=lwd[1], col=col.line)
         } else {
            h.o <- hist(x[isout],  plot=FALSE, breaks=h$breaks)
            h.i <- hist(x[!isout], plot=FALSE, breaks=h$breaks)
            brks.o <- h.o$breaks
            brks.i <- h.i$breaks
            nB.o <- length(brks.o)
            nB.i <- length(brks.i)
            y.o <- h.o$density
            y.i <- h.i$density
            z.o <- y.o / max(y.o, y.i)
            z.i <- y.i / max(y.o, y.i)
            rect(brks.i[-nB.i], 0, brks.i[-1], z.i, col=col.hist.i, border=border[1])
            rect(brks.o[-nB.o], 0, brks.o[-1], z.o, col=col.hist.o, border=border[1])
            res.o <- density(x[isout],  adjust=adjust[1], na.rm=TRUE)
            res.i <- density(x[!isout], adjust=adjust[1], na.rm=TRUE)
            res.o$y <- res.o$y / max(y.o, y.i)
            res.i$y <- res.i$y / max(y.o, y.i)
            if (lwd[1] > 0) {
               lines(res.i, lwd=lwd[1], col=col.line.i)
               lines(res.o, lwd=lwd[1], col=col.line.o)
            }
         }
         box()
      }

      ### draw scatterplot matrix

      X <- cbind(x$res[,het], x$res[,6:ncol(x$res)])
      lpairs(X, pch=pch, cex=cex, diag.panel=panel.hist, col=col.pnts, labels=labels, ...)

   }

   #########################################################################

}
