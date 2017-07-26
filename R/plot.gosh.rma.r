plot.gosh.rma <- function(x, het="I2", pch=16, cex, out, col, alpha, border,
xlim, ylim, xhist=TRUE, yhist=TRUE, hh=.3, breaks, adjust, lwd, labels, ...) {

   if (!inherits(x, "gosh.rma"))
      stop("Argument 'x' must be an object of class \"gosh.rma\".")

   het <- match.arg(het, c("QE", "I2", "H2", "tau2"))

   if (het == "tau2" && x$method == "FE")
      stop("Cannot plot 'tau2' for fixed-effects models.")

   if (missing(cex))
      cex <- 0.5

   if (missing(alpha))
      alpha <- nrow(x$res)^(-0.2)

   if (length(alpha) == 1)
      alpha <- c(alpha, 0.5, 0.9)

   if (length(alpha) == 2)
      alpha <- c(alpha[1], alpha[2], 0.9)

   if (missing(border))
      border <- "white"

   if (missing(breaks)) {
      nbrks <- round(sum(x$fit)^(1/2.5))
      nbrks[nbrks < 10] <- 10
      nbrks[nbrks > 100] <- 100
   } else {
      nbrks <- breaks
   }

   if (length(nbrks) == 1)
      nbrks <- c(nbrks, nbrks)

   missout <- ifelse(missing(out), TRUE, FALSE) ### need this for panel.hist()

   if (missout) {

      if (missing(col))
         col <- "black"

      col <- col2rgb(col) / 255
      col.points <- rgb(col[1], col[2], col[3], alpha[1])
      col.hist <- rgb(col[1], col[2], col[3], alpha[2])
      col.line <- rgb(col[1], col[2], col[3], alpha[3])

   } else {

      if (length(out) != 1)
         stop("Argument 'out' should only specify a single study.")

      if (out > x$k || out < 1)
         stop("Non-existing study chosen as potential outlier.")

      if (missing(col))
         col <- c("red", "blue")

      if (length(col) != 2)
         stop("Argument 'col' should specify two colors.")

      col.o <- col2rgb(col[1]) / 255
      col.i <- col2rgb(col[2]) / 255
      col.points.o <- rgb(col.o[1], col.o[2], col.o[3], alpha[1])
      col.points.i <- rgb(col.i[1], col.i[2], col.i[3], alpha[1])
      col.points <- ifelse(x$incl[,out], col.points.o, col.points.i)
      col.hist.o <- rgb(col.o[1], col.o[2], col.o[3], alpha[2])
      col.hist.i <- rgb(col.i[1], col.i[2], col.i[3], alpha[2])
      col.line.o <- rgb(col.o[1], col.o[2], col.o[3], alpha[3])
      col.line.i <- rgb(col.i[1], col.i[2], col.i[3], alpha[3])

   }

   if (length(hh) == 1)
      hh <- c(hh, hh)

   if (x$int.only && (any(hh < 0) | any(hh > 1)))
      stop("Invalid value(s) specified for 'hh' argument.")

   if (missing(adjust))
      adjust <- 1

   if (length(adjust) == 1)
      adjust <- c(adjust, adjust)

   if (missing(lwd))
      lwd <- 2

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
      par.mar.adj <- par.mar - c(0,-1,2.5,0.5)
      par.mar.adj[par.mar.adj < 0] <- 0
      on.exit(par(mar = par.mar))

      if (xhist & yhist)
         layout(mat=matrix(c(1,2,3,4), nrow=2, byrow=TRUE), widths=c(1-hh[2],hh[2]), heights=c(hh[1],1-hh[1]))
      if (xhist & !yhist)
         layout(mat=matrix(c(1,2), nrow=2, byrow=TRUE), heights=c(hh[1],1-hh[1]))
      if (!xhist & yhist)
         layout(mat=matrix(c(1,2), nrow=1, byrow=TRUE), widths=c(1-hh[2],hh[2]))

      if (missout) {

         breaks <- seq(min(x$res[,6], na.rm=TRUE), max(x$res[,6], na.rm=TRUE), length=nbrks[1])

         hx <- hist(x$res[,6], breaks=breaks, plot=FALSE)

         breaks <- seq(min(x$res[,het], na.rm=TRUE), max(x$res[,het], na.rm=TRUE), length=nbrks[2])

         hy <- hist(x$res[,het], breaks=breaks, plot=FALSE)

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
            rect(brks[-nB], 0, brks[-1], y, col=col.hist, border=border)
            if (lwd > 0)
               lines(d$x, d$y, lwd=lwd, col=col.line)
         }

      } else {

         isout <- x$incl[,out]

         breaks <- seq(min(x$res[,6], na.rm=TRUE), max(x$res[,6], na.rm=TRUE), length=nbrks[1])

         hx.o <- hist(x$res[isout,6],  breaks=breaks, plot=FALSE)
         hx.i <- hist(x$res[!isout,6], breaks=breaks, plot=FALSE)

         breaks <- seq(min(x$res[,het], na.rm=TRUE), max(x$res[,het], na.rm=TRUE), length=nbrks[2])

         hy.o <- hist(x$res[isout,het],  breaks=breaks, plot=FALSE)
         hy.i <- hist(x$res[!isout,het], breaks=breaks, plot=FALSE)

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
            rect(brks.i[-nB.i], 0, brks.i[-1], y.i, col=col.hist.i, border=border)
            rect(brks.o[-nB.o], 0, brks.o[-1], y.o, col=col.hist.o, border=border)
            if (lwd > 0) {
               lines(d.i$x, d.i$y, lwd=lwd, col=col.line.i)
               lines(d.o$x, d.o$y, lwd=lwd, col=col.line.o)
            }
         }

      }

      if (xhist & yhist)
         plot.new()

      par(mar = par.mar.adj)
      plot(x$res[,6], x$res[,het], xlim=xlim, ylim=ylim, pch=pch, cex=cex, col=col.points, bty="l", xlab=labels[1], ylab=labels[2], ...)

      if (missout) {

         if (yhist) {
            d <- density(x$res[,het], adjust=adjust[2], na.rm=TRUE)
            brks <- hy$breaks
            nB <- length(brks)
            y <- hy$density
            par(mar=c(par.mar.adj[1],0,par.mar.adj[3:4]))
            plot(NULL, xlim=c(0,max(hy$density,d$y)), ylim=ylim, xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
            rect(0, brks[-nB], y, brks[-1], col=col.hist, border=border)
            if (lwd > 0)
               lines(d$y, d$x, lwd=lwd, col=col.line)
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
            rect(0, brks.i[-nB.i], y.i, brks.i[-1], col=col.hist.i, border=border)
            rect(0, brks.o[-nB.o], y.o, brks.o[-1], col=col.hist.o, border=border)
            if (lwd > 0) {
               lines(d.i$y, d.i$x, lwd=lwd, col=col.line.i)
               lines(d.o$y, d.o$x, lwd=lwd, col=col.line.o)
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
         on.exit(par(usr))
         par(usr = c(usr[1:2], 0, 1.2 + hh[1]))
         breaks <- seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE), length=nbrks[1])
         if (missout) {
            h <- hist(x, plot=FALSE, breaks=breaks)
            brks <- h$breaks
            nB <- length(brks)
            y <- h$density
            z <- y / max(y)
            rect(brks[-nB], 0, brks[-1], z, col=col.hist, border=border)
            res <- density(x, adjust=adjust[1], na.rm=TRUE)
            res$y <- res$y / max(y)
            if (lwd > 0)
               lines(res, lwd=lwd, col=col.line)
         } else {
            h.o <- hist(x[isout],  plot=FALSE, breaks=breaks)
            h.i <- hist(x[!isout], plot=FALSE, breaks=breaks)
            brks.o <- h.o$breaks
            brks.i <- h.i$breaks
            nB.o <- length(brks.o)
            nB.i <- length(brks.i)
            y.o <- h.o$density
            y.i <- h.i$density
            z.o <- y.o / max(y.o)
            z.i <- y.i / max(y.i)
            rect(brks.i[-nB.i], 0, brks.i[-1], z.i, col=col.hist.i, border=border)
            rect(brks.o[-nB.o], 0, brks.o[-1], z.o, col=col.hist.o, border=border)
            res.o <- density(x[isout],  adjust=adjust[1], na.rm=TRUE)
            res.i <- density(x[!isout], adjust=adjust[1], na.rm=TRUE)
            res.o$y <- res.o$y / max(y.o)
            res.i$y <- res.i$y / max(y.i)
            if (lwd > 0) {
               lines(res.i, lwd=lwd, col=col.line.i)
               lines(res.o, lwd=lwd, col=col.line.o)
            }
         }
         box()
      }

      ### draw scatterplot matrix

      X <- cbind(x$res[,het], x$res[,6:ncol(x$res)])
      pairs(X, pch=pch, cex=cex, diag.panel=panel.hist, col=col.points, labels=labels, ...)

   }

   #########################################################################

}
