plot.gosh.rma <- function(x, het="I2", pch=19, cex, out, col, alpha, xhist=TRUE, yhist=TRUE, hh=.3, hcol, breaks, adjust, lwd, labels, ...) {

   if (!inherits(x, "gosh.rma"))
      stop("Argument 'x' must be an object of class \"gosh.rma\".")

   het <- match.arg(het, c("QE", "I2", "H2", "tau2"))

   if (het == "tau2" && x$method == "FE")
      stop("Cannot plot 'tau2' for fixed-effects models.")

   if (missing(cex))
      cex <- 0.5

   if (missing(alpha))
      alpha <- nrow(x$res)^(-.2)

   if (missing(out)) {

      if (missing(col))
         col <- rgb(0,0,0,alpha)

   } else {

      if (length(out) != 1)
         stop("Argument 'out' should only specify a single study.")

      if (out > x$k || out < 1)
         stop("Non-existing study chosen as potential outlier.")

      if (missing(col)) {
         col <- ifelse(x$incl[,out], rgb(1,0,0,alpha), rgb(0,0,1,alpha))
      } else {
         if (length(col) == 2)
            col <- ifelse(x$incl[,out], col[1], col[2])
      }
   }

   if (length(hh) == 1)
      hh <- c(hh, hh)

   if (any(hh < 0) | any(hh > 1))
      stop("Invalid value(s) specified for 'hh' argument.")

   if (missing(hcol))
      hcol <- c("gray70", "gray90")

   if (length(hcol) == 1)
      hcol <- c(hcol, hcol)

   if (missing(breaks))
      breaks <- "Sturges"

   if (missing(adjust))
      adjust <- 1

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
         labels <- c(labels, colnames(x$res)[-c(1:5)])
      }

   }

   #########################################################################

   if (x$int.only) {

      par.mar <- par("mar")
      par.mar.adj <- par.mar - c(0,-1,3.5,1.5)
      par.mar.adj[par.mar.adj < 0] <- 0
      on.exit(par(mar = par.mar))

      if (xhist & yhist)
         layout(mat=matrix(c(1,2,3,4), nrow=2, byrow=TRUE), widths=c(1-hh[2],hh[2]), heights=c(hh[1],1-hh[1]))
      if (xhist & !yhist)
         layout(mat=matrix(c(1,2), nrow=2, byrow=TRUE), heights=c(hh[1],1-hh[1]))
      if (!xhist & yhist)
         layout(mat=matrix(c(1,2), nrow=1, byrow=TRUE), widths=c(1-hh[2],hh[2]))

      if (xhist) {
         d <- density(x$res[,6], adjust=adjust, na.rm=TRUE)
         h <- hist(x$res[,6], breaks=breaks, plot=FALSE)
         brks <- h$breaks
         nB <- length(brks)
         y <- h$density
         par(mar=c(0,par.mar.adj[2:4]))
         plot(NULL, xlim = c(range(h$breaks)), ylim = c(0,max(h$density,d$y)), xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
         rect(brks[-nB], 0, brks[-1], y, col=hcol[1], border=hcol[2])
         if (lwd > 0)
            lines(d$x, d$y, lwd=lwd, col="black")
      }

      if (xhist & yhist)
         plot.new()

      par(mar = par.mar.adj)
      plot(x$res[,6], x$res[,het], pch=pch, cex=cex, col=col, bty="l", xlab=labels[1], ylab=labels[2], ...)

      if (yhist) {
         d <- density(x$res[,het], adjust=adjust, na.rm=TRUE)
         h <- hist(x$res[,het], breaks=breaks, plot=FALSE)
         brks <- h$breaks
         nB <- length(brks)
         y <- h$density
         par(mar=c(par.mar.adj[1],0,par.mar.adj[3:4]))
         plot(NULL, xlim = c(0,max(h$density,d$y)), ylim = c(range(h$breaks)), xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
         rect(0, brks[-nB], y, brks[-1], col=hcol[1], border=hcol[2])
         if (lwd > 0)
            lines(d$y, d$x, lwd=lwd, col="black")
      }

   } else {

      ### function for histograms with kernel density estimates on the diagonal

      panel.hist <- function(x, ...) {
         usr <- par("usr")
         on.exit(par(usr))
         par(usr = c(usr[1:2], 0, 1.5))
         h <- hist(x, plot=FALSE, breaks=breaks)
         breaks <- h$breaks
         nB <- length(breaks)
         y <- h$density
         z <- y / max(y)
         rect(breaks[-nB], 0, breaks[-1], z, col=hcol[1], border=hcol[2])
         res <- density(x, adjust=adjust, na.rm=TRUE)
         res$y <- res$y / max(y)
         if (lwd > 0)
            lines(res, lwd=lwd, col="black")
      }

      ### draw scatterplot matrix

      X <- cbind(x$res[,het], x$res[,6:ncol(x$res)])
      pairs(X, pch=pch, cex=cex, diag.panel=panel.hist, col=col, labels=labels, ...)

   }

   #########################################################################

}
