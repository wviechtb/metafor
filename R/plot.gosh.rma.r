plot.gosh.rma <- function(x, het="I2", pch=19, cex, out, col, alpha, breaks, adjust, lwd, labels, ...) {

   if (!is.element("gosh.rma", class(x)))
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

      labels <- c(labels, colnames(x$res)[-c(1:5)])

   }

   #########################################################################

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
      rect(breaks[-nB], 0, breaks[-1], z, col="gray60", border="gray90")
      res <- density(x, adjust=adjust, na.rm=TRUE)
      res$y <- res$y / max(y)
      lines(res, lwd=lwd, col="black")
   }

   ### draw scatterplot matrix

   X <- cbind(x$res[,het], x$res[,6:ncol(x$res)])
   pairs(X, pch=pch, cex=cex, diag.panel=panel.hist, col=col, labels=labels, ...)

}
