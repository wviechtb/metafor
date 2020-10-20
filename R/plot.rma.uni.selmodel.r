plot.rma.uni.selmodel <- function(x, xlim, n=1001, prec="max", scale=FALSE, rug=TRUE, add=FALSE, lwd=2, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(x, "rma.uni.selmodel"))
      stop(mstyle$stop("Argument 'x' must be an object of class \"rma.uni.selmodel\"."))

   if (missing(xlim))
      xlim <- c(x$pval.min, 1-x$pval.min)

   if (length(xlim) != 2L)
      stop(mstyle$stop("Argument 'xlim' should be a vector of length 2."))

   xlim <- sort(xlim)

   if (xlim[1] < 0 || xlim[2] > 1)
      stop(mstyle$stop("Values for 'xlim' should be between 0 and 1."))

   if (length(prec) != 1L)
      stop(mstyle$stop("Argument 'prec' should be of length 1."))

   if (is.character(prec)) {
      if (!is.element(prec, c("min", "max", "mean", "median")))
         stop(mstyle$stop("Unknown options specified for the 'prec' argument."))
      if (prec == "min")
         prec <- x$precis[["min"]]
      if (prec == "max")
         prec <- x$precis[["max"]]
      if (prec == "mean")
         prec <- x$precis[["mean"]]
      if (prec == "median")
         prec <- x$precis[["median"]]
   }

   delta <- x$delta
   steps <- x$steps
   xs <- seq(xlim[1], xlim[2], length=n)

   if (x$type == "stepfun") {
      xs <- unique(sort(c(xs, steps)))
      xs <- xs[xs >= xlim[1]]
      xs <- xs[xs <= xlim[2]]
      plot.type <- "S"
   } else {
      plot.type <- "l"
   }

   if (x$type == "beta")
      ys <- dbeta(xs, delta[1], delta[2])

   if (x$type == "halfnorm")
      ys <- ifelse(xs <= steps[1], 1, exp(-delta * prec * xs^2) / exp(-delta * prec * steps[1]^2))

   if (x$type == "negexp")
      ys <- ifelse(xs <= steps[1], 1, exp(-delta * prec * xs) / exp(-delta * prec * steps[1]))

   if (x$type == "logistic")
      ys <- ifelse(xs <= steps[1], 1, (2 * exp(-delta * prec * xs) / (1 + exp(-delta * prec * xs))) / (2 * exp(-delta * prec * steps[1]) / (1 + exp(-delta * prec * steps[1]))))

   if (x$type == "power")
      ys <- ifelse(xs <= steps[1], 1, (1-xs)^(prec*delta) / (1-steps[1])^(prec*delta))

   if (x$type == "negexppow")
      ys <- ifelse(xs <= steps[1], 1, exp(-delta[1] * prec * xs^(1/delta[2])) / exp(-delta[1] * prec * steps[1]^(1/delta[2])))

   if (x$type == "halfnorm2")
      ys <- ifelse(xs <= steps[1], 1, (delta[1] + exp(-delta[2] * prec * xs^2) / exp(-delta[2] * prec * steps[1]^2)) / (1 + delta[1]))

   if (x$type == "negexp2")
      ys <- ifelse(xs <= steps[1], 1, (delta[1] + exp(-delta[2] * prec * xs) / exp(-delta[2] * prec * steps[1])) / (1 + delta[1]))

   if (x$type == "logistic2")
      ys <- ifelse(xs <= steps[1], 1, (delta[1] + (2 * exp(-delta[2] * prec * xs) / (1 + exp(-delta[2] * prec * xs))) / (2 * exp(-delta[2] * prec * steps[1]) / (1 + exp(-delta[2] * prec * steps[1])))) / (1 + delta[1]))

   if (x$type == "power2")
      ys <- ifelse(xs <= steps[1], 1, (delta[1] + (1-xs)^(prec*delta[2]) / (1-steps[1])^(prec*delta[2])) / (1 + delta[1]))

   if (x$type == "stepfun")
      ys <- delta[sapply(xs, function(p) which(p <= x$steps)[1])] / prec


   if (scale) {

      is.inf.pos <- ys ==  Inf
      is.inf.neg <- ys == -Inf

      ys[is.infinite(ys)] <- NA

      rng <- max(ys, na.rm=TRUE) - min(ys, na.rm=TRUE)

      if (rng > .Machine$double.eps^0.5)
         ys <- (ys - min(ys, na.rm=TRUE)) / rng

      ys[is.inf.pos] <- 1
      ys[is.inf.neg] <- 0

   }

   if (add) {
      lines(xs, ys, type=plot.type, lwd=lwd, ...)
   } else {
      plot(xs, ys, type=plot.type, lwd=lwd, xlab="p-value", ylab="Relative Likelihood of Selection", ...)
   }

   if (rug && !add)
      rug(x$pvals, quiet=TRUE)

   sav <- data.frame(xs=xs, ys=ys)

   invisible(sav)

}
