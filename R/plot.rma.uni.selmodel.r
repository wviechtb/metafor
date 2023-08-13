plot.rma.uni.selmodel <- function(x, xlim, ylim, n=1000, prec="max", scale=FALSE,
   ci=FALSE, reps=1000, shade=TRUE, rug=TRUE, add=FALSE,
   lty=c("solid","dotted"), lwd=c(2,1), ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="rma.uni.selmodel")

   .start.plot(!add)

   if (is.element(x$type, c("trunc","truncest")))
      stop(mstyle$stop("Cannot draw the selection function for this type of selection model."))

   ### shade argument can either be a logical or a color

   if (is.logical(shade)) {
      shadecol <- .coladj(par("bg","fg"), dark=0.10, light=-0.10)
   }

   if (is.character(shade)) {
      shadecol <- shade
      shade <- TRUE
   }

   ddd <- list(...)

   lplot    <- function(..., seed) plot(...)
   llines   <- function(..., seed) lines(...)
   lrug     <- function(..., seed) rug(...)
   lpolygon <- function(..., seed) polygon(...)

   if (is.logical(ci))
      citype <- "boot"

   if (is.character(ci)) {
      citype <- tolower(ci)
      ci <- TRUE
   }

   if (!is.element(citype, c("boot", "wald")))
      stop(mstyle$stop("Unknown confidence interval type specified."))

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
   } else {
      if (is.numeric(prec) && !x$precspec)
         prec <- 1
   }

   delta <- x$delta
   steps <- x$steps

   ps <- seq(xlim[1], xlim[2], length.out=n)

   if (x$type == "stepfun") {
      ps <- unique(sort(c(ps, steps))) # make sure that the 'steps' values are part of 'ps'
      ps <- ps[ps >= xlim[1]]          # but only keep ps >= xlim[1]
      ps <- ps[ps <= xlim[2]]          #               ps <= xlim[2]
      plot.type <- "S"
   } else {
      plot.type <- "l"
   }

   wi.fun <- x$wi.fun

   ys <- wi.fun(ps, delta=delta, yi=x$yi, vi=x$vi, preci=prec, alternative=x$alternative, steps=x$steps)

   if (ci && citype == "boot" && all(is.na(x$vd)))
      ci <- FALSE

   if (ci && citype == "wald" && all(is.na(x$ci.lb.delta)) && all(is.na(x$ci.ub.delta)))
      ci <- FALSE

   if (ci && citype == "wald" && x$type != "stepfun" && sum(!x$delta.fix) >= 2L)
      stop(mstyle$stop("Cannot compute Wald-type confidence intervals for this selection model."))

   if (ci) {

      if (citype == "boot") {

         if (!is.null(ddd$seed))
            set.seed(ddd$seed)

         vd <- x$vd
         vd.na <- is.na(diag(vd))
         vd[vd.na,] <- 0
         vd[,vd.na] <- 0

         dsim <- .mvrnorm(reps, mu=delta, Sigma=vd)

         for (j in seq_len(ncol(dsim))) {
            dsim[,j] <- ifelse(dsim[,j] < x$delta.min[j], x$delta.min[j], dsim[,j])
            dsim[,j] <- ifelse(dsim[,j] > x$delta.max[j], x$delta.max[j], dsim[,j])
         }

         ys.ci <- lapply(ps, function(p) {
            ysim <- apply(dsim, 1, function(d) wi.fun(p, delta=d, yi=x$yi, vi=x$vi, preci=prec, alternative=x$alternative, steps=x$steps))
            quantile(ysim, probs=c(x$level/2, 1 - x$level/2))
         })
         ys.ci <- do.call(rbind, ys.ci)
         ys.lb <- ys.ci[,1]
         ys.ub <- ys.ci[,2]

      }

      if (citype == "wald") {

         ci.lb.delta <- x$ci.lb.delta
         ci.ub.delta <- x$ci.ub.delta

         if (x$type == "stepfun") {
            ci.lb.delta[x$delta.fix] <- delta[x$delta.fix]
            ci.ub.delta[x$delta.fix] <- delta[x$delta.fix]
         }

         ys.lb <- wi.fun(ps, delta=ci.lb.delta, yi=x$yi, vi=x$vi, preci=prec, alternative=x$alternative, steps=x$steps)
         ys.ub <- wi.fun(ps, delta=ci.ub.delta, yi=x$yi, vi=x$vi, preci=prec, alternative=x$alternative, steps=x$steps)

      }

   } else {

      ys.lb <- NA_real_
      ys.ub <- NA_real_

   }

   if (scale) {

      #is.inf.pos <- ys ==  Inf
      #is.inf.neg <- ys == -Inf

      ys[is.infinite(ys)] <- NA_real_

      rng.ys <- max(ys, na.rm=TRUE) - min(ys, na.rm=TRUE)
      min.ys <- min(ys, na.rm=TRUE)

      if (rng.ys > .Machine$double.eps^0.5) {
         ys <- (ys - min.ys) / rng.ys
         ys.lb <- (ys.lb - min.ys) / rng.ys
         ys.ub <- (ys.ub - min.ys) / rng.ys
      }

      #ys[is.inf.pos] <- 1
      #ys[is.inf.neg] <- 0

   }

   ys[ys < 0] <- 0
   ys.lb[ys.lb < 0] <- 0
   ys.ub[ys.ub < 0] <- 0

   if (missing(ylim)) {

      if (is.element(x$type, c("halfnorm", "negexp", "logistic", "power", "negexppow", "halfnorm2", "negexp2", "logistic2", "power2"))) {

         ylim <- c(0,1)

      } else {

         if (ci) {
            ylim <- c(min(c(ys.lb[is.finite(ys.lb)], ys[is.finite(ys)]), na.rm=TRUE), max(c(ys.ub[is.finite(ys.ub)], ys[is.finite(ys)]), na.rm=TRUE))
         } else {
            ylim <- range(ys[is.finite(ys)], na.rm=TRUE)
         }

      }

   } else {

      if (length(ylim) != 2L)
         stop(mstyle$stop("Argument 'ylim' should be a vector of length 2."))

      ylim <- sort(ylim)

   }

   if (!add)
      lplot(ps, ys, ylim=ylim, type="n", lwd=lwd, xlab="p-value", ylab="Relative Likelihood of Selection", ...)

   if (ci) {
      if (shade) {
         tmp <- approx(ps, ys.lb, n=10000, method="constant", f=1)
         ps.int.lb <- tmp$x
         ys.lb.int.lb <- tmp$y
         tmp <- approx(ps, ys.ub, n=10000, method="constant", f=1)
         ps.int.ub <- tmp$x
         ys.lb.int.ub <- tmp$y
         lpolygon(c(ps.int.lb,rev(ps.int.ub)), c(ys.lb.int.lb,rev(ys.lb.int.ub)), col=shadecol, border=NA)
         #lpolygon(c(ps,rev(ps)), c(ys.lb,rev(ys.ub)), col=shadecol, border=NA)
      }
      llines(ps, ys.lb, type=plot.type, lty=lty[2], lwd=lwd[2], ...)
      llines(ps, ys.ub, type=plot.type, lty=lty[2], lwd=lwd[2], ...)
   }

   if (rug && !add)
      lrug(x$pvals, quiet=TRUE)

   llines(ps, ys, type=plot.type, lty=lty[1], lwd=lwd[1], ...)

   sav <- data.frame(xs=ps, ys=ys, ys.lb=ys.lb, ys.ub=ys.ub)

   invisible(sav)

}
