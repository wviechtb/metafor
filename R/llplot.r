llplot <- function(measure, yi, vi, sei, ai, bi, ci, di, n1i, n2i, data, subset, drop00=TRUE,
xvals=1000, xlim, ylim, xlab, ylab, scale=TRUE,
lty, lwd, col, level=99.99, refline=0, ...) {

   #########################################################################

   mstyle <- .get.mstyle()

   ### data setup

   if (missing(measure))
      stop(mstyle$stop("Must specify an effect size or outcome measure via the 'measure' argument."))

   .chkclass(class(measure), notap="rma", type="Function")

   if (!is.element(measure, c("GEN", "OR")))
      stop(mstyle$stop("Currently only measure=\"GEN\" or measure=\"OR\" can be specified."))

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   if (measure == "OR" && !requireNamespace("BiasedUrn", quietly=TRUE))
      stop(mstyle$stop("Please install the 'BiasedUrn' package to use this function."))

   if (missing(xlab)) {
      if (measure == "GEN")
         xlab <- "Observed Outcome"
      if (measure == "OR")
         xlab <- "Log Odds Ratio"
   }

   if (missing(ylab)) {
      if (scale) {
         ylab <- "Scaled Likelihood"
      } else {
         ylab <- "Likelihood"
      }
   }

   level <- .level(level)

   ### get ... argument

   ddd <- list(...)

   ### set defaults or get onlyo1, addyi, and addvi arguments

   onlyo1 <- .chkddd(ddd$onlyo1, FALSE)
   addyi  <- .chkddd(ddd$addyi,  TRUE)
   addvi  <- .chkddd(ddd$addvi,  TRUE)

   .start.plot()

   #########################################################################

   ### check if data argument has been specified

   if (missing(data))
      data <- NULL

   if (is.null(data)) {
      data <- sys.frame(sys.parent())
   } else {
      if (!is.data.frame(data))
         data <- data.frame(data)
   }

   ### extract values, possibly from the data frame specified via data (arguments not specified are NULL)

   mf <- match.call()

   subset  <- .getx("subset", mf=mf, data=data)
   lty     <- .getx("lty",    mf=mf, data=data)
   lwd     <- .getx("lwd",    mf=mf, data=data)
   col     <- .getx("col",    mf=mf, data=data)

   if (measure == "GEN") {

      yi  <- .getx("yi",  mf=mf, data=data, checknumeric=TRUE)
      vi  <- .getx("vi",  mf=mf, data=data, checknumeric=TRUE)
      sei <- .getx("sei", mf=mf, data=data, checknumeric=TRUE)

      if (is.null(vi)) {
         if (is.null(sei)) {
            stop(mstyle$stop("Must specify 'vi' or 'sei' argument."))
         } else {
            vi <- sei^2
         }
      }

      if (!.all.specified(yi, vi))
         stop(mstyle$stop("Cannot construct plot. Check that all of the required information is specified\n  via the appropriate arguments (i.e., yi, vi)."))

      if (!.equal.length(yi, vi))
         stop(mstyle$stop("Supplied data vectors are not all of the same length."))

      k <- length(yi) ### number of outcomes before subsetting

      ### subsetting

      if (!is.null(subset)) {
         subset <- .chksubset(subset, k)
         yi <- .getsubset(yi, subset)
         vi <- .getsubset(vi, subset)
      }

   }

   if (measure == "OR") {

      ai  <- .getx("ai",  mf=mf, data=data, checknumeric=TRUE)
      bi  <- .getx("bi",  mf=mf, data=data, checknumeric=TRUE)
      ci  <- .getx("ci",  mf=mf, data=data, checknumeric=TRUE)
      di  <- .getx("di",  mf=mf, data=data, checknumeric=TRUE)
      n1i <- .getx("n1i", mf=mf, data=data, checknumeric=TRUE)
      n2i <- .getx("n2i", mf=mf, data=data, checknumeric=TRUE)

      if (!.equal.length(ai, bi, ci, di, n1i, n2i))
         stop(mstyle$stop("Supplied data vectors are not all of the same length."))

      n1i.inc <- n1i != ai + bi
      n2i.inc <- n2i != ci + di

      if (any(n1i.inc, na.rm=TRUE))
         stop(mstyle$stop("One or more 'n1i' values are not equal to 'ai + bi'."))
      if (any(n2i.inc, na.rm=TRUE))
         stop(mstyle$stop("One or more 'n2i' values are not equal to 'ci + di'."))

      bi <- replmiss(bi, n1i-ai)
      di <- replmiss(di, n2i-ci)

      if (!.all.specified(ai, bi, ci, di))
         stop(mstyle$stop("Cannot construct plot. Check that all of the required information is specified\n  via the appropriate arguments (i.e., ai, bi, ci, di or ai, n1i, ci, n2i)."))

      n1i <- ai + bi
      n2i <- ci + di

      if (any(c(ai > n1i, ci > n2i), na.rm=TRUE))
         stop(mstyle$stop("One or more event counts are larger than the corresponding group sizes."))

      if (any(c(ai, bi, ci, di) < 0, na.rm=TRUE))
         stop(mstyle$stop("One or more counts are negative."))

      if (any(c(n1i < 0, n2i < 0), na.rm=TRUE))
         stop(mstyle$stop("One or more group sizes are negative."))

      k <- length(ai) ### number of outcomes before subsetting

      ### note studies that have at least one zero cell

      id0 <- c(ai == 0L | bi == 0L | ci == 0L | di == 0L)
      id0[is.na(id0)] <- FALSE

      ### note studies that have no events or all events

      id00 <- c(ai == 0L & ci == 0L) | c(bi == 0L & di == 0L)
      id00[is.na(id00)] <- FALSE

      ### if drop00=TRUE, set counts to NA for studies that have no events (or all events) in both arms

      if (drop00) {
         ai[id00] <- NA_real_
         bi[id00] <- NA_real_
         ci[id00] <- NA_real_
         di[id00] <- NA_real_
      }

      ### subsetting

      if (!is.null(subset)) {
         subset <- .chksubset(subset, k)
         ai <- .getsubset(ai, subset)
         bi <- .getsubset(bi, subset)
         ci <- .getsubset(ci, subset)
         di <- .getsubset(di, subset)
      }

      dat <- .do.call(escalc, measure="OR", ai=ai, bi=bi, ci=ci, di=di, drop00=drop00, onlyo1=onlyo1, addyi=addyi, addvi=addvi)

      yi <- dat$yi ### one or more yi/vi pairs may be NA/NA
      vi <- dat$vi ### one or more yi/vi pairs may be NA/NA

   }

   #########################################################################

   ### study ids (1:k sequence before subsetting)

   ids <- seq_len(k)

   ### setting of lty, lwd, and col arguments (if a single value, repeat k times)
   ### if any of these arguments is not a single value, it must have the same length as the data before subsetting

   if (!is.null(lty)) {
      if (length(lty) == 1L) {
         lty <- rep(lty, k)
      } else {
         if (length(lty) != k)
            stop(mstyle$stop(paste0("Length of 'lty' argument (", length(lty), ") does not match length of data (", k, ").")))
      }
   }

   if (!is.null(lwd)) {
      if (length(lwd) == 1L) {
         lwd <- rep(lwd, k)
      } else {
         if (length(lwd) != k)
         stop(mstyle$stop(paste0("Length of 'lwd' argument (", length(lwd), ") does not match length of data (", k, ").")))
      }
   }

   if (!is.null(col)) {
      if (length(col) == 1L) {
         col <- rep(col, k)
      } else {
         if (length(col) != k)
            stop(mstyle$stop(paste0("Length of 'col' argument (", length(col), ") does not match length of data (", k, ").")))
      }
   }

   ### if a subset of studies is specified

   if (!is.null(subset)) {
      ids  <- .getsubset(ids,  subset)
      lty  <- .getsubset(lty,  subset)
      lwd  <- .getsubset(lwd,  subset)
      col  <- .getsubset(col,  subset)
      id0  <- .getsubset(id0,  subset)
      id00 <- .getsubset(id00, subset)
   }

   ### number of outcomes after subsetting

   k <- length(yi)

   ### check for NAs and act accordingly

   if (measure == "GEN") {
      has.na <- is.na(yi) | is.na(vi)
   }
   if (measure == "OR") {
      has.na <- is.na(ai) | is.na(bi) | is.na(ci) | is.na(di)
   }

   not.na <- !has.na

   if (any(has.na)) {

      if (na.act == "na.omit" || na.act == "na.exclude" || na.act == "na.pass") {
         yi   <- yi[not.na]
         vi   <- vi[not.na]
         if (measure == "OR") {
            ai   <- ai[not.na]
            bi   <- bi[not.na]
            ci   <- ci[not.na]
            di   <- di[not.na]
            id0  <- id0[not.na]
            id00 <- id00[not.na]
         }
         k    <- length(yi)
         ids  <- ids[not.na]
         lty  <- lty[not.na]
         lwd  <- lwd[not.na]
         col  <- col[not.na]
         warning(mstyle$warning(paste(sum(has.na), ifelse(sum(has.na) > 1, "studies", "study"), "with NAs omitted from plotting.")), call.=FALSE)
      }

      if (na.act == "na.fail")
         stop(mstyle$stop("Missing values in studies."))

   }

   ### at least one study left?

   if (k < 1)
      stop(mstyle$stop("Processing terminated since k = 0."))

   #########################################################################

   ### set default line types (id0 studies = dashed line, id00 studies = dotted line, all others = solid line)

   if (measure == "GEN") {
      if (is.null(lty))
         lty <- rep("solid", k)
   }

   if (measure == "OR") {
      if (is.null(lty))
         lty <- ifelse(id0 | id00, ifelse(id00, "dotted", "dashed"), "solid")
   }

   ### set default line widths (4.0 to 0.4 according to the rank of vi)

   if (is.null(lwd))
      lwd <- seq(from=4.0, to=0.4, length.out=k)[rank(vi)]

   ### set default line color (darker to lighter according to the rank of vi)

   if (is.null(col)) {
      col <- sapply(seq(from=0.8, to=0.2, length.out=k), function(x) .coladj(par("bg","fg"), dark=x, light=-x))
      col <- col[rank(vi)]
   }

   ### set x-axis limits

   ci.lb <- yi - qnorm(level/2, lower.tail=FALSE) * sqrt(vi)
   ci.ub <- yi + qnorm(level/2, lower.tail=FALSE) * sqrt(vi)

   if (missing(xlim)) {
      xlim <- c(min(ci.lb, na.rm=TRUE),max(ci.ub, na.rm=TRUE))
   } else {
      xlim <- sort(xlim)
   }

   xs <- seq(from=xlim[1], to=xlim[2], length.out=xvals)

   lls <- matrix(NA_real_, nrow=k, ncol=xvals)
   out <- matrix(TRUE, nrow=k, ncol=xvals)

   if (measure == "GEN") {

      for (i in seq_len(k)) {
         for (j in seq_len(xvals)) {
            lls[i,j] <- dnorm(yi[i], xs[j], sqrt(vi[i]))
            if (xs[j] >= ci.lb[i] & xs[j] <= ci.ub[i])
               out[i,j] <- FALSE
         }
      }

   }

   if (measure == "OR") {

      for (i in seq_len(k)) {
         for (j in seq_len(xvals)) {
            lls[i,j] <- .dnchgi(xs[j], ai=ai[i], bi=bi[i], ci=ci[i], di=di[i], random=FALSE, dnchgcalc="dFNCHypergeo", dnchgprec=1e-10)
            if (xs[j] >= ci.lb[i] & xs[j] <= ci.ub[i])
               out[i,j] <- FALSE
         }
      }

   }

   if (scale) {
      trapezoid <- function(x,y) sum(diff(x)*(y[-1]+y[-length(y)]))/2
      lls.sum <- rep(NA_real_, k)
      for (i in seq_len(k)) {
         lls.sum[i] <- trapezoid(xs[!is.na(lls[i,])], lls[i,!is.na(lls[i,])])
         lls[i,] <- lls[i,] / lls.sum[i]
      }
   }

   lls[out] <- NA_real_

   ### set y-axis limits

   if (missing(ylim)) {
      ylim <- c(0, max(lls, na.rm=TRUE))
   } else {
      ylim <- sort(ylim)
   }

   plot(NA, NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)

   for (i in seq_len(k)[order(1/vi)]) {
      lines(xs, lls[i,], lty=lty[i], lwd=lwd[i], col=col[i], ...)
   }

   if (is.numeric(refline))
      abline(v=refline, lty="solid", lwd=2, ...)

   invisible(lls)

}
