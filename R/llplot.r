llplot <- function(measure, yi, vi, sei, ai, bi, ci, di, n1i, n2i, data, subset, drop00=TRUE,
xvals=1000, xlim, ylim, xlab, ylab, scale=TRUE,
lty, lwd, col, level=99.99, refline=0, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   ### data setup

   if (missing(measure))
      stop(mstyle$stop("Must specify an effect size or outcome measure via the 'measure' argument."))

   if (inherits(measure, "rma"))
      stop(mstyle$stop("Function not applicable to 'rma' objects."))

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

   level <- ifelse(level == 0, 1, ifelse(level >= 1, (100-level)/100, ifelse(level > .5, 1-level, level)))

   ### get ... argument

   ddd <- list(...)

   ### set defaults or get onlyo1, addyi, and addvi arguments

   onlyo1 <- ifelse(is.null(ddd$onlyo1), FALSE, ddd$onlyo1)
   addyi  <- ifelse(is.null(ddd$addyi),  TRUE,  ddd$addyi)
   addvi  <- ifelse(is.null(ddd$addvi),  TRUE,  ddd$addvi)

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
   mf.subset <- mf[[match("subset", names(mf))]]
   mf.lty    <- mf[[match("lty",    names(mf))]]
   mf.lwd    <- mf[[match("lwd",    names(mf))]]
   mf.col    <- mf[[match("col",    names(mf))]]
   subset <- eval(mf.subset, data, enclos=sys.frame(sys.parent()))
   lty    <- eval(mf.lty,    data, enclos=sys.frame(sys.parent()))
   lwd    <- eval(mf.lwd,    data, enclos=sys.frame(sys.parent()))
   col    <- eval(mf.col,    data, enclos=sys.frame(sys.parent()))

   if (measure == "GEN") {

      mf.yi  <- mf[[match("yi",  names(mf))]]
      mf.vi  <- mf[[match("vi",  names(mf))]]
      mf.sei <- mf[[match("sei", names(mf))]]
      yi  <- eval(mf.yi,  data, enclos=sys.frame(sys.parent()))
      vi  <- eval(mf.vi,  data, enclos=sys.frame(sys.parent()))
      sei <- eval(mf.sei, data, enclos=sys.frame(sys.parent()))

      if (is.null(vi)) {
         if (is.null(sei)) {
            stop(mstyle$stop("Need to specify 'vi' or 'sei' argument."))
         } else {
            vi <- sei^2
         }
      }

      if (length(yi)==0L || length(vi)==0L)
         stop(mstyle$stop("Cannot extract outcomes. Check that all of the required \n  information is specified via the appropriate arguments."))

      if (length(yi) != length(vi))
         stop(mstyle$stop("Supplied data vectors are not all of the same length."))

      k <- length(yi) ### number of outcomes before subsetting

      ### subsetting

      if (!is.null(subset)) {
         yi <- yi[subset]
         vi <- vi[subset]
      }

   }

   if (measure == "OR") {

      mf.ai  <- mf[[match("ai",  names(mf))]]
      mf.bi  <- mf[[match("bi",  names(mf))]]
      mf.ci  <- mf[[match("ci",  names(mf))]]
      mf.di  <- mf[[match("di",  names(mf))]]
      mf.n1i <- mf[[match("n1i", names(mf))]]
      mf.n2i <- mf[[match("n2i", names(mf))]]
      ai  <- eval(mf.ai,  data, enclos=sys.frame(sys.parent()))
      bi  <- eval(mf.bi,  data, enclos=sys.frame(sys.parent()))
      ci  <- eval(mf.ci,  data, enclos=sys.frame(sys.parent()))
      di  <- eval(mf.di,  data, enclos=sys.frame(sys.parent()))
      n1i <- eval(mf.n1i, data, enclos=sys.frame(sys.parent()))
      n2i <- eval(mf.n2i, data, enclos=sys.frame(sys.parent()))
      if (is.null(bi)) bi <- n1i - ai
      if (is.null(di)) di <- n2i - ci

      if (length(ai)==0L || length(bi)==0L || length(ci)==0L || length(di)==0L)
         stop(mstyle$stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments."))

      if (!all(length(ai) == c(length(ai),length(bi),length(ci),length(di))))
         stop(mstyle$stop("Supplied data vectors are not all of the same length."))

      if (any(c(ai > n1i, ci > n2i), na.rm=TRUE))
         stop(mstyle$stop("One or more event counts are larger than the corresponding group sizes."))

      if (any(c(ai, bi, ci, di) < 0, na.rm=TRUE))
         stop(mstyle$stop("One or more counts are negative."))

      k <- length(ai) ### number of outcomes before subsetting

      ### note studies that have at least one zero cell

      id0 <- c(ai == 0L | bi == 0L | ci == 0L | di == 0L)
      id0[is.na(id0)] <- FALSE

      ### note studies that have no events or all events

      id00 <- c(ai == 0L & ci == 0L) | c(bi == 0L & di == 0L)
      id00[is.na(id00)] <- FALSE

      ### if drop00=TRUE, set counts to NA for studies that have no events (or all events) in both arms

      if (drop00) {
         ai[id00] <- NA
         bi[id00] <- NA
         ci[id00] <- NA
         di[id00] <- NA
      }

      ### subsetting

      if (!is.null(subset)) {
         ai <- ai[subset]
         bi <- bi[subset]
         ci <- ci[subset]
         di <- di[subset]
      }

      dat <- escalc(measure="OR", ai=ai, bi=bi, ci=ci, di=di, drop00=drop00, onlyo1=onlyo1, addyi=addyi, addvi=addvi)

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
      ids  <- ids[subset]
      lty  <- lty[subset]
      lwd  <- lwd[subset]
      col  <- col[subset]
      id0  <- id0[subset]
      id00 <- id00[subset]
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
         ai   <- ai[not.na]
         bi   <- bi[not.na]
         ci   <- ci[not.na]
         di   <- di[not.na]
         ids  <- ids[not.na]
         lty  <- lty[not.na]
         lwd  <- lwd[not.na]
         col  <- col[not.na]
         id0  <- id0[not.na]
         id00 <- id00[not.na]
         k    <- length(ai)
         warning(mstyle$warning("Studies with NAs omitted from plotting."), call.=FALSE)
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

   ### set default line color (gray0 to gray60 according to the rank of vi)

   if (is.null(col))
      col <- paste0("gray", round(seq(from=0, to=60, length.out=k))[rank(vi)])

   ### set x axis limits

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
      }
      #lls.sum <- rowSums(lls, na.rm=TRUE)
      lls <- apply(lls, 2, "/", lls.sum)
   }

   lls[out] <- NA

   ### set y axis limits

   if (missing(ylim)) {
      ylim <- c(0, max(lls, na.rm=TRUE))
   } else {
      ylim <- sort(ylim)
   }

   plot(NA, NA, xlim=c(xlim[1], xlim[2]), ylim=ylim, xlab=xlab, ylab=ylab, ...)

   for (i in seq_len(k)[order(1/vi)]) {
      lines(xs, lls[i,], lty=lty[i], lwd=lwd[i], col=col[i], ...)
   }

   if (is.numeric(refline))
      abline(v=refline, lty="solid", lwd=2, ...)

   invisible(lls)

}
