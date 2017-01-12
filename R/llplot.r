llplot <- function(measure="OR", ai, bi, ci, di, n1i, n2i, data, subset, drop00=TRUE,
xvals=1000, xlim, ylim, xlab, ylab, scale=TRUE,
lty, lwd, col, level=99.99, refline=0, ...) {

   #########################################################################

   ### data setup

   if (!is.element(measure, "OR"))
      stop("Currently only measure=\"OR\" can be specified.")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop("Unknown 'na.action' specified under options().")

   if (!requireNamespace("BiasedUrn", quietly=TRUE))
      stop("Please install the 'BiasedUrn' package to use this function.")

   if (missing(xlab))
      xlab <- "Log Odds Ratio"

   if (missing(ylab)) {
      if (scale) {
         ylab <- "Scaled Log Likelihood"
      } else {
         ylab <- "Log Likelihood"
      }
   }

   level <- ifelse(level > 1, (100-level)/100, ifelse(level > .5, 1-level, level))

   #########################################################################

   ### check if data argument has been specified

   if (missing(data))
      data <- NULL

   if (is.null(data)) {
      data <- sys.frame(sys.parent())
   } else {
      if (!is.data.frame(data)) {
         data <- data.frame(data)
      }
   }

   ### extract slab, subset, and mods values, possibly from the data frame specified via data (arguments not specified are NULL)

   mf <- match.call()
   mf.subset <- mf[[match("subset", names(mf))]]
   mf.lty    <- mf[[match("lty",   names(mf))]]
   mf.lwd    <- mf[[match("lwd",   names(mf))]]
   mf.col    <- mf[[match("col",   names(mf))]]
   subset <- eval(mf.subset, data, enclos=sys.frame(sys.parent()))
   lty    <- eval(mf.lty,   data, enclos=sys.frame(sys.parent()))
   lwd    <- eval(mf.lwd,   data, enclos=sys.frame(sys.parent()))
   col    <- eval(mf.col,   data, enclos=sys.frame(sys.parent()))

   mf.ai   <- mf[[match("ai",  names(mf))]]
   mf.bi   <- mf[[match("bi",  names(mf))]]
   mf.ci   <- mf[[match("ci",  names(mf))]]
   mf.di   <- mf[[match("di",  names(mf))]]
   mf.n1i  <- mf[[match("n1i", names(mf))]]
   mf.n2i  <- mf[[match("n2i", names(mf))]]
   ai      <- eval(mf.ai,  data, enclos=sys.frame(sys.parent()))
   bi      <- eval(mf.bi,  data, enclos=sys.frame(sys.parent()))
   ci      <- eval(mf.ci,  data, enclos=sys.frame(sys.parent()))
   di      <- eval(mf.di,  data, enclos=sys.frame(sys.parent()))
   n1i     <- eval(mf.n1i, data, enclos=sys.frame(sys.parent()))
   n2i     <- eval(mf.n2i, data, enclos=sys.frame(sys.parent()))
   if (is.null(bi)) bi <- n1i - ai
   if (is.null(di)) di <- n2i - ci

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

   dat <- escalc(measure="OR", ai=ai, bi=bi, ci=ci, di=di, drop00=drop00)

   yi <- dat$yi ### one or more yi/vi pairs may be NA/NA
   vi <- dat$vi ### one or more yi/vi pairs may be NA/NA

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
            stop("Length of 'lty' argument does not match data.")
      }
   }

   if (!is.null(lwd)) {
      if (length(lwd) == 1L) {
         lwd <- rep(lwd, k)
      } else {
         if (length(lwd) != k)
            stop("Length of 'lwd' argument does not match data.")
      }
   }

   if (!is.null(col)) {
      if (length(col) == 1L) {
         col <- rep(col, k)
      } else {
         if (length(col) != k)
            stop("Length of 'col' argument does not match data.")
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

   k <- length(ai)

   ### check for NAs and act accordingly

   aibicidi.na <- is.na(ai) | is.na(bi) | is.na(ci) | is.na(di)

   if (any(aibicidi.na)) {

      not.na <- !aibicidi.na

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
         warning("Studies with NAs omitted from plotting.")
      }

      if (na.act == "na.fail")
         stop("Missing values in studies.")

   } else {
      not.na <- rep(TRUE, k)
   }

   ### at least one study left?

   if (k < 1)
      stop("Processing terminated since k = 0.")

   #########################################################################

   ### set default line types (id0 studies = dashed line, id00 studies = dotted line, all others = solid line)

   if (is.null(lty))
      lty <- ifelse(id0 | id00, ifelse(id00, "dotted", "dashed"), "solid")

   ### set default line widths (4.0 to 0.2 according to the rank of vi)

   if (is.null(lwd))
      lwd <- seq(from=4.0, to=0.2, length=k)[rank(vi)]

   ### set default line color (gray0 to gray80 according to the rank of vi)

   if (is.null(col))
      col <- paste0("gray", round(seq(from=0, to=80, length=k))[rank(vi)])

   ### set x axis limits

   ci.lb <- yi - qnorm(level/2, lower.tail=FALSE) * sqrt(vi)
   ci.ub <- yi + qnorm(level/2, lower.tail=FALSE) * sqrt(vi)

   if (missing(xlim)) {
      xlim <- c(min(ci.lb, na.rm=TRUE),max(ci.ub, na.rm=TRUE))
   } else {
      xlim <- sort(xlim)
   }

   logORs <- seq(from=xlim[1], to=xlim[2], length.out=xvals)

   lls <- matrix(NA_real_, nrow=k, ncol=xvals)
   out <- matrix(TRUE, nrow=k, ncol=xvals)

   for (i in seq_len(k)) {
      for (j in seq_len(xvals)) {
         lls[i,j] <- .dnchgi(logORs[j], ai=ai[i], bi=bi[i], ci=ci[i], di=di[i], random=FALSE, dnchgcalc="dFNCHypergeo", dnchgprec=1e-10)
         if (logORs[j] >= ci.lb[i] & logORs[j] <= ci.ub[i])
            out[i,j] <- FALSE
      }
   }

   if (scale) {
      trapezoid <- function(x,y) sum(diff(x)*(y[-1]+y[-length(y)]))/2
      lls.sum <- rep(NA_real_, k)
      for (i in seq_len(k)) {
         lls.sum[i] <- trapezoid(logORs[!is.na(lls[i,])], lls[i,!is.na(lls[i,])])
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
      lines(logORs, lls[i,], lty=lty[i], lwd=lwd[i], col=col[i], ...)
   }

   if (is.numeric(refline))
      abline(v=refline, lty="solid", lwd=2, ...)

   invisible(lls)

}
