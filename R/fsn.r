fsn <- function(yi, vi, sei, data, type="Rosenthal", alpha=.05, target, weighted=FALSE, subset, digits, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   type <- match.arg(type, c("Rosenthal", "Orwin", "Rosenberg", "Binomial", "REM"))

   if (missing(target))
      target <- NULL

   ### set defaults for digits

   if (missing(digits)) {
      digits <- .set.digits(dmiss=TRUE)
   } else {
      digits <- .set.digits(digits, dmiss=FALSE)
   }

   ddd <- list(...)

   .chkdots(ddd, c("test", "verbose", "interval", "iters"))

   if (is.null(ddd$test)) {
      test <- "Stouffer"
   } else {
      test <- match.arg(ddd$test, c("Stouffer", "Fisher"))
   }

   if (is.null(ddd$verbose)) {
      verbose <- FALSE
   } else {
      verbose <- ddd$verbose
   }

   if (is.null(ddd$interval)) {
      interval <- c(0,1000)
   } else {
      interval <- ddd$interval
   }

   if (is.null(ddd$iters)) {
      iters <- 100000
   } else {
      iters <- ddd$iters
   }

   meanes  <- NA
   pval    <- NA
   rejrate <- NA

   #########################################################################

   ###### data setup

   ### check if data argument has been specified

   if (missing(data))
      data <- NULL

   if (is.null(data)) {
      data <- sys.frame(sys.parent())
   } else {
      if (!is.data.frame(data))
         data <- data.frame(data)
   }

   mf <- match.call()

   yi     <- .getx("yi",      mf=mf, data=data, checknumeric=TRUE)
   vi     <- .getx("vi",      mf=mf, data=data, checknumeric=TRUE)
   sei    <- .getx("sei",     mf=mf, data=data, checknumeric=TRUE)
   #weight <- .getx("weights", mf=mf, data=data, checknumeric=TRUE)
   subset <- .getx("subset",  mf=mf, data=data, checknumeric=TRUE)

   if (type %in% c("Rosenthal", "Rosenberg", "REM") || (type == "Orwin" && weighted)) {
      if (is.null(vi)) {
         if (is.null(sei)) {
            stop(mstyle$stop("Must specify 'vi' or 'sei' argument."))
         } else {
            vi <- sei^2
         }
      }
   } else {
      vi <- rep(0, length(yi))
   }

   ### check length of yi and vi

   if (length(yi) != length(vi))
      stop(mstyle$stop("Length of 'yi' and 'vi' (or 'sei') is not the same."))

   ### check 'vi' argument for potential misuse

   .chkviarg(mf$vi)

   ### if a subset of studies is specified

   if (!is.null(subset)) {
      subset <- .chksubset(subset, length(yi))
      yi <- .getsubset(yi, subset)
      vi <- .getsubset(vi, subset)
   }

   ### check for NAs in yi/vi and act accordingly

   yivi.na <- is.na(yi) | is.na(vi)

   if (any(yivi.na)) {

      not.na <- !yivi.na

      if (na.act == "na.omit" || na.act == "na.exclude" || na.act == "na.pass") {
         yi <- yi[not.na]
         vi <- vi[not.na]
      }

      if (na.act == "na.fail")
         stop(mstyle$stop("Missing values in results."))

   }

   #########################################################################

   if (type == "Rosenthal" && test == "Stouffer") {

      k      <- length(yi)
      zi     <- yi / sqrt(vi)
      z.avg  <- abs(sum(zi) / sqrt(k))
      pval   <- pnorm(z.avg, lower.tail=FALSE)
      fsnum  <- max(0, k * (z.avg / qnorm(alpha, lower.tail=FALSE))^2 - k)
      target <- NA

   }

   if (type == "Rosenthal" && test == "Fisher") {

      zi <- c(yi / sqrt(vi))
      pi <- pnorm(abs(zi), lower.tail=FALSE)
      pval <- .fsn.fisher(0, pi=pi, alpha=0)

      if (pval >= alpha) {
         fsnum <- 0
      } else {
         fsnum <- try(uniroot(.fsn.fisher, interval=interval, extendInt="upX", pi=pi, alpha=alpha)$root, silent=FALSE)
         if (inherits(fsnum, "try-error"))
            stop(mstyle$stop("Could not find fail-safe N using Fisher's method."))
      }
      target <- NA

   }

   if (type == "Orwin") {

      k <- length(yi)

      if (weighted) {
         wi <- 1/vi
         meanes <- .wmean(yi, wi)
      } else {
         meanes <- mean(yi)
      }

      if (is.null(target))
         target <- meanes / 2

      if (identical(target, 0)) {
         fsnum <- Inf
      } else {
         if (sign(target) != sign(meanes))
            target <- -1 * target
         fsnum <- max(0, k * (meanes - target) / target)
      }

   }

   if (type == "Rosenberg") {

      k      <- length(yi)
      wi     <- 1/vi
      meanes <- .wmean(yi, wi)
      zval   <- meanes / sqrt(1/sum(wi))
      w.p    <- (sum(wi*yi) / qnorm(alpha/2, lower.tail=FALSE))^2 - sum(wi)
      pval   <- 2*pnorm(abs(zval), lower.tail=FALSE)
      fsnum  <- max(0, k*w.p/sum(wi))
      target <- NA

   }

   if (type == "Binomial") {

      k    <- length(yi)
      kpos <- sum(yi > 0)
      pval <- binom.test(kpos, k)$p.value
      if (pval >= alpha) {
         fsnum <- 0
      } else {
         pvalnew <- pval
         fsnum <- 0
         while (pvalnew < alpha) {
            fsnum <- fsnum + 2
            pvalnew <- binom.test(kpos + fsnum/2, k + fsnum)$p.value
         }
      }
      target <- NA

   }

   if (type == "REM") {

      res <- .fsn.fitre(yi, vi)

      vnew   <- 1/mean(1/vi)
      tau2   <- res$tau2
      meanes <- res$est
      pval   <- res$pval

      if (is.null(target))
         target <- meanes / 2

      if (identical(target, 0)) {

         fsnum <- Inf

      } else {

         if (sign(target) != sign(meanes))
            target <- -1 * target

         diff.lo <- .fsn.re(0, yi=yi, vi=vi, vnew=vnew, tau2=tau2, target=target, alpha=alpha, iters=iters)

         if ((meanes > 0 && diff.lo < 0) || (meanes < 0 && diff.lo > 0)) {
            fsnum <- 0
         } else {
            fsnum <- try(uniroot(.fsn.re, interval=interval, tol=.001, extendInt=ifelse(meanes > 0,"downX","upX"), yi=yi, vi=vi, vnew=vnew, tau2=tau2, target=target, alpha=alpha, iters=iters, verbose=verbose)$root, silent=TRUE)
            if (inherits(fsnum, "try-error"))
               stop(mstyle$stop("Could not find fail-safe N based on a random-effects model."))
         }

         rejrate <- .fsn.fitnew(fsnum, yi, vi, vnew, tau2, alpha, iters)$rejrate

      }

   }

   if (is.finite(fsnum) && abs(fsnum - round(fsnum)) >= .Machine$double.eps^0.5) {
      fsnum <- ceiling(fsnum)
   } else {
      fsnum <- round(fsnum)
   }

   #########################################################################

   res <- list(type=type, fsnum=fsnum, alpha=alpha, pval=pval, meanes=meanes, target=target, rejrate=rejrate, digits=digits)

   class(res) <- "fsn"
   return(res)

}
