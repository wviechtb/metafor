fsn <- function(x, vi, sei, subset, data, type, alpha=.05, target,
                method, exact=FALSE, verbose=FALSE, digits, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   ### set defaults

   if (missing(target))
      target <- NULL

   ddd <- list(...)

   .chkdots(ddd, c("pool", "mumiss", "interval", "maxint", "tol", "maxiter", "tau2", "test", "weighted"))

   if (is.null(ddd$pool)) {
      pool <- "stouffer"
   } else {
      pool <- match.arg(tolower(ddd$pool), c("stouffer", "fisher"))
   }

   if (is.null(ddd$mumiss)) {
      mumiss <- 0
   } else {
      mumiss <- ddd$mumiss
   }

   # note: default interval set below; see [a] (based on k)

   if (is.null(ddd$maxint)) {
      maxint <- 10^7
   } else {
      maxint <- ddd$maxint
   }

   if (is.null(ddd$tol)) {
      tol <- .Machine$double.eps^0.25
   } else {
      tol <- ddd$tol
   }

   if (is.null(ddd$maxiter)) {
      maxiter <- 1000
   } else {
      maxiter <- ddd$maxiter
   }

   ### observed values (to be replaced as needed)

   est  <- NA_real_ # pooled estimate
   tau2 <- NA_real_ # tau^2 estimate
   pval <- NA_real_ # p-value

   ### defaults (to be replaced for type="General")

   est.fsn  <- NA_real_
   tau2.fsn <- NA_real_
   pval.fsn <- NA_real_
   ub.sign  <- ""

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

   mf <- match.call()

   x <- .getx("x", mf=mf, data=data)

   #########################################################################

   if (inherits(x, "rma")) {

      .chkclass(class(x), must="rma", notav=c("robust.rma", "rma.glmm", "rma.mv", "rma.ls", "rma.gen", "rma.uni.selmodel"))

      if (!x$int.only)
         stop(mstyle$stop("Method only applicable to models without moderators."))

      if (!missing(type) && type != "General")
         warning(mstyle$warning("Setting type='General' when using fsn() on a model object. "), call.=FALSE)

      type <- "General"

      if (!is.null(x$weights))
         stop(mstyle$stop("Cannot use function on models with custom weights."))

      if (!missing(vi) || !missing(sei) || !missing(subset))
         warning(mstyle$warning("Arguments 'vi', 'sei', and 'subset' ignored when 'x' is a model object."), call.=FALSE)

      yi <- x$yi
      vi <- x$vi

      ### set defaults for digits

      if (missing(digits)) {
         digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
      } else {
         digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
      }

   } else {

      if (!.is.vector(x))
         stop(mstyle$stop("Argument 'x' must be a vector or an 'rma' model object."))

      ### select/match type

      if (missing(type))
         type <- "Rosenthal"

      type.options <- c("rosenthal", "binomial", "orwin", "rosenberg", "general")

      type <- type.options[grep(tolower(type), type.options)[1]]

      if (is.na(type))
         stop(mstyle$stop("Unknown 'type' specified."))

      type <- paste0(toupper(substr(type, 1, 1)), substr(type, 2, nchar(type)))

      ### check if yi is numeric

      yi <- x

      if (!is.numeric(yi))
         stop(mstyle$stop("The object/variable specified for the 'x' argument is not numeric."))

      ### set defaults for digits

      if (missing(digits)) {
         digits <- .set.digits(dmiss=TRUE)
      } else {
         digits <- .set.digits(digits, dmiss=FALSE)
      }

      vi     <- .getx("vi",     mf=mf, data=data, checknumeric=TRUE)
      sei    <- .getx("sei",    mf=mf, data=data, checknumeric=TRUE)
      subset <- .getx("subset", mf=mf, data=data)

      if (is.null(vi)) {
         if (!is.null(sei))
            vi <- sei^2
      }

      if (type %in% c("Rosenthal", "Rosenberg", "General") && is.null(vi))
         stop(mstyle$stop("Must specify 'vi' or 'sei' argument."))

      ### ensure backwards compatibility with the 'weighted' argument when type="Orwin"

      if (type == "Orwin") {
         if (isTRUE(ddd$weighted) && is.null(vi)) # if weighted=TRUE, then check that the vi's are available
            stop(mstyle$stop("Must specify 'vi' or 'sei' argument."))
         if (isFALSE(ddd$weighted)) # if weighted=FALSE, then set vi <- 1 for unweighted
            vi <- 1
         if (is.null(ddd$weighted) && is.null(vi)) # if weighted is unspecified, set vi <- 1 if vi's are unspecified
            vi <- 1
      }

      ### allow easy setting of vi to a single value

      if (length(vi) == 1L)
         vi <- rep(vi, length(yi))

      ### check length of yi and vi

      if (length(yi) != length(vi))
         stop(mstyle$stop("Length of 'yi' and 'vi' (or 'sei') is not the same."))

      ### check 'vi' argument for potential misuse

      .chkviarg(mf$vi)

      #########################################################################

      ### if a subset of studies is specified

      if (!is.null(subset)) {
         subset <- .chksubset(subset, length(yi))
         yi <- .getsubset(yi, subset)
         vi <- .getsubset(vi, subset)
      }

      ### check for NAs and act accordingly

      has.na <- is.na(yi) | is.na(vi)

      if (any(has.na)) {

         not.na <- !has.na

         if (na.act == "na.omit" || na.act == "na.exclude" || na.act == "na.pass") {

            yi <- yi[not.na]
            vi <- vi[not.na]
            warning(mstyle$warning(paste(sum(has.na), ifelse(sum(has.na) > 1, "studies", "study"), "with NAs omitted.")), call.=FALSE)

         }

         if (na.act == "na.fail")
            stop(mstyle$stop("Missing values in data."))

      }

   }

   #########################################################################

   ### check for non-positive sampling variances

   if (any(vi <= 0))
      stop(mstyle$stop("Cannot use function when there are non-positive sampling variances in the data."))

   ### number of studies

   k <- length(yi)

   if (k == 1)
      stop(mstyle$stop("Stopped because k = 1."))

   ### set interval for uniroot() [a]

   if (is.null(ddd$interval)) {
      interval <- c(0, k*50)
   } else {
      interval <- ddd$interval
   }

   #########################################################################

   if (type == "Rosenthal" && pool == "stouffer") {

      zi     <- c(yi / sqrt(vi))
      z.avg  <- abs(sum(zi) / sqrt(k))
      pval   <- pnorm(z.avg, lower.tail=FALSE)
      fsnum  <- max(0, k * (z.avg / qnorm(alpha, lower.tail=FALSE))^2 - k)
      target <- NA_real_

   }

   if (type == "Rosenthal" && pool == "fisher") {

      zi <- c(yi / sqrt(vi))
      pi <- pnorm(abs(zi), lower.tail=FALSE)
      pval <- .fsn.fisher(0, pi=pi, alpha=0)

      if (pval >= alpha) {
         fsnum <- 0
      } else {
         fsnum <- try(uniroot(.fsn.fisher, interval=interval, extendInt="upX", tol=tol, maxiter=maxiter, pi=pi, alpha=alpha)$root, silent=FALSE)
         if (inherits(fsnum, "try-error"))
            stop(mstyle$stop("Could not find fail-safe N using Fisher's method for pooling p-values."))
      }
      target <- NA_real_

   }

   if (type == "Binomial") {

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
      target <- NA_real_

   }

   if (type == "Orwin") {

      wi <- 1 / vi
      est <- .wmean(yi, wi)

      if (is.null(target))
         target <- est / 2

      if (identical(target, 0)) {
         fsnum <- Inf
      } else {
         if (sign(target) != sign(est))
            target <- -1 * target
         fsnum <- max(0, k * (est - target) / target)
      }

   }

   if (type == "Rosenberg") {

      wi     <- 1 / vi
      est    <- .wmean(yi, wi)
      zval   <- est / sqrt(1/sum(wi))
      pval   <- 2*pnorm(abs(zval), lower.tail=FALSE)
      vt     <- 1 / mean(1/vi)
      #w.p   <- (sum(wi*yi) / qnorm(alpha/2, lower.tail=FALSE))^2 - sum(wi)
      #fsnum <- max(0, k*w.p/sum(wi))
      fsnum  <- max(0, ((sum(wi*yi) / qnorm(alpha/2, lower.tail=FALSE))^2 - sum(wi)) * vt)
      target <- NA_real_

   }

   if (type == "General") {

      if (missing(method)) {
         if (inherits(x, "rma")) {
            method <- x$method
         } else {
            method <- "REML"
         }
      }

      tau2fix <- NULL

      if (inherits(x, "rma") && x$tau2.fix)
         tau2fix <- x$tau2

      if (!is.null(ddd$tau2))
         tau2fix <- ddd$tau2

      test <- "z"

      if (inherits(x, "rma"))
         test <- x$test

      if (!is.null(ddd$test))
         test <- ddd$test

      if (test != "z")
         exact <- TRUE

      weighted <- TRUE

      if (inherits(x, "rma"))
         weighted <- x$weighted

      if (!is.null(ddd$weighted))
         weighted <- isTRUE(ddd$weighted)

      tmp <- try(rma(yi, vi, method=method, tau2=tau2fix, test=test, weighted=weighted, verbose=verbose), silent=!verbose)

      if (inherits(tmp, "try-error"))
         stop(mstyle$stop("Could not fit random-effects model (use verbose=TRUE for more info)."), call.=FALSE)

      vt   <- 1 / mean(1/vi)
      est  <- tmp$beta[1]
      tau2 <- tmp$tau2
      pval <- tmp$pval

      if (is.null(target)) {

         if (pval >= alpha) {

            fsnum <- 0

         } else {

            fsnum <- try(uniroot(.fsn.gen, interval=interval, extendInt="upX", tol=tol, maxiter=maxiter,
                                 yi=yi, vi=vi, vt=vt, est=est, tau2=tau2, tau2fix=tau2fix,
                                 test=test, weighted=weighted, target=target, alpha=alpha, exact=exact,
                                 method=method, mumiss=mumiss, upperint=max(interval), maxint=maxint, verbose=verbose)$root, silent=TRUE)

            if (inherits(fsnum, "try-error"))
               stop(mstyle$stop("Could not find fail-safe N based on a random-effects model (use verbose=TRUE for more info)."))

            if (fsnum > maxint)
               fsnum <- maxint

            tmp <- .fsn.gen(fsnum, yi=yi, vi=vi, vt=vt, est=est, tau2=tau2, tau2fix=tau2fix,
                            test=test, weighted=weighted, target=target, alpha=alpha, exact=exact,
                            method=method, mumiss=mumiss, upperint=max(interval), maxint=maxint, newest=TRUE)

         }

         target <- NA_real_

      } else {

         if (sign(target) != sign(est))
            target <- -1 * target

         if (identical(target, 0)) {
            fsnum <- Inf
         } else if (abs(target) >= abs(est)) {
            fsnum <- 0
         } else {

            fsnum <- try(uniroot(.fsn.gen, interval=interval, extendInt=ifelse(est > 0,"downX","upX"), tol=tol, maxiter=maxiter,
                                 yi=yi, vi=vi, vt=vt, est=est, tau2=tau2, tau2fix=tau2fix,
                                 test=test, weighted=weighted, target=target, alpha=alpha, exact=exact,
                                 method=method, mumiss=mumiss, upperint=max(interval), maxint=maxint, verbose=verbose)$root, silent=TRUE)

            if (inherits(fsnum, "try-error"))
               stop(mstyle$stop("Could not find fail-safe N based on a random-effects model (use verbose=TRUE for more info)."))

            if (fsnum > maxint)
               fsnum <- maxint

            tmp <- .fsn.gen(fsnum, yi=yi, vi=vi, vt=vt, est=est, tau2=tau2, tau2fix=tau2fix,
                            test=test, weighted=weighted, target=target, alpha=alpha, exact=exact,
                            method=method, mumiss=mumiss, upperint=max(interval), maxint=maxint, newest=TRUE)

         }

      }

      if (fsnum == 0) {
         est.fsn  <- est
         tau2.fsn <- tau2
         pval.fsn <- pval
      } else {
         est.fsn  <- tmp$est.fsn
         tau2.fsn <- tmp$tau2.fsn
         pval.fsn <- tmp$pval.fsn
      }

      if (fsnum >= maxint)
         ub.sign <- ">"

   }

   #########################################################################

   if (is.finite(fsnum) && abs(fsnum - round(fsnum)) >= .Machine$double.eps^0.5) {
      fsnum <- ceiling(fsnum)
   } else {
      fsnum <- round(fsnum)
   }

   res <- list(type=type, fsnum=fsnum, est=est, tau2=tau2, meanes=est, pval=pval, alpha=alpha, target=target,
               method=ifelse(type=="General", method, NA), est.fsn=est.fsn, tau2.fsn=tau2.fsn, pval.fsn=pval.fsn,
               ub.sign=ub.sign, digits=digits)

   class(res) <- "fsn"
   return(res)

}
