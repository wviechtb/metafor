tes.default <- function(x, vi, sei, subset,
   H0=0, alternative="two.sided", alpha=.05, theta, tau2,
   test, tes.alternative="greater", progbar=TRUE, tes.alpha=.10,
   digits, ...) {

   # allow multiple alpha values? plot for pval as a function of alpha?

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   alternative <- match.arg(alternative, c("two.sided", "greater", "less"))
   tes.alternative <- match.arg(tes.alternative, c("two.sided", "greater", "less"))

   if (missing(subset))
      subset <- NULL

   if (!is.null(subset))
      subset <- .setnafalse(subset)

   if (alpha <= 0 || alpha >= 1)
      stop(mstyle$stop("Value of 'alpha' needs to be > 0 and < 1."))

   if (tes.alpha <= 0 || tes.alpha >= 1)
      stop(mstyle$stop("Value of 'tes.alpha' needs to be > 0 and < 1."))

   if (alternative == "two.sided")
      crit <- qnorm(alpha/2, lower.tail=FALSE)
   if (alternative == "greater")
      crit <- qnorm(alpha, lower.tail=FALSE)
   if (alternative == "less")
      crit <- qnorm(alpha, lower.tail=TRUE)

   ### set defaults for digits

   if (missing(digits)) {
      digits <- .set.digits(dmiss=TRUE)
   } else {
      digits <- .set.digits(digits, dmiss=FALSE)
   }

   ddd <- list(...)

   .chkdots(ddd, c("correct", "rel.tol", "subdivisions", "tau2.lb", "find.lim"))

   if (!is.null(ddd$correct)) {
      correct <- ddd$correct
   } else {
      correct <- FALSE
   }

   if (!is.null(ddd$rel.tol)) {
      rel.tol <- ddd$rel.tol
   } else {
      rel.tol <- .Machine$double.eps^0.25
   }

   if (!is.null(ddd$subdivisions)) {
      subdivisions <- ddd$subdivisions
   } else {
      subdivisions <- 100L
   }

   if (!is.null(ddd$tau2.lb)) {
      tau2.lb <- ddd$tau2.lb
   } else {
      #tau2.lb <- 0.0001
      tau2.lb <- 0
   }

   if (!is.null(ddd$find.lim)) {
      find.lim <- ddd$find.lim
   } else {
      find.lim <- TRUE
   }

   #########################################################################

   ### check if sampling variances and/or standard errors are available

   if (missing(vi))
      vi <- NULL

   if (missing(sei))
      sei <- NULL

   if (is.null(vi)) {
      if (!is.null(sei))
         vi <- sei^2
   }

   if (is.null(vi))
      stop(mstyle$stop("Must specify 'vi' or 'sei' argument."))

   yi <- x

   ### check length of yi and vi

   if (length(yi) != length(vi))
      stop(mstyle$stop("Length of 'yi' and 'vi' (or 'sei') is not the same."))

   k.f <- length(yi)

   ### checks on H0

   if (length(H0) != 1L)
      stop(mstyle$stop("Argument 'H0' must specify a single value."))

   ### checks on theta

   if (missing(theta) || is.null(theta)) {

      single.theta <- TRUE
      est.theta <- TRUE
      theta <- rep(0, k.f)

   } else {

      if (length(theta) == 1L) {
         single.theta <- TRUE
         est.theta <- FALSE
         theta.1 <- theta
         theta <- rep(theta, k.f)
      } else {
         single.theta <- FALSE
         est.theta <- FALSE
      }

      if (length(theta) != k.f)
         stop(mstyle$stop("Length of 'theta' and 'yi' is not the same."))

   }

   #########################################################################

   ### if a subset of studies is specified

   if (!is.null(subset)) {
      yi <- yi[subset]
      vi <- vi[subset]
      theta <- theta[subset]
   }

   ### check for NAs and act accordingly

   has.na <- is.na(yi) | is.na(vi) | is.na(theta)

   if (any(has.na)) {

      not.na <- !has.na

      if (na.act == "na.omit" || na.act == "na.exclude" || na.act == "na.pass") {

         yi <- yi[not.na]
         vi <- vi[not.na]
         theta <- theta[not.na]
         warning(mstyle$warning("Studies with NAs omitted from test."), call.=FALSE)

      }

      if (na.act == "na.fail")
         stop(mstyle$stop("Missing values in data."))

   }

   #########################################################################

   k <- length(yi)

   if (k == 0L)
      stop(mstyle$stop("Stopped because k = 0."))

   sei <- sqrt(vi)
   zi  <- (yi - H0) / sei
   wi  <- 1 / vi

   if (est.theta) {
      theta.1 <- .wmean(yi, wi)
      theta <- rep(theta.1, k)
   }

   if (missing(tau2) || is.null(tau2) || tau2 <= tau2.lb) {

      if (alternative == "two.sided")
         pow <- pnorm(crit, mean=(theta-H0)/sei, 1, lower.tail=FALSE) + pnorm(-crit, mean=(theta-H0)/sei, 1, lower.tail=TRUE)
      if (alternative == "greater")
         pow <- pnorm(crit, mean=(theta-H0)/sei, 1, lower.tail=FALSE)
      if (alternative == "less")
         pow <- pnorm(crit, mean=(theta-H0)/sei, 1, lower.tail=TRUE)

   } else {

      tau <- sqrt(tau2)

      pow <- rep(NA_real_, k)

      for (i in seq_len(k)) {

         res <- try(integrate(.tes.intfun, lower=theta[i]-5*tau, upper=theta[i]+5*tau, theta=theta[i], tau=tau, sei=sei[i], H0=H0, alternative=alternative, crit=crit,
                              rel.tol=rel.tol, subdivisions=subdivisions, stop.on.error=FALSE), silent=TRUE)

         if (inherits(res, "try-error")) {
            stop(mstyle$stop(paste0("Could not integrate over density in study ", i, ".")))
         } else {
            pow[i] <- res$value
         }

      }

   }

   if (alternative == "two.sided")
      sig <- abs(zi) >= crit
   if (alternative == "greater")
      sig <- zi >= crit
   if (alternative == "less")
      sig <- zi <= crit

   E <- sum(pow)
   O <- sum(sig)

   if (tes.alternative == "two.sided")
      js <- 0:k
   if (tes.alternative == "greater")
      js <- O:k
   if (tes.alternative == "less")
      js <- 0:O

   if (missing(test) || is.null(test)) {
      tot <- sum(sapply(js, function(j) choose(k,j)))
      if (tot <= 10^6) {
         test <- "exact"
      } else {
         test <- "chi2"
      }
   } else {
      test <- match.arg(test, c("chi2", "binom", "exact"))
   }

   ### set defaults for progbar

   if (missing(progbar))
      progbar <- ifelse(test == "exact", TRUE, FALSE)

   if (test == "chi2") {
      res <- suppressWarnings(prop.test(O, k, p=E/k, alternative=tes.alternative, correct=correct))
      X2 <- unname(res$statistic)
      pval <- res$p.value
   }

   if (test == "binom") {
      res <- binom.test(O, k, p=E/k, alternative=tes.alternative)
      X2 <- NA
      pval <- binom.test(O, k, p=E/k, alternative=tes.alternative)$p.value
   }

   if (test == "exact") {

      X2 <- NA

      if (progbar)
         pbar <- pbapply::startpb(min=0, max=length(js))

      prj <- rep(NA_real_, length(js))
      id <- seq_len(k)

      for (j in 1:length(js)) {

         if (progbar)
            pbapply::setpb(pbar, j)

         if (js[j] == 0L) {
            prj[j] <- prod(1-pow)
         } else if (js[j] == k) {
            prj[j] <- prod(pow)
         } else {
            tmp <- try(suppressWarnings(sum(combn(k, js[j], FUN = function(i) {
               sel <- i
               not <- id[-i]
               prod(pow[sel])*prod(1-pow[not])
            }))), silent=TRUE)
            if (inherits(tmp, "try-error")) {
               if (progbar)
                  pbapply::closepb(pbar)
               stop(mstyle$stop(paste0("Number of combinations too large to do an exact test (use test=\"chi2\" or test=\"binomial\" instead).")))
            } else {
               prj[j] <- tmp
            }
         }

      }

      if (progbar)
         pbapply::closepb(pbar)

      if (tes.alternative == "two.sided")
         pval <- sum(prj[prj <= prj[O+1] + .Machine$double.eps^0.5])
      if (tes.alternative == "greater")
         pval <- sum(prj)
      if (tes.alternative == "less")
         pval <- sum(prj)

      pval[pval > 1] <- 1

   }

   theta.lim <- NULL

   if (find.lim && single.theta) {

      if (tes.alternative == "greater") {

         diff.H0 <- .tes.lim(H0, yi=yi, vi=vi, H0=H0, alternative=alternative, alpha=alpha, tau2=tau2, test=test, tes.alternative=tes.alternative, progbar=FALSE, tes.alpha=tes.alpha, correct=correct, rel.tol=rel.tol, subdivisions=subdivisions, tau2.lb=tau2.lb)
         if (diff.H0 >= 0) {
            theta.lim <- NA
         } else {
            if (theta.1 >= H0) {
               theta.lim <- try(uniroot(.tes.lim, interval=c(H0,theta.1), extendInt="upX", yi=yi, vi=vi, H0=H0, alternative=alternative, alpha=alpha, tau2=tau2, test=test, tes.alternative=tes.alternative, progbar=FALSE, tes.alpha=tes.alpha, correct=correct, rel.tol=rel.tol, subdivisions=subdivisions, tau2.lb=tau2.lb)$root, silent=TRUE)
            } else {
               theta.lim <- try(uniroot(.tes.lim, interval=c(theta.1,H0), extendInt="downX", yi=yi, vi=vi, H0=H0, alternative=alternative, alpha=alpha, tau2=tau2, test=test, tes.alternative=tes.alternative, progbar=FALSE, tes.alpha=tes.alpha, correct=correct, rel.tol=rel.tol, subdivisions=subdivisions, tau2.lb=tau2.lb)$root, silent=TRUE)
            }
            if (inherits(theta.lim, "try-error"))
               theta.lim <- NA
         }

      }

      if (tes.alternative == "less") {

         diff.H0 <- .tes.lim(H0, yi=yi, vi=vi, H0=H0, alternative=alternative, alpha=alpha, tau2=tau2, test=test, tes.alternative=tes.alternative, progbar=FALSE, tes.alpha=tes.alpha, correct=correct, rel.tol=rel.tol, subdivisions=subdivisions, tau2.lb=tau2.lb)
         if (diff.H0 <= 0) {
               theta.lim <- NA
            } else {
               if (theta.1 >= H0) {
                  theta.lim <- try(uniroot(.tes.lim, interval=c(H0,theta.1), extendInt="downX", yi=yi, vi=vi, H0=H0, alternative=alternative, alpha=alpha, tau2=tau2, test=test, tes.alternative=tes.alternative, progbar=FALSE, tes.alpha=tes.alpha, correct=correct, rel.tol=rel.tol, subdivisions=subdivisions, tau2.lb=tau2.lb)$root, silent=TRUE)
               } else {
                  theta.lim <- try(uniroot(.tes.lim, interval=c(theta.1,H0), extendInt="upX", yi=yi, vi=vi, H0=H0, alternative=alternative, alpha=alpha, tau2=tau2, test=test, tes.alternative=tes.alternative, progbar=FALSE, tes.alpha=tes.alpha, correct=correct, rel.tol=rel.tol, subdivisions=subdivisions, tau2.lb=tau2.lb)$root, silent=TRUE)
               }
               if (inherits(theta.lim, "try-error"))
                  theta.lim <- NA
         }

      }

      if (tes.alternative == "two.sided") {

         theta.lim.lb <- tes.default(x=yi, vi=vi, H0=H0, alternative=alternative, alpha=alpha, theta=theta.1, tau2=tau2, test=test, tes.alternative="greater", progbar=FALSE, tes.alpha=tes.alpha/2, correct=correct, rel.tol=rel.tol, subdivisions=subdivisions, tau2.lb=tau2.lb, find.lim=TRUE)$theta.lim
         theta.lim.ub <- tes.default(x=yi, vi=vi, H0=H0, alternative=alternative, alpha=alpha, theta=theta.1, tau2=tau2, test=test, tes.alternative="less", progbar=FALSE, tes.alpha=tes.alpha/2, correct=correct, rel.tol=rel.tol, subdivisions=subdivisions, tau2.lb=tau2.lb, find.lim=TRUE)$theta.lim
         theta.lim <- c(theta.lim.lb, theta.lim.ub)

      }

   }

   if (single.theta)
      theta <- theta.1

   res <- list(k=k, O=O, E=E, OEratio=O/E,
               test=test, X2=X2, pval=pval,
               power=pow, sig=sig, theta=theta,
               theta.lim=theta.lim, tes.alternative=tes.alternative, tes.alpha=tes.alpha,
               digits=digits)

   class(res) <- "tes"
   return(res)

}
