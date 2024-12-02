profile.rma.uni <- function(fitted,
   xlim, ylim, steps=20, lltol=1e-03, progbar=TRUE, parallel="no", ncpus=1, cl, plot=TRUE, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(fitted), must="rma.uni", notav=c("rma.ls", "rma.uni.selmodel", "rma.gen"))

   if (is.element(fitted$method, c("FE","EE","CE")))
      stop(mstyle$stop("Cannot profile tau^2 parameter for equal/fixed-effects models."))

   x <- fitted

   if (is.null(x$yi) || is.null(x$vi))
      stop(mstyle$stop("Information needed for profiling is not available in the model object."))

   if (anyNA(steps))
      stop(mstyle$stop("No missing values allowed in 'steps' argument."))

   if (length(steps) >= 2L) {
      if (missing(xlim))
         xlim <- range(steps)
      stepseq <- TRUE
   } else {
      if (steps < 2)
         stop(mstyle$stop("Argument 'steps' must be >= 2."))
      stepseq <- FALSE
   }

   parallel <- match.arg(parallel, c("no", "snow", "multicore"))

   if (parallel == "no" && ncpus > 1)
      parallel <- "snow"

   if (missing(cl))
      cl <- NULL

   if (!is.null(cl) && inherits(cl, "SOCKcluster")) {
      parallel <- "snow"
      ncpus <- length(cl)
   }

   if (parallel == "snow" && ncpus < 2)
      parallel <- "no"

   if (parallel == "snow" || parallel == "multicore") {

      if (!requireNamespace("parallel", quietly=TRUE))
         stop(mstyle$stop("Please install the 'parallel' package for parallel processing."))

      ncpus <- as.integer(ncpus)

      if (ncpus < 1L)
         stop(mstyle$stop("Argument 'ncpus' must be >= 1."))

   }

   if (!progbar) {
      pbo <- pbapply::pboptions(type="none")
      on.exit(pbapply::pboptions(pbo), add=TRUE)
   }

   ddd <- list(...)

   if (.isTRUE(ddd$time))
      time.start <- proc.time()

   pred <- isTRUE(ddd$pred)
   blup <- isTRUE(ddd$blup)

   newmods <- NULL

   if (pred) {

      if (!is.null(ddd$newmods))
         newmods <- ddd$newmods

      ### test if predict() works with the given newmods (and to get slab for [a])

      predres <- predict(x, newmods=newmods)

      if (length(predres$pred) == 0L)
         stop(mstyle$stop("Cannot compute predicted values."))

   }

   #########################################################################

   if (missing(xlim) || is.null(xlim)) {

      ### if the user has not specified xlim, set it automatically

      vc.ci <- try(suppressWarnings(confint(x)), silent=TRUE)

      if (inherits(vc.ci, "try-error")) {

         vc.lb <- NA_real_
         vc.ub <- NA_real_

      } else {

         ### min() and max() so the actual value is within the xlim bounds
         ### note: could still get NAs for the bounds if the CI is the empty set

         vc.lb <- min(x$tau2, vc.ci$random[1,2])
         vc.ub <- max(0.1, x$tau2, vc.ci$random[1,3]) # if CI is equal to null set, then this still gives vc.ub = 0.1

      }

      if (is.na(vc.lb) || is.na(vc.ub)) {

         ### if the CI method fails, try a Wald-type CI for tau^2

         vc.lb <- max(  0, x$tau2 - qnorm(0.995) * x$se.tau2)
         vc.ub <- max(0.1, x$tau2 + qnorm(0.995) * x$se.tau2)

      }

      if (is.na(vc.lb) || is.na(vc.ub)) {

         ### if this still results in NA bounds, use simple method

         vc.lb <- max(  0, x$tau2/4)
         vc.ub <- max(0.1, x$tau2*4)

      }

      ### if that fails, throw an error

      if (is.na(vc.lb) || is.na(vc.ub))
         stop(mstyle$stop("Cannot set 'xlim' automatically. Please set this argument manually."))

      xlim <- c(vc.lb, vc.ub)

      if (.isTRUE(ddd$sqrt))
         xlim <- sqrt(xlim)

   } else {

      if (length(xlim) != 2L)
         stop(mstyle$stop("Argument 'xlim' should be a vector of length 2."))

      xlim <- sort(xlim)

      ### note: if sqrt=TRUE, then xlim is assumed to be given in terms of tau

   }

   if (stepseq) {
      vcs <- steps
   } else {
      vcs <- seq(xlim[1], xlim[2], length.out=steps)
   }

   #return(vcs)

   if (length(vcs) <= 1L)
      stop(mstyle$stop("Cannot set 'xlim' automatically. Please set this argument manually."))

   if (!is.null(ddd[["code1"]]))
      eval(expr = parse(text = ddd[["code1"]]))

   ### if sqrt=TRUE, then the sequence of vcs are tau values, so square them for the actual profiling

   if (.isTRUE(ddd$sqrt))
      vcs <- vcs^2

   if (parallel == "no")
      res <- pbapply::pblapply(vcs, .profile.rma.uni, obj=x, parallel=parallel, profile=TRUE, pred=pred, blup=blup, newmods=newmods, code2=ddd$code2)

   if (parallel == "multicore")
      res <- pbapply::pblapply(vcs, .profile.rma.uni, obj=x, parallel=parallel, profile=TRUE, cl=ncpus, pred=pred, blup=blup, newmods=newmods, code2=ddd$code2)
      #res <- parallel::mclapply(vcs, .profile.rma.uni, obj=x, parallel=parallel, profile=TRUE, pred=pred, blup=blup, newmods=newmods, code2=ddd$code2, mc.cores=ncpus)

   if (parallel == "snow") {
      if (is.null(cl)) {
         cl <- parallel::makePSOCKcluster(ncpus)
         on.exit(parallel::stopCluster(cl), add=TRUE)
      }
      if (.isTRUE(ddd$LB)) {
         res <- parallel::parLapplyLB(cl, vcs, .profile.rma.uni, obj=x, parallel=parallel, profile=TRUE, pred=pred, blup=blup, newmods=newmods, code2=ddd$code2)
         #res <- parallel::clusterApplyLB(cl, vcs, .profile.rma.uni, obj=x, parallel=parallel, profile=TRUE, pred=pred, blup=blup, newmods=newmods, code2=ddd$code2)
         #res <- parallel::clusterMap(cl, .profile.rma.uni, vcs, MoreArgs=list(obj=x, parallel=parallel, profile=TRUE, pred=pred, blup=blup, newmods=newmods, code2=ddd$code2), .scheduling = "dynamic")
      } else {
         res <- pbapply::pblapply(vcs, .profile.rma.uni, obj=x, parallel=parallel, profile=TRUE, pred=pred, blup=blup, newmods=newmods, code2=ddd$code2, cl=cl)
         #res <- parallel::parLapply(cl, vcs, .profile.rma.uni, obj=x, parallel=parallel, profile=TRUE, pred=pred, blup=blup, newmods=newmods, code2=ddd$code2)
         #res <- parallel::clusterApply(cl, vcs, .profile.rma.uni, obj=x, parallel=parallel, profile=TRUE, pred=pred, blup=blup, newmods=newmods, code2=ddd$code2)
         #res <- parallel::clusterMap(cl, .profile.rma.uni, vcs, MoreArgs=list(obj=x, parallel=parallel, profile=TRUE, pred=pred, blup=blup, newmods=newmods, code2=ddd$code2))
      }
   }

   ### if sqrt=TRUE, then transform the tau^2 values back to tau values

   if (.isTRUE(ddd$sqrt)) {
      vcs <- sqrt(vcs)
      vc  <- sqrt(x$tau2)
   } else {
      vc <- x$tau2
   }

   lls   <- sapply(res, function(x) x$ll)
   beta  <- do.call(rbind, lapply(res, function(x) t(x$beta)))
   ci.lb <- do.call(rbind, lapply(res, function(x) t(x$ci.lb)))
   ci.ub <- do.call(rbind, lapply(res, function(x) t(x$ci.ub)))

   beta  <- data.frame(beta)
   ci.lb <- data.frame(ci.lb)
   ci.ub <- data.frame(ci.ub)
   names(beta)  <- rownames(x$beta)
   names(ci.lb) <- rownames(x$beta)
   names(ci.ub) <- rownames(x$beta)

   #########################################################################

   maxll <- c(logLik(x))

   if (x$method %in% c("ML","REML") && any(lls >= maxll + lltol, na.rm=TRUE))
      warning(mstyle$warning("At least one profiled log-likelihood value is larger than the log-likelihood of the fitted model."), call.=FALSE)

   if (all(is.na(lls)))
      warning(mstyle$warning("All model fits failed. Cannot draw profile likelihood plot."), call.=FALSE)

   if (.isTRUE(ddd$exp)) {
      lls <- exp(lls)
      maxll <- exp(maxll)
   }

   if (missing(ylim)) {

      if (any(is.finite(lls))) {
         if (xlim[1] <= vc && xlim[2] >= vc) {
            ylim <- range(c(maxll,lls[is.finite(lls)]), na.rm=TRUE)
         } else {
            ylim <- range(lls[is.finite(lls)], na.rm=TRUE)
         }
      } else {
         ylim <- rep(maxll, 2L)
      }

      if (!.isTRUE(ddd$exp))
         ylim <- ylim + c(-0.1, 0.1)

   } else {

      if (length(ylim) != 2L)
         stop(mstyle$stop("Argument 'ylim' should be a vector of length 2."))

      ylim <- sort(ylim)

   }

   if (.isTRUE(ddd$sqrt)) {
      xlab  <- expression(paste(tau, " Value"))
      title <- expression(paste("Profile Plot for ", tau))
   } else {
      xlab  <- expression(paste(tau^2, " Value"))
      title <- expression(paste("Profile Plot for ", tau^2))
   }

   sav <- list(tau2=vcs, ll=lls, beta=beta, ci.lb=ci.lb, ci.ub=ci.ub, comps=1, xlim=xlim, ylim=ylim, method=x$method, vc=vc, maxll=maxll, xlab=xlab, title=title, exp=ddd$exp, sqrt=ddd$sqrt)
   class(sav) <- "profile.rma"

   if (.isTRUE(ddd$sqrt))
      names(sav)[1] <- "tau"

   sav$I2 <- sapply(res, function(x) x$I2)
   sav$H2 <- sapply(res, function(x) x$H2)

   if (pred) {
      sav$pred       <- do.call(cbind, lapply(res, function(x) x$pred)) # use do.call(cbind, lapply()) instead of sapply() to always get a matrix, even when predicting a single value
      sav$pred.ci.lb <- do.call(cbind, lapply(res, function(x) x$pred.ci.lb))
      sav$pred.ci.ub <- do.call(cbind, lapply(res, function(x) x$pred.ci.ub))
      sav$pred.pi.lb <- do.call(cbind, lapply(res, function(x) x$pred.pi.lb))
      sav$pred.pi.ub <- do.call(cbind, lapply(res, function(x) x$pred.pi.ub))
      rownames(sav$pred) <- rownames(sav$pred.ci.lb) <- rownames(sav$pred.ci.ub) <- rownames(sav$pred.pi.lb) <- rownames(sav$pred.pi.ub) <- predres$slab # [a]
   }

   if (blup) {
      sav$blup       <- sapply(res, function(x) x$blup)
      sav$blup.se    <- sapply(res, function(x) x$blup.se)
      sav$blup.pi.lb <- sapply(res, function(x) x$blup.pi.lb)
      sav$blup.pi.ub <- sapply(res, function(x) x$blup.pi.ub)
      rownames(sav$blup) <- x$slab[x$not.na]
      rownames(sav$blup.se) <- x$slab[x$not.na]
      rownames(sav$blup.pi.lb) <- x$slab[x$not.na]
      rownames(sav$blup.pi.ub) <- x$slab[x$not.na]
   }

   #########################################################################

   if (plot)
      plot(sav, ...)

   #########################################################################

   if (.isTRUE(ddd$time)) {
      time.end <- proc.time()
      .print.time(unname(time.end - time.start)[3])
   }

   invisible(sav)

}
