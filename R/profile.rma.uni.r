profile.rma.uni <- function(fitted,
   xlim, ylim, steps=20, lltol=1e-04, progbar=TRUE, parallel="no", ncpus=1, cl=NULL, plot=TRUE, pch=19, cline=FALSE, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(fitted, "rma.uni"))
      stop(mstyle$stop("Argument 'fitted' must be an object of class \"rma.uni\"."))

   if (inherits(fitted, "rma.ls"))
      stop(mstyle$stop("Method not available for objects of class \"rma.ls\"."))

   if (steps < 2)
      stop(mstyle$stop("Argument 'steps' must be >= 2."))

   x <- fitted

   parallel <- match.arg(parallel, c("no", "snow", "multicore"))

   if (parallel == "no" && ncpus > 1)
      parallel <- "snow"

   if (!is.null(cl) && inherits(cl, "SOCKcluster")) {
      parallel <- "snow"
      ncpus <- length(cl)
   }

   if (parallel == "snow" || parallel == "multicore") {

      if (!requireNamespace("parallel", quietly=TRUE))
         stop(mstyle$stop("Please install the 'parallel' package for parallel processing."))

      ncpus <- as.integer(ncpus)

      if (ncpus < 1L)
         stop(mstyle$stop("Argument 'ncpus' must be >= 1."))

   }

   if (progbar && !requireNamespace("pbapply", quietly=TRUE))
      stop(mstyle$stop("Please install the 'pbapply' package to obtain a progress bar."))

   if (!progbar) {
      pbo <- pbapply::pboptions(type = "none")
      on.exit(pbapply::pboptions(pbo), add = TRUE)
   }

   ddd <- list(...)

   if (.isTRUE(ddd$time))
      time.start <- proc.time()

   #########################################################################

   if (missing(xlim)) {

      ### if the user has not specified xlim, get CI for tau^2 (suppress warnings)

      vc.ci <- try(suppressWarnings(confint(x)), silent=TRUE)

      if (inherits(vc.ci, "try-error")) {

         vc.lb <- NA
         vc.ub <- NA

      } else {

         ### min() and max() so the actual value is within the xlim bounds
         ### could still get NAs for the bounds if the CI is the empty set

         vc.lb <- min(x$tau2, vc.ci$random[1,2])
         vc.ub <- max(.1, x$tau2, vc.ci$random[1,3]) ### if CI is equal to null set, then this still gives vc.ub = .1

      }

      if (is.na(vc.lb) || is.na(vc.ub)) {

         ### if the CI method fails, try a Wald-type CI for tau^2

         vc.lb <- max( 0, x$tau2 - 1.96 * x$se.tau2)
         vc.ub <- max(.1, x$tau2 + 1.96 * x$se.tau2)

      }

      if (is.na(vc.lb) || is.na(vc.ub)) {

         ### if this still results in NA bounds, use simple method

         vc.lb <- max( 0, x$tau2/4)
         vc.ub <- max(.1, x$tau2*4)

      }

      ### if all of that fails, throw an error

      if (is.na(vc.lb) || is.na(vc.ub))
         stop(mstyle$stop("Cannot set 'xlim' automatically. Please set this argument manually."))

      xlim <- c(vc.lb, vc.ub)

   } else {

      if (length(xlim) != 2L)
         stop(mstyle$stop("Argument 'xlim' should be a vector of length 2."))

      xlim <- sort(xlim)

   }

   vcs <- seq(xlim[1], xlim[2], length.out=steps)
   #return(vcs)

   if (length(vcs) <= 1L)
      stop(mstyle$stop("Cannot set 'xlim' automatically. Please set this argument manually."))

   if (parallel == "no")
      res <- pbapply::pblapply(seq_len(steps), .profile.rma.uni, obj=x, parallel=parallel, profile=TRUE, vcs=vcs)

   if (parallel == "multicore")
      res <- pbapply::pblapply(seq_len(steps), .profile.rma.uni, obj=x, parallel=parallel, profile=TRUE, vcs=vcs, cl=ncpus)
      #res <- parallel::mclapply(seq_len(steps), .profile.rma.uni, obj=x, parallel=parallel, profile=TRUE, vcs=vcs, mc.cores=ncpus)

   if (parallel == "snow") {
      if (is.null(cl)) {
         cl <- parallel::makePSOCKcluster(ncpus)
         on.exit(parallel::stopCluster(cl))
      }
      if (.isTRUE(ddd$LB)) {
         res <- parallel::parLapplyLB(cl, seq_len(steps), .profile.rma.uni, obj=x, parallel=parallel, profile=TRUE, vcs=vcs)
         #res <- parallel::clusterApplyLB(cl, seq_len(steps), .profile.rma.uni, obj=x, parallel=parallel, profile=TRUE, vcs=vcs)
         #res <- parallel::clusterMap(cl, .profile.rma.uni, seq_len(steps), MoreArgs=list(obj=x, parallel=parallel, profile=TRUE, vcs=vcs), .scheduling = "dynamic")
      } else {
         res <- pbapply::pblapply(seq_len(steps), .profile.rma.uni, obj=x, parallel=parallel, profile=TRUE, vcs=vcs, cl=cl)
         #res <- parallel::parLapply(cl, seq_len(steps), .profile.rma.uni, obj=x, parallel=parallel, profile=TRUE, vcs=vcs)
         #res <- parallel::clusterApply(cl, seq_len(steps), .profile.rma.uni, obj=x, parallel=parallel, profile=TRUE, vcs=vcs)
         #res <- parallel::clusterMap(cl, .profile.rma.uni, seq_len(steps), MoreArgs=list(obj=x, parallel=parallel, profile=TRUE, vcs=vcs))
      }
   }

   lls <- sapply(res, function(x) x$ll)
   beta  <- do.call("rbind", lapply(res, function(x) t(x$beta)))
   ci.lb <- do.call("rbind", lapply(res, function(x) t(x$ci.lb)))
   ci.ub <- do.call("rbind", lapply(res, function(x) t(x$ci.ub)))

   #########################################################################

   if (x$method %in% c("ML", "REML") && any(lls >= logLik(x) + lltol, na.rm=TRUE))
      warning(mstyle$warning("At least one profiled log-likelihood value is larger than the log-likelihood of the fitted model."), call.=FALSE)

   if (all(is.na(lls)))
      warning(mstyle$warning("All model fits failed. Cannot draw profile likelihood plot."), call.=FALSE)

   beta  <- data.frame(beta)
   ci.lb <- data.frame(ci.lb)
   ci.ub <- data.frame(ci.ub)
   names(beta)  <- rownames(x$beta)
   names(ci.lb) <- rownames(x$beta)
   names(ci.ub) <- rownames(x$beta)

   if (missing(ylim)) {

      if (any(!is.na(lls))) {
         ylim <- range(lls, na.rm=TRUE)
      } else {
         ylim <- rep(logLik(x), 2)
      }
      ylim[1] <- ylim[1] - .1
      ylim[2] <- ylim[2] + .1

   } else {

      if (length(ylim) != 2L)
         stop(mstyle$stop("Argument 'ylim' should be a vector of length 2."))

      ylim <- sort(ylim)

   }

   xlab <- expression(paste(tau^2, " Value"))
   title <- expression(paste("Profile Plot for ", tau^2))

   sav <- list(tau2=vcs, ll=lls, beta=beta, ci.lb=ci.lb, ci.ub=ci.ub, comps=1, xlim=xlim, ylim=ylim, method=x$method, vc=x$tau2, maxll=logLik(x), xlab=xlab, title=title)
   class(sav) <- "profile.rma"

   #########################################################################

   if (plot)
      plot(sav, pch=pch, cline=cline, ...)

   #########################################################################

   if (.isTRUE(ddd$time)) {
      time.end <- proc.time()
      .print.time(unname(time.end - time.start)[3])
   }

   invisible(sav)

}
