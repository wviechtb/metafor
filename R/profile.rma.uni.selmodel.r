profile.rma.uni.selmodel <- function(fitted, tau2, delta,
   xlim, ylim, steps=20, lltol=1e-03, progbar=TRUE, parallel="no", ncpus=1, cl, plot=TRUE, pch=19, refline=TRUE, cline=FALSE, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(fitted), must="rma.uni.selmodel")

   if (steps < 2)
      stop(mstyle$stop("Argument 'steps' must be >= 2."))

   x <- fitted

   if (x$betaspec) ### TODO: consider allowing profiling over beta values as well
      stop(mstyle$stop("Cannot profile when one or more beta values were fixed."))

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

   #########################################################################

   ### check if user has not specified tau2 or delta argument

   if (missing(tau2) && missing(delta)) {

      mc <- match.call()

      ### total number of non-fixed components

      comps <- ifelse(!is.element(x$method, c("FE","EE","CE")) && !x$tau2.fix, 1, 0) + sum(!x$delta.fix)

      if (comps == 0)
         stop(mstyle$stop("No components in the model for which a profile likelihood can be constructed."))

      if (plot) {
         if (dev.cur() == 1) {
            par(mfrow=n2mfrow(comps))
            #on.exit(par(mfrow=c(1,1)), add=TRUE)
         }
      }

      sav <- list()
      j <- 0

      if (!is.element(x$method, c("FE","EE","CE")) && !x$tau2.fix) {
         j <- j + 1
         mc.vc <- mc
         mc.vc$tau2 <- 1
         mc.vc$time <- FALSE
         #mc.vc$fitted <- quote(x)
         if (progbar)
            cat(mstyle$verbose(paste("Profiling tau2\n")))
         sav[[j]] <- eval(mc.vc, envir=parent.frame())
      }

      if (any(!x$delta.fix)) {
         for (pos in seq_len(x$deltas)[!x$delta.fix]) {
            j <- j + 1
            mc.vc <- mc
            mc.vc$delta <- pos
            mc.vc$time <- FALSE
            #mc.vc$fitted <- quote(x)
            if (progbar)
               cat(mstyle$verbose(paste("Profiling delta =", pos, "\n")))
            sav[[j]] <- eval(mc.vc, envir=parent.frame())
         }
      }

      ### if there is just one component, turn the list of lists into a simple list

      if (comps == 1)
         sav <- sav[[1]]

      sav$comps <- comps

      if (.isTRUE(ddd$time)) {
         time.end <- proc.time()
         .print.time(unname(time.end - time.start)[3])
      }

      class(sav) <- "profile.rma"
      return(invisible(sav))

   }

   ### check if user has specified more than one of these arguments

   if (sum(!missing(tau2), !missing(delta)) > 1L)
      stop(mstyle$stop("Must specify only one of the 'tau2' or 'delta' arguments."))

   ### check if model actually contains (at least one) such a component and that it was actually estimated

   if (!missing(tau2) && (is.element(x$method, c("FE","EE","CE")) || x$tau2.fix))
      stop(mstyle$stop("Model does not contain an (estimated) 'tau2' component."))

   if (!missing(delta) && all(x$delta.fix))
      stop(mstyle$stop("Model does not contain any estimated 'delta' components."))

   ### check if user specified more than one tau2 or delta component

   if (!missing(tau2) && (length(tau2) > 1L))
      stop(mstyle$stop("Can only specify one 'tau2' component."))

   if (!missing(delta) && (length(delta) > 1L))
      stop(mstyle$stop("Can only specify one 'delta' component."))

   ### check if user specified a logical

   if (!missing(tau2) && is.logical(tau2) && isTRUE(tau2))
      tau2 <- 1

   if (!missing(delta) && is.logical(delta))
      stop(mstyle$stop("Must specify a number for the 'delta' component."))

   ### check if user specified a component that does not exist

   if (!missing(tau2) && (tau2 > 1 || tau2 <= 0))
      stop(mstyle$stop("No such 'tau2' component in the model."))

   if (!missing(delta) && (delta > x$deltas || delta <= 0))
      stop(mstyle$stop("No such 'delta' component in the model."))

   ### check if user specified a component that was fixed

   if (!missing(tau2) && x$tau2.fix)
      stop(mstyle$stop("Specified 'tau2' component was fixed."))

   if (!missing(delta) && x$delta.fix[delta])
      stop(mstyle$stop("Specified 'delta' component was fixed."))

   ### if everything is good so far, get value of the variance component and set 'comp'

   delta.pos <- NA_integer_

   if (!missing(tau2)) {
      vc <- x$tau2
      comp <- "tau2"
      tau2.pos <- 1
   }

   if (!missing(delta)) {
      vc <- x$delta[delta]
      comp <- "delta"
      delta.pos <- delta
   }

   #return(list(comp=comp, vc=vc))

   if (missing(xlim)) {

      ### if the user has not specified xlim, set it automatically

      if (comp == "tau2") {
         if (is.na(x$se.tau2)) {
            vc.lb <- max(0, vc/4)
            vc.ub <- min(max(.1, vc*4), x$tau2.max)
         } else {
            vc.lb <- max(0, vc - qnorm(.995) * x$se.tau2)
            vc.ub <- min(max(.1, vc + qnorm(.995) * x$se.tau2), x$tau2.max)
         }
      }
      if (comp == "delta") {
         if (is.na(x$se.delta[delta])) {
            vc.lb <- max(0, vc/4, x$delta.min[delta])
            vc.ub <- min(max(.1, vc*4), x$delta.max[delta])
         } else {
            vc.lb <- max(0, vc - qnorm(.995) * x$se.delta[delta], x$delta.min[delta])
            vc.ub <- min(max(.1, vc + qnorm(.995) * x$se.delta[delta]), x$delta.max[delta])
         }
      }

      ### if that fails, throw an error

      if (is.na(vc.lb) || is.na(vc.ub))
         stop(mstyle$stop("Cannot set 'xlim' automatically. Please set this argument manually."))

      xlim <- c(vc.lb, vc.ub)

   } else {

      if (length(xlim) != 2L)
         stop(mstyle$stop("Argument 'xlim' should be a vector of length 2."))

      xlim <- sort(xlim)

      if (comp == "tau2") {
         if (xlim[1] < 0)
            stop(mstyle$stop("Lower bound for profiling must be >= 0."))
      }
      if (comp == "delta") {
         if (xlim[1] < x$delta.min[delta])
            stop(mstyle$stop(paste0("Lower bound for profiling must be >= ", x$delta.min[delta], ".")))
         if (xlim[2] > x$delta.max[delta])
            stop(mstyle$stop(paste0("Upper bound for profiling must be <= ", x$delta.max[delta], ".")))
      }

   }

   vcs <- seq(xlim[1], xlim[2], length.out=steps)
   #return(vcs)

   if (length(vcs) <= 1L)
      stop(mstyle$stop("Cannot set 'xlim' automatically. Please set this argument manually."))

   if (parallel == "no")
      res <- pbapply::pblapply(vcs, .profile.rma.uni.selmodel, obj=x, comp=comp, delta.pos=delta.pos, parallel=parallel, profile=TRUE)

   if (parallel == "multicore")
      res <- pbapply::pblapply(vcs, .profile.rma.uni.selmodel, obj=x, comp=comp, delta.pos=delta.pos, parallel=parallel, profile=TRUE, cl=ncpus)
      #res <- parallel::mclapply(vcs, .profile.rma.uni.selmodel, obj=x, comp=comp, delta.pos=delta.pos, parallel=parallel, profile=TRUE, mc.cores=ncpus)

   if (parallel == "snow") {
      if (is.null(cl)) {
         cl <- parallel::makePSOCKcluster(ncpus)
         on.exit(parallel::stopCluster(cl), add=TRUE)
      }
      if (.isTRUE(ddd$LB)) {
         res <- parallel::parLapplyLB(cl, vcs, .profile.rma.uni.selmodel, obj=x, comp=comp, delta.pos=delta.pos, parallel=parallel, profile=TRUE)
         #res <- parallel::clusterApplyLB(cl, vcs, .profile.rma.uni.selmodel, obj=x, comp=comp, delta.pos=delta.pos, parallel=parallel, profile=TRUE)
      } else {
         res <- pbapply::pblapply(vcs, .profile.rma.uni.selmodel, obj=x, comp=comp, delta.pos=delta.pos, parallel=parallel, profile=TRUE, cl=cl)
         #res <- parallel::parLapply(cl, vcs, .profile.rma.uni.selmodel, obj=x, comp=comp, delta.pos=delta.pos, parallel=parallel, profile=TRUE)
         #res <- parallel::clusterApply(cl, vcs, .profile.rma.uni.selmodel, obj=x, comp=comp, delta.pos=delta.pos, parallel=parallel, profile=TRUE)
      }
   }

   lls <- sapply(res, function(x) x$ll)
   beta  <- do.call(rbind, lapply(res, function(x) t(x$beta)))
   ci.lb <- do.call(rbind, lapply(res, function(x) t(x$ci.lb)))
   ci.ub <- do.call(rbind, lapply(res, function(x) t(x$ci.ub)))

   #########################################################################

   if (any(lls >= logLik(x) + lltol, na.rm=TRUE))
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

      if (any(is.finite(lls))) {
         if (xlim[1] <= vc && xlim[2] >= vc) {
            ylim <- range(c(logLik(x),lls[is.finite(lls)]), na.rm=TRUE)
         } else {
            ylim <- range(lls[is.finite(lls)])
         }
      } else {
         ylim <- rep(logLik(x), 2L)
      }
      ylim <- ylim + c(-0.1, 0.1)

   } else {

      if (length(ylim) != 2L)
         stop(mstyle$stop("Argument 'ylim' should be a vector of length 2."))

      ylim <- sort(ylim)

   }

   if (comp == "tau2") {
      xlab <- expression(paste(tau^2, " Value"))
      title <- expression(paste("Profile Plot for ", tau^2))
   }
   if (comp == "delta") {
      if (x$deltas == 1L) {
         xlab <- expression(paste(delta, " Value"))
         title <- expression(paste("Profile Plot for ", delta))
      } else {
         xlab <- bquote(delta[.(delta)] ~ "Value")
         title <- bquote("Profile Plot for" ~ delta[.(delta)])
      }
   }

   sav <- list(vc=vcs, ll=lls, beta=beta, ci.lb=ci.lb, ci.ub=ci.ub, comps=1, ylim=ylim, method=x$method, vc=vc, maxll=logLik(x), xlab=xlab, title=title)
   names(sav)[1] <- switch(comp, tau2="tau2", delta="delta")
   class(sav) <- "profile.rma"

   #########################################################################

   if (plot)
      plot(sav, pch=pch, refline=refline, cline=cline, ...)

   #########################################################################

   if (.isTRUE(ddd$time)) {
      time.end <- proc.time()
      .print.time(unname(time.end - time.start)[3])
   }

   invisible(sav)

}
