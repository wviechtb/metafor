profile.rma.mv <- function(fitted, sigma2, tau2, rho, gamma2, phi,
   xlim, ylim, steps=20, lltol=1e-03, progbar=TRUE, parallel="no", ncpus=1, cl, plot=TRUE, pch=19, refline=TRUE, cline=FALSE, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(fitted), must="rma.mv")

   if (steps < 2)
      stop(mstyle$stop("Argument 'steps' must be >= 2."))

   x <- fitted

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
      on.exit(pbapply::pboptions(pbo))
   }

   ddd <- list(...)

   if (.isTRUE(ddd$time))
      time.start <- proc.time()

   if (!is.null(ddd$startmethod))
      warning(mstyle$warning("Argument 'startmethod' has been deprecated."))

   #########################################################################

   ### check if user has not specified one of the sigma2, tau2, rho, gamma2, or phi arguments

   if (missing(sigma2) && missing(tau2) && missing(rho) && missing(gamma2) && missing(phi)) {

      mc <- match.call()

      ### total number of non-fixed components

      comps <- ifelse(x$withS, sum(!x$vc.fix$sigma2), 0) + ifelse(x$withG, sum(!x$vc.fix$tau2) + sum(!x$vc.fix$rho), 0) + ifelse(x$withH, sum(!x$vc.fix$gamma2) + sum(!x$vc.fix$phi), 0)

      if (comps == 0)
         stop(mstyle$stop("No components in the model for which a profile likelihood can be constructed."))

      if (plot) {
         if (dev.cur() == 1) {
            par(mfrow=c(comps, 1))
            #on.exit(par(mfrow=c(1,1)), add=TRUE)
         }
      }

      sav <- list()
      j <- 0

      if (x$withS && any(!x$vc.fix$sigma2)) {
         for (pos in seq_len(x$sigma2s)[!x$vc.fix$sigma2]) {
            j <- j + 1
            mc.vc <- mc
            mc.vc$sigma2 <- pos
            mc.vc$time <- FALSE
            #mc.vc$fitted <- quote(x)
            if (progbar)
               cat(mstyle$verbose(paste("Profiling sigma2 =", pos, "\n")))
            sav[[j]] <- eval(mc.vc, envir=parent.frame())
         }
      }

      if (x$withG) {
         if (any(!x$vc.fix$tau2)) {
            for (pos in seq_len(x$tau2s)[!x$vc.fix$tau2]) {
               j <- j + 1
               mc.vc <- mc
               mc.vc$tau2 <- pos
               mc.vc$time <- FALSE
               #mc.vc$fitted <- quote(x)
               if (progbar)
                  cat(mstyle$verbose(paste("Profiling tau2 =", pos, "\n")))
               sav[[j]] <- eval(mc.vc, envir=parent.frame())
            }
         }
         if (any(!x$vc.fix$rho)) {
            for (pos in seq_len(x$rhos)[!x$vc.fix$rho]) {
               j <- j + 1
               mc.vc <- mc
               mc.vc$rho <- pos
               mc.vc$time <- FALSE
               #mc.vc$fitted <- quote(x)
               if (progbar)
                  cat(mstyle$verbose(paste("Profiling rho =", pos, "\n")))
               sav[[j]] <- eval(mc.vc, envir=parent.frame())
            }
         }
      }

      if (x$withH) {
         if (any(!x$vc.fix$gamma2)) {
            for (pos in seq_len(x$gamma2s)[!x$vc.fix$gamma2]) {
               j <- j + 1
               mc.vc <- mc
               mc.vc$gamma2 <- pos
               mc.vc$time <- FALSE
               #mc.vc$fitted <- quote(x)
               if (progbar)
                  cat(mstyle$verbose(paste("Profiling gamma2 =", pos, "\n")))
               sav[[j]] <- eval(mc.vc, envir=parent.frame())
            }
         }
         if (any(!x$vc.fix$phi)) {
            for (pos in seq_len(x$phis)[!x$vc.fix$phi]) {
               j <- j + 1
               mc.vc <- mc
               mc.vc$phi <- pos
               mc.vc$time <- FALSE
               #mc.vc$fitted <- quote(x)
               if (progbar)
                  cat(mstyle$verbose(paste("Profiling phi =", pos, "\n")))
               sav[[j]] <- eval(mc.vc, envir=parent.frame())
            }
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

   #if (missing(sigma2) && missing(tau2) && missing(rho) && missing(gamma2) && missing(phi))
   #   stop(mstyle$stop("Must specify one of the arguments 'sigma2', 'tau2', 'rho', 'gamma2', or 'phi'."))

   ### check if user has specified more than one of these arguments

   if (sum(!missing(sigma2), !missing(tau2), !missing(rho), !missing(gamma2), !missing(phi)) > 1L)
      stop(mstyle$stop("Must specify only one of the arguments 'sigma2', 'tau2', 'rho', 'gamma2', or 'phi'."))

   ### check if model actually contains (at least one) such a component and that it was actually estimated
   ### note: a component that is not in the model is NA; components that are fixed are TRUE

   if (!missing(sigma2) && (all(is.na(x$vc.fix$sigma2)) || all(x$vc.fix$sigma2)))
      stop(mstyle$stop("Model does not contain any (estimated) 'sigma2' components."))

   if (!missing(tau2) && (all(is.na(x$vc.fix$tau2)) || all(x$vc.fix$tau2)))
      stop(mstyle$stop("Model does not contain any (estimated) 'tau2' components."))

   if (!missing(rho) && c(all(is.na(x$vc.fix$rho)) || all(x$vc.fix$rho)))
      stop(mstyle$stop("Model does not contain any (estimated) 'rho' components."))

   if (!missing(gamma2) && (all(is.na(x$vc.fix$gamma2)) || all(x$vc.fix$gamma2)))
      stop(mstyle$stop("Model does not contain any (estimated) 'gamma2' components."))

   if (!missing(phi) && c(all(is.na(x$vc.fix$phi)) || all(x$vc.fix$phi)))
      stop(mstyle$stop("Model does not contain any (estimated) 'phi' components."))

   ### check if user specified more than one sigma2, tau2, rho, gamma2, or rho component

   if (!missing(sigma2) && (length(sigma2) > 1L))
      stop(mstyle$stop("Can only specify one 'sigma2' component."))

   if (!missing(tau2) && (length(tau2) > 1L))
      stop(mstyle$stop("Can only specify one 'tau2' component."))

   if (!missing(rho) && (length(rho) > 1L))
      stop(mstyle$stop("Can only specify one 'rho' component."))

   if (!missing(gamma2) && (length(gamma2) > 1L))
      stop(mstyle$stop("Can only specify one 'gamma2' component."))

   if (!missing(phi) && (length(phi) > 1L))
      stop(mstyle$stop("Can only specify one 'phi' component."))

   ### check if user specified a logical

   if (!missing(sigma2) && is.logical(sigma2))
      stop(mstyle$stop("Must specify the number for the 'sigma2' component."))

   if (!missing(tau2) && is.logical(tau2))
      stop(mstyle$stop("Must specify the number for the 'tau2' component."))

   if (!missing(rho) && is.logical(rho))
      stop(mstyle$stop("Must specify the number for the 'rho' component."))

   if (!missing(gamma2) && is.logical(gamma2))
      stop(mstyle$stop("Must specify the number for the 'gamma2' component."))

   if (!missing(phi) && is.logical(phi))
      stop(mstyle$stop("Must specify the number for the 'phi' component."))

   ### check if user specified a component that does not exist

   if (!missing(sigma2) && (sigma2 > length(x$vc.fix$sigma2) || sigma2 <= 0))
      stop(mstyle$stop("No such 'sigma2' component in the model."))

   if (!missing(tau2) && (tau2 > length(x$vc.fix$tau2) || tau2 <= 0))
      stop(mstyle$stop("No such 'tau2' component in the model."))

   if (!missing(rho) && (rho > length(x$vc.fix$rho) || rho <= 0))
      stop(mstyle$stop("No such 'rho' component in the model."))

   if (!missing(gamma2) && (gamma2 > length(x$vc.fix$gamma2) || gamma2 <= 0))
      stop(mstyle$stop("No such 'gamma2' component in the model."))

   if (!missing(phi) && (phi > length(x$vc.fix$phi) || phi <= 0))
      stop(mstyle$stop("No such 'phi' component in the model."))

   ### check if user specified a component that was fixed

   if (!missing(sigma2) && x$vc.fix$sigma2[sigma2])
      stop(mstyle$stop("Specified 'sigma2' component was fixed."))

   if (!missing(tau2) && x$vc.fix$tau2[tau2])
      stop(mstyle$stop("Specified 'tau2' component was fixed."))

   if (!missing(rho) && x$vc.fix$rho[rho])
      stop(mstyle$stop("Specified 'rho' component was fixed."))

   if (!missing(gamma2) && x$vc.fix$gamma2[gamma2])
      stop(mstyle$stop("Specified 'gamma2' component was fixed."))

   if (!missing(phi) && x$vc.fix$phi[phi])
      stop(mstyle$stop("Specified 'phi' component was fixed."))

   ### if everything is good so far, get value of the variance component and set 'comp'

   sigma2.pos <- NA
   tau2.pos   <- NA
   rho.pos    <- NA
   gamma2.pos <- NA
   phi.pos    <- NA

   if (!missing(sigma2)) {
      vc <- x$sigma2[sigma2]
      comp <- "sigma2"
      sigma2.pos <- sigma2
   }

   if (!missing(tau2)) {
      vc <- x$tau2[tau2]
      comp <- "tau2"
      tau2.pos <- tau2
   }

   if (!missing(rho)) {
      vc <- x$rho[rho]
      comp <- "rho"
      rho.pos <- rho
   }

   if (!missing(gamma2)) {
      vc <- x$gamma2[gamma2]
      comp <- "gamma2"
      gamma2.pos <- gamma2
   }

   if (!missing(phi)) {
      vc <- x$phi[phi]
      comp <- "phi"
      phi.pos <- phi
   }

   #return(list(comp=comp, vc=vc))

   if (missing(xlim)) {

      ### if the user has not specified xlim, set it automatically
      ### TODO: maybe try something based on CI later

      if (comp == "sigma2") {
         vc.lb <- max( 0, vc/4)
         vc.ub <- max(.1, vc*4)
      }
      if (comp == "tau2") {
         if (is.element(x$struct[1], c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD"))) {
            vc.lb <- max( 0, vc/2)
            vc.ub <- max(.1, vc*2)
         } else {
            vc.lb <- max( 0, vc/4)
            vc.ub <- max(.1, vc*4)
         }
      }
      if (comp == "gamma2") {
         if (is.element(x$struct[2], c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD"))) {
            vc.lb <- max( 0, vc/2)
            vc.ub <- max(.1, vc*2)
         } else {
            vc.lb <- max( 0, vc/4)
            vc.ub <- max(.1, vc*4)
         }
      }
      if (comp == "rho") {
         if (x$struct[1] == "CAR") {
            vc.lb <- max(0, vc-.5)
            vc.ub <- min(+.99999, vc+.5)
         }
         if (is.element(x$struct[1], c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYPL","PHYPD"))) {
            vc.lb <- vc/2
            vc.ub <- vc*2
         }
         if (!is.element(x$struct[1], c("CAR","SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYPL","PHYPD"))) {
            vc.lb <- max(-.99999, vc-.5)
            vc.ub <- min(+.99999, vc+.5)
         }
      }
      if (comp == "phi") {
         if (x$struct[2] == "CAR") {
            vc.lb <- max(0, vc-.5)
            vc.ub <- min(+.99999, vc+.5)
         }
         if (is.element(x$struct[2], c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYPL","PHYPD"))) {
            vc.lb <- vc/2
            vc.ub <- vc*2
         }
         if (!is.element(x$struct[2], c("CAR","SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYPL","PHYPD"))) {
            vc.lb <- max(-.99999, vc-.5)
            vc.ub <- min(+.99999, vc+.5)
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

      if (is.element(comp, c("sigma2", "tau2", "gamma2"))) {
         if (xlim[1] < 0)
            stop(mstyle$stop("Lower bound for profiling must be >= 0."))
      }
      if (comp == "rho") {
         if (is.element(x$struct[1], c("CAR","SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYPL","PHYPD")) && xlim[1] < 0)
            stop(mstyle$stop("Lower bound for profiling must be >= 0."))
         if (xlim[1] < -1)
            stop(mstyle$stop("Lower bound for profiling must be >= -1."))
         if (!is.element(x$struct[1], c("CAR","SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYPL","PHYPD")) && xlim[2] > 1)
            stop(mstyle$stop("Upper bound for profiling must be <= 1."))
      }
      if (comp == "phi") {
         if (is.element(x$struct[2], c("CAR","SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYPL","PHYPD")) && xlim[1] < 0)
            stop(mstyle$stop("Lower bound for profiling must be >= 0."))
         if (xlim[1] < -1)
            stop(mstyle$stop("Lower bound for profiling must be >= -1."))
         if (!is.element(x$struct[2], c("CAR","SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYPL","PHYPD")) && xlim[2] > 1)
            stop(mstyle$stop("Upper bound for profiling must be <= 1."))
      }

   }

   vcs <- seq(xlim[1], xlim[2], length.out=steps)
   #return(vcs)

   if (length(vcs) <= 1L)
      stop(mstyle$stop("Cannot set 'xlim' automatically. Please set this argument manually."))

   if (parallel == "no")
      res <- pbapply::pblapply(vcs, .profile.rma.mv, obj=x, comp=comp, sigma2.pos=sigma2.pos, tau2.pos=tau2.pos, rho.pos=rho.pos, gamma2.pos=gamma2.pos, phi.pos=phi.pos, parallel=parallel, profile=TRUE)

   if (parallel == "multicore")
      res <- pbapply::pblapply(vcs, .profile.rma.mv, obj=x, comp=comp, sigma2.pos=sigma2.pos, tau2.pos=tau2.pos, rho.pos=rho.pos, gamma2.pos=gamma2.pos, phi.pos=phi.pos, parallel=parallel, profile=TRUE, cl=ncpus)
      #res <- parallel::mclapply(vcs, .profile.rma.mv, obj=x, comp=comp, sigma2.pos=sigma2.pos, tau2.pos=tau2.pos, rho.pos=rho.pos, gamma2.pos=gamma2.pos, phi.pos=phi.pos, parallel=parallel, profile=TRUE, mc.cores=ncpus)

   if (parallel == "snow") {
      if (is.null(cl)) {
         cl <- parallel::makePSOCKcluster(ncpus)
         on.exit(parallel::stopCluster(cl), add=TRUE)
      }
      if (.isTRUE(ddd$LB)) {
         res <- parallel::parLapplyLB(cl, vcs, .profile.rma.mv, obj=x, comp=comp, sigma2.pos=sigma2.pos, tau2.pos=tau2.pos, rho.pos=rho.pos, gamma2.pos=gamma2.pos, phi.pos=phi.pos, parallel=parallel, profile=TRUE)
         #res <- parallel::clusterApplyLB(cl, vcs, .profile.rma.mv, obj=x, comp=comp, sigma2.pos=sigma2.pos, tau2.pos=tau2.pos, rho.pos=rho.pos, gamma2.pos=gamma2.pos, phi.pos=phi.pos, parallel=parallel, profile=TRUE)
      } else {
         res <- pbapply::pblapply(vcs, .profile.rma.mv, obj=x, comp=comp, sigma2.pos=sigma2.pos, tau2.pos=tau2.pos, rho.pos=rho.pos, gamma2.pos=gamma2.pos, phi.pos=phi.pos, parallel=parallel, profile=TRUE, cl=cl)
         #res <- parallel::parLapply(cl, vcs, .profile.rma.mv, obj=x, comp=comp, sigma2.pos=sigma2.pos, tau2.pos=tau2.pos, rho.pos=rho.pos, gamma2.pos=gamma2.pos, phi.pos=phi.pos, parallel=parallel, profile=TRUE)
         #res <- parallel::clusterApply(cl, vcs, .profile.rma.mv, obj=x, comp=comp, sigma2.pos=sigma2.pos, tau2.pos=tau2.pos, rho.pos=rho.pos, gamma2.pos=gamma2.pos, phi.pos=phi.pos, parallel=parallel, profile=TRUE)
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

      if (any(!is.na(lls))) {
         if (xlim[1] <= vc && xlim[2] >= vc) {
            ylim <- range(c(logLik(x),lls), na.rm=TRUE)
         } else {
            ylim <- range(lls, na.rm=TRUE)
         }
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

   if (comp == "sigma2") {
      if (x$sigma2s == 1L) {
         xlab <- expression(paste(sigma^2, " Value"))
         title <- expression(paste("Profile Plot for ", sigma^2))
      } else {
         xlab <- bquote(sigma[.(sigma2)]^2 ~ "Value")
         title <- bquote("Profile Plot for" ~ sigma[.(sigma2)]^2)
      }
   }
   if (comp == "tau2") {
      if (x$tau2s == 1L) {
         xlab <- expression(paste(tau^2, " Value"))
         title <- expression(paste("Profile Plot for ", tau^2))
      } else {
         xlab <- bquote(tau[.(tau2)]^2 ~ "Value")
         title <- bquote("Profile Plot for" ~ tau[.(tau2)]^2)
      }
   }
   if (comp == "rho") {
      if (x$rhos == 1L) {
         xlab <- expression(paste(rho, " Value"))
         title <- expression(paste("Profile Plot for ", rho))
      } else {
         xlab <- bquote(rho[.(rho)] ~ "Value")
         title <- bquote("Profile Plot for" ~ rho[.(rho)])
      }
   }
   if (comp == "gamma2") {
      if (x$gamma2s == 1L) {
         xlab <- expression(paste(gamma^2, " Value"))
         title <- expression(paste("Profile Plot for ", gamma^2))
      } else {
         xlab <- bquote(gamma[.(gamma2)]^2 ~ "Value")
         title <- bquote("Profile Plot for" ~ gamma[.(gamma2)]^2)
      }
   }
   if (comp == "phi") {
      if (x$phis == 1L) {
         xlab <- expression(paste(phi, " Value"))
         title <- expression(paste("Profile Plot for ", phi))
      } else {
         xlab <- bquote(phi[.(phi)] ~ "Value")
         title <- bquote("Profile Plot for" ~ phi[.(phi)])
      }
   }

   sav <- list(vc=vcs, ll=lls, beta=beta, ci.lb=ci.lb, ci.ub=ci.ub, comps=1, ylim=ylim, method=x$method, vc=vc, maxll=logLik(x), xlab=xlab, title=title)
   names(sav)[1] <- switch(comp, sigma2="sigma2", tau2="tau2", rho="rho", gamma2="gamma2", phi="phi")
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
