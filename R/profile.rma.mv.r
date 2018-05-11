profile.rma.mv <- function(fitted, sigma2, tau2, rho, gamma2, phi, xlim, ylim, steps=20, lltol=1e-06, startmethod="init", progbar=TRUE, parallel="no", ncpus=1, cl=NULL, plot=TRUE, pch=19, cline=FALSE, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(fitted, "rma.mv"))
      stop(mstyle$stop("Argument 'fitted' must be an object of class \"rma.mv\"."))

   if (steps < 2)
      stop(mstyle$stop("Argument 'steps' must be >= 2."))

   x <- fitted

   parallel <- match.arg(parallel, c("no", "snow", "multicore"))

   if (parallel == "no" && ncpus > 1)
      parallel <- "snow"

   #########################################################################

   ### check if user has specified one of the sigma2, tau2, rho, gamma2, or phi arguments

   if (missing(sigma2) && missing(tau2) && missing(rho) && missing(gamma2) && missing(phi)) {

      mc <- match.call()

      ### total number of non-fixed components

      comps <- ifelse(x$withS, sum(!x$vc.fix$sigma2), 0) + ifelse(x$withG, sum(!x$vc.fix$tau2) + sum(!x$vc.fix$rho), 0) + ifelse(x$withH, sum(!x$vc.fix$gamma2) + sum(!x$vc.fix$phi), 0)

      if (comps == 0)
         stop(mstyle$stop("No components to profile."))

      if (plot) {
         if (dev.cur() == 1) {
            par(mfrow=c(comps, 1))
            #on.exit(par(mfrow=c(1,1)))
         }
      }

      sav <- list()
      j <- 0

      if (x$withS && any(!x$vc.fix$sigma2)) {
         for (pos in seq_len(x$sigma2s)[!x$vc.fix$sigma2]) {
            j <- j + 1
            mc.vc <- mc
            mc.vc$sigma2 <- pos
            mc.vc$fitted <- quote(x)
            if (progbar)
               cat(mstyle$verbose(paste("Profiling sigma2 =", pos, "\n")))
            sav[[j]] <- eval(mc.vc)
         }
      }

      if (x$withG) {
         if (any(!x$vc.fix$tau2)) {
            for (pos in seq_len(x$tau2s)[!x$vc.fix$tau2]) {
               j <- j + 1
               mc.vc <- mc
               mc.vc$tau2 <- pos
               mc.vc$fitted <- quote(x)
               if (progbar)
                  cat(mstyle$verbose(paste("Profiling tau2 =", pos, "\n")))
               sav[[j]] <- eval(mc.vc)
            }
         }
         if (any(!x$vc.fix$rho)) {
            for (pos in seq_len(x$rhos)[!x$vc.fix$rho]) {
               j <- j + 1
               mc.vc <- mc
               mc.vc$rho <- pos
               mc.vc$fitted <- quote(x)
               if (progbar)
                  cat(mstyle$verbose(paste("Profiling rho =", pos, "\n")))
               sav[[j]] <- eval(mc.vc)
            }
         }
      }

      if (x$withH) {
         if (any(!x$vc.fix$gamma2)) {
            for (pos in seq_len(x$gamma2s)[!x$vc.fix$gamma2]) {
               j <- j + 1
               mc.vc <- mc
               mc.vc$gamma2 <- pos
               mc.vc$fitted <- quote(x)
               if (progbar)
                  cat(mstyle$verbose(paste("Profiling gamma2 =", pos, "\n")))
               sav[[j]] <- eval(mc.vc)
            }
         }
         if (any(!x$vc.fix$phi)) {
            for (pos in seq_len(x$phis)[!x$vc.fix$phi]) {
               j <- j + 1
               mc.vc <- mc
               mc.vc$phi <- pos
               mc.vc$fitted <- quote(x)
               if (progbar)
                  cat(mstyle$verbose(paste("Profiling phi =", pos, "\n")))
               sav[[j]] <- eval(mc.vc)
            }
         }
      }

      if (comps == 1)
         sav <- sav[[1]]

      sav$comps <- comps

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

      ### if the user has not specified xlim argument, set automatically
      ### TODO: maybe try something based on CI later

      if (is.element(comp, c("sigma2", "tau2", "gamma2"))) {
         #vc.lb <- max(.00001, log(vc)) ### old method
         #vc.ub <- max(.00001, exp(vc)) ### old method
         vc.lb <- max( 0, vc/4) ### new method
         vc.ub <- max(.1, vc*4) ### new method
      }
      if (comp == "rho") {
         if (x$struct[1] == "CAR") {
            vc.lb <- max(0, vc-.5)
            vc.ub <- min(+.99999, vc+.5)
         }
         if (is.element(x$struct[1], c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH"))) {
            vc.lb <- vc/4
            vc.ub <- vc*4
         }
         if (!is.element(x$struct[1], c("CAR","SPEXP","SPGAU","SPLIN","SPRAT","SPSPH"))) {
            vc.lb <- max(-.99999, vc-.5)
            vc.ub <- min(+.99999, vc+.5)
         }
      }
      if (comp == "phi") {
         if (x$struct[2] == "CAR") {
            vc.lb <- max(0, vc-.5)
            vc.ub <- min(+.99999, vc+.5)
         }
         if (is.element(x$struct[2], c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH"))) {
            vc.lb <- vc/4
            vc.ub <- vc*4
         }
         if (!is.element(x$struct[2], c("CAR","SPEXP","SPGAU","SPLIN","SPRAT","SPSPH"))) {
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
         if (is.element(x$struct[1], c("CAR","SPEXP","SPGAU","SPLIN","SPRAT","SPSPH")) && xlim[1] < 0)
            stop(mstyle$stop("Lower bound for profiling must be >= 0."))
         if (xlim[1] < -1)
            stop(mstyle$stop("Lower bound for profiling must be >= -1."))
         if (!is.element(x$struct[1], c("CAR","SPEXP","SPGAU","SPLIN","SPRAT","SPSPH")) && xlim[2] > 1)
            stop(mstyle$stop("Upper bound for profiling must be <= -1."))
      }
      if (comp == "phi") {
         if (is.element(x$struct[2], c("CAR","SPEXP","SPGAU","SPLIN","SPRAT","SPSPH")) && xlim[1] < 0)
            stop(mstyle$stop("Lower bound for profiling must be >= 0."))
         if (xlim[1] < -1)
            stop(mstyle$stop("Lower bound for profiling must be >= -1."))
         if (!is.element(x$struct[2], c("CAR","SPEXP","SPGAU","SPLIN","SPRAT","SPSPH")) && xlim[2] > 1)
            stop(mstyle$stop("Upper bound for profiling must be <= -1."))
      }

   }

   vcs <- seq(xlim[1], xlim[2], length=steps)

   if (length(vcs) <= 1)
      stop(mstyle$stop("Cannot set 'xlim' automatically. Please set this argument manually."))

   #return(vcs)

   ### extract control argument
   x.control <- x$control

   ### if not an empty list(), get position of sigma2.init, tau2.init, rho.init, gamma2.init, or phi.init arguments
   ### if these arguments were not specified, then the respective con.pos values are NA

   if (length(x.control) > 0) {
      con.pos.sigma2.init <- pmatch("sigma2.init", names(x.control))
      con.pos.tau2.init   <- pmatch("tau2.init",   names(x.control))
      con.pos.rho.init    <- pmatch("rho.init",    names(x.control))
      con.pos.gamma2.init <- pmatch("gamma2.init", names(x.control))
      con.pos.phi.init    <- pmatch("phi.init",    names(x.control))
   } else {
      con.pos.sigma2.init <- NA
      con.pos.tau2.init   <- NA
      con.pos.rho.init    <- NA
      con.pos.gamma2.init <- NA
      con.pos.phi.init    <- NA
   }

   if (parallel=="no") {

      lls   <- rep(NA_real_, length(vcs))
      beta  <- matrix(NA_real_, nrow=length(vcs), ncol=x$p)
      ci.lb <- matrix(NA_real_, nrow=length(vcs), ncol=x$p)
      ci.ub <- matrix(NA_real_, nrow=length(vcs), ncol=x$p)

      if (progbar)
         pbar <- txtProgressBar(min=0, max=steps, style=3)

      ### set any fixed components to their values

      sigma2.arg <- ifelse(x$vc.fix$sigma2, x$sigma2, NA)
      tau2.arg   <- ifelse(x$vc.fix$tau2, x$tau2, NA)
      rho.arg    <- ifelse(x$vc.fix$rho, x$rho, NA)
      gamma2.arg <- ifelse(x$vc.fix$gamma2, x$gamma2, NA)
      phi.arg    <- ifelse(x$vc.fix$phi, x$phi, NA)

      for (i in seq_along(vcs)) {

         if (comp == "sigma2")
            sigma2.arg[sigma2] <- vcs[i]

         if (comp == "tau2")
            tau2.arg[tau2] <- vcs[i]

         if (comp == "rho")
            rho.arg[rho] <- vcs[i]

         if (comp == "gamma2")
            gamma2.arg[gamma2] <- vcs[i]

         if (comp == "phi")
            phi.arg[phi] <- vcs[i]

         ### if this is the second model fit and the previous model fit was not a try-error

         if (startmethod == "prev" && i > 1 && !inherits(res, "try-error")) {

            ### set the sigma2.init argument to the estimated sigma2 value(s) from the previous model fit

            if (is.na(con.pos.sigma2.init)) {
               x.control$sigma2.init <- res$sigma2
            } else {
               x.control[[con.pos.sigma2.init]] <- res$sigma2
            }

            ### set the tau2.init argument to the estimated tau2 value(s) from the previous model fit

            if (is.na(con.pos.tau2.init)) {
               x.control$tau2.init <- res$tau2
            } else {
               x.control[[con.pos.tau2.init]] <- res$tau2
            }

            ### set the rho.init argument to the estimated rho value(s) from the previous model fit

            if (is.na(con.pos.rho.init)) {
               x.control$rho.init <- res$rho
            } else {
               x.control[[con.pos.rho.init]] <- res$rho
            }

            ### set the gamma2.init argument to the estimated gamma2 value(s) from the previous model fit

            if (is.na(con.pos.gamma2.init)) {
               x.control$gamma2.init <- res$gamma2
            } else {
               x.control[[con.pos.gamma2.init]] <- res$gamma2
            }

            ### set the phi.init argument to the estimated phi value(s) from the previous model fit

            if (is.na(con.pos.phi.init)) {
               x.control$phi.init <- res$phi
            } else {
               x.control[[con.pos.phi.init]] <- res$phi
            }

            res <- try(suppressWarnings(rma.mv(x$yi, V=x$V, W=x$W, mods=x$X, random=x$random, struct=x$struct, intercept=FALSE, data=x$mf.r, method=x$method, test=x$test, level=x$level, R=x$R, Rscale=x$Rscale, sigma2=sigma2.arg, tau2=tau2.arg, rho=rho.arg, gamma2=gamma2.arg, phi=phi.arg, sparse=x$sparse, dist=x$dist, control=x.control)), silent=TRUE)

         } else {

            res <- try(suppressWarnings(rma.mv(x$yi, V=x$V, W=x$W, mods=x$X, random=x$random, struct=x$struct, intercept=FALSE, data=x$mf.r, method=x$method, test=x$test, level=x$level, R=x$R, Rscale=x$Rscale, sigma2=sigma2.arg, tau2=tau2.arg, rho=rho.arg, gamma2=gamma2.arg, phi=phi.arg, sparse=x$sparse, dist=x$dist, control=x$control)), silent=TRUE)

         }

         if (progbar)
            setTxtProgressBar(pbar, i)

         if (inherits(res, "try-error"))
            next

         lls[i] <- c(logLik(res))
         beta[i,]  <- c(res$beta)
         ci.lb[i,] <- c(res$ci.lb)
         ci.ub[i,] <- c(res$ci.ub)

      }

      if (progbar)
         close(pbar)

   }

   if (parallel=="snow" || parallel == "multicore") {

      if (!requireNamespace("parallel", quietly=TRUE))
         stop(mstyle$stop("Please install the 'parallel' package for parallel processing."))

      ncpus <- as.integer(ncpus)

      if (ncpus < 1)
         stop(mstyle$stop("Argument 'ncpus' must be >= 1."))

      if (parallel == "multicore")
         res <- parallel::mclapply(vcs, .profile.rma.mv, obj=x, comp=comp, sigma2.pos=sigma2.pos, tau2.pos=tau2.pos, rho.pos=rho.pos, gamma2.pos=gamma2.pos, phi.pos=phi.pos, mc.cores=ncpus, parallel=parallel, profile=TRUE)

      if (parallel == "snow") {
         if (is.null(cl)) {
            cl <- parallel::makePSOCKcluster(ncpus)
            res <- parallel::parLapply(cl, vcs, .profile.rma.mv, obj=x, comp=comp, sigma2.pos=sigma2.pos, tau2.pos=tau2.pos, rho.pos=rho.pos, gamma2.pos=gamma2.pos, phi.pos=phi.pos, parallel=parallel, profile=TRUE)
            #res <- parallel::parLapplyLB(cl, vcs, .profile.rma.mv, obj=x, comp=comp, sigma2.pos=sigma2.pos, tau2.pos=tau2.pos, rho.pos=rho.pos, gamma2.pos=gamma2.pos, phi.pos=phi.pos, parallel=parallel, profile=TRUE)
            #res <- parallel::clusterApply(cl, vcs, .profile.rma.mv, obj=x, comp=comp, sigma2.pos=sigma2.pos, tau2.pos=tau2.pos, rho.pos=rho.pos, gamma2.pos=gamma2.pos, phi.pos=phi.pos, parallel=parallel, profile=TRUE)
            #res <- parallel::clusterApplyLB(cl, vcs, .profile.rma.mv, obj=x, comp=comp, sigma2.pos=sigma2.pos, tau2.pos=tau2.pos, rho.pos=rho.pos, gamma2.pos=gamma2.pos, phi.pos=phi.pos, parallel=parallel, profile=TRUE)
            parallel::stopCluster(cl)
         } else {
            res <- parallel::parLapply(cl, vcs, .profile.rma.mv, obj=x, comp=comp, sigma2.pos=sigma2.pos, tau2.pos=tau2.pos, rho.pos=rho.pos, gamma2.pos=gamma2.pos, phi.pos=phi.pos, parallel=parallel, profile=TRUE)
         }
      }

      lls <- sapply(res, function(z) z$ll)
      beta  <- do.call("rbind", lapply(res, function(z) t(z$beta)))
      ci.lb <- do.call("rbind", lapply(res, function(z) t(z$ci.lb)))
      ci.ub <- do.call("rbind", lapply(res, function(z) t(z$ci.ub)))

   }

   #########################################################################

   if (any(lls >= logLik(x) + lltol, na.rm=TRUE))
      warning(mstyle$warning("At least one profiled log-likelihood value is larger than the log-likelihood of the fitted model."))

   beta  <- data.frame(beta)
   ci.lb <- data.frame(ci.lb)
   ci.ub <- data.frame(ci.ub)
   names(beta)  <- rownames(x$beta)
   names(ci.lb) <- rownames(x$beta)
   names(ci.ub) <- rownames(x$beta)

   if (missing(ylim)) {

      ylim <- range(lls, na.rm=TRUE)
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
      plot(sav, pch=pch, cline=cline, ...)

   #########################################################################

   invisible(sav)

}
