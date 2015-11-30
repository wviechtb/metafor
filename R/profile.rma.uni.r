profile.rma.uni <- function(fitted, xlim, ylim, steps=20, progbar=TRUE, parallel="no", ncpus=1, cl=NULL, plot=TRUE, pch=19, ...) {

   if (!is.element("rma.uni", class(fitted)))
      stop("Argument 'fitted' must be an object of class \"rma.uni\".")

   if (steps < 2)
      stop("Argument 'steps' must be >= 2.")

   x <- fitted

   parallel <- match.arg(parallel, c("no", "snow", "multicore"))

   #########################################################################

   if (missing(xlim)) {

      ### if the user has not specified the xlim, get CI xlim for tau^2 (suppress warnings)

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

         ### if the profile method fails, try a Wald-type CI for tau^2

         vc.lb <- max(0, x$tau2 - 1.96 * x$se.tau2)
         vc.ub <- max(.1, x$tau2 + 1.96 * x$se.tau2)

      }

      if (is.na(vc.lb) || is.na(vc.ub)) {

         ### if this still results in NA bounds, use simple method

         #vc.lb <- max(0, log(x$tau2)) ### old method
         #vc.ub <- max(0, exp(x$tau2)) ### old method
         vc.lb <- max( 0, x$tau2/4) ### new method
         vc.ub <- max(.1, x$tau2*4) ### new method

      }

      ### if all of that fails, throw an error

      if (is.na(vc.lb) || is.na(vc.ub))
         stop("Cannot set 'xlim' automatically. Please set this argument manually.")

      xlim <- c(vc.lb, vc.ub)

   } else {

      if (length(xlim) != 2L)
         stop("Argument 'xlim' should be a vector of length 2.")

      xlim <- sort(xlim)

   }

   vcs <- seq(xlim[1], xlim[2], length=steps)

   if (length(vcs) <= 1)
      stop("Cannot set 'xlim' automatically. Please set this argument manually.")

   if (parallel=="no") {

      lls   <- rep(NA_real_, length(vcs))
      b     <- matrix(NA_real_, nrow=length(vcs), ncol=x$p)
      ci.lb <- matrix(NA_real_, nrow=length(vcs), ncol=x$p)
      ci.ub <- matrix(NA_real_, nrow=length(vcs), ncol=x$p)

      if (progbar)
         pbar <- txtProgressBar(min=0, max=steps, style=3)

      for (i in 1:length(vcs)) {

         res <- try(suppressWarnings(rma(x$yi, x$vi, weights=x$weights, mods=x$X, method=x$method, weighted=x$weighted, intercept=FALSE, knha=x$knha, level=x$level, control=x$control, tau2=vcs[i])), silent=TRUE)

         if (inherits(res, "try-error"))
            next

         lls[i] <- c(logLik(res))
         b[i,]  <- c(res$b)
         ci.lb[i,] <- c(res$ci.lb)
         ci.ub[i,] <- c(res$ci.ub)

         if (progbar)
            setTxtProgressBar(pbar, i)

      }

      if (progbar)
         close(pbar)

   }

   if (parallel=="snow" || parallel == "multicore") {

      if (!requireNamespace("parallel", quietly=TRUE))
         stop("Please install the 'parallel' package for parallel processing.")

      ncpus <- as.integer(ncpus)

      if (ncpus < 1)
         stop("Argument 'ncpus' must be >= 1.")

      if (parallel == "multicore")
         res <- parallel::mclapply(vcs, .profile.rma.uni, obj=x, mc.cores=ncpus, parallel=parallel)

      if (parallel == "snow") {
         if (is.null(cl)) {
            cl <- parallel::makePSOCKcluster(ncpus)
            res <- parallel::parLapply(cl, vcs, .profile.rma.uni, obj=x, parallel=parallel)
            parallel::stopCluster(cl)
         } else {
            res <- parallel::parLapply(cl, vcs, .profile.rma.uni, obj=x, parallel=parallel)
         }
      }

      lls <- sapply(res, function(z) z$ll)
      b   <- do.call("rbind", lapply(res, function(z) t(z$b)))
      ci.lb <- do.call("rbind", lapply(res, function(z) t(z$ci.lb)))
      ci.ub <- do.call("rbind", lapply(res, function(z) t(z$ci.ub)))

   }

   #########################################################################

   b <- data.frame(b)
   ci.lb <- data.frame(ci.lb)
   ci.ub <- data.frame(ci.ub)
   names(b) <- rownames(x$b)
   names(ci.lb) <- rownames(x$b)
   names(ci.ub) <- rownames(x$b)

   res <- list(tau2=vcs, ll=lls, b=b, ci.lb=ci.lb, ci.ub=ci.ub)

   class(res) <- "profile.rma.uni"

   #########################################################################

   if (plot) {

      if (missing(ylim)) {
         ylim <- range(lls, na.rm=TRUE)
         ylim[1] <- ylim[1] - .1
         ylim[2] <- ylim[2] + .1
      }

      xlab <- expression(paste(tau^2, " Value"))
      title <- expression(paste("Profile Plot for ", tau^2))

      plot(vcs, lls, type="o", xlab=xlab, ylab=paste(ifelse(x$method=="REML", "Restricted", ""), " Log-Likelihood", sep=""), main=title, bty="l", pch=pch, ylim=ylim, ...)
      abline(v=x$tau2,    lty="dotted")
      abline(h=logLik(x), lty="dotted")
      #abline(h=max(lls, na.rm=TRUE), lty="dotted")

   }

   #########################################################################

   invisible(res)

}
