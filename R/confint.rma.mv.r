confint.rma.mv <- function(object, parm, level, fixed=FALSE, sigma2, tau2, rho, gamma2, phi, digits, transf, targs, verbose=FALSE, control, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(object), must="rma.mv")

   if (!missing(parm))
      warning(mstyle$warning("Argument 'parm' (currently) ignored."), call.=FALSE)

   x <- object

   k <- x$k
   p <- x$p

   if (missing(level))
      level <- x$level

   if (missing(digits)) {
      digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
   }

   if (missing(transf))
      transf <- FALSE

   if (missing(targs))
      targs <- NULL

   if (missing(control))
      control <- list()

   level <- ifelse(level == 0, 1, ifelse(level >= 1, (100-level)/100, ifelse(level > .5, 1-level, level)))

   ddd <- list(...)

   .chkdots(ddd, c("time", "xlim"))

   if (.isTRUE(ddd$time))
      time.start <- proc.time()

   if (!is.null(ddd$xlim)) {
      if (length(ddd$xlim) != 2L)
         stop(mstyle$stop("Argument 'xlim' should be a vector of length 2."))
      control$vc.min <- ddd$xlim[1]
      control$vc.max <- ddd$xlim[2]
   }

   ### check if user has specified one of the sigma2, tau2, rho, gamma2, or phi arguments

   random <- !all(missing(sigma2), missing(tau2), missing(rho), missing(gamma2), missing(phi))

   if (!fixed && !random) {

      ### if both 'fixed' and 'random' are FALSE, obtain CIs for all variance/correlation components

      cl <- match.call()

      ### total number of non-fixed components

      comps <- ifelse(x$withS, sum(!x$vc.fix$sigma2), 0) +
               ifelse(x$withG, sum(!x$vc.fix$tau2) + sum(!x$vc.fix$rho), 0) +
               ifelse(x$withH, sum(!x$vc.fix$gamma2) + sum(!x$vc.fix$phi), 0)

      if (comps == 0)
         stop(mstyle$stop("No components for which a CI can be obtained."))

      res.all <- list()
      j <- 0

      if (x$withS && any(!x$vc.fix$sigma2)) {
         for (pos in seq_len(x$sigma2s)[!x$vc.fix$sigma2]) {
            j <- j + 1
            cl.vc <- cl
            cl.vc$sigma2 <- pos
            cl.vc$time <- FALSE
            #cl.vc$object <- quote(x)
            if (verbose)
               cat(mstyle$verbose(paste("\nObtaining CI for sigma2 =", pos, "\n")))
            res.all[[j]] <- eval(cl.vc, envir=parent.frame())
         }
      }

      if (x$withG) {
         if (any(!x$vc.fix$tau2)) {
            for (pos in seq_len(x$tau2s)[!x$vc.fix$tau2]) {
               j <- j + 1
               cl.vc <- cl
               cl.vc$tau2 <- pos
               cl.vc$time <- FALSE
               #cl.vc$object <- quote(x)
               if (verbose)
                  cat(mstyle$verbose(paste("\nObtaining CI for tau2 =", pos, "\n")))
               res.all[[j]] <- eval(cl.vc, envir=parent.frame())
            }
         }
         if (any(!x$vc.fix$rho)) {
            for (pos in seq_len(x$rhos)[!x$vc.fix$rho]) {
               j <- j + 1
               cl.vc <- cl
               cl.vc$rho <- pos
               cl.vc$time <- FALSE
               #cl.vc$object <- quote(x)
               if (verbose)
                  cat(mstyle$verbose(paste("\nObtaining CI for rho =", pos, "\n")))
               res.all[[j]] <- eval(cl.vc, envir=parent.frame())
            }
         }
      }

      if (x$withH) {
         if (any(!x$vc.fix$gamma2)) {
            for (pos in seq_len(x$gamma2s)[!x$vc.fix$gamma2]) {
               j <- j + 1
               cl.vc <- cl
               cl.vc$gamma2 <- pos
               cl.vc$time <- FALSE
               #cl.vc$object <- quote(x)
               if (verbose)
                  cat(mstyle$verbose(paste("\nObtaining CI for gamma2 =", pos, "\n")))
               res.all[[j]] <- eval(cl.vc, envir=parent.frame())
            }
         }
         if (any(!x$vc.fix$phi)) {
            for (pos in seq_len(x$phis)[!x$vc.fix$phi]) {
               j <- j + 1
               cl.vc <- cl
               cl.vc$phi <- pos
               cl.vc$time <- FALSE
               #cl.vc$object <- quote(x)
               if (verbose)
                  cat(mstyle$verbose(paste("\nObtaining CI for phi =", pos, "\n")))
               res.all[[j]] <- eval(cl.vc, envir=parent.frame())
            }
         }
      }

      if (.isTRUE(ddd$time)) {
         time.end <- proc.time()
         .print.time(unname(time.end - time.start)[3])
      }

      if (length(res.all) == 1L) {
         return(res.all[[1]])
      } else {
         res.all$digits <- digits
         class(res.all) <- "list.confint.rma"
         return(res.all)
      }

   }

   #########################################################################
   #########################################################################
   #########################################################################

   if (random) {

      type <- "PL"

      ######################################################################

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

      #return(list(comp=comp, vc=vc, sigma2.pos=sigma2.pos, tau2.pos=tau2.pos, rho.pos=rho.pos, gamma2.pos=gamma2.pos, phi.pos=phi.pos))

      ######################################################################

      ### set control parameters for uniroot() and possibly replace with user-defined values
      ### set vc.min and vc.max and possibly replace with user-defined values

      con <- list(tol=.Machine$double.eps^0.25, maxiter=1000, verbose=FALSE, eptries=10)

      if (is.element(comp, c("sigma2", "tau2", "gamma2"))) {
         con$vc.min <- 0
         con$vc.max <- max(ifelse(vc <= .Machine$double.eps^0.5, 10, max(10, vc*100)), con$vc.min)
      }
      if (comp == "rho") {
         if (is.element(x$struct[1], c("CS","HCS")))
            con$vc.min <- -1                                  ### this will fail most of the time but with retries, this may get closer to actual lower bound
            #con$vc.min <- min(-1/(x$g.nlevels.f[1] - 1), vc) ### this guarantees that cor matrix is semi-positive definite, but since V gets added, this is actually too strict
         if (is.element(x$struct[1], c("AR","HAR","CAR")))
            con$vc.min <- min(0, vc)                          ### negative autocorrelation parameters not considered (not even sensible for CAR)
         if (is.element(x$struct[1], c("UN","UNR","GEN")))
            con$vc.min <- -1                                  ### TODO: this will often fail! (but with retries, this should still work)
         con$vc.max <- 1
         if (is.element(x$struct[1], c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH"))) {
            con$vc.min <- 0                                   ### TODO: 0 basically always fails
            con$vc.max <- max(10, vc*10)
         }
         if (is.element(x$struct[1], c("PHYPL","PHYPD"))) {
            con$vc.min <- 0
            con$vc.max <- max(2, vc*2)
         }
      }
      if (comp == "phi") {
         if (is.element(x$struct[2], c("CS","HCS")))
            con$vc.min <- -1                                  ### this will fail most of the time but with retries, this may get closer to actual lower bound
            #con$vc.min <- min(-1/(x$h.nlevels.f[1] - 1), vc) ### this guarantees that cor matrix is semi-positive definite, but since V gets added, this is actually too strict
         if (is.element(x$struct[2], c("AR","HAR","CAR")))
            con$vc.min <- min(0, vc)                          ### negative autocorrelation parameters not considered (not even sensible for CAR)
         if (is.element(x$struct[2], c("UN","UNR","GEN")))
            con$vc.min <- -1                                  ### TODO: this will often fail! (but with retries, this should still work)
         con$vc.max <- 1
         if (is.element(x$struct[2], c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH"))) {
            con$vc.min <- 0                                   ### TODO: 0 basically always fails
            con$vc.max <- max(10, vc*10)
         }
         if (is.element(x$struct[2], c("PHYPL","PHYPD"))) {
            con$vc.min <- 0
            con$vc.max <- max(2, vc*2)
         }
      }

      con.pos <- pmatch(names(control), names(con))
      con[c(na.omit(con.pos))] <- control[!is.na(con.pos)]

      if (verbose)
         con$verbose <- verbose

      verbose <- con$verbose

      ######################################################################

      vc.lb <- NA
      vc.ub <- NA
      ci.null <- FALSE ### logical if CI is a null set
      lb.conv <- FALSE ### logical if search converged for lower bound (LB)
      ub.conv <- FALSE ### logical if search converged for upper bound (UB)
      lb.sign <- ""    ### for sign in case LB must be below vc.min ("<") or above vc.max (">")
      ub.sign <- ""    ### for sign in case UB must be below vc.min ("<") or above vc.max (">")

      ######################################################################
      ######################################################################
      ######################################################################

      ### Profile Likelihood method

      if (type == "PL") {

         if (con$vc.min > vc)
            stop(mstyle$stop("Lower bound of interval to be searched must be <= estimated value of component."))
         if (con$vc.max < vc)
            stop(mstyle$stop("Upper bound of interval to be searched must be >= estimated value of component."))

         objective <- qchisq(1-level, df=1)

         ###################################################################

         ### search for lower bound

         ### get diff value when setting component to vc.min; this value should be positive (i.e., discrepancy must be larger than critical value)
         ### if it is not, then the lower bound must be below vc.min

         epdiff <- abs(con$vc.min - vc) / con$eptries

         for (i in seq_len(con$eptries)) {

            res <- try(.profile.rma.mv(con$vc.min, obj=x, comp=comp, sigma2.pos=sigma2.pos, tau2.pos=tau2.pos, rho.pos=rho.pos, gamma2.pos=gamma2.pos, phi.pos=phi.pos, confint=TRUE, objective=objective, verbose=verbose), silent=TRUE)

            if (!inherits(res, "try-error")) {

               if (res < 0) {

                  vc.lb <- con$vc.min
                  lb.conv <- TRUE

                  if (is.element(comp, c("sigma2", "tau2", "gamma2")) && con$vc.min > 0)
                     lb.sign <- "<"

                  if (is.element(comp, c("rho", "phi")) && con$vc.min > -1)
                     lb.sign <- "<"

                  if (((comp == "rho" && is.element(x$struct[1], c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYPL","PHYPD"))) || (comp == "phi" && is.element(x$struct[2], c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYPL","PHYPD")))) && con$vc.min > 0)
                     lb.sign <- "<"

               } else {

                  res <- try(uniroot(.profile.rma.mv, interval=c(con$vc.min, vc), tol=con$tol, maxiter=con$maxiter, obj=x, comp=comp, sigma2.pos=sigma2.pos, tau2.pos=tau2.pos, rho.pos=rho.pos, gamma2.pos=gamma2.pos, phi.pos=phi.pos, confint=TRUE, objective=objective, verbose=verbose, check.conv=TRUE)$root, silent=TRUE)

                  ### check if uniroot method converged
                  if (!inherits(res, "try-error")) {
                     vc.lb <- res
                     lb.conv <- TRUE
                  }

               }

               break

            }

            con$vc.min <- con$vc.min + epdiff

         }

         if (verbose)
            cat("\n")

         ###################################################################

         ### search for upper bound

         ### get diff value when setting component to vc.max; this value should be positive (i.e., discrepancy must be larger than critical value)
         ### if it is not, then the upper bound must be above vc.max

         epdiff <- abs(con$vc.max - vc) / con$eptries

         for (i in seq_len(con$eptries)) {

            res <- try(.profile.rma.mv(con$vc.max, obj=x, comp=comp, sigma2.pos=sigma2.pos, tau2.pos=tau2.pos, rho.pos=rho.pos, gamma2.pos=gamma2.pos, phi.pos=phi.pos, confint=TRUE, objective=objective, verbose=verbose), silent=TRUE)

            if (!inherits(res, "try-error")) {

               if (res < 0) {

                  vc.ub <- con$vc.max
                  ub.conv <- TRUE

                  if (is.element(comp, c("sigma2", "tau2", "gamma2")))
                     ub.sign <- ">"

                  if (is.element(comp, c("rho", "phi")) && con$vc.max < 1)
                     ub.sign <- ">"

                  if ((comp == "rho" && is.element(x$struct[1], c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYPL","PHYPD"))) || (comp == "phi" && is.element(x$struct[2], c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYPL","PHYPD"))))
                     ub.sign <- ">"

               } else {

                  res <- try(uniroot(.profile.rma.mv, interval=c(vc, con$vc.max), tol=con$tol, maxiter=con$maxiter, obj=x, comp=comp, sigma2.pos=sigma2.pos, tau2.pos=tau2.pos, rho.pos=rho.pos, gamma2.pos=gamma2.pos, phi.pos=phi.pos, confint=TRUE, objective=objective, verbose=verbose, check.conv=TRUE)$root, silent=TRUE)

                  ### check if uniroot method converged
                  if (!inherits(res, "try-error")) {
                     vc.ub <- res
                     ub.conv <- TRUE
                  }

               }

               break

            }

            con$vc.max <- con$vc.max - epdiff

         }

         ###################################################################

      }

      ######################################################################
      ######################################################################
      ######################################################################

      if (!lb.conv)
         warning(mstyle$warning("Cannot obtain lower bound of profile likelihood CI due to convergence problems."), call.=FALSE)

      if (!ub.conv)
         warning(mstyle$warning("Cannot obtain upper bound of profile likelihood CI due to convergence problems."), call.=FALSE)

      ######################################################################

      vc <- c(vc, vc.lb, vc.ub)

      if (is.element(comp, c("sigma2", "tau2", "gamma2"))) {
         vcsqrt <- sqrt(ifelse(vc >= 0, vc, NA))
         res.random <- rbind(vc, vcsqrt)
         if (comp == "sigma2") {
            if (length(x$sigma2) == 1L) {
               rownames(res.random) <- c("sigma^2", "sigma")
            } else {
               rownames(res.random) <- paste0(c("sigma^2", "sigma"), ".", sigma2.pos)
            }
         }
         if (comp == "tau2") {
            if (length(x$tau2) == 1L) {
               rownames(res.random) <- c("tau^2", "tau")
            } else {
               rownames(res.random) <- paste0(c("tau^2", "tau"), ".", tau2.pos)
            }
         }
         if (comp == "gamma2") {
            if (length(x$gamma2) == 1L) {
               rownames(res.random) <- c("gamma^2", "gamma")
            } else {
               rownames(res.random) <- paste0(c("gamma^2", "gamma"), ".", gamma2.pos)
            }
         }
      } else {
         res.random <- rbind(vc)
         if (comp == "rho") {
            if (length(x$rho) == 1L) {
               rownames(res.random) <- "rho"
            } else {
               rownames(res.random) <- paste0("rho.", rho.pos)
            }
         }
         if (comp == "phi") {
            if (length(x$phi) == 1L) {
               rownames(res.random) <- "phi"
            } else {
               rownames(res.random) <- paste0("phi.", rho.pos)
            }
         }
      }

      colnames(res.random) <- c("estimate", "ci.lb", "ci.ub")

   }

   #########################################################################
   #########################################################################
   #########################################################################

   if (fixed) {

      if (x$test == "t") {
         crit <- sapply(seq_along(x$ddf), function(j) if (x$ddf[j] > 0) qt(level/2, df=x$ddf[j], lower.tail=FALSE) else NA)
      } else {
         crit <- qnorm(level/2, lower.tail=FALSE)
      }

      beta  <- c(x$beta)
      ci.lb <- c(beta - crit * x$se)
      ci.ub <- c(beta + crit * x$se)

      if (is.function(transf)) {
         if (is.null(targs)) {
            beta  <- sapply(beta, transf)
            ci.lb <- sapply(ci.lb, transf)
            ci.ub <- sapply(ci.ub, transf)
         } else {
            beta  <- sapply(beta, transf, targs)
            ci.lb <- sapply(ci.lb, transf, targs)
            ci.ub <- sapply(ci.ub, transf, targs)
         }
      }

      ### make sure order of intervals is always increasing

      tmp <- .psort(ci.lb, ci.ub)
      ci.lb <- tmp[,1]
      ci.ub <- tmp[,2]

      res.fixed <- cbind(estimate=beta, ci.lb=ci.lb, ci.ub=ci.ub)
      rownames(res.fixed) <- rownames(x$beta)

   }

   #########################################################################
   #########################################################################
   #########################################################################

   res <- list()

   if (fixed)
      res$fixed <- res.fixed

   if (random)
      res$random <- res.random

   res$digits <- digits

   if (random) {
      res$ci.null <- ci.null
      res$lb.sign <- lb.sign
      res$ub.sign <- ub.sign
      #res$vc.min <- con$vc.min
   }

   if (.isTRUE(ddd$time)) {
      time.end <- proc.time()
      .print.time(unname(time.end - time.start)[3])
   }

   class(res) <- "confint.rma"
   return(res)

}
