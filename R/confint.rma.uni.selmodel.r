confint.rma.uni.selmodel <- function(object, parm, level, fixed=FALSE, tau2, delta, digits, transf, targs, verbose=FALSE, control, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(object), must="rma.uni.selmodel")

   if (!missing(parm))
      warning(mstyle$warning("Argument 'parm' (currently) ignored."), call.=FALSE)

   x <- object

   if (x$betaspec) ### TODO: consider providing CIs also for this case
      stop(mstyle$stop("Cannot obtain confidence intervals when one or more beta values were fixed."))

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

   .chkdots(ddd, c("time"))

   if (.isTRUE(ddd$time))
      time.start <- proc.time()

   ### check if user has specified one of the tau2 or delta arguments

   random <- !all(missing(tau2), missing(delta))

   if (!fixed && !random) {

      ### if both 'fixed' and 'random' are FALSE, obtain CIs for tau2 and all selection model parameters

      cl <- match.call()

      ### total number of non-fixed components

      comps <- ifelse(x$method != "FE" && !x$tau2.fix, 1, 0) + sum(!x$delta.fix)

      if (comps == 0)
         stop(mstyle$stop("No components for which a CI can be obtained."))

      res.all <- list()
      j <- 0

      if (x$method != "FE" && !x$tau2.fix) {
         j <- j + 1
         cl.vc <- cl
         cl.vc$tau2 <- 1
         cl.vc$time <- FALSE
         #cl.vc$object <- quote(x)
         if (verbose)
            cat(mstyle$verbose(paste("\nObtaining CI for tau2\n")))
         res.all[[j]] <- eval(cl.vc, envir=parent.frame())
      }

      if (any(!x$delta.fix)) {
         for (pos in seq_len(x$deltas)[!x$delta.fix]) {
            j <- j + 1
            cl.vc <- cl
            cl.vc$delta <- pos
            cl.vc$time <- FALSE
            #cl.vc$object <- quote(x)
            if (verbose)
               cat(mstyle$verbose(paste("\nObtaining CI for delta =", pos, "\n")))
            res.all[[j]] <- eval(cl.vc, envir=parent.frame())
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

      if (sum(!missing(tau2), !missing(delta)) > 1L)
         stop(mstyle$stop("Must specify only one of the 'tau2' or 'delta' arguments."))

      ### check if model actually contains (at least one) such a component and that it was actually estimated

      if (!missing(tau2) && (x$method == "FE" || x$tau2.fix))
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
         stop(mstyle$stop("Must specify the number for the 'delta' component."))

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

      delta.pos <- NA

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

      #return(list(comp=comp, vc=vc, tau2.pos=tau2.pos, delta.pos=delta.pos))

      ######################################################################

      ### set control parameters for uniroot() and possibly replace with user-defined values
      ### set vc.min and vc.max and possibly replace with user-defined values

      con <- list(tol=.Machine$double.eps^0.25, maxiter=1000, verbose=FALSE, eptries=10)

      if (comp == "tau2") {
         con$vc.min <- 0
         con$vc.max <- min(max(ifelse(vc <= .Machine$double.eps^0.5, 10, max(10, vc*100)), con$vc.min), x$tau2.max)
      }

      if (comp == "delta") {
         con$vc.min <- max(0, x$delta.min[delta])
         con$vc.max <- min(max(ifelse(vc <= .Machine$double.eps^0.5, 10, max(10, vc*10)), con$vc.min), x$delta.max[delta])
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

            res <- try(.profile.rma.uni.selmodel(con$vc.min, obj=x, comp=comp, delta.pos=delta.pos, confint=TRUE, objective=objective, verbose=verbose), silent=TRUE)

            if (!inherits(res, "try-error")) {

               if (res < 0) {

                  vc.lb <- con$vc.min
                  lb.conv <- TRUE

                  if (comp == "tau2" && con$vc.min > 0)
                     lb.sign <- "<"

                  if (comp == "delta" && con$vc.min > 0)
                     lb.sign <- "<"

               } else {

                  res <- try(uniroot(.profile.rma.uni.selmodel, interval=c(con$vc.min, vc), tol=con$tol, maxiter=con$maxiter, obj=x, comp=comp, delta.pos=delta.pos, confint=TRUE, objective=objective, verbose=verbose, check.conv=TRUE)$root, silent=TRUE)

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

            res <- try(.profile.rma.uni.selmodel(con$vc.max, obj=x, comp=comp, delta.pos=delta.pos, confint=TRUE, objective=objective, verbose=verbose), silent=TRUE)

            if (!inherits(res, "try-error")) {

               if (res < 0) {

                  vc.ub <- con$vc.max
                  ub.conv <- TRUE

                  if (comp == "tau2")
                     ub.sign <- ">"

                  if (comp == "delta")
                     ub.sign <- ">"

               } else {

                  res <- try(uniroot(.profile.rma.uni.selmodel, interval=c(vc, con$vc.max), tol=con$tol, maxiter=con$maxiter, obj=x, comp=comp, delta.pos=delta.pos, confint=TRUE, objective=objective, verbose=verbose, check.conv=TRUE)$root, silent=TRUE)

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

      if (comp == "tau2") {
         vcsqrt <- sqrt(ifelse(vc >= 0, vc, NA))
         res.random <- rbind(vc, vcsqrt)
         rownames(res.random) <- c("tau^2", "tau")
      }
      if (comp == "delta") {
         res.random <- rbind(vc)
         if (x$deltas == 1L) {
            rownames(res.random) <- "delta"
         } else {
            rownames(res.random) <- paste0("delta.", delta.pos)
         }
      }

      colnames(res.random) <- c("estimate", "ci.lb", "ci.ub")

   }

   #########################################################################
   #########################################################################
   #########################################################################

   if (fixed) {

      if (is.element(x$test, c("t"))) {
         crit <- qt(level/2, df=x$dfs, lower.tail=FALSE)
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
