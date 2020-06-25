### Note: There is code below to obtain a profile likelihood CI for tau^2, but then
### I may have to introduce a 'type' argument to specify which type of CI to obtain.
### Actually, what would be most consistent is this:
### if method='ML/REML':    profile likelihood (PL) CI (based on the ML/REML likelihood)
### if method='EB/PM/PMM':  Q-profile (QP) CI
### if method='GENQ/GENQM': generalized Q-statistic (GENQ) CI (which also covers method='DL/HE' as special cases)
### if method='SJ':         method by Sidik & Jonkman (2005) (but this performs poorly, except if tau^2 is very large)
### if method='HS':         not sure since this is an ad-hoc estimator with no obvious underlying statistical principle
### Also could in principle compute Wald-type CIs (but those perform poorly except when k is very large).
### But it may be a bit late to change how the function works (right now, type="GENQ" if method="GENQ/GENQM" and type="QP" otherwise).

confint.rma.uni <- function(object, parm, level, fixed=FALSE, random=TRUE, digits, transf, targs, verbose=FALSE, control, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(object, "rma.uni"))
      stop(mstyle$stop("Argument 'object' must be an object of class \"rma.uni\"."))

   if (inherits(object, "robust.rma"))
      stop(mstyle$stop("Method not available for objects of class \"robust.rma\"."))

   if (inherits(object, "rma.ls"))
      stop(mstyle$stop("Method not available for objects of class \"rma.ls\"."))

   x <- object

   k <- x$k
   p <- x$p
   yi <- x$yi
   vi <- x$vi
   X <- x$X
   Y <- cbind(yi)
   weights <- x$weights

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

   if (!fixed && !random)
      stop(mstyle$stop("At least one of the arguments 'fixed' and 'random' must be TRUE."))

   if (x$method == "GENQ" || x$method == "GENQM") {
      type <- "GENQ"
   } else {
      type <- "QP"
   }

   #if (missing(type)) {
   #   if (x$method == "GENQ" || x$method == "GENQM") {
   #      type <- "GENQ"
   #   } else {
   #      type <- "QP"
   #   }
   #} else {
   #   type <- match.arg(type, c("QP", "GENQ", "PL"))
   #}

   level <- ifelse(level == 0, 1, ifelse(level >= 1, (100-level)/100, ifelse(level > .5, 1-level, level)))

   ddd <- list(...)

   if (.isTRUE(ddd$time))
      time.start <- proc.time()

   #########################################################################
   #########################################################################
   #########################################################################

   if (random) {

      if (k == 1)
         stop(mstyle$stop("Stopped because k = 1."))

      if (x$method == "FE")
         stop(mstyle$stop("Model does not contain a random-effects component."))

      if (x$tau2.fix)
         stop(mstyle$stop("Model does not contain an estimated random-effects component."))

      if (type == "GENQ" && !(is.element(x$method, c("GENQ","GENQM"))))
         stop(mstyle$stop("Model must be fitted with 'method=\"GENQ\" or 'method=\"GENQM\" to use this option."))

      ######################################################################

      ### set control parameters for uniroot() and possibly replace with user-defined values
      ### set tau2.min and tau2.max and possibly replace with user-defined values
      ### note: default tau2.min is smaller of 0 or tau2, since tau2 could in principle be negative
      ### note: default tau2.max must be larger than tau2 and tau2.min and really should be much larger (at least 100)

      tau2.min <- ifelse(is.null(x$control$tau2.min), min(0, x$tau2), x$control$tau2.min)
      tau2.max <- ifelse(is.null(x$control$tau2.max), max(100, x$tau2*10, tau2.min*10), x$control$tau2.max)

      ### user can in principle set non-sensical limits (i.e., tau2.min > tau2.max), but this is handled properly by the methods below

      con <- list(tol=.Machine$double.eps^0.25, maxiter=1000, tau2.min=tau2.min, tau2.max=tau2.max, verbose=FALSE)

      con[pmatch(names(control), names(con))] <- control

      if (verbose)
         con$verbose <- verbose

      verbose <- con$verbose

      #return(con)

      ######################################################################

      tau2.lb <- NA
      tau2.ub <- NA
      ci.null <- FALSE ### logical if CI is a null set
      lb.conv <- FALSE ### logical if search converged for lower bound (LB)
      ub.conv <- FALSE ### logical if search converged for upper bound (UB)
      lb.sign <- ""    ### for sign in case LB must be below tau2.min ("<") or above tau2.max (">")
      ub.sign <- ""    ### for sign in case UB must be below tau2.min ("<") or above tau2.max (">")

      ######################################################################

      ########################
      ### Q-profile method ###
      ########################

      if (type == "QP") {

         if (!x$allvipos)
            stop(mstyle$stop("Cannot compute CI for tau^2 when there are non-positive sampling variances in the data."))

         crit.u <- qchisq(level/2, k-p, lower.tail=FALSE) ### upper critical chi^2 value for df = k-p
         crit.l <- qchisq(level/2, k-p, lower.tail=TRUE)  ### lower critical chi^2 value for df = k-p

         QE.tau2.max <- .QE.func(con$tau2.max, Y=Y, vi=vi, X=X, k=k, objective=0)
         QE.tau2.min <- .QE.func(con$tau2.min, Y=Y, vi=vi, X=X, k=k, objective=0)

         #dfs <- 12; curve(dchisq(x, df=dfs), from=0, to=40, ylim=c(0,.1), xlab="", ylab=""); abline(v=qchisq(c(.025, .975), df=dfs)); text(qchisq(c(.025, .975), df=dfs)+1.6, .1, c("crit.l", "crit.u"))

         ###################################################################

         ### start search for upper bound

         if (QE.tau2.min < crit.l) {

            ### if QE.tau2.min is to the left of the crit.l, then both bounds are below tau2.min

            tau2.lb <- con$tau2.min
            tau2.ub <- con$tau2.min
            lb.sign <- "<"
            ub.sign <- "<"
            lb.conv <- TRUE
            ub.conv <- TRUE

            ### and if tau2.min <= 0, then the CI is equal to the null set

            if (con$tau2.min <= 0)
               ci.null <- TRUE

         } else {

            if (QE.tau2.max > crit.l) {

               ### if QE.tau2.max is to the right of crit.l, then upper bound > tau2.max, so set tau2.ub to >tau2.max

               tau2.ub <- con$tau2.max
               ub.sign <- ">"
               ub.conv <- TRUE

            } else {

               ### now QE.tau2.min is to the right of crit.l and QE.tau2.max is to the left of crit.l, so upper bound can be found

               res <- try(uniroot(.QE.func, interval=c(con$tau2.min, con$tau2.max), tol=con$tol, maxiter=con$maxiter, Y=Y, vi=vi, X=X, k=k, objective=crit.l, verbose=verbose, digits=digits)$root, silent=TRUE)

               ### check if uniroot method converged

               if (!inherits(res, "try-error")) {
                  tau2.ub <- res
                  ub.conv <- TRUE
               }

            }

         }

         ### end search for upper bound

         ###################################################################

         ### start search for lower bound

         if (QE.tau2.max > crit.u) {

            ### if QE.tau2.max is to the right of the crit.u, then both bounds are above tau2.max

            tau2.lb <- con$tau2.max
            tau2.ub <- con$tau2.max
            lb.sign <- ">"
            ub.sign <- ">"
            lb.conv <- TRUE
            ub.conv <- TRUE

         } else {

            if (QE.tau2.min < crit.u) {

               ### if QE.tau2.min is to the left of crit.u, then lower bound < tau2.min, so set tau2.lb to <tau2.min

               tau2.lb <- con$tau2.min
               lb.conv <- TRUE

               if (con$tau2.min > 0)
                  lb.sign <- "<"

            } else {

               ### now QE.tau2.min is to the right of crit.u and QE.tau2.max is to the left of crit.u, so lower bound can be found

               res <- try(uniroot(.QE.func, interval=c(con$tau2.min, con$tau2.max), tol=con$tol, maxiter=con$maxiter, Y=Y, vi=vi, X=X, k=k, objective=crit.u, verbose=verbose, digits=digits)$root, silent=TRUE)

               ### check if uniroot method converged

               if (!inherits(res, "try-error")) {
                  tau2.lb <- res
                  lb.conv <- TRUE
               }

            }

         }

         ### end search for lower bound

         ###################################################################

      }

      ######################################################################

      ###################
      ### GENQ method ###
      ###################

      if (type == "GENQ") {

         if (!requireNamespace("CompQuadForm", quietly=TRUE))
            stop(mstyle$stop("Please install the 'CompQuadForm' package when method='QGEN'."))

         A <- diag(weights, nrow=k, ncol=k)
         stXAX <- .invcalc(X=X, W=A, k=k)
         P <- A - A %*% X %*% stXAX %*% t(X) %*% A
         Q <- crossprod(Y,P) %*% Y

         ### note: .GENQ.func(tau2val, ..., Q=Q, level=0, getlower=TRUE) gives the area to the right of Q for a
         ### distribution with specified tau2val; and as we increase tau2val, so does the area to the right of Q

         GENQ.tau2.max <- .GENQ.func(con$tau2.max, P=P, vi=vi, Q=Q, level=0, k=k, p=p, getlower=TRUE)
         GENQ.tau2.min <- .GENQ.func(con$tau2.min, P=P, vi=vi, Q=Q, level=0, k=k, p=p, getlower=TRUE)

         ###################################################################

         ### start search for upper bound

         if (GENQ.tau2.min > 1 - level/2) {

            ### if GENQ.tau2.min is to the right of 1 - level/2, then both bounds are below tau2.min

            tau2.lb <- con$tau2.min
            tau2.ub <- con$tau2.min
            lb.sign <- "<"
            ub.sign <- "<"
            lb.conv <- TRUE
            ub.conv <- TRUE

            ### and if tau2.min = 0, then the CI is equal to the null set

            if (con$tau2.min <= 0)
               ci.null <- TRUE

         } else {

            if (GENQ.tau2.max < 1 - level/2) {

               ### if GENQ.tau2.max is to the left of 1 - level/2, then upper bound > tau2.max, so set tau2.ub to >tau2.max

               tau2.ub <- con$tau2.max
               ub.sign <- ">"
               ub.conv <- TRUE

            } else {

               ### now GENQ.tau2.min is to the left of 1 - level/2 and GENQ.tau2.max is to the right of 1 - level/2, so upper bound can be found

               res <- try(uniroot(.GENQ.func, c(con$tau2.min, con$tau2.max), P=P, vi=vi, Q=Q, level=level/2, k=k, p=p, getlower=FALSE, verbose=verbose, digits=digits)$root, silent=TRUE)

               ### check if uniroot method converged

               if (!inherits(res, "try-error")) {
                  tau2.ub <- res
                  ub.conv <- TRUE
               }

            }

         }

         ### end search for upper bound

         ###################################################################

         ### start search for lower bound

         if (GENQ.tau2.max < level/2) {

            ### if GENQ.tau2.max is to the left of level/2, then both bounds are abova tau2.max

            tau2.lb <- con$tau2.max
            tau2.ub <- con$tau2.max
            lb.sign <- ">"
            ub.sign <- ">"
            lb.conv <- TRUE
            ub.conv <- TRUE

         } else {

            if (GENQ.tau2.min > level/2) {

               ### if GENQ.tau2.min is to the right of level/2, then lower bound < tau2.min, so set tau2.lb to <tau2.min

               tau2.lb <- con$tau2.min
               lb.conv <- TRUE

               if (con$tau2.min > 0)
                  lb.sign <- "<"

            } else {

               ### now GENQ.tau2.max is to the right of level/2 and GENQ.tau2.min is to the left of level/2, so lower bound can be found

               res <- try(uniroot(.GENQ.func, c(con$tau2.min, con$tau2.max), P=P, vi=vi, Q=Q, level=level/2, k=k, p=p, getlower=TRUE, verbose=verbose, digits=digits)$root, silent=TRUE)

               ### check if uniroot method converged

               if (!inherits(res, "try-error")) {
                  tau2.lb <- res
                  lb.conv <- TRUE
               }

            }

         }

         ### end search for lower bound

         ###################################################################

      }

      ######################################################################

      #################
      ### PL method ###
      #################

      ### note: cannot actually use this at the moment

      if (type == "PL") {

         if (con$tau2.min > x$tau2)
            stop(mstyle$stop("Lower bound of interval to be searched must be <= actual value of component."))
         if (con$tau2.max < x$tau2)
            stop(mstyle$stop("Upper bound of interval to be searched must be >= actual value of component."))

         objective <- qchisq(1-level, df=1)

         ###################################################################

         ### start search for lower bound

         ### get diff value when setting component to tau2.min; this value should be positive (i.e., discrepancy must be larger than critical value)
         ### if it is not, then the lower bound must be below tau2.min

         res <- try(.profile.rma.uni(val = con$tau2.min, obj=x, CI=TRUE, objective=objective, verbose=verbose), silent=TRUE)

         if (!inherits(res, "try-error")) {

            if (res < 0) {

               tau2.lb <- con$tau2.min
               lb.conv <- TRUE

               if (con$tau2.min > 0)
                  lb.sign <- "<"

            } else {

               res <- try(uniroot(.profile.rma.uni, interval=c(con$tau2.min, x$tau2), tol=con$tol, maxiter=con$maxiter, obj=x, CI=TRUE, objective=objective, verbose=verbose, check.conv=TRUE)$root, silent=TRUE)

               ### check if uniroot method converged

               if (!inherits(res, "try-error")) {
                  tau2.lb <- res
                  lb.conv <- TRUE
               }

            }

         }

         ### end search for lower bound

         ###################################################################

         ### start search for upper bound

         ### get diff value when setting component to tau2.max; this value should be positive (i.e., discrepancy must be larger than critical value)
         ### if it is not, then the upper bound must be above tau2.max

         res <- try(.profile.rma.uni(val = con$tau2.max, obj=x, CI=TRUE, objective=objective, verbose=verbose), silent=TRUE)

         if (!inherits(res, "try-error")) {

            if (res < 0) {

               tau2.ub <- con$tau2.max
               ub.conv <- TRUE
               ub.sign <- ">"

            } else {

               res <- try(uniroot(.profile.rma.uni, interval=c(x$tau2, con$tau2.max), tol=con$tol, maxiter=con$maxiter, obj=x, CI=TRUE, objective=objective, verbose=verbose, check.conv=TRUE)$root, silent=TRUE)

               ### check if uniroot method converged

               if (!inherits(res, "try-error")) {
                  tau2.ub <- res
                  ub.conv <- TRUE
               }

            }

         }

         ### end search for upper bound

         ###################################################################

      }

      ######################################################################

      if (!lb.conv)
         warning(mstyle$warning("Error in iterative search for the lower bound."), call.=FALSE)

      if (!ub.conv)
         warning(mstyle$warning("Error in iterative search for the upper bound."), call.=FALSE)

      #if (lb.sign == "<" && con$tau2.min > 0)
      #   warning(mstyle$warning("Lower bound < tau2.min. Try decreasing tau2.min (via the 'control' argument)."), call.=FALSE)
      #if (ub.sign == ">")
      #   warning(mstyle$warning("Upper bound > tau2.max. Try increasing tau2.max (via the 'control' argument)."), call.=FALSE)

      ######################################################################

      I2.lb <- 100 * tau2.lb / (x$vt + tau2.lb)
      I2.ub <- 100 * tau2.ub / (x$vt + tau2.ub)
      H2.lb <- tau2.lb / x$vt + 1
      H2.ub <- tau2.ub / x$vt + 1

      tau2 <- c(x$tau2, tau2.lb, tau2.ub)
      tau  <- sqrt(c(ifelse(x$tau2 >= 0, x$tau2, NA), ifelse(tau2.lb >= 0, tau2.lb, NA), ifelse(tau2.ub >= 0, tau2.ub, NA)))
      I2   <- c(x$I2, I2.lb, I2.ub)
      H2   <- c(x$H2, H2.lb, H2.ub)

      res.random <- rbind("tau^2"=tau2, "tau"=tau, "I^2(%)"=I2, "H^2"=H2)
      colnames(res.random) <- c("estimate", "ci.lb", "ci.ub")

   }

   #########################################################################
   #########################################################################
   #########################################################################

   if (fixed) {

      if (is.element(x$test, c("knha","adhoc","t"))) {
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
      res$tau2.min <- con$tau2.min
   }

   if (.isTRUE(ddd$time)) {
      time.end <- proc.time()
      .print.time(unname(time.end - time.start)[3])
   }

   class(res) <- "confint.rma"
   return(res)

}
