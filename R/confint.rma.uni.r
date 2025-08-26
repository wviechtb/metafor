# What would be most consistent is this:
# if method='ML/REML':    profile likelihood (PL) CI (based on the ML/REML likelihood)
# if method='EB/PM/PMM':  Q-profile (QP) CI
# if method='GENQ/GENQM': generalized Q-statistic (GENQ) CI (which also covers method='DL/HE' as special cases)
# if method='SJ':         method by Sidik & Jonkman (2005) (but this performs poorly, except if tau^2 is very large)
# if method='HS':         not sure since this is an ad-hoc estimator with no obvious underlying statistical principle
# Also can compute Wald-type CIs (but those perform poorly except when k is very large).
# Too late to change how the function works (right now, type="GENQ" if method="GENQ/GENQM" and type="QP" otherwise).

confint.rma.uni <- function(object, parm, level, fixed=FALSE, random=TRUE, type, digits, transf, targs, verbose=FALSE, control, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(object), must="rma.uni", notav=c("robust.rma", "rma.ls", "rma.gen"))

   if (!missing(parm))
      warning(mstyle$warning("Argument 'parm' (currently) ignored."), call.=FALSE)

   x <- object

   k  <- x$k
   p  <- x$p
   yi <- x$yi
   vi <- x$vi
   X  <- x$X
   Y  <- cbind(yi)
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

   funlist <- lapply(list(transf.exp.int, transf.ilogit.int, transf.iprobit.int, transf.ztor.int, transf.iarcsin.int, transf.iahw.int, transf.iabt.int, transf.exp.mode, transf.ilogit.mode, transf.iprobit.mode, transf.ztor.mode, transf.iarcsin.mode, transf.iahw.mode, transf.iabt.mode), deparse)

   if (is.null(targs) && any(sapply(funlist, identical, deparse(transf))) && inherits(x, c("rma.uni","rma.glmm")) && length(x$tau2 == 1L))
      targs <- list(tau2=x$tau2)

   if (missing(control))
      control <- list()

   if (!fixed && !random)
      stop(mstyle$stop("At least one of the arguments 'fixed' and 'random' must be TRUE."))

   ddd <- list(...)

   .chkdots(ddd, c("time", "xlim", "extint"))

   if (.isTRUE(ddd$time))
      time.start <- proc.time()

   if (!is.null(ddd$xlim)) {
      if (length(ddd$xlim) == 1L)
         ddd$xlim <- c(0, ddd$xlim)
      if (length(ddd$xlim) != 2L)
         stop(mstyle$stop("Argument 'xlim' should be a vector of length 1 or 2."))
      control$tau2.min <- ddd$xlim[1]
      control$tau2.max <- ddd$xlim[2]
   }

   if (missing(type)) {
      if (x$method == "GENQ" || x$method == "GENQM") {
         type <- "genq"
      } else {
         type <- "qp"
      }
   } else {
      type <- tolower(type)
      if (!is.element(type, c("qp","genq","pl","ht","wald","wald.log","wald.sqrt")))
         stop(mstyle$stop("Unknown 'type' specified."))
   }

   level <- .level(level, stopon100=(type=="pl" && .isTRUE(ddd$extint)))

   #########################################################################
   #########################################################################
   #########################################################################

   if (random) {

      if (k == 1L)
         stop(mstyle$stop("Stopped because k = 1."))

      if (is.element(x$method, c("FE","EE","CE")))
         stop(mstyle$stop("Model does not contain a random-effects component."))

      if (x$tau2.fix)
         stop(mstyle$stop("Model does not contain an estimated random-effects component."))

      if (type == "genq" && !(is.element(x$method, c("GENQ","GENQM"))))
         stop(mstyle$stop("Model must be fitted with method=\"GENQ\" or method=\"GENQM\" to use this option."))

      ######################################################################

      ### set defaults for control parameters for uniroot() and replace with any user-defined values
      ### set tau2.min and tau2.max and replace with any user-defined values
      ### note: the default for tau2.min is the smaller of 0 and tau2, since tau2 could in principle be negative
      ### note: the default for tau2.max must be larger than tau2 and tau2.min and really should be much larger (at least 100)

      if (!is.null(x$control$tau2.min) && x$control$tau2.min == -min(x$vi))
         x$control$tau2.min <- x$control$tau2.min + 0.0001 # push tau2.min just a bit above -min(vi) to avoid division by zero

      tau2.min <- ifelse(is.null(x$control$tau2.min), min(0, x$tau2), x$control$tau2.min)
      tau2.max <- ifelse(is.null(x$control$tau2.max), max(100, x$tau2*10, tau2.min*10), x$control$tau2.max)

      ### user can in principle set non-sensical limits (i.e., tau2.min > tau2.max), but this is handled properly by the methods below

      con <- list(tol=.Machine$double.eps^0.25, maxiter=1000, tau2.min=tau2.min, tau2.max=tau2.max, verbose=FALSE)

      con.pos <- pmatch(names(control), names(con))
      con[c(na.omit(con.pos))] <- control[!is.na(con.pos)]

      if (verbose)
         con$verbose <- verbose

      verbose <- con$verbose

      #return(con)

      ######################################################################

      tau2.lb <- NA_real_
      tau2.ub <- NA_real_
      ci.null <- FALSE # logical if CI is a null set
      lb.conv <- FALSE # logical if search converged for lower bound (LB)
      ub.conv <- FALSE # logical if search converged for upper bound (UB)
      lb.sign <- ""    # for sign in case LB must be below tau2.min ("<") or above tau2.max (">")
      ub.sign <- ""    # for sign in case UB must be below tau2.min ("<") or above tau2.max (">")

      ######################################################################

      ########################
      ### Q-profile method ###
      ########################

      if (type == "qp") {

         if (!x$allvipos)
            stop(mstyle$stop("Cannot compute CI for tau^2 when there are non-positive sampling variances in the data."))

         crit.u <- qchisq(level/2, k-p, lower.tail=FALSE) # upper critical chi^2 value for df = k-p
         crit.l <- qchisq(level/2, k-p, lower.tail=TRUE)  # lower critical chi^2 value for df = k-p

         QE.tau2.max <- .QE.func(con$tau2.max, Y=Y, vi=vi, X=X, k=k, objective=0)
         QE.tau2.min <- try(.QE.func(con$tau2.min, Y=Y, vi=vi, X=X, k=k, objective=0), silent=TRUE)

         #dfs <- 12; curve(dchisq(x, df=dfs), from=0, to=40, ylim=c(0,0.1), xlab="", ylab=""); abline(v=qchisq(c(0.025, 0.975), df=dfs)); text(qchisq(c(0.025, 0.975), df=dfs)+1.6, 0.1, c("crit.l", "crit.u"))

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

      if (type == "genq") {

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

      if (type == "pl") {

         if (con$tau2.min > x$tau2)
            stop(mstyle$stop("Lower bound of interval to be searched must be <= actual value of component."))
         if (con$tau2.max < x$tau2)
            stop(mstyle$stop("Upper bound of interval to be searched must be >= actual value of component."))

         objective <- qchisq(1-level, df=1)

         ###################################################################

         ### start search for lower bound

         ### get diff value when setting component to tau2.min; this value should be positive (i.e., discrepancy must be larger than critical value)
         ### if it is not, then the lower bound must be below tau2.min

         res <- try(.profile.rma.uni(con$tau2.min, obj=x, confint=TRUE, objective=objective, verbose=verbose), silent=TRUE)

         if (!inherits(res, "try-error") && !is.na(res)) {

            if (res < 0) {

               tau2.lb <- con$tau2.min
               lb.conv <- TRUE

               if (con$tau2.min > 0)
                  lb.sign <- "<"

            } else {

               if (.isTRUE(ddd$extint)) {
                  res <- try(uniroot(.profile.rma.uni, interval=c(con$tau2.min, x$tau2), tol=con$tol, maxiter=con$maxiter, extendInt="downX", obj=x, confint=TRUE, objective=objective, verbose=verbose, check.conv=TRUE)$root, silent=TRUE)
               } else {
                  res <- try(uniroot(.profile.rma.uni, interval=c(con$tau2.min, x$tau2), tol=con$tol, maxiter=con$maxiter, obj=x, confint=TRUE, objective=objective, verbose=verbose, check.conv=TRUE)$root, silent=TRUE)
               }

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

         res <- try(.profile.rma.uni(con$tau2.max, obj=x, confint=TRUE, objective=objective, verbose=verbose), silent=TRUE)

         if (!inherits(res, "try-error") && !is.na(res)) {

            if (!.isTRUE(ddd$extint) && res < 0) {

               tau2.ub <- con$tau2.max
               ub.conv <- TRUE
               ub.sign <- ">"

            } else {

               if (.isTRUE(ddd$extint)) {
                  res <- try(uniroot(.profile.rma.uni, interval=c(x$tau2, con$tau2.max), tol=con$tol, maxiter=con$maxiter, extendInt="upX", obj=x, confint=TRUE, objective=objective, verbose=verbose, check.conv=TRUE)$root, silent=TRUE)
               } else {
                  res <- try(uniroot(.profile.rma.uni, interval=c(x$tau2, con$tau2.max), tol=con$tol, maxiter=con$maxiter, obj=x, confint=TRUE, objective=objective, verbose=verbose, check.conv=TRUE)$root, silent=TRUE)
               }

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

      #################
      ### HT method ###
      #################

      if (type == "ht") {

         if (!x$int.only)
            stop(mstyle$stop("Method only applicable to models without moderators."))
         #if (x$method != "DL")
         #   stop(mstyle$stop("Method only applicable when 'method=DL'."))
         if (x$k <= 2)
            stop(mstyle$stop("Method only applicable when k > 2."))

         if (x$QE > x$k) {
            se.lnH <- 1/2 * (log(x$QE) - log(x$k-1)) / (sqrt(2*x$QE) - sqrt(2*x$k-3))
         } else {
            se.lnH <- sqrt(1 / (2*(x$k-2)) * (1 - 1/(3*(x$k-2)^2))) # as in Higgins and Thompson (2002), p. 1549
            #se.lnH <- sqrt(1 / ((2*(x$k-2)) * (1 - 1/(3*(x$k-2)^2)))) # as in Borenstein et al. (2009), eq. 16.21
         }

         crit <- qnorm(level/2, lower.tail=FALSE)
         lb.conv <- TRUE
         ub.conv <- TRUE
         #H2.lb <- exp(log(sqrt(x$H2)) - crit * se.lnH)^2
         #H2.ub <- exp(log(sqrt(x$H2)) + crit * se.lnH)^2
         H2.lb <- exp(log(x$H2) - crit * 2*se.lnH) # note: SE[log(H^2)] = 2*SE[log(H)]
         H2.ub <- exp(log(x$H2) + crit * 2*se.lnH)
         I2.lb <- (H2.lb - 1) / H2.lb
         I2.ub <- (H2.ub - 1) / H2.ub
         tau2.lb <- max(0, I2.lb * x$vt / (1 - I2.lb))
         tau2.ub <- I2.ub * x$vt / (1 - I2.ub)

      }

      ######################################################################

      if (is.element(type, c("wald","wald.log","wald.sqrt"))) {
         crit <- qnorm(level/2, lower.tail=FALSE)
         lb.conv <- TRUE
         ub.conv <- TRUE
      }

      ###################
      ### Wald method ###
      ###################

      if (type == "wald") {

         tau2.lb <- x$tau2 - crit * x$se.tau2
         tau2.ub <- x$tau2 + crit * x$se.tau2
         tau2.lb <- max(ifelse(is.null(x$control$tau2.min), 0, x$control$tau2.min), tau2.lb)

      }

      #######################
      ### Wald.log method ###
      #######################

      if (type == "wald.log") {

         if (x$tau2 >= 0) {
            tau2.lb <- exp(log(x$tau2) - crit * x$se.tau2 / x$tau2)
            tau2.ub <- exp(log(x$tau2) + crit * x$se.tau2 / x$tau2)
            tau2.ub <- max(x$tau2, tau2.ub) # if tau2 is 0, then CI is 0 to tau2
         }

      }

      ########################
      ### Wald.sqrt method ###
      ########################

      if (type == "wald.sqrt") {

         if (x$tau2 >= 0) {
            tau2.lb <- (max(0, sqrt(x$tau2) - crit * x$se.tau2 / (2 * sqrt(x$tau2))))^2
            tau2.ub <- (sqrt(x$tau2) + crit * x$se.tau2 / (2 * sqrt(x$tau2)))^2
         }

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
      tau  <- sqrt(c(ifelse(x$tau2 >= 0, x$tau2, NA_real_), ifelse(tau2.lb >= 0, tau2.lb, NA_real_), ifelse(tau2.ub >= 0, tau2.ub, NA_real_)))
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
         crit <- qt(level/2, df=x$ddf, lower.tail=FALSE)
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
            if (!is.primitive(transf) && !is.null(targs) && length(formals(transf)) == 1L)
               stop(mstyle$stop("Function specified via 'transf' does not appear to have an argument for 'targs'."))
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
      res$tau2.min <- con$tau2.min
   }

   if (.isTRUE(ddd$time)) {
      time.end <- proc.time()
      .print.time(unname(time.end - time.start)[3])
   }

   class(res) <- "confint.rma"
   return(res)

}
