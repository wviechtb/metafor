confint.rma.ls <- function(object, parm, level, fixed=FALSE, alpha, digits, transf, targs, verbose=FALSE, control, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(object), must="rma.ls")

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

   ### check if user has specified alpha argument

   random <- !missing(alpha)

   if (!fixed && !random) {

      ### if both 'fixed' and 'random' are FALSE, obtain CIs for alpha parameters

      cl <- match.call()

      ### total number of non-fixed components

      comps <- sum(!x$alpha.fix)

      if (comps == 0)
         stop(mstyle$stop("No components for which a CI can be obtained."))

      res.all <- list()
      j <- 0

      if (any(!x$alpha.fix)) {
         for (pos in seq_len(x$alphas)[!x$alpha.fix]) {
            j <- j + 1
            cl.vc <- cl
            cl.vc$alpha <- pos
            cl.vc$time <- FALSE
            #cl.vc$object <- quote(x)
            if (verbose)
               cat(mstyle$verbose(paste("\nObtaining CI for alpha =", pos, "\n")))
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

      ### check if model actually contains (at least one) such a component and that it was actually estimated

      if (!missing(alpha) && all(x$alpha.fix))
         stop(mstyle$stop("Model does not contain any estimated 'alpha' components."))

      ### check if user specified more than one alpha component

      if (!missing(alpha) && (length(alpha) > 1L))
         stop(mstyle$stop("Can only specify one 'alpha' component."))

      ### check if user specified a logical

      if (!missing(alpha) && is.logical(alpha))
         stop(mstyle$stop("Must specify the number for the 'alpha' component."))

      ### check if user specified a component that does not exist

      if (!missing(alpha) && (alpha > x$alphas || alpha <= 0))
         stop(mstyle$stop("No such 'alpha' component in the model."))

      ### check if user specified a component that was fixed

      if (!missing(alpha) && x$alpha.fix[alpha])
         stop(mstyle$stop("Specified 'alpha' component was fixed."))

      ### if everything is good so far, get value of the variance component and set 'comp'

      alpha.pos <- NA

      if (!missing(alpha)) {
         vc <- x$alpha[alpha]
         comp <- "alpha"
         alpha.pos <- alpha
      }

      #return(list(comp=comp, vc=vc, alpha.pos=alpha.pos))

      ######################################################################

      ### set control parameters for uniroot() and possibly replace with user-defined values
      ### set vc.min and vc.max and possibly replace with user-defined values

      con <- list(tol=.Machine$double.eps^0.25, maxiter=1000, verbose=FALSE, eptries=10)

      if (comp == "alpha") {
         if (is.na(x$se.alpha[alpha])) {
            con$vc.min <- vc/4
            con$vc.max <- vc*4
         } else {
            con$vc.min <- vc - qnorm(.995) * x$se.alpha[alpha]
            con$vc.max <- vc + qnorm(.995) * x$se.alpha[alpha]
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

            res <- try(.profile.rma.ls(con$vc.min, obj=x, comp=comp, alpha.pos=alpha.pos, confint=TRUE, objective=objective, verbose=verbose), silent=TRUE)

            if (!inherits(res, "try-error") && !is.na(res)) {

               if (res < 0) {

                  vc.lb <- con$vc.min
                  lb.conv <- TRUE
                  lb.sign <- "<"

               } else {

                  res <- try(uniroot(.profile.rma.ls, interval=c(con$vc.min, vc), tol=con$tol, maxiter=con$maxiter, obj=x, comp=comp, alpha.pos=alpha.pos, confint=TRUE, objective=objective, verbose=verbose, check.conv=TRUE)$root, silent=TRUE)

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

            res <- try(.profile.rma.ls(con$vc.max, obj=x, comp=comp, alpha.pos=alpha.pos, confint=TRUE, objective=objective, verbose=verbose), silent=TRUE)

            if (!inherits(res, "try-error") && !is.na(res)) {

               if (res < 0) {

                  vc.ub <- con$vc.max
                  ub.conv <- TRUE
                  ub.sign <- ">"

               } else {

                  res <- try(uniroot(.profile.rma.ls, interval=c(vc, con$vc.max), tol=con$tol, maxiter=con$maxiter, obj=x, comp=comp, alpha.pos=alpha.pos, confint=TRUE, objective=objective, verbose=verbose, check.conv=TRUE)$root, silent=TRUE)

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

      if (comp == "alpha") {
         res.random <- rbind(vc)
         if (x$alphas == 1L) {
            rownames(res.random) <- "alpha"
         } else {
            rownames(res.random) <- paste0("alpha.", alpha.pos)
         }
      }

      colnames(res.random) <- c("estimate", "ci.lb", "ci.ub")

   }

   #########################################################################
   #########################################################################
   #########################################################################

   if (fixed) {

      if (x$test == "t") {
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
