hc.rma.uni <- function(object, digits, transf, targs, control, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(object, "rma.uni"))
      stop(mstyle$stop("Argument 'object' must be an object of class \"rma.uni\"."))

   if (inherits(object, "rma.ls"))
      stop(mstyle$stop("Method not available for objects of class \"rma.ls\"."))

   x <- object

   if (!x$int.only)
      stop(mstyle$stop("Method only applicable for models without moderators."))

   if (missing(digits))
      digits <- x$digits

   if (missing(transf))
      transf <- FALSE

   if (missing(targs))
      targs <- NULL

   yi <- x$yi
   vi <- x$vi
   k  <- length(yi)

   if (k == 1)
      stop(mstyle$stop("Stopped because k = 1."))

   if (!x$allvipos)
      stop(mstyle$stop("Cannot use method when one or more sampling variances are non-positive."))

   level <- ifelse(x$level == 0, 1, ifelse(x$level > 1, (100-x$level)/100, ifelse(x$level > .5, 1-x$level, x$level)))

   if (missing(control))
      control <- list()

   #########################################################################

   ### set control parameters for uniroot() and possibly replace with user-defined values
   con <- list(tol=.Machine$double.eps^0.25, maxiter=1000, verbose=FALSE)
   con[pmatch(names(control), names(con))] <- control

   #########################################################################

   ### original code by Henmi & Copas (2012), modified by Michael Dewey, small adjustments
   ### for consistency with other functions in the metafor package by Wolfgang Viechtbauer

   wi <- 1/vi ### fixed effects weights

   W1 <- sum(wi)
   W2 <- sum(wi^2) / W1
   W3 <- sum(wi^3) / W1
   W4 <- sum(wi^4) / W1

   ### fixed-effects estimate of theta
   beta <- sum(wi*yi) / W1

   ### Q statistic
   Q <- sum(wi * ((yi - beta)^2))

   ### DL estimate of tau^2
   tau2 <- max(0, (Q - (k-1)) / (W1 - W2))

   vb  <- (tau2 * W2 + 1) / W1 ### estimated Var of b
   se  <- sqrt(vb)             ### estimated SE of b
   VR  <- 1 + tau2 * W2        ### estimated Var of R
   SDR <- sqrt(VR)             ### estimated SD of R

   ### conditional mean of Q given R=r
   EQ <- function(r)
      (k - 1) + tau2 * (W1 - W2) + (tau2^2)*((1/VR^2) * (r^2) - 1/VR) * (W3 - W2^2)

   ### conditional variance of Q given R=r
   VQ <- function(r) {
      rsq <- r^2
      recipvr2 <- 1 / VR^2
      2 * (k - 1) + 4 * tau2 * (W1 - W2) +
      2 * tau2^2 * (W1*W2 - 2*W3 + W2^2) +
      4 * tau2^2 * (recipvr2 * rsq - 1/VR) * (W3 - W2^2) +
      4 * tau2^3 * (recipvr2 * rsq - 1/VR) * (W4 - 2*W2*W3 + W2^3) +
      2 * tau2^4 * (recipvr2 - 2 * (1/VR^3) * rsq) * (W3 - W2^2)^2
   }

   scale <- function(r){VQ(r)/EQ(r)}   ### scale parameter of the gamma distribution
   shape <- function(r){EQ(r)^2/VQ(r)} ### shape parameter of the gamma distribution

   ### inverse of f
   finv <- function(f)
      (W1/W2 - 1) * ((f^2) - 1) + (k - 1)

   ### equation to be solved
   eqn <- function(x) {
      integrand <- function(r) {
         pgamma(finv(r/x), scale=scale(SDR*r), shape=shape(SDR*r))*dnorm(r)
      }
      integral <- integrate(integrand, lower=x, upper=Inf)$value
      val <- integral - level / 2
      #cat(val, "\n")
      val
   }

   t0 <- try(uniroot(eqn, lower=0, upper=2, tol=con$tol, maxiter=con$maxiter))

   if (inherits(t0, "try-error"))
      stop(mstyle$stop("Error in uniroot()."))

   t0 <- t0$root
   u0 <- SDR * t0 ### (approximate) percentage point for the distribution of U

   #########################################################################

   ci.lb <- beta - u0 * se ### lower CI bound
   ci.ub <- beta + u0 * se ### upper CI bound

   beta.rma  <- x$beta
   se.rma    <- x$se
   ci.lb.rma <- x$ci.lb
   ci.ub.rma <- x$ci.ub

   ### if requested, apply transformation to yi's and CI bounds

   if (is.function(transf)) {
      if (is.null(targs)) {
         beta      <- sapply(beta, transf)
         beta.rma  <- sapply(beta.rma, transf)
         se        <- NA
         se.rma    <- NA
         ci.lb     <- sapply(ci.lb, transf)
         ci.ub     <- sapply(ci.ub, transf)
         ci.lb.rma <- sapply(ci.lb.rma, transf)
         ci.ub.rma <- sapply(ci.ub.rma, transf)
      } else {
         beta      <- sapply(beta, transf, targs)
         beta.rma  <- sapply(beta.rma, transf, targs)
         se        <- NA
         se.rma    <- NA
         ci.lb     <- sapply(ci.lb, transf, targs)
         ci.ub     <- sapply(ci.ub, transf, targs)
         ci.lb.rma <- sapply(ci.lb.rma, transf, targs)
         ci.ub.rma <- sapply(ci.ub.rma, transf, targs)
      }
   }

   ### make sure order of intervals is always increasing

   tmp <- .psort(ci.lb, ci.ub)
   ci.lb <- tmp[,1]
   ci.ub <- tmp[,2]

   tmp <- .psort(ci.lb.rma, ci.ub.rma)
   ci.lb.rma <- tmp[,1]
   ci.ub.rma <- tmp[,2]

   #########################################################################

   res <- list(beta=beta, se=se, ci.lb=ci.lb, ci.ub=ci.ub,
               beta.rma=beta.rma, se.rma=se.rma, ci.lb.rma=ci.lb.rma, ci.ub.rma=ci.ub.rma,
               method="DL", method.rma=x$method, tau2=tau2, tau2.rma=x$tau2, digits=digits)

   class(res) <- "hc.rma.uni"
   return(res)

}
