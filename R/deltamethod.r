deltamethod <- function(x, vcov, fun, order=1, level, H0=0, digits) {

   mstyle <- .get.mstyle()

   if (!requireNamespace("calculus", quietly=TRUE))
      stop(mstyle$stop("Please install the 'calculus' package to use this function."))

   if (missing(vcov))
      vcov <- NULL

   if (!is.function(fun))
      stop(mstyle$stop("Argument 'fun' must be a function."))

   if (!is.element(order, c(1,2)))
      stop(mstyle$stop("Argument 'order' must be equal to 1 or 2."))

   #########################################################################

   if (.is.vector(x)) {

      ### when x is a vector of coefficients

      coef <- x

      if (is.null(vcov))
         stop(mstyle$stop("Must specify the 'vcov' argument when 'x' is a vector."))

   } else {

      ### when x is not a vector (and then presumably a model object)

      coef <- try(coef(x))

      if (inherits(coef, "try-error"))
         stop(mstyle$stop("Cannot extract coefficients via coef() from 'x'."))

      if (!is.null(vcov))
         warning(mstyle$warning("Argument 'vcov' ignored when 'x' is a model object."))

      vcov <- try(vcov(x))

      if (inherits(vcov, "try-error"))
         stop(mstyle$stop("Cannot extract var-cov matrix via vcov() from 'x'."))

      if (is.list(coef) && names(coef)[1] == "beta")
         coef <- coef$beta

      if (is.list(vcov) && names(vcov)[1] == "beta")
         vcov <- vcov$beta

   }

   if (inherits(x, "rma")) {

      if (missing(digits)) {
         digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
      } else {
         digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
      }

      if (missing(level))
         level <- x$level

   } else {

      if (missing(digits))
         digits <- c(est=4, se=4, test=4, pval=4, ci=4)

      if (length(digits) == 1L)
         digits <- c(est=digits, se=digits, test=digits, pval=digits, ci=digits)

      if (missing(level))
         level <- 95

   }

   #########################################################################

   if (.is.vector(vcov) || nrow(vcov) == 1L || ncol(vcov) == 1L)
      vcov <- .diag(as.vector(vcov))

   if (!.is.square(vcov))
      stop(mstyle$stop("Argument 'vcov' must be a square matrix."))

   if (!is.null(dimnames(vcov)))
      vcov <- unname(vcov)

   if (!isSymmetric(vcov))
      stop(mstyle$stop("Argument 'vcov' must be a symmetric matrix."))

   p <- length(coef)
   pvcov <- nrow(vcov)

   if (p != pvcov)
      stop(mstyle$stop(paste0("Length of the 'coef' vector (", p, ") does not match the dimensions of 'vcov' (", pvcov, "x", pvcov, ").")))

   args <- formalArgs(fun)

   if (length(args) == 1L) {

      coef <- unname(coef)
      coef.transf <- try(fun(coef))

   } else {

      if (length(args) != p)
         stop(mstyle$stop(paste0("Number of function arguments (", length(args), ") does not match the number of coefficients (", p, ").")))

      names(coef) <- args
      coef.transf <- try(do.call(fun, args=as.list(coef)))

   }

   if (inherits(coef.transf, "try-error"))
      stop(mstyle$stop("Error when applying the function to the coefficient(s)."))

   if (!.is.vector(coef.transf))
      stop(mstyle$stop("Specified function does not return an atomic vector."))

   grad <- try(calculus::derivative(fun, var=coef, drop=FALSE))

   if (inherits(grad, "try-error"))
      stop(mstyle$stop("Error when computing the gradient."))

   if (ncol(grad) != p)
      stop(mstyle$stop(paste0("Length of the gradient (", ncol(grad), ") does not match the dimensions of 'vcov' (", pvcov, "x", pvcov, ").")))

   if (order == 2) {
      Hessian <- try(calculus::hessian(fun, var=coef, accuracy=4, drop=TRUE))
      if (inherits(Hessian, "try-error"))
         stop(mstyle$stop("Error when computing the Hessian."))
   }

   q <- length(coef.transf)

   if (length(H0) == 1L)
      H0 <- rep(H0, q)

   if (length(H0) != q)
      stop(mstyle$stop(paste0("Length of the 'H0' argument (", length(H0), ") does not match the number of transformed coefficients (", q, ").")))

   #########################################################################

   level <- .level(level)

   if (order == 1) {
      vcov.transf <- grad %*% vcov %*% t(grad)
   } else {
      vcov.transf <- grad %*% vcov %*% t(grad) + 1/2 * .tr(Hessian %*% vcov %*% vcov %*% Hessian)
   }

   rownames(vcov.transf) <- colnames(vcov.transf) <- names(coef.transf)

   crit <- qnorm(level/2, lower.tail=FALSE)
   se.transf <- sqrt(diag(vcov.transf))
   ci.lb <- coef.transf - crit * se.transf
   ci.ub <- coef.transf + crit * se.transf
   zval <- (coef.transf - H0) / se.transf
   pval <- 2*pnorm(abs(zval), lower.tail=FALSE)

   #########################################################################

   res <- list(tab = data.frame(coef=coef.transf, se=se.transf, zval=zval, pval=pval, ci.lb=ci.lb, ci.ub=ci.ub), vcov=vcov.transf, level=level, digits=digits, test="z")
   rownames(res$tab) <- names(coef.transf)
   class(res) <- "deltamethod"
   return(res)

}
