regtest.rma <- function(x, model="rma", predictor="sei", ret.fit=FALSE, digits, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="rma", notav=c("robust.rma", "rma.glmm", "rma.mv", "rma.ls", "rma.uni.selmodel"))

   model <- match.arg(model, c("lm", "rma"))
   predictor <- match.arg(predictor, c("sei", "vi", "ni", "ninv", "sqrtni", "sqrtninv"))

   if (missing(digits)) {
      digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
   }

   #########################################################################

   yi <- x$yi
   vi <- x$vi
   ni <- x$ni ### may be NULL
   p  <- x$p

   if (inherits(x, "rma.mh") || inherits(x, "rma.peto")) {
      weights <- NULL
      X <- cbind(rep(1,length(yi)))
   } else {
      weights <- x$weights
      X <- x$X
   }

   if (predictor == "sei")
      X <- cbind(X, sei=sqrt(vi))

   if (predictor == "vi")
      X <- cbind(X, vi=vi)

   if (is.element(predictor, c("ni", "ninv", "sqrtni", "sqrtninv"))) {

      if (is.null(ni)) {
         stop(mstyle$stop("No sample size information stored in model object."))

      } else {

         if (predictor == "ni")
            X <- cbind(X, ni=ni)
         if (predictor == "ninv")
            X <- cbind(X, ninv=1/ni)
         if (predictor == "sqrtni")
            X <- cbind(X, ni=sqrt(ni))
         if (predictor == "sqrtninv")
            X <- cbind(X, ni=1/sqrt(ni))

      }

   }

   ### check if X of full rank (if not, cannot carry out the test)

   tmp <- lm(yi ~ X - 1)
   coef.na <- is.na(coef(tmp))
   if (any(coef.na))
      stop(mstyle$stop("Model matrix no longer of full rank after addition of predictor. Cannot fit model."))

   if (model == "rma") {

      fit  <- rma.uni(yi, vi, weights=weights, mods=X, intercept=FALSE, method=x$method, weighted=x$weighted, test=x$test, level=x$level, tau2=ifelse(x$tau2.fix, x$tau2, NA), control=x$control, ...)
      zval <- fit$zval[p+1]
      pval <- fit$pval[p+1]
      dfs  <- fit$dfs

   } else {

      yi   <- c(yi) ### to remove attributes
      fit  <- lm(yi ~ X - 1, weights=1/vi)
      tmp  <- summary(fit)
      zval <- coef(tmp)[p+1,3]
      pval <- coef(tmp)[p+1,4]
      dfs  <- x$k - x$p - 1

   }

   if (predictor %in% c("sei", "vi", "ninv", "sqrtninv") && p == 1L && .is.intercept(X[,1])) {

      if (model=="lm") {
         est <- coef(tmp)[1,1]
         ci.lb <- est - qt(x$level/2, df=dfs, lower.tail=FALSE) * coef(tmp)[1,2]
         ci.ub <- est + qt(x$level/2, df=dfs, lower.tail=FALSE) * coef(tmp)[1,2]
      } else {
         est <- coef(fit)[1]
         ci.lb <- fit$ci.lb[1]
         ci.ub <- fit$ci.ub[1]
      }

   } else {

      est <- ci.lb <- ci.ub <- NULL

   }

   res <- list(model=model, predictor=predictor, zval=zval, pval=pval, dfs=dfs, method=x$method, digits=digits, ret.fit=ret.fit, fit=fit, est=est, ci.lb=ci.lb, ci.ub=ci.ub)

   class(res) <- "regtest"
   return(res)

}
