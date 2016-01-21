regtest.rma <- function(x, model="rma", predictor="sei", ret.fit=FALSE, ...) {

   #########################################################################

   if (!is.element("rma", class(x)))
      stop("Argument 'x' must be an object of class \"rma\".")

   if (is.element("robust.rma", class(x)))
      stop("Function not applicable to objects of class \"robust.rma\".")

   if (is.element("rma.glmm", class(x)))
      stop("Method not yet implemented for objects of class \"rma.glmm\". Sorry!")

   if (is.element("rma.mv", class(x)))
      stop("Method not yet implemented for objects of class \"rma.mv\". Sorry!")

   model <- match.arg(model, c("lm", "rma"))
   predictor <- match.arg(predictor, c("sei", "vi", "ni", "ninv", "sqrtni", "sqrtninv"))

   #########################################################################

   yi <- x$yi
   vi <- x$vi
   weights <- x$weights
   ni <- x$ni ### may be NULL
   X  <- x$X
   p  <- x$p

   if (predictor == "sei")
      X <- cbind(X, sei=sqrt(vi))

   if (predictor == "vi")
      X <- cbind(X, vi=vi)

   if (is.element(predictor, c("ni", "ninv", "sqrtni", "sqrtninv"))) {

      if (is.null(ni)) {
         stop("No sample size information stored in model object.")

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

   ### check if X is no longer of full rank (in which case we cannot carry out the test)

   tmp <- lm(yi ~ X - 1)
   coef.na <- is.na(coef(tmp))
   if (any(coef.na))
      stop("Model matrix no longer of full rank after addition of predictor. Cannot fit model.")

   if (model == "rma") {

      fit  <- rma.uni(yi, vi, weights=weights, mods=X, intercept=FALSE, method=x$method, weighted=x$weighted, knha=x$knha, tau2=ifelse(x$tau2.fix, x$tau2, NA), control=x$control, ...)
      zval <- fit$zval[p+1]
      pval <- fit$pval[p+1]
      dfs  <- fit$dfs

   } else {

      fit <- lm(yi ~ X - 1, weights=1/vi)

      fit  <- summary(fit)
      zval <- coef(fit)[p+1,3]
      pval <- coef(fit)[p+1,4]
      dfs  <- x$k - x$p - 1

   }

   res <- list(model=model, predictor=predictor, zval=zval, pval=pval, dfs=dfs, method=x$method, digits=x$digits, ret.fit=ret.fit, fit=fit)
   class(res) <- "regtest.rma"
   return(res)

}
