coef.matreg <- function(object, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(object), must="matreg")

   coefs <- c(object$tab$beta)
   names(coefs) <- rownames(object$tab)

   return(coefs)

}

vcov.matreg <- function(object, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(object), must="matreg")

   out <- object$vb
   return(out)

}

sigma.matreg <- function(object, REML=TRUE, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(object), must="matreg")

   if (object$test == "t") {

      if (REML) {
         sigma <- sqrt(object$sigma2.reml)
      } else {
         sigma <- sqrt(object$sigma2.ml)
      }

   } else {

      warning(mstyle$warning("Model does not contain a 'sigma' estimate."), call.=FALSE)
      return(invisible())

   }

   return(sigma)

}

confint.matreg <- function(object, parm, level, digits, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(object), must="matreg")

   if (!missing(parm))
      warning(mstyle$warning("Argument 'parm' (currently) ignored."), call.=FALSE)

   x <- object

   if (missing(level))
      level <- x$level

   if (missing(digits)) {
      digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
   }

   level <- .level(level)

   if (x$test=="t") {
      crit <- qt(level/2, df=x$df.residual, lower.tail=FALSE)
   } else {
      crit <- qnorm(level/2, lower.tail=FALSE)
   }

   beta  <- x$tab$beta
   ci.lb <- beta - crit * x$tab$se
   ci.ub <- beta + crit * x$tab$se

   res <- cbind(estimate=beta, ci.lb, ci.ub)
   res <- list(tab=res)
   rownames(res$tab) <- rownames(x$tab)

   res$digits <- digits

   class(res) <- "confint.matreg"
   return(res)

}

print.confint.matreg <- function(x, digits=x$digits, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(x), must="confint.matreg")

   digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)

   .space()

   tab <- cbind(fmtx(x$tab[,1,drop=FALSE], digits[["est"]]), fmtx(x$tab[,2:3,drop=FALSE], digits[["ci"]]))
   tab <- capture.output(print(tab, quote=FALSE, right=TRUE))
   .print.table(tab, mstyle)

   .space()

   invisible()

}

logLik.matreg <- function(object, REML=FALSE, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(object), must="matreg")

   if (object$test=="t") {

      if (REML) {
         val <- object$fit.stats["ll","REML"]
      } else {
         val <- object$fit.stats["ll","ML"]
      }

      attr(val, "nall") <- object$n
      attr(val, "nobs") <- ifelse(REML, object$df.residual, object$n)
      attr(val, "df")   <- object$parms

      class(val) <- "logLik"
      return(val)

   } else {

      warning(mstyle$warning("Cannot compute log-likelihood for this type of model object."), call.=FALSE)
      return(invisible())

   }

}

AIC.matreg <- function(object, ..., k=2, correct=FALSE, REML=FALSE) {

   mstyle <- .get.mstyle()

   .chkclass(class(object), must="matreg")

   if (object$test=="t") {

      if (missing(...)) {

         ### if there is just 'object'

         if (REML) {
            out <- ifelse(correct, object$fit.stats["AICc","REML"], object$fit.stats["AIC","REML"])
         } else {
            out <- ifelse(correct, object$fit.stats["AICc","ML"], object$fit.stats["AIC","ML"])
         }

      } else {

         ### if there is 'object' and additional objects via ...

         if (REML) {
            out <- sapply(list(object, ...), function(x) ifelse(correct, x$fit.stats["AICc","REML"], x$fit.stats["AIC","REML"]))
         } else {
            out <- sapply(list(object, ...), function(x) ifelse(correct, x$fit.stats["AICc","ML"], x$fit.stats["AIC","ML"]))
         }
         dfs <- sapply(list(object, ...), function(x) x$parms)

         out <- data.frame(df=dfs, AIC=out)

         if (correct)
            names(out)[2] <- "AICc"

         ### get names of objects; same idea as in stats:::AIC.default

         cl <- match.call()
         cl$k <- NULL
         cl$correct <- NULL
         cl$REML <- NULL
         rownames(out) <- as.character(cl[-1L])

      }

      return(out)

   } else {

      warning(mstyle$warning("Cannot compute AIC for this type of model object."), call.=FALSE)
      return(invisible())

   }

}

BIC.matreg <- function(object, ..., REML=FALSE) {

   mstyle <- .get.mstyle()

   .chkclass(class(object), must="matreg")

   if (object$test=="t") {

      if (missing(...)) {

         ### if there is just 'object'

         if (REML) {
            out <- object$fit.stats["BIC","REML"]
         } else {
            out <- object$fit.stats["BIC","ML"]
         }

      } else {

         ### if there is 'object' and additional objects via ...

         if (REML) {
            out <- sapply(list(object, ...), function(x) x$fit.stats["BIC","REML"])
         } else {
            out <- sapply(list(object, ...), function(x) x$fit.stats["BIC","ML"])
         }
         dfs <- sapply(list(object, ...), function(x) x$parms)

         out <- data.frame(df=dfs, BIC=out)

         ### get names of objects; same idea as in stats:::AIC.default

         cl <- match.call()
         cl$REML <- NULL
         rownames(out) <- as.character(cl[-1L])

      }

      return(out)

   } else {

      warning(mstyle$warning("Cannot compute BIC for this type of model object."), call.=FALSE)
      return(invisible())

   }

}
