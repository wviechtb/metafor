regtest.default <- function(x, vi, sei, ni, subset, model="rma", predictor="sei", ret.fit=FALSE, digits, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   if (missing(subset))
      subset <- NULL

   ### set defaults for digits

   if (missing(digits)) {
      digits <- .set.digits(dmiss=TRUE)
   } else {
      digits <- .set.digits(digits, dmiss=FALSE)
   }

   if (missing(ni))
      ni <- NULL

   model <- match.arg(model, c("lm", "rma"))
   predictor <- match.arg(predictor, c("sei", "vi", "ni", "ninv", "sqrtni", "sqrtninv"))

   ddd <- list(...)

   if (!is.null(ddd$level)) {
      level <- ifelse(level == 0, 1, ifelse(level >= 1, (100-level)/100, ifelse(level > .5, 1-level, level)))
   } else {
      level <- .05
   }

   #########################################################################

   ### check if sampling variances and/or standard errors are available

   if (missing(vi))
      vi <- NULL

   if (missing(sei))
      sei <- NULL

   if (is.null(vi)) {
      if (!is.null(sei))
         vi <- sei^2
   }

   if (is.null(vi))
      stop(mstyle$stop("Must specify 'vi' or 'sei' argument."))

   yi <- x

   ### check length of yi and vi

   if (length(yi) != length(vi))
      stop(mstyle$stop("Length of 'yi' and 'vi' (or 'sei') is not the same."))

   ### if ni has not been specified but is an attribute of yi, get it

   if (is.null(ni) && !is.null(attr(yi, "ni")))
      ni <- attr(yi, "ni")

   ### check length of yi and ni (only if ni is not NULL)
   ### if there is a mismatch, then ni cannot be trusted, so set it to NULL

   if (!is.null(ni) && length(ni) != length(yi))
      ni <- NULL

   ### if ni is now available, add it (back) as an attribute to yi

   if (!is.null(ni))
      attr(yi, "ni") <- ni

   ### if a subset of studies is specified

   if (!is.null(subset)) {
      subset <- .setnafalse(subset, k=length(yi))
      yi <- yi[subset]
      vi <- vi[subset]
      ni <- ni[subset]
   }

   ### check for NAs and act accordingly

   has.na <- is.na(yi) | is.na(vi) | (if (is.null(ni)) FALSE else is.na(ni))

   if (any(has.na)) {

      not.na <- !has.na

      if (na.act == "na.omit" || na.act == "na.exclude" || na.act == "na.pass") {

         yi <- yi[not.na]
         vi <- vi[not.na]
         ni <- ni[not.na]
         warning(mstyle$warning("Studies with NAs omitted from test."), call.=FALSE)

      }

      if (na.act == "na.fail")
         stop(mstyle$stop("Missing values in data."))

   }

   #########################################################################

   if (predictor == "sei")
      X <- cbind(1, sei=sqrt(vi))

   if (predictor == "vi")
      X <- cbind(1, vi=vi)

   if (is.element(predictor, c("ni", "ninv", "sqrtni", "sqrtninv"))) {

      if (is.null(ni)) {
         stop(mstyle$stop("Sample size information must be specified via 'ni' argument."))

      } else {

         if (predictor == "ni")
            X <- cbind(1, ni=ni)
         if (predictor == "ninv")
            X <- cbind(1, ninv=1/ni)
         if (predictor == "sqrtni")
            X <- cbind(1, ni=sqrt(ni))
         if (predictor == "sqrtninv")
            X <- cbind(1, ni=1/sqrt(ni))

      }

   }

   ### check if X of full rank (if not, cannot carry out the test)

   tmp <- lm(yi ~ X - 1)
   coef.na <- is.na(coef(tmp))
   if (any(coef.na))
      stop(mstyle$stop("Model matrix not of full. Cannot fit model."))

   if (model == "rma") {

      fit  <- rma.uni(yi, vi, mods=X, intercept=FALSE, ...)
      zval <- fit$zval[2]
      pval <- fit$pval[2]
      dfs  <- fit$dfs

   } else {

      yi   <- c(yi) ### to remove attributes
      fit  <- lm(yi ~ X - 1, weights=1/vi)
      tmp  <- summary(fit)
      zval <- coef(tmp)[2,3]
      pval <- coef(tmp)[2,4]
      dfs  <- length(yi) - 2

   }

   if (predictor %in% c("sei", "vi", "ninv", "sqrtninv")) {

      if (model=="lm") {
         est <- coef(tmp)[1,1]
         ci.lb <- est - qt(level/2, df=dfs, lower.tail=FALSE) * coef(tmp)[1,2]
         ci.ub <- est + qt(level/2, df=dfs, lower.tail=FALSE) * coef(tmp)[1,2]
      } else {
         est <- coef(fit)[1]
         ci.lb <- fit$ci.lb[1]
         ci.ub <- fit$ci.ub[1]
      }

   } else {

      est <- ci.lb <- ci.ub <- NULL

   }

   res <- list(model=model, predictor=predictor, zval=zval, pval=pval, dfs=dfs, method=fit$method, digits=digits, ret.fit=ret.fit, fit=fit, est=est, ci.lb=ci.lb, ci.ub=ci.ub)

   class(res) <- "regtest"
   return(res)

}
