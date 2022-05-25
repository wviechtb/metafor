regtest <- function(x, vi, sei, ni, subset, data, model="rma", predictor="sei", ret.fit=FALSE, digits, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   model <- match.arg(model, c("lm", "rma"))
   predictor <- match.arg(predictor, c("sei", "vi", "ni", "ninv", "sqrtni", "sqrtninv"))

   ddd <- list(...)

   .chkdots(ddd, c("level", "method", "test"))

   #########################################################################

   ### check if data argument has been specified

   if (missing(data))
      data <- NULL

   if (is.null(data)) {
      data <- sys.frame(sys.parent())
   } else {
      if (!is.data.frame(data))
         data <- data.frame(data)
   }

   mf <- match.call()

   x <- .getx("x", mf=mf, data=data)

   #########################################################################

   if (inherits(x, "rma")) {

      .chkclass(class(x), must="rma", notav=c("robust.rma", "rma.glmm", "rma.mv", "rma.ls", "rma.nn", "rma.uni.selmodel"))

      if (!missing(vi) || !missing(sei) || !missing(subset))
         warning(mstyle$warning("Arguments 'vi', 'sei', and 'subset' ignored when 'x' is a model object."), call.=FALSE)

      yi <- x$yi
      vi <- x$vi

      if (missing(ni)) {

         ni <- x$ni # may be NULL

      } else {

         ni <- .getx("ni", mf=mf, data=data, checknumeric=TRUE)

         if (!is.null(ni)) {

            if (length(ni) != x$k.all)
               stop(mstyle$stop(paste0("Length of variable specified via 'ni' (", length(ni), ") does not correspond to the size of the original dataset (", x$k.all, ").")))

            if (!is.null(x$subset))
               ni <- ni[x$subset]

            if (inherits(x, "rma.mh") || inherits(x, "rma.peto")) {
               ni <- ni[x$not.na.yivi]
            } else {
               ni <- ni[x$not.na]
            }

         }

      }

      k <- length(yi)

      ### set defaults for digits

      if (missing(digits)) {
         digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
      } else {
         digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
      }

      p <- x$p

      if (inherits(x, "rma.mh") || inherits(x, "rma.peto")) {
         X <- cbind(rep(1,k))
      } else {
         X <- x$X
      }

      if (!is.null(ddd$level)) {
         level <- .level(ddd$level)
      } else {
         level <- x$level
      }
      method   <- ifelse(is.null(ddd$method), x$method, ddd$method)
      test     <- ifelse(is.null(ddd$test), x$test, ddd$test)
      weights  <- x$weights
      weighted <- x$weighted
      tau2     <- ifelse(x$tau2.fix, x$tau2, NA)
      control  <- x$control

   } else {

      if (!.is.vector(x))
         stop(mstyle$stop("Argument 'x' must be a vector or an 'rma' model object."))

      yi <- x

      ### check if yi is numeric

      if (!is.numeric(yi))
         stop(mstyle$stop("The object/variable specified for the 'x' argument is not numeric."))

      ### set defaults for digits

      if (missing(digits)) {
         digits <- .set.digits(dmiss=TRUE)
      } else {
         digits <- .set.digits(digits, dmiss=FALSE)
      }

      if (!is.null(ddd$level)) {
         level <- .level(ddd$level)
      } else {
         level <- .05
      }

      k <- length(yi)

      vi     <- .getx("vi",     mf=mf, data=data, checknumeric=TRUE)
      sei    <- .getx("sei",    mf=mf, data=data, checknumeric=TRUE)
      ni     <- .getx("ni",     mf=mf, data=data, checknumeric=TRUE)
      subset <- .getx("subset", mf=mf, data=data)

      if (is.null(vi)) {
         if (!is.null(sei))
            vi <- sei^2
      }

      if (is.null(vi))
         stop(mstyle$stop("Must specify 'vi' or 'sei' argument."))

      ### check length of yi and vi

      if (length(vi) != k)
         stop(mstyle$stop("Length of 'yi' and 'vi' (or 'sei') is not the same."))

      ### check length of yi and ni

      if (!is.null(ni) && length(ni) != k)
         stop(mstyle$stop("Length of 'yi' and 'ni' is not the same."))

      ### check 'vi' argument for potential misuse

      .chkviarg(mf$vi)

      ### if ni has not been specified but is an attribute of yi, get it

      if (is.null(ni) && !is.null(attr(yi, "ni")))
         ni <- attr(yi, "ni")

      ### check length of yi and ni (only if ni is not NULL)
      ### if there is a mismatch, then ni cannot be trusted, so set it to NULL

      if (!is.null(ni) && length(ni) != k)
         ni <- NULL

      ### if ni is now available, add it (back) as an attribute to yi

      if (!is.null(ni))
         attr(yi, "ni") <- ni

      ### if a subset of studies is specified

      if (!is.null(subset)) {
         subset <- .setnafalse(subset, k=k)
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

      p <- 1L
      k <- length(yi)
      X <- cbind(rep(1,k))

      method   <- ifelse(is.null(ddd$method), "REML", ddd$method)
      test     <- ifelse(is.null(ddd$test), "z", ddd$test)
      weights  <- NULL
      weighted <- TRUE
      tau2     <- NA
      control  <- list()

   }

   #########################################################################

   if (predictor == "sei")
      X <- cbind(X, sei=sqrt(vi))

   if (predictor == "vi")
      X <- cbind(X, vi=vi)

   if (is.element(predictor, c("ni", "ninv", "sqrtni", "sqrtninv"))) {

      if (is.null(ni)) {

         stop(mstyle$stop("No sample size information available to use this predictor."))

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

      ddd$level  <- NULL
      ddd$method <- NULL
      ddd$test   <- NULL

      args <- list(yi=yi, vi=vi, weights=weights, mods=X, intercept=FALSE, method=method, weighted=weighted, test=test, level=level, tau2=tau2, control=control, ddd)
      fit  <- .do.call(rma.uni, args)
      zval <- fit$zval[p+1]
      pval <- fit$pval[p+1]
      ddf  <- fit$ddf

   } else {

      yi   <- c(yi) # remove attributes
      fit  <- lm(yi ~ X - 1, weights=1/vi)
      tmp  <- summary(fit)
      zval <- coef(tmp)[p+1,3]
      pval <- coef(tmp)[p+1,4]
      ddf  <- fit$df.residual

   }

   ### get the 'limit estimate'

   if (predictor %in% c("sei", "vi", "ninv", "sqrtninv") && p == 1L && .is.intercept(X[,1])) {

      if (model=="lm") {
         est <- coef(tmp)[1,1]
         ci.lb <- est - qt(level/2, df=ddf, lower.tail=FALSE) * coef(tmp)[1,2]
         ci.ub <- est + qt(level/2, df=ddf, lower.tail=FALSE) * coef(tmp)[1,2]
      } else {
         est <- coef(fit)[1]
         ci.lb <- fit$ci.lb[1]
         ci.ub <- fit$ci.ub[1]
      }

   } else {

      est <- ci.lb <- ci.ub <- NULL

   }

   res <- list(model=model, predictor=predictor, zval=zval, pval=pval, dfs=ddf, ddf=ddf, method=fit$method, digits=digits, ret.fit=ret.fit, fit=fit, est=est, ci.lb=ci.lb, ci.ub=ci.ub)

   class(res) <- "regtest"
   return(res)

}
