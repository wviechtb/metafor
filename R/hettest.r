hettest <- function(x, vi, sei, subset, data,
   method="REML", test="score", boot=TRUE, progbar=TRUE, digits, ...) {

   #########################################################################

   # control arguments for .hettest.esttau2i maxiter and threshold
   # what about the lrt for the random heteroscedastic heterogeneity model?
   # it's analogous to the lrt for tau^2 - probably doesn't work as well, so
   # one has to get this via anova(), not this function

   mstyle <- .get.mstyle()

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   test <- tolower(test)
   test <- match.arg(test, c("score", "lrt"))

   #########################################################################

   # check if the 'data' argument was specified

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

      on.exit(options(na.action=na.act), add=TRUE)

      .chkclass(class(x), must="rma", notav=c("rma.glmm", "rma.mv", "robust.rma", "rma.ls", "rma.gen", "rma.uni.selmodel"))

      if (is.null(x$yi) || is.null(x$vi))
         stop(mstyle$stop("Information needed to carry out the test is not available in the model object."))

      if (!x$int.only)
         stop(mstyle$stop("Test only applicable to models without moderators."))

      # set defaults for 'digits'

      if (missing(digits)) {
         digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
      } else {
         digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
      }

      hettest(c(x$yi), vi=x$vi, method=x$method, test=test, boot=boot, progbar=progbar, digits=digits, model=x, ...)

   } else {

      #########################################################################

      if (!.is.vector(x))
         stop(mstyle$stop("Argument 'x' must be a vector or an 'rma' model object."))

      yi <- x

      # check if 'yi' is numeric

      if (!is.numeric(yi))
         stop(mstyle$stop("The object/variable specified for the 'x' argument is not numeric."))

      # set defaults for 'digits'

      if (missing(digits)) {
         digits <- .set.digits(dmiss=TRUE)
      } else {
         digits <- .set.digits(digits, dmiss=FALSE)
      }

      vi     <- .getx("vi",     mf=mf, data=data, checknumeric=TRUE)
      sei    <- .getx("sei",    mf=mf, data=data, checknumeric=TRUE)
      subset <- .getx("subset", mf=mf, data=data)

      if (is.null(vi)) {
         if (!is.null(sei))
            vi <- sei^2
      }

      if (is.null(vi))
         stop(mstyle$stop("Must specify the 'vi' or 'sei' argument."))

      # check length of 'yi' and 'vi'

      if (length(yi) != length(vi))
         stop(mstyle$stop("Length of 'yi' and 'vi' (or 'sei') are not the same."))

      # check 'vi' argument for potential misuse

      .chkviarg(mf$vi)

      #########################################################################

      ddd <- list(...)

      .chkdots(ddd, c("model", "verbose", "seed", "mom"))

      verbose <- .chkddd(ddd$verbose, FALSE)
      mom     <- .chkddd(ddd$mom,     FALSE)

      #########################################################################

      k.f <- length(yi)

      # if a subset of studies is specified

      if (!is.null(subset)) {
         subset <- .chksubset(subset, k.f)
         yi <- .getsubset(yi, subset)
         vi <- .getsubset(vi, subset)
      }

      # check for NAs and act accordingly

      has.na <- is.na(yi) | is.na(vi)

      if (any(has.na)) {

         not.na <- !has.na

         if (na.act == "na.omit" || na.act == "na.exclude" || na.act == "na.pass") {

            yi <- yi[not.na]
            vi <- vi[not.na]
            warning(mstyle$warning(paste(sum(has.na), ifelse(sum(has.na) > 1, "studies", "study"), "with NAs omitted from test.")), call.=FALSE)

         }

         if (na.act == "na.fail")
            stop(mstyle$stop("Missing values in data."))

      }

      k <- length(yi)

      if (k <= 1L)
         stop(mstyle$stop("Stopped because k <= 1."))

      #########################################################################

      # set defaults for 'boot'

      if (is.numeric(boot) && length(boot) == 1L) {
         boot <- round(boot)
         boot[boot <= 1] <- 1
         iter <- boot
         boot <- TRUE
      } else {
         iter <- 1000
      }

      # set default for 'progbar'

      if (missing(progbar))
         progbar <- ifelse(boot, TRUE, FALSE)

      # check 'method' argument

      method <- method[1] # just in case

      if (!is.element(method, c("ML","REML")))
         warning(mstyle$warning("Test should be based on ML/REML estimation."), call.=FALSE)

      #########################################################################

      if (is.null(ddd$model)) {

         # fit the null model

         res0 <- try(rma(yi, vi, method=method, verbose=verbose), silent=!verbose)

         if (inherits(res0, "try-error"))
            stop(mstyle$stop("Could not fit the null model."))

      } else {

         res0 <- ddd$model

      }

      if (boot) {

         if (!is.null(ddd$seed))
            set.seed(ddd$seed)

         x2.boot <- rep(NA_real_, iter)
         sim <- simulate(res0, iter)
         sim[,1] <- yi

      }

      #########################################################################

      if (test == "score") {

         out <- .hettest.scoretest(yi, vi, mu_hat=coef(res0), tau2_hat=res0$tau2, method=method)

         x2   <- out$statistic
         df   <- out$df
         pval <- out$pval

         if (boot) {

            if (progbar)
               pbar <- pbapply::startpb(min=0, max=iter)

            for (b in seq_len(iter)) {

               if (progbar)
                  pbapply::setpb(pbar, b)

               tmp <- try(rma(sim[,b], vi, method=method), silent=TRUE)

               if (inherits(tmp, "try-error"))
                  next

               tmp <- .hettest.scoretest(sim[,b], vi, mu_hat=coef(tmp), tau2_hat=tmp$tau2, method=method)
               x2.boot[b] <- tmp$statistic

            }

            if (progbar)
               pbapply::closepb(pbar)

            pval <- mean(x2.boot >= x2, na.rm=TRUE)

         }

      }

      if (test == "lrt") {

         out <- .hettest.lrt(yi, vi, method=method, res0=res0, mom=mom)

         x2   <- out$statistic
         df   <- out$df
         pval <- out$pval

         if (boot) {

            if (progbar)
               pbar <- pbapply::startpb(min=0, max=iter)

            for (b in seq_len(iter)) {

               if (progbar)
                  pbapply::setpb(pbar, b)

               tmp <- try(rma(sim[,b], vi, method=method), silent=TRUE)

               if (inherits(tmp, "try-error"))
                  next

               tmp <- .hettest.lrt(sim[,b], vi, method=method, res0=tmp, mom=mom)
               x2.boot[b] <- tmp$statistic

            }

            if (progbar)
               pbapply::closepb(pbar)

            pval <- mean(x2.boot >= x2, na.rm=TRUE)

         }

      }

      res <- list(x2=x2, df=df, pval=pval, digits=digits)

      if (boot)
         res$x2.boot <- x2.boot

      if (test == "lrt")
         res$tau2i <- out$tau2i

      class(res) <- "hettest"
      return(res)

   }

}
