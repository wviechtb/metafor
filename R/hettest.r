hettest <- function(x, vi, sei, subset, data, method="REML", test="score", boot=TRUE, progbar=TRUE, digits, control, ...) {

   #########################################################################

   mstyle <- .get.mstyle()

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   test <- tolower(test)
   test <- match.arg(test, c("lrt", "wald", "score"))

   # check 'method' argument

   method <- method[1] # just in case

   if (!is.element(method, c("ML","REML")))
      stop(mstyle$warning("Test must be based on ML/REML estimation."), call.=FALSE)

   if (missing(control))
      control <- list()

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

      if (!missing(vi) || !missing(sei) || !missing(subset))
         warning(mstyle$warning("Arguments 'vi', 'sei', and 'subset' ignored when 'x' is a model object."), call.=FALSE)

      if (!is.null(mf$method))
         warning(mstyle$warning("Argument 'method' ignored when 'x' is a model object."), call.=FALSE)

      # set defaults for 'digits'

      if (missing(digits)) {
         digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
      } else {
         digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
      }

      if (is.element(x$method, c("ML","REML"))) {
         hettest(c(x$yi), vi=x$vi, method=x$method, test=test, boot=boot, progbar=progbar, digits=digits, model=x, ...)
      } else {
         warning(mstyle$warning(paste0("Model was not fitted with ML or REML estimation, but the test will be based on ", method, " estimation.")), call.=FALSE)
         hettest(c(x$yi), vi=x$vi, method=method, test=test, boot=boot, progbar=progbar, digits=digits, ...)
      }

   } else {

      #########################################################################

      if (!.is.vector(x))
         stop(mstyle$stop("Argument 'x' must be a vector with estimates or an 'rma' model object."))

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

      # set defaults for 'boot' (can be numeric to specify the number of iterations)

      if (is.numeric(boot)) {
         boot <- round(boot[1])
         boot[boot <= 1] <- 1
         iter <- boot
         boot <- TRUE
      } else {
         iter <- 1000
      }

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
         sim[,1] <- c(yi) # make sure that first bootstrap dataset corresponds to the actual data

      }

      #########################################################################

      # initial value(s) for tau2i in .hettest.esttau2i()

      tau2i.init <- rep(res0$tau2, k)
      tau2i.init[tau2i.init <= 0.01] <- 0.01

      # control settings

      con <- list(tau2i.init = tau2i.init,
                  threshold  = 10^-8,
                  maxiter    = 1000)

      con.pos <- pmatch(names(control), names(con))
      con[c(na.omit(con.pos))] <- control[!is.na(con.pos)]

      con$tau2i.init <- .expand1(con$tau2i.init, k)

      if (anyNA(con$tau2i.init))
         stop(mstyle$stop(paste0("No missing values allowed in 'tau2i.init'.")))

      if (any(con$tau2i.init < 0))
         stop(mstyle$stop("Value(s) of 'tau2i.init' must be >= 0."))

      if (length(con$tau2i.init) != k)
         stop(mstyle$stop(paste0("Length of the 'tau2i.init' argument (", length(con$tau2i.init), ") does not match the actual number of variance components (", k, ").")))

      #########################################################################

      if (test == "lrt") {

         out <- .hettest.lrt(yi, vi, method=method, res0=res0, mom=mom, tau2i.init=con$tau2i.init, threshold=con$threshold, maxiter=con$maxiter)

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

               tmp <- try(.hettest.lrt(sim[,b], vi, method=method, res0=tmp, mom=mom, tau2i.init=con$tau2i.init, threshold=con$threshold, maxiter=con$maxiter), silent=TRUE)
               if (inherits(tmp, "try-error"))
                  next
               x2.boot[b] <- tmp$statistic

            }

            if (progbar)
               pbapply::closepb(pbar)

            pval <- mean(x2.boot >= x2, na.rm=TRUE)

         }

      }

      if (test == "wald") {

         out <- .hettest.wald(yi, vi, method=method, res0=res0, mom=mom, tau2i.init=con$tau2i.init, threshold=con$threshold, maxiter=con$maxiter)

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

               tmp <- try(.hettest.wald(sim[,b], vi, method=method, res0=tmp, mom=mom, tau2i.init=con$tau2i.init, threshold=con$threshold, maxiter=con$maxiter), silent=TRUE)
               if (inherits(tmp, "try-error"))
                  next
               x2.boot[b] <- tmp$statistic

            }

            if (progbar)
               pbapply::closepb(pbar)

            pval <- mean(x2.boot >= x2, na.rm=TRUE)

         }

      }

      if (test == "score") {

         out <- .hettest.scoretest(yi, vi, mu_hat=coef(res0), tau2_hat=res0$tau2, method=method)

         x2   <- out$statistic
         df   <- out$df
         pval <- out$pval

         # the score test does not require that we estimate the tau^2_i values, but do this anyway

         if (mom) {
            hi <- hatvalues(res0)
            vari <- resid(res0)^2 / (1-hi)
            tau2i <- pmax(0, vari - vi)
         } else {
            tau2i <- try(.hettest.esttau2i(yi, vi, method=method, tau2i.init=con$tau2i.init, threshold=con$threshold, maxiter=con$maxiter), silent=TRUE)
            if (inherits(tau2i, "try-error"))
               tau2i <- NA
         }

         if (boot) {

            if (progbar)
               pbar <- pbapply::startpb(min=0, max=iter)

            for (b in seq_len(iter)) {

               if (progbar)
                  pbapply::setpb(pbar, b)

               tmp <- try(rma(sim[,b], vi, method=method), silent=TRUE)

               if (inherits(tmp, "try-error"))
                  next

               tmp <- try(.hettest.scoretest(sim[,b], vi, mu_hat=coef(tmp), tau2_hat=tmp$tau2, method=method), silent=TRUE)
               if (inherits(tmp, "try-error"))
                  next
               x2.boot[b] <- tmp$statistic

            }

            if (progbar)
               pbapply::closepb(pbar)

            pval <- mean(x2.boot >= x2, na.rm=TRUE)

         }

      }

      #########################################################################

      res <- list(x2=x2, df=df, pval=pval, method=method, test=test, boot=boot, iter=iter, digits=digits)

      if (boot)
         res$x2.boot <- x2.boot

      if (test %in% c("lrt","wald")) {
         res$tau2i <- out$tau2i
      } else {
         res$tau2i <- tau2i
      }

      if (test  == "wald")
         res$se.tau2i <- out$se.tau2i

      class(res) <- "hettest"
      return(res)

   }

}
