emmprep <- function(x, verbose=FALSE, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(x), must="rma")

   if (!requireNamespace("emmeans", quietly=TRUE))
      stop(mstyle$stop("Please install the 'emmeans' package to use this function."))

   if (any(x$coef.na))
      stop(mstyle$stop("Cannot use function when some redundant predictors were dropped from the model."))

   ### check if a formula is available

   formula <- formula(x)

   if (is.null(formula) && x$int.only)
      formula <- ~ 1

   if (is.null(formula))
      stop(mstyle$stop("Cannot use function when model was fitted without a formula specification."))

   if (verbose) {
      .space()
      cat("Extracted formula:  ~", paste(paste(formula)[-1], collapse=""), "\n")
   }

   ### get coefficients and corresponding var-cov matrix

   b  <- coef(x, type="beta")
   vb <- vcov(x, type="beta")

   ### change intrcpt to (Intercept)

   names(b)     <- sub("intrcpt", "(Intercept)", names(b))
   rownames(vb) <- sub("intrcpt", "(Intercept)", rownames(vb))
   colnames(vb) <- sub("intrcpt", "(Intercept)", colnames(vb))

   #########################################################################

   ddd <- list(...)

   ### get data and apply subsetting / removal of missings as needed

   if (is.null(ddd$data)) {

      dat <- x$data

      if (is.null(dat))
         stop(mstyle$stop("Cannot use function when the model object does not contain the original data."))

      if (!is.null(x$subset))
         dat <- dat[x$subset,,drop=FALSE]

      dat <- dat[x$not.na,,drop=FALSE]

   } else {

      dat <- ddd$data
      ddd$data <- NULL

   }

   ### set the degrees of freedom (use minimum value if there are multiple)

   if (is.null(ddd$df)) {

      if (is.na(x$ddf[1])) {
         ddf <- Inf
      } else {
         ddf <- min(x$ddf)
      }

   } else {

      ddf <- ddd$df
      ddd$df <- NULL

   }

   if (verbose && is.finite(ddf))
      cat("Degrees of freedom:", round(ddf, 2), "\n")

   ### set sigma for bias adjustment

   if (is.null(ddd$sigma)) {

      if (!inherits(x, c("rma.ls","rma.mv"))) {
         sigma <- sqrt(x$tau2)
      } else {
         sigma <- NA_real_
      }

   } else {

      sigma <- ddd$sigma
      ddd$sigma <- NULL

   }

   if (verbose && !is.na(sigma) && !is.element(x$method, c("FE","EE","CE")))
      cat("Value of tau^2:    ", round(sigma^2, 4), "\n")

   if (is.na(sigma))
      sigma <- 0

   ### create grid

   #out <- emmeans::qdrg(formula=formula, data=dat, coef=b, vcov=vb, df=ddf, sigma=sigma, ...)
   out <- do.call(emmeans::qdrg, c(list(formula=formula, data=dat, coef=b, vcov=vb, df=ddf, sigma=sigma), ddd))

   ### set (back)transformation

   if (is.null(ddd$tran)) {

      if (is.element(x$measure, c("RR","OR","MPORM","PETO","MPRR","MPOR","MPORC","MPPETO","IRR","ROM","D2OR","D2ORL","D2ORN","VR","CVR","PLN","IRLN","MNLN","SDLN","CVLN","ROMC","VRC","CVRC","REH","HR"))) {
         out@misc$tran <- "log"
         #out@misc$tran <- emmeans::make.tran("genlog", 0)
         #out <- update(out, emmeans::make.tran("genlog", 0))
         if (verbose) cat("Transformation:     log\n")
      }

      if (is.element(x$measure, c("PLO"))) {
         out@misc$tran <- "logit"
         if (verbose) cat("Transformation:     logit\n")
      }

      if (is.element(x$measure, c("PRZ"))) {
         out@misc$tran <- "probit"
         if (verbose) cat("Transformation:     probit\n")
      }

      if (is.element(x$measure, c("PAS"))) {
         out <- update(out, emmeans::make.tran("asin.sqrt", 1))
         if (verbose) cat("Transformation:     asin.sqrt\n")
      }

      if (is.element(x$measure, c("IRS"))) {
         out@misc$tran <- "sqrt"
         if (verbose) cat("Transformation:     sqrt\n")
      }

      if (is.element(x$measure, c("ZPHI","ZTET","ZPB","ZBIS","ZCOR","ZPCOR","ZSPCOR"))) {
         out@misc$tran$linkfun  <- transf.rtoz
         out@misc$tran$linkinv  <- transf.ztor
         out@misc$tran$mu.eta   <- function(eta) 1/cosh(eta)^2 # derivative of transf.ztor(eta) (= tanh(eta))
         out@misc$tran$valideta <- function(eta) all(is.finite(eta)) && all(abs(eta) <= 1)
         out@misc$tran$name     <- "r-to-z"
         if (verbose) cat("Transformation:     r-to-z\n")
      }

      if (is.element(x$measure, c("ZR2","ZR2F"))) {
         out@misc$tran$linkfun  <- transf.r2toz
         out@misc$tran$linkinv  <- transf.ztor2
         out@misc$tran$mu.eta   <- function(eta) 2*sinh(eta)/cosh(eta)^3 # derivative of transf.ztor2(eta) (= tanh(eta)^2)
         out@misc$tran$valideta <- function(eta) all(is.finite(eta)) && all(eta <= 1) && all(eta >= 0)
         out@misc$tran$name     <- "r-to-z"
         if (verbose) cat("Transformation:     r-to-z\n")
      }

      if (is.element(x$measure, c("AHW"))) {
         out@misc$tran$linkfun  <- transf.ahw
         out@misc$tran$linkinv  <- transf.iahw
         out@misc$tran$mu.eta   <- function(eta) 3*(1-eta)^2
         out@misc$tran$valideta <- function(eta) all(is.finite(eta)) && all(eta <= 1) && all(eta >= 0)
         out@misc$tran$name     <- "ahw"
         if (verbose) cat("Transformation:     ahw\n")
      }

      if (is.element(x$measure, c("ABT"))) {
         out@misc$tran$linkfun  <- transf.abt
         out@misc$tran$linkinv  <- transf.iabt
         out@misc$tran$mu.eta   <- function(eta) 1/(1-eta)
         out@misc$tran$valideta <- function(eta) all(is.finite(eta)) && all(eta <= 1) && all(eta >= 0)
         out@misc$tran$name     <- "abt"
         if (verbose) cat("Transformation:     abt\n")
      }

   } else {

      if (verbose) cat("Transformation:    ", ddd$tran, "\n")

   }

   if (verbose)
      .space()

   return(out)

}

############################################################################
