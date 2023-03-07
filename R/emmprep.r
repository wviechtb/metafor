emmprep <- function(x, verbose=FALSE, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

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

      if (is.element(x$measure, c("RR","OR","PETO","IRR","ROM","D2OR","D2ORL","D2ORN","CVR","VR","PLN","IRLN","SDLN","MNLN","CVLN","ROMC","CVRC","VRC","REH"))) {
         out@misc$tran <- "log"
         #out@misc$tran <- emmeans::make.tran("genlog", 0)
         #out <- update(out, emmeans::make.tran("genlog", 0))
         if (verbose) cat("Transformation:     log\n")
      }

      if (is.element(x$measure, c("PLO"))) {
         out@misc$tran <- "logit"
         if (verbose) cat("Transformation:     logit\n")
      }

      if (is.element(x$measure, c("PAS"))) {
         out <- update(out, emmeans::make.tran("asin.sqrt", 1))
         if (verbose) cat("Transformation:     asin.sqrt\n")
      }

      if (is.element(x$measure, c("IRS"))) {
         out@misc$tran <- "sqrt"
         if (verbose) cat("Transformation:     sqrt\n")
      }

      if (is.element(x$measure, c("ZCOR","ZPCOR"))) {
         out@misc$tran$linkfun  <- transf.rtoz
         out@misc$tran$linkinv  <- transf.ztor
         out@misc$tran$mu.eta   <- function(eta) 1/cosh(eta)^2
         out@misc$tran$valideta <- function(eta) all(is.finite(eta)) && all(abs(eta) <= 1)
         out@misc$tran$name     <- "r-to-z"
         if (verbose) cat("Transformation:     r-to-z\n")
      }

      #if (is.element(x$measure, c("AHW")))
      #   out@misc$tran <- "???"

      #if (is.element(x$measure, c("ABT")))
      #   out@misc$tran <- "???"

   } else {

      if (verbose) cat("Transformation:    ", ddd$tran, "\n")

   }

   if (verbose)
      .space()

   return(out)

}

############################################################################
