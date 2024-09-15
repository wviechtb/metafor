ranktest <- function(x, vi, sei, subset, data, digits, ...) {

   #########################################################################

   mstyle <- .get.mstyle()

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   ddd <- list(...)

   .chkdots(ddd, c("exact"))

   exact <- .chkddd(ddd$exact, TRUE)

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

   ############################################################################

   if (inherits(x, "rma")) {

      if (is.null(x$yi) || is.null(x$vi))
         stop(mstyle$stop("Information needed to carry out the test is not available in the model object."))

      if (!missing(vi) || !missing(sei) || !missing(subset))
         warning(mstyle$warning("Arguments 'vi', 'sei', and 'subset' ignored when 'x' is a model object."), call.=FALSE)

      if (!x$int.only)
         stop(mstyle$stop("Test only applicable to models without moderators."))

      yi <- x$yi
      vi <- x$vi

      ### set defaults for digits

      if (missing(digits)) {
         digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
      } else {
         digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
      }

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

      vi     <- .getx("vi",     mf=mf, data=data, checknumeric=TRUE)
      sei    <- .getx("sei",    mf=mf, data=data, checknumeric=TRUE)
      subset <- .getx("subset", mf=mf, data=data)

      if (is.null(vi)) {
         if (!is.null(sei))
            vi <- sei^2
      }

      if (is.null(vi))
         stop(mstyle$stop("Must specify 'vi' or 'sei' argument."))

      ### check length of yi and vi

      if (length(yi) != length(vi))
         stop(mstyle$stop("Length of 'yi' and 'vi' (or 'sei') is not the same."))

      ### check 'vi' argument for potential misuse

      .chkviarg(mf$vi)

      #########################################################################

      ### if a subset of studies is specified

      if (!is.null(subset)) {
         subset <- .chksubset(subset, length(yi))
         yi <- .getsubset(yi, subset)
         vi <- .getsubset(vi, subset)
      }

      ### check for NAs and act accordingly

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

   }

   #########################################################################

   wi <- 1/vi
   theta <- weighted.mean(yi, wi)
   vb <- 1 / sum(wi)

   vi.star <- vi - vb
   yi.star <- (yi - theta) / sqrt(vi.star)
   res <- cor.test(yi.star, vi, method="kendall", exact=exact)

   pval <- res$p.value
   tau  <- res$estimate

   res <- list(tau=tau, pval=pval, digits=digits)

   class(res) <- "ranktest"
   return(res)

}
