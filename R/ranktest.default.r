ranktest.default <- function(x, vi, sei, subset, digits, ...) {

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

   ddd <- list(...)

   .chkdots(ddd, c("exact"))

   if (is.null(ddd$exact)) {
      exact <- TRUE
   } else {
      exact <- ddd$exact
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

   #########################################################################

   ### if a subset of studies is specified

   if (!is.null(subset)) {
      yi <- yi[subset]
      vi <- vi[subset]
   }

   ### check for NAs and act accordingly

   has.na <- is.na(yi) | is.na(vi)

   if (any(has.na)) {

      not.na <- !has.na

      if (na.act == "na.omit" || na.act == "na.exclude" || na.act == "na.pass") {

         yi <- yi[not.na]
         vi <- vi[not.na]
         warning(mstyle$warning("Studies with NAs omitted from test."), call.=FALSE)

      }

      if (na.act == "na.fail")
         stop(mstyle$stop("Missing values in data."))

   }

   #########################################################################

   res  <- rma.uni(yi, vi, method="FE")
   beta <- c(res$beta)
   vb   <- c(res$vb)

   vi.star <- vi - vb
   yi.star <- (yi - beta) / sqrt(vi.star)
   res <- cor.test(yi.star, vi, method="kendall", exact=exact)

   pval <- res$p.value
   tau  <- res$estimate

   res <- list(tau=tau, pval=pval, digits=digits)

   class(res) <- "ranktest"
   return(res)

}
