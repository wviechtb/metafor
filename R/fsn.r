fsn <- function(yi, vi, sei, data, type="Rosenthal", alpha=.05, target, weighted=FALSE, subset, digits) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   type <- match.arg(type, c("Rosenthal", "Orwin", "Rosenberg"))

   if (missing(target))
      target <- NULL

   ### set defaults for digits

   if (missing(digits)) {
      digits <- .set.digits(dmiss=TRUE)
   } else {
      digits <- .set.digits(digits, dmiss=FALSE)
   }

   #########################################################################

   ###### data setup

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
   mf.yi      <- mf[[match("yi",      names(mf))]]
   mf.vi      <- mf[[match("vi",      names(mf))]]
   mf.sei     <- mf[[match("sei",     names(mf))]]
   #mf.weights <- mf[[match("weights", names(mf))]]
   mf.subset  <- mf[[match("subset",  names(mf))]]
   yi      <- eval(mf.yi,      data, enclos=sys.frame(sys.parent()))
   vi      <- eval(mf.vi,      data, enclos=sys.frame(sys.parent()))
   sei     <- eval(mf.sei,     data, enclos=sys.frame(sys.parent()))
   #weights <- eval(mf.weights, data, enclos=sys.frame(sys.parent()))
   subset  <- eval(mf.subset,  data, enclos=sys.frame(sys.parent()))

   if (type %in% c("Rosenthal", "Rosenberg") || (type == "Orwin" && weighted)) {
      if (is.null(vi)) {
         if (is.null(sei)) {
            stop(mstyle$stop("Need to specify 'vi' or 'sei' argument."))
         } else {
            vi <- sei^2
         }
      }
   } else {
      vi <- rep(0, length(yi))
   }

   ### check length of yi and vi

   if (length(yi) != length(vi))
      stop(mstyle$stop("Length of 'yi' and 'vi' (or 'sei') is not the same."))

   ### if a subset of studies is specified

   if (!is.null(subset)) {
      yi <- yi[subset]
      vi <- vi[subset]
   }

   ### check for NAs in yi/vi and act accordingly

   yivi.na <- is.na(yi) | is.na(vi)

   if (any(yivi.na)) {

      not.na <- !yivi.na

      if (na.act == "na.omit" || na.act == "na.exclude" || na.act == "na.pass") {
         yi <- yi[not.na]
         vi <- vi[not.na]
      }

      if (na.act == "na.fail")
         stop(mstyle$stop("Missing values in results."))

   }

   #########################################################################

   if (type == "Rosenthal") {

      k      <- length(yi)
      zi     <- yi / sqrt(vi)
      z.avg  <- abs(sum(zi) / sqrt(k))
      pval   <- pnorm(z.avg, lower.tail=FALSE)
      fsnum  <- max(0, k * (z.avg / qnorm(alpha, lower.tail=FALSE))^2 - k)
      meanes <- NA
      target <- NA

   }

   if (type == "Orwin") {

      k      <- length(yi)

      if (weighted) {
         wi <- 1/vi
         meanes <- sum(wi*yi)/sum(wi)
      } else {
         meanes <- mean(yi)
      }

      if (is.null(target))
         target <- meanes / 2

      if (identical(target, 0) || sign(target) != sign(meanes)) {
         fsnum <- Inf
      } else {
         fsnum <- max(0, k * (meanes - target) / target)
      }
      pval <- NA

   }

   if (type == "Rosenberg") {

      k      <- length(yi)
      wi     <- 1/vi
      meanes <- sum(wi*yi)/sum(wi)
      zval   <- meanes / sqrt(1/sum(wi))
      w.p    <- (sum(wi*yi) / qnorm(alpha/2, lower.tail=FALSE))^2 - sum(wi)
      pval   <- 2*pnorm(abs(zval), lower.tail=FALSE)
      fsnum  <- max(0, k*w.p/sum(wi))
      target <- NA

   }

   if (abs(fsnum - round(fsnum)) >= .Machine$double.eps^0.5) {
      fsnum <- ceiling(fsnum)
   } else {
      fsnum <- round(fsnum)
   }

   #########################################################################

   res <- list(type=type, fsnum=fsnum, alpha=alpha, pval=pval, meanes=meanes, target=target, digits=digits)

   class(res) <- "fsn"
   return(res)

}
