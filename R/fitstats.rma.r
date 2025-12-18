fitstats.rma <- function(object, ..., REML) {

   mstyle <- .get.mstyle()

   .chkclass(class(object), must="rma")

   # unless the 'REML' argument was specified, the 'method' of the first object
   # determines whether to show fit statistics based on the ML or REML likelihood

   if (missing(REML)) {
      if (object$method == "REML") {
         REML <- TRUE
      } else {
         REML <- FALSE
      }
   }

   if (missing(...)) {

      # if there is just 'object'

      if (REML) {
         out <- cbind(object$fit.stats$REML)
         colnames(out) <- "REML"
      } else {
         out <- cbind(object$fit.stats$ML)
         colnames(out) <- "ML"
      }

   } else {

      # if there is 'object' and additional objects via ...

      if (REML) {
         out <- sapply(list(object, ...), function(x) x$fit.stats$REML)
      } else {
         out <- sapply(list(object, ...), function(x) x$fit.stats$ML)
      }

      out <- data.frame(out)

      # get the names of the objects; same idea as in stats:::AIC.default

      cl <- match.call()
      cl$REML <- NULL
      names(out) <- as.character(cl[-1L])

      # check that all models were fitted to the same data

      chksums <- sapply(list(object, ...), function(x) x$chksumyi)

      if (any(chksums[1] != chksums))
         warning(mstyle$warning("Models not all fitted to the same data."), call.=FALSE)

   }

   rownames(out) <- c("logLik", "deviance", "AIC", "BIC", "AICc")
   return(out)

   #print(fmtx(out, object$digits[["fit"]]), quote=FALSE)
   #invisible(out)

}
