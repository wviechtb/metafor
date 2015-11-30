# Note: Works with "robust.rma" objects.

fitstats.rma <- function(object, ..., REML) {

   if (!is.element("rma", class(object)))
      stop("Argument 'object' must be an object of class \"rma\".")

   if (missing(REML)) {
      if (object$method == "REML") {
         REML <- TRUE
      } else {
         REML <- FALSE
      }
   }

   if (missing(...)) {

      ### if there is just 'object'

      if (REML) {
         out <- cbind(object$fit.stats$REML)
         colnames(out) <- "REML"
      } else {
         out <- cbind(object$fit.stats$ML)
         colnames(out) <- "ML"
      }

   } else {

      ### if there is 'object' and additional objects via ...

      if (REML) {
         out <- sapply(list(object, ...), function(x) x$fit.stats$REML)
      } else {
         out <- sapply(list(object, ...), function(x) x$fit.stats$ML)
      }

      out <- data.frame(out)

      ### get names of objects; same idea as in stats:::AIC.default

      cl <- match.call()
      cl$REML <- NULL
      names(out) <- as.character(cl[-1L])

      ### check that all models were fitted to the same data

      yis <- lapply(list(object, ...), function(x) as.vector(x$yi))
      if (!all(sapply(yis[-1], function(x) identical(x, yis[[1]]))))
         warning("Models not all fitted to the same data.")

   }

   rownames(out) <- c("logLik:", "deviance:", "AIC:", "BIC:", "AICc:")
   return(out)

}
