AIC.rma <- function(object, ..., k=2, correct=FALSE) {

   mstyle <- .get.mstyle()

   .chkclass(class(object), must="rma")

   if (missing(...)) {

      # if there is just 'object'

      if (object$method == "REML") {
         out <- ifelse(correct, object$fit.stats["AICc","REML"], object$fit.stats["AIC","REML"])
      } else {
         out <- ifelse(correct, object$fit.stats["AICc","ML"], object$fit.stats["AIC","ML"])
      }

   } else {

      # if there is 'object' and additional objects via ...

      if (object$method == "REML") {
         out <- sapply(list(object, ...), function(x) ifelse(correct, x$fit.stats["AICc","REML"], x$fit.stats["AIC","REML"]))
      } else {
         out <- sapply(list(object, ...), function(x) ifelse(correct, x$fit.stats["AICc","ML"], x$fit.stats["AIC","ML"]))
      }
      dfs <- sapply(list(object, ...), function(x) x$parms)

      out <- data.frame(df=dfs, AIC=out)

      if (correct)
         names(out)[2] <- "AICc"

      # get the names of the objects; same idea as in stats:::AIC.default

      cl <- match.call()
      cl$k <- NULL
      cl$correct <- NULL
      rownames(out) <- as.character(cl[-1L])

      # check that all models were fitted to the same data

      chksums <- sapply(list(object, ...), function(x) x$chksumyi)

      if (any(chksums[1] != chksums))
         warning(mstyle$warning("Models not all fitted to the same data."), call.=FALSE)

   }

   return(out)

}
