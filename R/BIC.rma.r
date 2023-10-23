BIC.rma <- function(object, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(object), must="rma")

   if (missing(...)) {

      ### if there is just 'object'

      if (object$method == "REML") {
         out <- object$fit.stats["BIC","REML"]
      } else {
         out <- object$fit.stats["BIC","ML"]
      }

   } else {

      ### if there is 'object' and additional objects via ...

      if (object$method == "REML") {
         out <- sapply(list(object, ...), function(x) x$fit.stats["BIC","REML"])
      } else {
         out <- sapply(list(object, ...), function(x) x$fit.stats["BIC","ML"])
      }
      dfs <- sapply(list(object, ...), function(x) x$parms)

      out <- data.frame(df=dfs, BIC=out)

      ### get names of objects; same idea as in stats:::AIC.default

      cl <- match.call()
      rownames(out) <- as.character(cl[-1L])

      ### check that all models were fitted to the same data

      yis <- lapply(list(object, ...), function(x) as.vector(x$yi))

      if (!all(sapply(yis[-1], function(x) identical(x, yis[[1]]))))
         warning(mstyle$warning("Models not all fitted to the same data."), call.=FALSE)

   }

   return(out)

}
