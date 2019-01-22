print.list.confint.rma <- function(x, digits, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(x, "list.confint.rma"))
      stop(mstyle$stop("Argument 'x' must be an object of class \"list.confint.rma\"."))

   if (missing(digits))
      digits <- x[[1]]$digits

   if (!exists(".rmspace"))
      cat("\n")

   len <- length(x)

   for (j in 1:len) {

      res.random <- formatC(x[[j]]$random, digits=digits, format="f")
      res.random[,2] <- paste0(x[[j]]$lb.sign, res.random[,2])
      res.random[,3] <- paste0(x[[j]]$ub.sign, res.random[,3])
      tmp <- capture.output(print(res.random, quote=FALSE, right=TRUE))
      .print.table(tmp, mstyle)

      if (j != len)
         cat("\n")

   }

   if (!exists(".rmspace"))
      cat("\n")

   invisible()

}
