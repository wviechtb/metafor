print.list.confint.rma <- function(x, digits=x$digits, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="list.confint.rma")

   digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)

   x$digits <- NULL # so length(x) is correct

   if (!exists(".rmspace"))
      cat("\n")

   len <- length(x)

   for (j in 1:len) {

      res.random <- .fcf(x[[j]]$random, digits[["var"]])
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
