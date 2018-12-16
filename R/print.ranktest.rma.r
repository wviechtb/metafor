print.ranktest.rma <- function(x, digits, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(x, "ranktest.rma"))
      stop(mstyle$stop("Argument 'x' must be an object of class \"ranktest.rma\"."))

   if (missing(digits))
      digits <- x$digits

   if (!exists(".rmspace"))
      cat("\n")

   cat(mstyle$section("Rank Correlation Test for Funnel Plot Asymmetry"))
   cat("\n\n")
   cat(mstyle$result(paste0("Kendall's tau = ", formatC(x$tau, digits=digits, format="f"), ", p ", .pval(x$pval, digits=digits, showeq=TRUE, sep=" "))))
   cat("\n")
   #cat("H0: true tau is equal to 0\n\n")

   if (!exists(".rmspace"))
      cat("\n")

   invisible()

}
