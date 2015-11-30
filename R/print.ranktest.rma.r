print.ranktest.rma <- function(x, digits, ...) {

   if (class(x) != "ranktest.rma")
      stop("Argument 'x' must be an object of class \"ranktest.rma\".")

   if (missing(digits))
      digits <- x$digits

   cat("\n")
   cat("Rank Correlation Test for Funnel Plot Asymmetry\n\n")
   cat("Kendall's tau = ", formatC(x$tau, digits=digits, format="f"), ", p ", .pval(x$pval, digits=digits, showeq=TRUE, sep=" "), "\n\n", sep="")
   #cat("H0: true tau is equal to 0\n\n")

   invisible()

}
