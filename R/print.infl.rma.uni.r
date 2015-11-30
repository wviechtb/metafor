print.infl.rma.uni <- function(x, digits, ...) {

   if (class(x) != "infl.rma.uni")
      stop("Argument 'x' must be an object of class \"infl.rma.uni\".")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop("Unknown 'na.action' specified under options().")

   if (missing(digits))
      digits <- x$digits

   #########################################################################

   inf <- round(x$inf, digits)
   dfbs <- round(x$dfbs, digits)

   inf$inf <- ifelse(!is.na(x$is.infl) & x$is.infl, "*", "")

   ### check for NAs and act accordingly

   any.na <- apply(is.na(inf), 1, any) | apply(is.na(dfbs), 1, any)

   if (any(any.na)) {

      #x$not.na <- !any.na ### this will omit even newly generated NAs (due to unfittable model or non-convergence)

      if (na.act == "na.omit") {
         inf  <- inf[x$not.na,]
         dfbs <- dfbs[x$not.na,,drop=FALSE] ### need drop=FALSE for intercept-only models
         out  <- list(inf=inf, dfbs=dfbs)
      }

      if (na.act == "na.exclude" || na.act == "na.pass") {
         out <- list(inf=inf, dfbs=dfbs)
      }

      if (na.act == "na.fail")
         stop("Missing values in results.")

   } else {
      out <- list(inf=inf, dfbs=dfbs)
   }

   #########################################################################

   ### if there is only one fixed effect, put output in a single table

   if (x$p == 1) {
      out <- cbind(inf[-ncol(inf)], dfbs, inf[ncol(inf)])
      colnames(out)[ncol(out)-1] <- "dfbs"
   }

   #########################################################################

   print(out)

}
