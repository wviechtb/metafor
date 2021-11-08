ranktest.rma <- function(x, digits, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="rma", notav=c("robust.rma", "rma.uni.selmodel"))

   if (missing(digits)) {
      digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
   }

   ddd <- list(...)

   .chkdots(ddd, c("exact"))

   if (is.null(ddd$exact)) {
      exact <- TRUE
   } else {
      exact <- ddd$exact
   }

   #########################################################################

   yi <- x$yi
   vi <- x$vi

   res  <- rma.uni(yi, vi, method="EE")
   beta <- c(res$beta)
   vb   <- c(res$vb)

   vi.star <- vi - vb
   yi.star <- (yi - beta) / sqrt(vi.star)
   res <- cor.test(yi.star, vi, method="kendall", exact=exact)

   pval <- res$p.value
   tau  <- res$estimate

   res <- list(tau=tau, pval=pval, digits=digits)

   class(res) <- "ranktest"
   return(res)

}
