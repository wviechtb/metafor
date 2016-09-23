ranktest.rma <- function(x, ...) {

   #########################################################################

   if (!inherits(x, "rma"))
      stop("Argument 'x' must be an object of class \"rma\".")

   if (inherits(x, "robust.rma"))
      stop("Function not applicable to objects of class \"robust.rma\".")

   #########################################################################

   yi <- x$yi
   vi <- x$vi

   res <- rma.uni(yi, vi, method="FE")
   b   <- c(res$b)
   vb  <- c(res$vb)

   vi.star <- vi - vb
   yi.star <- (yi - b) / sqrt(vi.star)
   res <- cor.test(yi.star, vi, method="kendall", exact=TRUE)

   pval <- res$p.value
   tau  <- res$estimate

   res <- list(tau=tau, pval=pval, digits=x$digits)

   class(res) <- "ranktest.rma"
   return(res)

}
