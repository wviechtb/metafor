.tes.intfun <- function(x, theta, tau, sei, H0, alternative, crit) {
   if (alternative == "two.sided")
      pow <- (pnorm(crit, mean=(x-H0)/sei, sd=1, lower.tail=FALSE) + pnorm(-crit, mean=(x-H0)/sei, sd=1, lower.tail=TRUE))
   if (alternative == "greater")
      pow <- pnorm(crit, mean=(x-H0)/sei, sd=1, lower.tail=FALSE)
   if (alternative == "less")
      pow <- pnorm(crit, mean=(x-H0)/sei, sd=1, lower.tail=TRUE)
   res <- pow * dnorm(x, theta, tau)
   return(res)
}

.tes.lim <- function(theta, yi, vi, H0, alternative, alpha, tau2, test, tes.alternative, progbar, tes.alpha, correct, rel.tol, subdivisions, tau2.lb) {
   pval <- tes(x=yi, vi=vi, H0=H0, alternative=alternative, alpha=alpha, theta=theta, tau2=tau2, test=test, tes.alternative=tes.alternative, progbar=progbar,
                       tes.alpha=tes.alpha, correct=correct, rel.tol=rel.tol, subdivisions=subdivisions, tau2.lb=tau2.lb, find.lim=FALSE)$pval
   #cat("theta = ", theta, " pval = ", pval, "\n")
   return(pval - tes.alpha)
}
