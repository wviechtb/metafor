.hettest.scoretest <- function(yi, vi, mu_hat, tau2_hat, method) {

   k <- length(yi)

   vari <- vi + tau2_hat
   wi <- 1 / vari

   if (method == "REML") {
      sumwi <- sum(wi)
      Si <- -1/2 * wi + 1/2 * wi^2 / sumwi + 1/2 * wi^2 * (yi - mu_hat)^2
      varSi <- -1 * (1/2 * wi^2 - wi^2 / sumwi + wi^4 / sumwi^2 - (wi^2 - wi^3 / sumwi))
      S <- sum(Si^2 / varSi)
   } else {
      Si <- -1/2 * wi + 1/2 * wi^2 * (yi - mu_hat)^2
      varSi <- 1 / (2*vari^2)
      S <- sum(Si^2 / varSi)
   }

   pval <- pchisq(S, df=k-1, lower.tail=FALSE)

   return(list(statistic=S, df=k-1, pval=pval))

}

.hettest.lrt <- function(yi, vi, method, res0, mom) {

   k <- length(yi)

   if (mom) {
      hi <- hatvalues(res0)
      vari <- resid(res0)^2 / (1-hi)
      tau2i <- pmax(0, vari - vi)
   } else {
      tau2i <- try(.hettest.esttau2i(yi, vi, method=method), silent=TRUE)
   }

   if (inherits(tau2i, "try-error"))
      stop(mstyle$stop("Could not estimate the tau^2_i values."))

   res1 <- rma(yi, vi=vi+tau2i, method=method, tau2=0)

   x2 <- c(-2 * (logLik(res0) - logLik(res1)))
   pval <- pchisq(x2, df=k-1, lower.tail=FALSE)
   #pval <- .pchisqmix(x2, df=k)

   return(list(statistic=x2, df=k-1, pval=pval, tau2i=tau2i))

}

.hettest.esttau2i <- function(yi, vi, method) {

   k <- length(yi)
   tau2i <- rep(0.2, k)
   diffs <- rep(1, k)
   conv <- 1
   iter <- 1
   while (any(diffs > 10^-5)) {
      if (iter > 500) {
         break
         conv <- 0
      }
      tau2i.old <- tau2i
      wi <- 1/(tau2i + vi)
      sumwi <- sum(wi)
      mu_hat <- sum(wi*yi) / sumwi
      if (method == "ML") {
         tau2i <- c(yi - mu_hat)^2 - vi
      } else {
         tau2i <- c(yi - mu_hat)^2 - vi + 1 / sumwi
      }
      tau2i[tau2i < 0] <- 0
      diffs <- abs(tau2i - tau2i.old)
      iter <- iter + 1
   }

   if (conv == 0)
      stop()

   return(tau2i)

}

.pchisqmix <- function(x2, df, tol = 1e-12) {
   j <- 1:df
   weights <- dbinom(j, df, 0.5)
   p <- sum(weights * pchisq(x2, df=j, lower.tail=FALSE))
   if (x2 < tol)
      p <- p + dbinom(0, df, 0.5)
   return(p)
}
