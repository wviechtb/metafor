.hettest.esttau2i <- function(yi, vi, method, res0, mom, tau2i.init, threshold, maxiter) {

   if (mom) {

      hi <- hatvalues(res0)
      vari <- resid(res0)^2 / (1-hi)
      tau2i <- pmax(0, vari - vi)

   } else {

      k <- length(yi)
      tau2i <- tau2i.init
      diffs <- rep(10, k)
      conv <- 1
      iter <- 1
      while (any(diffs > threshold)) {
         if (iter > maxiter) {
            conv <- 0
            break
         }
         tau2i.old <- tau2i
         wi <- 1 / (tau2i + vi)
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

   }

   return(tau2i)

}

.hettest.lrt <- function(yi, vi, method, res0, tau2i) {

   k <- length(yi)
   vari <- tau2i + vi
   wi <- 1 / vari
   mu_hat <- sum(wi*yi) / sum(wi)

   if (method == "ML") {
      llsat <- sum(dnorm(yi, mean=mu_hat, sd=sqrt(vari), log=TRUE))
   } else {
      llsat <- sum(dnorm(yi, mean=mu_hat, sd=sqrt(vari), log=TRUE)) + 1/2 * log(2*k*pi) - 1/2 * log(sum(wi))
   }

   x2 <- c(-2 * (logLik(res0) - llsat))
   pval <- pchisq(x2, df=k-1, lower.tail=FALSE)

   return(list(statistic=x2, df=k-1, pval=pval))

}

.hettest.wald <- function(yi, vi, method, tau2i) {

   k <- length(yi)
   vari <- tau2i + vi

   if (method == "ML") {
      V.tau2 <- diag(2 * vari^2)
   } else {
      wi <- 1 / vari
      S <- sum(wi)
      I <- outer(X = 1 / (sqrt(2) * S * vari^2), Y = 1 / (sqrt(2) * S * vari^2))
      #diag(I) <- 1 / (2 * vari^2) - 1 / (S * vari^3) + 1 / (2 * S^2 * vari^4)
      #diag(I) <- (vari - 1/S)^2 / (2 * vari^4) # same
      diag(I) <- 1/2 * (wi - 1 / (S * vari^2))^2 # same
      V.tau2 <- try(solve(I), silent=TRUE)
      if (inherits(V.tau2, "try-error"))
         return(list(statistic=NA_real_, df=k-1, pval=NA_real_, tau2i=tau2i, se.tau2i=rep(NA_real_, k)))
   }
   se.tau2i <- sqrt(diag(V.tau2))
   X <- cbind(1, -diag(k-1))
   tau2i.diff <- c(X %*% tau2i)
   V.diff <- X %*% V.tau2 %*% t(X)
   W.diff <- try(solve(V.diff), silent=TRUE)

   if (inherits(W.diff, "try-error"))
      return(list(statistic=NA_real_, df=k-1, pval=NA_real_, tau2i=tau2i, se.tau2i=se.tau2i))

   x2 <- c(tau2i.diff %*% W.diff %*% tau2i.diff)
   pval <- pchisq(x2, df=k-1, lower.tail=FALSE)

   return(list(statistic=x2, df=k-1, pval=pval, se.tau2i=se.tau2i))

}

.hettest.score <- function(yi, vi, method, res0) {

   k <- length(yi)
   mu_hat <- res0$b[1]
   tau2_hat <- res0$tau2
   vari <- tau2_hat + vi
   wi <- 1 / vari

   if (method == "ML") {
      #Ui <- 1/2 * wi^2 * (yi - mu_hat)^2 - 1/2 * wi
      #Iii <- 1 / (2*vari^2)
      #x2 <- sum(Ui^2 / Iii)
      x2 <- 1/2 * sum((wi * (yi - mu_hat)^2 - 1)^2) # same
   } else {
      S <- sum(wi)
      Ui <- 1/2 * wi^2 * (yi - mu_hat)^2 - 1/2 * wi + 1/2 * wi^2 / S
      #Iii <- -1 * (1/2 * wi^2 - wi^2 / S + wi^4 / S^2 - (wi^2 - wi^3 / S))
      #Iii <- 1/2 * (wi - 1 / (S * vari^2))^2
      #x2 <- sum(Ui^2 / Iii) # not the same!
      I <- outer(X = 1 / (sqrt(2) * S * vari^2), Y = 1 / (sqrt(2) * S * vari^2))
      diag(I) <- 1/2 * (wi - 1 / (S * vari^2))^2
      invI <- try(solve(I), silent=TRUE)
      if (inherits(invI, "try-error"))
         return(list(statistic=NA_real_, df=k-1, pval=NA_real_))
      x2 <- c(t(Ui) %*% invI %*% Ui)
   }

   pval <- pchisq(x2, df=k-1, lower.tail=FALSE)

   return(list(statistic=x2, df=k-1, pval=pval))

}

.hettest.ks <- function(x, cdf, ...) {
   k <- length(x)
   x <- sort(x)
   Fx <- cdf(x, ...)
   Dplus  <- max((1:k)/k - Fx)
   Dminus <- max(Fx - (0:(k-1))/k)
   Dval <- max(Dplus, Dminus)
   return(list(statistic=Dval, df=NA_integer_, pval=NA_real_))
}

.hettest.ad <- function(x, cdf, ...) {
   k <- length(x)
   x <- sort(x)
   Fx <- cdf(x, ...)
   eps <- .Machine$double.eps
   Fx <- pmin(pmax(Fx, eps), 1 - eps)
   A2 <- -k - mean((2*(1:k)-1) * (log(Fx) + log(1 - rev(Fx))))
   return(list(statistic=A2, df=NA_integer_, pval=NA_real_))
}
