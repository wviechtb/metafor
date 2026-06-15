.hettest.scoretest <- function(yi, vi, mu_hat, tau2_hat, method) {

   k <- length(yi)

   vari <- tau2_hat + vi
   wi <- 1 / vari

   if (method == "REML") {
      S <- sum(wi)
      Ui <- 1/2 * wi^2 * (yi - mu_hat)^2 - 1/2 * wi + 1/2 * wi^2 / S
      #Iii <- -1 * (1/2 * wi^2 - wi^2 / S + wi^4 / S^2 - (wi^2 - wi^3 / S))
      #Iii <- 1/2 * (wi - 1 / (S * vari^2))^2
      #x2 <- sum(Ui^2 / Iii) # not the same!
      I <- outer(X = 1 / (sqrt(2) * S * vari^2), Y = 1 / (sqrt(2) * S * vari^2))
      diag(I) <- 1/2 * (wi - 1 / (S * vari^2))^2
      invI <- try(solve(I), silent=TRUE)
      if (inherits(invI, "try-error")) {
         mstyle <- .get.mstyle()
         stop(mstyle$stop("Could not invert the informatrion matrix."), call.=FALSE)
      }
      x2 <- c(t(Ui) %*% invI %*% Ui)
   } else {
      #Ui <- 1/2 * wi^2 * (yi - mu_hat)^2 - 1/2 * wi
      #Iii <- 1 / (2*vari^2)
      #x2 <- sum(Ui^2 / Iii)
      x2 <- 1/2 * sum((wi * (yi - mu_hat)^2 - 1)^2) # same
   }

   pval <- pchisq(x2, df=k-1, lower.tail=FALSE)

   return(list(statistic=x2, df=k-1, pval=pval))

}

.hettest.lrt <- function(yi, vi, method, res0, mom, tau2i.init, threshold, maxiter) {

   k <- length(yi)

   if (mom) {
      hi <- hatvalues(res0)
      vari <- resid(res0)^2 / (1-hi)
      tau2i <- pmax(0, vari - vi)
   } else {
      tau2i <- try(.hettest.esttau2i(yi, vi, method=method, tau2i.init=tau2i.init, threshold=threshold, maxiter=maxiter), silent=TRUE)
   }

   if (inherits(tau2i, "try-error")) {
      mstyle <- .get.mstyle()
      stop(mstyle$stop("Could not estimate the tau^2_i values."), call.=FALSE)
   }

   res1 <- rma(yi, vi=vi+tau2i, method=method, tau2=0)

   x2 <- c(-2 * (logLik(res0) - logLik(res1)))
   pval <- pchisq(x2, df=k-1, lower.tail=FALSE)
   #pval <- .pchisqmix(x2, df=k)

   return(list(statistic=x2, df=k-1, pval=pval, tau2i=tau2i))

}

.hettest.wald <- function(yi, vi, method, res0, mom, tau2i.init, threshold, maxiter) {

   k <- length(yi)

   if (mom) {
      hi <- hatvalues(res0)
      vari <- resid(res0)^2 / (1-hi)
      tau2i <- pmax(0, vari - vi)
   } else {
      tau2i <- try(.hettest.esttau2i(yi, vi, method=method, tau2i.init=tau2i.init, threshold=threshold, maxiter=maxiter), silent=TRUE)
   }

   if (inherits(tau2i, "try-error")) {
      mstyle <- .get.mstyle()
      stop(mstyle$stop("Could not estimate the tau^2_i values."), call.=FALSE)
   }

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
         V.tau2 <- matrix(NA_real_, nrow=k, ncol=k)
   }
   se.tau2i <- sqrt(diag(V.tau2))
   X <- cbind(1, -diag(k-1))
   tau2i.diff <- c(X %*% tau2i)
   V.diff <- X %*% V.tau2 %*% t(X)
   W.diff <- try(solve(V.diff), silent=TRUE)
   if (inherits(W.diff, "try-error")) {
      x2 <- NA
   } else {
      x2 <- c(tau2i.diff %*% W.diff %*% tau2i.diff)
   }

   pval <- pchisq(x2, df=k-1, lower.tail=FALSE)

   return(list(statistic=x2, df=k-1, pval=pval, tau2i=tau2i, se.tau2i=se.tau2i))

}

.hettest.esttau2i <- function(yi, vi, method, tau2i.init, threshold, maxiter) {

   k <- length(yi)
   tau2i <- tau2i.init
   diffs <- rep(1, k)
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

   return(tau2i)

}

.ks.test <- function(x, cdf, ...) {
   k <- length(x)
   x <- sort(x)
   Fx <- cdf(x, ...)
   Dplus  <- max((1:k)/k - Fx)
   Dminus <- max(Fx - (0:(k-1))/k)
   max(Dplus, Dminus)
}

.ad.test <- function(x, cdf, ...) {
   k <- length(x)
   x <- sort(x)
   Fx <- cdf(x, ...)
   eps <- .Machine$double.eps
   Fx <- pmin(pmax(Fx, eps), 1 - eps)
   A2 <- -k - mean((2*(1:k)-1) * (log(Fx) + log(1 - rev(Fx))))
   A2
}
