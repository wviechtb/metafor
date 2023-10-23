############################################################################

### c(m) calculation function for bias correction of SMDs (mi = n1i + n2i - 2) or SMCC/SMCRs (mi = ni - 1)

.cmicalc <- function(mi, correct=TRUE) {

   ### this can overflow if mi is 'large' (on my machine, if mi >= 344)
   #cmi <- gamma(mi/2)/(sqrt(mi/2)*gamma((mi-1)/2))
   ### catch those cases and apply the approximate formula (which is accurate then)
   #is.na <- is.na(cmi)
   #cmi[is.na] <- 1 - 3/(4*mi[is.na] - 1)

   if (correct) {
      ### this avoids the problem with overflow altogether
      cmi <- ifelse(mi <= 1, NA_real_, exp(lgamma(mi/2) - log(sqrt(mi/2)) - lgamma((mi-1)/2)))
   } else {
      cmi <- rep(1, length(mi))
   }
   return(cmi)

}

############################################################################

### function to compute the tetrachoric correlation coefficient and its sampling variance

.rtet <- function(ai, bi, ci, di, maxcor=.9999) {

   mstyle <- .get.mstyle()

   if (!requireNamespace("mvtnorm", quietly=TRUE))
      stop(mstyle$stop("Please install the 'mvtnorm' package to compute this measure."), call.=FALSE)

   fn <- function(par, ai, bi, ci, di, maxcor, fixcut=FALSE) {

      rho <- par[1]
      cut.row <- par[2]
      cut.col <- par[3]

      ### truncate rho values outside of specified bounds
      if (abs(rho) > maxcor)
         rho <- sign(rho) * maxcor

      ### to substitute fixed cut values
      if (fixcut) {
         cut.row <- qnorm((ai+bi)/ni)
         cut.col <- qnorm((ai+ci)/ni)
      }

      #       │ ci | di    #   ci = lo X and hi Y    di = hi X and hi Y
      # var Y │----+----   #
      #       │ ai | bi    #   ai = lo X and lo Y    bi = hi X and lo Y
      #       ┼─────────
      #          var X
      #
      #      lo   hi
      #    +----+----+
      # lo | ai | bi |
      #    +----+----+ var Y
      # hi | ci | di |
      #    +----+----+
      #       var X

      R <- matrix(c(1,rho,rho,1), nrow=2, ncol=2)

      p.ai <- mvtnorm::pmvnorm(lower=c(-Inf,-Inf), upper=c(cut.col,cut.row), corr=R)
      p.bi <- mvtnorm::pmvnorm(lower=c(cut.col,-Inf), upper=c(+Inf,cut.row), corr=R)
      p.ci <- mvtnorm::pmvnorm(lower=c(-Inf,cut.row), upper=c(cut.col,+Inf), corr=R)
      p.di <- mvtnorm::pmvnorm(lower=c(cut.col,cut.row), upper=c(+Inf,+Inf), corr=R)

      ### in principle, should be able to compute these values with the following code, but this
      ### leads to more numerical instabilities when optimizing (possibly due to negative values)
      #p.y.lo <- pnorm(cut.row)
      #p.x.lo <- pnorm(cut.col)
      #p.ai <- mvtnorm::pmvnorm(lower=c(-Inf,-Inf), upper=c(cut.col,cut.row), corr=R)
      #p.bi <- p.y.lo - p.ai
      #p.ci <- p.x.lo - p.ai
      #p.di <- 1 - p.ai - p.bi - p.ci

      if (any(p.ai <= 0 || p.bi <= 0 || p.ci <= 0 || p.di <= 0)) {
         ll <- -Inf
      } else {
         ll <- ai*log(p.ai) + bi*log(p.bi) + ci*log(p.ci) + di*log(p.di)
      }

      return(-ll)

   }

   ni <- ai + bi + ci + di

   ### if one of the margins is equal to zero, then r_tet could in principle be equal to any value,
   ### but we define it here to be zero (presuming independence until evidence of dependence is found)
   ### but with infinite variance
   if ((ai + bi) == 0L || (ci + di) == 0L || (ai + ci) == 0L || (bi + di) == 0L)
      return(list(yi=0, vi=Inf))

   ### if bi and ci is zero, then r_tet must be +1 with zero variance
   if (bi == 0L && ci == 0L)
      return(list(yi=1, vi=0))

   ### if ai and di is zero, then r_tet must be -1 with zero variance
   if (ai == 0L && di == 0L)
      return(list(yi=-1, vi=0))

   ### cases where only one cell is equal to zero are handled further below

   ### in all other cases, first optimize over rho with cut values set to the sample values
   ### use suppressWarnings() to suppress "NA/Inf replaced by maximum positive value" warnings
   res <- try(suppressWarnings(optimize(fn, interval=c(-1,1), ai=ai, bi=bi, ci=ci, di=di, maxcor=maxcor, fixcut=TRUE)), silent=TRUE)

   ### check for non-convergence
   if (inherits(res, "try-error")) {
      warning(mstyle$warning("Could not estimate tetrachoric correlation coefficient."), call.=FALSE)
      return(list(yi=NA, vi=NA))
   }

   ### then use the value as the starting point and maximize over rho and the cut values
   ### (Nelder-Mead seems to do fine here; using L-BFGS-B doesn't seem to improve on this)
   res <- try(optim(par=c(res$minimum,qnorm((ai+bi)/ni),qnorm((ai+ci)/ni)), fn, ai=ai, bi=bi, ci=ci, di=di, maxcor=maxcor, fixcut=FALSE, hessian=TRUE), silent=TRUE)
   #res <- try(optim(par=c(res$minimum,qnorm((ai+bi)/ni),qnorm((ai+ci)/ni)), fn, method="L-BFGS-B", lower=c(-1,-Inf,-Inf), upper=c(1,Inf,Inf), ai=ai, bi=bi, ci=ci, di=di, maxcor=maxcor, fixcut=FALSE, hessian=TRUE), silent=TRUE)

   ### check for non-convergence
   if (inherits(res, "try-error")) {
      warning(mstyle$warning("Could not estimate tetrachoric correlation coefficient."), call.=FALSE)
      return(list(yi=NA, vi=NA))
   }

   ### take inverse of hessian and extract variance for estimate
   ### (using hessian() seems to lead to more problems, so stick with hessian from optim())
   vi <- try(chol2inv(chol(res$hessian))[1,1], silent=TRUE)
   #res$hessian <- try(chol2inv(chol(numDeriv::hessian(fn, x=res$par, ai=ai, bi=bi, ci=ci, di=di, maxcor=maxcor, fixcut=FALSE))), silent=TRUE)

   ### check for problems with computing the inverse
   if (inherits(vi, "try-error")) {
      warning(mstyle$warning("Could not estimate sampling variance of tetrachoric correlation coefficient."), call.=FALSE)
      vi <- NA
   }

   ### extract estimate
   yi <- res$par[1]

   ### but if bi or ci is zero, then r_tet must be +1
   if (bi == 0 || ci == 0)
      yi <- 1

   ### but if ai or di is zero, then r_tet must be -1
   if (ai == 0 || di == 0)
      yi <- -1

   ### note: what is the right variance when there is one zero cell?
   ### vi as estimated gets smaller as the table becomes more and more like
   ### a table with 0 diagonal/off-diagonal, which intuitively makes sense

   ### return estimate and sampling variance (and SE)
   return(list(yi=yi, vi=vi, sei=sqrt(vi)))

   ### Could consider implementing the Fisher scoring algorithm; first derivatives and
   ### elements of the information matrix are given in Tallis (1962). Could also consider
   ### estimating the variance from the inverse of the information matrix. But constructing
   ### the information matrix takes a bit of extra work and it is not clear to me how to
   ### handle estimated cell probabilities that go to zero here.

}

############################################################################

### function to calculate the Gaussian hypergeometric (Hypergeometric2F1) function

.Fcalc <- function(a, b, g, x) {

   mstyle <- .get.mstyle()

   if (!requireNamespace("gsl", quietly=TRUE))
      stop(mstyle$stop("Please install the 'gsl' package to use measure='UCOR'."), call.=FALSE)

   k.g <- length(g)
   k.x <- length(x)
   k   <- max(k.g, k.x)

   res <- rep(NA_real_, k)

   if (k.g == 1)
      g <- rep(g, k)
   if (k.x == 1)
      x <- rep(x, k)

   if (length(g) != length(x))
      stop(mstyle$stop("Length of 'g' and 'x' arguments is not the same."))

   for (i in seq_len(k)) {

      if (!is.na(g[i]) && !is.na(x[i]) && g[i] > (a+b)) {
         res[i] <- gsl::hyperg_2F1(a, b, g[i], x[i])
      } else {
         res[i] <- NA
      }

   }

   return(res)

}

############################################################################

### pdf of SMD (with or without bias correction)

.dsmd <- function(x, n1, n2, theta, correct=TRUE, xisg=FALSE, warn=FALSE) {

   nt <- n1 * n2 / (n1 + n2)
   m  <- n1 + n2 - 2

   cm <- .cmicalc(m)

   if (xisg)
      x <- x / cm

   if (!correct)
      cm <- 1

   if (warn) {
      res <- dt(x * sqrt(nt) / cm, df = m, ncp = sqrt(nt) * theta) * sqrt(nt) / cm
   } else {
      res <- suppressWarnings(dt(x * sqrt(nt) / cm, df = m, ncp = sqrt(nt) * theta) * sqrt(nt) / cm)
   }

   return(res)

}

#integrate(function(x) .dsmd(x, n1=4, n2=4, theta=.5), lower=-Inf, upper=Inf)
#integrate(function(x) x*.dsmd(x, n1=4, n2=4, theta=.5), lower=-Inf, upper=Inf)

### pdf of COR

.dcor <- function(x, n, rho) {

   x[x < -1] <- NA
   x[x >  1] <- NA

   ### only accurate for n >= 5
   n[n <= 4] <- NA

   ### calculate density
   res <- exp(log(n-2) + lgamma(n-1) + (n-1)/2 * log(1 - rho^2) + (n-4)/2 * log(1 - x^2) -
          1/2 * log(2*base::pi) - lgamma(n-1/2) - (n-3/2) * log(1 - rho*x)) *
          .Fcalc(1/2, 1/2, n-1/2, (rho*x + 1)/2)

   ### make sure that density is 0 for r = +-1
   res[abs(x) == 1] <- 0

   return(res)

}

#integrate(function(x) .dcor(x, n=5, rho=.8), lower=-1, upper=1)
#integrate(function(x) x*.dcor(x, n=5, rho=.8), lower=-1, upper=1) ### should not be rho due to bias!
#integrate(function(x) x*.Fcalc(1/2, 1/2, (5-2)/2, 1-x^2)*.dcor(x, n=5, rho=.8), lower=-1, upper=1) ### should be ~rho

### pdf of ZCOR

.dzcor <- function(x, n, rho, zrho) {

   ### only accurate for n >= 5
   n[n <= 4] <- NA

   ### if rho is missing, then back-transform zrho value(s)
   if (missing(rho))
      rho <- tanh(zrho)

   ### copy x to z and back-transform z values (so x = correlation)
   z <- x
   x <- tanh(z)

   ### calculate density
   res <- exp(log(n-2) + lgamma(n-1) + (n-1)/2 * log(1 - rho^2) + (n-4)/2 * log(1 - x^2) -
          1/2 * log(2*base::pi) - lgamma(n-1/2) - (n-3/2) * log(1 - rho*x) +
          log(4) + 2*z - 2*log(exp(2*z) + 1)) *
          .Fcalc(1/2, 1/2, n-1/2, (rho*x + 1)/2)

   ### make sure that density is 0 for r = +-1
   res[abs(x) == 1] <- 0

   return(res)

}

#integrate(function(x) .dzcor(x, n=5, rho=.8), lower=-100, upper=100)
#integrate(function(x) x*.dzcor(x, n=5, rho=.8), lower=-100, upper=100)

### pdf of ARAW

.daraw <- function(x, n, m, alpha) {
   res <- df((1-x)/(1-alpha), (n-1)*(m-1), (n-1)) / (1-alpha)
   res[alpha >=  1] <- 0
   res[alpha <= -1] <- 0
   return(res)
}

#integrate(function(x) .daraw(x, n=10, m=2, alpha=.8), lower=-Inf, upper=Inf)
#integrate(function(x) x*.daraw(x, n=10, m=2, alpha=.8), lower=-Inf, upper=Inf)

############################################################################

### function to convert p-values to t-statistics (need this to catch NULL
### since sign(NULL) and qt(NULL) throw errors)

.convp2t <- function(pval, df) {

   if (is.null(pval))
      return(NULL)

   df <- ifelse(df < 1, NA, df)
   pval <- ifelse(abs(pval) > 1, NA, pval)

   sign(pval) * qt(abs(pval)/2, df=df, lower.tail=FALSE)

}

### function to convert p-values to F-statistics (need this to catch NULL
### since qf(NULL) throws an error)

.convp2f <- function(pval, df1, df2) {

   if (is.null(pval))
      return(NULL)

   df1 <- ifelse(df1 < 1, NA, df1)
   df2 <- ifelse(df2 < 1, NA, df2)
   pval <- ifelse(pval < 0, NA, pval)
   pval <- ifelse(pval > 1, NA, pval)

   qf(pval, df1=df1, df2=df2, lower.tail=FALSE)

}

############################################################################
