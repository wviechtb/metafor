############################################################################

### density of non-central hypergeometric distribution (based on Liao and Rosen, 2001) from MCMCpack
### Liao, J. G. & Rosen, O. (2001). Fast and stable algorithms for computing and sampling from the
### noncentral hypergeometric distribution. The American Statistician, 55, 366-369.

.dnoncenhypergeom <- function (x=NA_real_, n1, n2, m1, psi) { ### x=ai, n1=ai+bi, n2=ci+di, m1=ai+ci, psi=ORi

   mstyle <- .get.mstyle()

   mode.compute <- function(n1, n2, m1, psi, ll, uu) {
      a <- psi - 1
      b <- -((m1 + n1 + 2) * psi + n2 - m1)
      c <- psi * (n1 + 1) * (m1 + 1)
      q <- b + sign(b) * sqrt(b * b - 4 * a * c)
      q <- -q/2
      mode <- trunc(c/q)
      if (uu >= mode && mode >= ll)
         return(mode)
      else return(trunc(q/a))
   }
   r.function <- function(n1, n2, m1, psi, i) {
      (n1 - i + 1) * (m1 - i + 1)/i/(n2 - m1 + i) * psi
   }
   ll <- max(0, m1 - n2)
   uu <- min(n1, m1)

   if (n1 < 0 | n2 < 0)
      stop(mstyle$stop("'n1' or 'n2' negative in dnoncenhypergeom()."), call.=FALSE)

   if (m1 < 0 | m1 > (n1 + n2))
      stop(mstyle$stop("'m1' out of range in dnoncenhypergeom()."))

   if (psi <= 0)
      stop(mstyle$stop("'psi' [odds ratio] negative in dnoncenhypergeom()."), call.=FALSE)

   if (!is.na(x) & (x < ll | x > uu))
      stop(mstyle$stop("'x' out of bounds in dnoncenhypergeom()."))

   if (!is.na(x) & length(x) > 1L)
      stop(mstyle$stop("'x' neither missing or scalar in dnoncenhypergeom()."), call.=FALSE)

   mode <- mode.compute(n1, n2, m1, psi, ll, uu)
   pi <- array(1, uu - ll + 1)
   shift <- 1 - ll

   if (mode < uu) {
      r1 <- r.function(n1, n2, m1, psi, (mode + 1):uu)
      pi[(mode + 1 + shift):(uu + shift)] <- cumprod(r1)
   }
   if (mode > ll) {
      r1 <- 1/r.function(n1, n2, m1, psi, mode:(ll + 1))
      pi[(mode - 1 + shift):(ll + shift)] <- cumprod(r1)
   }
   pi <- pi/sum(pi)
   if (is.na(x)) {
      return(cbind(ll:uu, pi))
   } else {
      return(pi[x + shift])
   }

}

############################################################################

### density of non-central hypergeometric distribution for fixed- and random/mixed-effects models

.dnchgi <- function(logOR, ai, bi, ci, di, mu.i, tau2, random, dnchgcalc, dnchgprec) {

   mstyle <- .get.mstyle()

   k <- length(logOR)
   dnchgi <- rep(NA_real_, k)

   ### beyond these values, the results from dFNCHypergeo (from BiasedUrn package) become unstable

   pow <- 12

   logOR[logOR < log(10^-pow)] <- log(10^-pow)
   logOR[logOR > log(10^pow)]  <- log(10^pow)

   for (i in seq_len(k)) {

      ORi <- exp(logOR[i])

      if (dnchgcalc == "dnoncenhypergeom") {
         res <- try(.dnoncenhypergeom(x=ai, n1=ai+bi, n2=ci+di, m1=ai+ci, psi=ORi))
      } else {
         res <- try(BiasedUrn::dFNCHypergeo(x=ai, m1=ai+bi, m2=ci+di, n=ai+ci, odds=ORi, precision=dnchgprec))
      }

      if (inherits(res, "try-error")) {
         stop(mstyle$stop(paste0("Could not compute density of non-central hypergeometric distribution in study ", i, ".")), call.=FALSE)
      } else {
         dnchgi[i] <- res
      }

   }

   if (random)
      dnchgi <- dnchgi * dnorm(logOR, mu.i, sqrt(tau2))

   return(dnchgi)

}

############################################################################

### joint density of k non-central hypergeometric distributions for fixed- and random/mixed-effects models

.dnchg <- function(parms, ai, bi, ci, di, X.fit, random, verbose=FALSE, digits, dnchgcalc, dnchgprec, intCtrl) {

   mstyle <- .get.mstyle()

   p    <- ncol(X.fit)
   k    <- length(ai)
   beta <- parms[seq_len(p)]                  ### first p elemenets in parms are the model coefficients
   tau2 <- ifelse(random, exp(parms[p+1]), 0) ### next value is tau^2 -- optimize over exp(tau^2) value or hold at 0 if random=FALSE
   mu.i <- X.fit %*% cbind(beta)

   lli  <- rep(NA_real_, k)

   if (!random) {

      for (i in seq_len(k)) {
         lli[i] <- log(.dnchgi(logOR=mu.i[i], ai=ai[i], bi=bi[i], ci=ci[i], di=di[i], random=random, dnchgcalc=dnchgcalc, dnchgprec=dnchgprec))
      }

      if (verbose)
         cat(mstyle$verbose(paste("ll =", fmtx(sum(lli), digits[["fit"]]), " ", fmtx(beta, digits[["est"]]), "\n")))

   }

   if (random) {

      for (i in seq_len(k)) {

         res <- try(integrate(.dnchgi, lower=intCtrl$lower, upper=intCtrl$upper, ai=ai[i], bi=bi[i], ci=ci[i], di=di[i], mu.i=mu.i[i], tau2=tau2, random=random, dnchgcalc=dnchgcalc, dnchgprec=dnchgprec, rel.tol=intCtrl$rel.tol, subdivisions=intCtrl$subdivisions, stop.on.error=FALSE), silent=!verbose)

         if (inherits(res, "try-error")) {
            stop(mstyle$stop(paste0("Could not integrate over density of non-central hypergeometric distribution in study ", i, ".")), call.=FALSE)
         } else {
            if (res$value > 0) {
               lli[i] <- log(res$value)
            } else {
               lli[i] <- -Inf
            }
         }

      }

      if (verbose)
         cat(mstyle$verbose(paste("ll = ", fmtx(sum(lli), digits[["fit"]]), " ", fmtx(tau2, digits[["var"]]), " ", fmtx(beta, digits[["est"]]), "\n")))

   }

   return(-sum(lli))

}

############################################################################
