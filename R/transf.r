transf.rtoz <- function(xi, ...)                   ### resulting value between -Inf (for -1) and +Inf (for +1)
   atanh(xi)

transf.ztor <- function(xi, ...)
   tanh(xi)

transf.logit <- function(xi, ...)                  ### resulting value between -Inf (for 0) and +Inf (for +1)
   qlogis(xi)

transf.ilogit <- function(xi, ...)
   plogis(xi)

transf.arcsin <- function(xi, ...)                 ### resulting value between 0 (for 0) and asin(1) = pi/2 (for +1)
   asin(sqrt(xi))

transf.iarcsin <- function(xi, ...) {
   zi <- sin(xi)^2
   zi[xi < 0] <- 0                                 ### if xi value is below 0 (e.g., CI bound), return 0
   zi[xi > asin(1)] <- 1                           ### if xi value is above maximum possible value, return 1
   return(c(zi))
}

transf.pft <- function(xi, ni, ...) {              ### Freeman-Tukey transformation for proportions
   xi <- xi*ni
   zi <- 1/2*(asin(sqrt(xi/(ni+1))) + asin(sqrt((xi+1)/(ni+1))))
   return(c(zi))
}

transf.ipft <- function(xi, ni, ...) {             ### inverse of Freeman-Tukey transformation for individual proportions
   zi <- suppressWarnings(1/2 * (1 - sign(cos(2*xi)) * sqrt(1 - (sin(2*xi)+(sin(2*xi)-1/sin(2*xi))/ni)^2)))
   zi <- ifelse(is.nan(zi), NA, zi)
   zi[xi > transf.pft(1,ni)] <- 1                  ### if xi is above upper limit, return 1
   zi[xi < transf.pft(0,ni)] <- 0                  ### if xi is below lower limit, return 0
   return(c(zi))
}

transf.ipft.hm <- function(xi, targs, ...) {       ### inverse of Freeman-Tukey transformation for a collection of proportions
   if (is.null(targs) || (is.list(targs) && is.null(targs$ni)))
      stop("Need to specify the sample sizes via the 'targs' argument.", call.=FALSE)
   if (is.list(targs)) {
      ni <- targs$ni
   } else {
      ni <- ni
   }
   nhm <- 1/(mean(1/ni, na.rm=TRUE))               ### calculate harmonic mean of the ni's
   zi <- suppressWarnings(1/2 * (1 - sign(cos(2*xi)) * sqrt(1 - (sin(2*xi)+(sin(2*xi)-1/sin(2*xi))/nhm)^2)))
   zi <- ifelse(is.nan(zi), NA, zi)                ### it may not be possible to calculate zi
   zi[xi > transf.pft(1,nhm)] <- 1                 ### if xi is above upper limit, return 1
   zi[xi < transf.pft(0,nhm)] <- 0                 ### if xi is below lower limit, return 0
   return(c(zi))
}

transf.isqrt <- function(xi, ...) {
   zi <- xi*xi
   zi[xi < 0] <- 0                                 ### if xi value is below 0 (e.g., CI bound), return 0
   return(c(zi))
}

transf.irft <- function(xi, ti, ...) {             ### Freeman-Tukey transformation for incidence rates
   zi <- 1/2*(sqrt(xi) + sqrt(xi + 1/ti))          ### xi is the incidence rate (not the number of events!)
   return(c(zi))
}

transf.iirft <- function(xi, ti, ...) {            ### inverse of Freeman-Tukey transformation for incidence rates (see Freeman-Tukey_incidence.r in code directory)
   #zi <- (1/ti - 2*xi^2 + ti*xi^4)/(4*xi^2*ti)    ### old version where transf.irft was not multiplied by 1/2
   zi <- (1/ti - 8*xi^2 + 16*ti*xi^4)/(16*xi^2*ti) ### xi is the incidence rate (not the number of events!)
   zi <- ifelse(is.nan(zi), NA, zi)
   zi[xi < transf.irft(0,ti)] <- 0                 ### if xi is below lower limit, return 0
   zi[zi <= .Machine$double.eps] <- 0              ### avoid finite precision errors in back-transformed values (transf.iirft(transf.irft(0, 1:200), 1:200))
   return(c(zi))
}

transf.ahw <- function(xi, ...) {                  ### resulting value between 0 (for alpha=0) and 1 (for alpha=1)
   #zi <- (1-xi)^(1/3)
   zi <- 1 - (1-xi)^(1/3)
   return(c(zi))
}

transf.iahw <- function(xi, ...) {
   #zi <- 1-xi^3
   zi <- 1 - (1-xi)^3
   zi <- ifelse(is.nan(zi), NA, zi)
   zi[xi > 1] <- 1                                 ### if xi is above upper limit, return 1
   zi[xi < 0] <- 0                                 ### if xi is below lower limit, return 0
   return(c(zi))
}

transf.abt <- function(xi, ...) {                  ### Bonett (2002) transformation of alphas (without bias correction)
#transf.abt <- function(xi, ni, ...) {             ### resulting value between 0 (for alpha=0) to Inf (for alpha=1)
   #zi <- log(1-xi) - log(ni/(ni-1))
   #zi <- log(1-xi)
   zi <- -log(1-xi)
   return(c(zi))
}

transf.iabt <- function(xi, ...) {                 ### inverse of Bonett (2002) transformation
#transf.iabt <- function(xi, ni, ...) {
   #zi <- 1 - exp(xi) * ni / (ni-1)
   #zi <- 1 - exp(xi)
   zi <- 1 - exp(-xi)
   zi <- ifelse(is.nan(zi), NA, zi)
   zi[xi < 0] <- 0                                 ### if xi is below lower limit, return 0
   return(c(zi))
}

transf.ztor.int <- function(xi, targs=NULL, ...) {

   if (is.null(targs$tau2))
      targs$tau2 <- 0
   if (is.null(targs$lower))
      targs$lower <- xi-5*sqrt(targs$tau2)
   if (is.null(targs$upper))
      targs$upper <- xi+5*sqrt(targs$tau2)

   toint <- function(zval, xi, tau2)
      tanh(zval) * dnorm(zval, mean=xi, sd=sqrt(tau2))

   cfunc <- function(xi, tau2, lower, upper)
      integrate(toint, lower=lower, upper=upper, xi=xi, tau2=tau2)$value

   if (targs$tau2 == 0) {
      zi <- transf.ztor(xi)
   } else {
      zi <- mapply(xi, FUN=cfunc, tau2=targs$tau2, lower=targs$lower, upper=targs$upper)
   }

   return(c(zi))

}

transf.exp.int <- function(xi, targs=NULL, ...) {

   if (is.null(targs$tau2))
      targs$tau2 <- 0
   if (is.null(targs$lower))
      targs$lower <- xi-5*sqrt(targs$tau2)
   if (is.null(targs$upper))
      targs$upper <- xi+5*sqrt(targs$tau2)

   toint <- function(zval, xi, tau2)
      exp(zval) * dnorm(zval, mean=xi, sd=sqrt(tau2))

   cfunc <- function(xi, tau2, lower, upper)
      integrate(toint, lower=lower, upper=upper, xi=xi, tau2=tau2)$value

   if (targs$tau2 == 0) {
      zi <- exp(xi)
   } else {
      zi <- mapply(xi, FUN=cfunc, tau2=targs$tau2, lower=targs$lower, upper=targs$upper)
   }

   return(c(zi))

}

transf.ilogit.int <- function(xi, targs=NULL, ...) {

   if (is.null(targs$tau2))
      targs$tau2 <- 0
   if (is.null(targs$lower))
      targs$lower <- xi-5*sqrt(targs$tau2)
   if (is.null(targs$upper))
      targs$upper <- xi+5*sqrt(targs$tau2)

   toint <- function(zval, xi, tau2)
      plogis(zval) * dnorm(zval, mean=xi, sd=sqrt(tau2))

   cfunc <- function(xi, tau2, lower, upper)
      integrate(toint, lower=lower, upper=upper, xi=xi, tau2=tau2)$value

   if (targs$tau2 == 0) {
      zi <- transf.ilogit(xi)
   } else {
      zi <- mapply(xi, FUN=cfunc, tau2=targs$tau2, lower=targs$lower, upper=targs$upper)
   }

   return(c(zi))

}

transf.dtou1 <- function(xi, ...) {
   u2i <- pnorm(abs(xi)/2)
   return((2*u2i - 1) / u2i)
}

transf.dtou2 <- function(xi, ...)
   pnorm(xi/2)

transf.dtou3 <- function(xi, ...)
   pnorm(xi)

transf.dtocles <- function(xi, ...)
   pnorm(xi/sqrt(2))

transf.dtorpb <- function(xi, n1i, n2i, ...) {
   if (missing(n1i) || missing(n2i)) {
      hi <- 4
   } else {
      if (length(n1i) != length(n2i))
         stop("Length of 'n1i' does not match length of 'n2i'.", call.=FALSE)
      if (length(n1i) != length(xi))
         stop("Length of 'n1i' and 'n2i' does not match length of 'xi'.", call.=FALSE)
      mi <- n1i + n2i - 2
      hi <- mi / n1i + mi / n2i
   }
   return(xi / sqrt(xi^2 + hi))
}

transf.dtobesd <- function(xi, ...) {
   rpbi <- xi / sqrt(xi^2 + 4)
   return(0.50 + rpbi/2)
}

transf.dtomd <- function(xi, targs=NULL, ...) {
   if (is.null(targs) || (is.list(targs) && is.null(targs$sd)))
      stop("Need to specify a standard deviation value via the 'targs' argument.", call.=FALSE)
   if (is.list(targs)) {
      sd <- targs$sd
   } else {
      sd <- targs
   }
   if (length(sd) != 1L)
      stop("Specify a single standard deviation value via the 'targs' argument.", call.=FALSE)
   return(xi * sd)
}

transf.logortord <- function(xi, pc, ...) {
   if (length(pc) == 1L)
      pc <- rep(pc, length(xi))
   if (length(xi) != length(pc))
      stop("Length of 'xi' does not match length of 'pc'.", call.=FALSE)
   if (any(pc < 0) || any(pc > 1))
      stop("The control group risk 'pc' must be between 0 and 1.", call.=FALSE)
   return(exp(xi)*pc / (1 - pc + pc * exp(xi)) - pc)
}

transf.logortorr <- function(xi, pc, ...) {
   if (length(pc) == 1L)
      pc <- rep(pc, length(xi))
   if (length(xi) != length(pc))
      stop("Length of 'xi' does not match length of 'pc'.", call.=FALSE)
   if (any(pc < 0) || any(pc > 1))
      stop("The control group risk 'pc' must be between 0 and 1.", call.=FALSE)
   return(exp(xi) / (pc * (exp(xi) - 1) + 1))
}
