############################################################################

.chktargsint <- function(targs) {

   if (length(targs) > 3L)
      stop("Length of 'targs' argument must be <= 3.", call.=FALSE)

   if (.is.vector(targs)) {
      if (is.null(names(targs))) {
         names(targs) <- c("tau2", "lower", "upper")[seq_along(targs)]
         targs <- as.list(targs)
      } else {
         targs <- list(tau2=unname(targs[startsWith(names(targs), "t")]), lower=unname(targs[startsWith(names(targs), "l")]), upper=unname(targs[startsWith(names(targs), "u")]))
         targs <- targs[lengths(targs) > 0L]
      }
   }

   if (any(lengths(targs) > 1L))
      stop("Elements of 'targs' arguments must be scalars.", call.=FALSE)

   if (is.null(targs$tau2))
      targs$tau2 <- 0

   return(targs)

}

############################################################################

transf.rtoz <- function(xi) { # resulting value between -Inf (for -1) and +Inf (for +1)
   xi[xi >  1] <-  1
   xi[xi < -1] <- -1
   atanh(xi) # same as 1/2 * log((1+xi)/(1-xi))
}

transf.ztor <- function(xi)
   tanh(xi) # same as (exp(2*xi)-1)/(exp(2*xi)+1)

transf.ztor.int <- function(xi, targs=NULL) {

   if (is.na(xi))
      return(NA_real_)

   targs <- .chktargsint(targs)

   if (is.null(targs$lower))
      targs$lower <- xi-10*sqrt(targs$tau2)
   if (is.null(targs$upper))
      targs$upper <- xi+10*sqrt(targs$tau2)

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

transf.r2toz <- function(xi) {
   xi[xi > 1] <- 1
   xi[xi < 0] <- 0
   atanh(sqrt(xi))
}

transf.ztor2 <- function(xi)
   tanh(xi)^2

############################################################################

transf.exp.int <- function(xi, targs=NULL) {

   if (is.na(xi))
      return(NA_real_)

   targs <- .chktargsint(targs)

   if (is.null(targs$lower))
      targs$lower <- xi-10*sqrt(targs$tau2)
   if (is.null(targs$upper))
      targs$upper <- xi+10*sqrt(targs$tau2)

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

############################################################################

transf.logit <- function(xi) # resulting value between -Inf (for 0) and +Inf (for +1)
   qlogis(xi)

transf.ilogit <- function(xi)
   plogis(xi)

transf.ilogit.int <- function(xi, targs=NULL) {

   if (is.na(xi))
      return(NA_real_)

   targs <- .chktargsint(targs)

   if (is.null(targs$lower))
      targs$lower <- xi-10*sqrt(targs$tau2)
   if (is.null(targs$upper))
      targs$upper <- xi+10*sqrt(targs$tau2)

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

############################################################################

transf.arcsin <- function(xi) # resulting value between 0 (for 0) and asin(1) = pi/2 (for 1)
   asin(sqrt(xi))

transf.iarcsin <- function(xi) {
   zi <- sin(xi)^2
   zi[xi < 0] <- 0       # if xi value is below 0 (e.g., CI bound), return 0
   zi[xi > asin(1)] <- 1 # if xi value is above maximum possible value, return 1
   return(c(zi))
}

# transf.iarcsin.int <- function(xi, targs=NULL) {
#
#   if (is.na(xi))
#      return(NA_real_)
#
#   targs <- .chktargsint(targs)
#
#   if (is.null(targs$lower))
#      targs$lower <- 0
#   if (is.null(targs$upper))
#      targs$upper <- asin(1)
#
#   toint <- function(zval, xi, tau2)
#      transf.iarcsin(zval) * dnorm(zval, mean=xi, sd=sqrt(tau2))
#
#   cfunc <- function(xi, tau2, lower, upper)
#      integrate(toint, lower=lower, upper=upper, xi=xi, tau2=tau2)$value
#
#   if (targs$tau2 == 0) {
#      zi <- transf.iarcsin(xi)
#   } else {
#      zi <- mapply(xi, FUN=cfunc, tau2=targs$tau2, lower=targs$lower, upper=targs$upper)
#   }
#
#   return(c(zi))
#
# }

############################################################################

transf.pft <- function(xi, ni) {                   # Freeman-Tukey transformation for proportions
   xi <- xi*ni
   zi <- 1/2*(asin(sqrt(xi/(ni+1))) + asin(sqrt((xi+1)/(ni+1))))
   return(c(zi))
}

transf.ipft <- function(xi, ni) {                  # inverse of Freeman-Tukey transformation for individual proportions
   zi <- suppressWarnings(1/2 * (1 - sign(cos(2*xi)) * sqrt(1 - (sin(2*xi)+(sin(2*xi)-1/sin(2*xi))/ni)^2)))
   zi <- ifelse(is.nan(zi), NA_real_, zi)
   zi[xi > transf.pft(1,ni)] <- 1                  # if xi is above upper limit, return 1
   zi[xi < transf.pft(0,ni)] <- 0                  # if xi is below lower limit, return 0
   return(c(zi))
}

transf.ipft.hm <- function(xi, targs) {            # inverse of Freeman-Tukey transformation for a collection of proportions
   if (is.null(targs) || (is.list(targs) && is.null(targs$ni)))
      stop("Must specify the sample sizes via the 'targs' argument.", call.=FALSE)
   if (is.list(targs)) {
      ni <- targs$ni
   } else {
      ni <- ni
   }
   nhm <- 1/(mean(1/ni, na.rm=TRUE))               # calculate harmonic mean of the ni's
   zi <- suppressWarnings(1/2 * (1 - sign(cos(2*xi)) * sqrt(1 - (sin(2*xi)+(sin(2*xi)-1/sin(2*xi))/nhm)^2)))
   zi <- ifelse(is.nan(zi), NA_real_, zi)          # it may not be possible to calculate zi
   zi[xi > transf.pft(1,nhm)] <- 1                 # if xi is above upper limit, return 1
   zi[xi < transf.pft(0,nhm)] <- 0                 # if xi is below lower limit, return 0
   return(c(zi))
}

############################################################################

transf.isqrt <- function(xi) {
   zi <- xi*xi
   zi[xi < 0] <- 0                                 # if xi value is below 0 (e.g., CI bound), return 0
   return(c(zi))
}

############################################################################

transf.irft <- function(xi, ti) {                  # Freeman-Tukey transformation for incidence rates
   zi <- 1/2*(sqrt(xi) + sqrt(xi + 1/ti))          # xi is the incidence rate (not the number of events!)
   return(c(zi))
}

transf.iirft <- function(xi, ti) {                 # inverse of Freeman-Tukey transformation for incidence rates (see Freeman-Tukey_incidence.r in code directory)
   #zi <- (1/ti - 2*xi^2 + ti*xi^4)/(4*xi^2*ti)    # old version where transf.irft was not multiplied by 1/2
   zi <- (1/ti - 8*xi^2 + 16*ti*xi^4)/(16*xi^2*ti) # xi is the incidence rate (not the number of events!)
   zi <- ifelse(is.nan(zi), NA_real_, zi)
   zi[xi < transf.irft(0,ti)] <- 0                 # if xi is below lower limit, return 0
   zi[zi <= .Machine$double.eps] <- 0              # avoid finite precision errors in back-transformed values (transf.iirft(transf.irft(0, 1:200), 1:200))
   return(c(zi))
}

############################################################################

transf.ahw <- function(xi) {              # resulting value between 0 (for alpha=0) and 1 (for alpha=1)
   #zi <- (1-xi)^(1/3)
   zi <- 1 - (1-xi)^(1/3)
   return(c(zi))
}

transf.iahw <- function(xi) {
   #zi <- 1-xi^3
   zi <- 1 - (1-xi)^3
   zi <- ifelse(is.nan(zi), NA_real_, zi)
   zi[xi > 1] <- 1                        # if xi is above upper limit, return 1
   zi[xi < 0] <- 0                        # if xi is below lower limit, return 0
   return(c(zi))
}

transf.abt <- function(xi) {              # Bonett (2002) transformation of alphas (without bias correction)
#transf.abt <- function(xi, ni) {         # resulting value between 0 (for alpha=0) to Inf (for alpha=1)
   #zi <- log(1-xi) - log(ni/(ni-1))
   #zi <- log(1-xi)
   zi <- -log(1-xi)
   return(c(zi))
}

transf.iabt <- function(xi) {             # inverse of Bonett (2002) transformation
#transf.iabt <- function(xi, ni) {
   #zi <- 1 - exp(xi) * ni / (ni-1)
   #zi <- 1 - exp(xi)
   zi <- 1 - exp(-xi)
   zi <- ifelse(is.nan(zi), NA_real_, zi)
   zi[xi < 0] <- 0                        # if xi is below lower limit, return 0
   return(c(zi))
}

############################################################################

transf.dtou1 <- function(xi) {
   u2i <- pnorm(abs(xi)/2)
   return((2*u2i - 1) / u2i)
}

transf.dtou2 <- function(xi)
   pnorm(xi/2)

transf.dtou3 <- function(xi)
   pnorm(xi)

transf.dtoovl <- function(xi)
   2*pnorm(-abs(xi)/2)

transf.dtocles <- function(xi) # note: this does not assume homoscedasticity
   pnorm(xi/sqrt(2))

transf.dtocliffd <- function(xi) # note: this does not assume homoscedasticity
   2 * pnorm(xi/sqrt(2)) - 1

transf.dtobesd <- function(xi) {
   rpbi <- xi / sqrt(xi^2 + 4)
   return(0.50 + rpbi/2)
}

transf.dtomd <- function(xi, targs=NULL) {
   if (is.null(targs) || (is.list(targs) && is.null(targs$sd)))
      stop("Must specify a standard deviation value via the 'targs' argument.", call.=FALSE)
   if (is.list(targs)) {
      sd <- targs$sd
   } else {
      sd <- targs
   }
   if (length(sd) != 1L)
      stop("Specify a single standard deviation value via the 'targs' argument.", call.=FALSE)
   return(xi * sd)
}

transf.dtorpb <- function(xi, n1i, n2i) {
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

transf.dtorbis <- function(xi, n1i, n2i) {
   if (missing(n1i) || missing(n2i)) {
      hi <- 4
      n1i <- 1
      n2i <- 1
   } else {
      if (length(n1i) != length(n2i))
         stop("Length of 'n1i' does not match length of 'n2i'.", call.=FALSE)
      if (length(n1i) != length(xi))
         stop("Length of 'n1i' and 'n2i' does not match length of 'xi'.", call.=FALSE)
      mi <- n1i + n2i - 2
      hi <- mi / n1i + mi / n2i
   }
   rpbi <- xi / sqrt(xi^2 + hi)
   pi <- n1i / (n1i + n2i)
   return(sqrt(pi*(1-pi)) / dnorm(qnorm(pi)) * rpbi)
}

transf.rpbtorbis <- function(xi, pi) {
   if (missing(pi)) {
      pi <- 0.5
   } else {
      pi <- .expand1(pi, length(xi))
      if (length(xi) != length(pi))
         stop("Length of 'xi' does not match length of 'pi'.", call.=FALSE)
   }
   if (any(pi < 0 | pi > 1, na.rm=TRUE))
      stop("One or more 'pi' values are < 0 or > 1.", call.=FALSE)
   return(sqrt(pi*(1-pi)) / dnorm(qnorm(pi)) * xi)
}

transf.rtorpb <- function(xi, pi) {
   if (missing(pi)) {
      pi <- 0.5
   } else {
      pi <- .expand1(pi, length(xi))
      if (length(xi) != length(pi))
         stop("Length of 'xi' does not match length of 'pi'.", call.=FALSE)
   }
   if (any(pi < 0 | pi > 1, na.rm=TRUE))
      stop("One or more 'pi' values are < 0 or > 1.", call.=FALSE)
   return(xi * dnorm(qnorm(pi)) / sqrt(pi*(1-pi)))
}

transf.rtod <- function(xi, n1i, n2i) {
   if (missing(n1i) || missing(n2i)) {
      hi <- 4
      n1i <- 1
      n2i <- 1
   } else {
      if (length(n1i) != length(n2i))
         stop("Length of 'n1i' does not match length of 'n2i'.", call.=FALSE)
      if (length(n1i) != length(xi))
         stop("Length of 'n1i' and 'n2i' does not match length of 'xi'.", call.=FALSE)
      mi <- n1i + n2i - 2
      hi <- mi / n1i + mi / n2i
   }
   if (any(c(n1i < 0, n2i < 0), na.rm=TRUE))
      stop("One or more values specified via the 'n1i' or 'n2i' arguments are negative.")
   pi <- n1i / (n1i + n2i)
   rpbi <- xi * dnorm(qnorm(pi)) / sqrt(pi*(1-pi))
   return(sqrt(hi) * rpbi / sqrt(1 - rpbi^2))
}

transf.rpbtod <- function(xi, n1i, n2i) {
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
   return(sqrt(hi) * xi / sqrt(1 - xi^2))
}

transf.lnortord <- function(xi, pc) {
   pc <- .expand1(pc, length(xi))
   if (length(xi) != length(pc))
      stop("Length of 'xi' does not match length of 'pc'.", call.=FALSE)
   if (any(pc < 0) || any(pc > 1))
      stop("The control group risk 'pc' must be between 0 and 1.", call.=FALSE)
   return(exp(xi)*pc / (1 - pc + pc * exp(xi)) - pc)
}

transf.lnortorr <- function(xi, pc) {
   pc <- .expand1(pc, length(xi))
   if (length(xi) != length(pc))
      stop("Length of 'xi' does not match length of 'pc'.", call.=FALSE)
   if (any(pc < 0) || any(pc > 1))
      stop("The control group risk 'pc' must be between 0 and 1.", call.=FALSE)
   return(exp(xi) / (pc * (exp(xi) - 1) + 1))
}

############################################################################

transf.lnortod.norm <- function(xi)
   xi / 1.65

transf.lnortod.logis <- function(xi)
   sqrt(3) / base::pi * xi

transf.dtolnor.norm <- function(xi)
   xi * 1.65

transf.dtolnor.logis <- function(xi)
   xi / sqrt(3) * base::pi

transf.lnortortet.pearson <- function(xi)
   cos(base::pi / (1 + sqrt(exp(xi))))

transf.lnortortet.digby <- function(xi)
   (exp(xi)^(3/4) - 1) / (exp(xi)^(3/4) + 1)

############################################################################
