### for profile(), confint(), and gosh()

.profile.rma.uni <- function(val, obj, parallel=FALSE, profile=FALSE, confint=FALSE, subset=FALSE, objective, model=0L, verbose=FALSE, outlist=NULL) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (parallel == "snow")
      library(metafor)

   if (profile || confint) {

      ### for profile and confint, fit model with tau2 fixed to 'val'

      res <- try(suppressWarnings(rma.uni(obj$yi, obj$vi, weights=obj$weights, mods=obj$X, intercept=FALSE,
                                          method=obj$method, weighted=obj$weighted, test=obj$test, level=obj$level,
                                          control=obj$control, tau2=val, skipr2=TRUE, outlist="minimal")), silent=TRUE)

   }

   if (profile) {

      if (inherits(res, "try-error")) {
         sav <- list(ll = NA, beta = matrix(NA, nrow=nrow(obj$beta), ncol=1), ci.lb = rep(NA, length(obj$ci.lb)), ci.ub = rep(NA, length(obj$ci.ub)))
      } else {
         sav <- list(ll = logLik(res), beta = res$beta, ci.lb = res$ci.lb, ci.ub = res$ci.ub)
      }

   }

   if (confint) {

      if (inherits(res, "try-error")) {

         if (verbose)
            cat(mstyle$verbose(paste("tau2 =", formatC(val, digits=obj$digits[["var"]], width=obj$digits[["var"]]+4, format="f"), "  LRT - objective = NA", "\n")))

         stop()

      } else {

         sav <- -2*(logLik(res) - logLik(obj)) - objective

         if (verbose)
            cat(mstyle$verbose(paste("tau2 =", formatC(val, digits=obj$digits[["var"]], width=obj$digits[["var"]]+4, format="f"), "  LRT - objective =", formatC(sav, digits=obj$digits[["test"]], width=obj$digits[["test"]]+4, format="f"), "\n")))

      }

   }

   if (subset) {

      ### for subset, fit model to subset as specified by 'val'

      if (model >= 1L) {

         # special cases for gosh() for FE and RE+DL models

         yi <- obj$yi[val]
         vi <- obj$vi[val]
         k <- length(yi)
         wi <- 1/vi
         sumwi <- sum(wi)
         est <- sum(wi*yi)/sumwi
         Q <- 0
         I2 <- 0
         H2 <- 1
         tau2 <- 0
         if (k > 1) {
            Q <- sum(wi * (yi - est)^2)
            I2 <- max(0, 100 * (Q - (k-1)) / Q)
            H2 <- Q / (k-1)
            if (model == 2L) {
               tau2 <- max(0, (Q - (k-1)) / (sumwi - sum(wi^2)/sumwi))
               wi <- 1 / (vi + tau2)
               est <- sum(wi*yi)/sum(wi)
            }
         }
         sav <- list(beta = est, k = k, QE = Q, I2 = I2, H2 = H2, tau2 = tau2)

      } else {

         sav <- try(suppressWarnings(rma.uni(obj$yi, obj$vi, weights=obj$weights, mods=obj$X, intercept=FALSE,
                                             method=obj$method, weighted=obj$weighted, test=obj$test, level=obj$level,
                                             control=obj$control, tau2=ifelse(obj$tau2.fix, obj$tau2, NA),
                                             subset=val, skipr2=TRUE, outlist=outlist)), silent=TRUE)

      }

   }

   return(sav)

}

.profile.rma.mv <- function(val, obj, comp, sigma2.pos, tau2.pos, rho.pos, gamma2.pos, phi.pos, parallel=FALSE, profile=FALSE, confint=FALSE, subset=FALSE, objective, verbose=FALSE) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (parallel == "snow")
      library(metafor)

   if (profile || confint) {

      ### for profile and confint, fit model with component fixed to 'val'

      ### set any fixed components to their values
      sigma2.arg <- ifelse(obj$vc.fix$sigma2, obj$sigma2, NA)
      tau2.arg   <- ifelse(obj$vc.fix$tau2, obj$tau2, NA)
      rho.arg    <- ifelse(obj$vc.fix$rho, obj$rho, NA)
      gamma2.arg <- ifelse(obj$vc.fix$gamma2, obj$gamma2, NA)
      phi.arg    <- ifelse(obj$vc.fix$phi, obj$phi, NA)

      if (comp == "sigma2")
         sigma2.arg[sigma2.pos] <- val

      if (comp == "tau2")
         tau2.arg[tau2.pos] <- val

      if (comp == "rho")
         rho.arg[rho.pos] <- val

      if (comp == "gamma2")
         gamma2.arg[gamma2.pos] <- val

      if (comp == "phi")
         phi.arg[phi.pos] <- val

      res <- try(suppressWarnings(rma.mv(obj$yi, V=obj$V, W=obj$W, mods=obj$X, random=obj$random, struct=obj$struct, intercept=FALSE,
                                         data=obj$mf.r, method=obj$method, test=obj$test, level=obj$level, R=obj$R, Rscale=obj$Rscale,
                                         sigma2=sigma2.arg, tau2=tau2.arg, rho=rho.arg, gamma2=gamma2.arg, phi=phi.arg, sparse=obj$sparse,
                                         dist=obj$dist, control=obj$control, outlist="minimal")), silent=TRUE)

   }

   if (profile) {

      if (inherits(res, "try-error")) {
         sav <- list(ll = NA, beta = matrix(NA, nrow=nrow(obj$beta), ncol=1), ci.lb = rep(NA, length(obj$ci.lb)), ci.ub = rep(NA, length(obj$ci.ub)))
      } else {
         sav <- list(ll = logLik(res), beta = res$beta, ci.lb = res$ci.lb, ci.ub = res$ci.ub)
      }

   }

   if (confint) {

      if (inherits(res, "try-error")) {

         if (verbose)
            cat(mstyle$verbose(paste("val =", formatC(val, digits=obj$digits[["var"]], width=obj$digits[["var"]]+4, format="f"), "  LRT - objective = NA", "\n")))

         stop()

      } else {

         sav <- -2*(logLik(res) - logLik(obj)) - objective

         if (verbose)
            cat(mstyle$verbose(paste("val =", formatC(val, digits=obj$digits[["var"]], width=obj$digits[["var"]]+4, format="f"), "  LRT - objective =", formatC(sav, digits=obj$digits[["fit"]], width=obj$digits[["fit"]]+4, format="f"), "\n")))

      }

   }

   return(sav)

}

.profile.rma.mh <- function(val, obj, parallel=FALSE, subset=FALSE, outlist=NULL) {

   if (parallel == "snow")
      library(metafor)

   if (subset) {

      ### for subset, fit model to subset as specified by 'val'

      if (is.element(obj$measure, c("RR","OR","RD"))) {
         sav <- try(suppressWarnings(rma.mh(ai=obj$ai, bi=obj$bi, ci=obj$ci, di=obj$di, measure=obj$measure,
                                            add=obj$add, to=obj$to, drop00=obj$drop00, correct=obj$correct,
                                            level=obj$level, subset=val, outlist=outlist)), silent=TRUE)
      } else {
         sav <- try(suppressWarnings(rma.mh(x1i=obj$x1i, x2i=obj$x2i, t1i=obj$t1i, t2i=obj$t2i, measure=obj$measure,
                                            add=obj$add, to=obj$to, drop00=obj$drop00, correct=obj$correct,
                                            level=obj$level, subset=val, outlist=outlist)), silent=TRUE)
      }

   }

   return(sav)

}

.profile.rma.peto <- function(val, obj, parallel=FALSE, subset=FALSE, outlist=NULL) {

   if (parallel == "snow")
      library(metafor)

   if (subset) {

      ### for subset, fit model to subset as specified by 'val'

      sav <- try(suppressWarnings(rma.peto(ai=obj$ai, bi=obj$bi, ci=obj$ci, di=obj$di,
                                           add=obj$add, to=obj$to, drop00=obj$drop00,
                                           level=obj$level, subset=val, outlist=outlist)), silent=TRUE)

   }

   return(sav)

}

.profile.rma.uni.selmodel <- function(val, obj, comp, delta.pos, parallel=FALSE, profile=FALSE, confint=FALSE, subset=FALSE, objective, verbose=FALSE) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (parallel == "snow")
      library(metafor)

   if (profile || confint) {

      ### for profile and confint, fit model with component fixed to 'val'

      ### set any fixed components to their values
      tau2.arg  <- ifelse(obj$method == "FE" || obj$tau2.fix, obj$tau2, NA)
      delta.arg <- ifelse(obj$delta.fix, obj$delta, NA)

      if (comp == "tau2")
         tau2.arg <- val

      if (comp == "delta")
         delta.arg[delta.pos] <- val

      ### reset steps to NA if !stepsspec (some types set steps=0 if steps was not specified)
      if (!obj$stepsspec)
         obj$steps <- NA

      res <- try(suppressWarnings(selmodel(obj, obj$type, alternative=obj$alternative, prec=obj$prec, scaleprec=obj$scaleprec, tau2=tau2.arg, delta=delta.arg, steps=obj$steps, verbose=FALSE, control=obj$control, skiphes=FALSE, skiphet=TRUE, defmap=obj$defmap, mapfun=obj$mapfun, mapinvfun=obj$mapinvfun)), silent=TRUE)

   }

   if (profile) {

      if (inherits(res, "try-error")) {
         sav <- list(ll = NA, beta = matrix(NA, nrow=nrow(obj$beta), ncol=1), ci.lb = rep(NA, length(obj$ci.lb)), ci.ub = rep(NA, length(obj$ci.ub)))
      } else {
         sav <- list(ll = logLik(res), beta = res$beta, ci.lb = res$ci.lb, ci.ub = res$ci.ub)
      }

   }

   if (confint) {

      if (inherits(res, "try-error")) {

         if (verbose)
            cat(mstyle$verbose(paste("val =", formatC(val, digits=obj$digits[["var"]], width=obj$digits[["var"]]+4, format="f"), "  LRT - objective = NA", "\n")))

         stop()

      } else {

         sav <- -2*(logLik(res) - logLik(obj)) - objective

         if (verbose)
            cat(mstyle$verbose(paste("val =", formatC(val, digits=obj$digits[["var"]], width=obj$digits[["var"]]+4, format="f"), "  LRT - objective =", formatC(sav, digits=obj$digits[["fit"]], width=obj$digits[["fit"]]+4, format="f"), "\n")))

      }

   }

   return(sav)

}
