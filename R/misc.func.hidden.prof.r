############################################################################

### for profile(), confint(), and gosh()

.profile.rma.uni <- function(val, obj, parallel=FALSE, profile=FALSE, confint=FALSE, subset=FALSE, vcs, objective, sel, FE=FALSE, verbose=FALSE) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (parallel == "snow")
      library(metafor)

   if (profile) {

      ### for profiling, fit model with tau2 fixed to the val-th value of 'vcs'

      res <- try(suppressWarnings(rma.uni(obj$yi, obj$vi, weights=obj$weights, mods=obj$X, intercept=FALSE, method=obj$method, weighted=obj$weighted, test=obj$test, level=obj$level, control=obj$control, tau2=vcs[val], skipr2=TRUE)), silent=TRUE)

      if (inherits(res, "try-error")) {
         sav <- list(ll = NA, beta = matrix(NA, nrow=nrow(obj$beta), ncol=1), ci.lb = rep(NA, length(obj$ci.lb)), ci.ub = rep(NA, length(obj$ci.ub)))
      } else {
         sav <- list(ll = logLik(res), beta = res$beta, ci.lb = res$ci.lb, ci.ub = res$ci.ub)
      }

   }

   if (confint) {

      ### for CI construction, fit model with tau2 fixed to 'val'

      res <- try(suppressWarnings(rma.uni(obj$yi, obj$vi, weights=obj$weights, mods=obj$X, intercept=FALSE, method=obj$method, weighted=obj$weighted, test=obj$test, level=obj$level, control=obj$control, tau2=val, skipr2=TRUE)), silent=TRUE)

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

      ### for subsetting, fit model to subset as specified in row 'val' of 'sel'

      if (FE) {

         if (parallel == "snow" || parallel == "multicore") {
            yi <- obj$yi[sel[val,]]
            vi <- obj$vi[sel[val,]]
         } else {
            yi <- obj$yi[sel]
            vi <- obj$vi[sel]
         }
         k <- length(yi)
         wi <- 1/vi
         est <- sum(wi*yi)/sum(wi)
         if (k > 1) {
            Q <- sum(wi * (yi - est)^2)
            I2 <- max(0, 100 * (Q - (k-1)) / Q)
            H2 <- Q / (k-1)
         } else {
            Q <- 0
            I2 <- 0
            H2 <- 1
         }
         tau2 <- 0
         if (parallel == "snow" || parallel == "multicore") {
            sav <- list(beta = est, het = c(k = k, QE = Q, I2 = I2, H2 = H2, tau2 = tau2))
         } else {
            sav <- list(beta = est, k = k, QE = Q, I2 = I2, H2 = H2, tau2 = tau2)
         }

      } else {

         res <- try(suppressWarnings(rma.uni(obj$yi, obj$vi, weights=obj$weights, mods=obj$X, intercept=FALSE, method=obj$method, weighted=obj$weighted, test=obj$test, level=obj$level, control=obj$control, tau2=ifelse(obj$tau2.fix, obj$tau2, NA), subset=sel[val,], skipr2=TRUE)), silent=TRUE)

         if (inherits(res, "try-error") || any(res$coef.na)) {
            sav <- list(beta = matrix(NA, nrow=nrow(obj$beta), ncol=1), het = rep(NA, 5))
         } else {
            sav <- list(beta = res$beta, het = c(res$k, res$QE, res$I2, res$H2, res$tau2))
         }

      }

   }

   return(sav)

}

.profile.rma.mv <- function(val, obj, comp, sigma2.pos, tau2.pos, rho.pos, gamma2.pos, phi.pos, parallel=FALSE, profile=FALSE, confint=FALSE, subset=FALSE, vcs, objective, verbose=FALSE) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (parallel == "snow")
      library(metafor)

   ### set any fixed components to their values
   sigma2.arg <- ifelse(obj$vc.fix$sigma2, obj$sigma2, NA)
   tau2.arg   <- ifelse(obj$vc.fix$tau2, obj$tau2, NA)
   rho.arg    <- ifelse(obj$vc.fix$rho, obj$rho, NA)
   gamma2.arg <- ifelse(obj$vc.fix$gamma2, obj$gamma2, NA)
   phi.arg    <- ifelse(obj$vc.fix$phi, obj$phi, NA)

   if (profile) {

      ### for profiling, fit model with component fixed to the val-th value of 'vcs'

      if (comp == "sigma2")
         sigma2.arg[sigma2.pos] <- vcs[val]

      if (comp == "tau2")
         tau2.arg[tau2.pos] <- vcs[val]

      if (comp == "rho")
         rho.arg[rho.pos] <- vcs[val]

      if (comp == "gamma2")
         gamma2.arg[gamma2.pos] <- vcs[val]

      if (comp == "phi")
         phi.arg[phi.pos] <- vcs[val]

      res <- try(suppressWarnings(rma.mv(obj$yi, V=obj$V, W=obj$W, mods=obj$X, random=obj$random, struct=obj$struct, intercept=FALSE, data=obj$mf.r, method=obj$method, test=obj$test, level=obj$level, R=obj$R, Rscale=obj$Rscale, sigma2=sigma2.arg, tau2=tau2.arg, rho=rho.arg, gamma2=gamma2.arg, phi=phi.arg, sparse=obj$sparse, dist=obj$dist, control=obj$control)), silent=TRUE)

      if (inherits(res, "try-error")) {
         sav <- list(ll = NA, beta = matrix(NA, nrow=nrow(obj$beta), ncol=1), ci.lb = rep(NA, length(obj$ci.lb)), ci.ub = rep(NA, length(obj$ci.ub)))
      } else {
         sav <- list(ll = logLik(res), beta = res$beta, ci.lb = res$ci.lb, ci.ub = res$ci.ub)
      }

   }

   if (confint) {

      ### for CI construction, fit model with component fixed to 'val'

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

      res <- try(suppressWarnings(rma.mv(obj$yi, V=obj$V, W=obj$W, mods=obj$X, random=obj$random, struct=obj$struct, intercept=FALSE, data=obj$mf.r, method=obj$method, test=obj$test, level=obj$level, R=obj$R, Rscale=obj$Rscale, sigma2=sigma2.arg, tau2=tau2.arg, rho=rho.arg, gamma2=gamma2.arg, phi=phi.arg, sparse=obj$sparse, dist=obj$dist, control=obj$control)), silent=TRUE)

      if (inherits(res, "try-error")) {

         if (verbose)
            cat(mstyle$verbose(paste("vc =", formatC(val, digits=obj$digits[["var"]], width=obj$digits[["var"]]+4, format="f"), "  LRT - objective = NA", "\n")))

         stop()

      } else {

         sav <- -2*(logLik(res) - logLik(obj)) - objective

         if (verbose)
            cat(mstyle$verbose(paste("vc =", formatC(val, digits=obj$digits[["var"]], width=obj$digits[["var"]]+4, format="f"), "  LRT - objective =", formatC(sav, digits=obj$digits[["fit"]], width=obj$digits[["fit"]]+4, format="f"), "\n")))

      }

   }

   return(sav)

}

.profile.rma.mh <- function(val, obj, parallel=FALSE, subset=FALSE, sel) {

   if (parallel == "snow")
      library(metafor)

   if (subset) {

      ### for subsetting, fit model to subset as specified in row 'val' of 'sel'

      if (is.element(obj$measure, c("RR","OR","RD"))) {
         res <- try(suppressWarnings(rma.mh(ai=obj$ai, bi=obj$bi, ci=obj$ci, di=obj$di, measure=obj$measure, add=obj$add, to=obj$to, drop00=obj$drop00, correct=obj$correct, subset=sel[val,])), silent=TRUE)
      } else {
         res <- try(suppressWarnings(rma.mh(x1i=obj$x1i, x2i=obj$x2i, t1i=obj$t1i, t2i=obj$t2i, measure=obj$measure, add=obj$add, to=obj$to, drop00=obj$drop00, correct=obj$correct, subset=sel[val,])), silent=TRUE)
      }

      if (inherits(res, "try-error")) {
         sav <- list(beta = NA, het = rep(NA, 5))
      } else {
         sav <- list(beta = res$beta, het = c(res$k, res$QE, res$I2, res$H2, res$tau2))
      }

   }

   return(sav)

}

.profile.rma.peto <- function(val, obj, parallel=FALSE, subset=FALSE, sel) {

   if (parallel == "snow")
      library(metafor)

   if (subset) {

      ### for subsetting, fit model to subset as specified in row 'val' of 'sel'

      res <- try(suppressWarnings(rma.peto(ai=obj$ai, bi=obj$bi, ci=obj$ci, di=obj$di, add=obj$add, to=obj$to, drop00=obj$drop00, subset=sel[val,])), silent=TRUE)

      if (inherits(res, "try-error")) {
         sav <- list(beta = NA, het = rep(NA, 5))
      } else {
         sav <- list(beta = res$beta, het = c(res$k, res$QE, res$I2, res$H2, res$tau2))
      }

   }

   return(sav)

}

############################################################################
