############################################################################

.compvif <- function(j, btt, vcov, xintercept, intercept, spec=NULL, colnames=NULL, obj=NULL, coef="beta", sim=FALSE) {

   x <- obj

   btt <- btt[[j]] # note: this might actually be att when computing (G)VIFs for the scale coefficients in location-scale model

   if (is.null(x)) {

      ### remove intercept (if there is one and intercept=FALSE) from vcov and adjust btt accordingly

      if (xintercept && !intercept) {
         vcov <- vcov[-1,-1,drop=FALSE]
         btt  <- btt - 1
         btt  <- btt[btt > 0]
      }

      rb   <- suppressWarnings(cov2cor(vcov))
      gvif <- det(rb[btt,btt,drop=FALSE]) * det(rb[-btt,-btt,drop=FALSE]) / det(rb)

   } else {

      ### if 'x' is not NULL, then reestimate the model for the computation of the (G)VIF

      if (xintercept && !intercept)
         btt <- btt[btt > 1]

      if (coef == "beta") {
         Xbtt <- x$X[,btt,drop=FALSE]
         Zbtt <- x$Z
         if (xintercept && !intercept && !identical(btt,1L))
            Xbtt <- cbind(1, Xbtt)
         outlist <- "vb=vb"
      } else {
         Xbtt <- x$X
         Zbtt <- x$Z[,btt,drop=FALSE]
         if (xintercept && !intercept && !identical(btt,1L))
            Zbtt <- cbind(1, Zbtt)
         outlist <- "va=va"
      }

      if (inherits(x, "rma.uni")) {
         if (inherits(x, "rma.ls")) {
            args <- list(yi=x$yi, vi=x$vi, weights=x$weights, mods=Xbtt, intercept=FALSE, scale=Zbtt, link=x$link, method=x$method, weighted=x$weighted,
                         test=x$test, level=x$level, alpha=ifelse(x$alpha.fix, x$alpha, NA), optbeta=x$optbeta, beta=ifelse(x$beta.fix, x$beta, NA), control=x$control, skiphes=FALSE, outlist=outlist)
         } else {
            args <- list(yi=x$yi, vi=x$vi, weights=x$weights, mods=Xbtt, intercept=FALSE, method=x$method, weighted=x$weighted,
                         test=x$test, level=x$level, tau2=ifelse(x$tau2.fix, x$tau2, NA), control=x$control, skipr2=TRUE, outlist=outlist)
         }
         tmp <- try(suppressWarnings(.do.call(rma.uni, args)), silent=TRUE)
      }

      if (inherits(x, "rma.mv")) {
         args <- list(yi=x$yi, V=x$V, W=x$W, mods=Xbtt, random=x$random, struct=x$struct, intercept=FALSE, data=x$mf.r, method=x$method, test=x$test, dfs=x$dfs, level=x$level, R=x$R, Rscale=x$Rscale,
                      sigma2=ifelse(x$vc.fix$sigma2, x$sigma2, NA), tau2=ifelse(x$vc.fix$tau2, x$tau2, NA), rho=ifelse(x$vc.fix$rho, x$rho, NA), gamma2=ifelse(x$vc.fix$gamma2, x$gamma2, NA), phi=ifelse(x$vc.fix$phi, x$phi, NA),
                      sparse=x$sparse, dist=x$dist, control=x$control, outlist=outlist)
         tmp <- try(suppressWarnings(.do.call(rma.mv, args)), silent=TRUE)
      }

      if (inherits(tmp, "try-error")) {

         if (sim) {
            return(NA_real_)
         } else {
            gvif <- NA_real_
         }

      } else {

         if (xintercept && !intercept) {
            gvif <- det(vcov(x, type=coef)[btt,btt,drop=FALSE]) / det(vcov(tmp, type=coef)[-1,-1,drop=FALSE])
         } else {
            gvif <- det(vcov(x, type=coef)[btt,btt,drop=FALSE]) / det(vcov(tmp, type=coef))
         }

      }

   }

   if (sim) {

      return(gvif)

   } else {

      m <- length(btt)
      gsif <- gvif^(1/(2*m))

      ### readjust btt if this was done earlier

      if (is.null(x) && xintercept && !intercept)
         btt <- btt + 1

      if (length(btt) == 1L) {
         coefname <- colnames[btt]
      } else {
         coefname <- ""
      }

      return(data.frame(spec=.format.btt(spec[[j]]), coefs=.format.btt(btt), coefname=coefname, m=m, vif=gvif, sif=gsif))

   }

}

############################################################################

.compvifsim <- function(l, obj, coef, btt=NULL, att=NULL, reestimate=FALSE, intercept=FALSE, parallel=FALSE, seed=NULL) {

   if (parallel == "snow")
      library(metafor)

   if (!is.null(seed))
      set.seed(seed+l)

   x <- obj

   if (coef == "beta") {

      if (reestimate) {
         outlist <- "nodata"
      } else {
         outlist <- "coef.na=coef.na, vb=vb"
      }

      if (inherits(x, "rma.uni")) {
         if (inherits(x, "rma.ls")) {
            args <- list(yi=x$yi, vi=x$vi, weights=x$weights, mods=apply(x$X, 2, sample), intercept=FALSE, scale=x$Z, link=x$link, method=x$method, weighted=x$weighted,
                         test=x$test, level=x$level, alpha=ifelse(x$alpha.fix, x$alpha, NA), optbeta=x$optbeta, beta=ifelse(x$beta.fix, x$beta, NA), control=x$control, skiphes=FALSE, outlist=outlist)
         } else {
            args <- list(yi=x$yi, vi=x$vi, weights=x$weights, mods=apply(x$X, 2, sample), intercept=FALSE, method=x$method, weighted=x$weighted,
                         test=x$test, level=x$level, tau2=ifelse(x$tau2.fix, x$tau2, NA), control=x$control, skipr2=TRUE, outlist=outlist)
         }
         tmp <- try(suppressWarnings(.do.call(rma.uni, args)), silent=TRUE)
      }

      if (inherits(x, "rma.mv")) {
         args <- list(yi=x$yi, V=x$V, W=x$W, mods=apply(x$X, 2, sample), random=x$random, struct=x$struct, intercept=FALSE, data=x$mf.r, method=x$method, test=x$test, dfs=x$dfs, level=x$level, R=x$R, Rscale=x$Rscale,
                      sigma2=ifelse(x$vc.fix$sigma2, x$sigma2, NA), tau2=ifelse(x$vc.fix$tau2, x$tau2, NA), rho=ifelse(x$vc.fix$rho, x$rho, NA), gamma2=ifelse(x$vc.fix$gamma2, x$gamma2, NA), phi=ifelse(x$vc.fix$phi, x$phi, NA),
                      sparse=x$sparse, dist=x$dist, control=x$control, outlist=outlist)
         tmp <- try(suppressWarnings(.do.call(rma.mv, args)), silent=TRUE)
      }

      if (inherits(tmp, "try-error") || any(tmp$coef.na))
         return(rep(NA_real_, length(btt)))

      vcov <- vcov(tmp, type="beta")

      obj <- if (reestimate) tmp else NULL

      vifs <- sapply(seq_along(btt), .compvif, btt=btt, vcov=vcov, xintercept=x$intercept, intercept=intercept, obj=obj, sim=TRUE)

   } else {

      if (reestimate) {
         outlist <- "nodata"
      } else {
         outlist <- "coef.na.Z=coef.na.Z, va=va"
      }

      args <- list(yi=x$yi, vi=x$vi, weights=x$weights, mods=x$X, intercept=FALSE, scale=apply(x$Z, 2, sample), link=x$link, method=x$method, weighted=x$weighted,
                   test=x$test, level=x$level, alpha=ifelse(x$alpha.fix, x$alpha, NA), optbeta=x$optbeta, beta=ifelse(x$beta.fix, x$beta, NA), control=x$control, skiphes=FALSE, outlist=outlist)
      tmp <- try(suppressWarnings(.do.call(rma.uni, args)), silent=TRUE)

      if (inherits(tmp, "try-error") || any(tmp$Z.coef.na))
         return(rep(NA_real_, length(att)))

      vcov <- vcov(tmp, type="alpha")

      obj <- if (reestimate) tmp else NULL

      vifs <- sapply(seq_along(att), .compvif, btt=att, vcov=vcov, xintercept=x$Z.intercept, intercept=intercept, obj=obj, sim=TRUE)

   }

   return(vifs)

}

############################################################################
