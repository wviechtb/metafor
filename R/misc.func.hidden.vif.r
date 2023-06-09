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
                         test=x$test, level=x$level, alpha=ifelse(x$alpha.fix, x$alpha, NA), optbeta=x$optbeta, beta=ifelse(x$beta.fix, x$beta, NA),
                         control=x$control, skiphes=FALSE, outlist=outlist)
         } else {
            args <- list(yi=x$yi, vi=x$vi, weights=x$weights, mods=Xbtt, intercept=FALSE, method=x$method, weighted=x$weighted,
                         test=x$test, level=x$level, tau2=ifelse(x$tau2.fix, x$tau2, NA), control=x$control, skipr2=TRUE, outlist=outlist)
         }
         tmp <- try(suppressWarnings(.do.call(rma.uni, args)), silent=TRUE)
      }

      if (inherits(x, "rma.mv")) {
         args <- list(yi=x$yi, V=x$V, W=x$W, mods=Xbtt, random=x$random, struct=x$struct, intercept=FALSE, data=x$mf.r, method=x$method,
                      test=x$test, dfs=x$dfs, level=x$level, R=x$R, Rscale=x$Rscale,
                      sigma2=ifelse(x$vc.fix$sigma2, x$sigma2, NA), tau2=ifelse(x$vc.fix$tau2, x$tau2, NA), rho=ifelse(x$vc.fix$rho, x$rho, NA),
                      gamma2=ifelse(x$vc.fix$gamma2, x$gamma2, NA), phi=ifelse(x$vc.fix$phi, x$phi, NA),
                      sparse=x$sparse, dist=x$dist, vccon=obj$vccon, control=x$control, outlist=outlist)
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

.compvifsim <- function(l, obj, coef, btt=NULL, att=NULL, reestimate=FALSE, intercept=FALSE, parallel=FALSE, seed=NULL, joinb=NULL, joina=NULL) {

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

      if (is.null(joinb)) {
         if (is.null(x$data) || is.null(x$formula.mods)) {
            Xperm <- apply(x$X, 2, sample)
         } else {
            #data <- x$data
            data <- get_all_vars(x$formula.mods, data=x$data) # only get variables that are actually needed
            if (!is.null(x$subset))
               data <- data[x$subset,,drop=FALSE]
            data <- data[x$not.na,,drop=FALSE]
            Xperm <- model.matrix(x$formula.mods, data=as.data.frame(lapply(data, sample)))
            #Xperm <- Xperm[,!x$coef.na,drop=FALSE]
         }
      } else {
         Xperm <- .permXvif(joinb, x$X)
      }

      if (inherits(x, "rma.uni")) {
         if (inherits(x, "rma.ls")) {
            args <- list(yi=x$yi, vi=x$vi, weights=x$weights, mods=Xperm, intercept=FALSE, scale=x$Z, link=x$link, method=x$method, weighted=x$weighted,
                         test=x$test, level=x$level, alpha=ifelse(x$alpha.fix, x$alpha, NA), optbeta=x$optbeta, beta=ifelse(x$beta.fix, x$beta, NA),
                         control=x$control, skiphes=FALSE, outlist=outlist)
         } else {
            args <- list(yi=x$yi, vi=x$vi, weights=x$weights, mods=Xperm, intercept=FALSE, method=x$method, weighted=x$weighted,
                         test=x$test, level=x$level, tau2=ifelse(x$tau2.fix, x$tau2, NA), control=x$control, skipr2=TRUE, outlist=outlist)
         }
         tmp <- try(suppressWarnings(.do.call(rma.uni, args)), silent=TRUE)
         #tmp <- try(.do.call(rma.uni, args))
      }

      if (inherits(x, "rma.mv")) {
         args <- list(yi=x$yi, V=x$V, W=x$W, mods=Xperm, random=x$random, struct=x$struct, intercept=FALSE, data=x$mf.r, method=x$method,
                      test=x$test, dfs=x$dfs, level=x$level, R=x$R, Rscale=x$Rscale,
                      sigma2=ifelse(x$vc.fix$sigma2, x$sigma2, NA), tau2=ifelse(x$vc.fix$tau2, x$tau2, NA), rho=ifelse(x$vc.fix$rho, x$rho, NA),
                      gamma2=ifelse(x$vc.fix$gamma2, x$gamma2, NA), phi=ifelse(x$vc.fix$phi, x$phi, NA),
                      sparse=x$sparse, dist=x$dist, vccon=obj$vccon, control=x$control, outlist=outlist)
         tmp <- try(suppressWarnings(.do.call(rma.mv, args)), silent=TRUE)
      }

      if (inherits(tmp, "try-error"))
         return(rep(NA_real_, length(btt)))

      if (any(tmp$coef.na))
         return(sapply(btt, function(x) if (any(which(tmp$coef.na) %in% x)) Inf else NA_real_))

      vcov <- vcov(tmp, type="beta")

      obj <- if (reestimate) tmp else NULL

      vifs <- sapply(seq_along(btt), .compvif, btt=btt, vcov=vcov, xintercept=x$intercept, intercept=intercept, obj=obj, sim=TRUE)

   } else {

      if (reestimate) {
         outlist <- "nodata"
      } else {
         outlist <- "coef.na.Z=coef.na.Z, va=va"
      }

      if (is.null(joina)) {
         if (is.null(x$data) || is.null(x$formula.scale)) {
            Zperm <- apply(x$Z, 2, sample)
         } else {
            #data <- x$data
            data <- get_all_vars(x$formula.scale, data=x$data) # only get variables that are actually needed
            if (!is.null(x$subset))
               data <- data[x$subset,,drop=FALSE]
            data <- data[x$not.na,,drop=FALSE]
            Zperm <- model.matrix(x$formula.scale, data=as.data.frame(lapply(data, sample)))
            #Zperm <- Zperm[,!x$coef.na.Z,drop=FALSE]
         }
      } else {
         Zperm <- .permXvif(joina, x$Z)
      }

      args <- list(yi=x$yi, vi=x$vi, weights=x$weights, mods=x$X, intercept=FALSE, scale=Zperm, link=x$link, method=x$method, weighted=x$weighted,
                   test=x$test, level=x$level, alpha=ifelse(x$alpha.fix, x$alpha, NA), optbeta=x$optbeta, beta=ifelse(x$beta.fix, x$beta, NA),
                   control=x$control, skiphes=FALSE, outlist=outlist)
      tmp <- try(suppressWarnings(.do.call(rma.uni, args)), silent=TRUE)
      #tmp <- try(.do.call(rma.uni, args))

      if (inherits(tmp, "try-error"))
         return(rep(NA_real_, length(att)))

      if (any(tmp$coef.na.Z))
         return(sapply(att, function(x) if (any(which(tmp$coef.na.Z) %in% x)) Inf else NA_real_))

      vcov <- vcov(tmp, type="alpha")

      obj <- if (reestimate) tmp else NULL

      vifs <- sapply(seq_along(att), .compvif, btt=att, vcov=vcov, xintercept=x$Z.intercept, intercept=intercept, obj=obj, sim=TRUE)

   }

   return(vifs)

}

.permXvif <- function(b, X) {
   ub <- unique(b)
   n <- nrow(X)
   for (j in 1:length(ub)) {
      pos <- which(ub[j] == b)
      X[,pos] <- X[sample(n),pos]
   }
   return(X)
}

############################################################################
