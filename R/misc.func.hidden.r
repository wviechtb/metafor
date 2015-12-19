############################################################################

### function to set default 'btt' value(s) or check specified 'btt' values

.set.btt <- function(btt, p, int.incl) {

   if (missing(btt) || is.null(btt)) {

      if (p > 1) {                        ### if the model matrix has more than one column
         if (int.incl) {
            btt <- seq.int(from=2, to=p)     ### and the model has an intercept term, test all coefficients except the intercept
         } else {
            btt <- seq_len(p)                ### and the model does not have an intercept term, test all coefficients
         }
      } else {
         btt <- 1                         ### if the model matrix has a single column, test that single coefficient
      }

   } else {

      btt <- btt[(btt >= 1) & (btt <= p)] ### weed out values below 1 or above p
      btt <- unique(round(btt))           ### round values and then take unique values

      if (length(btt) == 0L)              ### make sure that at least one valid value is left
         stop("Non-existent coefficients specified via 'btt'.")
   }

   return(btt)

}

############################################################################

### pairwise sorting of the elements of two vectors

.psort <- function(x,y) {

   ### t(apply(xy, 1, sort)) would be okay, but problematic if there are NAs;
   ### either they are removed completely (na.last=NA) or they are always put
   ### first/last (na.last=FALSE/TRUE); but we just want to leave the NAs in
   ### their position!

   if (is.null(x) || length(x) == 0) ### need to catch this
      return(NULL)

   if (missing(y)) {
      if (is.matrix(x)) {
         xy <- x
      } else {
         xy <- rbind(x) ### in case x is just a vector
      }
   } else {
      xy <- cbind(x,y)
   }

   n <- nrow(xy)

   for (i in 1:n) {
      if (anyNA(xy[i,]))
         next
      xy[i,] <- sort(xy[i,])
   }

   return(xy)

}

############################################################################

### c(m) calculation function for bias correction of SMDs (mi = n1i + n2i - 2) or SMCC/SMCRs (mi = ni - 1)

.cmicalc <- function(mi) {

   ### this can overflow if mi is 'large' (on my machine, if mi >= 344)
   #cmi <- gamma(mi/2)/(sqrt(mi/2)*gamma((mi-1)/2))
   ### catch those cases and apply the approximate formula (which is accurate then)
   #is.na <- is.na(cmi)
   #cmi[is.na] <- 1 - 3/(4*mi[is.na] - 1)

   ### this avoids the problem with overflow altogether
   cmi <- ifelse(mi <= 1, NA, exp(lgamma(mi/2) - log(sqrt(mi/2)) - lgamma((mi-1)/2)))
   return(cmi)

}

############################################################################

### function to obtain the trace of a matrix

.tr <- function(X)
   return(sum(diag(X)))

############################################################################

### function to test whether a vector is all equal to 1s (to find intercept(s) in a model matrix)

.is.int.func <- function(x, eps=1e-08)
   all(abs(x - 1) < eps)

############################################################################

### function to format p-values
### if showeq=FALSE, c(.001, .00001) becomes c("0.0010", "<.0001")
### if showeq=TRUE,  c(.001, .00001) becomes c("=0.0010", "<.0001")

.pval <- function(p, digits=4, showeq=FALSE, sep="") {

   cutoff  <- paste(c(".", rep(0,digits-1),1), collapse="")
   ncutoff <- as.numeric(cutoff)

   ifelse(is.na(p), paste0(ifelse(showeq, "=", ""), sep, NA),
                    ifelse(p >= ncutoff, paste0(ifelse(showeq, "=", ""), sep, formatC(p, digits=digits, format="f")),
                                         paste0("<", sep, cutoff)))

}

############################################################################

### function to calculate:
### solve(t(X) %*% W %*% X) = .invcalc(X=X, W=W, k=k)
### solve(t(X) %*% X)       = .invcalc(X=X, W=diag(k), k=k)
### without taking the actual inverse

.invcalc <- function(X, W, k) {

   sWX     <- sqrt(W) %*% X
   res.qrs <- qr.solve(sWX, diag(k))
   return(tcrossprod(res.qrs))

}

############################################################################

### function for confint.rma.uni() with Q-profile method and for the PM estimator

.QE.func <- function(tau2val, Y, vi, X, k, objective, verbose=FALSE, digits=4) {

   if (any(tau2val + vi < 0))
      stop("Some marginal variances are negative.")

   W     <- diag(1/(vi + tau2val), nrow=k, ncol=k)
   stXWX <- .invcalc(X=X, W=W, k=k)
   P     <- W - W %*% X %*% stXWX %*% crossprod(X,W)
   RSS   <- crossprod(Y,P) %*% Y

   if (verbose)
      cat("tau2 =", formatC(tau2val, digits=digits, width=digits+4, format="f"), " RSS - objective =", c(RSS - objective), "\n")

   return(RSS - objective)

}

############################################################################

### function for confint.rma.uni() with method="GENQ"

.GENQ.func <- function(tau2val, P, vi, Q, alpha, k, p, getlower, verbose=FALSE, digits=4) {

   S <- diag(sqrt(vi + tau2val), nrow=k, ncol=k)
   lambda <- Re(eigen(S %*% P %*% S, symmetric=TRUE, only.values=TRUE)$values)
   if (getlower) {
      res <- CompQuadForm::farebrother(Q, lambda[1:(k-p)])$res - alpha
   } else {
      res <- (1 - CompQuadForm::farebrother(Q, lambda[1:(k-p)])$res) - alpha
   }

   if (verbose)
      cat("tau2 =", formatC(tau2val, digits=digits, width=digits+4, format="f"), " objective =", res, "\n")

   return(res)

}

############################################################################

### function to construct vcov matrix for struct="UN" given vector of variances and correlations

.con.vcov.UN <- function(vars, cors) {
   dims <- length(vars)
   G <- matrix(1, nrow=dims, ncol=dims)
   G[upper.tri(G)] <- cors
   G[lower.tri(G)] <- t(G)[lower.tri(G)]
   H <- diag(sqrt(vars), nrow=dims, ncol=dims)
   return(H %*% G %*% H)
}

### function to construct G matrix for struct="UN" given vector of 'choled' variances and covariances

.con.vcov.UN.chol <- function(vars, covs) {
   dims <- length(vars)
   G <- matrix(0, nrow=dims, ncol=dims)
   G[upper.tri(G)] <- covs
   diag(G) <- vars
   return(crossprod(G))
}

############################################################################

### -1 times the log likelihood (regular or restricted) for rma.mv models

.ll.rma.mv <- function(par, reml, Y, M, X.fit, k, pX, # note: need to call model matrix X.fit due to hessian()
                       D.S, Z.G1, Z.G2, Z.H1, Z.H2,
                       sigma2.val, tau2.val, rho.val, gamma2.val, phi.val,
                       sigma2s, tau2s, rhos, gamma2s, phis,
                       withS, withG, withH,
                       struct, g.levels.r, h.levels.r,
                       tol, sparse, cholesky, posdefify, vctransf, verbose, very.verbose, digits, REMLf) {

   ### only NA values in sigma2.val, tau2.val, rho.val, gamma2.val, phi.val should be estimated; otherwise, replace with fixed values

   if (withS) {

      if (vctransf) {
         ### sigma2 is optimized in log-space, so exponentiate
         sigma2 <- ifelse(is.na(sigma2.val), exp(par[1:sigma2s]), sigma2.val)
      } else {
         ### for Hessian computation, leave as is
         sigma2 <- ifelse(is.na(sigma2.val), par[1:sigma2s], sigma2.val)
         sigma2[sigma2 < 0] <- 0
      }

      for (j in seq_len(sigma2s)) {
         M <- M + sigma2[j] * D.S[[j]]
      }

   }

   if (withG) {

      ### if cholesky=TRUE, back-transformation/substitution is done below; otherwise, back-transform and replace fixed values
      if (cholesky[1]) {
         tau2 <- par[(sigma2s+1):(sigma2s+tau2s)]
         rho  <- par[(sigma2s+tau2s+1):(sigma2s+tau2s+rhos)]
      } else {
         if (vctransf) {
            ### tau2 is optimized in log-space, so exponentiate
            tau2 <- ifelse(is.na(tau2.val), exp(par[(sigma2s+1):(sigma2s+tau2s)]), tau2.val)
            ### rho is optimized in r-to-z space, so use transf.ztor
            rho  <- ifelse(is.na(rho.val), transf.ztor(par[(sigma2s+tau2s+1):(sigma2s+tau2s+rhos)]), rho.val)
         } else {
            ### for Hessian computation, leave as is
            tau2 <- ifelse(is.na(tau2.val), par[(sigma2s+1):(sigma2s+tau2s)], tau2.val)
            rho  <- ifelse(is.na(rho.val), par[(sigma2s+tau2s+1):(sigma2s+tau2s+rhos)], rho.val)
            tau2[tau2 < 0] <- 0
            rho[rho >  1] <-  1
            rho[rho < -1] <- -1
         }
      }

      ncol.Z.G1 <- ncol(Z.G1)

      if (struct[1] == "CS") {
         G <- matrix(rho*tau2, nrow=ncol.Z.G1, ncol=ncol.Z.G1)
         diag(G) <- tau2
      }

      if (struct[1] == "HCS") {
         G <- matrix(rho, nrow=ncol.Z.G1, ncol=ncol.Z.G1)
         diag(G) <- 1
         G <- diag(sqrt(tau2), nrow=ncol.Z.G1, ncol=ncol.Z.G1) %*% G %*% diag(sqrt(tau2), nrow=ncol.Z.G1, ncol=ncol.Z.G1)
         diag(G) <- tau2
      }

      if (struct[1] == "UN") {
         if (cholesky[1]) {
            G <- .con.vcov.UN.chol(tau2, rho)
            tau2 <- diag(G)                 ### need this, so correct values are shown when verbose=TRUE
            rho <- cov2cor(G)[upper.tri(G)] ### need this, so correct values are shown when verbose=TRUE
            tau2[!is.na(tau2.val)] <- tau2.val[!is.na(tau2.val)] ### replace any fixed values
            rho[!is.na(rho.val)] <- rho.val[!is.na(rho.val)] ### replace any fixed values
         }
         G <- .con.vcov.UN(tau2, rho)
         if (posdefify) {
            G <- as.matrix(nearPD(G)$mat)   ### nearPD() in Matrix package
            tau2 <- diag(G)                 ### need this, so correct values are shown when verbose=TRUE
            rho <- cov2cor(G)[upper.tri(G)] ### need this, so correct values are shown when verbose=TRUE
         }
      }

      if (struct[1] == "ID" || struct[1] == "DIAG") {
         G <- diag(tau2, nrow=ncol.Z.G1, ncol=ncol.Z.G1)
      }

      if (struct[1] == "UNHO") {
         G <- matrix(1, nrow=ncol.Z.G1, ncol=ncol.Z.G1)
         G[upper.tri(G)] <- rho
         G[lower.tri(G)] <- t(G)[lower.tri(G)]
         G <- diag(sqrt(rep(tau2, ncol.Z.G1)), nrow=ncol.Z.G1, ncol=ncol.Z.G1) %*% G %*% diag(sqrt(rep(tau2, ncol.Z.G1)), nrow=ncol.Z.G1, ncol=ncol.Z.G1)
         if (posdefify) {
            G <- as.matrix(nearPD(G, keepDiag=TRUE)$mat) ### nearPD() in Matrix package
            tau2 <- G[1,1]                               ### need this, so correct values are shown when verbose=TRUE
            rho <- cov2cor(G)[upper.tri(G)]              ### need this, so correct values are shown when verbose=TRUE
         }
      }

      if (struct[1] == "AR") {
         if (ncol.Z.G1 > 1) {
            G <- toeplitz(ARMAacf(ar=rho, lag.max=ncol.Z.G1-1))
         } else {
            G <- diag(1)
         }
         G <- diag(sqrt(rep(tau2, ncol.Z.G1)), nrow=ncol.Z.G1, ncol=ncol.Z.G1) %*% G %*% diag(sqrt(rep(tau2, ncol.Z.G1)), nrow=ncol.Z.G1, ncol=ncol.Z.G1)
         diag(G) <- tau2
      }

      if (struct[1] == "HAR") {
         if (ncol.Z.G1 > 1) {
            G <- toeplitz(ARMAacf(ar=rho, lag.max=ncol.Z.G1-1))
         } else {
            G <- diag(1)
         }
         G <- diag(sqrt(tau2), nrow=ncol.Z.G1, ncol=ncol.Z.G1) %*% G %*% diag(sqrt(tau2), nrow=ncol.Z.G1, ncol=ncol.Z.G1)
         diag(G) <- tau2
      }

      ### set tau2 and corresponding rho value(s) to 0 for any levels that were removed
      if (any(g.levels.r)) {
         G[g.levels.r,] <- 0
         G[,g.levels.r] <- 0
      }

      if (sparse)
         G <- Matrix(G, sparse=TRUE)

      M <- M + (Z.G1 %*% G %*% t(Z.G1)) * tcrossprod(Z.G2)

   }

   if (withH) {

      ### if cholesky=TRUE, back-transformation/substitution is done below; otherwise, back-transform and replace fixed values
      if (cholesky[2]) {
         gamma2 <- par[(sigma2s+tau2s+rhos+1):(sigma2s+tau2s+rhos+gamma2s)]
         phi    <- par[(sigma2s+tau2s+rhos+gamma2s+1):(sigma2s+tau2s+rhos+gamma2s+phis)]
      } else {
         if (vctransf) {
            ### gamma2 is optimized in log-space, so exponentiate
            gamma2 <- ifelse(is.na(gamma2.val), exp(par[(sigma2s+tau2s+rhos+1):(sigma2s+tau2s+rhos+gamma2s)]), gamma2.val)
            ### phi is optimized in r-to-z space, so use transf.ztor
            phi    <- ifelse(is.na(phi.val), transf.ztor(par[(sigma2s+tau2s+rhos+gamma2s+1):(sigma2s+tau2s+rhos+gamma2s+phis)]), phi.val)
         } else {
            ### for Hessian computation, leave as is
            gamma2 <- ifelse(is.na(gamma2.val), par[(sigma2s+tau2s+rhos+1):(sigma2s+tau2s+rhos+gamma2s)], gamma2.val)
            phi    <- ifelse(is.na(phi.val), par[(sigma2s+tau2s+rhos+gamma2s+1):(sigma2s+tau2s+rhos+gamma2s+phis)], phi.val)
            gamma2[gamma2 < 0] <- 0
            phi[phi >  1] <-  1
            phi[phi < -1] <- -1
         }
      }

      ncol.Z.H1 <- ncol(Z.H1)

      if (struct[2] == "CS") {
         H <- matrix(phi*gamma2, nrow=ncol.Z.H1, ncol=ncol.Z.H1)
         diag(H) <- gamma2
      }

      if (struct[2] == "HCS") {
         H <- matrix(phi, nrow=ncol.Z.H1, ncol=ncol.Z.H1)
         diag(H) <- 1
         H <- diag(sqrt(gamma2), nrow=ncol.Z.H1, ncol=ncol.Z.H1) %*% H %*% diag(sqrt(gamma2), nrow=ncol.Z.H1, ncol=ncol.Z.H1)
         diag(H) <- gamma2
      }

      if (struct[2] == "UN") {
         if (cholesky[2]) {
            H <- .con.vcov.UN.chol(gamma2, phi)
            gamma2 <- diag(H)               ### need this, so correct values are shown when verbose=TRUE
            phi <- cov2cor(H)[upper.tri(H)] ### need this, so correct values are shown when verbose=TRUE
            gamma2[!is.na(gamma2.val)] <- gamma2.val[!is.na(gamma2.val)] ### replace any fixed values
            phi[!is.na(phi.val)] <- phi.val[!is.na(phi.val)] ### replace any fixed values
         }
         H <- .con.vcov.UN(gamma2, phi)
         if (posdefify) {
            H <- as.matrix(nearPD(H)$mat)   ### nearPD() in Matrix package
            gamma2 <- diag(H)               ### need this, so correct values are shown when verbose=TRUE
            phi <- cov2cor(H)[upper.tri(H)] ### need this, so correct values are shown when verbose=TRUE
         }
      }

      if (struct[2] == "ID" || struct[2] == "DIAG") {
         H <- diag(gamma2, nrow=ncol.Z.H1, ncol=ncol.Z.H1)
      }

      if (struct[2] == "UNHO") {
         H <- matrix(1, nrow=ncol.Z.H1, ncol=ncol.Z.H1)
         H[upper.tri(H)] <- phi
         H[lower.tri(H)] <- t(H)[lower.tri(H)]
         H <- diag(sqrt(rep(gamma2, ncol.Z.H1)), nrow=ncol.Z.H1, ncol=ncol.Z.H1) %*% H %*% diag(sqrt(rep(gamma2, ncol.Z.H1)), nrow=ncol.Z.H1, ncol=ncol.Z.H1)
         if (posdefify) {
            H <- as.matrix(nearPD(H, keepDiag=TRUE)$mat) ### nearPD() in Matrix package
            gamma2 <- H[1,1]                             ### need this, so correct values are shown when verbose=TRUE
            phi <- cov2cor(H)[upper.tri(H)]              ### need this, so correct values are shown when verbose=TRUE
         }
      }

      if (struct[2] == "AR") {
         if (ncol.Z.H1 > 1) {
            H <- toeplitz(ARMAacf(ar=phi, lag.max=ncol.Z.H1-1))
         } else {
            H <- diag(1)
         }
         H <- diag(sqrt(rep(gamma2, ncol.Z.H1)), nrow=ncol.Z.H1, ncol=ncol.Z.H1) %*% H %*% diag(sqrt(rep(gamma2, ncol.Z.H1)), nrow=ncol.Z.H1, ncol=ncol.Z.H1)
         diag(H) <- gamma2
      }

      if (struct[2] == "HAR") {
         if (ncol.Z.H1 > 1) {
            H <- toeplitz(ARMAacf(ar=phi, lag.max=ncol.Z.H1-1))
         } else {
            H <- diag(1)
         }
         H <- diag(sqrt(gamma2), nrow=ncol.Z.H1, ncol=ncol.Z.H1) %*% H %*% diag(sqrt(gamma2), nrow=ncol.Z.H1, ncol=ncol.Z.H1)
         diag(H) <- gamma2
      }

      ### set gamma2 and corresponding phi value(s) to 0 for any levels that were removed
      if (any(h.levels.r)) {
         H[h.levels.r,] <- 0
         H[,h.levels.r] <- 0
      }

      if (sparse)
         H <- Matrix(H, sparse=TRUE)

      M <- M + (Z.H1 %*% H %*% t(Z.H1)) * tcrossprod(Z.H2)

   }

   if (posdefify)
      M <- as.matrix(nearPD(M)$mat)

   if (verbose) {
      L <- try(chol(M), silent=!verbose)
   } else {
      L <- try(suppressWarnings(chol(M)), silent=!verbose)
   }

   if (inherits(L, "try-error")) {

      ### if M is not positive-definite, set the (restricted) log likelihood to -Inf
      ### this idea is based on: http://stats.stackexchange.com/q/11368/1934 (this is crude, but should
      ### move the parameter estimates away from values that create the non-positive-definite M matrix)

      llval <- -Inf

   } else {

      W <- chol2inv(L) ### need this for REML llval computation
      U <- chol(W)
      sX <- U %*% X.fit

      sY <- U %*% Y
      b  <- solve(crossprod(sX), crossprod(sX, sY))
      RSS.f <- sum(as.vector(sY - sX %*% b)^2)

      if (reml) {
         llval  <- -1/2 * (k-pX) * log(2*base::pi) + ifelse(REMLf, 1/2 * determinant(crossprod(X.fit), logarithm=TRUE)$modulus, 0) - 1/2 * determinant(M, logarithm=TRUE)$modulus - 1/2 * determinant(crossprod(X.fit,W) %*% X.fit, logarithm=TRUE)$modulus - 1/2 * RSS.f
      } else {
         llval  <- -1/2 * (k)    * log(2*base::pi)                                                                                 - 1/2 * determinant(M, logarithm=TRUE)$modulus                                                                           - 1/2 * RSS.f
      }

   }

   if ((vctransf && verbose) || (!vctransf && very.verbose)) {
      if (withS)
         cat("sigma2 =", ifelse(is.na(sigma2), NA, paste(formatC(sigma2, digits=digits, format="f", flag=" "), " ", sep="")), "  ", sep="")
      if (withG)
         cat("tau2 =",   ifelse(is.na(tau2),   NA, paste(formatC(tau2,   digits=digits, format="f", flag=" "), " ", sep="")), "  ", sep="")
      if (withG)
         cat("rho =",    ifelse(is.na(rho),    NA, paste(formatC(rho,    digits=digits, format="f", flag=" "), " ", sep="")), "  ", sep="")
      if (withH)
         cat("gamma2 =", ifelse(is.na(gamma2), NA, paste(formatC(gamma2, digits=digits, format="f", flag=" "), " ", sep="")), "  ", sep="")
      if (withH)
         cat("phi =",    ifelse(is.na(phi),    NA, paste(formatC(phi,    digits=digits, format="f", flag=" "), " ", sep="")), "  ", sep="")
      cat("  ll = ", ifelse(is.na(llval), NA, formatC(llval, digits=digits, format="f", flag=" ")), sep="", "\n")
   }

   return(-1 * llval)

}

############################################################################

### needed for multicore processing with profile()

.profile.rma.uni <- function(val, obj, parallel=FALSE, CI=FALSE, objective, verbose=FALSE) {

   if (parallel == "snow")
      library(metafor)

   res <- try(suppressWarnings(rma.uni(obj$yi, obj$vi, weights=obj$weights, mods=obj$X, method=obj$method, weighted=obj$weighted, intercept=FALSE, knha=obj$knha, level=obj$level, control=obj$control, tau2=val)), silent=TRUE)

   if (!CI) {

      if (inherits(res, "try-error")) {
         list(ll = NA, b = matrix(NA, nrow=nrow(obj$b), ncol=1), ci.lb = rep(NA, length(obj$ci.lb)), ci.ub = rep(NA, length(obj$ci.ub)))
      } else {
         list(ll = logLik(res), b = res$b, ci.lb = res$ci.lb, ci.ub = res$ci.ub)
      }

   } else {

      if (inherits(res, "try-error")) {

         if (verbose)
            cat("tau2 =", formatC(val, digits=obj$digits, width=obj$digits+4, format="f"), " LRT - objective = NA", "\n")

         stop()

      } else {

         difference <- -2*(logLik(res) - logLik(obj)) - objective

         if (verbose)
            cat("tau2 =", formatC(val, digits=obj$digits, width=obj$digits+4, format="f"), " LRT - objective =", difference, "\n")

         return(difference)

      }

   }

}

.profile.rma.mv <- function(val, obj, comp, sigma2.pos, tau2.pos, rho.pos, gamma2.pos, phi.pos, parallel=FALSE, CI=FALSE, objective, verbose=FALSE) {

   if (parallel == "snow")
      library(metafor)

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

   res <- try(suppressWarnings(rma.mv(obj$yi, obj$V, obj$W, mods=obj$X, random=obj$random, struct=obj$struct, intercept=FALSE, method=obj$method, tdist=obj$knha, level=obj$level, R=obj$R, Rscale=obj$Rscale, data=obj$mf.r, sigma2=sigma2.arg, tau2=tau2.arg, rho=rho.arg, gamma2=gamma2.arg, phi=phi.arg, control=obj$control)), silent=TRUE)

   if (!CI) {

      if (inherits(res, "try-error")) {
         list(ll = NA, b = matrix(NA, nrow=nrow(obj$b), ncol=1), ci.lb = rep(NA, length(obj$ci.lb)), ci.ub = rep(NA, length(obj$ci.ub)))
      } else {
         list(ll = logLik(res), b = res$b, ci.lb = res$ci.lb, ci.ub = res$ci.ub)
      }

   } else {

      if (inherits(res, "try-error")) {

         if (verbose)
            cat("vc =", formatC(val, digits=obj$digits, width=obj$digits+4, format="f"), " LRT - objective = NA", "\n")

         stop()

      } else {

         difference <- -2*(logLik(res) - logLik(obj)) - objective

         if (verbose)
            cat("vc =", formatC(val, digits=obj$digits, width=obj$digits+4, format="f"), " LRT - objective =", difference, "\n")

         return(difference)

      }

   }

}

############################################################################

### generate all possible permutations

.genperms <- function(k) {

   v <- seq_len(k)

   sub <- function(k, v) {
      if (k==1L) {
         matrix(v,1,k)
      } else {
         X  <-  NULL
         for(i in seq_len(k)) {
            X <- rbind(X, cbind(v[i], Recall(k-1, v[-i])))
         }
      X
      }
   }

   return(sub(k, v[seq_len(k)]))

}

### generate all possible unique permutations

.genuperms <- function(x) {

   z <- NULL

   sub <- function(x, y) {
      len.x <- length(x)
      if (len.x == 0L) {
         return(y)
      } else {
         prev.num <- 0
         for (i in seq_len(len.x)) {
            num <- x[i]
            if (num > prev.num) {
               prev.num <- num
               z <- rbind(z, Recall(x[-i], c(y,num)))
            }
         }
         return(z)
      }
   }

   return(sub(x, y=NULL))

}

############################################################################

### set axis label (for forest, funnel, and labbe functions)

.setlab <- function(measure, transf.char, atransf.char, gentype) {

   if (gentype == 1)
      lab <- "Observed Outcome"
   if (gentype == 2)
      lab <- "Overall Estimate" ### need this for forest.cumul.rma() function

   #########################################################################

   if (!is.null(measure)) {

      ######################################################################
      if (measure == "RR") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- "Log Relative Risk"
         } else {
            lab <- "Transformed Log Relative Risk"
            if (atransf.char == "exp" || atransf.char == "transf.exp.int")
               lab <- "Relative Risk (log scale)"
            if (transf.char == "exp" || transf.char == "transf.exp.int")
               lab <- "Relative Risk"
         }
      }
      if (is.element(measure, c("OR","PETO","D2OR","D2ORN","D2ORL"))) {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- "Log Odds Ratio"
         } else {
            lab <- "Transformed Log Odds Ratio"
            if (atransf.char == "exp" || atransf.char == "transf.exp.int")
               lab <- "Odds Ratio (log scale)"
            if (transf.char == "exp" || transf.char == "transf.exp.int")
               lab <- "Odds Ratio"
         }
      }
      if (measure == "RD") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- "Risk Difference"
         } else {
            lab <- "Transformed Risk Difference"
         }
      }
      if (measure == "AS") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- "Arcsine Transformed Risk Difference"
         } else {
            lab <- "Transformed Arcsine Transformed Risk Difference"
         }
      }
      if (measure == "PHI") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- "Phi Coefficient"
         } else {
            lab <- "Transformed Phi Coefficient"
         }
      }
      if (measure == "YUQ") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- "Yule's Q"
         } else {
            lab <- "Transformed Yule's Q"
         }
      }
      if (measure == "YUY") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- "Yule's Y"
         } else {
            lab <- "Transformed Yule's Y"
         }
      }
      ######################################################################
      if (measure == "IRR") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- "Log Incidence Rate Ratio"
         } else {
            lab <- "Transformed Log Incidence Relative Risk"
            if (atransf.char == "exp" || atransf.char == "transf.exp.int")
               lab <- "Incidence Rate Ratio (log scale)"
            if (transf.char == "exp" || transf.char == "transf.exp.int")
               lab <- "Incidence Rate Ratio"
         }
      }
      if (measure == "IRD") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- "Incidence Rate Difference"
         } else {
            lab <- "Transformed Incidence Rate Difference"
         }
      }
      if (measure == "IRSD") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- "Square-Root Transformed Incidence Rate Difference"
         } else {
            lab <- "Transformed Square-Root Transformed Incidence Rate Difference"
         }
      }
      ######################################################################
      if (measure == "MD") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- "Mean Difference"
         } else {
            lab <- "Transformed Mean Difference"
         }
      }
      if (is.element(measure, c("SMD","SMDH","PBIT","OR2D","OR2DN","OR2DL"))) {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- "Standardized Mean Difference"
         } else {
            lab <- "Transformed Standardized Mean Difference"
         }
      }
      if (measure == "ROM") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- "Log Ratio of Means"
         } else {
            lab <- "Transformed Log Ratio of Means"
            if (atransf.char == "exp" || atransf.char == "transf.exp.int")
               lab <- "Ratio of Means (log scale)"
            if (transf.char == "exp" || transf.char == "transf.exp.int")
               lab <- "Ratio of Means"
         }
      }
      if (measure == "RPB") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- "Point-Biserial Correlation"
         } else {
            lab <- "Transformed Point-Biserial Correlation"
         }
      }
      ######################################################################
      if (is.element(measure, c("COR","UCOR","RTET","RBIS"))) {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- "Correlation Coefficient"
         } else {
            lab <- "Transformed Correlation Coefficient"
         }
      }
      if (measure == "ZCOR") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- "Fisher's z Transformed Correlation Coefficient"
         } else {
            lab <- "Transformed Fisher's z Transformed Correlation Coefficient"
            if (atransf.char == "transf.ztor" || atransf.char == "transf.ztor.int")
               lab <- "Correlation Coefficient"
            if (transf.char == "transf.ztor" || transf.char == "transf.ztor.int")
               lab <- "Correlation Coefficient"
         }
      }
      ######################################################################
      if (measure == "PR") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- "Proportion"
         } else {
            lab <- "Transformed Proportion"
         }
      }
      if (measure == "PLN") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- "Log Proportion"
         } else {
            lab <- "Transformed Log Proportion"
            if (atransf.char == "exp" || atransf.char == "transf.exp.int")
               lab <- "Proportion (log scale)"
            if (transf.char == "exp" || transf.char == "transf.exp.int")
               lab <- "Proportion"
         }
      }
      if (measure == "PLO") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- "Log Odds"
         } else {
            lab <- "Transformed Log Odds"
            if (atransf.char == "transf.ilogit" || atransf.char == "transf.ilogit.int" || atransf.char == "plogis")
               lab <- "Proportion (logit scale)"
            if (transf.char == "transf.ilogit" || transf.char == "transf.ilogit.int" || transf.char == "plogis")
               lab <- "Proportion"
            if (atransf.char == "exp" || atransf.char == "transf.exp.int")
               lab <- "Odds (log scale)"
            if (transf.char == "exp" || transf.char == "transf.exp.int")
               lab <- "Odds"
         }
      }
      if (measure == "PAS") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- "Arcsine Transformed Proportion"
         } else {
            lab <- "Transformed Arcsine Transformed Proportion"
            if (atransf.char == "transf.iarcsin" || atransf.char == "transf.iarcsin.int")
               lab <- "Proportion (arcsine scale)"
            if (transf.char == "transf.iarcsin" || transf.char == "transf.iarcsin.int")
               lab <- "Proportion"
         }
      }
      if (measure == "PFT") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- "Double Arcsine Transformed Proportion"
         } else {
            lab <- "Transformed Double Arcsine Transformed Proportion"
            if (atransf.char == "transf.ift.hm")
               lab <- "Proportion"
            if (transf.char == "transf.ift.hm")
               lab <- "Proportion"
         }
      }
      ######################################################################
      if (measure == "IR") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- "Incidence Rate"
         } else {
            lab <- "Transformed Incidence Rate"
         }
      }
      if (measure == "IRLN") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- "Log Incidence Rate"
         } else {
            lab <- "Transformed Log Incidence Rate"
            if (atransf.char == "exp" || atransf.char == "transf.exp.int")
               lab <- "Incidence Rate (log scale)"
            if (transf.char == "exp" || transf.char == "transf.exp.int")
               lab <- "Incidence Rate"
         }
      }
      if (measure == "IRS") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- "Square-Root Transformed Incidence Rate"
         } else {
            lab <- "Transformed Square-Root Transformed Incidence Rate"
            if (atransf.char == "transf.isqrt" || atransf.char == "transf.isqrt.int")
               lab <- "Incidence Rate (square-root scale)"
            if (transf.char == "transf.isqrt" || transf.char == "transf.isqrt.int")
               lab <- "Incidence Rate"
         }
      }
      if (measure == "IRFT") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- "Freeman-Tukey Transformed Incidence Rate"
         } else {
            lab <- "Transformed Freeman-Tukey Transformed Incidence Rate"
         }
      }
      ######################################################################
      if (measure == "MN") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- "Mean"
         } else {
            lab <- "Transformed Mean"
         }
      }
      ######################################################################
      if (measure == "MC") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- "Mean Change"
         } else {
            lab <- "Transformed Mean Change"
         }
      }
      if (is.element(measure, c("SMCC","SMCR","SMCRH"))) {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- "Standardized Mean Change"
         } else {
            lab <- "Transformed Standardized Mean Change"
         }
      }
      if (measure == "ROMC") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- "Log Ratio of Means"
         } else {
            lab <- "Transformed Log Ratio of Means"
            if (atransf.char == "exp" || atransf.char == "transf.exp.int")
               lab <- "Ratio of Means (log scale)"
            if (transf.char == "exp" || transf.char == "transf.exp.int")
               lab <- "Ratio of Means"
         }
      }
      ######################################################################
      if (measure == "ARAW") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- "Coefficient alpha"
         } else {
            lab <- "Transformed Coefficient alpha"
         }
      }
      if (measure == "AHW") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- "Transformed Coefficient alpha"
         } else {
            lab <- "Transformed Coefficient alpha"
            if (atransf.char == "transf.iahw")
               lab <- "Coefficient alpha"
            if (transf.char == "transf.iahw")
               lab <- "Coefficient alpha"
         }
      }
      if (measure == "ABT") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- "Transformed Coefficient alpha"
         } else {
            lab <- "Transformed Coefficient alpha"
            if (atransf.char == "transf.iabt")
               lab <- "Coefficient alpha"
            if (transf.char == "transf.iabt")
               lab <- "Coefficient alpha"
         }
      }
      ######################################################################

   }

   return(lab)

}

############################################################################

### density of non-central hypergeometric distribution (based on Liao and Rosen, 2001) from MCMCpack
### Liao, J. G. & Rosen, O. (2001). Fast and stable algorithms for computing and sampling from the
### noncentral hypergeometric distribution. The American Statistician, 55, 366-369.

.dnoncenhypergeom <- function (x=NA, n1, n2, m1, psi) { ### x=ai, n1=ai+bi, n2=ci+di, m1=ai+ci, psi=ORi

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
      stop("'n1' or 'n2' negative in dnoncenhypergeom().\n")

   if (m1 < 0 | m1 > (n1 + n2))
      stop("'m1' out of range in dnoncenhypergeom().\n")

   if (psi <= 0)
      stop("'psi' [odds ratio] negative in dnoncenhypergeom().\n")

   if (!is.na(x) & (x < ll | x > uu))
      stop("'x' out of bounds in dnoncenhypergeom().\n")

   if (!is.na(x) & length(x) > 1)
      stop("'x' neither missing or scalar in dnoncenhypergeom().\n")

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
         stop(paste0("Could not compute density of non-central hypergeometric distribution in study ", i, "."))
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

.dnchg <- function(parms, ai, bi, ci, di, X.fit, random, verbose=FALSE, digits=4, dnchgcalc, dnchgprec, intCtrl) {

   p    <- ncol(X.fit)
   k    <- length(ai)
   b    <- parms[seq_len(p)]                  ### first p elemenets in parms are the model coefficients
   tau2 <- ifelse(random, exp(parms[p+1]), 0) ### next value is tau^2 -- optimize over exp(tau^2) value or hold at 0 if random=FALSE
   mu.i <- X.fit %*% cbind(b)

   lli  <- rep(NA_real_, k)

   if (!random) {

      for (i in seq_len(k)) {
         lli[i] <- log(.dnchgi(logOR=mu.i[i], ai=ai[i], bi=bi[i], ci=ci[i], di=di[i], random=random, dnchgcalc=dnchgcalc, dnchgprec=dnchgprec))
      }

      if (verbose)
         cat("ll =", formatC(sum(lli), digits=digits, format="f"), " ", formatC(b, digits=digits, format="f"), "\n")

   }

   if (random) {

      for (i in seq_len(k)) {

         res <- try(integrate(.dnchgi, lower=intCtrl$lower, upper=intCtrl$upper, ai=ai[i], bi=bi[i], ci=ci[i], di=di[i], mu.i=mu.i[i], tau2=tau2, random=random, dnchgcalc=dnchgcalc, dnchgprec=dnchgprec, rel.tol=intCtrl$rel.tol, subdivisions=intCtrl$subdivisions, stop.on.error=FALSE), silent=!verbose)

         if (inherits(res, "try-error")) {
            stop(paste0("Could not integrate over density of non-central hypergeometric distribution in study ", i, "."))
         } else {
            if (res$value > 0) {
               lli[i] <- log(res$value)
            } else {
               lli[i] <- -Inf
            }
         }

      }

      if (verbose)
         cat("ll = ", formatC(sum(lli), digits=digits, format="f"), " ", formatC(tau2, digits=digits, format="f"), " ", formatC(b, digits=digits, format="f"), "\n")

   }

   return(-sum(lli))

}

############################################################################

### -1 times the log likelihood (regular or restricted) for rma() location-scale models

.ll.rma.tau2 <- function(par, yi, vi, X, Z, reml, k, p, verbose, digits, REMLf) {

   #b.fe   <- par[1:p]
   #b.tau2 <- par[-c(1:p)]

   b.tau2 <- par

   ### compute predicted tau2 values and exponentiate (to force variances to be non-negative)
   tau2 <- exp(c(Z %*% b.tau2))

   ### compute weights
   wi <- 1/(vi + tau2)

   ### when using this, the optimization only pertains to the parameter(s) in 'b.tau2', as 'b.fe' is then fully
   ### determined by the current value(s) of 'b.tau2'; this is actually also how the standard RE/ME model is fitted;
   ### but is this really the best way of doing this? one could also optimize over b.fe and b.tau2 jointly
   W <- diag(wi, nrow=k, ncol=k)
   stXWX <- .invcalc(X=X, W=W, k=k)
   b.fe <- stXWX %*% crossprod(X,W) %*% as.matrix(yi)

   ### compute residual sum of squares
   RSS.f <- sum(wi*(yi - X %*% b.fe)^2)

   ### log-likelihood (could leave out additive constants)
   if (!reml) {
      llval <- -1/2 * (k)   * log(2*base::pi)                                                                             - 1/2 * sum(log(vi + tau2))                                                                   - 1/2 * RSS.f
   } else {
      llval <- -1/2 * (k-p) * log(2*base::pi) + ifelse(REMLf, 1/2 * determinant(crossprod(X), logarithm=TRUE)$modulus, 0) - 1/2 * sum(log(vi + tau2)) - 1/2 * determinant(crossprod(X,W) %*% X, logarithm=TRUE)$modulus - 1/2 * RSS.f
   }

   if (verbose) {
      cat("ll = ", ifelse(is.na(llval), NA, formatC(llval, digits=digits, format="f", flag=" ")), " ", sep="")
      cat("b.tau2 =", ifelse(is.na(b.tau2), NA, paste(formatC(b.tau2, digits=digits, format="f", flag=" "), " ", sep="")), "\n", sep="")
   }
   return(-1 * llval)

}

############################################################################

### function to compute the tetrachoric correlation coefficient and its sampling variance

.rtet <- function(ai, bi, ci, di, maxcor=.9999) {

   if (!requireNamespace("mvtnorm", quietly=TRUE))
      stop("Please install the 'mvtnorm' package to compute this measure.")

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

      #       │ ci : di    #   ci = lo X and hi Y    di = hi X and hi Y
      # var Y │∙∙∙∙:∙∙∙∙   #
      #       │ ai : bi    #   ai = lo X and lo Y    bi = hi X and lo Y
      #       ┼─────────
      #          var X
      #
      #      lo   hi
      #    ┌────┬────┐
      # lo │ ai │ bi │
      #    ├────┼────┤ var Y
      # hi │ ci │ di │
      #    └────┴────┘
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

   ### in all other cases, first optimize over rho with cut values set to sample values
   ### use suppressWarnings() to suppress "NA/Inf replaced by maximum positive value" warnings
   res <- try(suppressWarnings(optimize(fn, interval=c(-1,1), ai=ai, bi=bi, ci=ci, di=di, maxcor=maxcor, fixcut=TRUE)), silent=TRUE)

   ### check for non-convergence
   if (inherits(res, "try-error")) {
      warning("Could not estimate tetrachoric correlation coefficient.")
      return(list(yi=NA, vi=NA))
   }

   ### then use the value as the starting point and maximize over rho and the cut values
   ### (Nelder-Mead seems to do fine here; using L-BFGS-B doesn't seems to improve on this)
   res <- try(optim(par=c(res$minimum,qnorm((ai+bi)/ni),qnorm((ai+ci)/ni)), fn, ai=ai, bi=bi, ci=ci, di=di, maxcor=maxcor, fixcut=FALSE, hessian=TRUE), silent=TRUE)
   #res <- try(optim(par=c(res$minimum,qnorm((ai+bi)/ni),qnorm((ai+ci)/ni)), fn, method="L-BFGS-B", lower=c(-1,-Inf,-Inf), upper=c(1,Inf,Inf), ai=ai, bi=bi, ci=ci, di=di, maxcor=maxcor, fixcut=FALSE, hessian=TRUE), silent=TRUE)

   ### check for non-convergence
   if (inherits(res, "try-error")) {
      warning("Could not estimate tetrachoric correlation coefficient.")
      return(list(yi=NA, vi=NA))
   }

   ### take inverse of hessian and extract variance for estimate
   ### (using hessian() seems to lead to more problems, so stick with hessian from optim())
   vi <- try(chol2inv(chol(res$hessian))[1,1], silent=TRUE)
   #res$hessian <- try(chol2inv(chol(numDeriv::hessian(fn, x=res$par, ai=ai, bi=bi, ci=ci, di=di, maxcor=maxcor, fixcut=FALSE))), silent=TRUE)

   ### check for problems with computing the inverse
   if (inherits(vi, "try-error")) {
      warning("Could not estimate sampling variance of tetrachoric correlation coefficient.")
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
