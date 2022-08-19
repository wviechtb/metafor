permutest.rma.ls <- function(x, exact=FALSE, iter=1000, progbar=TRUE, digits, control, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="rma.ls")

   if (missing(digits)) {
      digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
   }

   ddd <- list(...)

   .chkdots(ddd, c("tol", "time", "seed", "verbose", "retpermdist", "permci", "skip.beta", "skip.alpha"))

   if (!is.null(ddd$tol)) # in case user specifies comptol in the old manner
      comptol <- ddd$tol

   if (.isTRUE(ddd$permci))
      warning(mstyle$warning("Permutation-based CIs for location-scale models not currently available."), call.=FALSE)

   if (.isTRUE(ddd$time))
      time.start <- proc.time()

   if (.isTRUE(ddd$skip.beta)) {
      skip.beta <- TRUE
   } else {
      skip.beta <- FALSE
   }

   if (.isTRUE(ddd$skip.alpha)) {
      skip.alpha <- TRUE
   } else {
      skip.alpha <- FALSE
   }

   iter <- round(iter)

   if (iter <= 1)
      stop(mstyle$stop("Argument 'iter' must be >= 2."))

   ### for intercept-only models, cannot run a permutation test

   if (x$Z.int.only) {
      skip.alpha <- TRUE
      warning(mstyle$warning("Cannot carry out a permutation test for an intercept-only scale model."), call.=FALSE)
   }

   if (skip.beta && skip.alpha)
      stop(mstyle$stop("Must run permutation test for at least one part of the model."))

   ### set control parameters and possibly replace with user-defined values

   if (missing(control))
      control <- list()

   con <- list(comptol=.Machine$double.eps^0.5, alternative="two.sided", p2defn="abs", stat="test")
   con.pos <- pmatch(names(control), names(con))
   con[c(na.omit(con.pos))] <- control[!is.na(con.pos)]

   con$alternative <- match.arg(con$alternative, c("two.sided", "less", "greater"))
   con$p2defn      <- match.arg(con$p2defn,      c("abs", "px2"))
   con$stat        <- match.arg(con$stat,        c("test", "coef"))

   if (exists("comptol", inherits=FALSE))
      con$comptol <- comptol

   #########################################################################
   #########################################################################
   #########################################################################

   if (!skip.beta) {

      ### calculate number of permutations for an exact permutation test

      if (x$int.only) {

         ### for intercept-only models, there are 2^k possible permutations of the signs

         X.exact.iter <- 2^x$k

      } else {

         ### for meta-regression models, there are k! possible permutations of the rows of the model matrix

         #X.exact.iter <- round(exp(lfactorial(x$k))) # note: without round(), not exactly an integer!

         ### however, when there are duplicated rows in the model matrix, the number of *unique* permutations
         ### is lower; the code below below determines the number of unique permutations

         ### order the X matrix

         X <- as.data.frame(x$X)[do.call(order, as.data.frame(x$X)),]

         ### determine groupings

         X.indices <- cumsum(c(TRUE, !duplicated(X)[-1]))

         ### this turns 1,1,1,2,2,3,4,4,4 into 1,1,1,4,4,6,7,7,7 so that the actual row numbers can be permuted

         X.indices <- rep(cumsum(rle(X.indices)$lengths) - (rle(X.indices)$lengths - 1), rle(X.indices)$lengths)

         ### determine exact number of unique permutations

         ind.table <- table(X.indices)

         X.exact.iter <- round(prod((max(ind.table)+1):x$k) / prod(factorial(ind.table[-which.max(ind.table)]))) # cancel largest value in numerator and denominator to reduce overflow problems
         #X.exact.iter <- round(factorial(x$k) / prod(factorial(ind.table)))       # definitional formula
         #X.exact.iter <- round(exp(lfactorial(x$k) - sum(lfactorial(ind.table)))) # using log of definitional formula and then round(exp())

         if (is.na(X.exact.iter))
            X.exact.iter <- Inf

      }

      ### if 'exact=TRUE' or if the number of iterations for an exact test are smaller
      ### than what is specified under 'iter', then carry out the exact test

      X.exact <- exact
      X.iter  <- iter

      if (X.exact || (X.exact.iter <= X.iter)) {
         X.exact <- TRUE
         X.iter <- X.exact.iter
      }

      if (X.iter == Inf)
         stop(mstyle$stop("Too many iterations required for an exact permutation test of the location model."))

      ######################################################################

      ### generate seed (needed when X.exact=FALSE)

      if (!X.exact) {

         seed <- as.integer(runif(1)*2e9)

         if (!is.null(ddd$seed)) {
            set.seed(ddd$seed)
         } else {
            set.seed(seed)
         }

      }

      ### elements that need to be returned

      outlist <- "beta=beta, zval=zval, QM=QM"

      ######################################################################

      if (progbar)
         cat(mstyle$verbose(paste0("Running ", X.iter, " iterations for an ", ifelse(X.exact, "exact", "approximate"), " permutation test of the location model.\n")))

      if (x$int.only) {

         ### permutation test for intercept-only model

         zval.perm <- try(rep(NA_real_, X.iter), silent=TRUE)

         if (inherits(zval.perm, "try-error"))
            stop(mstyle$stop("Number of iterations requested too large."))

         beta.perm <- try(rep(NA_real_, X.iter), silent=TRUE)

         if (inherits(beta.perm, "try-error"))
            stop(mstyle$stop("Number of iterations requested too large."))

         QM.perm <- try(rep(NA_real_, X.iter), silent=TRUE)

         if (inherits(QM.perm, "try-error"))
            stop(mstyle$stop("Number of iterations requested too large."))

         if (progbar)
            pbar <- pbapply::startpb(min=0, max=X.iter)

         if (X.exact) { # exact permutation test for intercept-only models

            signmat <- as.matrix(expand.grid(replicate(x$k, list(c(1,-1))), KEEP.OUT.ATTRS=FALSE))

            for (i in seq_len(X.iter)) {

               args <- list(yi=signmat[i,]*x$yi, vi=x$vi, weights=x$weights, intercept=TRUE, scale=x$Z, link=x$link, method=x$method, weighted=x$weighted, test=x$test, level=x$level, btt=1, alpha=ifelse(x$alpha.fix, x$alpha, NA), optbeta=x$optbeta, beta=ifelse(x$beta.fix, x$beta, NA), control=x$control, skiphes=TRUE, outlist=outlist)
               res <- try(suppressWarnings(.do.call(rma.uni, args)), silent=!isTRUE(ddd$verbose))

               if (inherits(res, "try-error"))
                  next

               beta.perm[i] <- res$beta[,1]
               zval.perm[i] <- res$zval
               QM.perm[i]   <- res$QM

               if (progbar)
                  pbapply::setpb(pbar, i)

            }

         } else { # approximate permutation test for intercept-only models

            i <- 1

            while (i <= X.iter) {

               signs <- sample(c(-1,1), x$k, replace=TRUE) # easier to understand (a tad slower for small k, but faster for larger k)
               #signs <- 2*rbinom(x$k,1,.5)-1

               args <- list(yi=signs*x$yi, vi=x$vi, weights=x$weights, intercept=TRUE, scale=x$Z, link=x$link, method=x$method, weighted=x$weighted, test=x$test, level=x$level, btt=1, alpha=ifelse(x$alpha.fix, x$alpha, NA), optbeta=x$optbeta, beta=ifelse(x$beta.fix, x$beta, NA), control=x$control, skiphes=TRUE, outlist=outlist)
               res <- try(suppressWarnings(.do.call(rma.uni, args)), silent=!isTRUE(ddd$verbose))

               if (inherits(res, "try-error"))
                  next

               beta.perm[i] <- res$beta[,1]
               zval.perm[i] <- res$zval
               QM.perm[i]   <- res$QM

               i <- i + 1

               if (progbar)
                  pbapply::setpb(pbar, i)

            }

         }

         ### the first random permutation is always the observed data (avoids possibility of p=0)

         if (!X.exact) {
            beta.perm[1] <- x$beta[,1]
            zval.perm[1] <- x$zval
            QM.perm[1]   <- x$QM
         }

         if (con$alternative == "two.sided") {

            if (con$p2defn == "abs") {

               ### absolute value definition of the two-sided p-value

               if (con$stat == "test") {
                  pval <- mean(abs(zval.perm) >= abs(x$zval) - con$comptol, na.rm=TRUE) # based on test statistic
               } else {
                  pval <- mean(abs(beta.perm) >= abs(c(x$beta)) - con$comptol, na.rm=TRUE) # based on coefficient
               }

            } else {

               ### two times the one-sided p-value definition of the two-sided p-value

               if (con$stat == "test") {
                  if (x$zval > median(zval.perm, na.rm=TRUE)) {
                     pval <- 2*mean(zval.perm >= x$zval - con$comptol, na.rm=TRUE) # based on test statistic
                  } else {
                     pval <- 2*mean(zval.perm <= x$zval + con$comptol, na.rm=TRUE)
                  }
               } else {
                  if (c(x$beta) > median(beta.perm, na.rm=TRUE)) {
                     pval <- 2*mean(beta.perm >= c(x$beta) - con$comptol, na.rm=TRUE) # based on coefficient
                  } else {
                     pval <- 2*mean(beta.perm <= c(x$beta) + con$comptol, na.rm=TRUE)
                  }
               }

            }

         }

         if (con$alternative == "less") {
            if (con$stat == "test") {
               pval <- mean(zval.perm <= x$zval + con$comptol, na.rm=TRUE) # based on test statistic
            } else {
               pval <- mean(beta.perm <= c(x$beta) + con$comptol, na.rm=TRUE) # based on coefficient
            }
         }

         if (con$alternative == "greater") {
            if (con$stat == "test") {
               pval <- mean(zval.perm >= x$zval - con$comptol, na.rm=TRUE) # based on test statistic
            } else {
               pval <- mean(beta.perm >= c(x$beta) - con$comptol, na.rm=TRUE) # based on coefficient
            }
         }

         pval[pval > 1] <- 1

         QMp <- mean(QM.perm >= x$QM - con$comptol, na.rm=TRUE)

      ######################################################################

      } else {

         ### permutation test for meta-regression model

         zval.perm <- try(suppressWarnings(matrix(NA_real_, nrow=X.iter, ncol=x$p)), silent=TRUE)

         if (inherits(zval.perm, "try-error"))
            stop(mstyle$stop("Number of iterations requested too large."))

         beta.perm <- try(suppressWarnings(matrix(NA_real_, nrow=X.iter, ncol=x$p)), silent=TRUE)

         if (inherits(beta.perm, "try-error"))
            stop(mstyle$stop("Number of iterations requested too large."))

         QM.perm <- try(rep(NA_real_, X.iter), silent=TRUE)

         if (inherits(QM.perm, "try-error"))
            stop(mstyle$stop("Number of iterations requested too large."))

         if (progbar)
            pbar <- pbapply::startpb(min=0, max=X.iter)

         if (X.exact) { # exact permutation test for meta-regression models

            #permmat <- .genperms(x$k)
            permmat <- .genuperms(X.indices) # use recursive algorithm to obtain all unique permutations

            for (i in seq_len(X.iter)) {

               args <- list(yi=x$yi, vi=x$vi, weights=x$weights, mods=cbind(X[permmat[i,],]), intercept=FALSE, scale=x$Z, link=x$link, method=x$method, weighted=x$weighted, test=x$test, level=x$level, btt=x$btt, alpha=ifelse(x$alpha.fix, x$alpha, NA), optbeta=x$optbeta, beta=ifelse(x$beta.fix, x$beta, NA), control=x$control, skiphes=FALSE, outlist=outlist)
               res <- try(suppressWarnings(.do.call(rma.uni, args)), silent=!isTRUE(ddd$verbose))

               if (inherits(res, "try-error"))
                  next

               beta.perm[i,] <- res$beta[,1]
               zval.perm[i,] <- res$zval
               QM.perm[i]    <- res$QM

               if (progbar)
                  pbapply::setpb(pbar, i)

            }

         } else { # approximate permutation test for meta-regression models

            i <- 1

            while (i <= X.iter) {

               args <- list(yi=x$yi, vi=x$vi, weights=x$weights, mods=cbind(X[sample(x$k),]), intercept=FALSE, scale=x$Z, link=x$link, method=x$method, weighted=x$weighted, test=x$test, level=x$level, btt=x$btt, alpha=ifelse(x$alpha.fix, x$alpha, NA), optbeta=x$optbeta, beta=ifelse(x$beta.fix, x$beta, NA), control=x$control, skiphes=FALSE, outlist=outlist)
               res <- try(suppressWarnings(.do.call(rma.uni, args)), silent=!isTRUE(ddd$verbose))

               if (inherits(res, "try-error"))
                  next

               beta.perm[i,] <- res$beta[,1]
               zval.perm[i,] <- res$zval
               QM.perm[i]    <- res$QM

               i <- i + 1

               if (progbar)
                  pbapply::setpb(pbar, i)

            }

         }

         ### the first random permutation is always the observed data (avoids possibility of p=0)

         if (!X.exact) {
            beta.perm[1,] <- x$beta[,1]
            zval.perm[1,] <- x$zval
            QM.perm[1]    <- x$QM
         }

         if (con$alternative == "two.sided") {

            if (con$p2defn == "abs") {

               ### absolute value definition of the two-sided p-value

               if (con$stat == "test") {
                  pval <- rowMeans(t(abs(zval.perm)) >= abs(x$zval) - con$comptol, na.rm=TRUE) # based on test statistics
               } else {
                  pval <- rowMeans(t(abs(beta.perm)) >= abs(c(x$beta)) - con$comptol, na.rm=TRUE) # based on coefficients
               }

            } else {

               ### two times the one-sided p-value definition of the two-sided p-value

               pval <- rep(NA_real_, x$p)

               if (con$stat == "test") {
                  for (j in seq_len(x$p)) {
                     if (x$zval[j] > median(zval.perm[,j], na.rm=TRUE)) {
                        pval[j] <- 2*mean(zval.perm[,j] >= x$zval[j] - con$comptol, na.rm=TRUE)
                     } else {
                        pval[j] <- 2*mean(zval.perm[,j] <= x$zval[j] + con$comptol, na.rm=TRUE)
                     }
                  }
               } else {
                  for (j in seq_len(x$p)) {
                     if (c(x$beta)[j] > median(beta.perm[,j], na.rm=TRUE)) {
                        pval[j] <- 2*mean(beta.perm[,j] >= c(x$beta)[j] - con$comptol, na.rm=TRUE)
                     } else {
                        pval[j] <- 2*mean(beta.perm[,j] <= c(x$beta)[j] + con$comptol, na.rm=TRUE)
                     }
                  }
               }

            }

         }

         if (con$alternative == "less") {
            if (con$stat == "test") {
               pval <- rowMeans(t(zval.perm) <= x$zval + con$comptol, na.rm=TRUE) # based on test statistics
            } else {
               pval <- rowMeans(t(beta.perm) <= c(x$beta) + con$comptol, na.rm=TRUE) # based on coefficients
            }
         }

         if (con$alternative == "greater") {
            if (con$stat == "test") {
               pval <- rowMeans(t(zval.perm) >= x$zval - con$comptol, na.rm=TRUE) # based on test statistics
            } else {
               pval <- rowMeans(t(beta.perm) >= c(x$beta) - con$comptol, na.rm=TRUE) # based on coefficients
            }
         }

         pval[pval > 1] <- 1

         QMp <- mean(QM.perm >= x$QM - con$comptol, na.rm=TRUE)

      }

      if (progbar)
         pbapply::closepb(pbar)

   } else {

      beta.perm <- NA
      zval.perm <- NA
      QM.perm <- NA
      pval <- x$pval
      QMp <- x$QMp
      X.exact.iter <- 0

   }

   #########################################################################
   #########################################################################
   #########################################################################

   if (!skip.alpha) {

      ### calculate number of permutations for an exact permutation test

      Z <- as.data.frame(x$Z)[do.call(order, as.data.frame(x$Z)),]
      Z.indices <- cumsum(c(TRUE, !duplicated(Z)[-1]))
      Z.indices <- rep(cumsum(rle(Z.indices)$lengths) - (rle(Z.indices)$lengths - 1), rle(Z.indices)$lengths)
      ind.table <- table(Z.indices)
      Z.exact.iter <- round(prod((max(ind.table)+1):x$k) / prod(factorial(ind.table[-which.max(ind.table)])))

      if (is.na(Z.exact.iter))
         Z.exact.iter <- Inf

      Z.exact <- exact
      Z.iter  <- iter

      if (Z.exact || (Z.exact.iter <= Z.iter)) {
         Z.exact <- TRUE
         Z.iter <- Z.exact.iter
      }

      if (Z.iter == Inf)
         stop(mstyle$stop("Too many iterations required for an exact permutation test of the scale model."))

      #########################################################################

      ### generate seed (needed when Z.exact=FALSE)

      if (!Z.exact) {

         seed <- as.integer(runif(1)*2e9)

         if (!is.null(ddd$seed)) {
            set.seed(ddd$seed)
         } else {
            set.seed(seed)
         }

      }

      ### elements that need to be returned

      outlist <- "alpha=alpha, zval.alpha=zval.alpha, QS=QS"

      #########################################################################

      if (progbar)
         cat(mstyle$verbose(paste0("Running ", Z.iter, " iterations for an ", ifelse(Z.exact, "exact", "approximate"), " permutation test of the scale model.\n")))

      ### permutation test for the scale model

      zval.alpha.perm <- try(suppressWarnings(matrix(NA_real_, nrow=Z.iter, ncol=x$q)), silent=TRUE)

      if (inherits(zval.alpha.perm, "try-error"))
         stop(mstyle$stop("Number of iterations requested too large."))

      alpha.perm <- try(suppressWarnings(matrix(NA_real_, nrow=Z.iter, ncol=x$q)), silent=TRUE)

      if (inherits(alpha.perm, "try-error"))
         stop(mstyle$stop("Number of iterations requested too large."))

      QS.perm <- try(rep(NA_real_, Z.iter), silent=TRUE)

      if (inherits(QS.perm, "try-error"))
         stop(mstyle$stop("Number of iterations requested too large."))

      if (progbar)
         pbar <- pbapply::startpb(min=0, max=Z.iter)

      if (Z.exact) { # exact permutation test for meta-regression models

         #permmat <- .genperms(x$k)
         permmat <- .genuperms(Z.indices) # use recursive algorithm to obtain all unique permutations

         for (i in seq_len(Z.iter)) {

            args <- list(yi=x$yi, vi=x$vi, weights=x$weights, mods=x$X, intercept=FALSE, scale=cbind(Z[permmat[i,],]), link=x$link, method=x$method, weighted=x$weighted, test=x$test, level=x$level, att=x$att, alpha=ifelse(x$alpha.fix, x$alpha, NA), optbeta=x$optbeta, beta=ifelse(x$beta.fix, x$beta, NA), control=x$control, skiphes=FALSE, outlist=outlist)
            res <- try(suppressWarnings(.do.call(rma.uni, args)), silent=!isTRUE(ddd$verbose))

            if (inherits(res, "try-error"))
               next

            alpha.perm[i,] <- res$alpha[,1]
            zval.alpha.perm[i,] <- res$zval.alpha
            QS.perm[i] <- res$QS

            if (progbar)
               pbapply::setpb(pbar, i)

         }

      } else { # approximate permutation test for meta-regression models

         i <- 1

         while (i <= Z.iter) {

            args <- list(yi=x$yi, vi=x$vi, weights=x$weights, mods=x$X, intercept=FALSE, scale=cbind(Z[sample(x$k),]), link=x$link, method=x$method, weighted=x$weighted, test=x$test, level=x$level, att=x$att, alpha=ifelse(x$alpha.fix, x$alpha, NA), optbeta=x$optbeta, beta=ifelse(x$beta.fix, x$beta, NA), control=x$control, skiphes=FALSE, outlist=outlist)
            res <- try(suppressWarnings(.do.call(rma.uni, args)), silent=!isTRUE(ddd$verbose))

            if (inherits(res, "try-error"))
               next

            alpha.perm[i,] <- res$alpha[,1]
            zval.alpha.perm[i,] <- res$zval.alpha
            QS.perm[i] <- res$QS

            i <- i + 1

            if (progbar)
               pbapply::setpb(pbar, i)

         }

      }

      ### the first random permutation is always the observed data (avoids possibility of p=0)

      if (!Z.exact) {
         alpha.perm[1,] <- x$alpha[,1]
         zval.alpha.perm[1,] <- x$zval.alpha
         QS.perm[1] <- x$QS
      }

      if (con$alternative == "two.sided") {

         if (con$p2defn == "abs") {

            ### absolute value definition of the two-sided p-value

            if (con$stat == "test") {
               pval.alpha <- rowMeans(t(abs(zval.alpha.perm)) >= abs(x$zval.alpha) - con$comptol, na.rm=TRUE) # based on test statistics
            } else {
               pval.alpha <- rowMeans(t(abs(alpha.perm)) >= abs(c(x$alpha)) - con$comptol, na.rm=TRUE) # based on coefficients
            }

         } else {

            ### two times the one-sided p-value definition of the two-sided p-value

            pval.alpha <- rep(NA_real_, x$q)

            if (con$stat == "test") {
               for (j in seq_len(x$q)) {
                  if (x$zval.alpha[j] > median(zval.alpha.perm[,j], na.rm=TRUE)) {
                     pval.alpha[j] <- 2*mean(zval.alpha.perm[,j] >= x$zval.alpha.[j] - con$comptol, na.rm=TRUE)
                  } else {
                     pval.alpha[j] <- 2*mean(zval.alpha.perm[,j] <= x$zval.alpha.[j] + con$comptol, na.rm=TRUE)
                  }
               }
            } else {
               for (j in seq_len(x$q)) {
                  if (c(x$alpha)[j] > median(alpha.perm[,j], na.rm=TRUE)) {
                     pval.alpha[j] <- 2*mean(alpha.perm[,j] >= c(x$alpha)[j] - con$comptol, na.rm=TRUE)
                  } else {
                     pval.alpha[j] <- 2*mean(alpha.perm[,j] <= c(x$alpha)[j] + con$comptol, na.rm=TRUE)
                  }
               }
            }

         }

      }

      if (con$alternative == "less") {
         if (con$stat == "test") {
            pval.alpha <- rowMeans(t(zval.alpha.perm) <= x$zval.alpha + con$comptol, na.rm=TRUE) # based on test statistics
         } else {
            pval.alpha <- rowMeans(t(alpha.perm) <= c(x$alpha) + con$comptol, na.rm=TRUE) # based on coefficients
         }
      }

      if (con$alternative == "greater") {
         if (con$stat == "test") {
            pval.alpha <- rowMeans(t(zval.alpha.perm) >= x$zval.alpha - con$comptol, na.rm=TRUE)    # based on test statistics
         } else {
            pval.alpha <- rowMeans(t(alpha.perm) >= c(x$alpha) - con$comptol, na.rm=TRUE) # based on coefficients
         }
      }

      pval.alpha[pval.alpha > 1] <- 1
      pval.alpha[x$alpha.fix] <- NA

      QSp <- mean(QS.perm >= x$QS - con$comptol, na.rm=TRUE)

      if (progbar)
         pbapply::closepb(pbar)

   } else {

      alpha.perm <- NA
      zval.alpha.perm <- NA
      QS.perm <- NA
      pval.alpha <- x$pval.alpha
      QSp <- NA
      Z.exact.iter <- 0

   }

   ############################################################################
   ############################################################################
   ############################################################################

   out <- list(pval=pval, QMdf=x$QMdf, QMp=QMp, beta=x$beta, se=x$se, zval=x$zval, ci.lb=x$ci.lb, ci.ub=x$ci.ub, QM=x$QM,
               pval.alpha=pval.alpha, QSdf=x$QSdf, QSp=QSp, alpha=x$alpha, se.alpha=x$se.alpha, zval.alpha=x$zval.alpha, ci.lb.alpha=x$ci.lb.alpha, ci.ub.alpha=x$ci.ub.alpha, QS=x$QS,
               k=x$k, p=x$p, btt=x$btt, m=x$m, test=x$test, dfs=x$dfs, ddf=x$ddf,
                      q=x$q, att=x$att, m.alpha=x$m.alpha, ddf.alpha=x$ddf.alpha,
               int.only=x$int.only, int.incl=x$int.incl, Z.int.only=x$Z.int.only, Z.int.incl=x$Z.int.incl,
               digits=digits, exact.iter=X.exact.iter, Z.exact.iter=Z.exact.iter,
               permci=FALSE, alternative=con$alternative, p2defn=con$p2defn, stat=con$stat)

   out$skip.beta <- skip.beta
   out$QM.perm <- QM.perm
   out$zval.perm <- data.frame(zval.perm)
   out$beta.perm <- data.frame(beta.perm)
   if (!skip.beta)
      names(out$zval.perm) <- names(out$beta.perm) <- colnames(x$X)

   out$skip.alpha <- skip.alpha
   out$QS.perm <- QS.perm
   out$zval.alpha.perm <- data.frame(zval.alpha.perm)
   out$alpha.perm <- data.frame(alpha.perm)
   if (!skip.alpha)
      names(out$zval.alpha.perm) <- names(out$alpha.perm) <- colnames(x$Z)

   if (.isTRUE(ddd$time)) {
      time.end <- proc.time()
      .print.time(unname(time.end - time.start)[3])
   }

   class(out) <- c("permutest.rma.ls", "permutest.rma.uni")
   return(out)

}
