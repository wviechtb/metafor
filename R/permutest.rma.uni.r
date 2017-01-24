permutest.rma.uni <- function(x, exact=FALSE, iter=1000, permci=FALSE, progbar=TRUE, retpermdist=FALSE, digits, tol, control, ...) {

   if (!inherits(x, "rma.uni"))
      stop("Argument 'x' must be an object of class \"rma.uni\".")

   if (inherits(x, "robust.rma"))
      stop("Method not yet implemented for objects of class \"robust.rma\". Sorry!")

   if (inherits(x, "rma.ls"))
      stop("Method not yet implemented for objects of class \"rma.ls\". Sorry!")

   if (missing(digits))
      digits <- x$digits

   ### default numerical tolerance

   if (missing(tol)) {
      comp.tol <- .Machine$double.eps^0.5
   } else {
      comp.tol <- tol
   }

   #########################################################################

   ### calculate number of permutations for an exact permutation test

   if (x$int.only) {

      ### for intercept-only models, 2^k possible permutations of the signs

      exact.iter <- 2^x$k

   } else {

      ### for meta-regression models, there are k! possible permutations of the rows of the model matrix

      #exact.iter <- round(exp(lfactorial(x$k))) ### note: without round(), not exactly an integer!

      ### however, when there are duplicated rows in the model matrix, the number of *unique* permutations
      ### is lower; the code below below determines the number of unique permutations

      ### order the X matrix

      X <- as.data.frame(x$X)[do.call(order, as.data.frame(x$X)),]

      ### determine groupings

      indices <- cumsum(c(TRUE, !duplicated(X)[-1]))

      ### this turns 1,1,1,2,2,3,4,4,4 into 1,1,1,4,4,6,7,7,7 so that the actual row numbers can be permutated

      indices <- rep(cumsum(rle(indices)$length) - (rle(indices)$length - 1), rle(indices)$length)

      ### determine exact number of unique permutations
      ### code below cancels largest ind.table value, which reduces overflow problems

      ind.table  <- table(indices)

      exact.iter <- round(prod((max(ind.table)+1):x$k) / prod(factorial(ind.table[-which.max(ind.table)]))) ### cancel largest value in numerator and denominator
      #exact.iter <- round(factorial(x$k) / prod(factorial(ind.table)))       ### definitional formula
      #exact.iter <- round(exp(lfactorial(x$k) - sum(lfactorial(ind.table)))) ### using log of definitional formula and then exp()

   }

   ### if 'exact=TRUE' or if the number of iterations for an exact test are smaller
   ### than what is specified under 'iter', then carry out the exact test

   if (exact || (exact.iter <= iter)) {
      exact <- TRUE
      iter <- exact.iter
   }

   if (iter == Inf)
      stop("Too many iterations required for exact permutation test.\n")

   #########################################################################

   ### generate seed (needed when exact=FALSE)

   if (!exact) {
      seed <- as.integer(runif(1)*2e9)
   } else {
      seed <- NA
   }

   ### set control parameters for uniroot() and possibly replace with user-defined values

   if (missing(control))
      control <- list()

   con <- list(tol=.Machine$double.eps^0.25, maxiter=100, alternative="two.sided", p2defn="abs", stat="test", cialt="one.sided", seed=seed, distfac=1)
   con[pmatch(names(control), names(con))] <- control

   if (!exact)
      set.seed(con$seed)

   #########################################################################

   if (progbar)
      cat("Running ", iter, " iterations for ", ifelse(exact, "exact", "approximate"), " permutation test.\n", sep="")

   if (x$int.only) {

      ### permutation test for intercept-only models

      zval.perm <- try(rep(NA_real_, iter), silent=TRUE)

      if (inherits(zval.perm, "try-error"))
         stop("Number of iterations requested too large.")

      beta.perm <- try(rep(NA_real_, iter), silent=TRUE)

      if (inherits(beta.perm, "try-error"))
         stop("Number of iterations requested too large.")

      QM.perm <- try(rep(NA_real_, iter), silent=TRUE)

      if (inherits(QM.perm, "try-error"))
         stop("Number of iterations requested too large.")

      if (progbar)
         pbar <- txtProgressBar(min=0, max=iter, style=3)

      if (exact) { ### exact permutation test for intercept-only models

         signmat <- as.matrix(expand.grid(replicate(x$k, list(c(1,-1))), KEEP.OUT.ATTRS=FALSE))

         for (i in seq_len(iter)) {

            res <- try(suppressWarnings(rma.uni(signmat[i,]*x$yi, x$vi, weights=x$weights, intercept=TRUE, method=x$method, weighted=x$weighted, test=x$test, btt=1, tau2=ifelse(x$tau2.fix, x$tau2, NA), control=x$control)), silent=FALSE)

            if (inherits(res, "try-error"))
               next

            beta.perm[i] <- coef(res)
            zval.perm[i] <- res$zval
            QM.perm[i]   <- res$QM

            if (progbar)
               setTxtProgressBar(pbar, i)

         }

      } else { ### approximate permutation test for intercept-only models

         i <- 1

         while (i <= iter) {

            signs <- sample(c(-1,1), x$k, replace=TRUE) ### easier to understand (a tad slower for small k, but faster for larger k)
            #signs <- 2*rbinom(x$k,1,.5)-1

            #if (i == 2) print(sample(c(-1,1), x$k, replace=TRUE))

            res <- try(suppressWarnings(rma.uni(signs*x$yi, x$vi, weights=x$weights, intercept=TRUE, method=x$method, weighted=x$weighted, test=x$test, btt=1, tau2=ifelse(x$tau2.fix, x$tau2, NA), control=x$control)), silent=FALSE)

            if (inherits(res, "try-error"))
               next

            beta.perm[i] <- coef(res)
            zval.perm[i] <- res$zval
            QM.perm[i]   <- res$QM

            i <- i + 1

            if (progbar)
               setTxtProgressBar(pbar, i)

         }

      }

      ### the first random permutation is always the observed data (avoids possibility of p=0)

      if (!exact) {
         beta.perm[1] <- coef(x)
         zval.perm[1] <- x$zval
         QM.perm[1]   <- x$QM
      }

      if (con$alternative == "two.sided") {

         if (con$p2defn == "abs") {

            ### absolute value definition of the two-sided p-value

            if (con$stat == "test") {
               pval <- mean(abs(zval.perm) >= abs(x$zval) - comp.tol, na.rm=TRUE) ### based on test statistic
            } else {
               pval <- mean(abs(beta.perm) >= abs(c(x$beta)) - comp.tol, na.rm=TRUE) ### based on coefficient
            }

         } else {

            ### two times the one-sided p-value definition of the two-sided p-value

            if (con$stat == "test") {
               if (x$zval > 0) {
                  pval <- 2*mean(zval.perm >= x$zval - comp.tol, na.rm=TRUE) ### based on test statistic
               } else {
                  pval <- 2*mean(zval.perm <= x$zval + comp.tol, na.rm=TRUE) ### based on coefficient
               }
            } else {
               if (c(x$beta) > 0) {
                  pval <- 2*mean(beta.perm >= c(x$beta) - comp.tol, na.rm=TRUE) ### based on test statistic
               } else {
                  pval <- 2*mean(beta.perm <= c(x$beta) + comp.tol, na.rm=TRUE) ### based on coefficient
               }
            }

         }

      }

      if (con$alternative == "less") {
         if (con$stat == "test") {
            pval <- mean(zval.perm <= x$zval + comp.tol, na.rm=TRUE) ### based on test statistic
         } else {
            pval <- mean(beta.perm <= c(x$beta) + comp.tol, na.rm=TRUE) ### based on coefficient
         }
      }

      if (con$alternative == "greater") {
         if (con$stat == "test") {
            pval <- mean(zval.perm >= x$zval - comp.tol, na.rm=TRUE) ### based on test statistic
         } else {
            pval <- mean(beta.perm >= c(x$beta) - comp.tol, na.rm=TRUE) ### based on coefficient
         }
      }

      pval[pval > 1] <- 1

      QMp <- mean(QM.perm >= x$QM - comp.tol, na.rm=TRUE)

   #########################################################################

   } else {

      ### permutation test for meta-regression models

      zval.perm <- try(suppressWarnings(matrix(NA_real_, nrow=iter, ncol=x$p)), silent=TRUE)

      if (inherits(zval.perm, "try-error"))
         stop("Number of iterations requested too large.")

      beta.perm <- try(suppressWarnings(matrix(NA_real_, nrow=iter, ncol=x$p)), silent=TRUE)

      if (inherits(beta.perm, "try-error"))
         stop("Number of iterations requested too large.")

      QM.perm <- try(rep(NA_real_, iter), silent=TRUE)

      if (inherits(QM.perm, "try-error"))
         stop("Number of iterations requested too large.")

      if (progbar)
         pbar <- txtProgressBar(min=0, max=iter, style=3)

      if (exact) { ### exact permutation test for meta-regression models

         #permmat <- .genperms(x$k)
         permmat <- .genuperms(indices) ### use recursive algorithm to obtain all unique permutations

         for (i in seq_len(iter)) {

            res <- try(suppressWarnings(rma.uni(x$yi, x$vi, weights=x$weights, mods=cbind(X[permmat[i,],]), intercept=FALSE, method=x$method, weighted=x$weighted, test=x$test, btt=x$btt, tau2=ifelse(x$tau2.fix, x$tau2, NA), control=x$control)), silent=FALSE)

            if (inherits(res, "try-error"))
               next

            beta.perm[i,] <- coef(res)
            zval.perm[i,] <- res$zval
            QM.perm[i]    <- res$QM

            if (progbar)
               setTxtProgressBar(pbar, i)

         }

      } else { ### approximate permutation test for meta-regression models

         i <- 1

         while (i <= iter) {

            res <- try(suppressWarnings(rma.uni(x$yi, x$vi, weights=x$weights, mods=cbind(X[sample(x$k),]), intercept=FALSE, method=x$method, weighted=x$weighted, test=x$test, btt=x$btt, tau2=ifelse(x$tau2.fix, x$tau2, NA), control=x$control)), silent=FALSE)

            if (inherits(res, "try-error"))
               next

            beta.perm[i,] <- coef(res)
            zval.perm[i,] <- res$zval
            QM.perm[i]    <- res$QM

            i <- i + 1

            if (progbar)
               setTxtProgressBar(pbar, i)

         }

      }

      ### the first random permutation is always the observed data (avoids possibility of p=0)

      if (!exact) {
         beta.perm[1,] <- coef(x)
         zval.perm[1,] <- x$zval
         QM.perm[1]    <- x$QM
      }

      if (con$alternative == "two.sided") {

         if (con$p2defn == "abs") {

            ### absolute value definition of the two-sided p-value

            if (con$stat == "test") {
               pval <- rowMeans(t(abs(zval.perm)) >= abs(x$zval) - comp.tol, na.rm=TRUE) ### based on test statistics
            } else {
               pval <- rowMeans(t(abs(beta.perm)) >= abs(c(x$beta)) - comp.tol, na.rm=TRUE) ### based on coefficients
            }

         } else {

            ### two times the one-sided p-value definition of the two-sided p-value

            pval <- rep(NA_real_, x$p)

            if (con$stat == "test") {
               for (j in seq_len(x$p)) {
                  if (x$zval[j] > 0) {
                     pval[j] <- 2*mean(zval.perm[,j] >= x$zval[j] - comp.tol, na.rm=TRUE)
                  } else {
                     pval[j] <- 2*mean(zval.perm[,j] <= x$zval[j] + comp.tol, na.rm=TRUE)
                  }
               }
            } else {
               for (j in seq_len(x$p)) {
                  if (c(x$beta)[j] > 0) {
                     pval[j] <- 2*mean(beta.perm[,j] >= c(x$beta)[j] - comp.tol, na.rm=TRUE)
                  } else {
                     pval[j] <- 2*mean(beta.perm[,j] <= c(x$beta)[j] + comp.tol, na.rm=TRUE)
                  }
               }
            }

         }

      }

      if (con$alternative == "less") {
         if (con$stat == "test") {
            pval <- rowMeans(t(zval.perm) <= x$zval + comp.tol, na.rm=TRUE) ### based on test statistics
         } else {
            pval <- rowMeans(t(beta.perm) <= c(x$beta) + comp.tol, na.rm=TRUE) ### based on coefficients
         }
      }

      if (con$alternative == "greater") {
         if (con$stat == "test") {
            pval <- rowMeans(t(zval.perm) >= x$zval - comp.tol, na.rm=TRUE) ### based on test statistics
         } else {
            pval <- rowMeans(t(beta.perm) >= c(x$beta) - comp.tol, na.rm=TRUE) ### based on coefficients
         }
      }

      pval[pval > 1] <- 1

      QMp <- mean(QM.perm >= x$QM - comp.tol, na.rm=TRUE)

   }

   if (progbar)
      close(pbar)

   #########################################################################

   ### permutation-based CI

   ci.lb <- x$ci.lb
   ci.ub <- x$ci.ub

   if ((is.logical(permci) && permci) || is.numeric(permci)) {

      level <- ifelse(x$level > 1, (100-x$level)/100, ifelse(x$level > .5, 1-x$level, x$level))

      ### check if it is even possible to reject at level

      if (1/iter > level / ifelse(con$cialt == "one.sided", 1, 2)) {

         permci <- FALSE
         warning("Cannot obtain ", x$level, "% permutation-based CI; number of permutations (", iter, ") too low.")

      } else {

         ### if permci is numeric, check if existing coefficients have been specified
         ### otherwise, CIs will be obtained for all model coefficients

         if (is.numeric(permci)) {
            coefs <- unique(round(permci))
            if (any(coefs > x$p) || any(coefs < 1))
               stop("Non-existent coefficients specified via 'permci'.")
            permci <- TRUE
         } else {
            coefs <- seq_len(x$p)
         }

         ci.lb <- rep(NA, x$p)
         ci.ub <- rep(NA, x$p)

         for (j in coefs) {

            if (progbar)
               cat("Searching for lower CI bound of coefficient ", j, ": \n", sep="")

            if (con$cialt == "one.sided") {
               con$alternative <- "greater"
            } else {
               con$alternative <- "two.sided"
            }

            #tmp <- try(uniroot(.permci, interval=c(x$ci.lb[j], coef(x)[j]), extendInt="upX", tol=con$tol, maxiter=con$maxiter, obj=x, j=j, exact=exact, iter=iter, progbar=progbar, comp.tol=comp.tol, level=level, digits=digits, control=con)$root, silent=TRUE)
            tmp <- try(uniroot(.permci, interval=c(x$ci.lb[j] - con$distfac*(coef(x)[j] - x$ci.lb[j]), coef(x)[j]), extendInt="no", tol=con$tol, maxiter=con$maxiter, obj=x, j=j, exact=exact, iter=iter, progbar=progbar, comp.tol=comp.tol, level=level, digits=digits, control=con)$root, silent=TRUE)

            if (inherits(tmp, "try-error")) {
               ci.lb[j] <- NA
            } else {
               ci.lb[j] <- tmp
            }

            if (progbar)
               cat("Searching for upper CI bound of coefficient ", j, ": \n", sep="")

            if (con$cialt == "one.sided") {
               con$alternative <- "less"
            } else {
               con$alternative <- "two.sided"
            }

            #tmp <- try(uniroot(.permci, interval=c(coef(x)[j], x$ci.ub[j]), extendInt="downX", tol=con$tol, maxiter=con$maxiter, obj=x, j=j, exact=exact, iter=iter, progbar=progbar, comp.tol=comp.tol, level=level, digits=digits, control=con)$root, silent=TRUE)
            tmp <- try(uniroot(.permci, interval=c(coef(x)[j], x$ci.ub[j] + con$distfac*(x$ci.ub[j] - coef(x)[j])), extendInt="no", tol=con$tol, maxiter=con$maxiter, obj=x, j=j, exact=exact, iter=iter, progbar=progbar, comp.tol=comp.tol, level=level, digits=digits, control=con)$root, silent=TRUE)

            if (inherits(tmp, "try-error")) {
               ci.ub[j] <- NA
            } else {
               ci.ub[j] <- tmp
            }

         }

      }

   }

   #########################################################################

   out <- list(pval=pval, QMp=QMp, beta=x$beta, se=x$se, zval=x$zval, ci.lb=ci.lb, ci.ub=ci.ub,
               QM=x$QM, k=x$k, p=x$p, btt=x$btt, m=x$m, test=x$test, dfs=x$dfs,
               int.only=x$int.only, digits=digits, exact.iter=exact.iter, permci=permci)

   if (retpermdist) {
      out$QM.perm   <- QM.perm
      out$zval.perm <- data.frame(zval.perm)
      out$beta.perm <- data.frame(beta.perm)
      names(out$zval.perm) <- names(out$beta.perm) <- colnames(x$X)
   }

   class(out) <- "permutest.rma.uni"
   return(out)

}
