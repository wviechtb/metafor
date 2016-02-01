permutest.rma.uni <- function(x, exact=FALSE, iter=1000, progbar=TRUE, retpermdist=FALSE, digits, tol, ...) {

   if (!inherits(x, "rma.uni"))
      stop("Argument 'x' must be an object of class \"rma.uni\".")

   if (missing(digits))
      digits <- x$digits

   if (missing(tol))
      tol <- .Machine$double.eps^0.5

   #########################################################################

   if (x$int.only) {

      exact.iter <- 2^x$k

   } else {

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

      ### total number of permutations
      #exact.iter <- round(exp(lfactorial(x$k))) ### note: without round(), not exactly an integer!

   }

   if (exact || (exact.iter <= iter)) {
      exact <- TRUE
      iter  <- exact.iter
   }

   if (iter == Inf)
      stop("Too many iterations required for exact permutation test.\n")

   if (progbar)
      cat("Running ", iter, " iterations for ", ifelse(exact, "exact", "approximate"), " permutation test.\n", sep="")

   #########################################################################

   if (x$int.only) {

      ### permutation test for intercept-only models

      zval.perm <- try(rep(NA_real_, iter), silent=TRUE)

      if (inherits(zval.perm, "try-error"))
         stop("Number of iterations requested too large.")

      QM.perm <- try(rep(NA_real_, iter), silent=TRUE)

      if (inherits(QM.perm, "try-error"))
         stop("Number of iterations requested too large.")

      if (progbar)
         pbar <- txtProgressBar(min=0, max=iter, style=3)

      if (exact) { ### exact permutation test for intercept-only models

         signmat <- as.matrix(expand.grid(replicate(x$k, list(c(1,-1))), KEEP.OUT.ATTRS=FALSE))

         for (i in seq_len(iter)) {

            res <- try(suppressWarnings(rma.uni(signmat[i,]*x$yi, x$vi, weights=x$weights, intercept=TRUE, method=x$method, weighted=x$weighted, knha=x$knha, btt=1, tau2=ifelse(x$tau2.fix, x$tau2, NA), control=x$control)), silent=FALSE)

            if (inherits(res, "try-error"))
               next

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
            res <- try(suppressWarnings(rma.uni(signs*x$yi, x$vi, weights=x$weights, intercept=TRUE, method=x$method, weighted=x$weighted, knha=x$knha, btt=1, tau2=ifelse(x$tau2.fix, x$tau2, NA), control=x$control)), silent=FALSE)

            if (inherits(res, "try-error"))
               next

            zval.perm[i] <- res$zval
            QM.perm[i]   <- res$QM
            i <- i + 1

            if (progbar)
               setTxtProgressBar(pbar, i)

         }

      }

      ### the first random permutation is always the observed test statistic (avoids possibility of p=0)

      if (!exact) {
         zval.perm[1] <- x$zval
         QM.perm[1]   <- x$QM
      }

      ### absolute value definition of the p-value

      #zval.orig.abs <- abs(x$zval)
      #pval <- mean(abs(zval.perm) >= zval.orig.abs, na.rm=TRUE)

      ### two times the one-sided p-value definition

      if (x$zval > 0) {
         pval <- 2*mean(zval.perm >= x$zval - tol, na.rm=TRUE)
      } else {
         pval <- 2*mean(zval.perm <= x$zval + tol, na.rm=TRUE)
      }
      pval[pval > 1] <- 1

      QMp <- mean(QM.perm >= x$QM - tol, na.rm=TRUE)

   #########################################################################

   } else {

      ### permutation test for models with moderators

      zval.perm <- try(suppressWarnings(matrix(NA_real_, nrow=iter, ncol=x$p)), silent=TRUE)

      if (inherits(zval.perm, "try-error"))
         stop("Number of iterations requested too large.")

      QM.perm <- try(rep(NA_real_, iter), silent=TRUE)

      if (inherits(QM.perm, "try-error"))
         stop("Number of iterations requested too large.")

      if (progbar)
         pbar <- txtProgressBar(min=0, max=iter, style=3)

      if (exact) { ### exact permutation test for models with moderators

         #permmat <- .genperms(x$k)
         permmat <- .genuperms(indices)

         for (i in seq_len(iter)) {

            res <- try(suppressWarnings(rma.uni(x$yi, x$vi, weights=x$weights, mods=cbind(X[permmat[i,],]), intercept=FALSE, method=x$method, weighted=x$weighted, knha=x$knha, btt=x$btt, tau2=ifelse(x$tau2.fix, x$tau2, NA), control=x$control)), silent=FALSE)

            if (inherits(res, "try-error"))
               next

            zval.perm[i,] <- res$zval
            QM.perm[i]    <- res$QM

            if (progbar)
               setTxtProgressBar(pbar, i)

         }

      } else { ### approximate permutation test for models with moderators

         i <- 1

         while (i <= iter) {

            res <- try(suppressWarnings(rma.uni(x$yi, x$vi, weights=x$weights, mods=cbind(X[sample(x$k),]), intercept=FALSE, method=x$method, weighted=x$weighted, knha=x$knha, btt=x$btt, tau2=ifelse(x$tau2.fix, x$tau2, NA), control=x$control)), silent=FALSE)

            if (inherits(res, "try-error"))
               next

            zval.perm[i,] <- res$zval
            QM.perm[i]    <- res$QM
            i <- i + 1

            if (progbar)
               setTxtProgressBar(pbar, i)

         }

      }

      ### the first random permutation is always the observed test statistic (avoids possibility of p=0)

      if (!exact) {
         zval.perm[1,] <- x$zval
         QM.perm[1]    <- x$QM
      }

      ### absolute value definition of the p-value

      #zval.orig.abs <- abs(x$zval)
      #pval <- rowMeans(t(abs(zval.perm)) >= zval.orig.abs, na.rm=TRUE)

      ### two times the one-sided p-value definition

      pval <- rep(NA_real_, x$p)
      for (j in seq_len(x$p)) {
         if (x$zval[j] > 0) {
            pval[j] <- 2*mean(zval.perm[,j] >= x$zval[j] - tol, na.rm=TRUE)
         } else {
            pval[j] <- 2*mean(zval.perm[,j] <= x$zval[j] + tol, na.rm=TRUE)
         }
      }
      pval[pval > 1] <- 1

      QMp <- mean(QM.perm >= x$QM - tol, na.rm=TRUE)

   }

   #########################################################################

   if (progbar)
      close(pbar)

   out <- list(pval=pval, QMp=QMp, b=x$b, se=x$se, zval=x$zval, ci.lb=x$ci.lb, ci.ub=x$ci.ub, QM=x$QM, k=x$k, p=x$p, btt=x$btt, m=x$m, knha=x$knha, dfs=x$dfs, int.only=x$int.only, digits=digits, exact.iter=exact.iter)

   if (retpermdist) {
      out$QM.perm   <- QM.perm
      out$zval.perm <- data.frame(zval.perm)
      names(out$zval.perm) <- colnames(x$X)
   }

   class(out) <- "permutest.rma.uni"
   return(out)

}
