permutest.rma.uni <- function(x, exact=FALSE, iter=1000, btt=x$btt, permci=FALSE, progbar=TRUE, digits, control, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(x), must="rma.uni", notav=c("robust.rma", "rma.ls", "rma.gen", "rma.uni.selmodel"))

   if (missing(digits)) {
      digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
   }

   if (is.null(x$yi) || is.null(x$vi))
      stop(mstyle$stop("Information needed to carry out permutation test is not available in the model object."))

   ddd <- list(...)

   .chkdots(ddd, c("tol", "time", "seed", "verbose", "fixed", "code1", "code2"))

   if (!is.null(ddd$tol)) # in case the user specified comptol in the old manner
      comptol <- ddd$tol

   fixed <- .chkddd(ddd$fixed, FALSE, .isTRUE(ddd$fixed))

   iter <- round(iter)

   if (iter <= 1)
      stop(mstyle$stop("Argument 'iter' must be >= 2."))

   if (.isTRUE(ddd$time))
      time.start <- proc.time()

   if (!missing(btt)) {

      btt <- .set.btt(btt, x$p, x$int.incl, colnames(x$X), fixed=fixed)

      args <- list(yi=x$yi, vi=x$vi, weights=x$weights, mods=x$X, intercept=FALSE, method=x$method, weighted=x$weighted,
                   test=x$test, level=x$level, btt=btt, tau2=ifelse(x$tau2.fix, x$tau2, NA), control=x$control, skipr2=TRUE)
      x <- try(suppressWarnings(.do.call(rma.uni, args)), silent=!isTRUE(ddd$verbose))

   }

   #########################################################################
   #########################################################################
   #########################################################################

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

   if (is.character(exact) && exact == "i")
      return(X.exact.iter)

   ### if 'exact=TRUE' or if the number of iterations for an exact test are smaller
   ### than what is specified under 'iter', then carry out the exact test

   X.exact <- exact
   X.iter  <- iter

   if (X.exact || (X.exact.iter <= X.iter)) {
      X.exact <- TRUE
      X.iter <- X.exact.iter
   }

   if (X.iter == Inf)
      stop(mstyle$stop("Too many iterations required for an exact permutation test."))

   #########################################################################

   ### generate seed (needed when X.exact=FALSE)

   if (!X.exact)
      seed <- as.integer(runif(1)*2e9)

   ### set defaults for control parameters and replace with any user-defined values

   if (missing(control))
      control <- list()

   con <- list(comptol=.Machine$double.eps^0.5, tol=.Machine$double.eps^0.25,
               maxiter=100, alternative="two.sided", p2defn="abs", stat="test",
               cialt="one.sided", distfac=1, extendInt="no")
   con.pos <- pmatch(names(control), names(con))
   con[c(na.omit(con.pos))] <- control[!is.na(con.pos)]

   con$alternative <- match.arg(con$alternative, c("two.sided", "less", "greater"))
   con$p2defn      <- match.arg(con$p2defn,      c("abs", "px2"))
   con$stat        <- match.arg(con$stat,        c("test", "coef"))

   if (exists("comptol", inherits=FALSE))
      con$comptol <- comptol

   if (!X.exact) {
      if (!is.null(ddd$seed)) {
         set.seed(ddd$seed)
      } else {
         set.seed(seed)
      }
   }

   ### elements that need to be returned

   outlist <- "beta=beta, zval=zval, QM=QM"

   #########################################################################

   if (progbar)
      cat(mstyle$verbose(paste0("Running ", X.iter, " iterations for an ", ifelse(X.exact, "exact", "approximate"), " permutation test.\n")))

   if (!is.null(ddd[["code1"]]))
      eval(expr = parse(text = ddd[["code1"]]))

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

            if (!is.null(ddd[["code2"]]))
               eval(expr = parse(text = ddd[["code2"]]))

            args <- list(yi=signmat[i,]*x$yi, vi=x$vi, weights=x$weights, intercept=TRUE, method=x$method, weighted=x$weighted,
                         test=x$test, level=x$level, btt=1, tau2=ifelse(x$tau2.fix, x$tau2, NA), control=x$control, skipr2=TRUE, outlist=outlist)
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

            if (!is.null(ddd[["code2"]]))
               eval(expr = parse(text = ddd[["code2"]]))

            signs <- sample(c(-1,1), x$k, replace=TRUE) # easier to understand (a tad slower for small k, but faster for larger k)
            #signs <- 2*rbinom(x$k,1,0.5)-1

            args <- list(yi=signs*x$yi, vi=x$vi, weights=x$weights, intercept=TRUE, method=x$method, weighted=x$weighted,
                         test=x$test, level=x$level, btt=1, tau2=ifelse(x$tau2.fix, x$tau2, NA), control=x$control, skipr2=TRUE, outlist=outlist)
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

   #########################################################################

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

            if (!is.null(ddd[["code2"]]))
               eval(expr = parse(text = ddd[["code2"]]))

            args <- list(yi=x$yi, vi=x$vi, weights=x$weights, mods=cbind(X[permmat[i,],]), intercept=FALSE, method=x$method, weighted=x$weighted,
                         test=x$test, level=x$level, btt=x$btt, tau2=ifelse(x$tau2.fix, x$tau2, NA), control=x$control, skipr2=TRUE, outlist=outlist)
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

            if (!is.null(ddd[["code2"]]))
               eval(expr = parse(text = ddd[["code2"]]))

            args <- list(yi=x$yi, vi=x$vi, weights=x$weights, mods=cbind(X[sample(x$k),]), intercept=FALSE, method=x$method, weighted=x$weighted,
                         test=x$test, level=x$level, btt=x$btt, tau2=ifelse(x$tau2.fix, x$tau2, NA), control=x$control, skipr2=TRUE, outlist=outlist)
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

   #########################################################################
   #########################################################################
   #########################################################################

   ### permutation-based CI

   ci.lb <- x$ci.lb
   ci.ub <- x$ci.ub

   if (.isTRUE(permci) || is.numeric(permci)) {

      level <- .level(x$level)

      ### check if it is even possible to reject at level

      if (1/X.iter > level / ifelse(con$cialt == "one.sided", 1, 2)) {

         permci <- FALSE
         warning(mstyle$warning(paste0("Cannot obtain ", 100*(1-x$level), "% permutation-based CI; number of permutations (", X.iter, ") too low.")), call.=FALSE)

      } else {

         ### if permci is numeric, check if existing coefficients have been specified
         ### otherwise, CIs will be obtained for all model coefficients

         if (is.numeric(permci)) {
            coefs <- unique(round(permci))
            if (any(coefs > x$p) || any(coefs < 1))
               stop(mstyle$stop("Non-existent coefficients specified via 'permci'."))
            permci <- TRUE
         } else {
            coefs <- seq_len(x$p)
         }

         ci.lb <- rep(NA_real_, x$p)
         ci.ub <- rep(NA_real_, x$p)

         for (j in coefs) {

            if (progbar)
               cat(mstyle$verbose(paste0("Searching for lower CI bound of coefficient ", j, ": \n")))

            if (con$cialt == "one.sided") {
               con$alternative <- "greater"
            } else {
               con$alternative <- "two.sided"
            }

            tmp <- try(uniroot(.permci, interval=c(x$ci.lb[j] - con$distfac*(x$beta[j,1] - x$ci.lb[j]), x$beta[j,1]), extendInt=ifelse(con$extendInt == "no", "no", "upX"), tol=con$tol, maxiter=con$maxiter, obj=x, j=j, exact=X.exact, iter=X.iter, progbar=progbar, level=level, digits=digits, control=con)$root, silent=TRUE)

            if (inherits(tmp, "try-error")) {
               ci.lb[j] <- NA_real_
            } else {
               ci.lb[j] <- tmp
            }

            if (progbar)
               cat(mstyle$verbose(paste0("Searching for upper CI bound of coefficient ", j, ": \n")))

            if (con$cialt == "one.sided") {
               con$alternative <- "less"
            } else {
               con$alternative <- "two.sided"
            }

            tmp <- try(uniroot(.permci, interval=c(x$beta[j,1], x$ci.ub[j] + con$distfac*(x$ci.ub[j] - x$beta[j,1])), extendInt=ifelse(con$extendInt == "no", "no", "downX"), tol=con$tol, maxiter=con$maxiter, obj=x, j=j, exact=X.exact, iter=X.iter, progbar=progbar, level=level, digits=digits, control=con)$root, silent=TRUE)

            if (inherits(tmp, "try-error")) {
               ci.ub[j] <- NA_real_
            } else {
               ci.ub[j] <- tmp
            }

         }

      }

   }

   #########################################################################

   out <- list(pval=pval, QMdf=x$QMdf, QMp=QMp, beta=x$beta, se=x$se, zval=x$zval, ci.lb=ci.lb, ci.ub=ci.ub, QM=x$QM,
               k=x$k, p=x$p, btt=x$btt, m=x$m, test=x$test, dfs=x$dfs, ddf=x$ddf,
               int.only=x$int.only, int.incl=x$int.incl,
               digits=digits, exact.iter=X.exact.iter,
               permci=permci, alternative=con$alternative, p2defn=con$p2defn, stat=con$stat)

   out$skip.beta <- FALSE
   out$QM.perm <- QM.perm
   out$zval.perm <- data.frame(zval.perm)
   out$beta.perm <- data.frame(beta.perm)
   names(out$zval.perm) <- names(out$beta.perm) <- colnames(x$X)

   if (.isTRUE(ddd$time)) {
      time.end <- proc.time()
      .print.time(unname(time.end - time.start)[3])
   }

   class(out) <- "permutest.rma.uni"
   return(out)

}
