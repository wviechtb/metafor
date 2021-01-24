anova.rma <- function(object, object2, btt, L, digits, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(object), must="rma", notap=c("rma.mh", "rma.peto"), notav="rma.glmm")

   if (missing(digits)) {
      digits <- .get.digits(xdigits=object$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=object$digits, dmiss=FALSE)
   }

   ddd <- list(...)

   .chkdots(ddd, c("test", "att", "K"))

   if (!is.null(ddd$att)) {
      if (!inherits(object, "rma.ls"))
         stop(mstyle$stop("Can only specify 'att' for location-scale models."))
      att <- ddd$att
   } else {
      att <- NULL
   }

   if (!is.null(ddd$K)) {
      if (!inherits(object, "rma.ls"))
         stop(mstyle$stop("Can only specify 'K' for location-scale models."))
      K <- ddd$K
   } else {
      K <- NULL
   }

   if (missing(object2)) {

      ### if only 'object' has been specified, can use function to test (sets) of coefficients via
      ### the 'btt' (or 'att') argument or one or more linear contrasts of the coefficients via the
      ### 'L' (or 'K') argument

      x <- object

      if (missing(L) && is.null(K)) {

         ### if 'L' (and 'K') has not been specified, then do a Wald-test via the 'btt' argument (can also use 'att' for location-scale models)

         if (inherits(object, "rma.ls") && !is.null(att)) {

            if (!missing(btt))
               stop(mstyle$stop("Can only specify either 'btt' or 'att', but not both."))

            ### set/check 'att' argument

            att <- .set.btt(att, x$q, x$Z.int.incl, colnames(x$Z))
            m <- length(att)

            QM <- try(as.vector(t(x$alpha)[att] %*% chol2inv(chol(x$vb.alpha[att,att])) %*% x$alpha[att]), silent=TRUE)

            if (inherits(QM, "try-error"))
               QM <- NA

            if (is.element(x$test, c("t"))) {
               QM <- QM / m
               QMp <- pf(QM, df1=m, df2=x$dfs.alpha, lower.tail=FALSE)
            } else {
               QMp <- pchisq(QM, df=m, lower.tail=FALSE)
            }

            res <- list(QM=QM, QMp=QMp, att=att, k=x$k, q=x$q, m=m, test=x$test, dfs=x$dfs.alpha, digits=digits, type="Wald.att")

         } else {

            ### set/check 'btt' argument

            btt <- .set.btt(btt, x$p, x$int.incl, colnames(x$X))
            m <- length(btt)

            QM <- try(as.vector(t(x$beta)[btt] %*% chol2inv(chol(x$vb[btt,btt])) %*% x$beta[btt]), silent=TRUE)

            if (inherits(QM, "try-error"))
               QM <- NA

            if (is.element(x$test, c("knha","adhoc","t"))) {
               QM  <- QM/m
               QMp <- pf(QM, df1=m, df2=x$dfs, lower.tail=FALSE)
            } else {
               QMp <- pchisq(QM, df=m, lower.tail=FALSE)
            }

            res <- list(QM=QM, QMp=QMp, btt=btt, k=x$k, p=x$p, m=m, test=x$test, dfs=x$dfs, digits=digits, type="Wald.btt")

         }

      } else {

         if (inherits(object, "rma.ls") && !is.null(K)) {

            ### if 'K' has been specified, then do Wald-type test(s) via 'K' argument

            if (!missing(L))
               stop(mstyle$stop("Can only specify either 'L' or 'K', but not both."))

            if (.is.vector(K))
               K <- rbind(K)

            if (is.data.frame(K))
               K <- as.matrix(K)

            if (is.character(K))
               stop(mstyle$stop("Argument 'K' must be a numeric vector/matrix."))

            ### if model has an intercept term and K has q-1 columns, assume user left out the intercept and add it automatically

            if (x$Z.int.incl && ncol(K) == (x$q-1))
               K <- cbind(1, K)

            ### if K has q+1 columns, assume that last column is the right-hand side
            ### leave this out for now; maybe add later or a 'rhs' argument (as linearHypothesis() from car package)

            #if (ncol(K) == (x$q+1)) {
            #   rhs <- K[,x$q+1]
            #   K <- K[,seq_len(x$q)]
            #}

            if (ncol(K) != x$q)
               stop(mstyle$stop(paste0("Length or number of columns of 'K' (", ncol(K), ") does not match number of scale coefficients (", x$q, ").")))

            m <- nrow(K)

            ### test of individual hypotheses

            Ka  <- K %*% x$alpha
            vKa <- K %*% x$vb.alpha %*% t(K)

            se <- sqrt(diag(vKa))
            zval <- c(Ka/se)

            if (is.element(x$test, c("t"))) {
               pval <- 2*pt(abs(zval), df=x$dfs.alpha, lower.tail=FALSE)
            } else {
               pval <- 2*pnorm(abs(zval), lower.tail=FALSE)
            }

            ### omnibus test of all hypotheses (only possible if 'K' is of full rank)

            QM  <- NA ### need this in case QM cannot be calculated below
            QMp <- NA ### need this in case QMp cannot be calculated below

            if (rankMatrix(K) == m) {

               QM <- try(as.vector(t(Ka) %*% chol2inv(chol(vKa)) %*% Ka), silent=TRUE)

               if (inherits(QM, "try-error"))
                  QM <- NA

               if (is.element(x$test, c("t"))) {
                  QM  <- QM/m
                  QMp <- pf(QM, df1=m, df2=x$dfs.alpha, lower.tail=FALSE)
               } else {
                  QMp <- pchisq(QM, df=m, lower.tail=FALSE)
               }

            }

            ### create a data frame with each row specifying the linear combination tested

            hyp <- rep("", m)
            for (j in seq_len(m)) {
               Lj <- round(K[j,], digits=digits[["est"]]) ### coefficients for the jth contrast
               sel <- Lj != 0 ### TRUE if coefficient is != 0
               hyp[j] <- paste(paste(Lj[sel], rownames(x$alpha)[sel], sep="*"), collapse=" + ") ### coefficient*variable + coefficient*variable ...
               hyp[j] <- gsub("1*", "", hyp[j], fixed=TRUE) ### turn '+1' into '+' and '-1' into '-'
               hyp[j] <- gsub("+ -", "- ", hyp[j], fixed=TRUE) ### turn '+ -' into '-'
            }
            hyp <- paste0(hyp, " = 0") ### add '= 0' at the right
            hyp <- data.frame(hyp, stringsAsFactors=FALSE)
            colnames(hyp) <- ""
            rownames(hyp) <- paste0(seq_len(m), ":") ### add '1:', '2:', ... as row names

            res <- list(QM=QM, QMp=QMp, hyp=hyp, Ka=Ka, se=se, zval=zval, pval=pval, k=x$k, q=x$q, m=m, test=x$test, dfs=x$dfs, digits=digits, type="Wald.Ka")

         } else {

            ### if 'L' has been specified, then do Wald-type test(s) via 'L' argument

            if (.is.vector(L))
               L <- rbind(L)

            if (is.data.frame(L))
               L <- as.matrix(L)

            if (is.character(L))
               stop(mstyle$stop("Argument 'L' must be a numeric vector/matrix."))

            ### if model has an intercept term and L has p-1 columns, assume user left out the intercept and add it automatically

            if (x$int.incl && ncol(L) == (x$p-1))
               L <- cbind(1, L)

            ### if L has p+1 columns, assume that last column is the right-hand side
            ### leave this out for now; maybe add later or a 'rhs' argument (as linearHypothesis() from car package)

            #if (ncol(L) == (x$p+1)) {
            #   rhs <- L[,x$p+1]
            #   L <- L[,seq_len(x$p)]
            #}

            if (ncol(L) != x$p)
               stop(mstyle$stop(paste0("Length or number of columns of 'L' (", ncol(L), ") does not match number of model coefficients (", x$p, ").")))

            m <- nrow(L)

            ### test of individual hypotheses

            Lb  <- L %*% x$beta
            vLb <- L %*% x$vb %*% t(L)

            se <- sqrt(diag(vLb))
            zval <- c(Lb/se)

            if (is.element(x$test, c("knha","adhoc","t"))) {
               pval <- 2*pt(abs(zval), df=x$dfs, lower.tail=FALSE)
            } else {
               pval <- 2*pnorm(abs(zval), lower.tail=FALSE)
            }

            ### omnibus test of all hypotheses (only possible if 'L' is of full rank)

            QM  <- NA ### need this in case QM cannot be calculated below
            QMp <- NA ### need this in case QMp cannot be calculated below

            if (rankMatrix(L) == m) {

               ### use try(), since this could fail: this could happen when the var-cov matrix of the
               ### fixed effects has been estimated using robust() -- 'vb' is then only guaranteed to
               ### be positive semidefinite, so for certain linear combinations, vLb could be singular
               ### (see Cameron & Miller, 2015, p. 326)

               QM <- try(as.vector(t(Lb) %*% chol2inv(chol(vLb)) %*% Lb), silent=TRUE)

               if (inherits(QM, "try-error"))
                  QM <- NA

               if (is.element(x$test, c("knha","adhoc","t"))) {
                  QM  <- QM/m
                  QMp <- pf(QM, df1=m, df2=x$dfs, lower.tail=FALSE)
               } else {
                  QMp <- pchisq(QM, df=m, lower.tail=FALSE)
               }

            }

            ### create a data frame with each row specifying the linear combination tested

            hyp <- rep("", m)
            for (j in seq_len(m)) {
               Lj <- round(L[j,], digits=digits[["est"]]) ### coefficients for the jth contrast
               sel <- Lj != 0 ### TRUE if coefficient is != 0
               hyp[j] <- paste(paste(Lj[sel], rownames(x$beta)[sel], sep="*"), collapse=" + ") ### coefficient*variable + coefficient*variable ...
               hyp[j] <- gsub("1*", "", hyp[j], fixed=TRUE) ### turn '+1' into '+' and '-1' into '-'
               hyp[j] <- gsub("+ -", "- ", hyp[j], fixed=TRUE) ### turn '+ -' into '-'
            }
            hyp <- paste0(hyp, " = 0") ### add '= 0' at the right
            hyp <- data.frame(hyp, stringsAsFactors=FALSE)
            colnames(hyp) <- ""
            rownames(hyp) <- paste0(seq_len(m), ":") ### add '1:', '2:', ... as row names

            res <- list(QM=QM, QMp=QMp, hyp=hyp, Lb=Lb, se=se, zval=zval, pval=pval, k=x$k, p=x$p, m=m, test=x$test, dfs=x$dfs, digits=digits, type="Wald.Lb")

         }

      }

   } else {

      ### if 'object' and 'object2' have been specified, can use function to
      ### do model comparisons via a likelihood ratio test (and fit indices)

      if (!inherits(object2, "rma"))
         stop(mstyle$stop("Argument 'object2' must be an object of class \"rma\"."))

      if (inherits(object2, c("rma.mh","rma.peto")))
         stop(mstyle$stop("Function not applicable for objects of class \"rma.mh\" or \"rma.peto\"."))

      if (inherits(object2, "rma.glmm"))
         stop(mstyle$stop("Method not available for objects of class \"rma.glmm\"."))

      if (!identical(class(object), class(object2)))
         stop(mstyle$stop("Class of 'object1' must be the same as class of 'object2'."))

      if (!is.null(ddd$test)) {
         test <- match.arg(ddd$test, c("LRT", "Wald"))
      } else {
         test <- "LRT"
      }

      ### assume 'object' is the full model and 'object2' the reduced model

      m.f <- object
      m.r <- object2

      ### number of parameters in the models

      p.f <- m.f$parms
      p.r <- m.r$parms

      ### check if they have the same number of parameters

      if (p.f == p.r)
         stop(mstyle$stop("Models have the same number of parameters. LRT not meaningful."))

      ### if p.f < p.r, then let 'object' be the reduced model and 'object2' the full model

      if (p.f < p.r) {
         m.f <- object2
         m.r <- object
         p.f <- m.f$parms
         p.r <- m.r$parms
      }

      ### check if models are based on the same data (TODO: also check for same weights?)

      if (inherits(object, "rma.uni")) {
         if (!(identical(as.vector(m.f$yi), as.vector(m.r$yi)) && identical(as.vector(m.f$vi), as.vector(m.r$vi)))) ### as.vector() to strip attributes/names
            stop(mstyle$stop("Observed outcomes and/or sampling variances not equal in the full and reduced model."))
      }

      if (inherits(object, "rma.mv")) {
         if (!(identical(as.vector(m.f$yi), as.vector(m.r$yi)) && identical(as.matrix(m.f$V), as.matrix(m.r$V)))) ### as.vector() to strip attributes/names, as.matrix() to make both V matrices non-sparse
            stop(mstyle$stop("Observed outcomes and/or sampling variances/covariances not equal in the full and reduced model."))
      }

      ### for Wald-type test, both models should be fitted using the same method

      if (test == "Wald" && (m.f$method != m.r$method))
         stop(mstyle$stop("Full and reduced model do not use the same 'method'."))

      ### for LRT, reduced model may use method="FE" and full model method="(RE)ML"
      ### which is fine, but the other way around doesn't really make sense

      if (m.f$method == "FE" && m.r$method != "FE")
         stop(mstyle$stop("Full model uses a fixed- and reduced model uses random/mixed-effects model."))

      ### could do even more checks for cases where the models are clearly not nested

      ######################################################################

      ### for 'rma.uni' objects, calculate pseudo R^2 value (based on the
      ### proportional reduction in tau^2) comparing full vs. reduced model

      if (inherits(object, "rma.uni") && !inherits(object, "rma.ls") && !inherits(object2, "rma.ls")) {
         if (m.f$method == "FE") {
            R2 <- NA
         } else if (identical(m.r$tau2,0)) {
            R2 <- 0
         } else {
            R2 <- 100 * max(0, (m.r$tau2 - m.f$tau2)/m.r$tau2)
         }
      } else {
         R2 <- NA
      }

      ### for 'rma.uni' objects, extract tau^2 estimates

      if (inherits(object, "rma.uni") && !inherits(object, "rma.ls") && !inherits(object2, "rma.ls")) {
         tau2.f <- m.f$tau2
         tau2.r <- m.r$tau2
      } else {
         tau2.f <- NA
         tau2.r <- NA
      }

      if (test == "LRT") {

         p.diff <- p.f - p.r

         if (m.f$method == "REML") {

            LRT <- m.r$fit.stats["dev","REML"] - m.f$fit.stats["dev","REML"]
            fit.stats.f <- t(m.f$fit.stats)["REML",] # to keep (row)names of fit.stats
            fit.stats.r <- t(m.r$fit.stats)["REML",] # to keep (row)names of fit.stats

            if (!identical(m.f$X, m.r$X))
               warning(mstyle$warning("Models with different fixed effects. REML comparisons are not meaningful."), call.=FALSE)

            ### in this case, one could consider just taking the ML deviances, but this
            ### is really ad-hoc; there is some theory in Welham & Thompson (1997) about
            ### LRTs for fixed effects when using REML estimation, but this involves
            ### additional work

         } else {

            LRT <- m.r$fit.stats["dev","ML"] - m.f$fit.stats["dev","ML"]
            fit.stats.f <- t(m.f$fit.stats)["ML",]
            fit.stats.r <- t(m.r$fit.stats)["ML",]

            ### in principle, a 'rma.uni' model may have been fitted with something besides ML
            ### estimation (for tau^2), so technically the deviance is then not really that as
            ### obtained with proper ML estimation; issue a warning about this?

         }

         ### set LRT to 0 if LRT < 0 (this should not happen, but could due to numerical issues)

         LRT[LRT < 0] <- 0

         pval <- pchisq(LRT, df=p.diff, lower.tail=FALSE)

         res <- list(fit.stats.f=fit.stats.f, fit.stats.r=fit.stats.r, p.f=p.f, p.r=p.r,
                     LRT=LRT, pval=pval, QE.f=m.f$QE, QE.r=m.r$QE,
                     tau2.f=tau2.f, tau2.r=tau2.r, R2=R2,
                     method=m.f$method, class.f=class(m.f), digits=digits, type="LRT")

      }

      if (test == "Wald") {

         btt <- setdiff(colnames(m.f$X), colnames(m.r$X))

         if (length(setdiff(colnames(m.r$X), colnames(m.f$X))) != 0L)
            stop(mstyle$stop("There are coefficients in the reduced model that are not in the full model."))

         btt <- charmatch(btt, colnames(m.f$X))

         if (any(is.na(btt)))
            stop(mstyle$stop("Cannot identify coefficients to test."))

         res <- anova(m.f, btt=btt)
         return(res)

      }


   }

   class(res) <- "anova.rma"
   return(res)

}
