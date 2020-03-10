anova.rma <- function(object, object2, btt, L, digits, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(object, "rma"))
      stop(mstyle$stop("Argument 'object' must be an object of class \"rma\"."))

   if (inherits(object, c("rma.mh", "rma.peto")))
      stop(mstyle$stop("Function not applicable for objects of class \"rma.mh\" or \"rma.peto\"."))

   if (inherits(object, "rma.glmm"))
      stop(mstyle$stop("Method not available for objects of class \"rma.glmm\"."))

   if (missing(digits)) {
      digits <- .get.digits(xdigits=object$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=object$digits, dmiss=FALSE)
   }

   if (missing(object2)) {

      ### if only 'object' has been specified, can use function to test (sets) of coefficients via
      ### the 'btt' argument or one or more linear contrasts of the coefficients via the 'L' argument

      x    <- object
      k    <- x$k
      p    <- x$p
      beta <- x$beta
      vb   <- x$vb

      if (missing(L)) {

         ### if 'L' has not been specified, then do Wald-test via 'btt' argument

         ### set/check 'btt' argument

         btt <- .set.btt(btt, p, x$int.incl, x$X)
         m <- length(btt) ### number of betas to test (m = p if all betas are tested)

         QM <- as.vector(t(beta)[btt] %*% chol2inv(chol(vb[btt,btt])) %*% beta[btt])

         if (is.element(x$test, c("knha","adhoc","t"))) {
            QM  <- QM/m
            QMp <- pf(QM, df1=m, df2=x$dfs, lower.tail=FALSE)
         } else {
            QMp <- pchisq(QM, df=m, lower.tail=FALSE)
         }

         res <- list(QM=QM, QMp=QMp, btt=btt, k=k, p=p, m=m, test=x$test, dfs=x$dfs, digits=digits, type="Wald.b")

      } else {

         ### if 'L' has been specified, then do Wald-type test(s) via 'L' argument

         if (.is.vector(L))
            L <- rbind(L)

         if (is.data.frame(L))
            L <- as.matrix(L)

         if (is.character(L))
            stop(mstyle$stop("Argument 'L' must be a numeric vector/matrix."))

         ### if model has an intercept term and L has p-1 columns, assume user left out the intercept and add it automatically

         if (x$int.incl && ncol(L) == (p-1))
            L <- cbind(1, L)

         ### if L has p+1 columns, assume that last column is the right-hand side
         ### leave this out for now; maybe add later or a 'rhs' argument (as linearHypothesis() from car package)

         #if (ncol(L) == (p+1)) {
         #   rhs <- L[,p+1]
         #   L <- L[,seq_len(p)]
         #}

         if (ncol(L) != p)
            stop(mstyle$stop("Length or number of columns of 'L' does not match number of model coefficients."))

         m <- nrow(L)

         ### test of individual hypotheses

         Lb  <- L %*% beta
         vLb <- L %*% vb %*% t(L)

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
            hyp[j] <- paste(paste(Lj[sel], rownames(beta)[sel], sep="*"), collapse=" + ") ### coefficient*variable + coefficient*variable ...
            hyp[j] <- gsub("1*", "", hyp[j], fixed=TRUE) ### turn '+1' into '+' and '-1' into '-'
            hyp[j] <- gsub("+ -", "- ", hyp[j], fixed=TRUE) ### turn '+ -' into '-'
         }
         hyp <- paste0(hyp, " = 0") ### add '= 0' at the right
         hyp <- data.frame(hyp, stringsAsFactors=FALSE)
         colnames(hyp) <- ""
         rownames(hyp) <- paste0(seq_len(m), ":") ### add '1:', '2:', ... as row names

         res <- list(QM=QM, QMp=QMp, hyp=hyp, Lb=Lb, se=se, zval=zval, pval=pval, k=k, p=p, m=m, test=x$test, dfs=x$dfs, digits=digits, type="Wald.L")

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

      ddd <- list(...)

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

      ### if they have the same number of parameters, they cannot be nested

      if (p.f == p.r)
         stop(mstyle$stop("Models have the same number of parameters. LRT not meaningful."))

      ### if p.f < p.r, then 'object' must be the reduced model and 'object2' the full model

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

      ### could do even more checks for cases where the models are clearly not nested ...

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
               warning(mstyle$warning("Models with different fixed effects. REML comparisons are not meaningful."))

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
