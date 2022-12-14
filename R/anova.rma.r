anova.rma <- function(object, object2, btt, X, att, Z, rhs, digits, refit=FALSE, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(object), must="rma", notap=c("rma.mh", "rma.peto"), notav="rma.glmm")

   if (missing(digits)) {
      digits <- .get.digits(xdigits=object$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=object$digits, dmiss=FALSE)
   }

   ddd <- list(...)

   .chkdots(ddd, c("test", "L", "verbose", "fixed"))

   if (!is.null(ddd$L))
      X <- ddd$L

   if (is.null(ddd$fixed)) {
      fixed <- FALSE
   } else {
      fixed <- .isTRUE(ddd$fixed)
   }

   if (!missing(att) && !inherits(object, "rma.ls"))
      stop(mstyle$stop("Can only specify 'att' for location-scale models."))

   if (!missing(Z) && !inherits(object, "rma.ls"))
      stop(mstyle$stop("Can only specify 'Z' for location-scale models."))

   if (missing(object2)) {

      ### if only 'object' has been specified, can use function to test one or multiple coefficients
      ### via the 'btt' (or 'att') argument or one or more linear contrasts of the coefficients via
      ### the 'X' (or 'Z') argument

      x <- object

      if (missing(X) && missing(Z)) {

         ### if 'X' (and 'Z') has not been specified, then do a Wald-test via the 'btt' argument (can also use 'att' for location-scale models)

         if (inherits(object, "rma.ls") && !missing(att)) {

            if (!missing(btt))
               stop(mstyle$stop("Can only specify either 'btt' or 'att', but not both."))

            ### set/check 'att' argument

            if (missing(att) || is.null(att)) {
               att <- x$att
            } else {
               if (is.character(att) && length(att) > 1L)
                  att <- as.list(att)
               if (is.list(att)) {
                  if (!missing(rhs))
                     stop(mstyle$stop("Cannot use 'rhs' argument when specifying a list for 'att'."))
                  sav <- lapply(att, function(attj) anova(x, att=attj, digits=digits, fixed=fixed))
                  names(sav) <- sapply(att, .format.btt)
                  class(sav) <- "list.anova.rma"
                  return(sav)
               }
               att <- .set.btt(att, x$q, x$Z.int.incl, colnames(x$Z), fixed=fixed)
            }

            m <- length(att)

            if (missing(rhs)) {
               rhs <- rep(0, m)
            } else {
               if (length(rhs) == 1L)
                  rhs <- rep(rhs, m)
               if (length(rhs) != m)
                  stop(mstyle$stop(paste0("Length of 'rhs' (", length(rhs), ") does not match the number of coefficients tested (", m, ").")))
            }

            x$alpha[att,] <- x$alpha[att,] - rhs

            QS <- try(as.vector(t(x$alpha)[att] %*% chol2inv(chol(x$va[att,att])) %*% x$alpha[att]), silent=TRUE)

            if (inherits(QS, "try-error"))
               QS <- NA

            if (is.element(x$test, c("knha","adhoc","t"))) {
               QS   <- QS / m
               QSdf <- c(m, x$QSdf[2])
               QSp  <- pf(QS, df1=QSdf[1], df2=QSdf[2], lower.tail=FALSE)
            } else {
               QSdf <- c(m, NA)
               QSp  <- pchisq(QS, df=QSdf[1], lower.tail=FALSE)
            }

            res <- list(QS=QS, QSdf=QSdf, QSp=QSp, att=att, k=x$k, q=x$q, m=m, test=x$test, digits=digits, type="Wald.att")

         } else {

            ### set/check 'btt' argument

            if (missing(btt) || is.null(btt)) {
               btt <- x$btt
            } else {
               if (is.character(btt) && length(btt) > 1L)
                  btt <- as.list(btt)
               if (is.list(btt)) {
                  if (!missing(rhs))
                     stop(mstyle$stop("Cannot use 'rhs' argument when specifying a list for 'btt'."))
                  sav <- lapply(btt, function(bttj) anova(x, btt=bttj, digits=digits, fixed=fixed))
                  names(sav) <- sapply(btt, .format.btt)
                  class(sav) <- "list.anova.rma"
                  return(sav)
               }
               btt <- .set.btt(btt, x$p, x$int.incl, colnames(x$X), fixed=fixed)
            }

            m <- length(btt)

            if (missing(rhs)) {
               rhs <- rep(0, m)
            } else {
               if (length(rhs) == 1L)
                  rhs <- rep(rhs, m)
               if (length(rhs) != m)
                  stop(mstyle$stop(paste0("Length of 'rhs' (", length(rhs), ") does not match the number of coefficients tested (", m, ").")))
            }

            x$b[btt,] <- x$beta[btt,] <- x$b[btt,] - rhs

            if (inherits(x, "robust.rma") && x$robumethod == "clubSandwich") {

               cs.wald <- try(clubSandwich::Wald_test(x, cluster=x$cluster, vcov=x$vb, test=x$wald_test, constraints=clubSandwich::constrain_zero(btt)), silent=!isTRUE(ddd$verbose))

               if (inherits(cs.wald, "try-error"))
                  stop(mstyle$stop("Could not obtain the cluster-robust Wald test (use verbose=TRUE for more details)."))

               QM   <- max(0, cs.wald$Fstat)
               QMdf <- c(cs.wald$df_num, cs.wald$df_denom)
               QMp  <- cs.wald$p_val

            } else {

               #QM <- try(as.vector(t((x$beta)[btt]-rhs) %*% chol2inv(chol(x$vb[btt,btt])) %*% (x$beta[btt]-rhs)), silent=TRUE)
               QM <- try(as.vector(t(x$beta)[btt] %*% chol2inv(chol(x$vb[btt,btt])) %*% x$beta[btt]), silent=TRUE)

               if (inherits(QM, "try-error"))
                  QM <- NA

               if (is.element(x$test, c("knha","adhoc","t"))) {
                  QM   <- QM / m
                  QMdf <- c(m, x$QMdf[2])
                  QMp  <- pf(QM, df1=QMdf[1], df2=QMdf[2], lower.tail=FALSE)
               } else {
                  QMdf <- c(m, NA)
                  QMp  <- pchisq(QM, df=QMdf[1], lower.tail=FALSE)
               }

            }

            res <- list(QM=QM, QMdf=QMdf, QMp=QMp, btt=btt, k=x$k, p=x$p, m=m, test=x$test, digits=digits, type="Wald.btt", class=class(x))

         }

      } else {

         if (inherits(object, "rma.ls") && !missing(Z)) {

            ### if 'Z' has been specified, then do Wald-type test(s) via 'Z' argument

            if (!missing(X))
               stop(mstyle$stop("Can only specify either 'X' or 'Z', but not both."))

            if (.is.vector(Z))
               Z <- rbind(Z)

            if (is.data.frame(Z))
               Z <- as.matrix(Z)

            if (is.character(Z))
               stop(mstyle$stop("Argument 'Z' must be a numeric vector/matrix."))

            ### if model has an intercept term and Z has q-1 columns, assume user left out the intercept and add it automatically

            if (x$Z.int.incl && ncol(Z) == (x$q-1))
               Z <- cbind(1, Z)

            if (ncol(Z) != x$q)
               stop(mstyle$stop(paste0("Length or number of columns of 'Z' (", ncol(Z), ") does not match the number of scale coefficients (", x$q, ").")))

            m <- nrow(Z)

            ### specification of the right-hand side

            if (missing(rhs)) {
               rhs <- rep(0, m)
            } else {
               if (length(rhs) == 1L)
                  rhs <- rep(rhs, m)
               if (length(rhs) != m)
                  stop(mstyle$stop(paste0("Length of 'rhs' (", length(rhs), ") does not match the number of linear combinations (", m, ").")))
            }

            ### test of individual hypotheses

            Za  <- Z %*% x$alpha - rhs
            vZa <- Z %*% x$va %*% t(Z)

            se <- sqrt(diag(vZa))
            zval <- c(Za/se)

            if (is.element(x$test, c("knha","adhoc","t"))) {
               pval <- if (x$ddf.alpha > 0) 2*pt(abs(zval), df=x$ddf.alpha, lower.tail=FALSE) else rep(NA,m)
            } else {
               pval <- 2*pnorm(abs(zval), lower.tail=FALSE)
            }

            ### omnibus test of all hypotheses (only possible if 'Z' is of full rank)

            QS  <- NA ### need this in case QS cannot be calculated below
            QSp <- NA ### need this in case QSp cannot be calculated below

            if (rankMatrix(Z) == m) {

               QS <- try(as.vector(t(Za) %*% chol2inv(chol(vZa)) %*% Za), silent=TRUE)

               if (inherits(QS, "try-error"))
                  QS <- NA

               if (is.element(x$test, c("knha","adhoc","t"))) {
                  QS   <- QS / m
                  QSdf <- c(m, x$QSdf[2])
                  QSp  <- if (QSdf[2] > 0) pf(QS, df1=QSdf[1], df2=QSdf[2], lower.tail=FALSE) else NA
               } else {
                  QSdf <- c(m, NA)
                  QSp  <- pchisq(QS, df=QSdf[1], lower.tail=FALSE)
               }

            }

            ### create a data frame with each row specifying the linear combination tested

            hyp <- rep("", m)
            for (j in seq_len(m)) {
               Zj <- round(Z[j,], digits[["est"]]) ### coefficients for the jth contrast
               sel <- Zj != 0 ### TRUE if coefficient is != 0
               hyp[j] <- paste(paste(Zj[sel], rownames(x$alpha)[sel], sep="*"), collapse=" + ") ### coefficient*variable + coefficient*variable ...
               hyp[j] <- gsub("1*", "", hyp[j], fixed=TRUE) ### turn '+1' into '+' and '-1' into '-'
               hyp[j] <- gsub("+ -", "- ", hyp[j], fixed=TRUE) ### turn '+ -' into '-'
            }
            if (identical(rhs, rep(0,m))) {
               hyp <- paste0(hyp, " = 0") ### add '= 0' at the right
            } else {
               if (length(unique(rhs)) == 1L) {
                  hyp <- paste0(hyp, " = ", round(rhs, digits=digits[["est"]])) ### add '= rhs' at the right
               } else {
                  hyp <- paste0(hyp, " = ", fmtx(rhs, digits=digits[["est"]])) ### add '= rhs' at the right
               }
            }
            hyp <- data.frame(hyp, stringsAsFactors=FALSE)
            colnames(hyp) <- ""
            rownames(hyp) <- paste0(seq_len(m), ":") ### add '1:', '2:', ... as row names

            res <- list(QS=QS, QSdf=QSdf, QSp=QSp, hyp=hyp, Za=Za, se=se, zval=zval, pval=pval, k=x$k, q=x$q, m=m, test=x$test, ddf=x$ddf.alpha, digits=digits, type="Wald.Za")

         } else {

            ### if 'X' has been specified, then do Wald-type test(s) via 'X' argument

            if (.is.vector(X))
               X <- rbind(X)

            if (is.data.frame(X))
               X <- as.matrix(X)

            if (is.character(X))
               stop(mstyle$stop("Argument 'X' must be a numeric vector/matrix."))

            ### if model has an intercept term and X has p-1 columns, assume user left out the intercept and add it automatically

            if (x$int.incl && ncol(X) == (x$p-1))
               X <- cbind(1, X)

            if (ncol(X) != x$p)
               stop(mstyle$stop(paste0("Length or number of columns of 'X' (", ncol(X), ") does not match the number of ", ifelse(inherits(object, "rma.ls"), "location", "model"), " coefficients (", x$p, ").")))

            m <- nrow(X)

            if (inherits(x, "robust.rma") && x$robumethod == "clubSandwich") {

               cs.lc <- try(clubSandwich::linear_contrast(x, cluster=x$cluster, vcov=x$vb, test=x$coef_test, contrasts=X, p_values=TRUE), silent=!isTRUE(ddd$verbose))

               if (inherits(cs.lc, "try-error"))
                  stop(mstyle$stop("Could not obtain the cluster-robust test(s) (use verbose=TRUE for more details)."))

               ddf <- cs.lc$df

               if (!missing(rhs))
                  warning(mstyle$warning("Cannot use 'rhs' argument for 'robust.rma' objects based on 'clubSandwich'."))

               rhs <- rep(0, m)

               Xb   <- cs.lc$Est
               se   <- cs.lc$SE
               zval <- c(Xb/se)
               pval <- cs.lc$p_val

               ### omnibus test of all hypotheses (only possible if 'X' is of full rank)

               QM   <- NA ### need this in case QMp cannot be calculated below
               QMdf <- NA ### need this in case X is not of full rank
               QMp  <- NA ### need this in case QMp cannot be calculated below

               if (rankMatrix(X) == m) {

                  cs.wald <- try(clubSandwich::Wald_test(x, cluster=x$cluster, vcov=x$vb, test=x$wald_test, constraints=X), silent=!isTRUE(ddd$verbose))

                  if (inherits(cs.wald, "try-error"))
                     stop(mstyle$stop("Could not obtain the cluster-robust omnibus Wald test (use verbose=TRUE for more details)."))

                  QM   <- max(0, cs.wald$Fstat)
                  QMdf <- c(cs.wald$df_num, cs.wald$df_denom)
                  QMp  <- cs.wald$p_val

               }

            } else {

               ### ddf calculation

               if (is.element(x$test, c("knha","adhoc","t"))) {
                  if (length(x$ddf) == 1L) {
                     ddf <- rep(x$ddf, m)
                  } else {
                     ddf <- rep(NA, m)
                     for (j in seq_len(m)) {
                        bn0 <- X[j,] != 0
                        ddf[j] <- min(x$ddf[bn0])
                     }
                  }
               } else {
                  ddf <- rep(NA, m)
               }

               ### specification of the right-hand side

               if (missing(rhs)) {
                  rhs <- rep(0, m)
               } else {
                  if (length(rhs) == 1L)
                     rhs <- rep(rhs, m)
                  if (length(rhs) != m)
                     stop(mstyle$stop(paste0("Length of 'rhs' (", length(rhs), ") does not match the number of linear combinations (", m, ").")))
               }

               ### test of individual hypotheses

               Xb  <- X %*% x$beta - rhs
               vXb <- X %*% x$vb %*% t(X)

               se <- sqrt(diag(vXb))
               zval <- c(Xb/se)

               if (is.element(x$test, c("knha","adhoc","t"))) {
                  pval <- sapply(seq_along(ddf), function(j) if (ddf[j] > 0) 2*pt(abs(zval[j]), df=ddf[j], lower.tail=FALSE) else NA)
               } else {
                  pval <- 2*pnorm(abs(zval), lower.tail=FALSE)
               }

               ### omnibus test of all hypotheses (only possible if 'X' is of full rank)

               QM   <- NA ### need this in case QMp cannot be calculated below
               QMdf <- NA ### need this in case X is not of full rank
               QMp  <- NA ### need this in case QMp cannot be calculated below

               if (rankMatrix(X) == m) {

                  ### use try(), since this could fail: this could happen when the var-cov matrix of the
                  ### fixed effects has been estimated using robust() -- 'vb' is then only guaranteed to
                  ### be positive semidefinite, so for certain linear combinations, vXb could be singular
                  ### (see Cameron & Miller, 2015, p. 326)

                  QM <- try(as.vector(t(Xb) %*% chol2inv(chol(vXb)) %*% Xb), silent=TRUE)

                  if (inherits(QM, "try-error"))
                     QM <- NA

                  if (is.element(x$test, c("knha","adhoc","t"))) {
                     QM   <- QM / m
                     QMdf <- c(m, min(ddf))
                     QMp  <- if (QMdf[2] > 0) pf(QM, df1=QMdf[1], df2=QMdf[2], lower.tail=FALSE) else NA
                  } else {
                     QMdf <- c(m, NA)
                     QMp  <- pchisq(QM, df=QMdf[1], lower.tail=FALSE)
                  }

               }

            }

            ### create a data frame with each row specifying the linear combination tested

            hyp <- rep("", m)
            for (j in seq_len(m)) {
               Xj <- round(X[j,], digits[["est"]]) ### coefficients for the jth contrast
               sel <- Xj != 0 ### TRUE if coefficient is != 0
               hyp[j] <- paste(paste(Xj[sel], rownames(x$beta)[sel], sep="*"), collapse=" + ") ### coefficient*variable + coefficient*variable ...
               hyp[j] <- gsub("1*", "", hyp[j], fixed=TRUE) ### turn '+1' into '+' and '-1' into '-'
               hyp[j] <- gsub("+ -", "- ", hyp[j], fixed=TRUE) ### turn '+ -' into '-'
            }
            if (identical(rhs, rep(0,m))) {
               hyp <- paste0(hyp, " = 0") ### add '= 0' at the right
            } else {
               if (length(unique(rhs)) == 1L) {
                  hyp <- paste0(hyp, " = ", round(rhs, digits=digits[["est"]])) ### add '= rhs' at the right
               } else {
                  hyp <- paste0(hyp, " = ", fmtx(rhs,  digits=digits[["est"]])) ### add '= rhs' at the right
               }
            }
            hyp <- data.frame(hyp, stringsAsFactors=FALSE)
            colnames(hyp) <- ""
            rownames(hyp) <- paste0(seq_len(m), ":") ### add '1:', '2:', ... as row names

            res <- list(QM=QM, QMdf=QMdf, QMp=QMp, hyp=hyp, Xb=Xb, se=se, zval=zval, pval=pval, k=x$k, p=x$p, m=m, test=x$test, ddf=ddf, digits=digits, type="Wald.Xb")

         }

      }

   } else {

      ### if 'object' and 'object2' have been specified, can use function to
      ### do model comparisons via a likelihood ratio test (and fit indices)

      if (!inherits(object2, "rma"))
         stop(mstyle$stop("Argument 'object2' must be an object of class \"rma\"."))

      if (inherits(object2, c("rma.mh","rma.peto")))
         stop(mstyle$stop("Function not applicable to objects of class \"rma.mh\" or \"rma.peto\"."))

      if (inherits(object2, "rma.glmm"))
         stop(mstyle$stop("Method not available for objects of class \"rma.glmm\"."))

      if (!identical(class(object), class(object2)))
         stop(mstyle$stop("Class of 'object' must be the same as class of 'object2'."))

      if (!is.null(ddd$test)) {
         test <- match.arg(ddd$test, c("LRT", "Wald"))
      } else {
         test <- "LRT"
      }

      ### assume 'object' is the full model and 'object2' the reduced model

      model.f <- object
      model.r <- object2

      ### number of parameters in the models

      parms.f <- model.f$parms
      parms.r <- model.r$parms

      ### check if they have the same number of parameters

      if (parms.f == parms.r)
         stop(mstyle$stop("Models have the same number of parameters. LRT not meaningful."))

      ### if parms.f < parms.r, then let 'object' be the reduced model and 'object2' the full model

      if (parms.f < parms.r) {
         model.f <- object2
         model.r <- object
         parms.f <- model.f$parms
         parms.r <- model.r$parms
      }

      ### check if models are based on the same data (TODO: also check for same weights?)
      ### note: using as.vector() to strip attributes/names, as.matrix() to make both V matrices non-sparse, and
      ###       isTRUE(all.equal()) because conversion to non-sparse can introduce some negligible discrepancies

      if (inherits(object, "rma.uni")) {
         if (!(identical(as.vector(model.f$yi), as.vector(model.r$yi)) && isTRUE(all.equal(as.vector(model.f$vi), as.vector(model.r$vi)))))
            stop(mstyle$stop("Observed outcomes and/or sampling variances not equal in the full and reduced model."))
      }

      if (inherits(object, "rma.mv")) {
         if (!(identical(as.vector(model.f$yi), as.vector(model.r$yi)) && isTRUE(all.equal(as.matrix(model.f$V), as.matrix(model.r$V)))))
            stop(mstyle$stop("Observed outcomes and/or sampling variances/covariances not equal in the full and reduced model."))
      }

      ### for Wald-type test, both models should be fitted using the same method

      if (test == "Wald" && (model.f$method != model.r$method))
         stop(mstyle$stop("Full and reduced model must use the same 'method' for the model fitting."))

      ### for LRTs, reduced model may use method="FE/EE/CE" and full model method="(RE)ML" but the other way around doesn't really make sense

      if (is.element(model.f$method, c("FE","EE","CE")) && !is.element(model.r$method, c("FE","EE","CE")))
         stop(mstyle$stop("Full model uses a fixed- and reduced model uses a random/mixed-effects model."))

      ### but have to check for a ML/REML mismatch

      if ((model.f$method == "ML" && model.r$method == "REML") || model.r$method == "ML" && model.f$method == "REML")
         stop(mstyle$stop(paste0("Mismatch between the use of ", model.f$method, " and ", model.r$method, " estimation in the full versus reduced model.")))

      ### for LRTs, using anything besides ML/REML is strictly speaking incorrect

      if (test == "LRT" && (!is.element(model.f$method, c("FE","EE","CE","ML","REML")) || !is.element(model.r$method, c("FE","EE","CE","ML","REML"))))
         warning(mstyle$warning("LRTs should be based on ML/REML estimation."))

      ### for LRTs based on REML estimation, check if fixed effects differ
      if (test == "LRT" && model.f$method == "REML" && (!identical(model.f$X, model.r$X))) {
         if (refit) {
            #message(mstyle$message("Refitting models with ML (instead of REML) estimation ..."))
            if (inherits(model.f, "rma.uni") && model.f$model == "rma.uni") {
               #model.f <- try(update(model.f, method="ML", data=model.f$data), silent=TRUE)
               args <- list(yi=model.f$yi, vi=model.f$vi, weights=model.f$weights, mods=model.f$X, intercept=FALSE, method="ML", weighted=model.f$weighted, test=model.f$test, level=model.f$level, tau2=ifelse(model.f$tau2.fix, model.f$tau2, NA), control=model.f$control, skipr2=TRUE)
               model.f <- try(suppressWarnings(.do.call(rma.uni, args)), silent=TRUE)
            } else {
               # note: this fails when building the docs with pkgdown; not sure why; the approach above at least works for 'rma.uni' objects and is more efficient as it skips the R^2 calculation
               model.f <- try(update(model.f, method="ML"), silent=TRUE)
            }
            if (inherits(model.f, "try-error"))
               stop(mstyle$stop("Refitting the full model with ML estimation failed."))
            if (inherits(model.r, "rma.uni") && model.r$model == "rma.uni") {
               #model.r <- try(update(model.r, method="ML", data=model.r$data), silent=TRUE)
               args <- list(yi=model.r$yi, vi=model.r$vi, weights=model.r$weights, mods=model.r$X, intercept=FALSE, method="ML", weighted=model.r$weighted, test=model.r$test, level=model.r$level, tau2=ifelse(model.r$tau2.fix, model.r$tau2, NA), control=model.r$control, skipr2=TRUE)
               model.r <- try(suppressWarnings(.do.call(rma.uni, args)), silent=TRUE)
            } else {
               model.r <- try(update(model.r, method="ML"), silent=TRUE)
            }
            if (inherits(model.r, "try-error"))
               stop(mstyle$stop("Refitting the reduced model with ML estimation failed."))
            parms.f <- model.f$parms
            parms.r <- model.r$parms
         } else {
            warning(mstyle$warning("REML comparisons not meaningful for models with different fixed effects\n(use 'refit=TRUE' to refit both models based on ML estimation)."), call.=FALSE)
         }
      }

      ### in this case, one could consider just taking the ML deviances, but this
      ### is really ad-hoc; there is some theory in Welham & Thompson (1997) about
      ### LRTs for fixed effects when using REML estimation, but this involves
      ### additional work

      ### could do even more checks for cases where the models are clearly not nested

      ######################################################################

      ### for 'rma.uni' objects, calculate pseudo R^2 value (based on the
      ### proportional reduction in tau^2) comparing full vs. reduced model

      if (inherits(object, "rma.uni") && !inherits(object, "rma.ls") && !inherits(object, "rma.gen")) {
         if (is.element(model.f$method, c("FE","EE","CE"))) {
            if (model.f$weighted) {
               if (is.null(model.f$weights)) {
                  lm.f <- lm(model.f$yi ~ model.f$X, weights=1/model.f$vi)
               } else {
                  lm.f <- lm(model.f$yi ~ model.f$X, weights=model.f$weights)
               }
            } else {
               lm.f <- lm(model.f$yi ~ model.f$X)
            }
            if (model.r$weighted) {
               if (is.null(model.r$weights)) {
                  lm.r <- lm(model.r$yi ~ model.r$X, weights=1/model.r$vi)
               } else {
                  lm.r <- lm(model.r$yi ~ model.r$X, weights=model.r$weights)
               }
            } else {
               lm.r <- lm(model.r$yi ~ model.r$X)
            }
            s2.f <- sigma(lm.f)^2
            s2.r <- sigma(lm.r)^2
            R2 <- 100 * max(0, (s2.r - s2.f) / s2.r)
         } else if (identical(model.r$tau2,0)) {
            R2 <- 0
         } else {
            R2 <- 100 * max(0, (model.r$tau2 - model.f$tau2) / model.r$tau2)
         }
      } else {
         R2 <- NA
      }

      ### for 'rma.uni' objects, extract tau^2 estimates

      if (inherits(object, "rma.uni") && !inherits(object, "rma.ls") && !inherits(object, "rma.gen")) {
         tau2.f <- model.f$tau2
         tau2.r <- model.r$tau2
      } else {
         tau2.f <- NA
         tau2.r <- NA
      }

      if (test == "LRT") {

         parms.diff <- parms.f - parms.r

         if (model.f$method == "REML") {

            LRT <- model.r$fit.stats["dev","REML"] - model.f$fit.stats["dev","REML"]
            fit.stats.f <- t(model.f$fit.stats)["REML",] # to keep (row)names of fit.stats
            fit.stats.r <- t(model.r$fit.stats)["REML",] # to keep (row)names of fit.stats

         } else {

            LRT <- model.r$fit.stats["dev","ML"] - model.f$fit.stats["dev","ML"]
            fit.stats.f <- t(model.f$fit.stats)["ML",]
            fit.stats.r <- t(model.r$fit.stats)["ML",]

         }

         ### set LRT to 0 if LRT < 0 (this should not happen, but could due to numerical issues)

         LRT[LRT < 0] <- 0

         pval <- pchisq(LRT, df=parms.diff, lower.tail=FALSE)

         res <- list(fit.stats.f=fit.stats.f, fit.stats.r=fit.stats.r, parms.f=parms.f, parms.r=parms.r,
                     LRT=LRT, pval=pval, QE.f=model.f$QE, QE.r=model.r$QE,
                     tau2.f=tau2.f, tau2.r=tau2.r, R2=R2,
                     method=model.f$method, class.f=class(model.f), digits=digits, type="LRT")

      }

      if (test == "Wald") {

         btt <- setdiff(colnames(model.f$X), colnames(model.r$X))

         if (length(btt) == 0L)
            stop(mstyle$stop("Full and reduced models appear to contain the same moderators."))

         if (length(setdiff(colnames(model.r$X), colnames(model.f$X))) != 0L)
            stop(mstyle$stop("There are coefficients in the reduced model that are not in the full model."))

         btt <- charmatch(btt, colnames(model.f$X))

         if (anyNA(btt))
            stop(mstyle$stop("Cannot identify coefficients to test."))

         res <- anova(model.f, btt=btt)
         return(res)

      }

   }

   class(res) <- "anova.rma"
   return(res)

}
