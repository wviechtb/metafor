predict.rma <- function(object, newmods, intercept, tau2.levels, gamma2.levels, hetvar,
                        addx=FALSE, level, adjust=FALSE, digits, transf, targs, vcov=FALSE, ...) {

   #########################################################################

   mstyle <- .get.mstyle()

   .chkclass(class(object), must="rma", notav="rma.ls")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   x <- object

   mf <- match.call()

   # so that pairmat() works when the model object is not specified

   if (any(grepl("pairmat(", as.character(mf), fixed=TRUE))) {
      try(assign("pairmat", object, envir=.metafor), silent=TRUE)
      on.exit(suppressWarnings(rm("pairmat", envir=.metafor)))
   }

   if (missing(newmods))
      newmods <- NULL

   if (missing(intercept)) {
      intercept <- x$intercept
      int.spec <- FALSE
   } else {
      int.spec <- TRUE
   }

   if (missing(tau2.levels))
      tau2.levels <- NULL

   if (missing(gamma2.levels))
      gamma2.levels <- NULL

   if (missing(hetvar)) {
      hetvar <- NULL
   } else {
      if (inherits(object, "rma.mv")) {
         if (!is.null(tau2.levels) || !is.null(gamma2.levels)) {
            tau2.levels <- NULL
            gamma2.levels <- NULL
            warning(mstyle$warning("Arguments 'tau2.levels' and 'gamma2.levels' ignored when specifying the 'hetvar' argument."), call.=FALSE)
         }
      }
      if (!is.numeric(hetvar))
         stop(mstyle$stop("Argument 'hetvar' must be a numeric vector."))
   }

   if (missing(level))
      level <- x$level

   if (missing(digits)) {
      digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
   }

   if (missing(transf))
      transf <- FALSE

   if (missing(targs))
      targs <- NULL

   funlist <- lapply(list(transf.exp.int, transf.ilogit.int, transf.iprobit.int, transf.ztor.int, transf.iarcsin.int, transf.iahw.int, transf.iabt.int, transf.dtocles.int, transf.exp.mode, transf.ilogit.mode, transf.iprobit.mode, transf.ztor.mode, transf.iarcsin.mode, transf.iahw.mode, transf.iabt.mode), deparse)

   if (is.null(targs) && any(sapply(funlist, identical, deparse(transf))) && inherits(x, c("rma.uni","rma.glmm")) && length(x$tau2 == 1L))
      targs <- list(tau2=x$tau2)

   level <- .level(level)

   if (!is.logical(adjust))
      stop(mstyle$stop("Argument 'adjust' must be a logical."))

   ddd <- list(...)

   .chkdots(ddd, c("pi.type", "predtype", "newvi", "verbose"))

   pi.type  <- .chkddd(ddd$pi.type, "default", tolower(ddd$pi.type))
   predtype <- .chkddd(ddd$predtype, pi.type, tolower(ddd$predtype))

   if (x$int.only && !is.null(newmods))
      stop(mstyle$stop("Cannot specify new moderator values for models without moderators."))

   rnames <- NULL

   #########################################################################

   # TODO: can this be simplified? (every time I sit down and stare at the mess below, it gives me a headache)

   if (is.null(newmods)) {

      # if no new moderator values are specified

      if (!inherits(object, "rma.mv") || (inherits(object, "rma.mv") && any(is.element(object$struct, c("GEN","GDIAG"))))) {

         # for rma.uni, rma.mh, rma.peto, and rma.glmm objects

         if (x$int.only) {                                                  # if intercept-only model predict only the intercept
            k.new <- 1L                                                     #
            X.new <- cbind(1)                                               #
         } else {                                                           # otherwise predict for all k.f studies (including studies with NAs)
            k.new <- x$k.f                                                  #
            X.new <- x$X.f                                                  #
         }                                                                  #

      } else {

         # for rma.mv objects

         if (x$int.only) {                                                       # if intercept-only model:
            if (!x$withG) {                                                      #   # if there is no G structure (and hence also no H structure)
               k.new <- 1L                                                       #   # then we just need to predict the intercept once
               X.new <- cbind(1)                                                 #
            }                                                                    #
            if (x$withG && x$withH) {                                            #   # if there is both a G and H structure
               if (is.null(tau2.levels) && is.null(gamma2.levels)) {             #      # and user has not specified tau2s.levels and gamma2.levels
                  k.new <- x$tau2s * x$gamma2s                                   #         # then we need to predict intercepts for all combinations of tau2 and gamma2 values
                  X.new <- cbind(rep(1,k.new))                                   #
                  if (x$tau2s == 1) {                                            #         # if there is only a single tau^2
                     tau2.levels <- rep(1,k.new)                                 #         # then tau2.levels should be 1 repeated k.new times
                  } else {                                                       #
                     tau2.levels <- rep(levels(x$mf.g.f$inner), each=x$gamma2s)  #       # otherwise repeat actual levels gamma2s times
                  }                                                              #
                  if (x$gamma2s == 1) {                                          #         # if there is only a single gamma^2 value
                     gamma2.levels <- rep(1,k.new)                               #         # then gamma2.levels should be 1 repeated k.new times
                  } else {                                                       #
                     gamma2.levels <- rep(levels(x$mf.h.f$inner), times=x$tau2s) #       # otherwise repeat actual levels tau2s times
                  }                                                              #
               }                                                                 #
               if ((!is.null(tau2.levels) && is.null(gamma2.levels)) ||          #   # if user specified only one of tau2.levels and gamma2.levels, throw an error
                   (is.null(tau2.levels) && !is.null(gamma2.levels)))            #
                  stop(mstyle$stop("Either specify both of 'tau2.levels' and 'gamma2.levels' or neither."))
               if (!is.null(tau2.levels) && !is.null(gamma2.levels)) {           #   # if user has specified both tau2s.levels and gamma2.levels
                  if (length(tau2.levels) != length(gamma2.levels))              #
                     stop(mstyle$stop("Length of 'tau2.levels' and 'gamma2.levels' are not the same."))
                  k.new <- length(tau2.levels)                                   #      # then we need to predict intercepts for those level combinations
                  X.new <- cbind(rep(1,k.new))                                   #
               }                                                                 #
            }                                                                    #
            if (x$withG && !x$withH) {                                           #   # if there is only a G structure (and no H structure)
               if (is.null(tau2.levels)) {                                       #      # and user has not specified tau2.levels
                  k.new <- x$tau2s                                               #         # then we need to predict intercepts for all tau2 values
                  X.new <- cbind(rep(1,k.new))                                   #
                  if (x$tau2s == 1) {                                            #
                     tau2.levels <- rep(1, k.new)                                #
                  } else {                                                       #
                     tau2.levels <- levels(x$mf.g.f$inner)                       #
                  }                                                              #
               } else {                                                          #      # and the user has specified tau2.levels
                  k.new <- length(tau2.levels)                                   #         # then we need to predict intercepts for those levels
                  X.new <- cbind(rep(1,k.new))                                   #
               }                                                                 #
               gamma2.levels <- rep(1, k.new)                                    #
            }                                                                    #
         } else {                                                                # if not an intercept-only model
            k.new <- x$k.f                                                       #   # then predict for all k.f studies (including studies with NAs)
            X.new <- x$X.f                                                       #
            if (!is.null(tau2.levels) || !is.null(gamma2.levels))                #
               warning(mstyle$warning("Arguments 'tau2.levels' and 'gamma2.levels' ignored when obtaining fitted values."), call.=FALSE)
            tau2.levels <- as.character(x$mf.g.f$inner)                          #
            gamma2.levels <- as.character(x$mf.h.f$inner)                        #
         }                                                                       #

      }

   } else {

      # if new moderator values have been specified

      if (!(.is.vector(newmods) || inherits(newmods, "matrix")))
         stop(mstyle$stop(paste0("Argument 'newmods' should be a vector or matrix, but is of class '", class(newmods)[1], "'.")))

      singlemod <- (NCOL(newmods) == 1L) && ((!x$int.incl && x$p == 1L) || (x$int.incl && x$p == 2L))

      if (singlemod) {                                                           # if single moderator (multiple k.new possible) (either without or with intercept in the model)
         k.new <- length(newmods)                                                # (but when specifying a matrix, it must be a column vector for this work)
         X.new <- cbind(c(newmods))                                              #
         if (.is.vector(newmods)) {                                              #
            rnames <- names(newmods)                                             #
         } else {                                                                #
            rnames <- rownames(newmods)                                          #
         }                                                                       #
      } else {                                                                   # in case the model has more than one predictor:
         if (.is.vector(newmods) || nrow(newmods) == 1L) {                       #   # if user gives one vector or one row matrix (only one k.new):
            k.new <- 1L                                                          #
            X.new <- rbind(newmods)                                              #
            if (inherits(newmods, "matrix"))                                     #
               rnames <- rownames(newmods)                                       #
         } else {                                                                #   # if user gives multiple rows and columns (multiple k.new):
            k.new <- nrow(newmods)                                               #
            X.new <- cbind(newmods)                                              #
            rnames <- rownames(newmods)                                          #
         }                                                                       #
         # allow matching of terms by names (note: only possible if all columns in X.new and x$X have colnames)
         if (!is.null(colnames(X.new)) && all(colnames(X.new) != "") && !is.null(colnames(x$X)) && all(colnames(x$X) != "")) {
            colnames.mod <- colnames(x$X)
            if (x$int.incl)
               colnames.mod <- colnames.mod[-1]
            pos <- sapply(colnames(X.new), function(colname) {
                     d <- c(adist(colname, colnames.mod, costs=c(ins=1, sub=Inf, del=Inf))) # compute edit distances with Inf costs for substitutions/deletions
                     if (all(is.infinite(d))) # if there is no match, then all elements are Inf
                        stop(mstyle$stop(paste0("Could not find variable '", colname, "' in the model.")))
                     d <- which(d == min(d)) # don't use which.min() since that only finds the first minimum
                     if (length(d) > 1L) # if there is no unique match, then there is more than one minimum
                        stop(mstyle$stop(paste0("Could not match up variable '", colname, "' uniquely to a variable in the model.")))
                     return(d)
                     })
            if (anyDuplicated(pos)) { # if the same name is used more than once, then there will be duplicated pos values
               dups <- paste(unique(colnames(X.new)[duplicated(pos)]), collapse=", ")
               stop(mstyle$stop(paste0("Found multiple matches for the same variable name (", dups, ").")))
            }
            if (length(pos) != length(colnames.mod)) {
               no.match <- colnames.mod[seq_along(colnames.mod)[-pos]]
               if (length(no.match) > 3L)
                  stop(mstyle$stop(paste0("Argument 'newmods' does not specify values for these variables: ", paste0(no.match[1:3], collapse=", "), ", ...")))
               if (length(no.match) > 1L)
                  stop(mstyle$stop(paste0("Argument 'newmods' does not specify values for these variables: ", paste0(no.match, collapse=", "))))
               if (length(no.match) == 1L)
                  stop(mstyle$stop(paste0("Argument 'newmods' does not specify values for this variable: ", no.match)))
            }
            X.new <- X.new[,order(pos),drop=FALSE]
            colnames(X.new) <- colnames.mod
         }
      }

      if (inherits(X.new[1,1], "character"))
         stop(mstyle$stop("Argument 'newmods' should only contain numeric variables."))

      # if the user has specified newmods and an intercept was included in the original model, add the intercept to X.new
      # but user can also decide to remove the intercept from the predictions with intercept=FALSE (but only do this when
      # newmods was not a matrix with p columns)

      if (!singlemod && ncol(X.new) == x$p) {

         if (int.spec)
            warning(mstyle$warning("Arguments 'intercept' ignored when 'newmods' includes 'p' columns."), call.=FALSE)

      } else {

         if (x$int.incl) {

            if (intercept) {
               X.new <- cbind(intrcpt=1, X.new)
            } else {
               X.new <- cbind(intrcpt=0, X.new)
            }

         }

      }

      if (ncol(X.new) != x$p)
         stop(mstyle$stop(paste0("Dimensions of 'newmods' (", ncol(X.new), ") do not the match dimensions of the model (", x$p, ").")))

   }

   if (is.null(X.new))
      stop(mstyle$stop("Matrix 'X.new' is NULL."))

   #return(list(k.new=k.new, tau2=x$tau2, gamma2=x$gamma2, tau2.levels=tau2.levels, gamma2.levels=gamma2.levels))

   #########################################################################

   # for rma.mv models with multiple tau^2 values, must use tau2.levels argument when using newmods to obtain prediction intervals

   if (inherits(object, "rma.mv") && x$withG) {

      if (x$tau2s > 1L) {

         if (is.null(tau2.levels)) {

            #warning(mstyle$warning("Must specify the 'tau2.levels' argument to obtain prediction intervals."), call.=FALSE)

         } else {

            # if tau2.levels argument is a character vector, check that specified tau^2 values actually exist
            if (!is.numeric(tau2.levels) && anyNA(pmatch(tau2.levels, x$g.levels.f[[1]], duplicates.ok=TRUE)))
               stop(mstyle$stop("Non-existing levels specified via 'tau2.levels' argument."))

            # if tau2.levels argument is numeric, check that specified tau^2 values actually exist
            if (is.numeric(tau2.levels)) {
               tau2.levels <- round(tau2.levels)
               if (any(tau2.levels < 1) || any(tau2.levels > x$g.nlevels.f[1]))
                  stop(mstyle$stop("Non-existing tau^2 values specified via 'tau2.levels' argument."))
            }

            # allow quick setting of all levels
            tau2.levels <- .expand1(tau2.levels, k.new)

            # check length of tau2.levels argument
            if (length(tau2.levels) != k.new)
               stop(mstyle$stop(paste0("Length of the 'tau2.levels' argument (", length(tau2.levels), ") does not match the number of predicted values (", k.new, ").")))

         }

      } else {

         tau2.levels <- rep(1, k.new)

      }

   }

   # for rma.mv models with multiple gamma^2 values, must use gamma.levels argument when using newmods to obtain prediction intervals

   if (inherits(object, "rma.mv") && x$withH) {

      if (x$gamma2s > 1L) {

         if (is.null(gamma2.levels)) {

            #warning(mstyle$warning("Must specify the 'gamma2.levels' argument to obtain prediction intervals."), call.=FALSE)

         } else {

            # if gamma2.levels argument is a character vector, check that specified gamma^2 values actually exist
            if (!is.numeric(gamma2.levels) && anyNA(pmatch(gamma2.levels, x$h.levels.f[[1]], duplicates.ok=TRUE)))
               stop(mstyle$stop("Non-existing levels specified via 'gamma2.levels' argument."))

            # if gamma2.levels argument is numeric, check that specified gamma^2 values actually exist
            if (is.numeric(gamma2.levels)) {
               gamma2.levels <- round(gamma2.levels)
               if (any(gamma2.levels < 1) || any(gamma2.levels > x$h.nlevels.f[1]))
                  stop(mstyle$stop("Non-existing gamma^2 values specified via 'gamma2.levels' argument."))
            }

            # allow quick setting of all levels
            gamma2.levels <- .expand1(gamma2.levels, k.new)

            # check length of gamma2.levels argument
            if (length(gamma2.levels) != k.new)
               stop(mstyle$stop(paste0("Length of the 'gamma2.levels' argument (", length(gamma2.levels), ") does not match the number of predicted values (", k.new, ").")))

         }

      } else {

         gamma2.levels <- rep(1, k.new)

      }

   }

   #########################################################################

   if (inherits(x, "robust.rma") && x$robumethod == "clubSandwich") {

      if (x$coef_test == "saddlepoint")
         stop(mstyle$stop("Cannot use method when saddlepoint correction was used."))

      cs.lc <- try(clubSandwich::linear_contrast(x, cluster=x$cluster, vcov=x$vb, test=x$coef_test, contrasts=X.new, p_values=FALSE, level=1-level), silent=!isTRUE(ddd$verbose))

      if (inherits(cs.lc, "try-error"))
         stop(mstyle$stop("Could not obtain the linear contrast(s) (use verbose=TRUE for more details)."))

      pred  <- cs.lc$Est
      se    <- cs.lc$SE
      vpred <- se^2
      ddf   <- cs.lc$df

      crit <- sapply(seq_along(ddf), function(j) if (ddf[j] > 0) qt(level/ifelse(adjust, 2*k.new, 2), df=ddf[j], lower.tail=FALSE) else NA_real_)

      #ci.lb <- cs.lc$CI_L
      #ci.ub <- cs.lc$CI_U
      ci.lb <- pred - crit * se
      ci.ub <- pred + crit * se

      x$test <- switch(x$coef_test, "z"="z", "naive-t"="t", "naive-tp"="t", "Satterthwaite"="t")

   } else {

      # ddf calculation for x$test %in% c("knha","adhoc","t") but also need this
      # for pi.ddf calculation when test="z" and predtype %in% c("riley","t")

      if (length(x$ddf) == 1L) {
         ddf <- rep(x$ddf, k.new)     # when test="z", x$ddf is NA, so this then results in a vector of NAs
      } else {
         ddf <- rep(NA_integer_, k.new)
         for (j in seq_len(k.new)) {
            bn0 <- X.new[j,] != 0     # determine which coefficients are involved in the linear contrast
            ddf[j] <- min(x$ddf[bn0]) # take the smallest ddf value for those coefficients
         }
      }
      ddf[is.na(ddf)] <- x$k - x$p    # when test="z", turn all NAs into the usual k-p dfs

      # predicted values, SEs, and confidence intervals

      pred  <- rep(NA_real_, k.new)
      vpred <- rep(NA_real_, k.new)

      for (i in seq_len(k.new)) {
         Xi.new   <- X.new[i,,drop=FALSE]
         pred[i]  <- Xi.new %*% x$beta
         vpred[i] <- Xi.new %*% tcrossprod(x$vb, Xi.new)
      }

      if (is.element(x$test, c("knha","adhoc","t"))) {
         crit <- sapply(seq_along(ddf), function(j) if (ddf[j] > 0) qt(level/ifelse(adjust, 2*k.new, 2), df=ddf[j], lower.tail=FALSE) else NA_real_)
      } else {
         crit <- qnorm(level/ifelse(adjust, 2*k.new, 2), lower.tail=FALSE)
      }

      vpred[vpred < 0] <- NA_real_
      se <- sqrt(vpred)
      ci.lb <- pred - crit * se
      ci.ub <- pred + crit * se

   }

   #########################################################################

   if (vcov)
      vcovpred <- symmpart(X.new %*% x$vb %*% t(X.new))

   if (predtype == "simple") {
      crit <- qnorm(level/ifelse(adjust, 2*k.new, 2), lower.tail=FALSE)
      vpred <- 0
   }

   pi.ddf <- ddf

   if (is.element(predtype, c("riley","t"))) {
      if (predtype == "riley")
         pi.ddf <- ddf - x$parms + x$p
      if (predtype == "t")
         pi.ddf <- ddf
      pi.ddf[pi.ddf < 1] <- 1
      crit <- sapply(seq_along(pi.ddf), function(j) if (pi.ddf[j] > 0) qt(level/ifelse(adjust, 2*k.new, 2), df=pi.ddf[j], lower.tail=FALSE) else NA_real_)
   }

   if (is.null(ddd$newvi)) {
      newvi <- 0
   } else {
      newvi <- ddd$newvi
      newvi <- .expand1(newvi, k.new)
      if (length(newvi) != k.new)
         stop(mstyle$stop(paste0("Length of the 'newvi' argument (", length(newvi), ") does not match the number of predicted values (", k.new, ").")))
   }

   #########################################################################

   # prediction intervals

   pi.se <- NULL

   if (is.null(hetvar)) {

      if (!inherits(object, "rma.mv")) {

         # for rma.uni, rma.mh, rma.peto, and rma.glmm objects (in rma.mh and rma.peto, tau2 = 0 by default and stored as such)

         pi.se <- sqrt(vpred + x$tau2 + newvi)
         pi.lb <- pred - crit * pi.se
         pi.ub <- pred + crit * pi.se

      } else {

         # for rma.mv objects

         if (!x$withG) {

            # if there is no G structure (and hence no H structure), there are no tau2 and gamma2 values, so just add the sum of all of the sigma2 values

            pi.se <- sqrt(vpred + sum(x$sigma2) + newvi)
            pi.lb <- pred - crit * pi.se
            pi.ub <- pred + crit * pi.se

         }

         if (x$withG && !x$withH) {

            # if there is a G structure but no H structure

            if (x$tau2s == 1L) {

               # if there is only a single tau^2 value, always add that (in addition to the sum of all of the sigma^2 values)

               pi.se <- sqrt(vpred + sum(x$sigma2) + x$tau2 + newvi)
               pi.lb <- pred - crit * pi.se
               pi.ub <- pred + crit * pi.se

            } else {

               if (is.null(tau2.levels)) {

                  # if user has not specified tau2.levels, cannot compute bounds

                  pi.lb <- rep(NA_real_, k.new)
                  pi.ub <- rep(NA_real_, k.new)
                  tau2.levels <- rep(NA, k.new)

               } else {

                  # if there are multiple tau^2 values, either let user define numerically which value(s) to use or
                  # match the position of the specified tau2.levels to the levels of the inner factor in the model

                  if (!is.numeric(tau2.levels))
                     tau2.levels <- pmatch(tau2.levels, x$g.levels.f[[1]], duplicates.ok=TRUE)

                  pi.se <- sqrt(vpred + sum(x$sigma2) + x$tau2[tau2.levels] + newvi)
                  pi.lb <- pred - crit * pi.se
                  pi.ub <- pred + crit * pi.se
                  tau2.levels <- x$g.levels.f[[1]][tau2.levels]

               }

            }

         }

         if (x$withG && x$withH) {

            # if there is a G structure and an H structure

            if (x$tau2s == 1L && x$gamma2s == 1L) {

               # if there is only a single tau^2 and gamma^2 value, always add that (in addition to the sum of all of the sigma^2 values)

               pi.se <- sqrt(vpred + sum(x$sigma2) + x$tau2 + x$gamma2 + newvi)
               pi.lb <- pred - crit * pi.se
               pi.ub <- pred + crit * pi.se

            } else {

               if (is.null(tau2.levels) || is.null(gamma2.levels)) {

                  # if user has not specified tau2.levels and gamma2.levels, cannot compute bounds

                  pi.lb <- rep(NA_real_, k.new)
                  pi.ub <- rep(NA_real_, k.new)
                  tau2.levels <- rep(NA, k.new)
                  gamma2.levels <- rep(NA, k.new)

               } else {

                  # if there are multiple tau^2 and/or gamma^2 values, either let user define numerically which value(s) to use or
                  # match the position of the specified tau2.levels and gamma2.levels to the levels of the inner factors in the model

                  if (!is.numeric(tau2.levels))
                     tau2.levels <- pmatch(tau2.levels, x$g.levels.f[[1]], duplicates.ok=TRUE)
                  if (!is.numeric(gamma2.levels))
                     gamma2.levels <- pmatch(gamma2.levels, x$h.levels.f[[1]], duplicates.ok=TRUE)

                  pi.se <- sqrt(vpred + sum(x$sigma2) + x$tau2[tau2.levels] + x$gamma2[gamma2.levels] + newvi)
                  pi.lb <- pred - crit * pi.se
                  pi.ub <- pred + crit * pi.se
                  tau2.levels <- x$g.levels.f[[1]][tau2.levels]
                  gamma2.levels <- x$h.levels.f[[1]][gamma2.levels]

               }

            }

         }

      }

   } else {

      hetvar <- .expand1(hetvar, k.new)
      if (length(hetvar) != k.new)
         stop(mstyle$stop(paste0("Length of the 'hetvar' argument (", length(newvi), ") does not match the number of predicted values (", k.new, ").")))

      pi.se <- sqrt(vpred + hetvar + newvi)
      pi.lb <- pred - crit * pi.se
      pi.ub <- pred + crit * pi.se

   }

   #########################################################################

   # apply transformation function if one has been specified

   if (is.function(transf)) {
      if (is.null(targs)) {
         pred  <- sapply(pred, transf)
         se    <- rep(NA_real_, k.new)
         ci.lb <- sapply(ci.lb, transf)
         ci.ub <- sapply(ci.ub, transf)
         pi.lb <- sapply(pi.lb, transf)
         pi.ub <- sapply(pi.ub, transf)
      } else {
         if (!is.primitive(transf) && !is.null(targs) && length(formals(transf)) == 1L)
            stop(mstyle$stop("Function specified via 'transf' does not appear to have an argument for 'targs'."))
         pred  <- sapply(pred, transf, targs)
         se    <- rep(NA_real_, k.new)
         ci.lb <- sapply(ci.lb, transf, targs)
         ci.ub <- sapply(ci.ub, transf, targs)
         pi.lb <- sapply(pi.lb, transf, targs)
         pi.ub <- sapply(pi.ub, transf, targs)
      }
      do.transf <- TRUE
   } else {
      do.transf <- FALSE
   }

   # make sure order of intervals is always increasing

   tmp <- .psort(ci.lb, ci.ub)
   ci.lb <- tmp[,1]
   ci.ub <- tmp[,2]

   tmp <- .psort(pi.lb, pi.ub)
   pi.lb <- tmp[,1]
   pi.ub <- tmp[,2]

   # use study labels from the object when the model has moderators and no new moderators have been specified
   # otherwise, just use consecutive numbers to label the predicted values

   if (is.null(newmods) && !x$int.only) {
      slab <- x$slab
   } else {
      slab <- seq_len(k.new)
      if (!is.null(rnames))
         slab <- rnames
   }

   # add row/colnames to vcovpred

   if (vcov)
      rownames(vcovpred) <- colnames(vcovpred) <- slab

   # but when predicting just a single value, use "" as study label

   if (k.new == 1L && is.null(rnames))
      slab <- ""

   # handle NAs

   not.na <- rep(TRUE, k.new)

   if (na.act == "na.omit") {
      if (is.null(newmods) && !x$int.only) {
         not.na <- x$not.na
      } else {
         not.na <- !is.na(pred)
      }
   }

   #if (na.act == "na.omit") {
   #   not.na <- !is.na(pred)
   #} else {
   #   not.na <- rep(TRUE, k.new)
   #}

   if (na.act == "na.fail" && any(!x$not.na))
      stop(mstyle$stop("Missing values in results."))

   out <- list(pred=pred[not.na], se=se[not.na], ci.lb=ci.lb[not.na], ci.ub=ci.ub[not.na], pi.lb=pi.lb[not.na], pi.ub=pi.ub[not.na], cr.lb=pi.lb[not.na], cr.ub=pi.ub[not.na])

   if (vcov)
      vcovpred <- vcovpred[not.na,not.na,drop=FALSE]

   if (na.act == "na.exclude" && is.null(newmods) && !x$int.only) {

      out <- lapply(out, function(val) ifelse(x$not.na, val, NA_real_))

      if (vcov) {
         vcovpred[!x$not.na,] <- NA_real_
         vcovpred[,!x$not.na] <- NA_real_
      }

   }

   # add tau2.levels values to list

   if (inherits(object, "rma.mv") && x$withG && x$tau2s > 1L)
      out$tau2.level <- tau2.levels

   # add gamma2.levels values to list

   if (inherits(object, "rma.mv") && x$withH && x$gamma2s > 1L)
      out$gamma2.level <- gamma2.levels

   # add X matrix to list

   if (addx) {
      out$X <- matrix(X.new[not.na,], ncol=x$p)
      colnames(out$X) <- colnames(x$X)
   }

   # add slab values to list

   out$slab <- slab[not.na]

   # add some additional info

   out$digits <- digits
   out$method <- x$method
   out$transf <- do.transf
   out$pred.type <- "location"

   if (x$test != "z")
      out$ddf <- ddf

   if ((x$test != "z" || is.element(predtype, c("riley","t"))) && predtype != "simple") {
      out$pi.dist <- "t"
      out$pi.ddf <- pi.ddf
   } else {
      out$pi.dist <- "norm"
   }

   out$pi.se <- pi.se

   # add some info to pi.lb

   attr(out$pi.lb, "level") <- level
   attr(out$pi.lb, "dist") <- out$pi.dist
   if (out$pi.dist == "t") {
      attr(out$pi.lb, "ddf") <- out$pi.ddf
   }
   attr(out$pi.lb, "se") <- pi.se

   # for rma.mv models with a GEN structure, remove PI bounds

   if (inherits(object, "rma.mv") && any(is.element(object$struct, c("GEN","GDIAG")))) {
      out$cr.lb <- NULL
      out$cr.ub <- NULL
      out$pi.lb <- NULL
      out$pi.ub <- NULL
      out$tau2.level <- NULL
      out$gamma2.level <- NULL
   }

   # for FE/EE/CE models, remove the PI bounds

   if (is.element(x$method, c("FE","EE","CE"))) {
      out$cr.lb <- NULL
      out$cr.ub <- NULL
      out$pi.lb <- NULL
      out$pi.ub <- NULL
   }

   # for certain transformations, remove the PI bounds

   funlist <- lapply(list(transf.exp.int,  transf.ilogit.int,  transf.iprobit.int,  transf.ztor.int,  transf.iarcsin.int,  transf.iahw.int,  transf.iabt.int, transf.dtocles.int,
                          transf.exp.mode, transf.ilogit.mode, transf.iprobit.mode, transf.ztor.mode, transf.iarcsin.mode, transf.iahw.mode, transf.iabt.mode), deparse)
   if (do.transf && any(sapply(funlist, identical, deparse(transf)))) {
      out$cr.lb <- NULL
      out$cr.ub <- NULL
      out$pi.lb <- NULL
      out$pi.ub <- NULL
   }

   class(out) <- c("predict.rma", "list.rma")

   if (vcov & !do.transf) {
      out <- list(pred=out)
      out$vcov <- vcovpred
   }

   return(out)

}
