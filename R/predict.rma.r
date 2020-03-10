predict.rma <- function(object, newmods, intercept, tau2.levels, gamma2.levels, addx=FALSE,
level, digits, transf, targs, vcov=FALSE, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(object, "rma"))
      stop(mstyle$stop("Argument 'object' must be an object of class \"rma\"."))

   if (inherits(object, "rma.ls"))
      stop(mstyle$stop("Method not available for objects of class \"rma.ls\"."))

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   x <- object

   if (missing(newmods))
      newmods <- NULL

   if (missing(intercept))
      intercept <- x$intercept

   if (missing(tau2.levels))
      tau2.levels <- NULL

   if (missing(gamma2.levels))
      gamma2.levels <- NULL

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

   level <- ifelse(level == 0, 1, ifelse(level >= 1, (100-level)/100, ifelse(level > .5, 1-level, level)))

   #########################################################################

   if (is.element(x$test, c("knha","adhoc","t"))) {
      crit <- qt(level/2, df=x$dfs, lower.tail=FALSE)
   } else {
      crit <- qnorm(level/2, lower.tail=FALSE)
   }

   if (x$int.only && !is.null(newmods))
      stop(mstyle$stop("Cannot specify new moderator values for models without moderators."))

   #########################################################################

   ### TODO: can this be simplified? (every time I sit down and stare at the mess below, it gives me a headache)

   if (is.null(newmods)) {

      ### if no new moderator values are specified

      if (!inherits(object, "rma.mv") || (inherits(object, "rma.mv") && any(object$struct=="GEN"))) {

         ### for rma.uni, rma.mh, rma.peto, and rma.glmm objects

         if (x$int.only) {                                                  # if intercept-only model predict only the intercept
            k.new <- 1                                                      #
            X.new <- cbind(1)                                               #
         } else {                                                           # otherwise predict for all k.f studies (including studies with NAs)
            k.new <- x$k.f                                                  #
            X.new <- x$X.f                                                  #
         }                                                                  #

      } else {

         ### for rma.mv objects

         if (x$int.only) {                                                    # if intercept-only model:
            if (!x$withG) {                                                   #   # if there is no G structure (and hence also no H structure)
               k.new <- 1                                                     #   # then we just need to predict the intercept once
               X.new <- cbind(1)                                              #
            }                                                                 #
            if (x$withG && x$withH) {                                         #   # if there is both a G and H structure
               if (is.null(tau2.levels) && is.null(gamma2.levels)) {          #      # and user has not specified tau2s.levels and gamma2.levels
                  k.new <- x$tau2s * x$gamma2s                                #         # then we need to predict intercepts for all combinations of tau2 and gamma2 values
                  X.new <- cbind(rep(1,k.new))                                #
                  if (x$tau2s == 1) {                                         #         # if there is only a single tau^2
                     tau2.levels <- rep(1,k.new)                              #         # then tau2.levels should be 1 repeated k.new times
                  } else {                                                    #
                     tau2.levels <- rep(levels(x$mf.g.f$inner), each=x$gamma2s)         # otherwise repeat actual levels gamma2s times
                  }                                                           #
                  if (x$gamma2s == 1) {                                       #         # if there is only a single gamma^2 value
                     gamma2.levels <- rep(1,k.new)                            #         # then gamma2.levels should be 1 repeated k.new times
                  } else {                                                    #
                     gamma2.levels <- rep(levels(x$mf.h.f$inner), times=x$tau2s)        # otherwise repeat actual levels tau2s times
                  }                                                           #
               }                                                              #
               if ((!is.null(tau2.levels) && is.null(gamma2.levels)) ||       #   # if user specifies only one of tau2.levels and gamma2.levels, throw an error
                   (is.null(tau2.levels) && !is.null(gamma2.levels)))
                  stop(mstyle$stop("Either specify both of 'tau2.levels' and 'gamma2.levels' or neither."))
               if (!is.null(tau2.levels) && !is.null(gamma2.levels)) {        #   # if user has specified both tau2s.levels and gamma2.levels
                  if (length(tau2.levels) != length(gamma2.levels))           #
                     stop(mstyle$stop("Length of 'tau2.levels' and 'gamma2.levels' is not the same."))
                  k.new <- length(tau2.levels)                                #      # then we need to predict intercepts for those level combinations
                  X.new <- cbind(rep(1,k.new))                                #
               }                                                              #
            }                                                                 #
            if (x$withG && !x$withH) {                                        #   # if there is only a G structure (and no H structure)
               if (is.null(tau2.levels)) {                                    #      # and user has not specified tau2.levels
                  k.new <- x$tau2s                                            #         # then we need to predict intercepts for all tau2 values
                  X.new <- cbind(rep(1,k.new))                                #
                  if (x$tau2s == 1) {                                         #
                     tau2.levels <- rep(1, k.new)                             #
                  } else {                                                    #
                     tau2.levels <- levels(x$mf.g.f$inner)                    #
                  }                                                           #
               } else {                                                       #      # and the user has specified tau2.levels
                  k.new <- length(tau2.levels)                                #         # then we need to predict intercepts for those levels
                  X.new <- cbind(rep(1,k.new))                                #
               }                                                              #
               gamma2.levels <- rep(1, k.new)
            }                                                                 #
         } else {                                                             # if not an intercept-only model
            k.new <- x$k.f                                                    #   # then predict for all k.f studies (including studies with NAs)
            X.new <- x$X.f                                                    #
            if (!is.null(tau2.levels) || !is.null(gamma2.levels))             #
               warning(mstyle$warning("Arguments 'tau2.levels' and 'gamma2.levels' ignored when obtaining fitted values."))
            tau2.levels <- as.character(x$mf.g.f$inner)                       #
            gamma2.levels <- as.character(x$mf.h.f$inner)                     #
         }                                                                    #

      }

   } else {

      ### if new moderator values have been specified

      if (!(.is.vector(newmods) || inherits(newmods, "matrix")))
         stop(mstyle$stop(paste0("Argument 'newmods' should be a vector or matrix, but is of class '", class(newmods), "'.")))

      if ((!x$int.incl && x$p == 1L) || (x$int.incl && x$p == 2L)) {
         k.new <- length(newmods)                               # if single moderator (multiple k.new possible) (either without or with intercept in the model)
         X.new <- cbind(c(newmods))                             #
      } else {                                                  # in case the model has more than one predictor:
         if (.is.vector(newmods) || nrow(newmods) == 1L) {      #   # if user gives one vector or one row matrix (only one k.new):
            k.new <- 1                                          #
            X.new <- rbind(newmods)                             #
         } else {                                               #   # if user gives multiple rows and columns (multiple k.new):
            k.new <- nrow(newmods)                              #
            X.new <- cbind(newmods)                             #
         }                                                      #
         ### allow matching of terms by names (note: only possible if all columns in X.new and x$X have colnames)
         if (!is.null(colnames(X.new)) && all(colnames(X.new) != "") && !is.null(colnames(x$X)) && all(colnames(x$X) != "")) {
            colnames.mod <- colnames(x$X)
            if (x$int.incl)
               colnames.mod <- colnames.mod[-1]
            pos <- sapply(colnames(X.new), function(colname) {
                     d <- c(adist(colname, colnames.mod, costs=c(ins=1, sub=Inf, del=Inf))) # compute edit distances with Inf costs for substitutions/deletions
                     if (all(is.infinite(d))) # if there is no match, then all elements are Inf
                        stop(mstyle$stop(paste0("Could not find variable '", colname, "' in the model.")), call. = FALSE)
                     d <- which(d == min(d)) # don't use which.min() since that only finds the first minimum
                     if (length(d) > 1L) # if there is no unique match, then there is more than one minimum
                        stop(mstyle$stop(paste0("Could not match up '", colname, "' uniquely to a variable in the model.")), call. = FALSE)
                     return(d)
                     })
            if (anyDuplicated(pos)) # if the same name is used more than once, then there will be duplicated pos values
               stop(mstyle$stop("Multiple matches for the same variable name."))
            colnames(X.new) <- colnames.mod[pos]
            pos <- sapply(colnames.mod, function(colname) {
                     d <- c(adist(colname, colnames(X.new), costs=c(ins=1, sub=Inf, del=Inf))) # compute edit distances with Inf costs for substitutions/deletions
                     d <- which(d == min(d)) # don't use which.min() since that only finds the first minimum
                     return(d)
                     })
            X.new <- X.new[,pos,drop=FALSE]
         }
      }

      if (inherits(X.new[1,1], "character"))
         stop(mstyle$stop(paste0("Argument 'newmods' should only contain numeric variables.")))

      ### if the user has specified newmods and an intercept was included in the original model, add the intercept to X.new
      ### but user can also decide to remove the intercept from the predictions with intercept=FALSE

      if (x$int.incl) {
         if (intercept) {
            X.new <- cbind(intrcpt=rep(1,k.new), X.new)
         } else {
            X.new <- cbind(intrcpt=rep(0,k.new), X.new)
         }
      }

      if (ncol(X.new) != x$p)
         stop(mstyle$stop("Dimensions of 'newmods' do not match dimensions of the model."))

   }

   #return(list(k.new=k.new, tau2=x$tau2, gamma2=x$gamma2, tau2.levels=tau2.levels, gamma2.levels=gamma2.levels))

   #########################################################################

   ### for rma.mv models with multiple tau^2 values, must use tau2.levels argument when using newmods to obtain credibility intervals

   if (inherits(object, "rma.mv") && x$withG) {

      if (x$tau2s > 1) {

         if (is.null(tau2.levels)) {

            #warning(mstyle$warning("Need to specify 'tau2.levels' argument to obtain credibility intervals."))

         } else {

            ### if tau2.levels argument is a character vector, check that specified tau^2 values actually exist
            if (!is.numeric(tau2.levels) && anyNA(pmatch(tau2.levels, x$g.levels.f[[1]], duplicates.ok=TRUE)))
               stop(mstyle$stop("Non-existing levels specified via 'tau2.levels' argument."))

            ### if tau2.levels argument is numeric, check that specified tau^2 values actually exist
            if (is.numeric(tau2.levels)) {
               tau2.levels <- round(tau2.levels)
               if (any(tau2.levels < 1) || any(tau2.levels > x$g.nlevels.f[1]))
                  stop(mstyle$stop("Non-existing tau^2 values specified via 'tau2.levels' argument."))
            }

            ### allow quick setting of all levels
            if (length(tau2.levels) == 1L)
               tau2.levels <- rep(tau2.levels, k.new)

            ### check length of tau2.levels argument
            if (length(tau2.levels) != k.new)
               stop(mstyle$stop("Length of 'tau2.levels' does not match number of predicted values."))

         }

      } else {

         tau2.levels <- rep(1, k.new)

      }

   }

   ### for rma.mv models with multiple gamma^2 values, must use gamma.levels argument when using newmods to obtain credibility intervals

   if (inherits(object, "rma.mv") && x$withH) {

      if (x$gamma2s > 1) {

         if (is.null(gamma2.levels)) {

            #warning(mstyle$warning("Need to specify 'gamma2.levels' argument to obtain credibility intervals."))

         } else {

            ### if gamma2.levels argument is a character vector, check that specified gamma^2 values actually exist
            if (!is.numeric(gamma2.levels) && anyNA(pmatch(gamma2.levels, x$h.levels.f[[1]], duplicates.ok=TRUE)))
               stop(mstyle$stop("Non-existing levels specified via 'gamma2.levels' argument."))

            ### if gamma2.levels argument is numeric, check that specified gamma^2 values actually exist
            if (is.numeric(gamma2.levels)) {
               gamma2.levels <- round(gamma2.levels)
               if (any(gamma2.levels < 1) || any(gamma2.levels > x$h.nlevels.f[1]))
                  stop(mstyle$stop("Non-existing gamma^2 values specified via 'gamma2.levels' argument."))
            }

            ### allow quick setting of all levels
            if (length(gamma2.levels) == 1L)
               gamma2.levels <- rep(gamma2.levels, k.new)

            ### check length of gamma2.levels argument
            if (length(gamma2.levels) != k.new)
               stop(mstyle$stop("Length of 'gamma2.levels' does not match number of predicted values."))

         }

      } else {

         gamma2.levels <- rep(1, k.new)

      }

   }

   #########################################################################

   ### predicted values, SEs, and confidence intervals

   pred  <- rep(NA_real_, k.new)
   vpred <- rep(NA_real_, k.new)

   for (i in seq_len(k.new)) {
      Xi.new   <- X.new[i,,drop=FALSE]
      pred[i]  <- Xi.new %*% x$beta
      vpred[i] <- Xi.new %*% tcrossprod(x$vb, Xi.new)
   }

   se <- sqrt(vpred)
   ci.lb <- pred - crit * se
   ci.ub <- pred + crit * se

   if (vcov)
      vcovpred <- X.new %*% x$vb %*% t(X.new)

   #########################################################################

   ### credibility/prediction intervals

   if (!inherits(object, "rma.mv")) {

      ### for rma.uni, rma.mh, rma.peto, and rma.glmm objects (in rma.mh and rma.peto, tau2 = 0 by default and stored as such)

      cr.lb <- pred - crit * sqrt(vpred + x$tau2)
      cr.ub <- pred + crit * sqrt(vpred + x$tau2)

   } else {

      ### for rma.mv objects

      if (!x$withG) {

         ### if there is no G structure (and hence no H structure), there are no tau2 and gamma2 values, so just add the sum of all of the sigma2 values

         cr.lb <- pred - crit * sqrt(vpred + sum(x$sigma2))
         cr.ub <- pred + crit * sqrt(vpred + sum(x$sigma2))

      }

      if (x$withG && !x$withH) {

         ### if there is a G structure but no H structure

         if (x$tau2s == 1) {

            ### if there is only a single tau^2 value, always add that (in addition to the sum of all of the sigma^2 values)

            cr.lb <- pred - crit * sqrt(vpred + sum(x$sigma2) + x$tau2)
            cr.ub <- pred + crit * sqrt(vpred + sum(x$sigma2) + x$tau2)

         } else {

            if (is.null(tau2.levels)) {

               ### if user has not specified tau2.levels, cannot compute bounds

               cr.lb <- rep(NA, k.new)
               cr.ub <- rep(NA, k.new)
               tau2.levels <- rep(NA, k.new)

            } else {

               ### if there are multiple tau^2 values, either let user define numerically which value(s) to use or
               ### match the position of the specified tau2.levels to the levels of the inner factor in the model

               if (!is.numeric(tau2.levels))
                  tau2.levels <- pmatch(tau2.levels, x$g.levels.f[[1]], duplicates.ok=TRUE)

               cr.lb <- pred - crit * sqrt(vpred + sum(x$sigma2) + x$tau2[tau2.levels])
               cr.ub <- pred + crit * sqrt(vpred + sum(x$sigma2) + x$tau2[tau2.levels])
               tau2.levels <- x$g.levels.f[[1]][tau2.levels]

            }

         }

      }

      if (x$withG && x$withH) {

         ### if there is a G structure and an H structure

         if (x$tau2s == 1 && x$gamma2s == 1) {

            ### if there is only a single tau^2 and gamma^2 value, always add that (in addition to the sum of all of the sigma^2 values)

            cr.lb <- pred - crit * sqrt(vpred + sum(x$sigma2) + x$tau2 + x$gamma2)
            cr.ub <- pred + crit * sqrt(vpred + sum(x$sigma2) + x$tau2 + x$gamma2)

         } else {

            if (is.null(tau2.levels) || is.null(gamma2.levels)) {

               ### if user has not specified tau2.levels and gamma2.levels, cannot compute bounds

               cr.lb <- rep(NA, k.new)
               cr.ub <- rep(NA, k.new)
               tau2.levels <- rep(NA, k.new)
               gamma2.levels <- rep(NA, k.new)

            } else {

               ### if there are multiple tau^2 and/or gamma^2 values, either let user define numerically which value(s) to use or
               ### match the position of the specified tau2.levels and gamma2.levels to the levels of the inner factors in the model

               if (!is.numeric(tau2.levels))
                  tau2.levels <- pmatch(tau2.levels, x$g.levels.f[[1]], duplicates.ok=TRUE)
               if (!is.numeric(gamma2.levels))
                  gamma2.levels <- pmatch(gamma2.levels, x$h.levels.f[[1]], duplicates.ok=TRUE)

               cr.lb <- pred - crit * sqrt(vpred + sum(x$sigma2) + x$tau2[tau2.levels] + x$gamma2[gamma2.levels])
               cr.ub <- pred + crit * sqrt(vpred + sum(x$sigma2) + x$tau2[tau2.levels] + x$gamma2[gamma2.levels])
               tau2.levels <- x$g.levels.f[[1]][tau2.levels]
               gamma2.levels <- x$h.levels.f[[1]][gamma2.levels]

            }

         }

      }

   }

   #########################################################################

   ### apply transformation function if one has been specified

   if (is.function(transf)) {
      if (is.null(targs)) {
         pred  <- sapply(pred, transf)
         se    <- rep(NA,k.new)
         ci.lb <- sapply(ci.lb, transf)
         ci.ub <- sapply(ci.ub, transf)
         cr.lb <- sapply(cr.lb, transf)
         cr.ub <- sapply(cr.ub, transf)
      } else {
         pred  <- sapply(pred, transf, targs)
         se    <- rep(NA,k.new)
         ci.lb <- sapply(ci.lb, transf, targs)
         ci.ub <- sapply(ci.ub, transf, targs)
         cr.lb <- sapply(cr.lb, transf, targs)
         cr.ub <- sapply(cr.ub, transf, targs)
      }
      transf <- TRUE
   }

   ### make sure order of intervals is always increasing

   tmp <- .psort(ci.lb, ci.ub)
   ci.lb <- tmp[,1]
   ci.ub <- tmp[,2]

   tmp <- .psort(cr.lb, cr.ub)
   cr.lb <- tmp[,1]
   cr.ub <- tmp[,2]

   ### use study labels from the object when the model has moderators and no new moderators have been specified
   ### otherwise, just use consecutive numbers to label the predicted values

   if (is.null(newmods) && !x$int.only) {
      slab <- x$slab
   } else {
      slab <- seq_len(k.new)
   }

   ### add row/colnames to vcovpred

   if (vcov)
      rownames(vcovpred) <- colnames(vcovpred) <- slab

   ### but when predicting just a single value, use "" as study label

   if (k.new == 1L)
      slab <- ""

   ### handle NAs

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

   out <- list(pred=pred[not.na], se=se[not.na], ci.lb=ci.lb[not.na], ci.ub=ci.ub[not.na], cr.lb=cr.lb[not.na], cr.ub=cr.ub[not.na])

   if (vcov)
      vcovpred <- vcovpred[not.na,not.na,drop=FALSE]

   if (na.act == "na.exclude" && is.null(newmods) && !x$int.only) {

      out <- lapply(out, function(val) ifelse(x$not.na, val, NA))

      if (vcov) {
         vcovpred[!x$not.na,] <- NA
         vcovpred[,!x$not.na] <- NA
      }

   }

   ### add tau2.levels values to list

   if (inherits(object, "rma.mv") && x$withG && x$tau2s > 1)
      out$tau2.level <- tau2.levels

   ### add gamma2.levels values to list

   if (inherits(object, "rma.mv") && x$withH && x$gamma2s > 1)
      out$gamma2.level <- gamma2.levels


   ### remove cr part for models with a GEN structure
   if (inherits(object, "rma.mv") && any(object$struct=="GEN")) {
      out$cr.lb <- NULL
      out$cr.ub <- NULL
      out$tau2.level <- NULL
      out$gamma2.level <- NULL
   }

   ### add X matrix to list

   if(addx)
      out$X <- matrix(X.new[not.na,], ncol=x$p)

   if (addx)
      colnames(out$X) <- colnames(x$X)

   ### add slab values to list

   out$slab <- slab[not.na]

   ### for FE models, remove the columns corresponding to the credibility interval bounds

   if (x$method == "FE") {
      out$cr.lb <- NULL
      out$cr.ub <- NULL
   }

   out$digits <- digits
   out$method <- x$method
   out$transf <- transf

   class(out) <- "list.rma"

   if (vcov & !transf) {
      out <- list(pred=out)
      out$vcov <- vcovpred
   }

   return(out)

}
