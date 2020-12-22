predict.rma.ls <- function(object, newmods, intercept, addx=FALSE, newscale, addz=FALSE,
level, digits, transf, targs, vcov=FALSE, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(object), must="rma.ls")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   x <- object

   if (missing(newmods))
      newmods <- NULL

   if (missing(intercept))
      intercept <- x$intercept

   if (missing(newscale))
      newscale <- NULL

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

   ddd <- list(...)

   .chkdots(ddd, c("pi.type"))

   if (is.null(ddd$pi.type)) {
      pi.type <- "default"
   } else {
      pi.type <- ddd$pi.type
   }

   if (x$int.only && !is.null(newmods))
      stop(mstyle$stop("Cannot specify new moderator values for models without moderators."))

   if (!(x$Z.int.only && identical(newscale, 1)) && x$Z.int.only && !is.null(newscale))
      stop(mstyle$stop("Cannot specify new scale values for models without scale variables."))

   #########################################################################

   if (!is.null(newmods)) {

      ### if newmods has been specified

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
                        stop(mstyle$stop(paste0("Could not match up variable '", colname, "' uniquely to a variable in the model.")), call. = FALSE)
                     return(d)
                     })
            if (anyDuplicated(pos)) { # if the same name is used more than once, then there will be duplicated pos values
               dups <- paste(unique(colnames(X.new)[duplicated(pos)]), collapse=", ")
               stop(mstyle$stop(paste0("Found multiple matches for the same variable name (", dups, ").")))
            }
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
            X.new <- cbind(intrcpt=1, X.new)
         } else {
            X.new <- cbind(intrcpt=0, X.new)
         }
      }

      if (ncol(X.new) != x$p)
         stop(mstyle$stop(paste0("Dimensions of 'newmods' (", ncol(X.new), ") do not match dimensions of the model (", x$p, ").")))

   }

   if (!is.null(newscale)) {

      if (!(.is.vector(newscale) || inherits(newscale, "matrix")))
         stop(mstyle$stop(paste0("Argument 'newscale' should be a vector or matrix, but is of class '", class(newscale), "'.")))

      if ((!x$Z.int.incl && x$q == 1L) || (x$Z.int.incl && x$q == 2L)) {
         Z.k.new <- length(newscale)                            # if single moderator (multiple k.new possible) (either without or with intercept in the model)
         Z.new <- cbind(c(newscale))                            #
      } else {                                                  # in case the model has more than one predictor:
         if (.is.vector(newscale) || nrow(newscale) == 1L) {    #   # if user gives one vector or one row matrix (only one k.new):
            Z.k.new <- 1                                        #
            Z.new <- rbind(newscale)                            #
         } else {                                               #   # if user gives multiple rows and columns (multiple k.new):
            Z.k.new <- nrow(newscale)                           #
            Z.new <- cbind(newscale)                            #
         }                                                      #
         ### allow matching of terms by names (note: only possible if all columns in Z.new and x$Z have colnames)
         if (!is.null(colnames(Z.new)) && all(colnames(Z.new) != "") && !is.null(colnames(x$Z)) && all(colnames(x$Z) != "")) {
            colnames.mod <- colnames(x$Z)
            if (x$Z.int.incl)
               colnames.mod <- colnames.mod[-1]
            pos <- sapply(colnames(Z.new), function(colname) {
                     d <- c(adist(colname, colnames.mod, costs=c(ins=1, sub=Inf, del=Inf))) # compute edit distances with Inf costs for substitutions/deletions
                     if (all(is.infinite(d))) # if there is no match, then all elements are Inf
                        stop(mstyle$stop(paste0("Could not find variable '", colname, "' from 'newscale' in the model.")), call. = FALSE)
                     d <- which(d == min(d)) # don't use which.min() since that only finds the first minimum
                     if (length(d) > 1L) # if there is no unique match, then there is more than one minimum
                        stop(mstyle$stop(paste0("Could not match up variable '", colname, "' from 'newscale' uniquely to a variable in the model.")), call. = FALSE)
                     return(d)
                     })
            if (anyDuplicated(pos)) { # if the same name is used more than once, then there will be duplicated pos values
               dups <- paste(unique(colnames(Z.new)[duplicated(pos)]), collapse=", ")
               stop(mstyle$stop(paste0("Found multiple matches for the same variable name (", dups, ") in 'newscale'.")))
            }
            colnames(Z.new) <- colnames.mod[pos]
            pos <- sapply(colnames.mod, function(colname) {
                     d <- c(adist(colname, colnames(Z.new), costs=c(ins=1, sub=Inf, del=Inf))) # compute edit distances with Inf costs for substitutions/deletions
                     d <- which(d == min(d)) # don't use which.min() since that only finds the first minimum
                     return(d)
                     })
            Z.new <- Z.new[,pos,drop=FALSE]
         }
      }

      if (inherits(Z.new[1,1], "character"))
         stop(mstyle$stop(paste0("Argument 'newscale' should only contain numeric variables.")))

      ### if the user has specified newscale and an intercept was included in the original model, add the intercept to Z.new
      ### one special case: when the scale model is an intercept-only model, one can set newscale=1 to obtain the predicted
      ### intercept (which can be converted to tau^2 with transf=exp when using a log link)

      if (!(x$Z.int.only && dim(Z.new) == c(1L,1L) && Z.new == 1) && x$Z.int.incl)
         Z.new <- cbind(intrcpt=1, Z.new)

      if (ncol(Z.new) != x$q)
         stop(mstyle$stop(paste0("Dimensions of 'newscale' (", ncol(Z.new), ") do not match dimensions of the scale model (", x$q, ").")))

   }

   # four possibilities:
   # 1) newmods not specified, newscale not specified: get the fitted values of the studies and ci/pi bounds thereof
   # 2) newmods     specified, newscale not specified: get the predicted mu values for these newmods values and ci bounds thereof
   #                                                      (note: cannot compute pi bounds, since the tau^2 values cannot be predicted)
   # 3) newmods not specified, newscale     specified: get the predicted alpha values and ci bounds thereof
   #                                                      (transf=exp to obtain predicted tau^2 values when using the default log link)
   # 4) newmods     specified, newscale     specified: get the predicted mu values for these newmods values and ci/pi bounds thereof

   pred.mui <- TRUE

   if (is.null(newmods)) {

      if (is.null(newscale)) {

         k.new  <- x$k.f
         X.new  <- x$X.f
         Z.new  <- x$Z.f
         tau2.f <- x$tau2.f

      } else {

         k.new    <- Z.k.new
         addx     <- FALSE
         pred.mui <- FALSE

      }

   } else {

      if (is.null(newscale)) {

         Z.new  <- matrix(NA, nrow=k.new, ncol=x$q)
         tau2.f <- rep(NA, k.new)
         addz   <- FALSE

      } else {

         tau2.f <- rep(NA_real_, Z.k.new)

         for (i in seq_len(Z.k.new)) {
            Zi.new    <- Z.new[i,,drop=FALSE]
            tau2.f[i] <- Zi.new %*% x$alpha
         }

         if (x$link == "log") {
            tau2.f <- exp(tau2.f)
         } else {
            if (any(tau2.f < 0)) {
               warning(mstyle$warning(paste0("Negative predicted 'tau2' values constrained to 0.")), call.=FALSE)
               tau2.f[tau2.f < 0] <- 0
            }
         }

         if (length(tau2.f) == 1L) {
            Z.new  <- Z.new[rep(1,k.new),]
            tau2.f <- rep(tau2.f, k.new)
         }

         if (length(tau2.f) != k.new)
            stop(mstyle$stop(paste0("Dimensions of 'newmods' (", k.new, ") do not match dimensions of newscale (", length(tau2.f), ").")))

      }

   }

   #return(list(k.new=k.new, tau2=x$tau2, gamma2=x$gamma2, tau2.levels=tau2.levels, gamma2.levels=gamma2.levels))

   #########################################################################

   ### predicted values, SEs, and confidence intervals

   pred  <- rep(NA_real_, k.new)
   vpred <- rep(NA_real_, k.new)

   if (pred.mui) {

      for (i in seq_len(k.new)) {
         Xi.new   <- X.new[i,,drop=FALSE]
         pred[i]  <- Xi.new %*% x$beta
         vpred[i] <- Xi.new %*% tcrossprod(x$vb, Xi.new)
      }

      if (is.element(x$test, c("t"))) {
         crit <- qt(level/2, df=x$dfs, lower.tail=FALSE)
      } else {
         crit <- qnorm(level/2, lower.tail=FALSE)
      }

   } else {

      for (i in seq_len(k.new)) {
         Zi.new   <- Z.new[i,,drop=FALSE]
         pred[i]  <- Zi.new %*% x$alpha
         vpred[i] <- Zi.new %*% tcrossprod(x$vb.alpha, Zi.new)
      }

      if (is.element(x$test, c("t"))) {
         crit <- qt(level/2, df=x$dfs.alpha, lower.tail=FALSE)
      } else {
         crit <- qnorm(level/2, lower.tail=FALSE)
      }

   }

   se <- sqrt(vpred)
   ci.lb <- pred - crit * se
   ci.ub <- pred + crit * se

   #########################################################################

   if (pred.mui) {

      if (vcov)
         vcovpred <- X.new %*% x$vb %*% t(X.new)

      if (pi.type == "simple")
         crit <- qnorm(level/2, lower.tail=FALSE)

      if (pi.type == "riley") {
         dfs <- x$k - x$p - 1
         if (dfs <= 0) {
            crit <- Inf
         } else {
            crit <- qt(level/2, df=x$k-x$p-1, lower.tail=FALSE)
         }
      }

      if (pi.type == "simple")
         vpred <- 0

      ### prediction intervals

      pi.lb <- pred - crit * sqrt(vpred + tau2.f)
      pi.ub <- pred + crit * sqrt(vpred + tau2.f)

   } else {

      if (vcov)
         vcovpred <- Z.new %*% x$vb.alpha %*% t(Z.new)

      pi.lb <- NA
      pi.ub <- NA

   }

   #########################################################################

   ### apply transformation function if one has been specified

   if (is.function(transf)) {
      if (is.null(targs)) {
         pred  <- sapply(pred, transf)
         se    <- rep(NA,k.new)
         ci.lb <- sapply(ci.lb, transf)
         ci.ub <- sapply(ci.ub, transf)
         pi.lb <- sapply(pi.lb, transf)
         pi.ub <- sapply(pi.ub, transf)
      } else {
         pred  <- sapply(pred, transf, targs)
         se    <- rep(NA,k.new)
         ci.lb <- sapply(ci.lb, transf, targs)
         ci.ub <- sapply(ci.ub, transf, targs)
         pi.lb <- sapply(pi.lb, transf, targs)
         pi.ub <- sapply(pi.ub, transf, targs)
      }
      do.transf <- TRUE
   } else {
      do.transf <- FALSE
   }

   ### make sure order of intervals is always increasing

   tmp <- .psort(ci.lb, ci.ub)
   ci.lb <- tmp[,1]
   ci.ub <- tmp[,2]

   tmp <- .psort(pi.lb, pi.ub)
   pi.lb <- tmp[,1]
   pi.ub <- tmp[,2]

   ### when predicting tau^2 values, set negative tau^2 values and CI bounds to 0

   if (!pred.mui && x$link=="identity" && !is.function(transf)) {
      if (any(pred < 0))
         warning(mstyle$warning(paste0("Negative predicted 'tau2' values constrained to 0.")), call.=FALSE)
      pred[pred < 0]   <- 0
      ci.lb[ci.lb < 0] <- 0
      ci.ub[ci.ub < 0] <- 0
   }

   ### use study labels from the object when the model has moderators and no new moderators have been specified

   if (pred.mui) {
      if (is.null(newmods)) {
         slab <- x$slab
      } else {
         slab <- seq_len(k.new)
      }
   } else {
      if (is.null(newscale)) {
         slab <- x$slab
      } else {
         slab <- seq_len(k.new)
      }
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
      if (pred.mui) {
         if (is.null(newmods)) {
            not.na <- x$not.na
         } else {
            not.na <- !is.na(pred)
         }
      } else {
         if (is.null(newscale)) {
            not.na <- x$not.na
         } else {
            not.na <- !is.na(pred)
         }
      }
   }

   if (na.act == "na.fail" && any(!x$not.na))
      stop(mstyle$stop("Missing values in results."))

   out <- list(pred=pred[not.na], se=se[not.na], ci.lb=ci.lb[not.na], ci.ub=ci.ub[not.na], pi.lb=pi.lb[not.na], pi.ub=pi.ub[not.na], cr.lb=pi.lb[not.na], cr.ub=pi.ub[not.na])

   if (vcov)
      vcovpred <- vcovpred[not.na,not.na,drop=FALSE]

   if (na.act == "na.exclude" && is.null(newmods)) {

      out <- lapply(out, function(val) ifelse(x$not.na, val, NA))

      if (vcov) {
         vcovpred[!x$not.na,] <- NA
         vcovpred[,!x$not.na] <- NA
      }

   }

   ### add X matrix to list

   if (addx) {
      out$X <- matrix(X.new[not.na,], ncol=x$p)
      colnames(out$X) <- colnames(x$X)
   }

   ### add Z matrix to list

   if (addz) {
      out$Z <- matrix(Z.new[not.na,], ncol=x$q)
      colnames(out$Z) <- colnames(x$Z)
   }

   ### add slab values to list

   out$slab <- slab[not.na]

   ### for FE models, remove the columns corresponding to the prediction interval bounds

   if (x$method == "FE" || !pred.mui) {
      out$cr.lb <- NULL
      out$cr.ub <- NULL
      out$pi.lb <- NULL
      out$pi.ub <- NULL
   }

   out$digits <- digits
   out$method <- x$method
   out$transf <- do.transf

   class(out) <- "list.rma"

   if (vcov & !do.transf) {
      out <- list(pred=out)
      out$vcov <- vcovpred
   }

   return(out)

}
