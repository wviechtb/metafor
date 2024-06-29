predict.rma.ls <- function(object, newmods, intercept, addx=FALSE, newscale, addz=FALSE,
level, digits, transf, targs, vcov=FALSE, ...) {

   #########################################################################

   mstyle <- .get.mstyle()

   .chkclass(class(object), must="rma.ls")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   x <- object

   mf <- match.call()

   if (any(grepl("pairwise(", as.character(mf), fixed=TRUE))) {
      try(assign("pairwise", object, envir=.metafor), silent=TRUE)
      on.exit(suppressWarnings(rm("pairwise", envir=.metafor)))
   }

   if (missing(newmods))
      newmods <- NULL

   if (missing(intercept)) {
      intercept <- x$intercept
      int.spec <- FALSE
   } else {
      int.spec <- TRUE
   }

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

   level <- .level(level)

   ddd <- list(...)

   .chkdots(ddd, c("pi.type", "newvi"))

   pi.type <- .chkddd(ddd$pi.type, "default", tolower(ddd$pi.type))

   if (!is.null(newmods) && x$int.only && !(x$int.only && identical(newmods, 1)))
      stop(mstyle$stop("Cannot specify new moderator values for models without moderators."))

   if (!is.null(newscale) && x$Z.int.only && !(x$Z.int.only && identical(newscale, 1)))
      stop(mstyle$stop("Cannot specify new scale values for models without scale variables."))

   rnames <- NULL

   #########################################################################

   if (!is.null(newmods)) {

      ### if newmods has been specified

      if (!(.is.vector(newmods) || inherits(newmods, "matrix")))
         stop(mstyle$stop(paste0("Argument 'newmods' should be a vector or matrix, but is of class '", class(newmods), "'.")))

      singlemod <- (NCOL(newmods) == 1L) && ((!x$int.incl && x$p == 1L) || (x$int.incl && x$p == 2L))

      if (singlemod) {                                           # if single moderator (multiple k.new possible) (either without or with intercept in the model)
         k.new <- length(newmods)                                # (but when specifying a matrix, it must be a column vector for this work)
         X.new <- cbind(c(newmods))                              #
         if (.is.vector(newmods)) {                              #
            rnames <- names(newmods)                             #
         } else {                                                #
            rnames <- rownames(newmods)                          #
         }                                                       #
      } else {                                                   # in case the model has more than one predictor:
         if (.is.vector(newmods) || nrow(newmods) == 1L) {       #   # if user gives one vector or one row matrix (only one k.new):
            k.new <- 1                                           #
            X.new <- rbind(newmods)                              #
            if (inherits(newmods, "matrix"))                     #
               rnames <- rownames(newmods)                       #
         } else {                                                #   # if user gives multiple rows and columns (multiple k.new):
            k.new <- nrow(newmods)                               #
            X.new <- cbind(newmods)                              #
            rnames <- rownames(newmods)                          #
         }                                                       #
         ### allow matching of terms by names (note: only possible if all columns in X.new and x$X have colnames)
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
         stop(mstyle$stop(paste0("Argument 'newmods' should only contain numeric variables.")))

      ### if the user has specified newmods and an intercept was included in the original model, add the intercept to X.new
      ### but user can also decide to remove the intercept from the predictions with intercept=FALSE
      ### one special case: when the location model is an intercept-only model, one can set newmods=1 to obtain the predicted intercept

      if (!singlemod && ncol(X.new) == x$p) {

         if (int.spec)
            warning(mstyle$warning("Arguments 'intercept' ignored when 'newmods' includes 'p' columns."), call.=FALSE)

      } else {

         if (x$int.incl && !(x$int.only && ncol(X.new) == 1L && nrow(X.new) == 1L && X.new[1,1] == 1)) {

            if (intercept) {
               X.new <- cbind(intrcpt=1, X.new)
            } else {
               X.new <- cbind(intrcpt=0, X.new)
            }

         }

      }

      if (ncol(X.new) != x$p)
         stop(mstyle$stop(paste0("Dimensions of 'newmods' (", ncol(X.new), ") do not match the dimensions of the model (", x$p, ").")))

   }

   if (!is.null(newscale)) {

      if (!(.is.vector(newscale) || inherits(newscale, "matrix")))
         stop(mstyle$stop(paste0("Argument 'newscale' should be a vector or matrix, but is of class '", class(newscale), "'.")))

      singlescale <- (NCOL(newscale) == 1L) && ((!x$Z.int.incl && x$q == 1L) || (x$Z.int.incl && x$q == 2L))

      if (singlescale) {                                         # if single moderator (multiple k.new possible) (either without or with intercept in the model)
         Z.k.new <- length(newscale)                             #
         Z.new <- cbind(c(newscale))                             #
         if (is.null(rnames)) {                                  #
            if (.is.vector(newscale)) {                          #
               rnames <- names(newscale)                         #
            } else {                                             #
               rnames <- rownames(newscale)                      #
            }                                                    #
         }                                                       #
      } else {                                                   # in case the model has more than one predictor:
         if (.is.vector(newscale) || nrow(newscale) == 1L) {     #   # if user gives one vector or one row matrix (only one k.new):
            Z.k.new <- 1                                         #
            Z.new <- rbind(newscale)                             #
            if (is.null(rnames) && inherits(newscale, "matrix")) #
               rnames <- rownames(newscale)                      #
         } else {                                                #   # if user gives multiple rows and columns (multiple k.new):
            Z.k.new <- nrow(newscale)                            #
            Z.new <- cbind(newscale)                             #
            if (is.null(rnames))                                 #
               rnames <- rownames(newscale)                      #
         }                                                       #
         ### allow matching of terms by names (note: only possible if all columns in Z.new and x$Z have colnames)
         if (!is.null(colnames(Z.new)) && all(colnames(Z.new) != "") && !is.null(colnames(x$Z)) && all(colnames(x$Z) != "")) {
            colnames.mod <- colnames(x$Z)
            if (x$Z.int.incl)
               colnames.mod <- colnames.mod[-1]
            pos <- sapply(colnames(Z.new), function(colname) {
                     d <- c(adist(colname, colnames.mod, costs=c(ins=1, sub=Inf, del=Inf))) # compute edit distances with Inf costs for substitutions/deletions
                     if (all(is.infinite(d))) # if there is no match, then all elements are Inf
                        stop(mstyle$stop(paste0("Could not find variable '", colname, "' from 'newscale' in the model.")))
                     d <- which(d == min(d)) # don't use which.min() since that only finds the first minimum
                     if (length(d) > 1L) # if there is no unique match, then there is more than one minimum
                        stop(mstyle$stop(paste0("Could not match up variable '", colname, "' from 'newscale' uniquely to a variable in the model.")))
                     return(d)
                     })
            if (anyDuplicated(pos)) { # if the same name is used more than once, then there will be duplicated pos values
               dups <- paste(unique(colnames(Z.new)[duplicated(pos)]), collapse=", ")
               stop(mstyle$stop(paste0("Found multiple matches for the same variable name (", dups, ") in 'newscale'.")))
            }
            if (length(pos) != length(colnames.mod)) {
               no.match <- colnames.mod[seq_along(colnames.mod)[-pos]]
               if (length(no.match) > 3L)
                  stop(mstyle$stop(paste0("Argument 'newscale' does not specify values for these variables: ", paste0(no.match[1:3], collapse=", "), ", ...")))
               if (length(no.match) > 1L)
                  stop(mstyle$stop(paste0("Argument 'newscale' does not specify values for these variables: ", paste0(no.match, collapse=", "))))
               if (length(no.match) == 1L)
                  stop(mstyle$stop(paste0("Argument 'newscale' does not specify values for this variable: ", no.match)))
            }
            Z.new <- Z.new[,order(pos),drop=FALSE]
            colnames(Z.new) <- colnames.mod
         }
      }

      if (inherits(Z.new[1,1], "character"))
         stop(mstyle$stop(paste0("Argument 'newscale' should only contain numeric variables.")))

      ### if the user has specified newscale and an intercept was included in the original model, add the intercept to Z.new
      ### but user can also decide to remove the intercept from the predictions with intercept=FALSE (only when predicting log(tau^2))
      ### one special case: when the scale model is an intercept-only model, one can set newscale=1 to obtain the predicted intercept
      ### (which can be converted to tau^2 with transf=exp when using a log link)

      if (!singlescale && ncol(Z.new) == x$q) {

         if (int.spec)
            warning(mstyle$warning("Arguments 'intercept' ignored when 'newscale' includes 'q' columns."), call.=FALSE)

      } else {

         if (x$Z.int.incl && !(x$Z.int.only && ncol(Z.new) == 1L && nrow(Z.new) == 1L && Z.new[1,1] == 1)) {
            if (is.null(newmods)) {
               if (intercept) {
                  Z.new <- cbind(intrcpt=1, Z.new)
               } else {
                  Z.new <- cbind(intrcpt=0, Z.new)
               }
            } else {
               Z.new <- cbind(intrcpt=1, Z.new)
            }
         }

      }

      if (ncol(Z.new) != x$q)
         stop(mstyle$stop(paste0("Dimensions of 'newscale' (", ncol(Z.new), ") do not match the dimensions of the scale model (", x$q, ").")))

   }

   # four possibilities for location-scale models:
   # 1) newmods not specified, newscale not specified: get the fitted values of the studies and ci/pi bounds thereof
   # 2) newmods     specified, newscale not specified: get the predicted mu values for these newmods values and ci bounds thereof
   #                                                   (note: cannot compute pi bounds, since the tau^2 values cannot be predicted)
   # 3) newmods not specified, newscale     specified: get the predicted log(tau^2) (or tau^2) values and ci bounds thereof
   #                                                   (transf=exp to obtain predicted tau^2 values when using the default log link)
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

         Z.new  <- matrix(NA_real_, nrow=k.new, ncol=x$q)
         tau2.f <- rep(NA_real_, k.new)
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

         if (k.new == 1L && Z.k.new > 1L) {
            X.new <- X.new[rep(1,Z.k.new),,drop=FALSE]
            k.new <- Z.k.new
         }

         if (length(tau2.f) == 1L && k.new > 1L) {
            Z.new  <- Z.new[rep(1,k.new),,drop=FALSE]
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

      ddf <- ifelse(is.na(x$ddf), x$k - x$p, x$ddf)

      for (i in seq_len(k.new)) {
         Xi.new   <- X.new[i,,drop=FALSE]
         pred[i]  <- Xi.new %*% x$beta
         vpred[i] <- Xi.new %*% tcrossprod(x$vb, Xi.new)
      }

      if (is.element(x$test, c("knha","adhoc","t"))) {
         crit <- if (ddf > 0) qt(level/2, df=ddf, lower.tail=FALSE) else NA_real_
      } else {
         crit <- qnorm(level/2, lower.tail=FALSE)
      }

   } else {

      ddf <- ifelse(is.na(x$ddf.alpha), x$k - x$q, x$ddf.alpha)

      for (i in seq_len(k.new)) {
         Zi.new   <- Z.new[i,,drop=FALSE]
         pred[i]  <- Zi.new %*% x$alpha
         vpred[i] <- Zi.new %*% tcrossprod(x$va, Zi.new)
      }

      if (is.element(x$test, c("knha","adhoc","t"))) {
         crit <- if (ddf > 0) qt(level/2, df=ddf, lower.tail=FALSE) else NA_real_
      } else {
         crit <- qnorm(level/2, lower.tail=FALSE)
      }

   }

   vpred[vpred < 0] <- NA_real_
   se <- sqrt(vpred)
   ci.lb <- pred - crit * se
   ci.ub <- pred + crit * se

   #########################################################################

   if (pred.mui) {

      if (vcov)
         vcovpred <- symmpart(X.new %*% x$vb %*% t(X.new))

      if (pi.type == "simple") {
         crit <- qnorm(level/2, lower.tail=FALSE)
         vpred <- 0
      }

      pi.ddf <- ddf

      if (is.element(pi.type, c("riley","t"))) {
         if (pi.type == "riley")
            pi.ddf <- x$k - x$p - x$q
         if (pi.type == "t")
            pi.ddf <- x$k - x$p
         pi.ddf[pi.ddf < 1] <- 1
         crit <- qt(level/2, df=pi.ddf, lower.tail=FALSE)
      }

      if (is.null(ddd$newvi)) {
         newvi <- 0
      } else {
         newvi <- ddd$newvi
         newvi <- .expand1(newvi, k.new)
         if (length(newvi) != k.new)
            stop(mstyle$stop(paste0("Length of 'newvi' argument (", length(newvi), ") does not match the number of predicted values (", k.new, ").")))
      }

      ### prediction intervals

      pi.lb <- pred - crit * sqrt(vpred + tau2.f + newvi)
      pi.ub <- pred + crit * sqrt(vpred + tau2.f + newvi)

   } else {

      if (vcov)
         vcovpred <- symmpart(Z.new %*% x$va %*% t(Z.new))

      pi.lb <- NA_real_
      pi.ub <- NA_real_

   }

   #########################################################################

   ### apply transformation function if one has been specified

   if (is.function(transf)) {
      if (is.null(targs)) {
         pred  <- sapply(pred, transf)
         se    <- rep(NA_real_, k.new)
         ci.lb <- sapply(ci.lb, transf)
         ci.ub <- sapply(ci.ub, transf)
         pi.lb <- sapply(pi.lb, transf)
         pi.ub <- sapply(pi.ub, transf)
      } else {
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
         if (!is.null(rnames))
            slab <- rnames
      }
   } else {
      if (is.null(newscale)) {
         slab <- x$slab
      } else {
         slab <- seq_len(k.new)
         if (!is.null(rnames))
            slab <- rnames
      }
   }

   ### add row/colnames to vcovpred

   if (vcov)
      rownames(vcovpred) <- colnames(vcovpred) <- slab

   ### but when predicting just a single value, use "" as study label

   if (k.new == 1L && is.null(rnames))
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

      out <- lapply(out, function(val) ifelse(x$not.na, val, NA_real_))

      if (vcov) {
         vcovpred[!x$not.na,] <- NA_real_
         vcovpred[,!x$not.na] <- NA_real_
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

   ### for FE/EE/CE models, remove the columns corresponding to the prediction interval bounds

   if (is.element(x$method, c("FE","EE","CE")) || !pred.mui) {
      out$cr.lb <- NULL
      out$cr.ub <- NULL
      out$pi.lb <- NULL
      out$pi.ub <- NULL
   }

   out$digits <- digits
   out$method <- x$method
   out$transf <- do.transf
   out$pred.type <- ifelse(pred.mui, "location", "scale")

   if (x$test != "z")
      out$ddf <- ddf

   if (pred.mui && (x$test != "z" || is.element(pi.type, c("riley","t"))) && pi.type != "simple")
      out$pi.ddf <- pi.ddf

   class(out) <- c("predict.rma", "list.rma")

   if (vcov & !do.transf) {
      out <- list(pred=out)
      out$vcov <- vcovpred
   }

   return(out)

}
