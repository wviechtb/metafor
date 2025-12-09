predict.matreg <- function(object, newmods, intercept, addx=FALSE, level, adjust=FALSE, digits, transf, targs, vcov=FALSE, ...) {

   #########################################################################

   mstyle <- .get.mstyle()

   .chkclass(class(object), must="matreg")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   x <- object

   mf <- match.call()

   if (missing(newmods))
      stop(mstyle$stop("Argument 'newmods' must be specified."))

   if (missing(intercept)) {
      intercept <- x$intercept
      int.spec <- FALSE
   } else {
      int.spec <- TRUE
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

   level <- .level(level)

   if (!is.logical(adjust))
      stop(mstyle$stop("Argument 'adjust' must be a logical."))

   rnames <- NULL

   p <- nrow(x$tab) # number of coefficients (also counts the intercept even if it is NA)

   #########################################################################

   if (!(.is.vector(newmods) || inherits(newmods, "matrix")))
      stop(mstyle$stop(paste0("Argument 'newmods' should be a vector or matrix, but is of class '", class(newmods)[1], "'.")))

   singlemod <- (NCOL(newmods) == 1L) && ((!x$intercept && p == 1L) || (x$intercept && p == 2L))

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
      # allow matching of terms by names (note: only possible if all columns in X.new have colnames)
      if (!is.null(colnames(X.new)) && all(colnames(X.new) != "")) {
         colnames.mod <- rownames(x$tab)
         if (x$intercept)
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

   if (!singlemod && ncol(X.new) == p) {

      if (int.spec)
         warning(mstyle$warning("Arguments 'intercept' ignored when 'newmods' includes 'p' columns."), call.=FALSE)

   } else {

      if (x$intercept) {

         if (intercept) {
            X.new <- cbind(intrcpt=1, X.new)
         } else {
            X.new <- cbind(intrcpt=0, X.new)
         }

      }

   }

   if (ncol(X.new) != p)
      stop(mstyle$stop(paste0("Dimensions of 'newmods' (", ncol(X.new), ") do not the match dimensions of the model (", p, ").")))

   if (is.null(X.new))
      stop(mstyle$stop("Matrix 'X.new' is NULL."))

   #########################################################################

   # predicted values, SEs, and confidence intervals

   pred  <- rep(NA_real_, k.new)
   vpred <- rep(NA_real_, k.new)

   for (i in seq_len(k.new)) {
      Xi.new   <- X.new[i,,drop=FALSE]
      beta <- x$tab$beta
      beta[Xi.new == 0] <- 0
      pred[i]  <- Xi.new %*% cbind(beta)
      vb <- x$vb
      vb[Xi.new == 0,] <- 0
      vb[,Xi.new == 0] <- 0
      vpred[i] <- Xi.new %*% tcrossprod(vb, Xi.new)
   }

   if (x$test == "t") {
      crit <- qt(level/ifelse(adjust, 2*k.new, 2), df=x$df.residual, lower.tail=FALSE)
   } else {
      crit <- qnorm(level/ifelse(adjust, 2*k.new, 2), lower.tail=FALSE)
   }

   vpred[vpred < 0] <- NA_real_
   se <- sqrt(vpred)
   ci.lb <- pred - crit * se
   ci.ub <- pred + crit * se

   if (vcov) {
      vcovpred <- symmpart(X.new %*% x$vb %*% t(X.new))
      diag(vcovpred) <- se^2
   }

   #########################################################################

   # apply transformation function if one has been specified

   if (is.function(transf)) {
      if (is.null(targs)) {
         pred  <- sapply(pred, transf)
         se    <- rep(NA_real_, k.new)
         ci.lb <- sapply(ci.lb, transf)
         ci.ub <- sapply(ci.ub, transf)
      } else {
         if (!is.primitive(transf) && !is.null(targs) && length(formals(transf)) == 1L)
            stop(mstyle$stop("Function specified via 'transf' does not appear to have an argument for 'targs'."))
         pred  <- sapply(pred, transf, targs)
         se    <- rep(NA_real_, k.new)
         ci.lb <- sapply(ci.lb, transf, targs)
         ci.ub <- sapply(ci.ub, transf, targs)
      }
      do.transf <- TRUE
   } else {
      do.transf <- FALSE
   }

   # make sure order of intervals is always increasing

   tmp <- .psort(ci.lb, ci.ub)
   ci.lb <- tmp[,1]
   ci.ub <- tmp[,2]

   slab <- seq_len(k.new)

   if (!is.null(rnames))
      slab <- rnames

   # add row/colnames to vcovpred

   if (vcov)
      rownames(vcovpred) <- colnames(vcovpred) <- slab

   # but when predicting just a single value, use "" as study label

   if (k.new == 1L && is.null(rnames))
      slab <- ""

   # handle NAs

   not.na <- rep(TRUE, k.new)

   if (na.act == "na.omit")
      not.na <- !is.na(pred)

   if (na.act == "na.fail" && any(!x$not.na))
      stop(mstyle$stop("Missing values in results."))

   out <- list(pred=pred[not.na], se=se[not.na], ci.lb=ci.lb[not.na], ci.ub=ci.ub[not.na])

   if (vcov)
      vcovpred <- vcovpred[not.na,not.na,drop=FALSE]

   # add X matrix to list

   if (addx) {
      out$X <- matrix(X.new[not.na,], ncol=p)
      colnames(out$X) <- rownames(x$tab)
   }

   # add slab values to list

   out$slab <- slab[not.na]

   # add some additional info

   out$digits <- digits
   out$transf <- do.transf

   if (x$test != "z")
      out$ddf <- x$df.residual

   class(out) <- c("predict.matreg", "list.rma")

   if (vcov & !do.transf) {
      out <- list(pred=out)
      out$vcov <- vcovpred
   }

   return(out)

}
