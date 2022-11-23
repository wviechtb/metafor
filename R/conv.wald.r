conv.wald <- function(data, out, ci.lb, ci.ub, zval, pval, level=95, transf, include, checkci=TRUE, verbose=FALSE, ...) {

   # add a ni argument to add to attributes(x$yi)?
   # what about allowing t-distribution / dfs?

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (missing(data))
      stop(mstyle$stop("Must specify 'data' argument."))

   .chkclass(class(data), must="escalc")

   x <- data

   if (missing(transf))
      transf <- FALSE

   ### get ... argument and check for extra/superfluous arguments

   ddd <- list(...)

   .chkdots(ddd, c("cifac"))

   cifac <- ifelse(is.null(ddd$cifac), 0.1, ddd$cifac)

   #########################################################################

   if (!is.null(attr(x, "yi.names"))) { # if yi.names attributes is available
      yi.name <- attr(x, "yi.names")[1] # take the first entry to be the yi variable
   } else {                             # if not, see if 'yi' is in the object and assume that is the yi variable
      if (!is.element("yi", names(x)))
         stop(mstyle$stop("Cannot determine name of the 'yi' variable."))
      yi.name <- "yi"
   }
   if (!is.null(attr(x, "vi.names"))) { # if vi.names attributes is available
      vi.name <- attr(x, "vi.names")[1] # take the first entry to be the vi variable
   } else {                             # if not, see if 'vi' is in the object and assume that is the vi variable
      if (!is.element("vi", names(x)))
         stop(mstyle$stop("Cannot determine name of the 'vi' variable."))
      vi.name <- "vi"
   }

   #########################################################################

   k <- nrow(x)

   mf <- match.call()

   out     <- .getx("out",     mf=mf, data=data, checknumeric=TRUE)
   ci.lb   <- .getx("ci.lb",   mf=mf, data=data, checknumeric=TRUE)
   ci.ub   <- .getx("ci.ub",   mf=mf, data=data, checknumeric=TRUE)
   zval    <- .getx("zval",    mf=mf, data=data, checknumeric=TRUE)
   pval    <- .getx("pval",    mf=mf, data=data, checknumeric=TRUE)
   level   <- .getx("level",   mf=mf, data=data, checknumeric=TRUE)
   include <- .getx("include", mf=mf, data=data)

   if (is.null(level))
      level <- 95

   if (is.null(out))
      out <- rep(NA, k)
   if (is.null(zval))
      zval <- rep(NA, k)
   if (is.null(pval))
      pval <- rep(NA, k)
   if (is.null(ci.lb))
      ci.lb <- rep(NA, k)
   if (is.null(ci.ub))
      ci.ub <- rep(NA, k)

   ### if include is NULL, set to TRUE vector

   if (is.null(include))
      include <- rep(TRUE, k)

   ### turn numeric include vector into logical vector

   include <- .chksubset(include, k, stoponk0=FALSE)

   ### set inputs to NA for rows not to be included

   out[!include]   <- NA
   ci.lb[!include] <- NA
   ci.ub[!include] <- NA
   zval[!include]  <- NA
   pval[!include]  <- NA

   ### check p-values

   if (any(pval < 0, na.rm=TRUE) || any(pval > 1, na.rm=TRUE))
      stop(mstyle$stop("One or more p-values are < 0 or > 1."))

   ### if level is a single value, expand to the appropriate length

   if (length(level) == 1L)
      level <- rep(level, k)

   if (length(level) != k)
      stop(mstyle$stop(paste0("Length of the 'level' argument (", length(level), ") does not correspond to the size of the dataset (", k, ").")))

   if (!.equal.length(out, ci.lb, ci.ub, zval, pval, level))
      stop(mstyle$stop("Supplied data vectors are not all of the same length."))

   level <- .level(level, allow.vector=TRUE)
   crit <- qnorm(level/2, lower.tail=FALSE)

   ### apply transformation function if one has been specified

   if (is.function(transf)) {
      out   <- sapply(out,   transf)
      ci.lb <- sapply(ci.lb, transf)
      ci.ub <- sapply(ci.ub, transf)
   }

   ### replace missing x$yi values

   if (verbose) {
      k.repl <- sum(is.na(x[[yi.name]]) & !is.na(out))
      if (k.repl > 0L)
         message(mstyle$message("Replacing ", k.repl, " missing 'yi' value", ifelse(k.repl > 1, "s", ""), " via 'out'."))
   }

   x[[yi.name]] <- replmiss(x[[yi.name]], out)

   ### convert Wald-type CIs to sampling variances

   vi <- ((ci.ub-ci.lb)/(2*crit))^2

   ### check if yi is about halfway between CI bounds

   if (checkci) {

      # |-------------+-------------|
      # lb           yi            ub
      #               |---| (ub+lb)/2
      #
      # if the difference is more than 10% of the CI range, then flag this row

      diffs <- abs((ci.ub+ci.lb)/2 - x[[yi.name]]) / (ci.ub - ci.lb)
      #x$diffs <- diffs
      diffslarge <- diffs > cifac
      diffslarge[!is.na(x[[vi.name]])] <- NA # when x$vi is not missing, ignore diffslarge

      if (any(diffslarge, na.rm=TRUE)) {
         diffslarge <- which(diffslarge)
         if (length(diffslarge) > 5) {
            diffslarge <- paste0(paste0(head(diffslarge, 5), collapse=", "), ", ...")
         } else {
            diffslarge <- paste0(diffslarge, collapse=", ")
         }
         warning(mstyle$warning("The observed outcome does not appear to be halfway between '(ci.lb, ci.ub)' in row(s): ", diffslarge), call.=FALSE)
      }

   }

   ### replace missing x$vi values

   if (verbose) {
      k.repl <- sum(is.na(x[[vi.name]]) & !is.na(vi))
      if (k.repl > 0L)
         message(mstyle$message("Replacing ", k.repl, " missing 'vi' value", ifelse(k.repl > 1, "s", ""), " via '(ci.lb, ci.ub)'."))
   }

   x[[vi.name]] <- replmiss(x[[vi.name]], vi)

   ### convert two-sided p-values to Wald-type test statistics and replace missing zval values

   zval <- replmiss(zval, qnorm(pval/2, lower.tail=FALSE))

   ### convert Wald-type test statistics to sampling variances and replace missing x$vi values

   if (verbose) {
      k.repl <- sum(is.na(x[[vi.name]]) & !is.na(x[[yi.name]]) & !is.na(zval))
      if (k.repl > 0L)
         message(mstyle$message("Replacing ", k.repl, " missing 'vi' value", ifelse(k.repl > 1, "s", ""), " via 'zval' or 'pval'."))
   }

   x[[vi.name]] <- replmiss(x[[vi.name]], (x[[yi.name]] / zval)^2)

   return(x)

   #########################################################################

}
