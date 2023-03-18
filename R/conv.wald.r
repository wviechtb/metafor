conv.wald <- function(out, ci.lb, ci.ub, zval, pval, n, data, include,
                      level=95, transf, check=TRUE, var.names, append=TRUE, replace="ifna", ...) {

   # TODO: allow t-distribution based CIs/tests (then also need dfs argument)?

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (missing(out) && missing(ci.lb) && missing(ci.ub) && missing(zval) && missing(pval))
      stop(mstyle$stop("Must specify at least some of these arguments: 'out', 'ci.lb', 'ci.ub', 'zval', 'pval'."))

   replace <- match.arg(replace, c("ifna","all"))

   ### get ... argument and check for extra/superfluous arguments

   ddd <- list(...)

   .chkdots(ddd, c("cifac"))

   cifac <- ifelse(is.null(ddd$cifac), 0.1, ddd$cifac)

   #########################################################################

   if (missing(data))
      data <- NULL

   has.data <- !is.null(data)

   if (is.null(data)) {
      data <- sys.frame(sys.parent())
   } else {
      if (!is.data.frame(data))
         data <- data.frame(data)
   }

   x <- data

   ### checks on var.names argument

   if (missing(var.names)) {

      if (inherits(x, "escalc")) {

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

      } else {

         yi.name <- "yi"
         vi.name <- "vi"

      }

   } else {

      if (length(var.names) != 2L)
         stop(mstyle$stop("Argument 'var.names' must be of length 2."))

      if (any(var.names != make.names(var.names, unique=TRUE))) {
         var.names <- make.names(var.names, unique=TRUE)
         warning(mstyle$warning(paste0("Argument 'var.names' does not contain syntactically valid variable names.\nVariable names adjusted to: var.names = c('", var.names[1], "','", var.names[2], "').")), call.=FALSE)
      }

      yi.name <- var.names[1]
      vi.name <- var.names[2]

   }

   if (missing(transf))
      transf <- FALSE

   #########################################################################

   mf <- match.call()

   out     <- .getx("out",     mf=mf, data=x, checknumeric=TRUE)
   ci.lb   <- .getx("ci.lb",   mf=mf, data=x, checknumeric=TRUE)
   ci.ub   <- .getx("ci.ub",   mf=mf, data=x, checknumeric=TRUE)
   zval    <- .getx("zval",    mf=mf, data=x, checknumeric=TRUE)
   pval    <- .getx("pval",    mf=mf, data=x, checknumeric=TRUE)
   n       <- .getx("n",       mf=mf, data=x, checknumeric=TRUE)
   level   <- .getx("level",   mf=mf, data=x, checknumeric=TRUE, default=95)
   include <- .getx("include", mf=mf, data=x)

   if (!.equal.length(out, ci.lb, ci.ub, zval, pval, n))
      stop(mstyle$stop("Supplied data vectors are not all of the same length."))

   k <- max(length(out), length(ci.lb), length(ci.ub), length(zval), length(pval), length(n))

   if (is.null(out))
      out <- rep(NA_real_, k)
   if (is.null(ci.lb))
      ci.lb <- rep(NA_real_, k)
   if (is.null(ci.ub))
      ci.ub <- rep(NA_real_, k)
   if (is.null(zval))
      zval <- rep(NA_real_, k)
   if (is.null(pval))
      pval <- rep(NA_real_, k)
   if (is.null(n))
      n <- rep(NA_real_, k)

   ### if include is NULL, set to TRUE vector

   if (is.null(include))
      include <- rep(TRUE, k)

   ### turn numeric include vector into logical vector

   include <- .chksubset(include, k, stoponk0=FALSE)

   ### set inputs to NA for rows not to be included

   out[!include]   <- NA_real_
   ci.lb[!include] <- NA_real_
   ci.ub[!include] <- NA_real_
   zval[!include]  <- NA_real_
   pval[!include]  <- NA_real_
   n[!include]     <- NA_real_

   ### check p-values

   if (any(pval < 0, na.rm=TRUE) || any(pval > 1, na.rm=TRUE))
      stop(mstyle$stop("One or more p-values are < 0 or > 1."))

   ### if level is a single value, expand to the appropriate length

   if (length(level) == 1L)
      level <- rep(level, k)

   if (length(level) != k)
      stop(mstyle$stop(paste0("Length of the 'level' argument (", length(level), ") does not correspond to the size of the dataset (", k, ").")))

   level <- .level(level, allow.vector=TRUE)
   crit <- qnorm(level/2, lower.tail=FALSE)

   ### apply transformation function if one has been specified

   if (is.function(transf)) {
      out   <- sapply(out,   transf)
      ci.lb <- sapply(ci.lb, transf)
      ci.ub <- sapply(ci.ub, transf)
   }

   ### set up data frame if 'data' was not specified

   if (!has.data) {
      x <- data.frame(rep(NA_real_, k), rep(NA_real_, k))
      names(x) <- c(yi.name, vi.name)
   }

   #########################################################################

   ### replace missing x$yi values

   if (replace=="ifna") {
      x[[yi.name]] <- replmiss(x[[yi.name]], out)
   } else {
      x[[yi.name]][!is.na(out)] <- out[!is.na(out)]
   }

   ### replace missing ni attribute values (or add 'ni' attribute if at least one value is not missing)

   if (!is.null(attributes(x[[yi.name]])$ni)) {
      attributes(x[[yi.name]])$ni <- replmiss(attributes(x[[yi.name]])$ni, n)
   } else {
      if (any(!is.na(n)))
         attr(x[[yi.name]], "ni") <- n
   }

   #########################################################################

   ### convert Wald-type CIs to sampling variances

   vi <- ((ci.ub-ci.lb)/(2*crit))^2

   ### check if yi is about halfway between CI bounds

   if (check) {

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

   ### convert two-sided p-values to Wald-type test statistics and replace missing zval values

   zval <- replmiss(zval, qnorm(pval/2, lower.tail=FALSE))

   ### convert Wald-type test statistics to sampling variances and replace missing vi values

   vi <- replmiss(vi, (x[[yi.name]] / zval)^2)

   ### note: if both (ci.lb,ci.ub) and zval/pval is available, then this favors
   ### the back-calculation based on (ci.lb,ci.ub) which seems reasonable

   ### TODO: could consider checking if the back-calculated vi's differs in this case
   ### (or if x$vi is already available)

   ### replace missing x$vi values

   if (replace=="ifna") {
      x[[vi.name]] <- replmiss(x[[vi.name]], vi)
   } else {
      x[[vi.name]][!is.na(vi)] <- vi[!is.na(vi)]
   }

   #########################################################################

   measure <- attr(x[[yi.name]], "measure")

   if (is.null(measure))
      measure <- "GEN"

   x <- escalc(measure=measure, data=x, yi=x[[yi.name]], vi=x[[vi.name]], var.names=c(yi.name,vi.name))

   if (!append)
      x <- x[,c(yi.name, vi.name)]

   return(x)

   #########################################################################

}
