conv.delta <- function(yi, vi, data, include, transf, var.names, append=TRUE, replace="ifna", ...) {

   # TODO: be able to specify ni to be added to ni attribute

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (missing(yi) || missing(vi))
      stop(mstyle$stop("Must specify 'yi' and 'vi' arguments."))

   if (missing(transf))
      stop(mstyle$stop("Must specify 'transf' argument."))

   replace <- match.arg(replace, c("ifna","all"))

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

   #########################################################################

   mf <- match.call()

   yi      <- .getx("yi",      mf=mf, data=x, checknumeric=TRUE)
   vi      <- .getx("vi",      mf=mf, data=x, checknumeric=TRUE)
   include <- .getx("include", mf=mf, data=x)

   ### check length of yi and vi

   if (length(yi) != length(vi))
      stop(mstyle$stop("Length of 'yi' and 'vi' is not the same."))

   ### check 'vi' argument for potential misuse

   .chkviarg(mf$vi)

   k <- length(yi)

   ### if include is NULL, set to TRUE vector

   if (is.null(include))
      include <- rep(TRUE, k)

   ### turn numeric include vector into logical vector

   include <- .chksubset(include, k, stoponk0=FALSE)

   ### set inputs to NA for rows not to be included

   yi[!include] <- NA_real_
   vi[!include] <- NA_real_

   ### get ... args

   args <- names(sapply(mf[-1], deparse))
   rmargs <- c("yi", "vi", "data", "include", "transf", "var.names", "append", "replace")
   dotargs <- args[!args %in% rmargs]

   dotarglist <- list()

   for (i in seq_along(dotargs)) {
      dotarglist[[i]] <- .getx(dotargs[i], mf=mf, data=x, checknumeric=TRUE)
      if (length(dotarglist[[i]]) == 1L)
         dotarglist[[i]] <- rep(dotarglist[[i]], k)
      names(dotarglist)[i] <- dotargs[i]
   }

   #print(dotarglist)

   argmatch <- pmatch(names(dotarglist), table=c("func","method","side"), duplicates.ok=TRUE)

   if (!all(is.na(argmatch)))
      stop(mstyle$stop("One or more arguments in ... (partially) match an argument from numDeriv::grad()."))

   #########################################################################

   #ddd <- list(c(yi), ...)
   #yi.t  <- unlist(.mapply(FUN=transf, dots=ddd, MoreArgs=NULL))
   #deriv <- unlist(.mapply(FUN=.compgrad, dots=ddd, MoreArgs=list(func=transf)))
   #vi.t  <- vi * deriv^2

   #dat <- data.frame(yi=yi.t, vi=vi.t)
   #return(dat)

   yi.t  <- rep(NA_real_, k)
   vi.t  <- rep(NA_real_, k)
   deriv <- rep(NA_real_, k)

   for (i in 1:k) {

      args <- c(yi[i], as.list(sapply(dotarglist, `[`, i)))
      #print(args)

      tmp <- try(suppressWarnings(do.call(transf, args)), silent=TRUE)
      #tmp <- try(do.call(transf, args), silent=FALSE)

      if (inherits(tmp, "try-error")) {
         yi.t[i] <- NA
      } else {
         yi.t[i] <- tmp
      }

      args <- c(args, func=transf)
      #print(args)

      tmp <- try(suppressWarnings(do.call(numDeriv::grad, args)), silent=TRUE)
      #tmp <- try(do.call(numDeriv::grad, args))

      if (inherits(tmp, "try-error")) {
         vi.t[i] <- NA
      } else {
         vi.t[i] <- vi[i] * tmp^2
      }

      #tmp <- try(suppressWarnings(numDeriv::grad(func=transf, yi[i])), silent=TRUE)

      #if (inherits(tmp, "try-error")) {
      #   deriv[i] <- NA
      #} else {
      #   deriv[i] <- tmp
      #}

      #vi.t[i] <- vi[i] * deriv[i]^2

   }

   #########################################################################

   ### set up data frame if 'data' was not specified

   if (!has.data) {
      x <- data.frame(rep(NA_real_, k), rep(NA_real_, k))
      names(x) <- c(yi.name, vi.name)
   }

   ### replace missing x$yi values

   if (replace=="ifna") {
      x[[yi.name]] <- replmiss(x[[yi.name]], yi.t)
   } else {
      x[[yi.name]][!is.na(yi.t)] <- yi.t[!is.na(yi.t)]
   }

   ### replace missing x$vi values

   if (replace=="ifna") {
      x[[vi.name]] <- replmiss(x[[vi.name]], vi.t)
   } else {
      x[[vi.name]][!is.na(vi.t)] <- vi.t[!is.na(vi.t)]
   }

   #if (!inherits(x, "escalc"))
   x <- escalc(data=x, yi=x[[yi.name]], vi=x[[vi.name]], var.names=c(yi.name,vi.name))

   if (!append)
      x <- x[,c(yi.name, vi.name)]

   return(x)

}
