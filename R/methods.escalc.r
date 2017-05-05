############################################################################

print.escalc <- function(x, digits, ...) {

   if (!inherits(x, "escalc"))
      stop("Argument 'x' must be an object of class \"escalc\".")

   attr(x, "class") <- NULL

   if (missing(digits))
      digits <- attr(x, "digits")

   if (is.null(digits))
      digits <- 4

   ### get positions of the variable names in the object
   ### note: if the object no longer contains a particular variable, match() returns NA;
   ### use na.omit(), so that length() is then zero (as needed for if statements below)

   yi.pos    <- na.omit(match(attr(x, "yi.names"),    names(x)))
   vi.pos    <- na.omit(match(attr(x, "vi.names"),    names(x)))
   sei.pos   <- na.omit(match(attr(x, "sei.names"),   names(x)))
   zi.pos    <- na.omit(match(attr(x, "zi.names"),    names(x)))
   ci.lb.pos <- na.omit(match(attr(x, "ci.lb.names"), names(x)))
   ci.ub.pos <- na.omit(match(attr(x, "ci.ub.names"), names(x)))

   ### get rownames attribute so we can back-assign it

   rnames <- attr(x, "row.names")
   x <- data.frame(x)
   rownames(x) <- rnames

   ### round variables according to the digits argument

   if (length(yi.pos) > 0)
      x[yi.pos] <- apply(x[yi.pos], 2, formatC, digits=digits, format="f")

   if (length(vi.pos) > 0)
      x[vi.pos] <- apply(x[vi.pos], 2, formatC, digits=digits, format="f")

   if (length(sei.pos) > 0)
      x[sei.pos] <- apply(x[sei.pos], 2, formatC, digits=digits, format="f")

   if (length(zi.pos) > 0)
      x[zi.pos] <- apply(x[zi.pos], 2, formatC, digits=digits, format="f")

   if (length(ci.lb.pos) > 0)
      x[ci.lb.pos] <- apply(x[ci.lb.pos], 2, formatC, digits=digits, format="f")

   if (length(ci.ub.pos) > 0)
      x[ci.ub.pos] <- apply(x[ci.ub.pos], 2, formatC, digits=digits, format="f")

   print(x, ...)

}

############################################################################

summary.escalc <- function(object, out.names=c("sei","zi","ci.lb","ci.ub"), var.names,
H0=0, append=TRUE, replace=TRUE, level=95, digits, transf, ...) {

   if (!inherits(object, "escalc"))
      stop("Argument 'object' must be an object of class \"escalc\".")

   x <- object

   level <- ifelse(level > 1, (100-level)/100, ifelse(level > .5, 1-level, level))

   crit <- qnorm(level/2, lower.tail=FALSE)

   if (length(out.names) != 4)
      stop("Argument 'out.names' must be of length 4.")

   if (any(out.names != make.names(out.names, unique=TRUE))) {
      out.names <- make.names(out.names, unique=TRUE)
      warning(paste0("Argument 'out.names' does not contain syntactically valid variable names.\n  Variable names adjusted to: out.names = c('", out.names[1], "', '", out.names[2], "', '", out.names[3], "', '", out.names[4], "')."))
   }

   if (missing(transf))
      transf <- FALSE

   #########################################################################

   ### figure out names of yi and vi variables (if possible) and extract the values (if possible)

   if (missing(var.names)) {               ### if var.names not specified, take from object if possible

      if (!is.null(attr(x, "yi.names"))) { ### if yi.names attributes is available
         yi.name <- attr(x, "yi.names")[1] ### take the first entry to be the yi variable
      } else {                             ### if not, see if 'yi' is in the object and assume that is the yi variable
         if (!is.element("yi", names(x)))
            stop("Cannot determine name of the 'yi' variable.")
         yi.name <- "yi"
      }
      if (!is.null(attr(x, "vi.names"))) { ### if vi.names attributes is available
         vi.name <- attr(x, "vi.names")[1] ### take the first entry to be the vi variable
      } else {                             ### if not, see if 'vi' is in the object and assume that is the vi variable
         if (!is.element("vi", names(x)))
            stop("Cannot determine name of the 'vi' variable.")
         vi.name <- "vi"
      }

   } else {

      if (length(var.names) != 2)
         stop("Argument 'var.names' must be of length 2.")

      if (any(var.names != make.names(var.names, unique=TRUE))) {
         var.names <- make.names(var.names, unique=TRUE)
         warning(paste0("Argument 'var.names' does not contain syntactically valid variable names.\n  Variable names adjusted to: var.names = c('", var.names[1], "', '", var.names[2], "')."))
      }

      yi.name <- var.names[1]
      vi.name <- var.names[2]

   }

   yi <- x[[yi.name]]
   vi <- x[[vi.name]]

   if (is.null(yi) || is.null(vi))
      stop(paste0("Cannot find variables '", yi.name, "' and/or '", vi.name, "' in the data frame."))

   #########################################################################

   k <- length(yi)

   if (length(H0) == 1)
      H0 <- rep(H0, k)

   ### compute sei, zi, and lower/upper CI bounds; when applying a transformation, compute the transformed outcome and CI bounds

   sei <- sqrt(vi)
   zi  <- c(yi - H0) / sei
   if (is.function(transf)) {
      ci.lb <- mapply(transf, yi - crit * sei, ...)
      ci.ub <- mapply(transf, yi + crit * sei, ...)
      yi    <- mapply(transf, yi, ...)
      attr(x, "transf") <- TRUE
   } else {
      ci.lb <- yi - crit * sei
      ci.ub <- yi + crit * sei
      attr(x, "transf") <- FALSE
   }

   ### make sure order of intervals is always increasing

   tmp <- .psort(ci.lb, ci.ub)
   ci.lb <- tmp[,1]
   ci.ub <- tmp[,2]

   x[[yi.name]] <- yi
   x[[vi.name]] <- vi

   #return(cbind(yi, vi, sei, zi, ci.lb, ci.ub))

   ### put together dataset

   if (append) {

      ### if user wants to append

      dat <- data.frame(x)

      if (replace) {

         ### and wants to replace all values

         dat[[out.names[1]]] <- sei   ### if variable does not exists in dat, it will be added
         dat[[out.names[2]]] <- zi    ### if variable does not exists in dat, it will be added
         dat[[out.names[3]]] <- ci.lb ### if variable does not exists in dat, it will be added
         dat[[out.names[4]]] <- ci.ub ### if variable does not exists in dat, it will be added

      } else {

         ### and only wants to replace any NA values

         if (is.element(out.names[1], names(dat))) { ### if sei variable is in data frame, replace NA values with newly calculated values
            is.na.sei <- is.na(dat[[out.names[1]]])
            dat[[out.names[1]]][is.na.sei] <- sei[is.na.sei]
         } else {
            dat[[out.names[1]]] <- sei               ### if sei variable does not exist in dat, just add as new variable
         }

         if (is.element(out.names[2], names(dat))) { ### if zi variable is in data frame, replace NA values with newly calculated values
            is.na.zi <- is.na(dat[[out.names[2]]])
            dat[[out.names[2]]][is.na.zi] <- zi[is.na.zi]
         } else {
            dat[[out.names[2]]] <- zi                ### if zi variable does not exist in dat, just add as new variable
         }

         if (is.element(out.names[3], names(dat))) { ### if ci.lb variable is in data frame, replace NA values with newly calculated values
            is.na.ci.lb <- is.na(dat[[out.names[3]]])
            dat[[out.names[3]]][is.na.ci.lb] <- ci.lb[is.na.ci.lb]
         } else {
            dat[[out.names[3]]] <- ci.lb             ### if ci.lb variable does not exist in dat, just add as new variable
         }

         if (is.element(out.names[4], names(dat))) { ### if ci.ub variable is in data frame, replace NA values with newly calculated values
            is.na.ci.ub <- is.na(dat[[out.names[4]]])
            dat[[out.names[4]]][is.na.ci.ub] <- ci.ub[is.na.ci.ub]
         } else {
            dat[[out.names[4]]] <- ci.ub             ### if ci.ub variable does not exist in dat, just add as new variable
         }

      }

   } else {

      ### if user does not want to append

      dat <- data.frame(yi, vi, sei, zi, ci.lb, ci.ub)
      names(dat) <- c(yi.name, vi.name, out.names)

   }

   ### update existing digits attribute if digits is specified

   if (!missing(digits)) {
      attr(dat, "digits") <- digits
   } else {
      attr(dat, "digits") <- attr(x, "digits")
   }

   if (is.null(attr(dat, "digits"))) ### in case x no longer has a 'digits' attribute
      attr(dat, "digits") <- 4

   ### update existing var.names attribute if var.names is specified
   ### and make sure all other yi.names and vi.names are added back in

   if (!missing(var.names)) {
      attr(dat, "yi.names") <- unique(c(var.names[1], attr(object, "yi.names")))
   } else {
      attr(dat, "yi.names") <- unique(c(yi.name, attr(object, "yi.names")))
   }

   if (!missing(var.names)) {
      attr(dat, "vi.names") <- unique(c(var.names[2], attr(object, "vi.names")))
   } else {
      attr(dat, "vi.names") <- unique(c(vi.name, attr(object, "vi.names")))
   }

   ### add 'sei.names', 'zi.names', 'ci.lb.names', and 'ci.ub.names' to the first position of the corresponding attributes
   ### note: if "xyz" is not an attribute of the object, attr(object, "xyz") returns NULL, so this works fine

   attr(dat, "sei.names")   <- unique(c(out.names[1], attr(object, "sei.names")))
   attr(dat, "zi.names")    <- unique(c(out.names[2], attr(object, "zi.names")))
   attr(dat, "ci.lb.names") <- unique(c(out.names[3], attr(object, "ci.lb.names")))
   attr(dat, "ci.ub.names") <- unique(c(out.names[4], attr(object, "ci.ub.names")))

   ### TODO: clean up attribute elements that are no longer actually part of the object

   class(dat) <- c("escalc", "data.frame")
   return(dat)

}

############################################################################

"[.escalc" <- function(x, i, ...) {

   dat <- NextMethod("[")

   ### add measure attribute back to each yi variable (if that variable is still part of dat)

   yi.names <- attr(x, "yi.names")
   yi.names <- yi.names[is.element(yi.names, names(dat))]

   for (l in seq_along(yi.names)) {

      #eval(parse(text=paste0("attr(dat$", yi.names[l], ", 'measure') <- attr(x$", yi.names[l], ", 'measure')")))
      attr(dat[[yi.names[l]]], "measure") <- attr(x[[yi.names[l]]], "measure")

      ### if selecting rows, also subset ni and slab attributes and add them back to each yi variable

      if (!missing(i)) {
         attr(dat[[yi.names[l]]], "ni")   <- attr(x[[yi.names[l]]], "ni")[i]
         attr(dat[[yi.names[l]]], "slab") <- attr(x[[yi.names[l]]], "slab")[i]
      }

   }

   ### add var.names and out.names attributes back to object (but only if they exist and only keep variables still in the dataset)

   all.names <- c("yi.names", "vi.names", "sei.names", "zi.names", "ci.lb.names", "ci.ub.names")

   for (l in seq_along(all.names)) {
      if (any(is.element(attr(x, all.names[l]), names(dat)))) ### check if any of the variables still exist in the dataset
         attr(dat, all.names[l]) <- attr(x, all.names[l])[is.element(attr(x, all.names[l]), names(dat))]
   }

   ### add digits attribute back to object (but not to vectors)

   if (!is.null(attr(x, "digits")) && !is.null(dim(dat)))
      attr(dat, "digits") <- attr(x, "digits")

   return(dat)

}

############################################################################

cbind.escalc <- function (..., deparse.level=1) {

   dat <- data.frame(..., check.names = FALSE)

   allargs <- list(...)

   ### for each element, extract the 'var.names' and 'out.names' attributes and add entire set back to the object

   yi.names    <- NULL
   vi.names    <- NULL
   sei.names   <- NULL
   zi.names    <- NULL
   ci.lb.names <- NULL
   ci.ub.names <- NULL
   digits      <- NULL

   for (arg in allargs) {
      yi.names    <- c(attr(arg, "yi.names"),    yi.names)
      vi.names    <- c(attr(arg, "vi.names"),    vi.names)
      sei.names   <- c(attr(arg, "sei.names"),   sei.names)
      zi.names    <- c(attr(arg, "zi.names"),    zi.names)
      ci.lb.names <- c(attr(arg, "ci.lb.names"), ci.lb.names)
      ci.ub.names <- c(attr(arg, "ci.ub.names"), ci.ub.names)
      digits      <- c(attr(arg, "digits"),      digits)
   }

   ### but only keep unique variable names

   attr(dat, "yi.names")    <- unique(yi.names)
   attr(dat, "vi.names")    <- unique(vi.names)
   attr(dat, "sei.names")   <- unique(sei.names)
   attr(dat, "zi.names")    <- unique(zi.names)
   attr(dat, "ci.lb.names") <- unique(ci.lb.names)
   attr(dat, "ci.ub.names") <- unique(ci.ub.names)

   ### add 'digits' attribute back (use the first value)

   attr(dat, "digits") <- digits[1]

   class(dat) <- c("escalc", "data.frame")
   return(dat)

}

############################################################################

rbind.escalc <- function (..., deparse.level=1) {

   dat <- rbind.data.frame(..., deparse.level = deparse.level)

   allargs <- list(...)

   yi.names <- attr(dat, "yi.names")
   yi.names <- yi.names[is.element(yi.names, names(dat))]

   for (i in seq_along(yi.names)) {

      yi.pos <- which(names(allargs[[1]]) == yi.names[i])

      ### get 'ni' attribute from all arguments
      ni <- lapply(allargs, function(x) attr(x[[yi.pos]], "ni"))

      ### if none of them are missing, then combine and add back to variable
      ### otherwise remove 'ni' attribute, since it won't be of the right length
      if (all(sapply(ni, function(x) !is.null(x)))) {
         attr(dat[[yi.pos]], "ni") <- unlist(ni)
      } else {
         attr(dat[[yi.pos]], "ni") <- NULL
      }

      ### get 'slab' attribute from all arguments
      slab <- lapply(allargs, function(x) attr(x[[yi.pos]], "slab"))

      ### if none of them are missing, then combine and add back to variable (and make sure they are unique)
      ### otherwise remove 'slab' attribute, since it won't be of the right length
      if (all(sapply(slab, function(x) !is.null(x)))) {
         attr(dat[[yi.pos]], "slab") <- .make.unique(unlist(slab))
      } else {
         attr(dat[[yi.pos]], "slab") <- NULL
      }

   }

   return(dat)

}

############################################################################
