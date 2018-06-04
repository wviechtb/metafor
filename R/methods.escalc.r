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
