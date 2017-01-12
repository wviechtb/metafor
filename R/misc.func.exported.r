############################################################################

"[.list.rma" <- function(x, i, ...) { ### removed j argument (see below), so can only select rows, not columns

   out <- x

   attr(out, "class") <- NULL

   slab.pos <- which(names(out) == "slab")

   if (!missing(i)) ### for X element
      out[seq_len(slab.pos-1)] <- lapply(out[seq_len(slab.pos-1)], function(r) if (class(r) == "matrix") r[i,] else r[i])

   ### catch cases where user selects values outside 1:k

   if (length(out[[1]]) == 0L)
      return(NULL)

   #out <- out[j] ### this causes all kinds of problems, so left out for now

   out$slab <- x$slab[i]

   ### slab can only contain NAs if user selects values outside 1:k

   if (anyNA(out$slab))
      return(NULL)

   out$digits <- x$digits
   out$transf <- x$transf
   out$method <- x$method

   class(out) <- "list.rma"
   return(out)

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

rbind.escalc <- function (..., deparse.level=1) {

   dat <- rbind.data.frame(..., deparse.level = deparse.level)

   allargs <- list(...)

   yi.names <- attr(dat, "yi.names")
   yi.names <- yi.names[is.element(yi.names, names(dat))]

   for (i in seq_along(yi.names)) {

      ### get 'ni' attribute from all arguments
      ni <- lapply(allargs, function(x) attr(x[[yi.names[i]]], "ni"))

      ### if none of them are missing, then combine and add back to variable
      if (all(sapply(ni, function(x) !is.null(x))))
         attr(dat[[yi.names[i]]], "ni") <- unlist(ni)

      ### get 'slab' attribute from all arguments
      slab <- lapply(allargs, function(x) attr(x[[yi.names[i]]], "slab"))

      ### if none of them are missing, then combine and add back to variable (and make sure they are unique)
      if (all(sapply(slab, function(x) !is.null(x))))
         attr(dat[[yi.names[i]]], "slab") <- make.unique(as.character(unlist(slab)))

   }

   return(dat)

}

############################################################################

replmiss <- function(x, y) {

   if (length(y) == 0L)
      y <- NA

   ### catch cases where x is of length 0

   if (length(x) == 0L)
      return(x)

   ### in case user specifies a constant to use for replacement

   if (length(y) == 1L)
      y <- rep(y, length(x))

   ### check that x and y are of the same length

   if (length(x) != length(y))
      stop("Length of 'x' and 'y' do not match.")

   #x <- ifelse(is.na(x), y, x) ### this is quite a bit slower than the following
   x[is.na(x)] <- y[is.na(x)]

   return(x)

}

############################################################################
