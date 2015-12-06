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

"[.escalc" <- function(x, ...) {

   dat <- NextMethod("[")

   ### add measure attribute back to each yi variable (if that variable is still part of dat)

   yi.names <- attr(x, "yi.names")

   if (!is.null(yi.names)) {
      for (i in 1:length(yi.names)) {
         if (!is.element(yi.names[i], names(dat))) {
            next
         } else {
            eval(parse(text=paste0("attr(dat$", yi.names[i], ", 'measure') <- attr(x$", yi.names[i], ", 'measure')")))
         }
      }
   }

   ### note: ni attribute is lost upon subsetting, but it's better that way (since it would also need
   ### to be subsetted, as otherwise it could end up having a different length than yi itself)

   ### same goes for slab attribute

   ### add var.names and out.names attributes back to object

   attr(dat, "yi.names")    <- attr(x, "yi.names")
   attr(dat, "vi.names")    <- attr(x, "vi.names")
   attr(dat, "sei.names")   <- attr(x, "sei.names")
   attr(dat, "zi.names")    <- attr(x, "zi.names")
   attr(dat, "ci.lb.names") <- attr(x, "ci.lb.names")
   attr(dat, "ci.ub.names") <- attr(x, "ci.ub.names")

   ### add digits attribute back to object

   if (!is.null(attr(x, "digits")))
      attr(dat, "digits") <- attr(x, "digits")

   ### TODO: clean up attribute elements that are no longer actually part of the object

   return(dat)

}

############################################################################

cbind.escalc <- function (..., deparse.level=1) {

   dat <- data.frame(..., check.names = FALSE)

   arguments <- list(...)

   ### for each element, extract the 'var.names' and 'out.names' attributes and add that entire set back to the object

   yi.names    <- NULL
   vi.names    <- NULL
   sei.names   <- NULL
   zi.names    <- NULL
   ci.lb.names <- NULL
   ci.ub.names <- NULL
   digits      <- NULL

   for (arg in arguments) {
      yi.names    <- c(attr(arg, "yi.names"),    yi.names)
      vi.names    <- c(attr(arg, "vi.names"),    vi.names)
      sei.names   <- c(attr(arg, "sei.names"),   sei.names)
      zi.names    <- c(attr(arg, "zi.names"),    zi.names)
      ci.lb.names <- c(attr(arg, "ci.lb.names"), ci.lb.names)
      ci.ub.names <- c(attr(arg, "ci.ub.names"), ci.ub.names)
      digits      <- c(attr(arg, "digits"), digits)
   }

   attr(dat, "yi.names")    <- unique(yi.names)
   attr(dat, "vi.names")    <- unique(vi.names)
   attr(dat, "sei.names")   <- unique(sei.names)
   attr(dat, "zi.names")    <- unique(zi.names)
   attr(dat, "ci.lb.names") <- unique(ci.lb.names)
   attr(dat, "ci.ub.names") <- unique(ci.ub.names)

   # var.name.list <- list()
   # digits.list   <- list()

   # i <- 0

   # for (arg in arguments) {
      # i <- i + 1
      # var.name.list[[i]] <- attr(arg, "var.names")
      # digits.list[[i]]   <- attr(arg, "digits")
   # }

   # if (length(var.name.list) > 0)
      # attr(dat, "var.names") <- var.name.list[[1]]

   # if (length(digits.list) > 0)
      # attr(dat, "digits") <- digits.list[[1]]

   ### add 'digits' attribute back (use the first value)

   attr(dat, "digits") <- digits[1]

   ### TODO: clean up attribute elements that are no longer actually part of the object

   class(dat) <- c("escalc", "data.frame")
   return(dat)

}

############################################################################

replmiss <- function(x, y) {

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
