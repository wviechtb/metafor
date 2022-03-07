summary.escalc <- function(object, out.names=c("sei","zi","pval","ci.lb","ci.ub"), var.names,
H0=0, append=TRUE, replace=TRUE, level=95, olim, digits, transf, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(object), must="escalc")

   x <- object

   level <- .level(level)

   crit <- qnorm(level/2, lower.tail=FALSE)

   if (length(out.names) != 5L)
      stop(mstyle$stop("Argument 'out.names' must be of length 5."))

   if (any(out.names != make.names(out.names, unique=TRUE))) {
      out.names <- make.names(out.names, unique=TRUE)
      warning(mstyle$warning(paste0("Argument 'out.names' does not contain syntactically valid variable names.\nVariable names adjusted to: out.names = c('", out.names[1], "', '", out.names[2], "', '", out.names[3], "', '", out.names[4], "', '", out.names[5], "').")), call.=FALSE)
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
            stop(mstyle$stop("Cannot determine name of the 'yi' variable."))
         yi.name <- "yi"
      }
      if (!is.null(attr(x, "vi.names"))) { ### if vi.names attributes is available
         vi.name <- attr(x, "vi.names")[1] ### take the first entry to be the vi variable
      } else {                             ### if not, see if 'vi' is in the object and assume that is the vi variable
         if (!is.element("vi", names(x)))
            stop(mstyle$stop("Cannot determine name of the 'vi' variable."))
         vi.name <- "vi"
      }

   } else {

      if (length(var.names) != 2L)
         stop(mstyle$stop("Argument 'var.names' must be of length 2."))

      if (any(var.names != make.names(var.names, unique=TRUE))) {
         var.names <- make.names(var.names, unique=TRUE)
         warning(mstyle$warning(paste0("Argument 'var.names' does not contain syntactically valid variable names.\nVariable names adjusted to: var.names = c('", var.names[1], "', '", var.names[2], "').")), call.=FALSE)
      }

      yi.name <- var.names[1]
      vi.name <- var.names[2]

   }

   yi <- x[[yi.name]]
   vi <- x[[vi.name]]

   if (is.null(yi) || is.null(vi))
      stop(mstyle$stop(paste0("Cannot find variables '", yi.name, "' and/or '", vi.name, "' in the data frame.")))

   #########################################################################

   k <- length(yi)

   if (length(H0) == 1L)
      H0 <- rep(H0, k)

   ### compute sei, zi, and lower/upper CI bounds; when applying a transformation, compute the transformed outcome and CI bounds

   sei  <- sqrt(vi)
   zi   <- c(yi - H0) / sei
   pval <- 2*pnorm(abs(zi), lower.tail=FALSE)
   if (is.function(transf)) {
      ci.lb <- mapply(transf, yi - crit * sei, ...)
      ci.ub <- mapply(transf, yi + crit * sei, ...)
      yi    <- mapply(transf, yi, ...)
      attr(x, "transf") <- TRUE
      vi   <- NULL
      sei  <- NULL
      zi   <- NULL
      pval <- NULL
   } else {
      ci.lb <- yi - crit * sei
      ci.ub <- yi + crit * sei
      attr(x, "transf") <- FALSE
   }

   ### make sure order of intervals is always increasing

   tmp <- .psort(ci.lb, ci.ub)
   ci.lb <- tmp[,1]
   ci.ub <- tmp[,2]

   ### apply observation/outcome limits if specified

   if (!missing(olim)) {
      if (length(olim) != 2L)
         stop(mstyle$stop("Argument 'olim' must be of length 2."))
      olim <- sort(olim)
      yi[yi < olim[1]] <- olim[1] # note: zi and pval are based on unconstrained yi
      yi[yi > olim[2]] <- olim[2]
      ci.lb[ci.lb < olim[1]] <- olim[1]
      ci.ub[ci.ub > olim[2]] <- olim[2]
   }

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
         dat[[out.names[3]]] <- pval  ### if variable does not exists in dat, it will be added
         dat[[out.names[4]]] <- ci.lb ### if variable does not exists in dat, it will be added
         dat[[out.names[5]]] <- ci.ub ### if variable does not exists in dat, it will be added

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

         if (is.element(out.names[3], names(dat))) { ### if pval variable is in data frame, replace NA values with newly calculated values
            is.na.pval <- is.na(dat[[out.names[3]]])
            dat[[out.names[3]]][is.na.pval] <- pval[is.na.pval]
         } else {
            dat[[out.names[3]]] <- pval             ### if pval variable does not exist in dat, just add as new variable
         }

         if (is.element(out.names[4], names(dat))) { ### if ci.lb variable is in data frame, replace NA values with newly calculated values
            is.na.ci.lb <- is.na(dat[[out.names[4]]])
            dat[[out.names[4]]][is.na.ci.lb] <- ci.lb[is.na.ci.lb]
         } else {
            dat[[out.names[4]]] <- ci.lb             ### if ci.lb variable does not exist in dat, just add as new variable
         }

         if (is.element(out.names[5], names(dat))) { ### if ci.ub variable is in data frame, replace NA values with newly calculated values
            is.na.ci.ub <- is.na(dat[[out.names[5]]])
            dat[[out.names[5]]][is.na.ci.ub] <- ci.ub[is.na.ci.ub]
         } else {
            dat[[out.names[5]]] <- ci.ub             ### if ci.ub variable does not exist in dat, just add as new variable
         }

      }

   } else {

      ### if user does not want to append

      if (is.function(transf)) {
         dat <- data.frame(yi, ci.lb, ci.ub)
         names(dat) <- c(yi.name, out.names[4:5])
      } else {
         dat <- data.frame(yi, vi, sei, zi, pval, ci.lb, ci.ub)
         names(dat) <- c(yi.name, vi.name, out.names)
      }
   }

   ### update existing digits attribute if digits is specified

   if (!missing(digits)) {
      attr(dat, "digits") <- .get.digits(digits=digits, xdigits=attr(x, "digits"), dmiss=FALSE)
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

   ### add 'sei.names', 'zi.names', 'pval.names', 'ci.lb.names', and 'ci.ub.names' to the first position of the corresponding attributes
   ### note: if "xyz" is not an attribute of the object, attr(object, "xyz") returns NULL, so this works fine

   attr(dat, "sei.names")   <- unique(c(out.names[1], attr(object, "sei.names")))
   attr(dat, "zi.names")    <- unique(c(out.names[2], attr(object, "zi.names")))
   attr(dat, "pval.names")  <- unique(c(out.names[3], attr(object, "pval.names")))
   attr(dat, "ci.lb.names") <- unique(c(out.names[4], attr(object, "ci.lb.names")))
   attr(dat, "ci.ub.names") <- unique(c(out.names[5], attr(object, "ci.ub.names")))

   ### TODO: clean up attribute elements that are no longer actually part of the object

   class(dat) <- c("escalc", "data.frame")
   return(dat)

}
