############################################################################

### function to set default 'btt' value(s) or check specified 'btt' values

.set.btt <- function(btt, p, int.incl, Xnames, fixed=FALSE) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (missing(btt) || is.null(btt)) {

      if (p > 1L) {                       ### if the model matrix has more than one column
         if (int.incl) {
            btt <- seq.int(from=2, to=p)     ### and the model has an intercept term, test all coefficients except the intercept
         } else {
            btt <- seq_len(p)                ### and the model does not have an intercept term, test all coefficients
         }
      } else {
         btt <- 1L                        ### if the model matrix has a single column, test that single coefficient
      }

   } else {

      if (is.character(btt)) {

         btt <- grep(btt, Xnames, fixed=fixed)

         if (length(btt) == 0L)
            stop(mstyle$stop("Cannot identify coefficient(s) corresponding to the specified 'btt' string."), call.=FALSE)

      } else {

         ### round, take unique values, sort, and turn into integer(s)
         btt <- as.integer(sort(unique(round(btt))))

         ### check for mix of positive and negative values
         if (any(btt < 0) && any(btt > 0))
            stop(mstyle$stop("Cannot mix positive and negative 'btt' values."), call.=FALSE)

         ### keep/remove from 1:p vector as specified
         btt <- seq_len(p)[btt]

         ### (1:5)[5:6] yields c(5, NA) so remove NAs if this happens
         btt <- btt[!is.na(btt)]

         ### make sure that at least one valid value is left
         if (length(btt) == 0L)
            stop(mstyle$stop("Non-existent coefficients specified via 'btt'."), call.=FALSE)

      }

   }

   return(btt)

}

### function to format 'btt' value(s) for printing

.format.btt <- function(btt) {

   sav <- c()

   if (length(btt) > 1L) {

      btt <- sort(btt)

      while (length(btt) > 0L) {

         x <- rle(diff(btt))

         if (x$values[1] == 1 && length(x$values) != 0L) {
            sav <- c(sav, c(btt[1], ":", btt[x$lengths[1] + 1]))
            btt <- btt[-c(1:(x$lengths[1] + 1))]
            #sav <- c(sav, ", ") # this adds a space between multiple a:b sets
            sav <- c(sav, ",")
         } else {
            sav <- c(sav, btt[1], ",")
            btt <- btt[-1]
         }

      }

      sav <- paste0(sav[-length(sav)], collapse="")

   } else {

      sav <- paste0(btt)

   }

   return(sav)

}

############################################################################

### pairwise sorting of the elements of two vectors

.psort <- function(x,y) {

   ### t(apply(xy, 1, sort)) would be okay, but problematic if there are NAs;
   ### either they are removed completely (na.last=NA) or they are always put
   ### first/last (na.last=FALSE/TRUE); but we just want to leave the NAs in
   ### their position!

   if (is.null(x) || length(x) == 0L) ### need to catch this
      return(NULL)

   if (missing(y)) {
      if (is.matrix(x)) {
         xy <- x
      } else {
         xy <- rbind(x) ### in case x is just a vector
      }
   } else {
      xy <- cbind(x,y)
   }

   n <- nrow(xy)

   for (i in seq_len(n)) {
      if (anyNA(xy[i,]))
         next
      xy[i,] <- sort(xy[i,])
   }

   colnames(xy) <- NULL

   return(xy)

}

############################################################################

### function to take the square root of a vector of numbers, giving NA for negative numbers (without a warning)

.sqrt <- function(x)
   sapply(x, function(x) if (is.na(x) || x < 0) NA_real_ else sqrt(x))

### function to obtain the trace of a matrix

.tr <- function(X)
   return(sum(diag(X)))

### function to check if a matrix is square

.is.square <- function(X)
   NROW(X) == NCOL(X)

### use NROW/NCOL to better deal with scalars; compare:
### (V <- list(matrix(1, nrow=2, ncol=2), 3, c(1,4), cbind(c(2,1)))); sapply(V, function(x) nrow(x) == ncol(x)); sapply(V, function(x) NROW(x) == NCOL(x))

### function to test whether a vector is all equal to 1s (e.g., to find intercept(s) in a model matrix)

.is.intercept <- function(x, eps=1e-08)
   return(all(abs(x - 1) < eps))

### function to test whether a vector is a dummy variable (i.e., consists of only 0s and 1s)

.is.dummy <- function(x, eps=1e-08)
   return(all(abs(x) < eps | abs(x - 1) < eps))
   #return(all(sapply(x, identical, 0) | sapply(x, identical, 1)))

### function to test whether something is a vector (in the sense of being atomic, not a matrix, and not NULL)

.is.vector <- function(x)
   is.atomic(x) && !is.matrix(x) && !is.null(x)

### function to test if a string is an integer and to return the integer if so (otherwise return NA)

.is.stringint <- function(x) {
   is.int <- grepl("^[0-9]+L?$", x)
   if (is.int) {
      x <- sub("L", "", x, fixed=TRUE)
      x <- as.integer(x)
   } else {
      x <- NA
   }
   return(x)
}

### function to test if x is a matrix and that also covers Matrix objects

.is.matrix <- function(x)
   is.matrix(x) || inherits(x, "Matrix")

### function to test if x is numeric but also allow a (vector of) NA

.is.numeric <- function(x) {
   if (all(is.na(x)))
      return(TRUE)
   is.numeric(x)
}

############################################################################

### function to format p-values (no longer used; use fmtp() instead)
### if showeq=FALSE, c(.001, .00001) becomes c("0.0010", "<.0001")
### if showeq=TRUE,  c(.001, .00001) becomes c("=0.0010", "<.0001")
### if add0=FALSE, "<.0001"; if add0=TRUE, "<0.0001"

.pval <- function(p, digits=4, showeq=FALSE, sep="", add0=FALSE) {

   digits  <- max(digits, 1)
   cutoff  <- paste(c(".", rep(0,digits-1),1), collapse="")
   ncutoff <- as.numeric(cutoff)

   ifelse(is.na(p), paste0(ifelse(showeq, "=", ""), sep, "NA"),
                    ifelse(p >= ncutoff, paste0(ifelse(showeq, "=", ""), sep, formatC(p, digits=digits, format="f")),
                                         paste0("<", sep, ifelse(add0, "0", ""), cutoff)))

}

### function to format/round values in general (no longer used; use fmtx() instead)

.fcf <- function(x, digits) {

   if (all(is.na(x))) { # since formatC(NA, format="f", digits=2) fails
      rep("NA", length(x))
   } else {
      trimws(formatC(x, format="f", digits=digits))
   }

}

### function to handle 'level' argument

.level <- function(level, allow.vector=FALSE) {

   if (!allow.vector && length(level) != 1L) {
      mstyle <- .get.mstyle("crayon" %in% .packages())
      stop(mstyle$stop("Argument 'level' must specify a single value."), call.=FALSE)
   }

   if (!is.numeric(level)) {
      mstyle <- .get.mstyle("crayon" %in% .packages())
      stop(mstyle$stop("The 'level' argument must be numeric."), call.=FALSE)
   }

   ifelse(level == 0, 1, ifelse(level >= 1, (100-level)/100, ifelse(level > .5, 1-level, level)))

}

############################################################################

### function to print a named (character) vector right aligned with
### a gap of two spaces between adjacent values and no padding

.print.vector <- function(x, minfoot=NA, print.gap=2) {

   empty.last.colname <- colnames(x)[length(colnames(x))] == ""

   if (is.null(names(x)))
      names(x) <- seq_along(x)

   gap <- paste0(rep(" ", print.gap), collapse="")

   len.n   <- nchar(names(x))
   len.x   <- nchar(x, keepNA=FALSE)
   len.max <- pmax(len.n, len.x)
   #format  <- sapply(len.max, function(x) paste("%", x, "s", sep=""))

   #row.n <- paste(sprintf(format, names(x)), collapse=gap) # sprintf("%3s", "\u00b9") isn't right
   #row.x <- paste(sprintf(format, x), collapse=gap)

   #f <- function(x, n)
   #   paste0(paste0(rep(" ", n-nchar(x)), collapse=""), x, collapse="")
   #row.n <- paste(mapply(f, names(x), len.max), collapse=gap)
   #row.x <- paste(mapply(f, unname(x), len.max), collapse=gap)

   if (is.na(minfoot)) {
      row.n <- paste(mapply(formatC, names(x), width=len.max), collapse=gap) # formatC("\u00b9", width=3) works
      row.x <- paste(mapply(formatC, x, width=len.max), collapse=gap)
   } else {
      row.n <- mapply(formatC, names(x), width=len.max)
      row.n[minfoot] <- paste0(" ", row.n[minfoot])
      row.n <- paste(row.n, collapse=gap)
      row.x <- mapply(formatC, x, width=len.max)
      if (empty.last.colname) {
         row.x[length(row.x)] <- paste0(" ", row.x[length(row.x)])
      } else {
         row.x[length(row.x)] <- paste0(row.x[length(row.x)], " ")
      }
      row.x <- paste(row.x, collapse=gap)
   }

   cat(row.n, "\n", row.x, "\n", sep="")

}

.addfootsym <- function(x, cols, footsym) {
   nc <- length(cols)
   if (length(footsym) == 1L)
      footsym <- rep(footsym, nc)
   if (length(footsym) != nc)
      stop(paste0("Length of 'cols' not the same as length of 'footsym' in .addfootsym()."), call.=FALSE)
   for (i in seq_along(cols)) {
      colnames(x)[cols[i]] <- paste0(colnames(x)[cols[i]], footsym[i])
      x[[cols[i]]] <- paste0(x[[cols[i]]], " ")
   }
   return(x)
}

############################################################################

.space <- function(x=TRUE) {
   no.rmspace <- !exists(".rmspace")
   if (no.rmspace && x)
      cat("\n")
   if (!no.rmspace && !x)
      cat("\n")
}

.get.footsym <- function() {

   if (exists(".footsym")) {
      fs <- get(".footsym")
   } else {
      fs <- c("\u00b9", "1)", "\u00b2", "2)", "\u00b3", "3)")
   }

   return(fs)

}

# .footsym <- c("\u00b9", "\u00b9\u207e", "\u00b2", "\u00b2\u207e", "\u00b3", "\u00b3\u207e")

############################################################################

### function that prints the model fitting time

.print.time <- function(x) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   hours   <- floor(x/60/60)
   minutes <- floor(x/60) - hours*60
   seconds <- round(x - minutes*60 - hours*60*60, ifelse(x > 60, 0, 2))

   cat("\n")
   cat(mstyle$message(paste("Processing time:", hours, ifelse(hours == 0 || hours > 1, "hours,", "hour,"), minutes, ifelse(minutes == 0 || minutes > 1, "minutes,", "minute,"), seconds, ifelse(x < 60 || seconds == 0 || seconds > 1, "seconds", "second"))))
   cat("\n")

}

############################################################################

### function like make.unique(), but starts at .1 for the first instance
### of a repeated element

.make.unique <- function(x) {

   if (is.null(x))
      return(NULL)

   x <- as.character(x)
   ux <- unique(x)

   for (i in seq_along(ux)) {
      xiTF <- x == ux[i]
      xi <- x[xiTF]
      if (length(xi) == 1L)
         next
      x[xiTF] <- paste(xi, seq_along(xi), sep=".")
   }

   return(x)

}

############################################################################

### function to check if extra/superfluous arguments are specified via ...

.chkdots <- function(ddd, okargs) {

   for (i in seq_along(okargs))
      ddd[okargs[i]] <- NULL

   if (length(ddd) > 0L) {
      mstyle <- .get.mstyle("crayon" %in% .packages())
      warning(mstyle$warning(paste0("Extra argument", ifelse(length(ddd) > 1L, "s ", " "), "(", paste0("'", names(ddd), "'", collapse=", "), ") disregarded.")), call.=FALSE)
   }

}

############################################################################

.getx <- function(x, mf, data, enclos=sys.frame(sys.parent(n=2)), checknull=TRUE, checknumeric=FALSE, default) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   mf.getx <- match.call()
   dname <- deparse1(mf.getx[[match("data", names(mf.getx))]])
   dname <- deparse1(mf[[match(dname, names(mf))]])

   mf.x <- mf[[match(x, names(mf))]]

   if (!is.null(dname) && dname %in% names(data) && grepl("$", deparse1(mf.x), fixed=TRUE) || grepl("[[", deparse1(mf.x), fixed=TRUE))
      data <- NULL

   out <- try(eval(mf.x, data, enclos), silent=TRUE) # NULL if x was not specified

   if (inherits(out, "try-error") || is.function(out))
      stop(mstyle$stop(paste0("Cannot find the object/variable ('", deparse(mf.x), "') specified for the '", x, "' argument.")), call.=FALSE)

   # note: is.function() check catches case where 'vi' is the utils::vi() function and other shenanigans

   # check if x is actually one of the elements in the call

   spec <- x %in% names(mf)

   # out could be NULL if it is not a specified argument; if so, apply default if there is one

   if (is.null(out) && !spec && !missing(default))
      out <- default

   if (checknull) {

      # when using something like fun(dat$blah) and blah doesn't exist in dat, then get NULL

      if (spec && is.null(out)) {
         mf.txt <- deparse(mf.x)
         if (mf.txt == "NULL") {
            mf.txt <- " "
         } else {
            mf.txt <- paste0(" ('", mf.txt, "') ")
         }
         stop(mstyle$stop(paste0(deparse(mf)[1], ":\nThe object/variable", mf.txt, "specified for the '", x, "' argument is NULL.")), call.=FALSE)
      }

   }

   if (checknumeric && !is.null(out) && !is.list(out) && !.is.numeric(out[1])) # using [1] so is.numeric(Matrix(1:3)[1]) works
      stop(mstyle$stop(paste0("The object/variable specified for the '", x, "' argument is not numeric.")), call.=FALSE)

   return(out)

}

.getfromenv <- function(what, element, envir=.metafor, default=NULL) {

   x <- try(get(what, envir=envir, inherits=FALSE), silent=TRUE)
   if (inherits(x, "try-error")) {
      return(default)
   } else {
      if (missing(element)) {
         return(x)
      } else {
         x <- x[[element]]
         if (is.null(x)) {
            return(default)
         } else {
            return(x)
         }
      }
   }

}

### a version of do.call() that allows for the arguments to be passed via ... (i.e., can either be a list or not) and removes NULL arguments

.do.call <- function(fun, ...) {
   if (is.list(..1) && ...length() == 1L) {
      args <- c(...)
   } else {
      args <- list(...)
   }
   args <- args[!sapply(args, is.null)]
   do.call(fun, args)
}

############################################################################

.chkclass <- function(class, must, notap, notav, type="Method") {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   obj <- as.character(match.call()[2])
   obj <- substr(obj, 7, nchar(obj)-1)

   if (!missing(must) && !is.element(must, class))
      stop(mstyle$stop(paste0("Argument '", obj, "' must be an object of class \"", must, "\".")), call.=FALSE)

   if (!missing(notap) && any(is.element(notap, class)))
      stop(mstyle$stop(paste0(type, " not applicable to objects of class \"", class[1], "\".")), call.=FALSE)
      #stop(mstyle$stop(paste0("Method not applicable to objects of class \"", paste0(class, collapse=", "), "\".")), call.=FALSE)

   if (!missing(notav) && any(is.element(notav, class)))
      stop(mstyle$stop(paste0(type, " not available for objects of class \"", class[1], "\".")), call.=FALSE)
      #stop(mstyle$stop(paste0("Method not available for objects of class \"", paste0(class, collapse=", "), "\".")), call.=FALSE)

}

############################################################################

.chkviarg <- function(x) {

   runvicheck <- .getfromenv("runvicheck", default=TRUE)

   if (runvicheck) {

      x <- deparse(x)
      xl <- tolower(x)

      ok <- TRUE

      # starts with 'se' or 'std'
      if (any(grepl("^se", xl)))
         ok <- FALSE
      if (any(grepl("^std", xl)))
         ok <- FALSE
      # ends with 'se' or 'std'
      if (any(grepl("se$", xl)))
         ok <- FALSE
      if (any(grepl("std$", xl)))
         ok <- FALSE
      # catch cases where vi=<data frame>$se and vi=<data frame>$std
      if (any(grepl("^[[:alpha:]][[:alnum:]_.]*\\$se", xl)))
         ok <- FALSE
      if (any(grepl("^[[:alpha:]][[:alnum:]_.]*\\$std", xl)))
         ok <- FALSE

      # but if ^, *, or ( appears, don't issue a warning
      if (any(grepl("^", xl, fixed=TRUE)))
         ok <- TRUE
      if (any(grepl("*", xl, fixed=TRUE)))
         ok <- TRUE
      if (any(grepl("(", xl, fixed=TRUE)))
         ok <- TRUE

      if (!ok) {
         mstyle <- .get.mstyle("crayon" %in% .packages())
         warning(mstyle$warning(paste0("The 'vi' argument is for specifying the sampling variances,\nbut '", x, "' sounds like this variable may contain standard\nerrors (maybe use 'sei=", x, "' instead?).")), call.=FALSE)
         try(assign("runvicheck", FALSE, envir=.metafor), silent=TRUE)
      }

   }

}

############################################################################

### check that the lengths of all non-zero length elements given via ... are equal to each other

.equal.length <- function(...) {

   ddd <- list(...)
   ks <- lengths(ddd)   # get the length of each element in ddd
   if (all(ks == 0L)) { # if all elements have length 0 (are NULL), return TRUE
      return(TRUE)
   } else {
      ks <- ks[ks > 0L]  # keep the non-zero lengths
      return(length(unique(ks)) == 1L) # check that they are all identical
   }

}

### check that all elements given via ... are not of length 0 (are not NULL)

.all.specified <- function(...) {

   ddd <- list(...)
   #all(!sapply(ddd, is.null))
   not0 <- lengths(ddd) != 0L
   all(not0)

}

############################################################################

### set axis label (for forest, funnel, and labbe functions)

.setlab <- function(measure, transf.char, atransf.char, gentype, short=FALSE) {

   if (gentype == 1)
      lab <- "Observed Outcome"
   if (gentype == 2)
      lab <- "Overall Estimate" # for forest.cumul.rma() function
   if (gentype == 3)
      lab <- "Estimate"         # for header

   #########################################################################

   if (!is.null(measure)) {

      ######################################################################
      if (is.element(measure, c("RR","MPRR"))) {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "Log[RR]", "Log Risk Ratio")
         } else {
            lab <- ifelse(short, lab, "Transformed Log Risk Ratio")
            funlist <- lapply(list(exp, transf.exp.int), deparse)
            if (any(sapply(funlist, identical, atransf.char)))
               lab <- ifelse(short, "Risk Ratio", "Risk Ratio (log scale)")
            if (any(sapply(funlist, identical, transf.char)))
               lab <- ifelse(short, "Risk Ratio", "Risk Ratio")
         }
      }
      if (is.element(measure, c("OR","PETO","D2OR","D2ORN","D2ORL","MPOR","MPORC","MPPETO","MPORM"))) {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "Log[OR]", "Log Odds Ratio")
         } else {
            lab <- ifelse(short, lab, "Transformed Log Odds Ratio")
            funlist <- lapply(list(exp, transf.exp.int), deparse)
            if (any(sapply(funlist, identical, atransf.char)))
               lab <- ifelse(short, "Odds Ratio", "Odds Ratio (log scale)")
            if (any(sapply(funlist, identical, transf.char)))
               lab <- ifelse(short, "Odds Ratio", "Odds Ratio")
         }
      }
      if (is.element(measure, c("RD","MPRD"))) {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "Risk Difference", "Risk Difference")
         } else {
            lab <- ifelse(short, lab, "Transformed Risk Difference")
         }
      }
      if (measure == "AS") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "Arcsine RD", "Arcsine Transformed Risk Difference")
         } else {
            lab <- ifelse(short, lab, "Transformed Arcsine Transformed Risk Difference")
         }
      }
      if (measure == "PHI") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "Phi", "Phi Coefficient")
         } else {
            lab <- ifelse(short, lab, "Transformed Phi Coefficient")
         }
      }
      if (measure == "ZPHI") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, expression('Fisher\'s ' * z[phi]), "Fisher's z Transformed Phi Coefficient")
         } else {
            lab <- ifelse(short, lab, "Transformed Fisher's z Transformed Phi Coefficient")
            funlist <- lapply(list(transf.ztor, transf.ztor.int, tanh), deparse)
            if (any(sapply(funlist, identical, atransf.char)))
               lab <- ifelse(short, "Phi", "Phi Coefficient")
            if (any(sapply(funlist, identical, transf.char)))
               lab <- ifelse(short, "Phi", "Phi Coefficient")
         }
      }
      if (measure == "YUQ") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "Yule's Q", "Yule's Q")
         } else {
            lab <- ifelse(short, lab, "Transformed Yule's Q")
         }
      }
      if (measure == "YUY") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "Yule's Y", "Yule's Y")
         } else {
            lab <- ifelse(short, lab, "Transformed Yule's Y")
         }
      }
      ######################################################################
      if (measure == "IRR") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "Log[IRR]", "Log Incidence Rate Ratio")
         } else {
            lab <- ifelse(short, lab, "Transformed Log Incidence Rate Ratio")
            funlist <- lapply(list(exp, transf.exp.int), deparse)
            if (any(sapply(funlist, identical, atransf.char)))
               lab <- ifelse(short, "Rate Ratio", "Incidence Rate Ratio (log scale)")
            if (any(sapply(funlist, identical, transf.char)))
               lab <- ifelse(short, "Rate Ratio", "Incidence Rate Ratio")
         }
      }
      if (measure == "IRD") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "IRD", "Incidence Rate Difference")
         } else {
            lab <- ifelse(short, lab, "Transformed Incidence Rate Difference")
         }
      }
      if (measure == "IRSD") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "IRSD", "Square Root Transformed Incidence Rate Difference")
         } else {
            lab <- ifelse(short, lab, "Transformed Square Root Transformed Incidence Rate Difference")
         }
      }
      ######################################################################
      if (measure == "MD") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "MD", "Mean Difference")
         } else {
            lab <- ifelse(short, lab, "Transformed Mean Difference")
         }
      }
      if (is.element(measure, c("SMD","SMDH","SMD1","SMD1H","PBIT","OR2D","OR2DN","OR2DL"))) {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "SMD", "Standardized Mean Difference")
         } else {
            lab <- ifelse(short, lab, "Transformed Standardized Mean Difference")
         }
      }
      if (measure == "ROM") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "Log[RoM]", "Log Ratio of Means")
         } else {
            lab <- ifelse(short, lab, "Transformed Log Ratio of Means")
            funlist <- lapply(list(exp, transf.exp.int), deparse)
            if (any(sapply(funlist, identical, atransf.char)))
               lab <- ifelse(short, "Ratio of Means", "Ratio of Means (log scale)")
            if (any(sapply(funlist, identical, transf.char)))
               lab <- ifelse(short, "Ratio of Means", "Ratio of Means")
         }
      }
      if (measure == "RPB") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "Correlation", "Point-Biserial Correlation Coefficient")
         } else {
            lab <- ifelse(short, lab, "Transformed Point-Biserial Correlation Coefficient")
         }
      }
      if (measure == "ZPB") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, expression('Fisher\'s ' * z[phi]), "Fisher's z Transformed Point-Biserial Correlation Coefficient")
         } else {
            lab <- ifelse(short, lab, "Transformed Fisher's z Transformed Point-Biserial Correlation Coefficient")
            funlist <- lapply(list(transf.ztor, transf.ztor.int, tanh), deparse)
            if (any(sapply(funlist, identical, atransf.char)))
               lab <- ifelse(short, "Correlation", "Point-Biserial Correlation Coefficient")
            if (any(sapply(funlist, identical, transf.char)))
               lab <- ifelse(short, "Correlation", "Point-Biserial Correlation Coefficient")
         }
      }
      if (measure == "CVR") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "Log[CVR]", "Log Coefficient of Variation Ratio")
         } else {
            lab <- ifelse(short, lab, "Transformed Log Coefficient of Variation Ratio")
            funlist <- lapply(list(exp, transf.exp.int), deparse)
            if (any(sapply(funlist, identical, atransf.char)))
               lab <- ifelse(short, "CVR", "Coefficient of Variation Ratio (log scale)")
            if (any(sapply(funlist, identical, transf.char)))
               lab <- ifelse(short, "CVR", "Coefficient of Variation Ratio")
         }
      }
      if (measure == "VR") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "Log[VR]", "Log Variability Ratio")
         } else {
            lab <- ifelse(short, lab, "Transformed Log Variability Ratio")
            funlist <- lapply(list(exp, transf.exp.int), deparse)
            if (any(sapply(funlist, identical, atransf.char)))
               lab <- ifelse(short, "VR", "Variability Ratio (log scale)")
            if (any(sapply(funlist, identical, transf.char)))
               lab <- ifelse(short, "VR", "Variability Ratio")
         }
      }
      ######################################################################
      if (is.element(measure, c("COR","UCOR","RTET","RBIS"))) {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "Correlation", "Correlation Coefficient")
         } else {
            lab <- ifelse(short, lab, "Transformed Correlation Coefficient")
         }
      }
      if (is.element(measure, c("ZCOR","ZTET","ZBIS"))) {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, expression('Fisher\'s ' * z[r]), "Fisher's z Transformed Correlation Coefficient")
         } else {
            lab <- ifelse(short, lab, "Transformed Fisher's z Transformed Correlation Coefficient")
            funlist <- lapply(list(transf.ztor, transf.ztor.int, tanh), deparse)
            if (any(sapply(funlist, identical, atransf.char)))
               lab <- ifelse(short, "Correlation", "Correlation Coefficient")
            if (any(sapply(funlist, identical, transf.char)))
               lab <- ifelse(short, "Correlation", "Correlation Coefficient")
         }
      }
      ######################################################################
      if (measure == "PCOR") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "Correlation", "Partial Correlation Coefficient")
         } else {
            lab <- ifelse(short, lab, "Transformed Partial Correlation Coefficient")
         }
      }
      if (measure == "ZPCOR") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, expression('Fisher\'s ' * z[r]), "Fisher's z Transformed Partial Correlation Coefficient")
         } else {
            lab <- ifelse(short, lab, "Transformed Fisher's z Transformed Partial Correlation Coefficient")
            funlist <- lapply(list(transf.ztor, transf.ztor.int, tanh), deparse)
            if (any(sapply(funlist, identical, atransf.char)))
               lab <- ifelse(short, "Correlation", "Partial Correlation Coefficient")
            if (any(sapply(funlist, identical, transf.char)))
               lab <- ifelse(short, "Correlation", "Partial Correlation Coefficient")
         }
      }
      if (measure == "SPCOR") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "Correlation", "Semi-Partial Correlation Coefficient")
         } else {
            lab <- ifelse(short, lab, "Transformed Semi-Partial Correlation Coefficient")
         }
      }
      if (measure == "ZSPCOR") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, expression('Fisher\'s ' * z[r]), "Fisher's z Transformed Semi-Partial Correlation Coefficient")
         } else {
            lab <- ifelse(short, lab, "Transformed Fisher's z Transformed Semi-Partial Correlation Coefficient")
            funlist <- lapply(list(transf.ztor, transf.ztor.int, tanh), deparse)
            if (any(sapply(funlist, identical, atransf.char)))
               lab <- ifelse(short, "Correlation", "Semi-Partial Correlation Coefficient")
            if (any(sapply(funlist, identical, transf.char)))
               lab <- ifelse(short, "Correlation", "Semi-Partial Correlation Coefficient")
         }
      }
      ######################################################################
      if (measure == "PR") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "Proportion", "Proportion")
         } else {
            lab <- ifelse(short, lab, "Transformed Proportion")
         }
      }
      if (measure == "PLN") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "Log[Pr]", "Log Proportion")
         } else {
            lab <- ifelse(short, lab, "Transformed Log Proportion")
            funlist <- lapply(list(exp, transf.exp.int), deparse)
            if (any(sapply(funlist, identical, atransf.char)))
               lab <- ifelse(short, "Proportion", "Proportion (log scale)")
            if (any(sapply(funlist, identical, transf.char)))
               lab <- ifelse(short, "Proportion", "Proportion")
         }
      }
      if (measure == "PLO") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "Log[Odds]", "Log Odds")
         } else {
            lab <- ifelse(short, lab, "Transformed Log Odds")
            funlist <- lapply(list(transf.ilogit, transf.ilogit.int, plogis), deparse)
            if (any(sapply(funlist, identical, atransf.char)))
               lab <- ifelse(short, "Proportion", "Proportion (logit scale)")
            if (any(sapply(funlist, identical, transf.char)))
               lab <- ifelse(short, "Proportion", "Proportion")
            funlist <- lapply(list(exp, transf.exp.int), deparse)
            if (any(sapply(funlist, identical, atransf.char)))
               lab <- ifelse(short, "Odds", "Odds (log scale)")
            if (any(sapply(funlist, identical, transf.char)))
               lab <- ifelse(short, "Odds", "Odds")
         }
      }
      if (measure == "PAS") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, expression(arcsin(sqrt(p))), "Arcsine Transformed Proportion")
         } else {
            lab <- ifelse(short, lab, "Transformed Arcsine Transformed Proportion")
            funlist <- lapply(list(transf.iarcsin), deparse)
            if (any(sapply(funlist, identical, atransf.char)))
               lab <- ifelse(short, "Proportion", "Proportion (arcsine scale)")
            if (any(sapply(funlist, identical, transf.char)))
               lab <- ifelse(short, "Proportion", "Proportion")
         }
      }
      if (measure == "PFT") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "PFT", "Double Arcsine Transformed Proportion")
         } else {
            lab <- ifelse(short, lab, "Transformed Double Arcsine Transformed Proportion")
            funlist <- lapply(list(transf.ipft.hm), deparse)
            if (any(sapply(funlist, identical, atransf.char)))
               lab <- ifelse(short, "Proportion", "Proportion")
            if (any(sapply(funlist, identical, transf.char)))
               lab <- ifelse(short, "Proportion", "Proportion")
         }
      }
      ######################################################################
      if (measure == "IR") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "Rate", "Incidence Rate")
         } else {
            lab <- ifelse(short, lab, "Transformed Incidence Rate")
         }
      }
      if (measure == "IRLN") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "Log[IR]", "Log Incidence Rate")
         } else {
            lab <- ifelse(short, lab, "Transformed Log Incidence Rate")
            funlist <- lapply(list(exp, transf.exp.int), deparse)
            if (any(sapply(funlist, identical, atransf.char)))
               lab <- ifelse(short, "Rate", "Incidence Rate (log scale)")
            if (any(sapply(funlist, identical, transf.char)))
               lab <- ifelse(short, "Rate", "Incidence Rate")
         }
      }
      if (measure == "IRS") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "Sqrt[IR]", "Square Root Transformed Incidence Rate")
         } else {
            lab <- ifelse(short, lab, "Transformed Square Root Transformed Incidence Rate")
            funlist <- lapply(list(transf.isqrt, atransf.char), deparse)
            if (any(sapply(funlist, identical, atransf.char)))
               lab <- ifelse(short, "Rate", "Incidence Rate (square root scale)")
            if (any(sapply(funlist, identical, transf.char)))
               lab <- ifelse(short, "Rate", "Incidence Rate")
         }
      }
      if (measure == "IRFT") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "IRFT", "Freeman-Tukey Transformed Incidence Rate")
         } else {
            lab <- ifelse(short, lab, "Transformed Freeman-Tukey Transformed Incidence Rate")
         }
      }
      ######################################################################
      if (measure == "MN") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "Mean", "Mean")
         } else {
            lab <- ifelse(short, lab, "Transformed Mean")
         }
      }
      if (measure == "MNLN") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "Log[Mean]", "Log Mean")
         } else {
            lab <- ifelse(short, lab, "Transformed Log Mean")
            funlist <- lapply(list(exp, transf.exp.int), deparse)
            if (any(sapply(funlist, identical, atransf.char)))
               lab <- ifelse(short, "Mean", "Mean (log scale)")
            if (any(sapply(funlist, identical, transf.char)))
               lab <- ifelse(short, "Mean", "Mean")
         }
      }
      if (measure == "SMN") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "SM", "Standardized Mean")
         } else {
            lab <- ifelse(short, lab, "Transformed Standardized Mean")
         }
      }
      if (measure == "CVLN") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "Log[CV]", "Log Coefficient of Variation")
         } else {
            lab <- ifelse(short, lab, "Transformed Log Coefficient of Variation")
            funlist <- lapply(list(exp, transf.exp.int), deparse)
            if (any(sapply(funlist, identical, atransf.char)))
               lab <- ifelse(short, "CV", "Coefficient of Variation (log scale)")
            if (any(sapply(funlist, identical, transf.char)))
               lab <- ifelse(short, "CV", "Coefficient of Variation")
         }
      }
      if (measure == "SDLN") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "Log[SD]", "Log Standard Deviation")
         } else {
            lab <- ifelse(short, lab, "Transformed Log Standard Deviation")
            funlist <- lapply(list(exp, transf.exp.int), deparse)
            if (any(sapply(funlist, identical, atransf.char)))
               lab <- ifelse(short, "SD", "Standard Deviation (log scale)")
            if (any(sapply(funlist, identical, transf.char)))
               lab <- ifelse(short, "SD", "Standard Deviation")
         }
      }
      ######################################################################
      if (measure == "MC") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "Mean Change", "Mean Change")
         } else {
            lab <- ifelse(short, lab, "Transformed Mean Change")
         }
      }
      if (is.element(measure, c("SMCC","SMCR","SMCRH"))) {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "SMC", "Standardized Mean Change")
         } else {
            lab <- ifelse(short, lab, "Transformed Standardized Mean Change")
         }
      }
      if (measure == "ROMC") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "Log[RoM]", "Log Ratio of Means")
         } else {
            lab <- ifelse(short, lab, "Transformed Log Ratio of Means")
            funlist <- lapply(list(exp, transf.exp.int), deparse)
            if (any(sapply(funlist, identical, atransf.char)))
               lab <- ifelse(short, "Ratio of Means", "Ratio of Means (log scale)")
            if (any(sapply(funlist, identical, transf.char)))
               lab <- ifelse(short, "Ratio of Means", "Ratio of Means")
         }
      }
      if (measure == "CVRC") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "Log[CVR]", "Log Coefficient of Variation Ratio")
         } else {
            lab <- ifelse(short, lab, "Transformed Log Coefficient of Variation Ratio")
            funlist <- lapply(list(exp, transf.exp.int), deparse)
            if (any(sapply(funlist, identical, atransf.char)))
               lab <- ifelse(short, "CVR", "Coefficient of Variation Ratio (log scale)")
            if (any(sapply(funlist, identical, transf.char)))
               lab <- ifelse(short, "CVR", "Coefficient of Variation Ratio")
         }
      }
      if (measure == "VRC") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "Log[VR]", "Log Variability Ratio")
         } else {
            lab <- ifelse(short, lab, "Transformed Log Variability Ratio")
            funlist <- lapply(list(exp, transf.exp.int), deparse)
            if (any(sapply(funlist, identical, atransf.char)))
               lab <- ifelse(short, "VR", "Variability Ratio (log scale)")
            if (any(sapply(funlist, identical, transf.char)))
               lab <- ifelse(short, "VR", "Variability Ratio")
         }
      }
      ######################################################################
      if (measure == "ARAW") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "Alpha", "Cronbach's alpha")
         } else {
            lab <- ifelse(short, lab, "Transformed Cronbach's alpha")
         }
      }
      if (measure == "AHW") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, expression('Alpha'[HW]), "Transformed Cronbach's alpha")
         } else {
            lab <- ifelse(short, lab, "Transformed Cronbach's alpha")
            funlist <- lapply(list(transf.iahw), deparse)
            if (any(sapply(funlist, identical, atransf.char)))
               lab <- ifelse(short, "Alpha", "Cronbach's alpha")
            if (any(sapply(funlist, identical, transf.char)))
               lab <- ifelse(short, "Alpha", "Cronbach's alpha")
         }
      }
      if (measure == "ABT") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, expression('Alpha'[B]), "Transformed Cronbach's alpha")
         } else {
            lab <- ifelse(short, lab, "Transformed Cronbach's alpha")
            funlist <- lapply(list(transf.iabt), deparse)
            if (any(sapply(funlist, identical, atransf.char)))
               lab <- ifelse(short, "Alpha", "Cronbach's alpha")
            if (any(sapply(funlist, identical, transf.char)))
               lab <- ifelse(short, "Alpha", "Cronbach's alpha")
         }
      }
      ######################################################################
      if (measure == "REH") {
         if (identical(transf.char, "FALSE") && identical(atransf.char, "FALSE")) {
            lab <- ifelse(short, "Log[REH]", "Log Relative Excess Heterozygosity")
         } else {
            lab <- ifelse(short, lab, "Transformed Log Relative Excess Heterozygosity")
            funlist <- lapply(list(exp, transf.exp.int), deparse)
            if (any(sapply(funlist, identical, atransf.char)))
               lab <- ifelse(short, "REH", "Relative Excess Heterozygosity (log scale)")
            if (any(sapply(funlist, identical, transf.char)))
               lab <- ifelse(short, "REH", "Relative Excess Heterozygosity")
         }
      }
      ######################################################################

   }

   return(lab)

}

############################################################################

### stuff related to colored/styled output

.get.mstyle <- function(withcrayon) {

   if (withcrayon) {

      if (exists(".mstyle")) {
         .mstyle <- get(".mstyle")
         if (!is.list(.mstyle))
            .mstyle <- list(.mstyle)
      } else {
         .mstyle <- list()
      }

      if (is.null(.mstyle$section)) {
         section <- crayon::bold
      } else {
         section <- .mstyle$section
      }
      if (is.null(.mstyle$header)) {
         header <- crayon::underline
      } else {
         header <- .mstyle$header
      }
      if (is.null(.mstyle$body1)) {
         body1 <- crayon::reset
      } else {
         body1 <- .mstyle$body1
      }
      if (is.null(.mstyle$body2)) {
         body2 <- crayon::reset
      } else {
         body2 <- .mstyle$body2
      }
      if (is.null(.mstyle$na)) {
         na <- crayon::reset
      } else {
         na <- .mstyle$na
      }
      if (is.null(.mstyle$text)) {
         text <- crayon::reset
      } else {
         text <- .mstyle$text
      }
      if (is.null(.mstyle$result)) {
         result <- crayon::reset
      } else {
         result <- .mstyle$result
      }
      if (is.null(.mstyle$stop)) {
         stop <- crayon::combine_styles(crayon::red, crayon::bold)
      } else {
         stop <- .mstyle$stop
      }
      if (is.null(.mstyle$warning)) {
         warning <- crayon::yellow
      } else {
         warning <- .mstyle$warning
      }
      if (is.null(.mstyle$message)) {
         message <- crayon::green
      } else {
         message <- .mstyle$message
      }
      if (is.null(.mstyle$verbose)) {
         verbose <- crayon::cyan
      } else {
         verbose <- .mstyle$verbose
      }
      if (is.null(.mstyle$legend)) {
         legend <- crayon::silver
         #legend <- crayon::make_style("gray90")
      } else {
         legend <- .mstyle$legend
      }

   } else {

      tmp <- function(...) paste0(...)
      section <- tmp
      header  <- tmp
      body1   <- tmp
      body2   <- tmp
      na      <- tmp
      text    <- tmp
      result  <- tmp
      stop    <- tmp
      warning <- tmp
      message <- tmp
      verbose <- tmp
      legend  <- tmp

   }

   return(list(section=section, header=header, body1=body1, body2=body2, na=na, text=text, result=result, stop=stop, warning=warning, message=message, verbose=verbose, legend=legend))

}

.print.output <- function(x, mstyle) {

   if (missing(mstyle)) {
      for (i in seq_along(x)) {
         cat(x[i], "\n")
      }
   } else {
      for (i in seq_along(x)) {
         cat(mstyle(x[i]), "\n")
      }
   }

}

.is.even <- function(x) x %% 2 == 0

.print.table <- function(x, mstyle) {

   is.header <- !grepl(" [-0-9]", x)
   #is.header <- !grepl("^\\s*[0-9]", x)
   has.header <- any(is.header)

   for (i in seq_along(x)) {
      if (is.header[i]) {
         #x[i] <- trimws(x[i], which="right")
         x[i] <- mstyle$header(x[i])
      } else {
         x[i] <- gsub("NA", mstyle$na("NA"), x[i], fixed=TRUE)
         if (.is.even(i-has.header)) {
            x[i] <- mstyle$body2(x[i])
         } else {
            x[i] <- mstyle$body1(x[i])
         }
      }
      cat(x[i], "\n")
   }

}

#.set.mstyle.1 <- str2lang(".mstyle <- list(section=make_style(\"gray90\")$bold, header=make_style(\"skyblue1\")$bold$underline, body=make_style(\"skyblue2\"), text=make_style(\"slateblue3\"), result=make_style(\"slateblue1\"))")
#eval(metafor:::.set.mstyle.1)

############################################################################

.set.digits <- function(digits, dmiss) {

   res <- c(est=4, se=4, test=4, pval=4, ci=4, var=4, sevar=4, fit=4, het=4)

   if (exists(".digits")) {
      .digits <- get(".digits")
      if (is.null(names(.digits)) && length(.digits) == 1L) {
         # if .digits is a single unnamed scalar, set all digit values to that value
         res <- c(est=.digits, se=.digits, test=.digits, pval=.digits, ci=.digits, var=.digits, sevar=.digits, fit=.digits, het=.digits)
      } else if (any(names(.digits) != "") && any(names(.digits) == "")) {
         # if .digits has (at least) one unnamed element, use it to set all unnamed elements to that digits value
         pos <- pmatch(names(.digits), names(res))
         res[c(na.omit(pos))] <- .digits[!is.na(pos)]
         otherval <- .digits[names(.digits) == ""][1]
         res[(1:9)[-c(na.omit(pos))]] <- otherval
      } else {
         pos <- pmatch(names(.digits), names(res))
         res[c(na.omit(pos))] <- .digits[!is.na(pos)]
      }
   }

   if (!dmiss) {
      if (is.null(names(digits))) {
         res <- c(est=digits[[1]], se=digits[[1]], test=digits[[1]], pval=digits[[1]], ci=digits[[1]], var=digits[[1]], sevar=digits[[1]], fit=digits[[1]], het=digits[[1]])
      } else {
         pos <- pmatch(names(digits), names(res))
         res[c(na.omit(pos))] <- digits[!is.na(pos)]
      }
   }

   ### p-values are always given to at least 2 digits
   if (res["pval"] <= 1)
      res["pval"] <- 2

   res

}

.get.digits <- function(digits, xdigits, dmiss) {

   res <- xdigits

   if (exists(".digits")) {
      .digits <- get(".digits")
      pos <- pmatch(names(.digits), names(res))
      res[c(na.omit(pos))] <- .digits[!is.na(pos)]
   }

   if (!dmiss) {
      if (is.null(names(digits))) {
         res <- c(est=digits[[1]], se=digits[[1]], test=digits[[1]], pval=digits[[1]], ci=digits[[1]], var=digits[[1]], sevar=digits[[1]], fit=digits[[1]], het=digits[[1]])
      } else {
         pos <- pmatch(names(digits), names(res))
         res[c(na.omit(pos))] <- digits[!is.na(pos)]
      }
   }

   ### so we can still print objects created with older metafor versions (where xdigit is just an unnamed scalar)
   if (length(res) == 1L && is.null(names(res)))
      res <- c(est=res[[1]], se=res[[1]], test=res[[1]], pval=res[[1]], ci=res[[1]], var=res[[1]], sevar=res[[1]], fit=res[[1]], het=res[[1]])

   ### p-values are always given to at least 2 digits
   if (!is.null(res["pval"]) && res["pval"] <= 1)
      res["pval"] <- 2

   res

}

############################################################################

### check if x is logical and TRUE/FALSE (NAs and NULL always evaluate as FALSE)

.isTRUE <- function(x)
   !is.null(x) && is.logical(x) && !is.na(x) && x

.isFALSE <- function(x)
   !is.null(x) && is.logical(x) && !is.na(x) && !x

# not sure anymore why I implemented these; c(isTRUE(NULL), isTRUE(NA), isFALSE(NULL), isFALSE(NA)) are all FALSE

############################################################################

### shorten a character vector so that elements remain distinguishable

.shorten <- function(x, minlen) {

   y <- x

   x <- c(na.omit(x))

   n <- length(unique(x))

   maxlen <- max(nchar(unique(x)))

   for (l in seq_len(maxlen)) {
      tab <- table(x, substr(x, 1, l))
      if (nrow(tab) == n && ncol(tab) == n && sum(tab[upper.tri(tab)]) == 0 && sum(tab[lower.tri(tab)]) == 0)
         break
   }

   if (!missing(minlen) && l < minlen) {
      if (minlen > maxlen)
         minlen <- maxlen
      l <- minlen
   }

   return(substr(y, 1, l))

}

############################################################################

### simplified version of what mvtnorm::rmvnorm() does

.mvrnorm <- function(n, mu, Sigma) {

   p <- nrow(Sigma)
   eS <- eigen(Sigma, symmetric = TRUE)
   eval <- eS$values
   evec <- eS$vectors

   Y <- matrix(rnorm(p * n), nrow = n, byrow = TRUE) %*% t(evec %*% (t(evec) * sqrt(pmax(eval, 0))))
   Y <- sweep(Y, 2, mu, "+")

   return(Y)

}

############################################################################

### check subset argument (if logical, make sure it's of the right length and set NAs to FALSE; if
### numeric, remove NAs and 0's and check that values are not beyond k)

.chksubset <- function(x, k, stoponk0=TRUE) {

   if (is.null(x)) # if x is NULL, return x (i.e., NULL)
      return(x)

   mstyle <- .get.mstyle("crayon" %in% .packages())

   argname <- deparse(substitute(x))

   if (length(x) == 0L)
      stop(mstyle$stop(paste0("Argument '", argname, "' is of length 0.")), call.=FALSE)

   if (is.logical(x)) {
      if (length(x) != k)
         stop(mstyle$stop(paste0("Length of the '", argname, "' argument (", length(x), ") is not of length k = ", k, ".")), call.=FALSE)
      #x <- x[seq_len(k)]     # keep only elements 1:k from x
      if (anyNA(x))           # if x includes any NA elements
         x[is.na(x)] <- FALSE # set NA elements to FALSE
   }

   if (is.numeric(x)) {
      if (anyNA(x))             # if x includes any NA elements
         x <- x[!is.na(x)]      # remove them
      x <- as.integer(round(x))
      x <- x[x != 0L]           # also remove any 0's
      if (any(x > 0L) && any(x < 0L))
         stop(mstyle$stop(paste0("Cannot mix positive and negative values in '", argname, "' argument.")), call.=FALSE)
      if (all(x > 0L)) {
         if (any(x > k))
            stop(mstyle$stop(paste0("Argument '", argname, "' includes values larger than k = ", k, ".")), call.=FALSE)
         x <- is.element(seq_len(k), x)
      } else {
         if (any(x < -k))
            stop(mstyle$stop(paste0("Argument '", argname, "' includes values larger than k = ", k, ".")), call.=FALSE)
         x <- !is.element(seq_len(k), abs(x))
      }
   }

   if (stoponk0 && !any(x))
      stop(mstyle$stop(paste0("Stopped because k = 0 after subsetting.")), call.=FALSE)

   return(x)

}

### get subset function that works for matrices and data frames (selecting rows by default but rows
### and columns when col=TRUE) and vectors and also checks that x is of the same length as subset

.getsubset <- function(x, subset, col=FALSE, drop=FALSE) {

   if (is.null(x) || is.null(subset)) # if x or subset is NULL, return x
      return(x)

   mstyle <- .get.mstyle("crayon" %in% .packages())

   xname <- deparse(substitute(x))

   k <- length(subset)

   if (.is.matrix(x) || is.data.frame(x)) {
      if (nrow(x) != k)
         stop(mstyle$stop(paste0("Element '", xname, "' is not of length ", k, ".")), call.=FALSE)
      if (col) {
         x <- x[subset,subset,drop=drop]
      } else {
         x <- x[subset,,drop=drop]
      }
   } else {
      if (length(x) != k)
         stop(mstyle$stop(paste0("Element '", xname, "' is not of length ", k, ".")), call.=FALSE)
      x <- x[subset]
   }

   return(x)

}

############################################################################

# function to compute a weighted mean (this one works a bit different than
# stats:::weighted.mean.default)

.wmean <- function (x, w, na.rm=FALSE) {
   if (na.rm) {
      i <- !(is.na(x) | is.na(w)) # only include x if x and w are not missing
      x <- x[i]
      w <- w[i]
   }
   sum(x*w) / sum(w)
}

############################################################################

.tes.intfun <- function(x, theta, tau, sei, H0, alternative, crit) {
   if (alternative == "two.sided")
      pow <- (pnorm(crit, mean=(x-H0)/sei, sd=1, lower.tail=FALSE) + pnorm(-crit, mean=(x-H0)/sei, sd=1, lower.tail=TRUE))
   if (alternative == "greater")
      pow <- pnorm(crit, mean=(x-H0)/sei, sd=1, lower.tail=FALSE)
   if (alternative == "less")
      pow <- pnorm(crit, mean=(x-H0)/sei, sd=1, lower.tail=TRUE)
   res <- pow * dnorm(x, theta, tau)
   return(res)
}

.tes.lim <- function(theta, yi, vi, H0, alternative, alpha, tau2, test, tes.alternative, progbar, tes.alpha, correct, rel.tol, subdivisions, tau2.lb) {
   pval <- tes(x=yi, vi=vi, H0=H0, alternative=alternative, alpha=alpha, theta=theta, tau2=tau2, test=test, tes.alternative=tes.alternative, progbar=progbar,
                       tes.alpha=tes.alpha, correct=correct, rel.tol=rel.tol, subdivisions=subdivisions, tau2.lb=tau2.lb, find.lim=FALSE)$pval
   #cat("theta = ", theta, " pval = ", pval, "\n")
   return(pval - tes.alpha)
}

############################################################################

.fsn.fisher <- function(fsnum, pi, alpha) {
   k <- length(pi)
   X2 <- -2*sum(log(c(pi, rep(0.5, fsnum))))
   return(pchisq(X2, df=2*(k+fsnum), lower.tail=FALSE) - alpha)
}

.fsn.fitre <- function(yi, vi) {

   k     <- length(yi)
   wi    <- 1/vi
   sumwi <- sum(wi)
   est   <- sum(wi*yi)/sumwi
   Q     <- sum(wi * (yi - est)^2)
   tau2  <- max(0, (Q - (k-1)) / (sumwi - sum(wi^2)/sumwi))
   wi    <- 1 / (vi + tau2)
   sumwi <- sum(wi)
   est   <- sum(wi*yi)/sumwi
   se    <- sqrt(1 / sumwi)
   zval  <- est / se
   pval  <- 2*pnorm(abs(zval), lower.tail=FALSE)

   return(list(est=est, se=se, zval=zval, pval=pval, tau2=tau2))

}

.fsn.fitnew <- function(new, yi, vi, vnew, tau2, alpha, iters) {

   new <- ceiling(new)

   mus   <- rep(NA_real_, iters)
   pvals <- rep(NA_real_, iters)

   for (j in seq_len(iters)) {
      yinew <- c(yi, rnorm(new, 0, sqrt(vnew+tau2)))
      vinew <- c(vi, rep(vnew, new))
      tmp <- .fsn.fitre(yinew, vinew)
      mus[j] <- tmp$est
      pvals[j] <- tmp$pval
   }

   return(list(mean = mean(mus), rejrate = mean(pvals <= alpha)))

}

.fsn.re <- function(fsnum, yi, vi, vnew, tau2, target, alpha, iters, verbose=FALSE) {

   fsnum <- ceiling(fsnum)
   tmp <- .fsn.fitnew(fsnum, yi, vi, vnew, tau2, alpha, iters)
   est <- tmp$mean
   diff <- est - target
   if (verbose)
      cat("fsnum =", formatC(fsnum, width=4, format="d"), "  est =", fmtx(est), "  target =", fmtx(target), "  diff =", fmtx(diff, flag=" "), "\n")
   return(diff)

}

############################################################################

.chkopt <- function(optimizer, optcontrol) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   ### set NLOPT_LN_BOBYQA as the default algorithm for nloptr optimizer
   ### and by default use a relative convergence criterion of 1e-8 on the function value

   if (optimizer == "nloptr" && !is.element("algorithm", names(optcontrol)))
      optcontrol$algorithm <- "NLOPT_LN_BOBYQA"

   if (optimizer == "nloptr" && !is.element("ftol_rel", names(optcontrol)))
      optcontrol$ftol_rel <- 1e-8

   ### for mads, set trace=FALSE and tol=1e-6 by default

   if (optimizer == "mads" && !is.element("trace", names(optcontrol)))
      optcontrol$trace <- FALSE

   if (optimizer == "mads" && !is.element("tol", names(optcontrol)))
      optcontrol$tol <- 1e-6

   ### for subplex, set reltol=1e-8 by default (the default in subplex() is .Machine$double.eps)

   if (optimizer == "subplex" && !is.element("reltol", names(optcontrol)))
      optcontrol$reltol <- 1e-8

   ### for BBoptim, set trace=FALSE by default

   if (optimizer == "BBoptim" && !is.element("trace", names(optcontrol)))
      optcontrol$trace <- FALSE

   ### for solnp, set trace=FALSE by default

   if (optimizer == "solnp" && !is.element("trace", names(optcontrol)))
      optcontrol$trace <- FALSE

   ### check that the required packages are installed

   if (is.element(optimizer, c("uobyqa","newuoa","bobyqa"))) {
      if (!requireNamespace("minqa", quietly=TRUE))
         stop(mstyle$stop("Please install the 'minqa' package to use this optimizer."), call.=FALSE)
   }

   if (is.element(optimizer, c("nloptr","ucminf","lbfgsb3c","subplex","optimParallel"))) {
      if (!requireNamespace(optimizer, quietly=TRUE))
         stop(mstyle$stop(paste0("Please install the '", optimizer, "' package to use this optimizer.")), call.=FALSE)
   }

   if (is.element(optimizer, c("hjk","nmk","mads"))) {
      if (!requireNamespace("dfoptim", quietly=TRUE))
         stop(mstyle$stop("Please install the 'dfoptim' package to use this optimizer."), call.=FALSE)
   }

   if (optimizer == "BBoptim") {
      if (!requireNamespace("BB", quietly=TRUE))
         stop(mstyle$stop("Please install the 'BB' package to use this optimizer."), call.=FALSE)
   }

   if (optimizer == "solnp") {
      if (!requireNamespace("Rsolnp", quietly=TRUE))
         stop(mstyle$stop("Please install the 'Rsolnp' package to use this optimizer."), call.=FALSE)
   }

   if (optimizer == "constrOptim.nl") {
      if (!requireNamespace("alabama", quietly=TRUE))
         stop(mstyle$stop("Please install the 'alabama' package to use this optimizer."), call.=FALSE)
   }

   if (optimizer == "Rcgmin") {
      if (!requireNamespace("Rcgmin", quietly=TRUE))
         stop(mstyle$stop("Please install the 'Rcgmin' package to use this optimizer."), call.=FALSE)
   }

   if (optimizer == "Rvmmin") {
      if (!requireNamespace("Rvmmin", quietly=TRUE))
         stop(mstyle$stop("Please install the 'Rvmmin' package to use this optimizer."), call.=FALSE)
   }

   #########################################################################

   if (is.element(optimizer, c("optim","constrOptim"))) {
      par.arg <- "par"
      ctrl.arg <- ", control=optcontrol"
   }

   if (optimizer == "nlminb") {
      par.arg <- "start"
      ctrl.arg <- ", control=optcontrol"
   }

   if (is.element(optimizer, c("uobyqa","newuoa","bobyqa"))) {
      par.arg <- "par"
      optimizer <- paste0("minqa::", optimizer) # need to use this since loading nloptr masks bobyqa() and newuoa() functions
      ctrl.arg <- ", control=optcontrol"
   }

   if (optimizer == "nloptr") {
      par.arg <- "x0"
      optimizer <- paste0("nloptr::nloptr") # need to use this due to requireNamespace()
      ctrl.arg <- ", opts=optcontrol"
   }

   if (optimizer == "nlm") {
      par.arg <- "p" # because of this, must use argument name pX for p (number of columns in X matrix)
      ctrl.arg <- paste(names(optcontrol), unlist(optcontrol), sep="=", collapse=", ")
      if (nchar(ctrl.arg) != 0L)
         ctrl.arg <- paste0(", ", ctrl.arg)
   }

   if (is.element(optimizer, c("hjk","nmk","mads"))) {
      par.arg <- "par"
      optimizer <- paste0("dfoptim::", optimizer) # need to use this so that the optimizers can be found
      ctrl.arg <- ", control=optcontrol"
   }

   if (is.element(optimizer, c("ucminf","lbfgsb3c","subplex"))) {
      par.arg <- "par"
      optimizer <- paste0(optimizer, "::", optimizer) # need to use this due to requireNamespace()
      ctrl.arg <- ", control=optcontrol"
   }

   if (optimizer == "BBoptim") {
      par.arg <- "par"
      optimizer <- "BB::BBoptim"
      ctrl.arg <- ", quiet=TRUE, control=optcontrol"
   }

   if (optimizer == "solnp") {
      par.arg <- "pars"
      optimizer <- "Rsolnp::solnp"
      ctrl.arg <- ", control=optcontrol"
   }

   if (optimizer == "constrOptim.nl") {
      par.arg <- "par"
      optimizer <- "alabama::constrOptim.nl"
      if ("control.outer" %in% names(optcontrol)) {
         # can specify 'control.outer' to be passed to constrOptim.nl(), but when using
         # the 'method' argument, must escape " or use ' for this to work; for example:
         # control=list(optimizer="constrOptim.nl", control.outer=list(method="'Nelder-Mead'"))
         control.outer <- paste0("control.outer=list(", paste(names(optcontrol$control.outer), unlist(optcontrol$control.outer), sep="=", collapse=", "), ")")
         ctrl.arg <- paste0(", control.optim=optcontrol, ", control.outer)
         optcontrol$control.outer <- NULL
      } else {
         ctrl.arg <- ", control.optim=optcontrol, control.outer=list(trace=FALSE)"
      }
   }

   if (optimizer == "Rcgmin") {
      par.arg <- "par"
      optimizer <- "Rcgmin::Rcgmin"
      #ctrl.arg <- ", gr='grnd', control=optcontrol"
      ctrl.arg <- ", control=optcontrol"
   }

   if (optimizer == "Rvmmin") {
      par.arg <- "par"
      optimizer <- "Rvmmin::Rvmmin"
      #ctrl.arg <- ", gr='grnd', control=optcontrol"
      ctrl.arg <- ", control=optcontrol"
   }

   if (optimizer == "optimParallel") {
      par.arg <- "par"
      optimizer <- "optimParallel::optimParallel"
      ctrl.arg <- ", control=optcontrol, parallel=parallel"
   }

   return(list(optimizer=optimizer, optcontrol=optcontrol, par.arg=par.arg, ctrl.arg=ctrl.arg))

}

.chkconv <- function(optimizer, opt.res, optcontrol, fun, verbose) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (optimizer == "optimParallel::optimParallel" && verbose) {
      tmp <- capture.output(print(opt.res$loginfo))
      .print.output(tmp, mstyle$verbose)
   }

   ### convergence checks

   if (inherits(opt.res, "try-error"))
      stop(mstyle$stop(paste0("Error during the optimization. Use verbose=TRUE and see help(", fun, ") for more details on the optimization routines.")), call.=FALSE)

   if (optimizer == "lbfgsb3c::lbfgsb3c" && is.null(opt.res$convergence)) # special provision for lbfgsb3c in case 'convergence' is missing
      opt.res$convergence <- -99

   if (is.element(optimizer, c("optim","constrOptim","nlminb","dfoptim::hjk","dfoptim::nmk","lbfgsb3c::lbfgsb3c","subplex::subplex","BB::BBoptim","Rsolnp::solnp","alabama::constrOptim.nl","Rcgmin::Rcgmin","Rvmmin:Rvmmin","optimParallel::optimParallel")) && opt.res$convergence != 0)
      stop(mstyle$stop(paste0("Optimizer (", optimizer, ") did not achieve convergence (convergence = ", opt.res$convergence, ").")), call.=FALSE)

   if (is.element(optimizer, c("dfoptim::mads")) && opt.res$convergence > optcontrol$tol)
      stop(mstyle$stop(paste0("Optimizer (", optimizer, ") did not achieve convergence (convergence = ", opt.res$convergence, ").")), call.=FALSE)

   if (is.element(optimizer, c("minqa::uobyqa","minqa::newuoa","minqa::bobyqa")) && opt.res$ierr != 0)
      stop(mstyle$stop(paste0("Optimizer (", optimizer, ") did not achieve convergence (ierr = ", opt.res$ierr, ").")), call.=FALSE)

   if (optimizer=="nloptr::nloptr" && !(opt.res$status >= 1 && opt.res$status <= 4))
      stop(mstyle$stop(paste0("Optimizer (", optimizer, ") did not achieve convergence (status = ", opt.res$status, ").")), call.=FALSE)

   if (optimizer=="ucminf::ucminf" && !(opt.res$convergence == 1 || opt.res$convergence == 2))
      stop(mstyle$stop(paste0("Optimizer (", optimizer, ") did not achieve convergence (convergence = ", opt.res$convergence, ").")), call.=FALSE)

   if (verbose > 2) {
      cat("\n")
      tmp <- capture.output(print(opt.res))
      .print.output(tmp, mstyle$verbose)
   }

   ### copy estimated values to 'par'

   if (optimizer=="nloptr::nloptr")
      opt.res$par <- opt.res$solution
   if (optimizer=="nlm")
      opt.res$par <- opt.res$estimate
   if (optimizer=="Rsolnp::solnp")
      opt.res$par <- opt.res$pars

   return(opt.res$par)

}

############################################################################

.coltail <- function(h, val, tail="upper", mult=1, col, border, freq, ...) {

   h$counts  <- h$counts  * mult
   h$density <- h$density * mult

   if (tail == "lower") {

      above <- which(h$breaks > val)
      if (length(above) > 0L) {
         pos <- above[1]
         h$breaks[pos] <- val
      }
      sel <- h$breaks <= val
      if (sum(sel) >= 2L) {
         h$breaks  <- h$breaks[sel]
         h$counts  <- h$counts[sel[-1]]
         h$density <- h$density[sel[-1]]
         h$mids    <- h$mids[sel[-1]]
         lines(h, col=col, border=border, freq=freq, ...)
      }

   } else {

      below <- which(h$breaks < val)
      if (length(below) > 0L) {
         pos <- below[length(below)]
         h$breaks[pos] <- val
      }
      sel <- h$breaks >= val
      if (sum(sel) >= 2L) {
         len <- length(below)
         h$breaks  <- h$breaks[sel]
         h$counts  <- h$counts[sel[-len]]
         h$density <- h$density[sel[-len]]
         h$mids    <- h$mids[sel[-len]]
         lines(h, col=col, border=border, freq=freq, ...)
      }

   }

}

############################################################################
