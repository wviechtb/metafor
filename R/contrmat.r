contrmat <- function(data, grp1, grp2, last, shorten=FALSE, minlen=2, append=TRUE) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!is.data.frame(data))
      data <- data.frame(data)

   ### get variable names

   varnames <- names(data)

   ### number of variables

   nvars <- length(varnames)

   ############################################################################

   ### checks on 'grp1' argument

   if (length(grp1) != 1L)
      stop(mstyle$stop("Argument 'grp1' must of length 1."))

   if (!(is.character(grp1) | is.numeric(grp1)))
      stop(mstyle$stop("Argument 'grp1' must either be a character string or a number."))

   if (is.character(grp1)) {

      grp1.pos <- charmatch(grp1, varnames)

      if (is.na(grp1.pos))
         stop(mstyle$stop("Argument 'grp1' must be the name of a variable in the data frame."))

      if (grp1.pos == 0L)
         stop(mstyle$stop("No ambiguous match found for variable name specified via 'grp1' argument."))

   } else {

      grp1.pos <- round(grp1)

      if (grp1.pos < 1 | grp1.pos > nvars)
         stop(mstyle$stop("Specified position of 'grp1' variable does not exist in the data frame."))

   }

   ### get grp1 variable

   grp1 <- data[[grp1.pos]]

   ### make sure there are no missing values in grp1 variable

   if (anyNA(grp1))
      stop(mstyle$stop("Variable specified via 'grp1' argument should not contain missing values."))

   ############################################################################

   ### checks on 'grp2' argument

   if (length(grp2) != 1L)
      stop(mstyle$stop("Argument 'grp2' must of length 1."))

   if (!(is.character(grp2) | is.numeric(grp2)))
      stop(mstyle$stop("Argument 'grp2' must either be a character string or a number."))

   if (is.character(grp2)) {

      grp2.pos <- charmatch(grp2, varnames)

      if (is.na(grp2.pos))
         stop(mstyle$stop("Argument 'grp2' must be the name of a variable in the data frame."))

      if (grp2.pos == 0L)
         stop(mstyle$stop("No ambiguous match found for variable name specified via 'grp2' argument."))

   } else {

      grp2.pos <- round(grp2)

      if (grp2.pos < 1 | grp2.pos > nvars)
         stop(mstyle$stop("Specified position of 'grp2' variable does not exist in the data frame."))

   }

   ### get grp2 variable

   grp2 <- data[[grp2.pos]]

   ### make sure there are no missing values in grp2 variable

   if (anyNA(grp2))
      stop(mstyle$stop("Variable specified via 'grp2' argument should not contain missing values."))

   ############################################################################

   ### get all levels (of grp1 and grp2)

   if (is.factor(grp1) && is.factor(grp2) && identical(levels(grp1), levels(grp2))) {
      lvls <- levels(grp1)
   } else {
      lvls <- sort(unique(c(levels(factor(grp1)), levels(factor(grp2)))))
   }

   ############################################################################

   ### checks on 'last' argument

   ### if last is not specified, place most common grp2 group last

   if (missing(last))
      last <- names(sort(table(grp2), decreasing=TRUE)[1])

   if (length(last) != 1L)
      stop(mstyle$stop("Argument 'last' must be of length one."))

   ### if last is set to NA, leave last unchanged

   if (is.na(last))
      last <- tail(lvls, 1)

   last.pos <- charmatch(last, lvls)

   if (is.na(last.pos))
      stop(mstyle$stop("Could not find specified group in 'grp1' or 'grp2' variables."))

   if (last.pos == 0L)
      stop(mstyle$stop("No ambiguous match found for group specified via 'last' argument."))

   last <- lvls[last.pos]

   ### reorder levels so that the reference level is always last

   lvls <- c(lvls[-last.pos], lvls[last.pos])

   ############################################################################

   ### turn grp1 and grp2 into factors with all levels

   grp1 <- factor(grp1, levels=lvls)
   grp2 <- factor(grp2, levels=lvls)

   ### create contrast matrix

   X <- model.matrix(~ grp1 - 1, contrasts.arg = list(grp1 = "contr.treatment")) - model.matrix(~ grp2 - 1, contrasts.arg = list(grp2 = "contr.treatment"))
   attr(X, "assign") <- NULL
   attr(X, "contrasts") <- NULL

   ### shorten variables names (if shorten=TRUE)

   if (shorten)
      lvls <- .shorten(lvls, minlen=minlen)

   ### add variable names

   colnames(X) <- lvls

   ### append to original data if requested

   if (append)
      X <- cbind(data, X)

   ############################################################################

   return(X)

}
