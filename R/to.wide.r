to.wide <- function(data, study, group, ref, grpvars, postfix=c(".1",".2"), addid=TRUE, addcomp=TRUE, adddesign=TRUE, minlen=2) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!(is.data.frame(data)))
      data <- data.frame(data)

   ### get variable names

   varnames <- names(data)

   ### number of variables

   nvars <- length(varnames)

   ############################################################################

   ### checks on 'study' argument

   if (length(study) != 1L)
      stop(mstyle$stop("Argument 'study' must of length 1."))

   if (!(is.character(study) | is.numeric(study)))
      stop(mstyle$stop("Argument 'study' must either be a character string or a scalar."))

   if (is.character(study)) {

      study.pos <- charmatch(study, varnames)

      if (is.na(study.pos))
         stop(mstyle$stop("Argument 'study' must be the name of a variable in the data frame."))

      if (study.pos == 0L)
         stop(mstyle$stop("No ambiguous match found for variable name specified via 'study' argument."))

   } else {

      study.pos <- round(study)

      if (study.pos < 1 | study.pos > nvars)
         stop(mstyle$stop("Specified position of 'study' variable does not exist in the data frame."))

   }

   ### get study variable

   study <- data[[study.pos]]

   ### make sure there are no missing values in study variable

   if (anyNA(study))
      stop(mstyle$stop("Variable specified via 'study' argument should not contain missing values."))

   ############################################################################

   ### checks on 'group' argument

   if (length(group) != 1L)
      stop(mstyle$stop("Argument 'group' must of length 1."))

   if (!(is.character(group) | is.numeric(group)))
      stop(mstyle$stop("Argument 'group' must either be a character string or a scalar."))

   if (is.character(group)) {

      group.pos <- charmatch(group, varnames)

      if (is.na(group.pos))
         stop(mstyle$stop("Argument 'group' must be the name of a variable in the data frame."))

      if (group.pos == 0L)
         stop(mstyle$stop("No ambiguous match found for variable name specified via 'group' argument."))

   } else {

      group.pos <- round(group)

      if (group.pos < 1 | group.pos > nvars)
         stop(mstyle$stop("Specified position of 'group' variable does not exist in the data frame."))

   }

   ### get group variable

   group <- data[[group.pos]]

   ### make sure there are no missing values in group variable

   if (anyNA(group))
      stop(mstyle$stop("Variable specified via 'group' argument should not contain missing values."))

   ### get levels of the group variable

   if (is.factor(group)) {
      lvls <- levels(group)
   } else {
      lvls <- sort(unique(group))
   }

   ############################################################################

   ### checks on 'ref' argument

   ### if ref is not specified, use first group as the reference group

   if (missing(ref))
      ref <- lvls[1]

   if (length(ref) != 1L)
      stop(mstyle$stop("Argument 'ref' must be of length one."))

   ref.pos <- charmatch(ref, lvls)

   if (is.na(ref.pos))
      stop(mstyle$stop("Could not find specified reference level in 'group' variable."))

   if (ref.pos == 0L)
      stop(mstyle$stop("No ambiguous match found for reference level specified via 'ref' argument."))

   ############################################################################

   ### reorder levels and data so that the reference level is always last

   lvls <- c(lvls[-ref.pos], lvls[ref.pos])
   data <- data[order(study, factor(group, levels=lvls)),]

   ### get study and group variables again

   study <- data[[study.pos]]
   group <- data[[group.pos]]

   ############################################################################

   ### checks on 'grpvars' argument

   if (!(is.character(grpvars) | is.numeric(grpvars)))
      stop(mstyle$stop("Argument 'grpvars' must either be a string or numeric vector."))

   if (is.character(grpvars)) {

      grpvars.pos <- unique(charmatch(grpvars, varnames))

      if (any(is.na(grpvars.pos)))
         stop(mstyle$stop("Argument 'grpvars' must be the names of variables in the data frame."))

      if (any(grpvars.pos == 0L))
         stop(mstyle$stop("One or multiple ambiguous matches for variable names specified via 'grpvars' argument."))

   } else {

      grpvars.pos <- unique(round(grpvars))

      if (any(grpvars.pos < 1) | any(grpvars.pos > nvars))
         stop(mstyle$stop("Specified positions of 'grpvars' variables do not exist in the data frame."))

   }

   ### in case the group variable is not specified as part of the group variables, add it

   if (!(group.pos %in% grpvars.pos))
      grpvars.pos <- c(group.pos, grpvars.pos)

   ### and make sure that group.pos is always in the first position of grpvars.pos

   grpvars.pos <- unique(c(group.pos, grpvars.pos))

   ############################################################################

   ### restructure data set into wide format

   restruct <- function(x) {
      if (nrow(x) > 1L) {
         cbind(x[-nrow(x),], x[rep(nrow(x),nrow(x)-1),grpvars.pos])
      } else {
         # to handle one-arm studies
         unname(c(x, rep(NA, length(grpvars.pos))))
      }
   }

   dat <- lapply(split(data, study), restruct)
   dat <- do.call(rbind, dat)

   ### add postfix to outcome variable names

   names(dat)[grpvars.pos] <- paste0(names(dat)[grpvars.pos], postfix[1])
   names(dat)[(nvars+1):ncol(dat)] <- paste0(names(dat)[(nvars+1):ncol(dat)], postfix[2])

   ### fix row names

   rownames(dat) <- seq_len(nrow(dat))

   ############################################################################

   ### generate comp variable

   grps <- .shorten(as.character(data[[group.pos]]), minlen=minlen)

   restruct <- function(x) {
      if (length(x) > 1L) {
         paste0(x[-length(x)], "-", x[length(x)])
      } else {
         NA
      }
   }

   comp <- unlist(sapply(split(grps, study), restruct))

   ### generate design variable

   restruct <- function(x) {
      if (length(x) > 1L) {
         rep(paste0(x, collapse="-"), length(x)-1)
      } else {
         NA
      }
   }

   design <- unlist(sapply(split(grps, study), restruct))

   ############################################################################

   ### add row id to dataset

   if (addid) {
      dat$id <- 1:nrow(dat)
      ### make sure that row id variable is always the first variable in the dataset
      #id.pos <- which(names(dat) == "id")
      #dat <- dat[c(id.pos, seq_along(names(dat))[-id.pos])]
   }

   ### add comp variable to dataset

   if (addcomp)
      dat$comp <- comp

   ### add design variable to dataset

   if (adddesign)
      dat$design <- design

   ############################################################################

   return(dat)

}
