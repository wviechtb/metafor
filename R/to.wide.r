to.wide <- function(data, study, grp, ref, grpvars, postfix=c(".1",".2"),
addid=TRUE, addcomp=TRUE, adddesign=TRUE, minlen=2, var.names=c("id","comp","design")) {

   mstyle <- .get.mstyle()

   if (missing(data))
      stop(mstyle$stop("Argument 'data' must be specified."))

   if (!is.data.frame(data))
      data <- data.frame(data)

   # get variable names

   varnames <- names(data)

   # number of variables

   nvars <- length(varnames)

   # checks on the 'var.names' argument

   if (length(var.names) != 3L)
      stop(mstyle$stop("Argument 'var.names' must of length 3."))

   if (!inherits(var.names, "character"))
      stop(mstyle$stop("Argument 'var.names' must of vector with character strings."))

   if (any(var.names != make.names(var.names, unique=TRUE))) {
      var.names <- make.names(var.names, unique=TRUE)
      warning(mstyle$warning(paste0("Argument 'var.names' does not contain syntactically valid variable names.\nVariable names adjusted to: var.names = c('", var.names[1], "','", var.names[2], "','", var.names[3], "').")), call.=FALSE)
   }

   ############################################################################

   # checks on the 'study' argument

   if (missing(study))
      stop(mstyle$stop("Argument 'study' must be specified."))

   if (length(study) != 1L)
      stop(mstyle$stop("Argument 'study' must of length 1."))

   if (!(is.character(study) | is.numeric(study)))
      stop(mstyle$stop("Argument 'study' must either be a character string or a scalar."))

   if (is.character(study)) {

      study.pos <- charmatch(study, varnames)

      if (is.na(study.pos) || study.pos == 0L)
         stop(mstyle$stop("Could not find or uniquely identify variable specified via the 'study' argument."))

   } else {

      study.pos <- round(study)

      if (study.pos < 1 | study.pos > nvars)
         stop(mstyle$stop("Specified position of 'study' variable does not exist in the data frame."))

   }

   # get study variable

   study <- data[[study.pos]]

   # make sure there are no missing values in study variable

   if (anyNA(study))
      stop(mstyle$stop("Variable specified via 'study' argument should not contain missing values."))

   ############################################################################

   # checks on the 'grp' argument

   if (missing(grp))
      stop(mstyle$stop("Argument 'grp' must be specified."))

   if (length(grp) != 1L)
      stop(mstyle$stop("Argument 'grp' must of length 1."))

   if (!(is.character(grp) || is.numeric(grp)))
      stop(mstyle$stop("Argument 'grp' must either be a character string or a scalar."))

   if (is.character(grp)) {

      grp.pos <- charmatch(grp, varnames)

      if (is.na(grp.pos) || grp.pos == 0L)
         stop(mstyle$stop("Could not find or uniquely identify variable specified via the 'grp' argument."))

   } else {

      grp.pos <- round(grp)

      if (grp.pos < 1 | grp.pos > nvars)
         stop(mstyle$stop("Specified position of 'grp' variable does not exist in the data frame."))

   }

   # get grp variable

   grp <- data[[grp.pos]]

   # make sure there are no missing values in group variable

   if (anyNA(grp))
      stop(mstyle$stop("Variable specified via 'grp' argument should not contain missing values."))

   # get levels of the group variable

   if (is.factor(grp)) {
      lvls <- levels(grp)
   } else {
      lvls <- sort(unique(grp))
   }

   ############################################################################

   # checks on the 'ref' argument

   # if ref is not specified, use the most common group as the reference group

   if (missing(ref))
      ref <- names(sort(table(grp), decreasing=TRUE)[1])

   if (length(ref) != 1L)
      stop(mstyle$stop("Argument 'ref' must be of length one."))

   ref.pos <- charmatch(ref, lvls)

   if (is.na(ref.pos) || ref.pos == 0L)
      stop(mstyle$stop("Could not find or uniquely identify reference group specified via the 'ref' argument."))

   ############################################################################

   # reorder levels and data so that the reference level is always last

   lvls <- c(lvls[-ref.pos], lvls[ref.pos])
   data <- data[order(study, factor(grp, levels=lvls)),]

   # get study and group variables again

   study <- data[[study.pos]]
   grp   <- data[[grp.pos]]

   ############################################################################

   # checks on the 'grpvars' argument

   if (!(is.character(grpvars) || is.numeric(grpvars)))
      stop(mstyle$stop("Argument 'grpvars' must either be a string or numeric vector."))

   if (is.character(grpvars)) {

      grpvars.pos <- unique(charmatch(grpvars, varnames))

      if (anyNA(grpvars.pos) || any(grpvars.pos == 0L))
         stop(mstyle$stop("Could not find or uniquely identify variable(s) specified via the 'grpvars' argument."))

   } else {

      grpvars.pos <- unique(round(grpvars))

      if (any(grpvars.pos < 1) | any(grpvars.pos > nvars))
         stop(mstyle$stop("Specified positions of 'grpvars' variables do not exist in the data frame."))

   }

   # in case the group variable is not specified as part of the group variables, add it

   if (!(grp.pos %in% grpvars.pos))
      grpvars.pos <- c(grp.pos, grpvars.pos)

   # and make sure that grp.pos is always in the first position of grpvars.pos

   grpvars.pos <- union(grp.pos, grpvars.pos)

   ############################################################################

   # restructure data set into wide format

   restruct <- function(x) {
      if (nrow(x) > 1L) {
         cbind(x[-nrow(x),], x[rep(nrow(x),nrow(x)-1L),grpvars.pos])
      } else {
         # to handle one-arm studies
         unname(c(x, rep(NA, length(grpvars.pos))))
      }
   }

   dat <- lapply(split(data, study), restruct)
   dat <- do.call(rbind, dat)

   # add postfix to outcome variable names

   names(dat)[grpvars.pos] <- paste0(names(dat)[grpvars.pos], postfix[1])
   names(dat)[(nvars+1):ncol(dat)] <- paste0(names(dat)[(nvars+1):ncol(dat)], postfix[2])

   # fix row names

   rownames(dat) <- seq_len(nrow(dat))

   ############################################################################

   # generate comp variable

   grps <- .shorten(as.character(data[[grp.pos]]), minlen=minlen)

   restruct <- function(x) {
      if (length(x) > 1L) {
         paste0(x[-length(x)], "-", x[length(x)])
      } else {
         NA
      }
   }

   comp <- unlist(sapply(split(grps, study), restruct))

   # generate design variable

   restruct <- function(x) {
      if (length(x) > 1L) {
         rep(paste0(x, collapse="-"), length(x)-1L)
      } else {
         NA
      }
   }

   design <- unlist(sapply(split(grps, study), restruct))

   ############################################################################

   # add row id to dataset

   if (addid) {

      dat[[var.names[1]]] <- seq_len(nrow(dat))

      # make sure that row id variable is always the first variable in the dataset
      #id.pos <- which(names(dat) == "id")
      #dat <- dat[c(id.pos, seq_along(names(dat))[-id.pos])]

   }

   # add comp variable to dataset

   if (addcomp)
      dat[[var.names[2]]] <- comp

   # add design variable to dataset

   if (adddesign)
      dat[[var.names[3]]] <- design

   ############################################################################

   return(dat)

}
