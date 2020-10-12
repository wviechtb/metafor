############################################################################

### function to set default 'btt' value(s) or check specified 'btt' values

.set.btt <- function(btt, p, int.incl, Xnames) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (missing(btt) || is.null(btt)) {

      if (p > 1) {                        ### if the model matrix has more than one column
         if (int.incl) {
            btt <- seq.int(from=2, to=p)     ### and the model has an intercept term, test all coefficients except the intercept
         } else {
            btt <- seq_len(p)                ### and the model does not have an intercept term, test all coefficients
         }
      } else {
         btt <- 1                         ### if the model matrix has a single column, test that single coefficient
      }

   } else {

      if (is.character(btt)) {

         btt <- grep(btt, Xnames)

         if (length(btt) == 0L)
            stop(mstyle$stop("Cannot identify coefficient(s) corresponding to the specified 'btt' string."))

      } else {

         ### round, take unique values, and sort
         btt <- sort(unique(round(btt)))

         ### check for mix of positive and negative values
         if (any(btt < 0) && any(btt > 0))
            stop(mstyle$stop("Cannot mix positive and negative 'btt' values."))

         ### keep/remove from 1:p vector as specified
         btt <- seq_len(p)[btt]

         ### (1:5)[5:6] yields c(5, NA) so remove NAs if this happens
         btt <- btt[!is.na(btt)]

         ### make sure that at least one valid value is left
         if (length(btt) == 0L)
            stop(mstyle$stop("Non-existent coefficients specified via 'btt'."))

      }

   }

   return(btt)

}

### function to format 'btt' values for printing

.format.btt <- function(btt) {

   sav <- c()

   if (length(btt) > 1L) {

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

############################################################################

### function to format p-values
### if showeq=FALSE, c(.001, .00001) becomes c("0.0010", "<.0001")
### if showeq=TRUE,  c(.001, .00001) becomes c("=0.0010", "<.0001")
### if add0=FALSE, "<.0001"; if add0=TRUE, "<0.0001"

.pval <- function(p, digits=4, showeq=FALSE, sep="", add0=FALSE) {

   digits  <- max(digits, 1)
   cutoff  <- paste(c(".", rep(0,digits-1),1), collapse="")
   ncutoff <- as.numeric(cutoff)

   ifelse(is.na(p), paste0(ifelse(showeq, "=", ""), sep, NA),
                    ifelse(p >= ncutoff, paste0(ifelse(showeq, "=", ""), sep, formatC(p, digits=digits, format="f")),
                                         paste0("<", sep, ifelse(add0, "0", ""), cutoff)))

}

### function to format/round values in general

.fcf <- function(x, digits) {

   if (all(is.na(x))) { # since formatC(NA, format="f", digits=2) fails
      x
   } else {
      trimws(formatC(x, format="f", digits=digits))
   }

}

############################################################################

### function to print a named (character) vector right aligned with
### a gap of two spaces between adjacent values and no padding

.print.vector <- function(x) {

   if (is.null(names(x)))
      names(x) <- seq_along(x)

   len.n   <- nchar(names(x))
   len.x   <- nchar(x, keepNA=FALSE)
   len.max <- pmax(len.n, len.x)
   format  <- sapply(len.max, function(x) paste("%", x, "s", sep=""))

   row.n <- paste(sprintf(format, names(x)), collapse="  ")
   row.x <- paste(sprintf(format, x), collapse="  ")

   cat(row.n, "\n", row.x, "\n", sep="")

}

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

   mstyle <- .get.mstyle("crayon" %in% .packages())

   for (i in seq_along(okargs))
      ddd[okargs[i]] <- NULL

   if (length(ddd) > 0L)
      warning(mstyle$warning(paste0("Extra argument", ifelse(length(ddd) > 1L, "s ", " "), "(", paste0("'", names(ddd), "'", collapse=", "), ") disregarded.")), call.=FALSE)

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
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "Log[RR]", "Log Risk Ratio")
         } else {
            lab <- ifelse(short, lab, "Transformed Log Risk Ratio")
            if (atransf.char == "exp" || atransf.char == "transf.exp.int")
               lab <- ifelse(short, "Risk Ratio", "Risk Ratio (log scale)")
            if (transf.char == "exp" || transf.char == "transf.exp.int")
               lab <- ifelse(short, "Risk Ratio", "Risk Ratio")
         }
      }
      if (is.element(measure, c("OR","PETO","D2OR","D2ORN","D2ORL","MPOR","MPORC","MPPETO"))) {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "Log[OR]", "Log Odds Ratio")
         } else {
            lab <- ifelse(short, lab, "Transformed Log Odds Ratio")
            if (atransf.char == "exp" || atransf.char == "transf.exp.int")
               lab <- ifelse(short, "Odds Ratio", "Odds Ratio (log scale)")
            if (transf.char == "exp" || transf.char == "transf.exp.int")
               lab <- ifelse(short, "Odds Ratio", "Odds Ratio")
         }
      }
      if (is.element(measure, c("RD","MPRD"))) {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "Risk Difference", "Risk Difference")
         } else {
            lab <- ifelse(short, lab, "Transformed Risk Difference")
         }
      }
      if (measure == "AS") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "Arcsine RD", "Arcsine Transformed Risk Difference")
         } else {
            lab <- ifelse(short, lab, "Transformed Arcsine Transformed Risk Difference")
         }
      }
      if (measure == "PHI") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "Phi", "Phi Coefficient")
         } else {
            lab <- ifelse(short, lab, "Transformed Phi Coefficient")
         }
      }
      if (measure == "YUQ") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "Yule's Q", "Yule's Q")
         } else {
            lab <- ifelse(short, lab, "Transformed Yule's Q")
         }
      }
      if (measure == "YUY") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "Yule's Y", "Yule's Y")
         } else {
            lab <- ifelse(short, lab, "Transformed Yule's Y")
         }
      }
      ######################################################################
      if (measure == "IRR") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "Log[IRR]", "Log Incidence Rate Ratio")
         } else {
            lab <- ifelse(short, lab, "Transformed Log Incidence Rate Ratio")
            if (atransf.char == "exp" || atransf.char == "transf.exp.int")
               lab <- ifelse(short, "Rate Ratio", "Incidence Rate Ratio (log scale)")
            if (transf.char == "exp" || transf.char == "transf.exp.int")
               lab <- ifelse(short, "Rate Ratio", "Incidence Rate Ratio")
         }
      }
      if (measure == "IRD") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "IRD", "Incidence Rate Difference")
         } else {
            lab <- ifelse(short, lab, "Transformed Incidence Rate Difference")
         }
      }
      if (measure == "IRSD") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "IRSD", "Square Root Transformed Incidence Rate Difference")
         } else {
            lab <- ifelse(short, lab, "Transformed Square Root Transformed Incidence Rate Difference")
         }
      }
      ######################################################################
      if (measure == "MD") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "MD", "Mean Difference")
         } else {
            lab <- ifelse(short, lab, "Transformed Mean Difference")
         }
      }
      if (is.element(measure, c("SMD","SMDH","PBIT","OR2D","OR2DN","OR2DL","SMD1"))) {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "SMD", "Standardized Mean Difference")
         } else {
            lab <- ifelse(short, lab, "Transformed Standardized Mean Difference")
         }
      }
      if (measure == "ROM") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "Log[RoM]", "Log Ratio of Means")
         } else {
            lab <- ifelse(short, lab, "Transformed Log Ratio of Means")
            if (atransf.char == "exp" || atransf.char == "transf.exp.int")
               lab <- ifelse(short, "Ratio of Means", "Ratio of Means (log scale)")
            if (transf.char == "exp" || transf.char == "transf.exp.int")
               lab <- ifelse(short, "Ratio of Means", "Ratio of Means")
         }
      }
      if (measure == "RPB") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "Correlation", "Point-Biserial Correlation")
         } else {
            lab <- ifelse(short, lab, "Transformed Point-Biserial Correlation")
         }
      }
      if (measure == "CVR") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "Log[CVR]", "Log Coefficient of Variation Ratio")
         } else {
            lab <- ifelse(short, lab, "Transformed Log Coefficient of Variation Ratio")
            if (atransf.char == "exp" || atransf.char == "transf.exp.int")
               lab <- ifelse(short, "CVR", "Coefficient of Variation Ratio (log scale)")
            if (transf.char == "exp" || transf.char == "transf.exp.int")
               lab <- ifelse(short, "CVR", "Coefficient of Variation Ratio")
         }
      }
      if (measure == "VR") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "Log[VR]", "Log Variability Ratio")
         } else {
            lab <- ifelse(short, lab, "Transformed Log Variability Ratio")
            if (atransf.char == "exp" || atransf.char == "transf.exp.int")
               lab <- ifelse(short, "VR", "Variability Ratio (log scale)")
            if (transf.char == "exp" || transf.char == "transf.exp.int")
               lab <- ifelse(short, "VR", "Variability Ratio")
         }
      }
      ######################################################################
      if (is.element(measure, c("COR","UCOR","RTET","RBIS"))) {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "Correlation", "Correlation Coefficient")
         } else {
            lab <- ifelse(short, lab, "Transformed Correlation Coefficient")
         }
      }
      if (measure == "ZCOR") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, expression('Fisher\'s ' * z[r]), "Fisher's z Transformed Correlation Coefficient")
         } else {
            lab <- ifelse(short, lab, "Transformed Fisher's z Transformed Correlation Coefficient")
            if (atransf.char == "transf.ztor" || atransf.char == "transf.ztor.int")
               lab <- ifelse(short, "Correlation", "Correlation Coefficient")
            if (transf.char == "transf.ztor" || transf.char == "transf.ztor.int")
               lab <- ifelse(short, "Correlation", "Correlation Coefficient")
         }
      }
      ######################################################################
      if (measure == "PCOR") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "Correlation", "Partial Correlation Coefficient")
         } else {
            lab <- ifelse(short, lab, "Transformed Partial Correlation Coefficient")
         }
      }
      if (measure == "ZPCOR") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, expression('Fisher\'s ' * z[r]), "Fisher's z Transformed Partial Correlation Coefficient")
         } else {
            lab <- ifelse(short, lab, "Transformed Fisher's z Transformed Partial Correlation Coefficient")
            if (atransf.char == "transf.ztor" || atransf.char == "transf.ztor.int")
               lab <- ifelse(short, "Correlation", "Partial Correlation Coefficient")
            if (transf.char == "transf.ztor" || transf.char == "transf.ztor.int")
               lab <- ifelse(short, "Correlation", "Partial Correlation Coefficient")
         }
      }
      if (measure == "SPCOR") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "Correlation", "Semi-Partial Correlation Coefficient")
         } else {
            lab <- ifelse(short, lab, "Transformed Semi-Partial Correlation Coefficient")
         }
      }
      ######################################################################
      if (measure == "PR") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "Proportion", "Proportion")
         } else {
            lab <- ifelse(short, lab, "Transformed Proportion")
         }
      }
      if (measure == "PLN") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "Log[Pr]", "Log Proportion")
         } else {
            lab <- ifelse(short, lab, "Transformed Log Proportion")
            if (atransf.char == "exp" || atransf.char == "transf.exp.int")
               lab <- ifelse(short, "Proportion", "Proportion (log scale)")
            if (transf.char == "exp" || transf.char == "transf.exp.int")
               lab <- ifelse(short, "Proportion", "Proportion")
         }
      }
      if (measure == "PLO") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "Log[Odds]", "Log Odds")
         } else {
            lab <- ifelse(short, lab, "Transformed Log Odds")
            if (atransf.char == "transf.ilogit" || atransf.char == "transf.ilogit.int" || atransf.char == "plogis")
               lab <- ifelse(short, "Proportion", "Proportion (logit scale)")
            if (transf.char == "transf.ilogit" || transf.char == "transf.ilogit.int" || transf.char == "plogis")
               lab <- ifelse(short, "Proportion", "Proportion")
            if (atransf.char == "exp" || atransf.char == "transf.exp.int")
               lab <- ifelse(short, "Odds", "Odds (log scale)")
            if (transf.char == "exp" || transf.char == "transf.exp.int")
               lab <- ifelse(short, "Odds", "Odds")
         }
      }
      if (measure == "PAS") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, expression(arcsin(sqrt(p))), "Arcsine Transformed Proportion")
         } else {
            lab <- ifelse(short, lab, "Transformed Arcsine Transformed Proportion")
            if (atransf.char == "transf.iarcsin" || atransf.char == "transf.iarcsin.int")
               lab <- ifelse(short, "Proportion", "Proportion (arcsine scale)")
            if (transf.char == "transf.iarcsin" || transf.char == "transf.iarcsin.int")
               lab <- ifelse(short, "Proportion", "Proportion")
         }
      }
      if (measure == "PFT") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "PFT", "Double Arcsine Transformed Proportion")
         } else {
            lab <- ifelse(short, lab, "Transformed Double Arcsine Transformed Proportion")
            if (atransf.char == "transf.ipft.hm")
               lab <- ifelse(short, "Proportion", "Proportion")
            if (transf.char == "transf.ipft.hm")
               lab <- ifelse(short, "Proportion", "Proportion")
         }
      }
      ######################################################################
      if (measure == "IR") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "Rate", "Incidence Rate")
         } else {
            lab <- ifelse(short, lab, "Transformed Incidence Rate")
         }
      }
      if (measure == "IRLN") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "Log[IR]", "Log Incidence Rate")
         } else {
            lab <- ifelse(short, lab, "Transformed Log Incidence Rate")
            if (atransf.char == "exp" || atransf.char == "transf.exp.int")
               lab <- ifelse(short, "Rate", "Incidence Rate (log scale)")
            if (transf.char == "exp" || transf.char == "transf.exp.int")
               lab <- ifelse(short, "Rate", "Incidence Rate")
         }
      }
      if (measure == "IRS") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "Sqrt[IR]", "Square Root Transformed Incidence Rate")
         } else {
            lab <- ifelse(short, lab, "Transformed Square Root Transformed Incidence Rate")
            if (atransf.char == "transf.isqrt" || atransf.char == "transf.isqrt.int")
               lab <- ifelse(short, "Rate", "Incidence Rate (square root scale)")
            if (transf.char == "transf.isqrt" || transf.char == "transf.isqrt.int")
               lab <- ifelse(short, "Rate", "Incidence Rate")
         }
      }
      if (measure == "IRFT") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "IRFT", "Freeman-Tukey Transformed Incidence Rate")
         } else {
            lab <- ifelse(short, lab, "Transformed Freeman-Tukey Transformed Incidence Rate")
         }
      }
      ######################################################################
      if (measure == "MN") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "Mean", "Mean")
         } else {
            lab <- ifelse(short, lab, "Transformed Mean")
         }
      }
      if (measure == "MNLN") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "Log[Mean]", "Log Mean")
         } else {
            lab <- ifelse(short, lab, "Transformed Log Mean")
            if (atransf.char == "exp" || atransf.char == "transf.exp.int")
               lab <- ifelse(short, "Mean", "Mean (log scale)")
            if (transf.char == "exp" || transf.char == "transf.exp.int")
               lab <- ifelse(short, "Mean", "Mean")
         }
      }
      if (measure == "CVLN") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "Log[CV]", "Log Coefficient of Variation")
         } else {
            lab <- ifelse(short, lab, "Transformed Log Coefficient of Variation")
            if (atransf.char == "exp" || atransf.char == "transf.exp.int")
               lab <- ifelse(short, "CV", "Coefficient of Variation (log scale)")
            if (transf.char == "exp" || transf.char == "transf.exp.int")
               lab <- ifelse(short, "CV", "Coefficient of Variation")
         }
      }
      if (measure == "SDLN") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "Log[SD]", "Log Standard Deviation")
         } else {
            lab <- ifelse(short, lab, "Transformed Log Standard Deviation")
            if (atransf.char == "exp" || atransf.char == "transf.exp.int")
               lab <- ifelse(short, "SD", "Standard Deviation (log scale)")
            if (transf.char == "exp" || transf.char == "transf.exp.int")
               lab <- ifelse(short, "SD", "Standard Deviation")
         }
      }
      ######################################################################
      if (measure == "MC") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "Mean Change", "Mean Change")
         } else {
            lab <- ifelse(short, lab, "Transformed Mean Change")
         }
      }
      if (is.element(measure, c("SMCC","SMCR","SMCRH"))) {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "SMC", "Standardized Mean Change")
         } else {
            lab <- ifelse(short, lab, "Transformed Standardized Mean Change")
         }
      }
      if (measure == "ROMC") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "Log[RoM]", "Log Ratio of Means")
         } else {
            lab <- ifelse(short, lab, "Transformed Log Ratio of Means")
            if (atransf.char == "exp" || atransf.char == "transf.exp.int")
               lab <- ifelse(short, "Ratio of Means", "Ratio of Means (log scale)")
            if (transf.char == "exp" || transf.char == "transf.exp.int")
               lab <- ifelse(short, "Ratio of Means", "Ratio of Means")
         }
      }
      if (measure == "CVRC") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "Log[CVR]", "Log Coefficient of Variation Ratio")
         } else {
            lab <- ifelse(short, lab, "Transformed Log Coefficient of Variation Ratio")
            if (atransf.char == "exp" || atransf.char == "transf.exp.int")
               lab <- ifelse(short, "CVR", "Coefficient of Variation Ratio (log scale)")
            if (transf.char == "exp" || transf.char == "transf.exp.int")
               lab <- ifelse(short, "CVR", "Coefficient of Variation Ratio")
         }
      }
      if (measure == "VRC") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "Log[VR]", "Log Variability Ratio")
         } else {
            lab <- ifelse(short, lab, "Transformed Log Variability Ratio")
            if (atransf.char == "exp" || atransf.char == "transf.exp.int")
               lab <- ifelse(short, "VR", "Variability Ratio (log scale)")
            if (transf.char == "exp" || transf.char == "transf.exp.int")
               lab <- ifelse(short, "VR", "Variability Ratio")
         }
      }
      ######################################################################
      if (measure == "ARAW") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, "Alpha", "Cronbach's alpha")
         } else {
            lab <- ifelse(short, lab, "Transformed Cronbach's alpha")
         }
      }
      if (measure == "AHW") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, expression('Alpha'[HW]), "Transformed Cronbach's alpha")
         } else {
            lab <- ifelse(short, lab, "Transformed Cronbach's alpha")
            if (atransf.char == "transf.iahw")
               lab <- ifelse(short, "Alpha", "Cronbach's alpha")
            if (transf.char == "transf.iahw")
               lab <- ifelse(short, "Alpha", "Cronbach's alpha")
         }
      }
      if (measure == "ABT") {
         if (transf.char == "FALSE" && atransf.char == "FALSE") {
            lab <- ifelse(short, expression('Alpha'[B]), "Transformed Cronbach's alpha")
         } else {
            lab <- ifelse(short, lab, "Transformed Cronbach's alpha")
            if (atransf.char == "transf.iabt")
               lab <- ifelse(short, "Alpha", "Cronbach's alpha")
            if (transf.char == "transf.iabt")
               lab <- ifelse(short, "Alpha", "Cronbach's alpha")
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
      if (is.null(.mstyle$body)) {
         body <- crayon::reset
      } else {
         body <- .mstyle$body
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
      } else {
         legend <- .mstyle$legend
      }

   } else {

      tmp <- function(...) paste0(...)
      section <- tmp
      header  <- tmp
      body    <- tmp
      text    <- tmp
      result  <- tmp
      stop    <- tmp
      warning <- tmp
      message <- tmp
      verbose <- tmp
      legend  <- tmp

   }

   return(list(section=section, header=header, body=body, text=text, result=result, stop=stop, warning=warning, message=message, verbose=verbose, legend=legend))

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

.print.table <- function(x, mstyle) {

   is.header <- !grepl(" [-0-9]", x)

   for (i in seq_along(x)) {
      if (is.header[i]) {
         x[i] <- trimws(x[i], which="right")
         x[i] <- mstyle$header(x[i])
      } else {
         x[i] <- mstyle$body(x[i])
      }
      cat(x[i], "\n")
   }

}

#.set.mstyle.1 <- parse(text=".mstyle <- list(section=make_style(\"gray90\")$bold, header=make_style(\"skyblue1\")$bold$underline, body=make_style(\"skyblue2\"), text=make_style(\"slateblue3\"), result=make_style(\"slateblue1\"))")
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

   ### so we can still print objects created with older metafor versions (where xdigit will be just an unnamed scalar)
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

############################################################################

### to register getfit method for 'rma.uni' and 'rma.mv' objects: eval(metafor:::.glmulti)

.glmulti <- parse(text="

if (!(\"glmulti\" %in% .packages()))
   stop(\"Need to load the 'glmulti' package first to use this code.\")

setOldClass(\"rma.uni\")

setMethod(\"getfit\", \"rma.uni\", function(object, ...) {
   if (object$test==\"z\") {
      cbind(estimate=coef(object), se=sqrt(diag(vcov(object))), df=Inf)
   } else {
      cbind(estimate=coef(object), se=sqrt(diag(vcov(object))), df=object$k-object$p)
   }
})

setOldClass(\"rma.mv\")

setMethod(\"getfit\", \"rma.mv\", function(object, ...) {
   if (object$test==\"z\") {
      cbind(estimate=coef(object), se=sqrt(diag(vcov(object))), df=Inf)
   } else {
      cbind(estimate=coef(object), se=sqrt(diag(vcov(object))), df=object$k-object$p)
   }
})

setOldClass(\"rma.glmm\")

setMethod(\"getfit\", \"rma.glmm\", function(object, ...) {
   if (object$test==\"z\") {
      cbind(estimate=coef(object), se=sqrt(diag(vcov(object))), df=Inf)
   } else {
      cbind(estimate=coef(object), se=sqrt(diag(vcov(object))), df=object$k-object$p)
   }
})

")

### helper functions to make MuMIn work together with metafor

.MuMIn <- parse(text="

makeArgs.rma <- function (obj, termNames, comb, opt, ...) {
   ret <- MuMIn:::makeArgs.default(obj, termNames, comb, opt)
   names(ret)[1L] <- \"mods\"
   ret
}

coefTable.rma <- function (model, ...) {
  MuMIn:::.makeCoefTable(model$b, model$se, coefNames = rownames(model$b))
}

")

### helper functions to make mice work together with metafor (note: no longer
### needed, as there are glance and tidy methods for rma objects in broom now)

.mice <- parse(text="

glance.rma <- function (x, ...)
   data.frame(df.residual=df.residual(x))

tidy.rma <- function (x, ...) {
   ret <- coef(summary(x))
   colnames(ret)[2] <- \"std.error\"
   ret$term <- rownames(ret)
   return(ret)
}

")

############################################################################

### shorten a string vector so that elements remain distinguishable

.shorten <- function(x, minlen) {

   y <- x

   x <- c(na.omit(x))

   n <- length(unique(x))

   maxlen <- max(nchar(unique(x)))

   for (l in 1:maxlen) {
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

.tes.intfun <- function(x, theta, tau, sei, H0, alternative, crit) {
   if (alternative == "two.sided")
      pow <- (pnorm(crit, mean=(x-H0)/sei, 1, lower.tail=FALSE) + pnorm(-crit, mean=(x-H0)/sei, 1, lower.tail=TRUE))
   if (alternative == "greater")
      pow <- pnorm(crit, mean=(x-H0)/sei, 1, lower.tail=FALSE)
   if (alternative == "less")
      pow <- pnorm(crit, mean=(x-H0)/sei, 1, lower.tail=TRUE)
   res <- pow * dnorm(x, theta, tau)
   return(res)
}

.tes.lim <- function(theta, yi, vi, H0, alternative, alpha, tau2, test, tes.alternative, progbar, tes.alpha, correct, rel.tol, subdivisions, tau2.lb) {
   pval <- tes.default(x=yi, vi=vi, H0=H0, alternative=alternative, alpha=alpha, theta=theta, tau2=tau2, test=test, tes.alternative=tes.alternative, progbar=progbar,
                       tes.alpha=tes.alpha, correct=correct, rel.tol=rel.tol, subdivisions=subdivisions, tau2.lb=tau2.lb, find.lim=FALSE)$pval
   #cat("theta = ", theta, " pval = ", pval, "\n")
   return(pval - tes.alpha)
}

############################################################################
