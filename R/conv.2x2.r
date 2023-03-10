conv.2x2 <- function(ori, ri, x2i, ni, n1i, n2i, correct=TRUE, data, include,
                    var.names=c("ai","bi","ci","di"), append=TRUE, replace="ifna") {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   replace <- match.arg(replace, c("ifna","all"))

   #########################################################################

   if (missing(data))
      data <- NULL

   has.data <- !is.null(data)

   if (is.null(data)) {
      data <- sys.frame(sys.parent())
   } else {
      if (!is.data.frame(data))
         data <- data.frame(data)
   }

   ### checks on var.names argument

   if (length(var.names) != 4L)
      stop(mstyle$stop("Argument 'var.names' must be of length 4."))

   if (any(var.names != make.names(var.names, unique=TRUE))) {
      var.names <- make.names(var.names, unique=TRUE)
      warning(mstyle$warning(paste0("Argument 'var.names' does not contain syntactically valid variable names.\nVariable names adjusted to: var.names = c('", var.names[1], "','", var.names[2], "','", var.names[3], "','", var.names[2], "').")), call.=FALSE)
   }

   #########################################################################

   mf <- match.call()

   ori     <- .getx("ori",     mf=mf, data=data, checknumeric=TRUE)
   ri      <- .getx("ri",      mf=mf, data=data, checknumeric=TRUE)
   x2i     <- .getx("x2i",     mf=mf, data=data, checknumeric=TRUE)
   ni      <- .getx("ni",      mf=mf, data=data, checknumeric=TRUE)
   n1i     <- .getx("n1i",     mf=mf, data=data, checknumeric=TRUE)
   n2i     <- .getx("n2i",     mf=mf, data=data, checknumeric=TRUE)
   correct <- .getx("correct", mf=mf, data=data, default=TRUE)
   include <- .getx("include", mf=mf, data=data)

   if (!.equal.length(ori, ri, x2i, ni, n1i, n2i))
      stop(mstyle$stop("Supplied data vectors are not all of the same length."))

   k <- max(length(ori), length(ri), length(x2i), length(ni), length(n1i), length(n2i))

   if (is.null(ori))
      ori <- rep(NA_real_, k)
   if (is.null(ri))
      ri <- rep(NA_real_, k)
   if (is.null(x2i))
      x2i <- rep(NA_real_, k)
   if (is.null(ni))
      ni <- rep(NA_real_, k)
   if (is.null(n1i))
      n1i <- rep(NA_real_, k)
   if (is.null(n2i))
      n2i <- rep(NA_real_, k)

   ### handle correct argument

   if (length(correct) == 1L)
      correct <- rep(correct, k)

   if (length(correct) != k)
      stop(mstyle$stop(paste0("Length of 'correct' argument (", length(correct), ") does not match length of data (", k, ").")))

   correct[is.na(correct)] <- TRUE

   ### if include is NULL, set to TRUE vector

   if (is.null(include))
      include <- rep(TRUE, k)

   ### turn numeric include vector into logical vector

   include <- .chksubset(include, k, stoponk0=FALSE)

   ### set inputs to NA for rows not to be included

   ori[!include] <- NA_real_
   ri[!include]  <- NA_real_
   x2i[!include] <- NA_real_
   ni[!include]  <- NA_real_
   n1i[!include] <- NA_real_
   n2i[!include] <- NA_real_

   ### round ni, n1i, and n2i

   ni  <- round(ni)
   n1i <- round(n1i)
   n2i <- round(n2i)

   ### checks on values

   if (any(c(ni < 0, n1i < 0, n2i < 0), na.rm=TRUE))
      stop(mstyle$stop("One or more sample sizes or marginal counts are negative."))

   if (any(c(n1i > ni, n2i > ni), na.rm=TRUE))
      stop(mstyle$stop("One or more marginal counts are larger than the sample sizes."))

   if (any(abs(ri) > 1, na.rm=TRUE))
      stop(mstyle$stop("One or more phi coefficients are > 1 or < -1."))

   ### compute marginal proportions for the two variables

   p1i <- n1i / ni
   p2i <- n2i / ni

   #########################################################################

   p11i <- rep(NA_real_, k)

   for (i in seq_len(k)) {

      if (is.na(ni[i]) || is.na(n1i[i]) || is.na(n2i[i]))
         next

      if (!is.na(ori[i])) {

         p1. <- p1i[i]
         p2. <- 1-p1i[i]
         p.1 <- p2i[i]
         p.2 <- 1-p2i[i]

         x <- ori[i] * (p1. + p.1) + p2. - p.1
         y <- sqrt(x^2 - 4 * p1. * p.1 * ori[i] * (ori[i]-1))

         p11i[i] <- (x - y) / (2 * (ori[i] - 1))

      }

      # note: when x2i=0, then sign(0) = 0 and hence ri is automatically 0, which is correct
      # (i.e., we do not want to use the continuity correction in this case)
      if (is.na(ri[i]) && !is.na(x2i[i])) {
         if (correct[i]) {
            ri[i] <- sign(x2i[i]) * (sqrt(abs(x2i[i])/ni[i]) + ni[i] / (2*sqrt(n1i[i]*(ni[i]-n1i[i])*n2i[i]*(ni[i]-n2i[i]))))
         } else {
            ri[i] <- sign(x2i[i]) * sqrt(abs(x2i[i])/ni[i])
         }
      }

      if (is.na(p11i[i]) && !is.na(ri[i]))
         p11i[i] <- p1i[i]*p2i[i] + ri[i] * sqrt(p1i[i]*(1-p1i[i])*p2i[i]*(1-p2i[i]))

   }

   ai <- round(ni * p11i)
   bi <- n1i - ai
   ci <- n2i - ai
   di <- ni - ai - bi - ci

   #print(matrix(c(ai,bi,ci,di), nrow=2, byrow=TRUE))

   ### check for negative cell frequencies

   hasneg <- (ai < 0) | (bi < 0) | (ci < 0) | (di < 0)

   if (any(hasneg, na.rm=TRUE)) {

      warning(mstyle$warning(paste0("There are negative cell frequencies in table", ifelse(sum(hasneg, na.rm=TRUE) > 1, "s ", " "),
                                    paste0(which(hasneg), collapse=","), ".")), call.=FALSE)

      ai[hasneg] <- NA_real_
      bi[hasneg] <- NA_real_
      ci[hasneg] <- NA_real_
      di[hasneg] <- NA_real_

   }

   #########################################################################

   if (has.data && append) {

      if (is.element(var.names[1], names(data))) {
         if (replace=="ifna") {
            data[[var.names[1]]] <- replmiss(data[[var.names[1]]], ai)
         } else {
            data[[var.names[1]]][!is.na(ai)] <- ai[!is.na(ai)]
         }
      } else {
         data <- cbind(data, ai)
         names(data)[length(names(data))] <- var.names[1]
      }

      if (is.element(var.names[2], names(data))) {
         if (replace=="ifna") {
            data[[var.names[2]]] <- replmiss(data[[var.names[2]]], bi)
         } else {
            data[[var.names[2]]][!is.na(bi)] <- bi[!is.na(bi)]
         }
      } else {
         data <- cbind(data, bi)
         names(data)[length(names(data))] <- var.names[2]
      }

      if (is.element(var.names[3], names(data))) {
         if (replace=="ifna") {
            data[[var.names[3]]] <- replmiss(data[[var.names[3]]], ci)
         } else {
            data[[var.names[3]]][!is.na(ci)] <- ai[!is.na(ci)]
         }
      } else {
         data <- cbind(data, ci)
         names(data)[length(names(data))] <- var.names[3]
      }

      if (is.element(var.names[4], names(data))) {
         if (replace=="ifna") {
            data[[var.names[4]]] <- replmiss(data[[var.names[4]]], di)
         } else {
            data[[var.names[4]]][!is.na(di)] <- ai[!is.na(di)]
         }
      } else {
         data <- cbind(data, di)
         names(data)[length(names(data))] <- var.names[4]
      }

   } else {

      data <- data.frame(ai, bi, ci, di)
      names(data) <- var.names

   }

   return(data)

}
