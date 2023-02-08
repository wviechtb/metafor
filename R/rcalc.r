rcalc <- function(x, ni, data, rtoz=FALSE, nfun="min", sparse=FALSE, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!(inherits(x, "formula") || inherits(x, "matrix") || inherits(x, "list")))
      stop(mstyle$stop("Argument 'x' must be either a formula, a matrix, or a list of matrices."))

   if (missing(ni))
      stop(mstyle$stop("Argument 'ni' must be specified."))

   if (is.character(nfun))
      nfun <- match.arg(nfun, c("min", "harmonic", "mean"))

   ### get ... argument and check for extra/superfluous arguments

   ddd <- list(...)

   .chkdots(ddd, c("upper", "simplify", "rowid", "vnames", "noid"))

   if (is.null(ddd$upper)) {
      upper <- TRUE
   } else {
      upper <- ddd$upper
   }

   if (is.null(ddd$simplify)) {
      simplify <- TRUE
   } else {
      simplify <- ddd$simplify
   }

   na.act <- getOption("na.action")
   on.exit(options(na.action=na.act), add=TRUE)

   ############################################################################

   ### in case x is a formula, process it

   if (inherits(x, "formula")) {

      if (missing(data))
         stop(mstyle$stop("Must specify 'data' argument when 'x' is a formula."))

      if (!is.data.frame(data))
         data <- data.frame(data)

      ### extract ni
      mf <- match.call()
      ni <- .getx("ni", mf=mf, data=data, checknumeric=TRUE)

      ### get all variables from data
      options(na.action = "na.pass")
      dat <- get_all_vars(x, data=data)
      options(na.action = na.act)

      ### if no study id has been specified, assume it is a single study
      if (ncol(dat) == 3L) {
         dat[[4]] <- 1
         noid <- TRUE
      } else {
         noid <- FALSE
      }

      vnames <- names(dat)

      ### check that there are really 4 variables
      if (ncol(dat) != 4L)
         stop(mstyle$stop(paste0("Formula should contain 4 variables, but contains ", ncol(dat), " variables.")))

      ### check that there are no missings in the variable identifiers
      if (anyNA(c(dat[[2]],dat[[3]])))
         stop(mstyle$stop("No missing values allowed in variable identifiers."))

      id <- dat[[4]]

      ### check that ni has the same length as there are rows in 'data'
      if (length(ni) != nrow(data))
         stop(mstyle$stop("Argument 'ni' must be of the same length as the data frame specified via 'data'."))

      ### check that there are no missings in the study identifier
      if (anyNA(id))
         stop(mstyle$stop("No missing values allowed in study identifier."))

      ### need these to correctly sort 'dat' and 'V' back into the original order at the end
      ### (and need to order within rows, so that matching works correctly)
      id.var1 <- paste0(id, ".", as.character(dat[[2]]))
      id.var2 <- paste0(id, ".", as.character(dat[[3]]))
      id.var1.id.var2 <- .psort(id.var1, id.var2)
      id.var1 <- id.var1.id.var2[,1]
      id.var2 <- id.var1.id.var2[,2]
      rowid <- paste0(id.var1, ".", id.var2)

      dat <- split(dat, id)
      ni  <- split(ni, id)

      Rlist <- list()
      nmi <- rep(NA, length(ni))

      for (i in seq_along(dat)) {

         if (any(ni[[i]] <= 0, na.rm=TRUE))
            stop(mstyle$stop(paste0("One or more sample sizes are <= 0 in study ", dat[[i]][[4]][[1]], ".")))

         if (is.function(nfun)) {
            nfunnmi <- nfun(ni[[i]])
            if (length(nfunnmi) != 1L)
               stop(mstyle$stop("Function specified via 'nfun' does not return a single value."))
            nmi[i] <- nfunnmi
         } else {
            if (nfun == "min")
               nmi[i] <- min(ni[[i]], na.rm=TRUE)
            if (nfun == "harmonic")
               nmi[i] <- 1 / mean(1/ni[[i]], na.rm=TRUE)
            if (nfun == "mean")
               nmi[i] <- mean(ni[[i]], na.rm=TRUE)
         }

         var1 <- as.character(dat[[i]][[2]])
         var2 <- as.character(dat[[i]][[3]])

         var1.var2 <- paste0(var1, ".", var2)

         var1.var2.eq <- var1 == var2
         if (any(var1.var2.eq))
            stop(mstyle$stop(paste0("Identical var1-var2 pair", ifelse(sum(var1.var2.eq) >= 2L, "s", ""), " (", paste0(var1.var2[var1.var2.eq], collapse=", "), ") in study ", dat[[i]][[4]][[1]], ".")))

         var1.var2.dup <- duplicated(var1.var2)
         if (any(var1.var2.dup))
            stop(mstyle$stop(paste0("Duplicated var1-var2 pair", ifelse(sum(var1.var2.dup) >= 2L, "s", ""), " (", paste0(var1.var2[var1.var2.dup], collapse=", "), ") in study ", dat[[i]][[4]][[1]], ".")))

         ri <- dat[[i]][[1]]

         if (any(abs(ri) > 1, na.rm=TRUE))
            stop(mstyle$stop(paste0("One or more correlations are > 1 or < -1 in study ", dat[[i]][[4]][[1]], ".")))

         vars <- sort(unique(c(var1, var2)))

         Ri <- matrix(NA, nrow=length(vars), ncol=length(vars))
         diag(Ri) <- 1
         rownames(Ri) <- colnames(Ri) <- vars

         for (j in seq_along(var1)) {
            Ri[var1[j],var2[j]] <- Ri[var2[j],var1[j]] <- ri[j]
         }

         Rlist[[i]] <- Ri

      }

      names(Rlist) <- names(dat)

      return(rcalc(Rlist, ni=nmi, simplify=simplify, rtoz=rtoz, sparse=sparse, rowid=rowid, vnames=vnames, noid=noid))

   }

   ############################################################################

   ### in case x is a list, need to loop through elements

   if (is.list(x)) {

      k <- length(x)

      if (length(x) != length(ni))
         stop(mstyle$stop("Argument 'ni' must be of the same length as there are elements in 'x'."))

      res <- list()

      for (i in seq_len(k)) {
         res[[i]] <- rcalc(x[[i]], ni[i], upper=upper, rtoz=rtoz, ...)
      }

      if (is.null(names(x)))
         names(x) <- seq_len(k)

      if (simplify) {

         ki  <- sapply(res, function(x) NROW(x$dat))
         dat <- cbind(id=rep(names(x), times=ki), do.call(rbind, lapply(res, "[[", "dat")))
         if (sparse) {
            V <- bdiag(lapply(res, "[[", "V"))
         } else {
            V <- bldiag(lapply(res, "[[", "V"))
         }
         rownames(V) <- colnames(V) <- unlist(lapply(res, function(x) rownames(x$V)))

         if (!is.null(ddd$rowid)) {
            rowid <- match(ddd$rowid, paste0(dat[[1]], ".", as.character(dat[[2]]), ".", dat[[1]], ".", as.character(dat[[3]])))
            dat <- dat[rowid,]
            V <- V[rowid,rowid]
         }

         if (!is.null(ddd$vnames)) {
            names(dat)[1:3] <- ddd$vnames[c(4,2,3)]
            names(dat)[4]   <- paste0(ddd$vnames[2], ".", ddd$vnames[3])
         }

         if (!is.null(ddd$noid) && ddd$noid) {
            dat[[1]] <- NULL
         }

         rownames(dat) <- seq_len(nrow(dat))
         return(list(dat=dat, V=V))

      } else {

         names(res) <- names(x)
         return(res)

      }

   }

   ############################################################################

   ### check if x is square matrix

   if (!is.matrix(x))
      stop(mstyle$stop("Argument 'x' must be a matrix."))

   if (dim(x)[1] != dim(x)[2])
      stop(mstyle$stop("Argument 'x' must be a square matrix."))

   ### set default dimension names

   dimsx  <- nrow(x)
   dnames <- paste0("x", seq_len(dimsx))

   ### in case x has dimension names, use those

   if (!is.null(rownames(x)))
      dnames <- rownames(x)
   if (!is.null(colnames(x)))
      dnames <- colnames(x)

   ### in case x is a 1x1 (or 0x0) matrix, return nothing

   if (dimsx <= 1L)
      return(list(dat=NULL, V=NULL))

   ### make x symmetric, depending on whether we use upper or lower part

   if (upper) {
      x[lower.tri(x)] <- t(x)[lower.tri(x)]
   } else {
      x[upper.tri(x)] <- t(x)[upper.tri(x)]
   }

   ### check if x is symmetric (can be skipped since x must now be symmetric)

   #if (!isSymmetric(x))
   #   stop(mstyle$stop("Argument 'x' must be a symmetric matrix."))

   ### stack upper/lower triangular part of x into a column vector (this is always done column-wise!)

   if (upper) {
      ri <- cbind(x[upper.tri(x)])
   } else {
      ri <- cbind(x[lower.tri(x)])
   }

   ### check that correlations are in [-1,1]

   if (any(abs(ri) > 1, na.rm=TRUE))
      stop(mstyle$stop("One or more correlations are > 1 or < -1."))

   ### check that sample sizes are > 0

   if (isTRUE(ni <= 0))
      stop(mstyle$stop("One or more sample sizes are <= 0."))

   ### apply r-to-z transformation if requested

   if (rtoz)
      ri <- 1/2 * log((1 + ri)/(1 - ri))

   ### I and J are matrices with 1:dimsx for rows and columns, respectively

   I <- matrix(seq_len(dimsx), nrow=dimsx, ncol=dimsx)
   J <- matrix(seq_len(dimsx), nrow=dimsx, ncol=dimsx, byrow=TRUE)

   ### get upper/lower triangular elements of I and J

   if (upper) {
      I <- I[upper.tri(I)]
      J <- J[upper.tri(J)]
   } else {
      I <- I[lower.tri(I)]
      J <- J[lower.tri(J)]
   }

   ### dimensions in V (must be dimsx*(dimsx-1)/2)

   dimsV <- length(ri)

   ### set up V matrix

   V <- matrix(NA, nrow=dimsV, ncol=dimsV)

   for (ro in seq_len(dimsV)) {
      for (co in seq_len(dimsV)) {

         i <- I[ro]
         j <- J[ro]
         k <- I[co]
         l <- J[co]

         ### Olkin & Finn (1995), equation 5, page 157

         V[ro,co] <- 1/2 * x[i,j]*x[k,l] * (x[i,k]^2 + x[i,l]^2 + x[j,k]^2 + x[j,l]^2) +
                     x[i,k]*x[j,l] + x[i,l]*x[j,k] -
                     (x[i,j]*x[i,k]*x[i,l] + x[j,i]*x[j,k]*x[j,l] + x[k,i]*x[k,j]*x[k,l] + x[l,i]*x[l,j]*x[l,k])

         ### Steiger (1980), equation 2, page 245 (provides the same result)

         #V[ro,co] <- 1/2 * ((x[i,k] - x[i,j]*x[j,k]) * (x[j,l] - x[j,k]*x[k,l]) +
         #                   (x[i,l] - x[i,k]*x[k,l]) * (x[j,k] - x[j,i]*x[i,k]) +
         #                   (x[i,k] - x[i,l]*x[l,k]) * (x[j,l] - x[j,i]*x[i,l]) +
         #                   (x[i,l] - x[i,j]*x[j,l]) * (x[j,k] - x[j,l]*x[l,k]))

         ### Steiger (1980), equation 11, page 247 for r-to-z transformed values

         if (rtoz)
            V[ro,co] <- V[ro,co] / ((1 - x[i,j]^2) * (1 - x[k,l]^2))

      }
   }

   ### divide V by (n-1) for raw correlations and by (n-3) for r-to-z transformed correlations

   if (isTRUE(ni >= 5)) {
      if (rtoz) {
         V <- V/(ni-3)
      } else {
         V <- V/(ni-1)
      }
   } else {
      V <- NA*V
   }

   ### create matrix with var1 and var2 names and sort rowwise

   dmat <- cbind(dnames[I], dnames[J])
   dmat <- t(apply(dmat, 1, sort))

   ### set row/column names for V

   var1.var2 <- paste0(dmat[,1], ".", dmat[,2])
   rownames(V) <- colnames(V) <- var1.var2

   #return(list(dat=data.frame(var1=dmat[,1], var2=dmat[,2], var1.var2=var1.var2, yi=ri, vi=unname(diag(V)), ni=ni, stringsAsFactors=FALSE), V=V))
   return(list(dat=data.frame(var1=dmat[,1], var2=dmat[,2], var1.var2=var1.var2, yi=ri, ni=ni, stringsAsFactors=FALSE), V=V))

}
