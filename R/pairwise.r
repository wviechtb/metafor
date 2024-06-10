pairwise <- function(x, btt, btt2, ...) {

   mstyle <- .get.mstyle()

   if (missing(x)) {
      x <- .getfromenv("pairwise", envir=.metafor)
   }# else {
   #   if (is.atomic(x)) {
   #      btt <- x
   #      x <- .getfromenv("pairwise", envir=.metafor)
   #   }
   #}

   if (is.null(x))
      stop(mstyle$stop("Need to specify 'x' argument."), call.=FALSE)

   .chkclass(class(x), must="rma")

   if (x$int.only)
      stop(mstyle$stop("Cannot construct contrast matrices for intercept-only models."))

   if (missing(btt) || is.null(btt))
      stop(mstyle$stop("Need to specify 'btt' argument."), call.=FALSE)

   ddd <- list(...)

   .chkdots(ddd, c("fixed"))

   fixed <- .chkddd(ddd$fixed, FALSE, .isTRUE(ddd$fixed))

   #########################################################################

   btt <- .set.btt(btt, x$p, x$int.incl, colnames(x$X), fixed=fixed)

   p <- length(btt)

   if (p == 1L)
      stop(mstyle$stop("Need to specify multiple coefficients via argument 'btt' for pairwise comparisons."), call.=FALSE)

   names <- rownames(x$beta)
   connames <- rep("", p*(p-1)/2)

   X <- matrix(0, nrow=p*(p-1)/2, ncol=x$p)
   row <- 0

   for (i in 1:(p-1)) {
      btti <- btt[i]
      for (j in (i+1):p) {
         bttj <- btt[j]
         row <- row + 1
         X[row,btti] <- -1
         X[row,bttj] <- +1
         connames[row] <- paste0(names[btti], "-", names[bttj])
      }
   }

   rownames(X) <- connames

   #########################################################################

   ### in case btt2 is specified, add these coefficients to X

   if (!missing(btt2)) {

      btt <- .set.btt(btt2, x$p, x$int.incl, colnames(x$X), fixed=fixed)

      p <- length(btt)

      Xadd <- matrix(0, nrow=p, ncol=x$p)

      for (i in 1:p) {
         Xadd[i,btt[i]] <- 1
      }

      rownames(Xadd) <- names[btt]

      X <- rbind(Xadd, X)

   }

   #########################################################################

   return(X)

}
