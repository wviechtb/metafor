trimfill.rma.uni <- function(x, side, estimator="L0", maxiter=100, verbose=FALSE, ilim, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="rma.uni", notav=c("robust.rma", "rma.ls", "rma.nn", "rma.uni.selmodel"))

   if (!x$int.only)
      stop(mstyle$stop("Trim-and-fill method only applicable for models without moderators."))

   if (missing(side))
      side <- NULL

   estimator <- match.arg(estimator, c("L0", "R0", "Q0"))

   if (x$k == 1)
      stop(mstyle$stop("Stopped because k = 1."))

   #########################################################################

   yi <- x$yi
   vi <- x$vi
   wi <- x$weights
   ni <- x$ni

   ### determine side (if none is specified)

   if (is.null(side)) {
      args <- list(yi=yi, vi=vi, weights=wi, mods=sqrt(vi), method=x$method, weighted=x$weighted, outlist="beta=beta", ...)
      res <- suppressWarnings(.do.call(rma.uni, args))
      ### TODO: add check in case there are problems with fitting the model
      if (res$beta[2] < 0) {
         side <- "right"
      } else {
         side <- "left"
      }
   } else {
      side <- match.arg(side, c("left", "right"))
   }

   ### flip data if examining right side

   if (side == "right")
      yi <- -1*yi

   ### sort data by increasing yi

   ix <- sort(yi, index.return=TRUE)$ix
   yi <- yi[ix]
   vi <- vi[ix]
   wi <- wi[ix]
   ni <- ni[ix]

   #########################################################################

   k <- length(yi)

   k0.sav <- -1
   k0     <-  0 ### estimated number of missing studies
   iter   <-  0 ### iteration counter

   if (verbose)
      cat("\n")

   while (abs(k0 - k0.sav) > 0) {

      k0.sav <- k0 ### save current value of k0

      iter <- iter + 1

      if (iter > maxiter)
         stop(mstyle$stop("Trim and fill algorithm did not converge."))

      ### truncated data

      yi.t <- yi[seq_len(k-k0)]
      vi.t <- vi[seq_len(k-k0)]
      wi.t <- wi[seq_len(k-k0)]

      args <- list(yi=yi.t, vi=vi.t, weights=wi.t, method=x$method, weighted=x$weighted, outlist="beta=beta", ...)
      res <- suppressWarnings(.do.call(rma.uni, args))

      ### intercept estimate based on truncated data

      beta <- c(res$beta)

      yi.c     <- yi - beta                            ### centered values
      yi.c.r   <- rank(abs(yi.c), ties.method="first") ### ranked absolute centered values
      yi.c.r.s <- sign(yi.c) * yi.c.r                  ### signed ranked centered values

      ### estimate the number of missing studies with the R0 estimator

      if (estimator == "R0") {
         k0 <- (k - max(-1*yi.c.r.s[yi.c.r.s < 0])) - 1
         se.k0 <- sqrt(2*max(0,k0) + 2)
      }

      ### estimate the number of missing studies with the L0 estimator

      if (estimator == "L0") {
         Sr <- sum(yi.c.r.s[yi.c.r.s > 0])
         k0 <- (4*Sr - k*(k+1)) / (2*k - 1)
         varSr <- 1/24 * (k*(k+1)*(2*k+1) + 10*k0^3 + 27*k0^2 + 17*k0 - 18*k*k0^2 - 18*k*k0 + 6*k^2*k0)
         se.k0 <- 4*sqrt(varSr) / (2*k - 1)
      }

      ### estimate the number of missing studies with the Q0 estimator

      if (estimator == "Q0") {
         Sr <- sum(yi.c.r.s[yi.c.r.s > 0])
         k0 <- k - 1/2 - sqrt(2*k^2 - 4*Sr + 1/4)
         varSr <- 1/24 * (k*(k+1)*(2*k+1) + 10*k0^3 + 27*k0^2 + 17*k0 - 18*k*k0^2 - 18*k*k0 + 6*k^2*k0)
         se.k0 <- 2*sqrt(varSr) / sqrt((k-1/2)^2 - k0*(2*k - k0 -1))
      }

      ### round k0 and make sure that k0 is non-negative

      k0    <- max(0, round(k0))
      se.k0 <- max(0, se.k0)

      if (verbose)
         cat(mstyle$verbose(paste0("Iteration: ", formatC(iter, format="f", digits=0, width=nchar(maxiter), flag="-"),
                                   "  missing = ", formatC(k0, format="f", digits=0, width=nchar(k), flag="-"),
                                   "  beta = ",    formatC(ifelse(side == "right", -1*beta, beta), format="f", digits=x$digits[["est"]]), "\n")))

   }

   #########################################################################

   ### if estimated number of missing studies is > 0

   if (k0 > 0) {

      ### flip data back if side is right

      if (side == "right") {
         yi.c <- -1 * (yi.c - beta)
      } else {
         yi.c <- yi.c - beta
      }

      ### create filled-in data set

      yi.fill <- c(x$yi.f, -1*yi.c[(k-k0+1):k])

      ### apply limits if specified

      if (!missing(ilim)) {
         ilim <- sort(ilim)
         if (length(ilim) != 2L)
            stop(mstyle$stop("Argument 'ilim' must be of length 2."))
         yi.fill[yi.fill < ilim[1]] <- ilim[1]
         yi.fill[yi.fill > ilim[2]] <- ilim[2]
      }

      vi.fill <- c(x$vi.f, vi[(k-k0+1):k])
      wi.fill <- c(x$weights.f, wi[(k-k0+1):k])
      ni.fill <- c(x$ni.f, ni[(k-k0+1):k])

      ### add measure attribute to the yi.fill vector

      attr(yi.fill, "measure") <- x$measure

      ### fit model with imputed data

      args <- list(yi=yi.fill, vi=vi.fill, weights=wi.fill, ni=ni.fill, method=x$method, weighted=x$weighted, digits=x$digits, ...)
      res <- suppressWarnings(.do.call(rma.uni, args))

      ### fill, ids, and slab are of length 'k.f + k0' (i.e., subsetted but with NAs)

      res$fill <- c(rep(FALSE,x$k.f), rep(TRUE,k0))
      res$ids  <- c(x$ids, (max(x$ids)+1):(max(x$ids)+k0))

      if (x$slab.null) {
         res$slab <- c(paste("Study", x$ids), paste("Filled", seq_len(k0)))
      } else {
         res$slab <- c(x$slab, paste("Filled", seq_len(k0)))
      }
      res$slab.null <- FALSE

   } else {

      ### in case 0 studies are imputed

      res      <- x
      res$fill <- rep(FALSE,k)

   }

   res$k0     <- k0
   res$se.k0  <- se.k0
   res$side   <- side
   res$k0.est <- estimator
   res$k.all  <- x$k.all + k0

   if (estimator == "R0") {
      m <- -1:(k0-1)
      res$p.k0 <- 1 - sum(choose(0+m+1, m+1) * 0.5^(0+m+2))
   } else {
      res$p.k0 <- NA
   }

   class(res) <- c("rma.uni.trimfill", class(res))
   return(res)

}
