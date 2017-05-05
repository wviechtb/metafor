simulate.rma <- function (object, nsim = 1, seed = NULL, ...) {

   if (!inherits(object, "rma"))
      stop("Argument 'object' must be an object of class \"rma\".")

   if (inherits(object, "rma.glmm"))
      stop("Method not yet implemented for objects of class \"rma.glmm\". Sorry!")
   if (inherits(object, "rma.mh"))
      stop("Method not yet implemented for objects of class \"rma.mh\". Sorry!")
   if (inherits(object, "rma.peto"))
      stop("Method not yet implemented for objects of class \"rma.peto\". Sorry!")

   ### as in stats:::simulate.lm
   if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      runif(1)
   if (is.null(seed)) {
      RNGstate <- get(".Random.seed", envir = .GlobalEnv)
   } else {
      R.seed <- get(".Random.seed", envir = .GlobalEnv)
      set.seed(seed)
      RNGstate <- structure(seed, kind = as.list(RNGkind()))
      on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
   }

   ### fitted values
   ftd <- c(object$X %*% object$beta)

   ### for rma.uni objects, just need rnorm() (note: this also covers rma.ls objects)

   if (inherits(object, "rma.uni")) {

      val <- replicate(nsim, rnorm(object$k, mean=ftd, sd=sqrt(object$vi + object$tau2)))

   }

   ### for rma.mv objects, need mvrnorm() from MASS

   if (inherits(object, "rma.mv")) {

      if (!requireNamespace("MASS", quietly=TRUE))
         stop("Please install the 'MASS' package to simulate from this model.")

      val <- replicate(nsim, MASS::mvrnorm(object$k, mu=ftd, Sigma=object$M))

   }

   val <- as.data.frame(val)

   colnames(val) <- paste0("sim_", seq_len(nsim))
   rownames(val) <- object$slab[object$not.na]

   attr(val, "seed") <- RNGstate

   return(val)

}
