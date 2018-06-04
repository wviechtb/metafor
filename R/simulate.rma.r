simulate.rma <- function (object, nsim = 1, seed = NULL, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(object, "rma"))
      stop(mstyle$stop("Argument 'object' must be an object of class \"rma\"."))

   if (inherits(object, "rma.glmm"))
      stop(mstyle$stop("Method not available for objects of class \"rma.glmm\"."))

   if (inherits(object, "rma.mh"))
      stop(mstyle$stop("Method not available for objects of class \"rma.mh\"."))

   if (inherits(object, "rma.peto"))
      stop(mstyle$stop("Method not available for objects of class \"rma.peto\"."))

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

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
         stop(mstyle$stop("Please install the 'MASS' package to simulate from this model."))

      val <- replicate(nsim, MASS::mvrnorm(1, mu=ftd, Sigma=object$M))

   }

   res <- matrix(NA_real_, nrow=object$k.f, ncol=nsim)
   res[object$not.na,] <- val
   res <- as.data.frame(res)

   rownames(res) <- object$slab
   colnames(res) <- paste0("sim_", seq_len(nsim))

   if (na.act == "na.omit")
      res <- res[object$not.na,,drop=FALSE]

   attr(res, "seed") <- RNGstate

   return(res)

}
