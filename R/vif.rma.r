vif.rma <- function(x, btt, att, table=FALSE, reestimate=FALSE, sim=FALSE, progbar=TRUE,
                    seed=NULL, parallel="no", ncpus=1, cl, digits, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="rma")

   # allow vif() for 'rma.glmm', 'robust.rma', and 'rma.uni.selmodel' objects based on the same principle (but not sim/reestimate)

   if (missing(digits)) {
      digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
   }

   ### determine for which types of coefficients (G)VIFs will be computed

   vif.loc <- !x$int.only

   if (inherits(x, "rma.ls") && !x$Z.int.only) {
      vif.scale <- TRUE
   } else {
      vif.scale <- FALSE
   }

   if (!vif.loc && !vif.scale)
      stop(mstyle$stop("VIFs not applicable to intercept-only models."))

   if (!is.null(seed))
      set.seed(seed)

   ddd <- list(...)

   .chkdots(ddd, c("fixed", "intercept", "time", "LB", "joinb", "joina"))

   if (is.null(ddd$fixed)) {
      fixed <- FALSE
   } else {
      fixed <- .isTRUE(ddd$fixed)
   }

   if (is.null(ddd$intercept)) {
      intercept <- FALSE
   } else {
      intercept <- .isTRUE(ddd$intercept)
   }

   if (is.null(ddd$joinb)) {
      joinb <- NULL
   } else {
      joinb <- ddd$joinb
   }

   if (is.null(ddd$joina)) {
      joina <- NULL
   } else {
      joina <- ddd$joina
   }

   if (.isTRUE(ddd$time))
      time.start <- proc.time()

   ### process 'sim' argument (if TRUE, set sim to 1000, otherwise use given value)

   if (is.logical(sim)) {
      sim <- ifelse(isTRUE(sim), 1000, 0)
   } else {
      sim <- round(sim)
      if (sim <= 1)
         stop(mstyle$stop("Argument 'sim' must be >= 2."))
   }

   ### do not allow sim and reestimate for 'rma.glmm', 'robust.rma', and 'rma.uni.selmodel' objects

   if (sim >= 2 && inherits(x, "rma.glmm"))
      stop(mstyle$stop("Cannot use 'sim' with models of class 'rma.glmm'."))

   if (sim >= 2 && inherits(x, "robust.rma"))
      stop(mstyle$stop("Cannot use 'sim' with models of class 'robust.rma'."))

   if (sim >= 2 && inherits(x, "rma.uni.selmodel"))
      stop(mstyle$stop("Cannot use 'sim' with models of class 'rma.uni.selmodel'."))

   if (reestimate && inherits(x, "rma.glmm"))
      stop(mstyle$stop("Cannot use 'restimate=TRUE' with models of class 'rma.glmm'."))

   if (reestimate && inherits(x, "robust.rma"))
      stop(mstyle$stop("Cannot use 'restimate=TRUE' with models of class 'robust.rma'."))

   if (reestimate && inherits(x, "rma.uni.selmodel"))
      stop(mstyle$stop("Cannot use 'restimate=TRUE' with models of class 'rma.uni.selmodel'."))

   ### check if btt/att have been specified

   bttmiss <- missing(btt) || is.null(btt)
   attmiss <- missing(att) || is.null(att)

   if (!attmiss && !inherits(x, "rma.ls"))
      stop(mstyle$stop("Argument 'att' only relevant for location-scale models."))

   ### handle parallel (and related) arguments

   parallel <- match.arg(parallel, c("no", "snow", "multicore"))

   if (parallel == "no" && ncpus > 1)
      parallel <- "snow"

   if (missing(cl))
      cl <- NULL

   if (!is.null(cl) && inherits(cl, "SOCKcluster")) {
      parallel <- "snow"
      ncpus <- length(cl)
   }

   if (parallel == "snow" && ncpus < 2)
      parallel <- "no"

   if (sim <= 1)
      parallel <- "no"

   if (parallel == "snow" || parallel == "multicore") {

      if (!requireNamespace("parallel", quietly=TRUE))
         stop(mstyle$stop("Please install the 'parallel' package for parallel processing."))

      ncpus <- as.integer(ncpus)

      if (ncpus < 1L)
         stop(mstyle$stop("Argument 'ncpus' must be >= 1."))

   }

   if (parallel == "snow") {

      if (is.null(cl)) {
         cl <- parallel::makePSOCKcluster(ncpus)
         on.exit(parallel::stopCluster(cl), add=TRUE)
      }

   }

   if (!progbar) {
      pbo <- pbapply::pboptions(type="none")
      on.exit(pbapply::pboptions(pbo), add=TRUE)
   }


   #########################################################################

   if (vif.loc) {

      ### process/set btt argument

      if (bttmiss) {
         if (x$intercept && !intercept) {
            btt <- as.list(2:x$p)
         } else {
            btt <- as.list(seq_len(x$p))
         }
      }

      if (is.character(btt)) # turn btt=c("foo","bar") into list("foo","bar")
         btt <- as.list(btt)

      if (!is.list(btt))
         btt <- list(btt)

      spec <- btt
      btt <- lapply(btt, .set.btt, x$p, x$int.incl, colnames(x$X), fixed=fixed)

      if (x$intercept && !intercept && any(sapply(btt, function(bttj) length(bttj) == 1L && bttj == 1L)))
         stop(mstyle$stop("Cannot compute VIF(s) for the specified 'btt' argument."))

      ### get var-cov matrix of the fixed effects (location coefficients)

      vcov <- vcov(x, type="beta")

      ### compute (G)VIF for each element in the btt list

      obj <- if (reestimate) x else NULL

      res <- list()
      res$vif <- lapply(seq_along(btt), .compvif, btt=btt, vcov=vcov, xintercept=x$intercept, intercept=intercept, spec=spec, colnames=colnames(x$X), obj=obj, sim=FALSE)

      ### add coefficient names

      if (bttmiss) {
         names(res$vif) <- sapply(res$vif, function(x) x$coefname)
      } else {
         names(res$vif) <- sapply(res$vif, function(x) x$coefs)
      }

      ### add (G)VIFs as vector

      res$vifs <- sapply(res$vif, function(x) x$vif)

      ### add coefficient table if requested

      if (table && bttmiss) {
         res$table <- coef(summary(x), type="beta")
         res$test  <- x$test
      }

      res$bttspec <- !bttmiss
      res$digits  <- digits
      class(res)  <- "vif.rma"

      ######################################################################

      ### if sim >= 2, simulate corresponding (G)VIFs under independence

      sim.loc <- sim

      ### but skip this if all (G)VIFs are equal to 1

      if (all(sapply(res$vif, function(x) x$vif) == 1, na.rm=TRUE))
         sim.loc <- 0

      if (sim >= 2 && any(x$coef.na)) {
         warning(mstyle$warning("Cannot use 'sim' when some redundant predictors were dropped from the model."))
         sim.loc <- 0
      }

      if (sim.loc >= 2) {

         if (parallel == "no")
            vif.sim <- pbapply::pblapply(seq_len(sim.loc), .compvifsim, obj=x, coef="beta", btt=btt, att=NULL, reestimate=reestimate, intercept=intercept, parallel=parallel, seed=seed, joinb=joinb)

         if (parallel == "multicore")
            vif.sim <- pbapply::pblapply(seq_len(sim.loc), .compvifsim, obj=x, coef="beta", btt=btt, att=NULL, reestimate=reestimate, intercept=intercept, parallel=parallel, seed=seed, joinb=joinb, cl=ncpus)

         if (parallel == "snow") {

            if (.isTRUE(ddd$LB)) {
               vif.sim <- parallel::parLapplyLB(cl, seq_len(sim.loc), .compvifsim, obj=x, coef="beta", btt=btt, att=NULL, reestimate=reestimate, intercept=intercept, parallel=parallel, seed=seed, joinb=joinb)
            } else {
               vif.sim <- pbapply::pblapply(seq_len(sim.loc), .compvifsim, obj=x, coef="beta", btt=btt, att=NULL, reestimate=reestimate, intercept=intercept, parallel=parallel, seed=seed, joinb=joinb, cl=cl)
            }

         }

         vif.sim <- do.call(rbind, vif.sim)
         rownames(vif.sim) <- seq_len(sim.loc)
         colnames(vif.sim) <- seq_along(btt)

         if (!is.null(joinb) || is.null(x$data) || is.null(x$formula.mods)) {
            attr(vif.sim, "type") <- "X"
         } else {
            attr(vif.sim, "type") <- "data"
         }

         res$sim <- vif.sim

         vifs <- sapply(res$vif, function(x) x$vif)
         res$prop <- apply(vifs >= t(vif.sim), 1, mean, na.rm=TRUE)

      }

      ######################################################################

   } else {

      res <- NULL

   }

   #########################################################################

   if (vif.scale) {

      res.loc <- res

      ### process/set att argument

      if (attmiss) {
         if (x$Z.intercept && !intercept) {
            att <- as.list(2:x$q)
         } else {
            att <- as.list(seq_len(x$q))
         }
      }

      if (is.character(att))
         att <- as.list(att)

      if (!is.list(att))
         att <- list(att)

      spec <- att
      att <- lapply(att, .set.btt, x$q, x$Z.int.incl, colnames(x$Z), fixed=fixed)

      if (x$Z.intercept && !intercept && any(sapply(att, function(attj) length(attj) == 1L && attj == 1L)))
         stop(mstyle$stop("Cannot compute VIF(s) for the specified 'att' argument."))

      ### get var-cov matrix of the fixed effects (scale coefficients)

      vcov <- vcov(x, type="alpha")

      ### compute (G)VIF for each element in the att list

      obj <- if (reestimate) x else NULL

      res.scale <- list()
      res.scale$vif <- lapply(seq_along(att), .compvif, btt=att, vcov=vcov, xintercept=x$Z.intercept, intercept=intercept, spec=spec, colnames=colnames(x$Z), obj=obj, coef="alpha", sim=FALSE)

      ### add coefficient names

      if (attmiss) {
         names(res.scale$vif) <- sapply(res.scale$vif, function(x) x$coefname)
      } else {
         names(res.scale$vif) <- sapply(res.scale$vif, function(x) x$coefs)
      }

      ### add (G)VIFs as vector

      res.scale$vifs <- sapply(res.scale$vif, function(x) x$vif)

      ### add coefficient table if requested

      if (table && attmiss) {
         res.scale$table <- coef(summary(x), type="alpha")
         res.scale$test  <- x$test
      }

      res.scale$attspec <- !attmiss
      res.scale$digits  <- digits
      class(res.scale)  <- "vif.rma"

      ######################################################################

      ### if sim >= 2, simulate corresponding (G)VIFs under independence

      sim.scale <- sim

      ### but skip this if all (G)VIFs are equal to 1

      if (all(sapply(res.scale$vif, function(x) x$vif) == 1, na.rm=TRUE))
         sim.scale <- 0

      if (sim >= 2 && any(x$coef.na.Z)) {
         warning(mstyle$warning("Cannot use 'sim' when some redundant predictors were dropped from the model."))
         sim.scale <- 0
      }

      if (sim.scale >= 2) {

         if (parallel == "no")
            vif.sim <- pbapply::pblapply(seq_len(sim.scale), .compvifsim, obj=x, coef="alpha", btt=NULL, att=att, reestimate=reestimate, intercept=intercept, parallel=parallel, seed=seed, joina=joina)

         if (parallel == "multicore")
            vif.sim <- pbapply::pblapply(seq_len(sim.scale), .compvifsim, obj=x, coef="alpha", btt=NULL, att=att, reestimate=reestimate, intercept=intercept, parallel=parallel, seed=seed, joina=joina, cl=ncpus)

         if (parallel == "snow") {

            if (.isTRUE(ddd$LB)) {
               vif.sim <- parallel::parLapplyLB(cl, seq_len(sim.scale), .compvifsim, obj=x, coef="alpha", btt=NULL, att=att, reestimate=reestimate, intercept=intercept, parallel=parallel, seed=seed, joina=joina)
            } else {
               vif.sim <- pbapply::pblapply(seq_len(sim.scale), .compvifsim, obj=x, coef="alpha", btt=NULL, att=att, reestimate=reestimate, intercept=intercept, parallel=parallel, seed=seed, joina=joina, cl=cl)
            }

         }

         vif.sim <- do.call(rbind, vif.sim)
         rownames(vif.sim) <- seq_len(sim.scale)
         colnames(vif.sim) <- seq_along(att)

         if (!is.null(joina) || is.null(x$data) || is.null(x$formula.scale)) {
            attr(vif.sim, "type") <- "X"
         } else {
            attr(vif.sim, "type") <- "data"
         }

         res.scale$sim <- vif.sim

         vifs <- sapply(res.scale$vif, function(x) x$vif)
         res.scale$prop <- apply(vifs >= t(vif.sim), 1, mean, na.rm=TRUE)

      }

      ######################################################################

      if (vif.loc) {
         res <- list(beta=res.loc, alpha=res.scale)
         class(res) <- "vif.rma"
      } else {
         res <- res.scale
      }

   }

   #########################################################################

   if (.isTRUE(ddd$time)) {
      time.end <- proc.time()
      .print.time(unname(time.end - time.start)[3])
   }

   return(res)

}
