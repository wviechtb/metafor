gosh.rma <- function(x, subsets, progbar=TRUE, parallel="no", ncpus=1, cl=NULL, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(x, "rma"))
      stop(mstyle$stop("Argument 'x' must be an object of class \"rma\"."))

   if (inherits(x, "rma.glmm"))
      stop(mstyle$stop("Method not available for objects of class \"rma.glmm\"."))

   if (inherits(x, "rma.mv"))
      stop(mstyle$stop("Method not available for objects of class \"rma.mv\"."))

   if (inherits(x, "robust.rma"))
      stop(mstyle$stop("Method not available for objects of class \"robust.rma\"."))

   if (inherits(x, "rma.ls"))
      stop(mstyle$stop("Method not available for objects of class \"rma.ls\"."))

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   if (x$k == 1)
      stop(mstyle$stop("Stopped because k = 1."))

   parallel <- match.arg(parallel, c("no", "snow", "multicore"))

   if (parallel == "no" && ncpus > 1)
      parallel <- "snow"

   if (!is.null(cl) && inherits(cl, "SOCKcluster")) {
      parallel <- "snow"
      ncpus <- length(cl)
   }

   if (parallel == "snow" && ncpus < 2)
      parallel <- "no"

   if (parallel == "snow" || parallel == "multicore") {

      if (!requireNamespace("parallel", quietly=TRUE))
         stop(mstyle$stop("Please install the 'parallel' package for parallel processing."))

      ncpus <- as.integer(ncpus)

      if (ncpus < 1L)
         stop(mstyle$stop("Argument 'ncpus' must be >= 1."))

   }

   if (!progbar) {
      pbo <- pbapply::pboptions(type="none")
      on.exit(pbapply::pboptions(pbo))
   }

   ddd <- list(...)

   .chkdots(ddd, c("seed", "time", "LB"))

   if (.isTRUE(ddd$time))
      time.start <- proc.time()

   ### total number of possible subsets

   N.tot <- sum(choose(x$k, x$p:x$k))

   ### if 'subsets' is missing, include all possible subsets if N.tot is <= 10^6
   ### and otherwise include 10^6 random subsets; if the user specifies 'subsets'
   ### and N.tot <= subsets, then again include all possible subsets

   if (missing(subsets)) {
      if (N.tot <= 10^6) {
         exact <- TRUE
      } else {
         exact <- FALSE
         N.tot <- 10^6
      }
   } else {
      if (N.tot <= subsets) {
         exact <- TRUE
      } else {
         exact <- FALSE
         N.tot <- subsets
      }
   }

   if (N.tot == Inf)
      stop(mstyle$stop("Too many iterations required for all combinations."))

   if (progbar)
      message(paste0("Fitting ", N.tot, " models (based on ", ifelse(exact, "all possible", "random"), " subsets)."))

   #########################################################################

   ### generate inclusion matrix (either exact or at random)

   if (exact) {

      incl <- as.matrix(expand.grid(replicate(x$k, list(c(FALSE,TRUE))), KEEP.OUT.ATTRS=FALSE))
      incl <- incl[rowSums(incl) >= x$p,,drop=FALSE]

      ### slower, but does not generate rows that need to be filtered out (as above)
      #incl <- lapply(x$p:x$k, function(m) apply(combn(x$k,m), 2, function(l) 1:x$k %in% l))
      #incl <- t(do.call(cbind, incl))

   } else {

      if (!is.null(ddd$seed))
         set.seed(ddd$seed)

      j <- sample(x$p:x$k, N.tot, replace=TRUE, prob=dbinom(x$p:x$k, x$k, 0.5))
      incl <- t(sapply(j, function(m) seq_len(x$k) %in% sample(x$k, m)))

   }

   colnames(incl) <- seq_len(x$k)

   ### check if model is a standard FE model or a standard RE model with the DL estimators

   model <- 0L
   if (x$method=="FE" && x$weighted && is.null(x$weights) && x$int.only)
      model <- 1L
   if (x$method=="DL" && x$weighted && is.null(x$weights) && x$int.only)
      model <- 2L

   #########################################################################

   outlist <- "beta=beta, k=k, QE=QE, I2=I2, H2=H2, tau2=tau2, coef.na=coef.na"

   if (parallel == "no") {

      if (inherits(x, "rma.uni"))
         res <- pbapply::pbapply(incl, 1, .profile.rma.uni, obj=x, parallel=parallel, subset=TRUE, model=model, outlist=outlist)

      if (inherits(x, "rma.mh"))
         res <- pbapply::pbapply(incl, 1, .profile.rma.mh, obj=x, parallel=parallel, subset=TRUE, outlist=outlist)

      if (inherits(x, "rma.peto"))
         res <- pbapply::pbapply(incl, 1, .profile.rma.peto, obj=x, parallel=parallel, subset=TRUE, outlist=outlist)

   }

   if (parallel == "multicore") {

      if (inherits(x, "rma.uni"))
         res <- pbapply::pbapply(incl, 1, .profile.rma.uni, obj=x, parallel=parallel, subset=TRUE, model=model, outlist=outlist, cl=ncpus)
         #res <- parallel::mclapply(asplit(incl, 1), .profile.rma.uni, obj=x, mc.cores=ncpus, parallel=parallel, subset=TRUE, model=model, outlist=outlist)

      if (inherits(x, "rma.mh"))
         res <- pbapply::pbapply(incl, 1, .profile.rma.mh, obj=x, parallel=parallel, subset=TRUE, outlist=outlist, cl=ncpus)
         #res <- parallel::mclapply(asplit(incl, 1), .profile.rma.mh, obj=x, mc.cores=ncpus, parallel=parallel, subset=TRUE, outlist=outlist)

      if (inherits(x, "rma.peto"))
         res <- pbapply::pbapply(incl, 1, .profile.rma.peto, obj=x, parallel=parallel, subset=TRUE, outlist=outlist, cl=ncpus)
         #res <- parallel::mclapply(asplit(incl, 1), .profile.rma.peto, obj=x, mc.cores=ncpus, parallel=parallel, subset=TRUE, outlist=outlist)

   }

   if (parallel == "snow") {

      if (is.null(cl)) {
         cl <- parallel::makePSOCKcluster(ncpus)
         on.exit(parallel::stopCluster(cl), add=TRUE)
      }

      if (inherits(x, "rma.uni")) {
         if (.isTRUE(ddd$LB)) {
            res <- parallel::parLapplyLB(cl, asplit(incl, 1), .profile.rma.uni, obj=x, parallel=parallel, subset=TRUE, model=model, outlist=outlist)
         } else {
            res <- pbapply::pbapply(incl, 1, .profile.rma.uni, obj=x, parallel=parallel, subset=TRUE, model=model, outlist=outlist, cl=cl)
            #res <- parallel::parLapply(cl, asplit(incl, 1), .profile.rma.uni, obj=x, parallel=parallel, subset=TRUE, model=model, outlist=outlist)
         }
      }

      if (inherits(x, "rma.mh")) {
         if (.isTRUE(ddd$LB)) {
            res <- parallel::parLapplyLB(cl, asplit(incl, 1), .profile.rma.mh, obj=x, parallel=parallel, subset=TRUE, outlist=outlist)
         } else {
            res <- pbapply::pbapply(incl, 1, .profile.rma.mh, obj=x, parallel=parallel, subset=TRUE, outlist=outlist, cl=cl)
            #res <- parallel::parLapply(cl, asplit(incl, 1), .profile.rma.mh, obj=x, parallel=parallel, subset=TRUE, outlist=outlist)
         }
      }

      if (inherits(x, "rma.peto")) {
         if (.isTRUE(ddd$LB)) {
            res <- parallel::parLapplyLB(cl, asplit(incl, 1), .profile.rma.peto, obj=x, parallel=parallel, subset=TRUE, outlist=outlist)
         } else {
            res <- pbapply::pbapply(incl, 1, .profile.rma.peto, obj=x, parallel=parallel, subset=TRUE, outlist=outlist, cl=cl)
            #res <- parallel::parLapply(cl, asplit(incl, 1), .profile.rma.peto, obj=x, parallel=parallel, subset=TRUE, outlist=outlist)
         }
      }

   }

   beta <- do.call("rbind", lapply(res, function(x) if (inherits(x, "try-error") || any(x$coef.na)) NA else t(x$beta)))
   het  <- do.call("rbind", lapply(res, function(x) if (inherits(x, "try-error") || any(x$coef.na)) NA else c(x$k, x$QE, x$I2, x$H2, x$tau2)))

   if (all(is.na(het)))
      stop(mstyle$stop("All model fits failed."))

   #########################################################################

   ### in case a model fit was skipped, this guarantees that we still get
   ### a value for k in the first column of the het matrix for each model

   het[,1] <- rowSums(incl)

   ### set column names

   colnames(het) <- c("k", "QE", "I2", "H2", "tau2")

   if (x$int.only) {
      colnames(beta) <- "estimate"
   } else {
      colnames(beta) <- colnames(x$X)
   }

   ### combine het and beta objects and order incl and res by k

   res <- data.frame(het, beta)
   incl <- incl[order(res$k),,drop=FALSE]
   res <- res[order(res$k),,drop=FALSE]

   ### fix rownames

   rownames(res)  <- seq_len(nrow(res))
   rownames(incl) <- seq_len(nrow(incl))

   ### was model fitted successfully / all values are not NA?

   fit <- apply(res, 1, function(x) all(!is.na(x)))

   ### print processing time

   if (.isTRUE(ddd$time)) {
      time.end <- proc.time()
      .print.time(unname(time.end - time.start)[3])
   }

   ### list to return

   out <- list(res=res, incl=incl, fit=fit, k=x$k, int.only=x$int.only, method=x$method, measure=x$measure, digits=x$digits)

   class(out) <- "gosh.rma"
   return(out)

}
