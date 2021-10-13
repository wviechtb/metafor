cooks.distance.rma.mv <- function(model, progbar=FALSE, cluster, reestimate=TRUE, parallel="no", ncpus=1, cl, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(model), must="rma.mv")

   #if (inherits(model, "robust.rma")) ### can compute Cook's distance also for 'robust.rma' objects
   #   stop(mstyle$stop("Method not available for objects of class \"robust.rma\"."))

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   x <- model

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

   misscluster <- ifelse(missing(cluster), TRUE, FALSE)

   if (misscluster)
      cluster <- seq_len(x$k.all)

   ddd <- list(...)

   .chkdots(ddd, c("btt", "time", "LB"))

   btt <- .set.btt(ddd$btt, x$p, int.incl=FALSE, Xnames=colnames(x$X))
   m <- length(btt)

   if (.isTRUE(ddd$time))
      time.start <- proc.time()

   #########################################################################

   ### process cluster variable

   # note: cluster variable is assumed to be of the same length as the size of
   # the original dataset passed to the model fitting function and so we apply
   # the same subsetting and removing of missings (if necessary) as was done
   # during model fitting

   if (!is.null(x$subset))
      cluster <- cluster[x$subset]

   cluster.f <- cluster

   cluster <- cluster[x$not.na]

   ### checks on cluster variable

   if (anyNA(cluster.f))
      stop(mstyle$stop("No missing values allowed in 'cluster' variable."))

   if (length(cluster.f) == 0L)
      stop(mstyle$stop(paste0("Cannot find 'cluster' variable (or it has zero length).")))

   if (length(cluster) != x$k)
      stop(mstyle$stop(paste0("Length of variable specified via 'cluster' (", length(cluster), ") does not match length of data (", x$k, ").")))

   ### cluster ids and number of clusters

   ids <- unique(cluster)
   n <- length(ids)

   #########################################################################

   ### calculate inverse of variance-covariance matrix under the full model

   svb <- chol2inv(chol(x$vb[btt,btt,drop=FALSE]))

   if (parallel == "no")
      res <- pbapply::pblapply(seq_len(n), .cooks.distance.rma.mv, obj=x, parallel=parallel, svb=svb, cluster=cluster, ids=ids, reestimate=reestimate, btt=btt)

   if (parallel == "multicore")
      res <- pbapply::pblapply(seq_len(n), .cooks.distance.rma.mv, obj=x, parallel=parallel, svb=svb, cluster=cluster, ids=ids, reestimate=reestimate, btt=btt, cl=ncpus)
      #res <- parallel::mclapply(seq_len(n), .cooks.distance.rma.mv, obj=x, parallel=parallel, svb=svb, cluster=cluster, ids=ids, reestimate=reestimate, btt=btt, mc.cores=ncpus)

   if (parallel == "snow") {
      if (is.null(cl)) {
         cl <- parallel::makePSOCKcluster(ncpus)
         on.exit(parallel::stopCluster(cl), add=TRUE)
      }
      if (.isTRUE(ddd$LB)) {
         res <- parallel::parLapplyLB(cl, seq_len(n), .cooks.distance.rma.mv, obj=x, parallel=parallel, svb=svb, cluster=cluster, ids=ids, reestimate=reestimate, btt=btt)
         #res <- parallel::clusterApplyLB(cl, seq_len(n), .cooks.distance.rma.mv, obj=x, parallel=parallel, svb=svb, cluster=cluster, ids=ids, reestimate=reestimate, btt=btt)
      } else {
         res <- pbapply::pblapply(seq_len(n), .cooks.distance.rma.mv, obj=x, parallel=parallel, svb=svb, cluster=cluster, ids=ids, reestimate=reestimate, btt=btt, cl=cl)
         #res <- parallel::parLapply(cl, seq_len(n), .cooks.distance.rma.mv, obj=x, parallel=parallel, svb=svb, cluster=cluster, ids=ids, reestimate=reestimate, btt=btt)
         #res <- parallel::clusterApply(cl, seq_len(n), .cooks.distance.rma.mv, obj=x, parallel=parallel, svb=svb, cluster=cluster, ids=ids, reestimate=reestimate, btt=btt)
      }
   }

   cook.d <- sapply(res, function(x) x$cook.d)

   #########################################################################

   if (na.act == "na.omit") {
      out <- cook.d
      if (misscluster) {
         names(out) <- x$slab[x$not.na]
      } else {
         names(out) <- ids
         out <- out[order(ids)]
      }
   }

   if (na.act == "na.exclude" || na.act == "na.pass") {

      ids.f <- unique(cluster.f)

      out <- rep(NA_real_, length(ids.f))
      out[match(ids, ids.f)] <- cook.d

      if (misscluster) {
         names(out) <- x$slab
      } else {
         names(out) <- ids.f
         out <- out[order(ids.f)]
      }

   }

   if (na.act == "na.fail" && any(!x$not.na))
      stop(mstyle$stop("Missing values in results."))

   if (.isTRUE(ddd$time)) {
      time.end <- proc.time()
      .print.time(unname(time.end - time.start)[3])
   }

   return(out)

}
