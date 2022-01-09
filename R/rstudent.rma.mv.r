rstudent.rma.mv <- function(model, digits, progbar=FALSE, cluster, reestimate=TRUE, parallel="no", ncpus=1, cl, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(model), must="rma.mv", notav="robust.rma")

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
      on.exit(pbapply::pboptions(pbo), add=TRUE)
   }

   if (missing(digits)) {
      digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
   }

   misscluster <- ifelse(missing(cluster), TRUE, FALSE)

   if (misscluster) {
      cluster <- seq_len(x$k.all)
   } else {
      mf <- match.call()
      cluster <- .getx("cluster", mf=mf, data=x$data)
   }

   ddd <- list(...)

   .chkdots(ddd, c("time", "LB"))

   if (.isTRUE(ddd$time))
      time.start <- proc.time()

   #########################################################################

   ### process cluster variable

   # note: cluster variable is assumed to be of the same length as the size of
   # the original dataset passed to the model fitting function and so we apply
   # the same subsetting and removing of missings (if necessary) as was done
   # during model fitting

   if (length(cluster) != x$k.all)
      stop(mstyle$stop(paste0("Length of variable specified via 'cluster' (", length(cluster), ") does not match length of data (", x$k.all, ").")))

   if (!is.null(x$subset))
      cluster <- cluster[x$subset]

   cluster.f <- cluster

   cluster <- cluster[x$not.na]

   ### checks on cluster variable

   if (anyNA(cluster.f))
      stop(mstyle$stop("No missing values allowed in 'cluster' variable."))

   if (length(cluster.f) == 0L)
      stop(mstyle$stop(paste0("Cannot find 'cluster' variable (or it has zero length).")))

   ### cluster ids and number of clusters

   ids <- unique(cluster)
   n <- length(ids)

   #########################################################################

   if (parallel == "no")
      res <- pbapply::pblapply(seq_len(n), .rstudent.rma.mv, obj=x, parallel=parallel, cluster=cluster, ids=ids, reestimate=reestimate)

   if (parallel == "multicore")
      res <- pbapply::pblapply(seq_len(n), .rstudent.rma.mv, obj=x, parallel=parallel, cluster=cluster, ids=ids, reestimate=reestimate, cl=ncpus)
      #res <- parallel::mclapply(seq_len(n), .rstudent.rma.mv, obj=x, parallel=parallel, cluster=cluster, ids=ids, reestimate=reestimate, mc.cores=ncpus)

   if (parallel == "snow") {
      if (is.null(cl)) {
         cl <- parallel::makePSOCKcluster(ncpus)
         on.exit(parallel::stopCluster(cl), add=TRUE)
      }
      if (.isTRUE(ddd$LB)) {
         res <- parallel::parLapplyLB(cl, seq_len(n), .rstudent.rma.mv, obj=x, parallel=parallel, cluster=cluster, ids=ids, reestimate=reestimate)
         #res <- parallel::clusterApplyLB(cl, seq_len(n), .rstudent.rma.mv, obj=x, parallel=parallel, cluster=cluster, ids=ids, reestimate=reestimate)
      } else {
         res <- pbapply::pblapply(seq_len(n), .rstudent.rma.mv, obj=x, parallel=parallel, cluster=cluster, ids=ids, reestimate=reestimate, cl=cl)
         #res <- parallel::parLapply(cl, seq_len(n), .rstudent.rma.mv, obj=x, parallel=parallel, cluster=cluster, ids=ids, reestimate=reestimate)
         #res <- parallel::clusterApply(cl, seq_len(n), .rstudent.rma.mv, obj=x, parallel=parallel, cluster=cluster, ids=ids, reestimate=reestimate)
      }
   }

   delresid   <- rep(NA_real_, x$k)
   sedelresid <- rep(NA_real_, x$k)

   pos <- unlist(sapply(res, function(x) x$pos))

   delresid[pos]   <- unlist(sapply(res, function(x) x$delresid))
   sedelresid[pos] <- unlist(sapply(res, function(x) x$sedelresid))

   X2   <- sapply(res, function(x) x$X2)
   k.id <- sapply(res, function(x) x$k.id)

   #########################################################################

   delresid[abs(delresid) < 100 * .Machine$double.eps] <- 0

   resid   <- rep(NA_real_, x$k.f)
   seresid <- rep(NA_real_, x$k.f)

   resid[x$not.na]   <- delresid
   seresid[x$not.na] <- sedelresid

   stresid <- resid / seresid

   #########################################################################

   if (na.act == "na.omit") {
      out <- list(resid=resid[x$not.na], se=seresid[x$not.na], z=stresid[x$not.na])
      if (!misscluster)
         out$cluster <- cluster.f[x$not.na]
      out$slab <- x$slab[x$not.na]
   }

   if (na.act == "na.exclude" || na.act == "na.pass") {
      out <- list(resid=resid, se=seresid, z=stresid)
      if (!misscluster)
         out$cluster <- cluster.f
      out$slab <- x$slab
   }

   if (na.act == "na.fail" && any(!x$not.na))
      stop(mstyle$stop("Missing values in results."))

   if (.isTRUE(ddd$time)) {
      time.end <- proc.time()
      .print.time(unname(time.end - time.start)[3])
   }

   if (misscluster) {

      out$digits <- digits

      class(out) <- "list.rma"
      return(out)

   } else {

      out <- list(out)

      if (na.act == "na.omit") {
         out[[2]] <- list(X2=X2[order(ids)], k=k.id[order(ids)], slab=ids[order(ids)])
      }

      if (na.act == "na.exclude" || na.act == "na.pass") {

         ids.f <- unique(cluster.f)

         X2.f <- rep(NA_real_, length(ids.f))
         X2.f[match(ids, ids.f)] <- X2

         k.id.f <- sapply(ids.f, function(id) sum((id == cluster.f) & x$not.na))

         out[[2]] <- list(X2=X2.f[order(ids.f)], k=k.id.f[order(ids.f)], slab=ids.f[order(ids.f)])

      }

      out[[1]]$digits <- digits
      out[[2]]$digits <- digits

      names(out) <- c("obs", "cluster")

      class(out[[1]]) <- "list.rma"
      class(out[[2]]) <- "list.rma"
      attr(out[[1]], ".rmspace") <- TRUE
      attr(out[[2]], ".rmspace") <- TRUE
      return(out)

   }

}
