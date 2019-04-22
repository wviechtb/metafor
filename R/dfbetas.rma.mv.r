dfbetas.rma.mv <- function(model, progbar=FALSE, cluster, reestimate=TRUE, parallel="no", ncpus=1, cl=NULL, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(model, "rma.mv"))
      stop(mstyle$stop("Argument 'model' must be an object of class \"rma.mv\"."))

   if (inherits(model, "robust.rma"))
      stop(mstyle$stop("Method not available for objects of class \"robust.rma\"."))

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   x <- model

   parallel <- match.arg(parallel, c("no", "snow", "multicore"))

   if (parallel == "no" && ncpus > 1)
      parallel <- "snow"

   misscluster <- ifelse(missing(cluster), TRUE, FALSE)

   if (misscluster)
      cluster <- seq_len(x$k.all)

   #########################################################################

   ### process cluster variable
   ### note: cluster variable is assumed to be of the same length as the original data passed to the model fitting function
   ###       so we have to apply the same subsetting (if necessary) and removing of missings as done during model fitting

   if (!is.null(x$subset))
      cluster <- cluster[x$subset]

   cluster.f <- cluster

   cluster <- cluster[x$not.na]

   ### checks on cluster variable

   if (anyNA(cluster.f))
      stop(mstyle$stop("No missing values allowed in 'cluster' variable."))

   if (length(cluster) != x$k)
      stop(mstyle$stop("Length of variable specified via 'cluster' does not match length of data."))

   ### cluster ids and number of clusters

   ids <- unique(cluster)
   n <- length(ids)

   #########################################################################

   if (parallel=="no") {

      dfbs <- matrix(NA_real_, nrow=n, ncol=x$p)

      if (progbar)
         pbar <- txtProgressBar(min=0, max=n, style=3)

      for (i in seq_len(n)) {

         if (progbar)
            setTxtProgressBar(pbar, i)

         incl <- cluster %in% ids[i]

         if (reestimate) {

            ### set initial values to estimates from full model

            control             <- x$control
            control$sigma2.init <- x$sigma2
            control$tau2.init   <- x$tau2
            control$rho.init    <- x$rho
            control$gamma2.init <- x$gamma2
            control$phi.init    <- x$phi

            ### fit model without data from ith cluster

            res <- try(suppressWarnings(rma.mv(x$yi, V=x$V, W=x$W, mods=x$X, random=x$random, struct=x$struct, intercept=FALSE, data=x$mf.r, method=x$method, test=x$test, level=x$level, R=x$R, Rscale=x$Rscale, sigma2=ifelse(x$vc.fix$sigma2, x$sigma2, NA), tau2=ifelse(x$vc.fix$tau2, x$tau2, NA), rho=ifelse(x$vc.fix$rho, x$rho, NA), gamma2=ifelse(x$vc.fix$gamma2, x$gamma2, NA), phi=ifelse(x$vc.fix$phi, x$phi, NA), sparse=x$sparse, dist=x$dist, control=control, subset=!incl)), silent=TRUE)

         } else {

            ### set values of variance/correlation components to those from the 'full' model

            res <- try(suppressWarnings(rma.mv(x$yi, V=x$V, W=x$W, mods=x$X, random=x$random, struct=x$struct, intercept=FALSE, data=x$mf.r, method=x$method, test=x$test, level=x$level, R=x$R, Rscale=x$Rscale, sigma2=x$sigma2, tau2=x$tau2, rho=x$rho, gamma2=x$gamma2, phi=x$phi, sparse=x$sparse, dist=x$dist, control=x$control, subset=!incl)), silent=TRUE)

         }

         if (inherits(res, "try-error"))
            next

         ### removing a cluster could lead to a model coefficient becoming inestimable

         if (any(res$coef.na))
            next

         ### fit model based on all data but with var/cor components fixed to those from res

         tmp <- try(suppressWarnings(rma.mv(x$yi, V=x$V, W=x$W, mods=x$X, random=x$random, struct=x$struct, intercept=FALSE, data=x$mf.r, method=x$method, test=x$test, level=x$level, R=x$R, Rscale=x$Rscale, sigma2=res$sigma2, tau2=res$tau2, rho=res$rho, gamma2=res$gamma2, phi=res$phi, sparse=x$sparse, dist=x$dist, control=x$control)), silent=TRUE)

         ### compute dfbeta value(s)

         dfb <- x$beta - res$beta

         ### compute dfbetas

         dfbs[i,] <- dfb / sqrt(diag(tmp$vb))

      }

      if (progbar)
         close(pbar)

   }

   if (parallel=="snow" || parallel == "multicore") {

      if (!requireNamespace("parallel", quietly=TRUE))
         stop(mstyle$stop("Please install the 'parallel' package for parallel processing."))

      ncpus <- as.integer(ncpus)

      if (ncpus < 1)
         stop(mstyle$stop("Argument 'ncpus' must be >= 1."))

      if (parallel == "multicore")
         res <- parallel::mclapply(seq_len(n), .dfbetas.rma.mv, obj=x, mc.cores=ncpus, parallel=parallel, cluster=cluster, ids=ids, reestimate=reestimate)

      if (parallel == "snow") {
         if (is.null(cl)) {
            cl <- parallel::makePSOCKcluster(ncpus)
            on.exit(parallel::stopCluster(cl))
         }
         res <- parallel::parLapply(cl, seq_len(n), .dfbetas.rma.mv, obj=x, parallel=parallel, cluster=cluster, ids=ids, reestimate=reestimate)
      }

      dfbs <- lapply(res, function(x) x$dfbs)
      dfbs <- do.call(rbind, dfbs)

   }

   #########################################################################

   if (na.act == "na.omit") {
      out <- dfbs
      if (misscluster) {
         rownames(out) <- x$slab[x$not.na]
      } else {
         rownames(out) <- ids
         out <- out[order(ids),]
      }
   }

   if (na.act == "na.exclude" || na.act == "na.pass") {

      ids.f <- unique(cluster.f)

      out <- matrix(NA_real_, nrow=length(ids.f), ncol=x$p)
      out[match(ids, ids.f),] <- dfbs

      if (misscluster) {
         rownames(out) <- x$slab
      } else {
         rownames(out) <- ids.f
         out <- out[order(ids.f),]
      }

   }

   if (na.act == "na.fail" && any(!x$not.na))
      stop(mstyle$stop("Missing values in results."))

   colnames(out) <- rownames(x$beta)

   out <- data.frame(out)

   return(out)

}
