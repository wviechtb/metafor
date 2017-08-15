rstudent.rma.mv <- function(model, digits, progbar=FALSE, cluster, reestimate=TRUE, parallel="no", ncpus=1, cl=NULL, ...) {

   if (!inherits(model, "rma.mv"))
      stop("Argument 'model' must be an object of class \"rma.mv\".")

   if (inherits(model, "robust.rma"))
      stop("Method not yet implemented for objects of class \"robust.rma\". Sorry!")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop("Unknown 'na.action' specified under options().")

   x <- model

   parallel <- match.arg(parallel, c("no", "snow", "multicore"))

   if (missing(digits))
      digits <- x$digits

   misscluster <- ifelse(missing(cluster), TRUE, FALSE)

   if (misscluster)
      cluster <- 1:x$k.all

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
      stop("No missing values allowed in 'cluster' variable.")

   if (length(cluster) != x$k)
      stop("Length of variable specified via 'cluster' does not match length of data.")

   ### cluster ids and number of clusters

   ids <- unique(cluster)
   n <- length(ids)

   #########################################################################

   if (parallel=="no") {

      delresid   <- rep(NA_real_, x$k)
      sedelresid <- rep(NA_real_, x$k)
      X2   <- rep(NA_integer_, n)
      k.id <- rep(NA_integer_, n)

      if (progbar)
         pbar <- txtProgressBar(min=0, max=n, style=3)

      for (i in seq_len(n)) {

         if (progbar)
            setTxtProgressBar(pbar, i)

         incl <- cluster %in% ids[i]

         k.id[i] <- sum(incl)

         if (reestimate) {

            ### set initial values to estimates from full model

            control             <- x$control
            control$sigma2.init <- x$sigma2
            control$tau2.init   <- x$tau2
            control$rho.init    <- x$rho
            control$gamma2.init <- x$gamma2
            control$phi.init    <- x$phi

            ### fit model without data from ith cluster

            res <- try(suppressWarnings(rma.mv(x$yi, V=x$V, W=x$W, mods=x$X, random=x$random, struct=x$struct, intercept=FALSE, data=x$mf.r, method=x$method, test=x$test, level=x$level, R=x$R, Rscale=x$Rscale, sigma2=ifelse(x$vc.fix$sigma2, x$sigma2, NA), tau2=ifelse(x$vc.fix$tau2, x$tau2, NA), rho=ifelse(x$vc.fix$rho, x$rho, NA), gamma2=ifelse(x$vc.fix$gamma2, x$gamma2, NA), phi=ifelse(x$vc.fix$phi, x$phi, NA), sparse=x$sparse, control=control, subset=!incl)), silent=TRUE)

         } else {

            ### set values of variance/correlation components to those from the 'full' model

            res <- try(suppressWarnings(rma.mv(x$yi, V=x$V, W=x$W, mods=x$X, random=x$random, struct=x$struct, intercept=FALSE, data=x$mf.r, method=x$method, test=x$test, level=x$level, R=x$R, Rscale=x$Rscale, sigma2=x$sigma2, tau2=x$tau2, rho=x$rho, gamma2=x$gamma2, phi=x$phi, sparse=x$sparse, control=x$control, subset=!incl)), silent=TRUE)

         }

         if (inherits(res, "try-error"))
            next

         ### removing a cluster could lead to a model coefficient becoming inestimable

         if (any(res$coef.na))
            next

         ### fit model based on all data but with var/cor components fixed to those from res

         tmp <- try(suppressWarnings(rma.mv(x$yi, V=x$V, W=x$W, mods=x$X, random=x$random, struct=x$struct, intercept=FALSE, data=x$mf.r, method=x$method, test=x$test, level=x$level, R=x$R, Rscale=x$Rscale, sigma2=res$sigma2, tau2=res$tau2, rho=res$rho, gamma2=res$gamma2, phi=res$phi, sparse=x$sparse, control=x$control)), silent=TRUE)

         Xi <- x$X[incl,,drop=FALSE]
         delpred  <- Xi %*% res$beta
         vdelpred <- Xi %*% res$vb %*% t(Xi)
         delresid[incl] <- x$yi[incl] - delpred
         sedelresid[incl] <- sqrt(diag(tmp$M[incl,incl,drop=FALSE] + vdelpred))

         sve <- try(chol2inv(chol(tmp$M[incl,incl,drop=FALSE] + vdelpred)), silent=TRUE)
         #sve <- try(solve(tmp$M[incl,incl,drop=FALSE] + vdelpred), silent=TRUE)

         if (inherits(sve, "try-error"))
            next

         X2[i] <- rbind(delresid[incl]) %*% sve %*% cbind(delresid[incl])

      }

      if (progbar)
         close(pbar)

   }

   if (parallel=="snow" || parallel == "multicore") {

      if (!requireNamespace("parallel", quietly=TRUE))
         stop("Please install the 'parallel' package for parallel processing.")

      ncpus <- as.integer(ncpus)

      if (ncpus < 1)
         stop("Argument 'ncpus' must be >= 1.")

      if (parallel == "multicore")
         res <- parallel::mclapply(seq_len(n), .rstudent.rma.mv, obj=x, mc.cores=ncpus, parallel=parallel, cluster=cluster, ids=ids, reestimate=reestimate)

      if (parallel == "snow") {
         if (is.null(cl)) {
            cl <- parallel::makePSOCKcluster(ncpus)
            res <- parallel::parLapply(cl, seq_len(n), .rstudent.rma.mv, obj=x, parallel=parallel, cluster=cluster, ids=ids, reestimate=reestimate)
            parallel::stopCluster(cl)
         } else {
            res <- parallel::parLapply(cl, seq_len(n), .rstudent.rma.mv, obj=x, parallel=parallel, cluster=cluster, ids=ids, reestimate=reestimate)
         }
      }

      delresid   <- rep(NA_real_, x$k)
      sedelresid <- rep(NA_real_, x$k)

      pos <- unlist(sapply(res, function(z) z$pos))

      delresid[pos]   <- unlist(sapply(res, function(z) z$delresid))
      sedelresid[pos] <- unlist(sapply(res, function(z) z$sedelresid))

      X2   <- sapply(res, function(z) z$X2)
      k.id <- sapply(res, function(z) z$k.id)

   }

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
      stop("Missing values in results.")

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
      return(out)

   }

}
