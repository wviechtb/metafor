rstandard.rma.mv <- function(model, digits, cluster, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(model, "rma.mv"))
      stop(mstyle$stop("Argument 'model' must be an object of class \"rma.mv\"."))

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   x <- model

   if (missing(digits))
      digits <- x$digits

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

   #########################################################################

   options(na.action="na.omit")
   H <- hatvalues(x, type="matrix")
   options(na.action = na.act)

   #########################################################################

   ImH <- diag(x$k) - H
   #ei <- ImH %*% cbind(x$yi)
   ei <- c(x$yi - x$X %*% x$beta)

   ei[abs(ei) < 100 * .Machine$double.eps] <- 0
   #ei[abs(ei) < 100 * .Machine$double.eps * median(abs(ei), na.rm=TRUE)] <- 0 ### see lm.influence

   if (inherits(x, "robust.rma")) {
      ve <- ImH %*% tcrossprod(x$meat,ImH)
   } else {
      ve <- ImH %*% tcrossprod(x$M,ImH)
   }

   #ve <- x$M + x$X %*% x$vb %*% t(x$X) - 2*H%*%x$M
   sei <- sqrt(diag(ve))

   #########################################################################

   if (!misscluster) {

      ### cluster ids and number of clusters

      ids <- unique(cluster)
      n <- length(ids)

      X2 <- rep(NA_real_, n)
      k.id <- rep(NA_integer_, n)

      for (i in seq_len(n)) {

         incl <- cluster %in% ids[i]
         k.id[i] <- sum(incl)

         sve <- try(chol2inv(chol(as.matrix(ve[incl,incl,drop=FALSE]))), silent=TRUE)
         #sve <- try(solve(ve[incl,incl,drop=FALSE]), silent=TRUE)

         if (inherits(sve, "try-error"))
            next

         X2[i] <- rbind(ei[incl]) %*% sve %*% cbind(ei[incl])

      }

   }

   #########################################################################

   resid   <- rep(NA_real_, x$k.f)
   seresid <- rep(NA_real_, x$k.f)
   stresid <- rep(NA_real_, x$k.f)

   resid[x$not.na]   <- ei
   seresid[x$not.na] <- sei
   stresid[x$not.na] <- ei / sei

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
