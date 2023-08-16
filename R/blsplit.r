blsplit <- function(x, cluster, fun, args, sort=FALSE) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (missing(cluster))
      stop(mstyle$stop("Must specify 'cluster' variable."))

   if (!is.matrix(x) && !inherits(x, "dgCMatrix"))
      stop(mstyle$stop("Argument 'x' must be a matrix."))

   if (!isSymmetric(x))
      stop(mstyle$stop("Argument 'x' must be a symmetric matrix."))

   k <- nrow(x)

   if (length(cluster) != k)
      stop(mstyle$stop(paste0("Length of variable specified via 'cluster' (", length(cluster), ") does not correspond to the dimensions of the matrix (", k, "x", k, ").")))

   res <- list()

   clusters <- unique(cluster)

   if (sort)
      clusters <- sort(clusters)

   for (i in seq_along(clusters)) {
      res[[i]] <- x[cluster == clusters[i], cluster == clusters[i], drop=FALSE]
   }

   names(res) <- clusters

   if (!missing(fun)) {
      if (missing(args)) {
         res <- lapply(res, fun)
      } else {
         args <- as.list(args)
         for (i in 1:length(res)) {
            res[[i]] <- do.call(fun, c(unname(res[i]), args))
         }
      }
   }

   return(res)

}
