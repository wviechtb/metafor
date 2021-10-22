aggregate.escalc <- function(x, cluster, time, obs, V, struct="CS", rho, phi, weighted=TRUE, checkpd=TRUE, fun, na.rm=TRUE, subset, select, digits, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="escalc")

   if (any(!is.element(struct, c("ID","CS","CAR","CS+CAR","CS*CAR"))))
      stop(mstyle$stop("Unknown 'struct' specified."))

   if (missing(cluster))
      stop(mstyle$stop("Must specify 'cluster' variable."))

   if (length(na.rm) == 1L)
      na.rm <- c(na.rm, na.rm)

   k <- nrow(x)

   #########################################################################

   ### extract V, cluster, time, and subset variables

   mf <- match.call()

   V       <- .getx("V",       mf=mf, data=x)
   cluster <- .getx("cluster", mf=mf, data=x)
   time    <- .getx("time",    mf=mf, data=x)
   obs     <- .getx("obs",     mf=mf, data=x)
   subset  <- .getx("subset",  mf=mf, data=x)

   #########################################################################

   ### checks on cluster variable

   if (anyNA(cluster))
      stop(mstyle$stop("No missing values allowed in 'cluster' variable."))

   if (length(cluster) != k)
      stop(mstyle$stop(paste0("Length of variable specified via 'cluster' (", length(cluster), ") does not match length of data (", k, ").")))

   ucluster <- unique(cluster)
   n <- length(ucluster)

   #########################################################################

   ### get vi variable

   if (!is.null(attr(x, "vi.names"))) { # if vi.names attributes is available
      vi.name <- attr(x, "vi.names")[1] # take the first entry to be the vi variable
   } else {                             # if not, see if 'vi' is in the object and assume that is the vi variable
      if (!is.element("vi", names(x)))
         stop(mstyle$stop("Cannot determine name of the 'vi' variable."))
      vi.name <- "vi"
   }

   if (is.null(x[[vi.name]]))
      stop(mstyle$stop("Cannot find 'vi' variable in data frame."))

   #########################################################################

   if (is.null(V)) {

      ### if V is not specified

      vi <- x[[vi.name]]

      ### construct V matrix based on the specified structure

      if (struct=="ID")
         R <- diag(1, nrow=k, ncol=k)

      if (is.element(struct, c("CS","CS+CAR","CS*CAR"))) {

         if (missing(rho))
            stop(mstyle$stop("Must specify 'rho' for this var-cov structure."))

         if (length(rho) == 1L)
            rho <- rep(rho, n)

         if (length(rho) != n)
            stop(mstyle$stop(paste0("Length of 'rho' (", length(rho), ") does not match the number of clusters (", n, ").")))

         if (any(rho > 1) || any(rho < -1))
            stop(mstyle$stop("Value(s) of 'rho' must be in [-1,1]."))

      }

      if (is.element(struct, c("CAR","CS+CAR","CS*CAR"))) {

         if (missing(phi))
            stop(mstyle$stop("Must specify 'phi' for this var-cov structure."))

         if (length(phi) == 1L)
            phi <- rep(phi, n)

         if (length(phi) != n)
            stop(mstyle$stop(paste0("Length of 'phi' (", length(phi), ") does not match the number of clusters (", n, ").")))

         if (any(phi > 1) || any(phi < 0))
            stop(mstyle$stop("Value(s) of 'phi' must be in [0,1]."))

         ### checks on time variable

         if (!is.element("time", names(mf)))
            stop(mstyle$stop("Must specify 'time' variable for this var-cov structure."))

         if (length(time) != k)
            stop(mstyle$stop(paste0("Length of variable specified via 'time' (", length(time), ") does not match length of data (", k, ").")))

         if (struct == "CS*CAR") {

            ### checks on obs variable

         if (!is.element("obs", names(mf)))
            stop(mstyle$stop("Must specify 'obs' variable for this var-cov structure."))

            if (length(obs) != k)
               stop(mstyle$stop(paste0("Length of variable specified via 'obs' (", length(obs), ") does not match length of data (", k, ").")))

         }

      }

      if (struct=="CS") {

         R <- matrix(0, nrow=k, ncol=k)
         for (i in 1:n) {
            R[cluster == ucluster[i], cluster == ucluster[i]] <- rho[i]
         }

      }

      if (struct == "CAR") {

         R <- matrix(0, nrow=k, ncol=k)
         for (i in 1:n) {
            R[cluster == ucluster[i], cluster == ucluster[i]] <- outer(time[cluster == ucluster[i]], time[cluster == ucluster[i]], function(x,y) phi[i]^(abs(x-y)))
         }

      }

      if (struct == "CS+CAR") {

         R <- matrix(0, nrow=k, ncol=k)
         for (i in 1:n) {
            R[cluster == ucluster[i], cluster == ucluster[i]] <- rho[i] + (1 - rho[i]) * outer(time[cluster == ucluster[i]], time[cluster == ucluster[i]], function(x,y) phi[i]^(abs(x-y)))
         }

      }

      if (struct == "CS*CAR") {

         R <- matrix(0, nrow=k, ncol=k)
         for (i in 1:n) {
            R[cluster == ucluster[i], cluster == ucluster[i]] <- outer(obs[cluster == ucluster[i]], obs[cluster == ucluster[i]], function(x,y) ifelse(x==y, 1, rho[i])) * outer(time[cluster == ucluster[i]], time[cluster == ucluster[i]], function(x,y) phi[i]^(abs(x-y)))
         }

      }

      diag(R) <- 1
      S <- diag(sqrt(as.vector(vi)), nrow=k, ncol=k)
      V <- S %*% R %*% S

   } else {

      ### if V is specified

      if (.is.vector(V)) {

         if (length(V) == 1L)
            V <- rep(V, k)

         if (length(V) != k)
            stop(mstyle$stop(paste0("Length of 'V' (", length(V), ") does not match length of data frame (", k, ").")))

         V <- diag(as.vector(V), nrow=k, ncol=k)

      }

      if (is.data.frame(V))
         V <- as.matrix(V)

      if (!is.null(dimnames(V)))
         V <- unname(V)

      if (!.is.square(V))
         stop(mstyle$stop("'V' must be a square matrix."))

      if (!isSymmetric(V))
         stop(mstyle$stop("'V' must be a symmetric matrix."))

      if (nrow(V) != k)
         stop(mstyle$stop(paste0("Dimensions of 'V' (", nrow(V), "x", ncol(V), ") do not match length of data frame (", k, ").")))

      ### check that covariances are really 0 for estimates belonging to different clusters
      ### note: if na.rm[1] is FALSE, there may be missings in V, so skip check in those clusters

      for (i in 1:n) {
         if (any(abs(V[cluster == ucluster[i], cluster != ucluster[i]]) >= .Machine$double.eps, na.rm=TRUE))
            warning(mstyle$warning(paste0("Estimates in cluster '", ucluster[i], "' appear to have non-zero covariances with estimates belonging to different clusters.")), call.=FALSE)
      }

   }

   if (!is.null(attr(x, "yi.names"))) { # if yi.names attributes is available
      yi.name <- attr(x, "yi.names")[1] # take the first entry to be the yi variable
   } else {                             # if not, see if 'yi' is in the object and assume that is the yi variable
      if (!is.element("yi", names(x)))
         stop(mstyle$stop("Cannot determine name of the 'yi' variable."))
      yi.name <- "yi"
   }

   if (is.null(x[[yi.name]]))
      stop(mstyle$stop("Cannot find 'yi' variable in data frame."))

   ### note: there may be multiple yi/vi pairs; only first will be used

   yi <- as.vector(x[[yi.name]])

   ### if 'subset' is not null, apply subset

   if (!is.null(subset)) {

      subset <- .setnafalse(subset, k=k)

      x  <- x[subset,,drop=FALSE]
      yi <- yi[subset]
      V  <- V[subset,subset,drop=FALSE]
      cluster <- cluster[subset]

      k <- nrow(x)
      ucluster <- unique(cluster)
      n <- length(ucluster)

      if (k == 0L)
         stop(mstyle$stop("Processing terminated since k == 0."))

   }

   ### remove missings in yi/vi/V if na.rm[1] is TRUE

   if (na.rm[1]) {

      has.na <- is.na(yi) | .anyNAv(V)
      not.na <- !has.na

      if (any(has.na)) {

         x  <- x[not.na,]
         yi <- yi[not.na]
         V  <- V[not.na,not.na,drop=FALSE]
         cluster <- cluster[not.na]

      }

      k <- nrow(x)
      ucluster <- unique(cluster)
      n <- length(ucluster)

      if (k == 0L)
         stop(mstyle$stop("Processing terminated since k == 0."))

   }

   ### check that 'V' is positive definite (in each cluster)

   if (checkpd) {

      all.pd <- TRUE

      for (i in 1:n) {

         Vi <- V[cluster == ucluster[i], cluster == ucluster[i]]

         if (!anyNA(Vi) && any(eigen(Vi, symmetric=TRUE, only.values=TRUE)$values <= .Machine$double.eps)) {
            all.pd <- FALSE
            warning(mstyle$warning(paste0("'V' appears to be not positive definite in cluster ", ucluster[i], ".")), call.=FALSE)
         }

      }

      if (!all.pd)
         stop(mstyle$stop("Cannot aggregate estimates with a non-positive-definite 'V' matrix."))

   }

   ### compute aggregated estimates and corresponding sampling variances

   yi.agg <- rep(NA_real_, n)
   vi.agg <- rep(NA_real_, n)

   for (i in 1:n) {

      Vi <- V[cluster == ucluster[i], cluster == ucluster[i]]

      if (weighted) {

         Wi <- try(chol2inv(chol(Vi)), silent=TRUE)

         if (inherits(Wi, "try-error"))
            stop(mstyle$stop(paste0("Cannot take inverse of 'V' in cluster ", ucluster[i], ".")))

         sumWi <- sum(Wi)
         yi.agg[i] <- sum(Wi %*% cbind(yi[cluster == ucluster[i]])) / sumWi
         vi.agg[i] <- 1 / sumWi

      } else {

         ki <- sum(cluster == ucluster[i])
         yi.agg[i] <- sum(yi[cluster == ucluster[i]]) / ki
         vi.agg[i] <- sum(Vi) / ki^2

      }

   }

   if (!missing(fun)) {

      if (!is.list(fun) || length(fun) != 3 || any(sapply(fun, function(f) !is.function(f))))
         stop(mstyle$stop("Argument 'fun' must be a list of functions of length 3."))

      fun1 <- fun[[1]]
      fun2 <- fun[[2]]
      fun3 <- fun[[3]]

   } else {

      fun1 <- function(x) {
         m <- mean(x, na.rm=na.rm[2])
         if (is.nan(m)) NA else m
      }
      fun2 <- fun1
      fun3 <- function(x) {
         if (na.rm[2]) {
            tab <- table(na.omit(x))
            #tab <- table(x, useNA=ifelse(na.rm[2], "no", "ifany"))
         } else {
            tab <- table(x, useNA="ifany")
         }
         val <- tail(names(sort(tab)), 1)
         if (is.null(val)) NA else val
      }

   }

   ### turn 'cluster' into a factor with the desired levels, such that split() will give the same order

   fcluster <- factor(cluster, levels=ucluster)

   xsplit <- split(x, fcluster)

   xagg <- lapply(xsplit, function(xi) {
      tmp <- lapply(xi, function(xij) {
         if (inherits(xij, c("numeric","integer"))) {
            fun1(xij)
         } else if (inherits(xij, c("logical"))) {
            fun2(xij)
         } else {
            fun3(xij)
         }
      })
      as.data.frame(tmp)
   })

   xagg <- do.call(rbind, xagg)

   ### turn variables that were factors back into factors

   facs <- sapply(x, is.factor)

   if (any(facs)) {
      for (j in which(facs)) {
         xagg[[j]] <- factor(xagg[[j]])
      }
   }

   ### put yi.agg and vi.agg into the aggregate data at their respective positions

   xagg[which(names(xagg) == yi.name)] <- yi.agg
   xagg[which(names(xagg) == vi.name)] <- vi.agg

   ### add back some attributes

   measure <- attr(x[[yi.name]], "measure")

   if (is.null(measure))
      measure <- "GEN"

   attr(xagg[[yi.name]], "measure") <- measure

   attr(xagg, "yi.names") <- yi.name
   attr(xagg, "vi.names") <- vi.name

   if (!missing(digits)) {
      attr(xagg, "digits") <- .get.digits(digits=digits, xdigits=attr(x, "digits"), dmiss=FALSE)
   } else {
      attr(xagg, "digits") <- attr(x, "digits")
   }

   if (is.null(attr(xagg, "digits"))) ### in case x no longer has a 'digits' attribute
      attr(xagg, "digits") <- 4

   class(xagg) <- c("escalc", "data.frame")

   ### if 'select' is not missing, select variables to include in the output

   if (!missing(select)) {
      nl <- as.list(seq_along(x))
      names(nl) <- names(x)
      sel <- eval(substitute(select), nl, parent.frame())
      xagg <- xagg[,sel,drop=FALSE]
   }

   rownames(xagg) <- NULL

   return(xagg)

}
