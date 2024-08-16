vcalc <- function(vi, cluster, subgroup, obs, type, time1, time2, grp1, grp2, w1, w2,
data, rho, phi, rvars, checkpd=TRUE, nearpd=FALSE, sparse=FALSE, ..., new = TRUE) {

   mstyle <- .get.mstyle()

   ############################################################################

   if (missing(vi))
      stop(mstyle$stop("Must specify 'vi' variable."))

   if (missing(cluster))
      stop(mstyle$stop("Must specify 'cluster' variable."))

   ### get ... argument and check for extra/superfluous arguments

   ddd <- list(...)

   .chkdots(ddd, c("nearPD", "retdat"))

   if (.isTRUE(ddd$nearPD))
      nearpd <- TRUE

   ### check if data argument has been specified

   if (missing(data))
      data <- NULL

   if (is.null(data) && !missing(rvars))
      stop(mstyle$stop("Must specify 'data' argument when using 'rvars'."))

   if (is.null(data)) {
      data <- sys.frame(sys.parent())
   } else {
      if (!is.data.frame(data))
         data <- data.frame(data)
   }

   subgroup.spec <- !missing(subgroup)
   type.spec     <- !missing(type)
   obs.spec      <- !missing(obs)
   grp1.spec     <- !missing(grp1)
   grp2.spec     <- !missing(grp2)
   time1.spec    <- !missing(time1)
   time2.spec    <- !missing(time2)
   w1.spec       <- !missing(w1)
   w2.spec       <- !missing(w2)

   mf <- match.call()

   vi       <- .getx("vi",       mf=mf, data=data, checknumeric=TRUE)
   cluster  <- .getx("cluster",  mf=mf, data=data)
   subgroup <- .getx("subgroup", mf=mf, data=data)
   type     <- .getx("type",     mf=mf, data=data)
   obs      <- .getx("obs",      mf=mf, data=data)
   grp1     <- .getx("grp1",     mf=mf, data=data)
   grp2     <- .getx("grp2",     mf=mf, data=data)
   time1    <- .getx("time1",    mf=mf, data=data, checknumeric=TRUE)
   time2    <- .getx("time2",    mf=mf, data=data, checknumeric=TRUE)
   w1       <- .getx("w1",       mf=mf, data=data, checknumeric=TRUE)
   w2       <- .getx("w2",       mf=mf, data=data, checknumeric=TRUE)

   ############################################################################

   # to be able to quickly set vi to a constant (e.g., 1) for all rows

   if (length(vi) == 1L && length(cluster) > 1L)
      vi <- rep(vi, length(cluster))

   k <- length(vi)

   if (k == 1L)
      stop(mstyle$stop("Processing terminated since k = 1.")) # could also do: return(matrix(vi, nrow=1, ncol=1))

   #########################################################################

   ### checks on cluster variable

   if (anyNA(cluster))
      stop(mstyle$stop("No missing values allowed in 'cluster' variable."))

   if (length(cluster) != k)
      stop(mstyle$stop(paste0("Length of variable specified via 'cluster' (", length(cluster), ") does not match length of 'vi' (", k, ").")))

   ### checks on subgroup variable

   if (subgroup.spec) {

      if (anyNA(subgroup))
         stop(mstyle$stop("No missing values allowed in 'subgroup' variable."))

      if (length(subgroup) != k)
         stop(mstyle$stop(paste0("Length of variable specified via 'subgroup' (", length(subgroup), ") does not match length of 'vi' (", k, ").")))

      cluster <- paste0(cluster, ".", subgroup)

   }

   ucluster <- unique(cluster)
   n <- length(ucluster)

   #########################################################################

   if (missing(rvars)) {

      ############################################################################

      ### process type variable

      if (type.spec) {
         if (missing(rho))
            stop(mstyle$stop("Must specify 'rho' when 'type' is specified."))
      } else {
         type <- rep(1, k)
      }

      if (anyNA(type))
         stop(mstyle$stop("No missing values allowed in 'type' variable."))

      if (length(type) != k)
         stop(mstyle$stop(paste0("Length of variable specified via 'type' (", length(type), ") does not match length of 'vi' (", k, ").")))

      ### process obs variable

      if (obs.spec) {
         if (missing(rho))
            stop(mstyle$stop("Must specify 'rho' when 'obs' is specified."))
      } else {
         #obs <- ave(cluster, cluster, FUN=seq_along)
         obs <- rep(1, k)
      }

      if (anyNA(obs))
         stop(mstyle$stop("No missing values allowed in 'obs' variable."))

      if (length(obs) != k)
         stop(mstyle$stop(paste0("Length of variable specified via 'obs' (", length(obs), ") does not match length of 'vi' (", k, ").")))

      ### process grp1 and grp2 variables

      #if ((grp1.spec && !grp2.spec) || (!grp1.spec && grp2.spec))
      #   stop(mstyle$stop("Either specify both 'grp1' and 'grp2' or neither."))

      if ((grp2.spec && !grp1.spec))
         stop(mstyle$stop("Either specify only 'grp1', both 'grp1' and 'grp2', or neither."))

      if (!grp1.spec)
         grp1 <- rep(1, k)

      if (!grp2.spec)
         grp2 <- rep(2, k)

      if (anyNA(grp1))
         stop(mstyle$stop("No missing values allowed in 'grp1' variable."))

      if (anyNA(grp2))
         stop(mstyle$stop("No missing values allowed in 'grp2' variable."))

      if (length(grp1) != k)
         stop(mstyle$stop(paste0("Length of variable specified via 'grp1' (", length(grp1), ") does not match length of 'vi' (", k, ").")))

      if (length(grp2) != k)
         stop(mstyle$stop(paste0("Length of variable specified via 'grp2' (", length(grp2), ") does not match length of 'vi' (", k, ").")))

      ### process time1 and time2 variables

      if ((time2.spec && !time1.spec))
         stop(mstyle$stop("Either specify only 'time1', both 'time1' and 'time2', or neither."))

      if (time2.spec && !grp2.spec)
         stop(mstyle$stop("Must specify 'grp2' when 'time2' is specified."))

      if (!time1.spec)
         time1 <- rep(1, k)

      if (!time2.spec)
         time2 <- time1

      if (time1.spec || time2.spec) {
         if (missing(phi))
            stop(mstyle$stop("Must specify 'phi' when 'time1' and/or 'time2' is specified."))
      } else {
         phi <- 1
      }

      if (abs(phi) > 1)
         stop(mstyle$stop("Value of argument 'phi' must be in [-1,1]."))

      if (anyNA(time1))
         stop(mstyle$stop("No missing values allowed in 'time1' variable."))

      if (anyNA(time2))
         stop(mstyle$stop("No missing values allowed in 'time2' variable."))

      if (length(time1) != k)
         stop(mstyle$stop(paste0("Length of variable specified via 'time1' (", length(time1), ") does not match length of 'vi' (", k, ").")))

      if (length(time2) != k)
         stop(mstyle$stop(paste0("Length of variable specified via 'time2' (", length(time2), ") does not match length of 'vi' (", k, ").")))

      if (!is.numeric(time1))
         stop(mstyle$stop("Variable 'time1' must be a numeric variable."))

      if (!is.numeric(time2))
         stop(mstyle$stop("Variable 'time2' must be a numeric variable."))

      ### process w1 and w2 variables

      if ((w2.spec && !w1.spec))
         stop(mstyle$stop("Either specify only 'w1', both 'w1' and 'w2', or neither."))

      if (w2.spec && !grp2.spec)
         stop(mstyle$stop("Must specify 'grp2' when 'w2' is specified."))

      if (!w1.spec)
         w1 <- rep(1, k)

      if (!w2.spec)
         w2 <- w1

      if (anyNA(w1))
         stop(mstyle$stop("No missing values allowed in 'w1' variable."))

      if (anyNA(w2))
         stop(mstyle$stop("No missing values allowed in 'w2' variable."))

      if (length(w1) != k)
         stop(mstyle$stop(paste0("Length of variable specified via 'w1' (", length(w1), ") does not match length of 'vi' (", k, ").")))

      if (length(w2) != k)
         stop(mstyle$stop(paste0("Length of variable specified via 'w2' (", length(w2), ") does not match length of 'vi' (", k, ").")))

      if (!is.numeric(w1))
         stop(mstyle$stop("Variable 'w1' must be a numeric variable."))

      if (!is.numeric(w2))
         stop(mstyle$stop("Variable 'w2' must be a numeric variable."))

      ############################################################################

      ### process/create rho

      if (!missing(rho) && !(.is.vector(rho) || is.matrix(rho) || is.list(rho)))
         stop(mstyle$stop("Argument 'rho' must either be a vector, a matrix, or a list."))

      if (type.spec) {
         if (obs.spec) {
            # both type and obs are specified
            if (.is.vector(rho)) {
               if (length(rho) != 2L)
                  stop(mstyle$stop("When 'type' and 'obs' are both specified, 'rho' must specify both the within- and between-construct correlations."))
               rho <- as.list(rho)
            } else {
               if (is.matrix(rho)) {
                  stop(mstyle$stop("When 'type' and 'obs' are both specified, 'rho' must specify both the within- and between-construct correlations."))
               } else {
                  if (length(rho) != 2L)
                     stop(mstyle$stop("When 'type' and 'obs' are both specified and 'rho' is a list, then it must have two elements."))
               }
            }
         } else {
            # only type is specified
            if (.is.vector(rho)) {
               if (length(rho) != 1L)
                  stop(mstyle$stop("When only 'type' is specified, 'rho' must be a scalar."))
               rho <- list(0, rho)
            } else {
               if (is.matrix(rho)) {
                  rho <- list(0, rho)
               } else {
                  if (length(rho) != 1L)
                     stop(mstyle$stop("When only 'type' is specified, 'rho' must have a single list element."))
                  rho <- list(0, rho[[1]])
               }
            }
         }
      } else {
         if (obs.spec) {
            # only obs is specified
            if (.is.vector(rho)) {
               if (length(rho) != 1L)
                  stop(mstyle$stop("When only 'obs' is specified, 'rho' must be a scalar."))
               rho <- list(rho, 0)
            } else {
               if (is.matrix(rho)) {
                  rho <- list(rho, 0)
               } else {
                  if (length(rho) != 1L)
                     stop(mstyle$stop("When only 'obs' is specified, 'rho' must have a single list element."))
                  rho <- list(rho[[1]], 0)
               }
            }
         } else {
            # neither type nor obs is specified
            rho <- list(0, 0)
         }
      }

      if (length(rho[[1]]) == 1L) {
         rho[[1]] <- matrix(rho[[1]], nrow=length(unique(obs)), ncol=length(unique(obs)))
         diag(rho[[1]]) <- 1
         rownames(rho[[1]]) <- colnames(rho[[1]]) <- unique(obs)
      }
      if (length(rho[[2]]) == 1L) {
         rho[[2]] <- matrix(rho[[2]], nrow=length(unique(type)), ncol=length(unique(type)))
         diag(rho[[2]]) <- 1
         rownames(rho[[2]]) <- colnames(rho[[2]]) <- unique(type)
      }

      if (any(!sapply(rho, .is.square)))
         stop(mstyle$stop("All matrices specified via 'rho' argument must be square matrices."))

      if (any(abs(rho[[1]]) > 1) || any(abs(rho[[2]]) > 1))
         stop(mstyle$stop("All correlations specified via 'rho' must be in [-1,1]."))

      if (is.null(dimnames(rho[[1]])) || is.null(dimnames(rho[[2]])))
         stop(mstyle$stop("Any matrices specified via 'rho' must have dimension names."))

      if (is.null(rownames(rho[[1]])))
         rownames(rho[[1]]) <- colnames(rho[[1]])
      if (is.null(rownames(rho[[2]])))
         rownames(rho[[2]]) <- colnames(rho[[2]])
      if (is.null(colnames(rho[[1]])))
         colnames(rho[[1]]) <- rownames(rho[[1]])
      if (is.null(colnames(rho[[2]])))
         colnames(rho[[2]]) <- rownames(rho[[2]])

      if (!all(unique(obs) %in% rownames(rho[[1]])))
         stop("There are 'obs' values with no corresponding row/column in the correlation matrix.")
      if (!all(unique(type) %in% rownames(rho[[2]])))
         stop("There are 'type' values with no corresponding row/column in the correlation matrix.")

      #return(rho)

      ############################################################################

      #### turn obs and type into character variables to that [obs[i],obs[j]] and [type[i],type[j]] below work correctly

      obs  <- as.character(obs)
      type <- as.character(type)

      ### construct R matrix

      if (sparse) {
         R <- Matrix(0, nrow=k, ncol=k)
      } else {
         R <- matrix(0, nrow=k, ncol=k)
      }

      if (new) {
        cluster_set <- unique(cluster) 
        
        for (cl in cluster_set) {
          cl_i <- which(cl == cluster)
          k_c <- length(cl_i)
          R_c <- matrix(0, nrow = k_c, ncol = k_c)
          diag(R_c) <- 1
          if (k_c > 1L) {
            for (i in 2:k_c) {
              for (j in 1:i) {
                ci <- cl_i[i]
                cj <- cl_i[j]
                R_c[i,j] <- ifelse(type[ci]==type[cj], ifelse(obs[ci]==obs[cj], 1, rho[[1]][obs[ci],obs[cj]]), rho[[2]][type[ci],type[cj]]) *
                  (ifelse(grp1[ci]==grp1[cj], ifelse(time1[ci]==time1[cj], 1, phi^abs(time1[ci]-time1[cj])), 0) * sqrt(1/w1[ci] * 1/w1[cj]) -
                     ifelse(grp1[ci]==grp2[cj], ifelse(time1[ci]==time2[cj], 1, phi^abs(time1[ci]-time2[cj])), 0) * sqrt(1/w1[ci] * 1/w2[cj]) -
                     ifelse(grp2[ci]==grp1[cj], ifelse(time2[ci]==time1[cj], 1, phi^abs(time2[ci]-time1[cj])), 0) * sqrt(1/w2[ci] * 1/w1[cj]) +
                     ifelse(grp2[ci]==grp2[cj], ifelse(time2[ci]==time2[cj], 1, phi^abs(time2[ci]-time2[cj])), 0) * sqrt(1/w2[ci] * 1/w2[cj])) /
                  (sqrt(1/w1[ci] + 1/w2[ci] - 2*ifelse(grp1[ci]==grp2[ci], ifelse(time1[ci]==time2[ci], 1, phi^abs(time1[ci]-time2[ci])), 0) * sqrt(1/w1[ci] * 1/w2[ci])) *
                     sqrt(1/w1[cj] + 1/w2[cj] - 2*ifelse(grp1[cj]==grp2[cj], ifelse(time1[cj]==time2[cj], 1, phi^abs(time1[cj]-time2[cj])), 0) * sqrt(1/w1[cj] * 1/w2[cj])))
              }
            }
          }
          R_c[upper.tri(R_c)] <- t(R_c)[upper.tri(R_c)]
          R[cl_i, cl_i] <- R_c
        }
        
      } else {
        diag(R) <- 1
        
        for (i in 2:k) {
          for (j in 1:i) {
            if (cluster[i] == cluster[j]) {
              
              R[i,j] <- ifelse(type[i]==type[j], ifelse(obs[i]==obs[j], 1, rho[[1]][obs[i],obs[j]]), rho[[2]][type[i],type[j]]) *
                (ifelse(grp1[i]==grp1[j], ifelse(time1[i]==time1[j], 1, phi^abs(time1[i]-time1[j])), 0) * sqrt(1/w1[i] * 1/w1[j]) -
                   ifelse(grp1[i]==grp2[j], ifelse(time1[i]==time2[j], 1, phi^abs(time1[i]-time2[j])), 0) * sqrt(1/w1[i] * 1/w2[j]) -
                   ifelse(grp2[i]==grp1[j], ifelse(time2[i]==time1[j], 1, phi^abs(time2[i]-time1[j])), 0) * sqrt(1/w2[i] * 1/w1[j]) +
                   ifelse(grp2[i]==grp2[j], ifelse(time2[i]==time2[j], 1, phi^abs(time2[i]-time2[j])), 0) * sqrt(1/w2[i] * 1/w2[j])) /
                (sqrt(1/w1[i] + 1/w2[i] - 2*ifelse(grp1[i]==grp2[i], ifelse(time1[i]==time2[i], 1, phi^abs(time1[i]-time2[i])), 0) * sqrt(1/w1[i] * 1/w2[i])) *
                   sqrt(1/w1[j] + 1/w2[j] - 2*ifelse(grp1[j]==grp2[j], ifelse(time1[j]==time2[j], 1, phi^abs(time1[j]-time2[j])), 0) * sqrt(1/w1[j] * 1/w2[j])))
              
            }
          }
        }
        
        R[upper.tri(R)] <- t(R)[upper.tri(R)]
        
      }

   } else {

      ### when rvars are specified

      ### warn user if non-relevant arguments have been specified

      not.miss <- c(type.spec, obs.spec, grp1.spec, grp2.spec, time1.spec, time2.spec, w1.spec, w2.spec, !missing(rho), !missing(phi))

      if (any(not.miss)) {
         args <- c("type", "obs", "grp1", "grp2", "time1", "time2", "w1", "w2", "rho", "phi")
         warning(mstyle$warning("Argument", ifelse(sum(not.miss) > 1, "s", ""), " '", paste0(args[not.miss], collapse=","), "' ignored for when 'rvars' is specified."), call.=FALSE)
      }

      ### get position of rvars in data

      nl <- as.list(seq_along(data))
      names(nl) <- names(data)
      rvars <- try(eval(substitute(rvars), envir=nl, enclos=NULL), silent=TRUE)

      if (inherits(rvars, "try-error"))
         stop(mstyle$stop("Could not find all variables specified via 'rvars' in 'data'."))

      ### get rvars from data

      has.colon <- grepl(":", deparse1(substitute(rvars)), fixed=TRUE)

      if (has.colon && length(rvars) == 2L) {
         rvars <- data[seq(from = rvars[1], to = rvars[2])]
      } else {
         rvars <- data[rvars]
      }

      ### check that number of rvars makes sense given the k per cluster

      k.cluster <- tapply(cluster, cluster, length)

      if (max(k.cluster) > length(rvars))
         stop(mstyle$stop(paste0("There ", ifelse(length(rvars) == 1L, "is 1 variable ", paste0("are ", length(rvars), " variables ")), "specified via 'rvars', but there are clusters with more rows.")))

      if (max(k.cluster) != length(rvars))
         stop(mstyle$stop(paste0("There ", ifelse(length(rvars) == 1L, "is 1 variable ", paste0("are ", length(rvars), " variables ")), "specified via 'rvars', but no cluster with this many rows.")))

      ### construct R matrix based on rvars

      R <- list()

      for (i in seq_len(n)) {
         x <- rvars[cluster == ucluster[i],]
         x <- x[seq_len(nrow(x))]
         if (anyNA(x[lower.tri(x, diag=TRUE)]))
            warning(mstyle$warning(paste0("There are missing values in 'rvals' for cluster ", ucluster[i], ".")), call.=FALSE)
         x[upper.tri(x)] <- t(x)[upper.tri(x)]
         R[[i]] <- as.matrix(x)
      }

      names(R) <- ucluster

      #R <- lapply(split(rvars, cluster), function(x) {
      #   k <- nrow(x)
      #   x <- x[seq_len(k)]
      #   x[upper.tri(x)] <- t(x)[upper.tri(x)]
      #   as.matrix(x)
      #   })

      #R <- bldiag(R, order=cluster)

      R <- bldiag(R)
      R <- Matrix(R, sparse=TRUE)

   }

   #return(R)

   ############################################################################

   ### check that 'R' is positive definite in each cluster

   if (checkpd || nearpd) {

      for (i in seq_len(n)) {

         Ri <- R[cluster == ucluster[i], cluster == ucluster[i]]

         if (!anyNA(Ri) && !.chkpd(Ri)) {

            if (nearpd) {
               Ri <- try(as.matrix(nearPD(Ri, corr=TRUE)$mat), silent=TRUE)
               if (inherits(Ri, "try-error")) {
                  warning(mstyle$warning(paste0("Using nearPD() failed in cluster ", ucluster[i], ".")), call.=FALSE)
               } else {
                  if (!anyNA(Ri) && !.chkpd(Ri))
                     warning(mstyle$warning(paste0("The var-cov matrix still appears to be not positive definite in cluster ", ucluster[i], " even after nearPD().")), call.=FALSE)
                  R[cluster == ucluster[i], cluster == ucluster[i]] <- Ri
               }
            } else {
               warning(mstyle$warning(paste0("The var-cov matrix appears to be not positive definite in cluster ", ucluster[i], ".")), call.=FALSE)
            }

         }

      }

   }

   ############################################################################

   ### turn R into V

   vi <- as.vector(vi)

   S <- Diagonal(k, sqrt(vi))
   V <- S %*% R %*% S

   if (!sparse)
      V <- as.matrix(V)

   if (.isTRUE(ddd$retdat))
      V <- data.frame(cluster, type, obs, grp1, grp2, time1, time2, w1, w2, vi, V=V)

   return(V)

}
