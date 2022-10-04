############################################################################

### function to test for missings in a var-cov matrix

.anyNAv <- function(x) {
   k <- nrow(x)
   not.na <- not.na.diag <- !is.na(diag(x))
   for (i in seq_len(k)[not.na.diag]) {
      not.na[i] <- !anyNA(x[i, seq_len(k)[not.na.diag]])
   }
   return(!not.na)
}

### function to test each row for any missings in the lower triangular part of a matrix

#.anyNAv <- function(x)
#   return(sapply(seq_len(nrow(x)), FUN=function(i) anyNA(x[i,seq_len(i)])))

### function above is faster (and does not require making a copy of the object)

#.anyNAv <- function(X) {
#   X[upper.tri(X)] <- 0
#   return(apply(is.na(X), 1, any))
#}

############################################################################

.process.G.aftersub <- function(mf.g, struct, formula, tau2, rho, isG, k, sparse, verbose) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (verbose > 1)
      message(mstyle$message(paste0("Processing '", paste0(formula, collapse=""), "' term (#1) ...")))

   ### number of variables in model frame

   nvars <- ncol(mf.g)

   ### check that the number of variables is correct for the chosen structure

   if (is.element(struct, c("CS","HCS","UN","UNR","AR","HAR","CAR","ID","DIAG","PHYBM","PHYPL","PHYPD")) && nvars != 2)
      stop(mstyle$stop(paste0("Only a single inner variable allowed for an '~ inner | outer' term when 'struct=\"", struct, "\"'.")), call.=FALSE)

   ### get variables names in mf.g

   g.names <- names(mf.g) ### names for inner and outer factors/variables

   ### check that inner variable is a factor (or character variable) for structures that require this (no longer required)

   #if (is.element(struct, c("CS","HCS","UN","UNR","ID","DIAG")) && !is.factor(mf.g[[1]]) && !is.character(mf.g[[1]]))
   #   stop(mstyle$stop(paste0("Inner variable in '~ inner | outer' term must be a factor or character variable when 'struct=\"", struct, "\"'.")), call.=FALSE)

   ### for struct="CAR", check that inner term is numeric and get the unique numeric values

   if (is.element(struct, c("CAR"))) {
      if (!is.numeric(mf.g[[1]]))
         stop(mstyle$stop("Inner variable in '~ inner | outer' term must be numeric for 'struct=\"CAR\"'."), call.=FALSE)
      g.values <- sort(unique(round(mf.g[[1]], digits=8))) ### aweful hack to avoid floating points issues
   } else {
      g.values <- NULL
   }

   ### turn each variable in mf.g into a factor (not for SP/PHY structures or GEN)
   ### if a variable was a factor to begin with, this drops any unused levels, but order of existing levels is preserved

   if (is.element(struct, c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD","GEN","GDIAG"))) {
      mf.g <- data.frame(mf.g[-nvars], outer=factor(mf.g[[nvars]]))
   } else {
      mf.g <- data.frame(inner=factor(mf.g[[1]]), outer=factor(mf.g[[2]]))
   }

   ### check if there are any NAs anywhere in mf.g

   if (anyNA(mf.g))
      stop(mstyle$stop("No NAs allowed in variables specified via the 'random' argument."), call.=FALSE)

   ### get number of levels of each variable in mf.g (vector with two values, for the inner and outer factor)

   #g.nlevels <- c(nlevels(mf.g[[1]]), nlevels(mf.g[[2]])) ### works only for factors
   if (is.element(struct, c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD","GEN","GDIAG"))) {
      g.nlevels <- c(length(unique(apply(mf.g[-nvars], 1, paste, collapse=" + "))), length(unique(mf.g[[nvars]])))
   } else {
      g.nlevels <- c(length(unique(mf.g[[1]])), length(unique(mf.g[[2]])))
   }

   ### get levels of each variable in mf.g

   #g.levels <- list(levels(mf.g[[1]]), levels(mf.g[[2]])) ### works only for factors
   if (is.element(struct, c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD","GEN","GDIAG"))) {
      g.levels <- list(sort(unique(apply(mf.g[-nvars], 1, paste, collapse=" + "))), sort(unique((mf.g[[nvars]]))))
   } else {
      #g.levels <- list(sort(unique(as.character(mf.g[[1]]))), sort(unique(as.character(mf.g[[2]]))))
      g.levels <- list(as.character(sort(unique(mf.g[[1]]))), as.character(sort(unique(mf.g[[2]]))))
   }

   ### determine appropriate number of tau2 and rho values (note: this is done *after* subsetting)
   ### note: if g.nlevels[1] is 1, then technically there is no correlation, but we still need one
   ### rho for the optimization function (this rho is fixed to 0 further in the rma.mv() function)

   if (is.element(struct, c("CS","ID","AR","CAR","SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD"))) {
      tau2s <- 1
      rhos  <- 1
   }
   if (is.element(struct, c("HCS","DIAG","HAR"))) {
      tau2s <- g.nlevels[1]
      rhos  <- 1
   }
   if (struct == "UN") {
      tau2s <- g.nlevels[1]
      rhos  <- ifelse(g.nlevels[1] > 1, g.nlevels[1]*(g.nlevels[1]-1)/2, 1)
   }
   if (struct == "UNR") {
      tau2s <- 1
      rhos  <- ifelse(g.nlevels[1] > 1, g.nlevels[1]*(g.nlevels[1]-1)/2, 1)
   }
   if (struct == "GEN") {
      p     <- nvars - 1
      tau2s <- p
      rhos  <- ifelse(p > 1, p*(p-1)/2, 1)
   }
   if (struct == "GDIAG") {
      p     <- nvars - 1
      tau2s <- p
      rhos  <- 1
   }

   ### set default value(s) for tau2 if it is unspecified

   if (is.null(tau2))
      tau2 <- rep(NA_real_, tau2s)

   ### set default value(s) for rho argument if it is unspecified

   if (is.null(rho))
      rho <- rep(NA_real_, rhos)

   ### allow quickly setting all tau2 values to a fixed value

   if (length(tau2) == 1L)
      tau2 <- rep(tau2, tau2s)

   ### allow quickly setting all rho values to a fixed value

   if (length(rho) == 1L)
      rho <- rep(rho, rhos)

   ### check if tau2 and rho are of correct length

   if (length(tau2) != tau2s)
      stop(mstyle$stop(paste0("Length of ", ifelse(isG, 'tau2', 'gamma2'), " argument (", length(tau2), ") does not match actual number of variance components (", tau2s, ").")), call.=FALSE)
   if (length(rho) != rhos)
      stop(mstyle$stop(paste0("Length of ", ifelse(isG, 'rho', 'phi'), " argument (", length(rho), ") does not match actual number of correlations (", rhos, ").")), call.=FALSE)

   ### checks on any fixed values of tau2 and rho arguments

   if (any(tau2 < 0, na.rm=TRUE))
      stop(mstyle$stop(paste0("Specified value(s) of ", ifelse(isG, 'tau2', 'gamma2'), " must be >= 0.")), call.=FALSE)
   if (is.element(struct, c("CAR")) && any(rho > 1 | rho < 0, na.rm=TRUE))
      stop(mstyle$stop(paste0("Specified value(s) of ", ifelse(isG, 'rho', 'phi'), " must be in [0,1].")), call.=FALSE)
   if (is.element(struct, c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYPL","PHYPD")) && any(rho < 0, na.rm=TRUE))
      stop(mstyle$stop(paste0("Specified value(s) of ", ifelse(isG, 'rho', 'phi'), " must be >= 0.")), call.=FALSE)
   if (!is.element(struct, c("CAR","SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD")) && any(rho > 1 | rho < -1, na.rm=TRUE))
      stop(mstyle$stop(paste0("Specified value(s) of ", ifelse(isG, 'rho', 'phi'), " must be in [-1,1].")), call.=FALSE)

   ### create model matrix for inner and outer factors of mf.g

   if (is.element(struct, c("CS","HCS","UN","UNR","AR","HAR","CAR","ID","DIAG"))) {

      if (g.nlevels[1] == 1) {
         Z.G1 <- cbind(rep(1,k))
      } else {
         if (sparse) {
            Z.G1 <- sparse.model.matrix(~ mf.g[[1]] - 1)
         } else {
            Z.G1 <- model.matrix(~ mf.g[[1]] - 1)
         }
      }

   }

   if (is.element(struct, c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD"))) {

      if (sparse) {
         Z.G1 <- Diagonal(k)
      } else {
         Z.G1 <- diag(1, nrow=k, ncol=k)
      }

   }

   if (is.element(struct, c("GEN","GDIAG"))) {

      if (sparse) {
         Z.G1 <- Matrix(as.matrix(mf.g[-nvars]), sparse=TRUE)
      } else {
         Z.G1 <- as.matrix(mf.g[-nvars])
      }

   }

   if (g.nlevels[2] == 1) {
      Z.G2 <- cbind(rep(1,k))
   } else {
      if (sparse) {
         Z.G2 <- sparse.model.matrix(~ mf.g[[nvars]] - 1)
      } else {
         Z.G2 <- model.matrix(~ mf.g[[nvars]] - 1)
      }
   }

   attr(Z.G1, "assign")    <- NULL
   attr(Z.G1, "contrasts") <- NULL
   attr(Z.G2, "assign")    <- NULL
   attr(Z.G2, "contrasts") <- NULL

   return(list(mf.g=mf.g, g.names=g.names, g.nlevels=g.nlevels,
               g.levels=g.levels, g.values=g.values,
               tau2s=tau2s, rhos=rhos, tau2=tau2, rho=rho, Z.G1=Z.G1, Z.G2=Z.G2))

}

############################################################################

.process.G.afterrmna <- function(mf.g, g.nlevels, g.levels, g.values, struct, formula, tau2, rho, Z.G1, Z.G2, isG, sparse, distspec, verbose) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (verbose > 1)
      message(mstyle$message(paste0("Processing '", paste0(formula, collapse=""), "' term (#2) ...")))

   ### number of variables in model frame

   nvars <- ncol(mf.g)

   ### copy g.nlevels and g.levels

   g.nlevels.f <- g.nlevels
   g.levels.f  <- g.levels

   ### redo: turn each variable in mf.g into a factor (not for SP structures or GEN)
   ### (reevaluates the levels present, but order of existing levels is preserved)

   if (is.element(struct, c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD","GEN","GDIAG"))) {
      mf.g <- data.frame(mf.g[-nvars], outer=factor(mf.g[[nvars]]))
   } else {
      mf.g <- data.frame(inner=factor(mf.g[[1]]), outer=factor(mf.g[[2]]))
   }

   ### redo: get number of levels of each variable in mf.g (vector with two values, for the inner and outer factor)

   #g.nlevels <- c(nlevels(mf.g[[1]]), nlevels(mf.g[[2]])) ### works only for factors
   if (is.element(struct, c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD","GEN","GDIAG"))) {
      g.nlevels <- c(length(unique(apply(mf.g[-nvars], 1, paste, collapse=" + "))), length(unique(mf.g[[nvars]])))
   } else {
      g.nlevels <- c(length(unique(mf.g[[1]])), length(unique(mf.g[[2]])))
   }

   ### redo: get levels of each variable in mf.g

   #g.levels <- list(levels(mf.g[[1]]), levels(mf.g[[2]])) ### works only for factors
   if (is.element(struct, c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD","GEN","GDIAG"))) {
      g.levels <- list(sort(unique(apply(mf.g[-nvars], 1, paste, collapse=" + "))), sort(unique((mf.g[[nvars]]))))
   } else {
      #g.levels <- list(sort(unique(as.character(mf.g[[1]]))), sort(unique(as.character(mf.g[[2]]))))
      g.levels <- list(as.character(sort(unique(mf.g[[1]]))), as.character(sort(unique(mf.g[[2]]))))
   }

   ### determine which levels of the inner factor were removed

   g.levels.r <- !is.element(g.levels.f[[1]], g.levels[[1]])

   ### warn if any levels were removed (not for "AR","CAR","SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","GEN","GDIAG")

   if (any(g.levels.r) && !is.element(struct, c("AR","CAR","SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","GEN","GDIAG")))
      warning(mstyle$warning(paste0("One or more levels of inner factor (i.e., ", paste(g.levels.f[[1]][g.levels.r], collapse=", "), ") removed due to NAs.")), call.=FALSE)

   ### for "ID", "DIAG", and "GDIAG", fix rho to 0

   if (is.element(struct, c("ID","DIAG","GDIAG")))
      rho <- 0

   ### if there is only a single arm for "CS","HCS","AR","HAR","CAR" (either to begin with or after removing NAs), then fix rho to 0

   if (g.nlevels[1] == 1 && is.element(struct, c("CS","HCS","AR","HAR","CAR")) && is.na(rho)) {
      rho <- 0
      warning(mstyle$warning(paste0("Inner factor has only a single level, so fixed value of ", ifelse(isG, 'rho', 'phi'), " to 0.")), call.=FALSE)
   }

   ### if there is only a single arm for SP/PHY structures or GEN/GDIAG (either to begin with or after removing NAs), cannot fit model

   if (g.nlevels[1] == 1 && is.element(struct, c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD","GEN","GDIAG")))
      stop(mstyle$stop("Cannot fit model since inner term only has a single level."), call.=FALSE)

   ### k per level of the inner factor

   if (is.element(struct, c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD","GEN","GDIAG"))) {
      g.levels.k <- table(factor(apply(mf.g[-nvars], 1, paste, collapse=" + "), levels=g.levels.f[[1]]))
   } else {
      g.levels.k <- table(factor(mf.g[[1]], levels=g.levels.f[[1]]))
   }

   ### for "HCS","UN","DIAG","HAR": if a particular level of the inner factor only occurs once, then set corresponding tau2 value to 0 (if not already fixed)
   ### note: no longer done; variance component should still be (weakly) identifiable

   #if (is.element(struct, c("HCS","UN","DIAG","HAR"))) {
   #   if (any(is.na(tau2) & g.levels.k == 1)) {
   #      tau2[is.na(tau2) & g.levels.k == 1] <- 0
   #      warning(mstyle$warning("Inner factor has k=1 for one or more levels. Corresponding 'tau2' value(s) fixed to 0."), call.=FALSE)
   #   }
   #}

   ### check if each study has only a single arm (could be different arms!)
   ### for "CS","HCS","AR","HAR","CAR" must then fix rho to 0 (if not already fixed)
   ### for SP/PHY structures cannot fit model; for GEN rho may still be (weakly) identifiable

   if (g.nlevels[2] == nrow(mf.g)) {
      if (is.element(struct, c("CS","HCS","AR","HAR","CAR")) && is.na(rho)) {
         rho <- 0
         warning(mstyle$warning(paste0("Each level of the outer factor contains only a single level of the inner factor, so fixed value of ", ifelse(isG, 'rho', 'phi'), " to 0.")), call.=FALSE)
      }
      if (is.element(struct, c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD")))
         stop(mstyle$stop("Cannot fit model since each level of the outer factor contains only a single level of the inner term."), call.=FALSE)
   }

   g.levels.comb.k <- NULL

   if (!is.element(struct, c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD","GEN","GDIAG"))) {

      ### create matrix where each row (= study) indicates how often each arm occurred
      ### then turn this into a list (with each element equal to a row (= study))

      g.levels.comb.k <- crossprod(Z.G2, Z.G1)
      g.levels.comb.k <- split(g.levels.comb.k, seq_len(nrow(g.levels.comb.k)))

      ### create matrix for each element (= study) that indicates which combinations occurred
      ### sum up all matrices (numbers indicate in how many studies each combination occurred)
      ### take upper triangle part that corresponds to the arm combinations (in order of rho)

      g.levels.comb.k <- lapply(g.levels.comb.k, function(x) outer(x,x, FUN="&"))
      g.levels.comb.k <- Reduce("+", g.levels.comb.k)
      g.levels.comb.k <- g.levels.comb.k[lower.tri(g.levels.comb.k)]

      ### UN/UNR: if a particular combination of arms never occurs in any of the studies, then must fix the corresponding rho to 0 (if not already fixed)
      ### this also takes care of the case where each study has only a single arm

      if (is.element(struct, c("UN","UNR")) && any(g.levels.comb.k == 0 & is.na(rho))) {
         rho[g.levels.comb.k == 0] <- 0
         warning(mstyle$warning(paste0("Some combinations of the levels of the inner factor never occurred. Corresponding ", ifelse(isG, 'rho', 'phi'), " value(s) fixed to 0.")), call.=FALSE)
      }

      ### if there was only a single arm for "UN" or "UNR" to begin with, then fix rho to 0
      ### (technically there is then no rho at all to begin with, but rhos was still set to 1 earlier for the optimization routine)
      ### (if there is a single arm after removing NAs, then this is dealt with below by setting tau2 and rho values to 0)

      if (is.element(struct, c("UN","UNR")) && g.nlevels.f[1] == 1 && is.na(rho)) {
         rho <- 0
         warning(mstyle$warning(paste0("Inner factor has only a single level, so fixed value of ", ifelse(isG, 'rho', 'phi'), " to 0.")), call.=FALSE)
      }

   }

   ### construct G matrix for the various structures

   if (struct == "CS") {
      G <- matrix(rho*tau2, nrow=g.nlevels.f[1], ncol=g.nlevels.f[1])
      diag(G) <- tau2
   }

   if (struct == "HCS") {
      G <- matrix(rho, nrow=g.nlevels.f[1], ncol=g.nlevels.f[1])
      diag(G) <- 1
      G <- diag(sqrt(tau2), nrow=g.nlevels.f[1], ncol=g.nlevels.f[1]) %*% G %*% diag(sqrt(tau2), nrow=g.nlevels.f[1], ncol=g.nlevels.f[1])
      diag(G) <- tau2
   }

   if (is.element(struct, c("UN","GEN"))) {
      G <- .con.vcov.UN(tau2, rho)
   }

   if (struct == "UNR") {
      G <- .con.vcov.UNR(tau2, rho)
   }

   if (is.element(struct, c("GDIAG"))) {
      G <- diag(tau2, nrow=length(tau2), ncol=length(tau2))
   }

   if (is.element(struct, c("ID","DIAG"))) {
      G <- diag(tau2, nrow=g.nlevels.f[1], ncol=g.nlevels.f[1])
   }

   if (struct == "AR") {
      if (is.na(rho)) {
         G <- matrix(NA_real_, nrow=g.nlevels.f[1], ncol=g.nlevels.f[1])
      } else {
         ### is g.nlevels.f[1] == 1 even possible here?
         if (g.nlevels.f[1] > 1) {
            G <- toeplitz(ARMAacf(ar=rho, lag.max=g.nlevels.f[1]-1))
         } else {
            G <- diag(1)
         }
      }
      G <- diag(sqrt(rep(tau2, g.nlevels.f[1])), nrow=g.nlevels.f[1], ncol=g.nlevels.f[1]) %*% G %*% diag(sqrt(rep(tau2, g.nlevels.f[1])), nrow=g.nlevels.f[1], ncol=g.nlevels.f[1])
      diag(G) <- tau2
   }

   if (struct == "HAR") {
      if (is.na(rho)) {
         G <- matrix(NA_real_, nrow=g.nlevels.f[1], ncol=g.nlevels.f[1])
      } else {
         ### is g.nlevels.f[1] == 1 even possible here?
         if (g.nlevels.f[1] > 1) {
            G <- toeplitz(ARMAacf(ar=rho, lag.max=g.nlevels.f[1]-1))
         } else {
            G <- diag(1)
         }
      }
      G <- diag(sqrt(tau2), nrow=g.nlevels.f[1], ncol=g.nlevels.f[1]) %*% G %*% diag(sqrt(tau2), nrow=g.nlevels.f[1], ncol=g.nlevels.f[1])
      diag(G) <- tau2
   }

   if (struct == "CAR") {
      if (is.na(rho)) {
         G <- matrix(NA_real_, nrow=g.nlevels.f[1], ncol=g.nlevels.f[1])
      } else {
         ### is g.nlevels.f[1] == 1 even possible here?
         if (g.nlevels.f[1] > 1) {
            G <- outer(g.values, g.values, function(x,y) rho^(abs(x-y)))
         } else {
            G <- diag(1)
         }
      }
      G <- diag(sqrt(rep(tau2, g.nlevels.f[1])), nrow=g.nlevels.f[1], ncol=g.nlevels.f[1]) %*% G %*% diag(sqrt(rep(tau2, g.nlevels.f[1])), nrow=g.nlevels.f[1], ncol=g.nlevels.f[1])
      diag(G) <- tau2
   }

   if (is.element(struct, c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD"))) {

      ### remove the '| outer' part from the formula and add '- 1'
      formula <- as.formula(paste0(strsplit(paste0(formula, collapse=""), "|", fixed=TRUE)[[1]][1], "- 1", collapse=""))

      ### create distance matrix

      if (.is.matrix(distspec)) {

         if (anyNA(distspec))
            stop(mstyle$stop("No missing values allowed in matrices specified via 'dist'."), call.=FALSE)
         if (!.is.square(distspec))
            stop(mstyle$stop("Distance matrices specified via 'dist' must be square matrices."), call.=FALSE)
         if (!isSymmetric(unname(distspec)))
            stop(mstyle$stop("Distance matrices specified via 'dist' must be symmetric matrices."), call.=FALSE)
         if (is.null(rownames(distspec)))
            rownames(distspec) <- colnames(distspec)
         if (is.null(colnames(distspec)))
            colnames(distspec) <- rownames(distspec)
         if (length(colnames(distspec)) != length(unique(colnames(distspec))))
            stop(mstyle$stop("Distance matrices specified via 'dist' must have unique dimension names."), call.=FALSE)
         if (any(!is.element(as.character(mf.g[[1]]), colnames(distspec))))
            stop(mstyle$stop(paste0("There are levels in '", colnames(mf.g)[1], "' for which there are no matching rows/columns in the corresponding 'dist' matrix.")), call.=FALSE)

         if (is.element(struct, c("PHYBM","PHYPL","PHYPD")) && !all.equal(min(distspec), 0))
            warning(mstyle$warning("Minimum value in the distance matrix is not 0."), call.=FALSE)

         if (is.element(struct, c("PHYBM","PHYPL","PHYPD")) && !all.equal(max(distspec), 2))
            warning(mstyle$warning("Maximum value in the distance matrix is not 2."), call.=FALSE)

         Dmat <- distspec[as.character(mf.g[[1]]), as.character(mf.g[[1]])]

      } else {

         if (is.element(struct, c("PHYBM","PHYPL","PHYPD")))
            stop(mstyle$stop("Must supply distance matrix via 'dist' for phylogenetic correlation structures."), call.=FALSE)

         Cmat <- model.matrix(formula, data=mf.g[-nvars])

         if (is.function(distspec)) {
            Dmat <- distspec(Cmat)
         } else {
            if (is.element(distspec, c("euclidean", "maximum", "manhattan")))
               Dmat <- as.matrix(dist(Cmat, method=distspec))
            if (distspec == "gcd")
               Dmat <- sp::spDists(Cmat, longlat=TRUE)
         }

      }

      if (sparse)
         Dmat <- Matrix(Dmat, sparse=TRUE)

   } else {

      Dmat <- NULL

   }

   if (struct == "SPEXP") {
      Rmat <- exp(-Dmat/rho)
      G <- tau2 * Rmat * tcrossprod(Z.G2)
   }

   if (struct == "SPGAU") {
      Rmat <- exp(-Dmat^2/rho^2)
      G <- tau2 * Rmat * tcrossprod(Z.G2)
   }

   if (struct == "SPLIN") {
      Rmat <- (1 - Dmat/rho) * I(Dmat < rho)
      G <- tau2 * Rmat * tcrossprod(Z.G2)
   }

   if (struct == "SPRAT") {
      Rmat <- 1 - (Dmat/rho)^2 / (1 + (Dmat/rho)^2)
      G <- tau2 * Rmat * tcrossprod(Z.G2)
   }

   if (struct == "SPSPH") {
      Rmat <- (1 - 3/2*Dmat/rho + 1/2*(Dmat/rho)^3) * I(Dmat < rho)
      G <- tau2 * Rmat * tcrossprod(Z.G2)
   }

   if (struct == "PHYBM") {
      rho <- max(Dmat)
      Rmat <- 1 - Dmat/rho
      G <- tau2 * Rmat * tcrossprod(Z.G2)
   }

   if (struct == "PHYPL") {
      Rmat <- rho * (1 - Dmat/max(Dmat))
      diag(Rmat) <- 1
      Rmat[Dmat == 0] <- 1
      G <- tau2 * Rmat * tcrossprod(Z.G2)
   }

   if (struct == "PHYPD") {
      Rmat <- 1 - Dmat/max(Dmat)
      G <- tau2 * Rmat^rho * tcrossprod(Z.G2)
   }

   ### for spatial and phylogeny structures, compute a much more sensible initial value for rho

   if (is.element(struct, c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD"))) {
      if (struct == "PHYBM")
         rho.init <- max(Dmat)
      if (struct == "PHYPL")
         rho.init <- 0.5
      if (struct == "PHYPD")
         rho.init <- 1
      if (!is.element(struct, c("PHYBM","PHYPL","PHYPD")))
         rho.init <- unname(suppressMessages(quantile(Dmat[lower.tri(Dmat)], 0.25))) # suppressMessages() to avoid '<sparse>[ <logic> ] : .M.sub.i.logical() maybe inefficient' messages when sparse=TRUE
   } else {
      rho.init <- NULL
   }

   ### for "CS","AR","CAR","ID" set tau2 value to 0 for any levels that were removed

   if (any(g.levels.r) && is.element(struct, c("CS","AR","CAR","ID"))) {
      G[g.levels.r,] <- 0
      G[,g.levels.r] <- 0
   }

   ### for "HCS","HAR","DIAG" set tau2 value(s) to 0 for any levels that were removed

   if (any(g.levels.r) && is.element(struct, c("HCS","HAR","DIAG"))) {
      G[g.levels.r,] <- 0
      G[,g.levels.r] <- 0
      tau2[g.levels.r] <- 0
      warning(mstyle$warning(paste0("Fixed ", ifelse(isG, 'tau2', 'gamma2'), " to 0 for removed level(s).")), call.=FALSE)
   }

   ### for "UN", set tau2 value(s) and corresponding rho(s) to 0 for any levels that were removed

   if (any(g.levels.r) && struct == "UN") {
      G[g.levels.r,] <- 0
      G[,g.levels.r] <- 0
      tau2[g.levels.r] <- 0
      rho <- G[lower.tri(G)]
      warning(mstyle$warning(paste0("Fixed ", ifelse(isG, 'tau2', 'gamma2'), " and corresponding ", ifelse(isG, 'rho', 'phi'), " value(s) to 0 for removed level(s).")), call.=FALSE)
   }

   ### for "UNR", set rho(s) to 0 corresponding to any levels that were removed

   if (any(g.levels.r) && struct == "UNR") {
      G[g.levels.r,] <- 0
      G[,g.levels.r] <- 0
      diag(G) <- tau2 ### don't really need this
      rho <- G[lower.tri(G)]
      warning(mstyle$warning(paste0("Fixed ", ifelse(isG, 'rho', 'phi'), " value(s) to 0 for removed level(s).")), call.=FALSE)
   }

   ### special handling for the bivariate model:
   ### if tau2 (for "CS","AR","CAR","UNR") or either tau2.1 or tau2.2 (for "HCS","UN","HAR") is fixed to 0, then rho must be fixed to 0

   if (g.nlevels.f[1] == 2) {
      if (is.element(struct, c("CS","AR","CAR","UNR")) && !is.na(tau2) && tau2 == 0)
         rho <- 0
      if (is.element(struct, c("HCS","UN","HAR")) && ((!is.na(tau2[1]) && tau2[1] == 0) || (!is.na(tau2[2]) && tau2[2] == 0)))
         rho <- 0
   }

   return(list(mf.g=mf.g, g.nlevels=g.nlevels, g.nlevels.f=g.nlevels.f,
               g.levels=g.levels, g.levels.f=g.levels.f, g.levels.r=g.levels.r, g.levels.k=g.levels.k, g.levels.comb.k=g.levels.comb.k,
               tau2=tau2, rho=rho, G=G, Dmat=Dmat, rho.init=rho.init))

}

############################################################################

### function to construct var-cov matrix for "UN" and "GEN" structures given vector of variances and correlations

.con.vcov.UN <- function(vars, cors, vccov=FALSE) {
   dims <- length(vars)
   if (vccov) {
      G <- matrix(0, nrow=dims, ncol=dims)
      G[lower.tri(G)] <- cors
      G[upper.tri(G)] <- t(G)[upper.tri(G)]
      diag(G) <- vars
      return(G)
   } else {
      R <- matrix(1, nrow=dims, ncol=dims)
      R[lower.tri(R)] <- cors
      R[upper.tri(R)] <- t(R)[upper.tri(R)]
      S <- diag(sqrt(vars), nrow=dims, ncol=dims)
      return(S %*% R %*% S)
   }
}

### function to construct var-cov matrix for "UN" and "GEN" structures given vector of 'choled' variances and covariances

.con.vcov.UN.chol <- function(vars, covs) {
   dims <- length(vars)
   G <- matrix(0, nrow=dims, ncol=dims)
   G[lower.tri(G)] <- covs
   diag(G) <- vars
   return(crossprod(G))
}

### function to construct var-cov matrix for "UNR" structure given the variance and correlations

.con.vcov.UNR <- function(var, cors) {
   dims <- round((1 + sqrt(1 + 8*length(cors)))/2)
   G <- matrix(1, nrow=dims, ncol=dims)
   G[lower.tri(G)] <- cors
   G[upper.tri(G)] <- t(G)[upper.tri(G)]
   return(var * G)
}

### function to construct var-cov matrix for "UNR" structure given the variance and vector of 'choled' correlations

.con.vcov.UNR.chol <- function(var, cors) {
   dims <- round((1 + sqrt(1 + 8*length(cors)))/2)
   G <- matrix(0, nrow=dims, ncol=dims)
   G[lower.tri(G)] <- cors
   diag(G) <- 1
   return(var * crossprod(G))
}

############################################################################

### function to construct var-cov matrix (G or H) for '~ inner | outer' terms

.con.E <- function(v, r, v.val, r.val, Z1, Z2, levels.r, values, Dmat, struct, cholesky, vctransf, vccov, nearpd, sparse) {

   ### if cholesky=TRUE, back-transformation/substitution is done below; otherwise, back-transform and replace fixed values

   if (!cholesky) {
      if (vctransf) {
         v <- ifelse(is.na(v.val), exp(v), v.val)           ### variances are optimized in log space, so exponentiate
         if (struct == "CAR")
            r <- ifelse(is.na(r.val), plogis(r), r.val)     ### CAR correlation is optimized in qlogis space, so use plogis
         if (is.element(struct, c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD")))
            r <- ifelse(is.na(r.val), exp(r), r.val)        ### spatial and phylogenetic 'correlation' parameter is optimized in log space, so exponentiate
         if (!is.element(struct, c("CAR","SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD")))
            r <- ifelse(is.na(r.val), tanh(r), r.val)       ### other correlations are optimized in atanh space, so use tanh
      } else {
         ### for Hessian computation, can choose to leave as is
         v <- ifelse(is.na(v.val), v, v.val)
         r <- ifelse(is.na(r.val), r, r.val)
         v[v < 0] <- 0
         if (struct == "CAR") {
            r[r < 0] <- 0
            r[r > 1] <- 1
         }
         if (is.element(struct, c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD"))) {
            r[r < 0] <- 0
         }
         if (!is.element(struct, c("CAR","SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD")) && !vccov) {
            r[r < -1] <- -1
            r[r > 1] <- 1
         }
      }
      v <- ifelse(v <= .Machine$double.eps*10, 0, v) ### don't do this with Cholesky factorization, since values can be negative
   }

   ncol.Z1 <- ncol(Z1)

   if (struct == "CS") {
      E <- matrix(r*v, nrow=ncol.Z1, ncol=ncol.Z1)
      diag(E) <- v
   }

   if (struct == "HCS") {
      E <- matrix(r, nrow=ncol.Z1, ncol=ncol.Z1)
      diag(E) <- 1
      E <- diag(sqrt(v), nrow=ncol.Z1, ncol=ncol.Z1) %*% E %*% diag(sqrt(v), nrow=ncol.Z1, ncol=ncol.Z1)
      diag(E) <- v
   }

   if (is.element(struct, c("UN","GEN"))) {
      if (cholesky) {
         E <- .con.vcov.UN.chol(v, r)
         v <- diag(E)                  ### need this, so correct values are shown when verbose=TRUE
         r <- cov2cor(E)[lower.tri(E)] ### need this, so correct values are shown when verbose=TRUE
         v[!is.na(v.val)] <- v.val[!is.na(v.val)] ### replace any fixed values
         r[!is.na(r.val)] <- r.val[!is.na(r.val)] ### replace any fixed values
      }
      E <- .con.vcov.UN(v, r, vccov)
      if (nearpd) {
         E <- as.matrix(nearPD(E)$mat) ### nearPD() in Matrix package
         v <- diag(E)                  ### need this, so correct values are shown when verbose=TRUE
         r <- cov2cor(E)[lower.tri(E)] ### need this, so correct values are shown when verbose=TRUE
      }
   }

   if (struct == "UNR") {
      if (cholesky) {
         E <- .con.vcov.UNR.chol(v, r)
         v <- diag(E)[1,1]             ### need this, so correct values are shown when verbose=TRUE
         r <- cov2cor(E)[lower.tri(E)] ### need this, so correct values are shown when verbose=TRUE
         v[!is.na(v.val)] <- v.val[!is.na(v.val)] ### replace any fixed values
         r[!is.na(r.val)] <- r.val[!is.na(r.val)] ### replace any fixed values
      }
      E <- .con.vcov.UNR(v, r)
      if (nearpd) {
         E <- as.matrix(nearPD(E, keepDiag=TRUE)$mat) ### nearPD() in Matrix package
         v <- E[1,1]                   ### need this, so correct values are shown when verbose=TRUE
         r <- cov2cor(E)[lower.tri(E)] ### need this, so correct values are shown when verbose=TRUE
      }
   }

   if (struct == "GDIAG") {
      E <- diag(v, nrow=length(v), ncol=length(v))
   }

   if (is.element(struct, c("ID","DIAG")))
      E <- diag(v, nrow=ncol.Z1, ncol=ncol.Z1)

   if (struct == "AR") {
      if (ncol.Z1 > 1) {
         E <- toeplitz(ARMAacf(ar=r, lag.max=ncol.Z1-1))
      } else {
         E <- diag(1)
      }
      E <- diag(sqrt(rep(v, ncol.Z1)), nrow=ncol.Z1, ncol=ncol.Z1) %*% E %*% diag(sqrt(rep(v, ncol.Z1)), nrow=ncol.Z1, ncol=ncol.Z1)
      diag(E) <- v
   }

   if (struct == "HAR") {
      if (ncol.Z1 > 1) {
         E <- toeplitz(ARMAacf(ar=r, lag.max=ncol.Z1-1))
      } else {
         E <- diag(1)
      }
      E <- diag(sqrt(v), nrow=ncol.Z1, ncol=ncol.Z1) %*% E %*% diag(sqrt(v), nrow=ncol.Z1, ncol=ncol.Z1)
      diag(E) <- v
   }

   if (struct == "CAR") {
      if (ncol.Z1 > 1) {
         E <- outer(values, values, function(x,y) r^(abs(x-y)))
      } else {
         E <- diag(1)
      }
      E <- diag(sqrt(rep(v, ncol.Z1)), nrow=ncol.Z1, ncol=ncol.Z1) %*% E %*% diag(sqrt(rep(v, ncol.Z1)), nrow=ncol.Z1, ncol=ncol.Z1)
      diag(E) <- v
   }

   if (struct == "SPEXP")
      E <- v * exp(-Dmat/r) * tcrossprod(Z2)

   if (struct == "SPGAU")
      E <- v * exp(-Dmat^2/r^2) * tcrossprod(Z2)

   if (struct == "SPLIN")
      E <- v * ((1 - Dmat/r) * I(Dmat < r)) * tcrossprod(Z2)

   if (struct == "SPRAT")
      E <- v * (1 - (Dmat/r)^2 / (1 + (Dmat/r)^2)) * tcrossprod(Z2)

   if (struct == "SPSPH")
      E <- v * ((1 - 3/2*Dmat/r + 1/2*(Dmat/r)^3) * I(Dmat < r)) * tcrossprod(Z2)

   if (struct == "PHYBM") {
      r <- max(Dmat)
      E <- 1 - Dmat/r
      E <- v * E * tcrossprod(Z2)
   }

   if (struct == "PHYPL") {
      E <- r * (1 - Dmat/max(Dmat))
      diag(E) <- 1
      E[Dmat == 0] <- 1
      E <- v * E * tcrossprod(Z2)
   }

   if (struct == "PHYPD") {
      E <- 1 - Dmat/max(Dmat)
      E <- v * E^r * tcrossprod(Z2)
   }

   ### set variance and corresponding correlation value(s) to 0 for any levels that were removed

   if (!is.element(struct, c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","PHYBM","PHYPL","PHYPD","GEN","GDIAG")) && any(levels.r)) {
      E[levels.r,] <- 0
      E[,levels.r] <- 0
   }

   if (sparse)
      E <- Matrix(E, sparse=TRUE)

   return(list(v=v, r=r, E=E))

}

############################################################################

### -1 times the log likelihood (regular or restricted) for rma.mv models

.ll.rma.mv <- function(par, reml, Y, M, A, X, k, pX, # note: pX due to nlm(); M=V to begin with
                       D.S, Z.G1, Z.G2, Z.H1, Z.H2, g.Dmat, h.Dmat,
                       sigma2.val, tau2.val, rho.val, gamma2.val, phi.val,
                       sigma2s, tau2s, rhos, gamma2s, phis,
                       withS, withG, withH,
                       struct, g.levels.r, h.levels.r, g.values, h.values,
                       sparse, cholesky, nearpd, vctransf, vccov,
                       verbose, digits, REMLf, dofit=FALSE, hessian=FALSE) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   ### only NA values in sigma2.val, tau2.val, rho.val, gamma2.val, phi.val should be estimated; otherwise, replace with fixed values

   if (withS) {

      if (vctransf) {
         sigma2 <- ifelse(is.na(sigma2.val), exp(par[seq_len(sigma2s)]), sigma2.val) ### sigma2 is optimized in log space, so exponentiate
      } else {
         sigma2 <- ifelse(is.na(sigma2.val), par[seq_len(sigma2s)], sigma2.val)      ### for Hessian computation, can choose to leave as is
         sigma2[sigma2 < 0] <- 0
      }

      #if (any(is.nan(sigma2)))
      #   return(Inf)

      ### set really small sigma2 values equal to 0 (anything below .Machine$double.eps*10 is essentially 0)
      sigma2 <- ifelse(sigma2 <= .Machine$double.eps*10, 0, sigma2)

      for (j in seq_len(sigma2s)) {
         M <- M + sigma2[j] * D.S[[j]]
      }

   }

   if (withG) {

      vars <- par[(sigma2s+1):(sigma2s+tau2s)]
      cors <- par[(sigma2s+tau2s+1):(sigma2s+tau2s+rhos)]

      resG <- .con.E(v=vars, r=cors,
                     v.val=tau2.val, r.val=rho.val, Z1=Z.G1, Z2=Z.G2, levels.r=g.levels.r, values=g.values, Dmat=g.Dmat,
                     struct=struct[1], cholesky=cholesky[1], vctransf=vctransf, vccov=vccov, nearpd=nearpd, sparse=sparse)
      tau2 <- resG$v
      rho  <- resG$r
      G    <- resG$E

      M <- M + (Z.G1 %*% G %*% t(Z.G1)) * tcrossprod(Z.G2)

   }

   if (withH) {

      vars <- par[(sigma2s+tau2s+rhos+1):(sigma2s+tau2s+rhos+gamma2s)]
      cors <- par[(sigma2s+tau2s+rhos+gamma2s+1):(sigma2s+tau2s+rhos+gamma2s+phis)]

      resH <- .con.E(v=vars, r=cors,
                     v.val=gamma2.val, r.val=phi.val, Z1=Z.H1, Z2=Z.H2, levels.r=h.levels.r, values=h.values, Dmat=h.Dmat,
                     struct=struct[2], cholesky=cholesky[2], vctransf=vctransf, vccov=vccov, nearpd=nearpd, sparse=sparse)
      gamma2 <- resH$v
      phi    <- resH$r
      H      <- resH$E

      M <- M + (Z.H1 %*% H %*% t(Z.H1)) * tcrossprod(Z.H2)

   }

   ### put estimates so far into .metafor environment

   if (!hessian) {
      pars <- list(sigma2 = if (withS) sigma2 else NULL,
                   tau2   = if (withG) tau2 else NULL,
                   rho    = if (withG) rho else NULL,
                   gamma2 = if (withH) gamma2 else NULL,
                   phi    = if (withH) phi else NULL)
      try(assign("rma.mv", pars, envir=.metafor), silent=TRUE)
   }

   ### note: if M is sparse, then using nearPD() could blow up

   if (nearpd)
      M <- as.matrix(nearPD(M)$mat)

   if (verbose > 1) {
      W <- try(chol2inv(chol(M)), silent=FALSE)
   } else {
      W <- try(suppressWarnings(chol2inv(chol(M))), silent=TRUE)
   }

   ### note: need W for REML llval computation

   if (inherits(W, "try-error")) {

      ### if M is not positive-definite, set the (restricted) log likelihood to -Inf
      ### this idea is based on: https://stats.stackexchange.com/q/11368/1934 (this is crude, but should
      ### move the parameter estimates away from values that create the non-positive-definite M matrix)

      if (dofit) {
         stop(mstyle$stop("Final variance-covariance matrix not positive definite."), call.=FALSE)
      } else {
         llval <- -Inf
      }

   } else {

      if (verbose > 1) {
         U <- try(chol(W), silent=FALSE)
      } else {
         U <- try(suppressWarnings(chol(W)), silent=TRUE)
      }

      ### Y ~ N(Xbeta, M), so UY ~ N(UXbeta, UMU) where UMU = I
      ### return(U %*% M %*% U)

      if (inherits(U, "try-error")) {

         if (dofit) {
            stop(mstyle$stop("Cannot fit model based on estimated marginal variance-covariance matrix."), call.=FALSE)
         } else {
            llval <- -Inf
         }

      } else {

         if (!dofit || is.null(A)) {

            sX   <- U %*% X
            sY   <- U %*% Y
            beta <- solve(crossprod(sX), crossprod(sX, sY))
            RSS  <- sum(as.vector(sY - sX %*% beta)^2)
            if (dofit)
               vb <- matrix(solve(crossprod(sX)), nrow=pX, ncol=pX)

         } else {

            stXAX <- chol2inv(chol(as.matrix(t(X) %*% A %*% X)))
            #stXAX <- tcrossprod(qr.solve(sX, diag(k)))
            beta  <- matrix(stXAX %*% crossprod(X,A) %*% Y, ncol=1)
            RSS   <- as.vector(t(Y - X %*% beta) %*% W %*% (Y - X %*% beta))
            vb    <- matrix(stXAX %*% t(X) %*% A %*% M %*% A %*% X %*% stXAX, nrow=pX, ncol=pX)

         }

         llvals <- c(NA_real_, NA_real_)

         if (dofit || !reml)
            llvals[1]  <- -1/2 * (k) * log(2*base::pi) - 1/2 * determinant(M, logarithm=TRUE)$modulus - 1/2 * RSS

         if (dofit || reml)
            llvals[2]  <- -1/2 * (k-pX) * log(2*base::pi) + ifelse(REMLf, 1/2 * determinant(crossprod(X), logarithm=TRUE)$modulus, 0) +
                          -1/2 * determinant(M, logarithm=TRUE)$modulus - 1/2 * determinant(crossprod(X,W) %*% X, logarithm=TRUE)$modulus - 1/2 * RSS

         if (dofit) {

            res <- list(beta=beta, vb=vb, M=M, llvals=llvals)

            if (withS)
               res$sigma2 <- sigma2

            if (withG) {
               res$G <- G
               res$tau2 <- tau2
               res$rho <- rho
            }

            if (withH) {
               res$H <- H
               res$gamma2 <- gamma2
               res$phi <- phi
            }

            return(res)

         } else {

            llval <- ifelse(reml, llvals[2], llvals[1])

         }

      }

   }

   if ((vctransf && verbose) || (!vctransf && (verbose > 1))) {

      if (!hessian) {

         iteration <- .getfromenv("iteration", default=NULL)

         if (!is.null(iteration)) {
            #cat(mstyle$verbose(paste0("Iteration ", iteration, "\t")))
            cat(mstyle$verbose(paste0("Iteration ", formatC(iteration, width=5, flag="-", format="f", digits=0), " ")))
            try(assign("iteration", iteration+1, envir=.metafor), silent=TRUE)
         }

      }

      cat(mstyle$verbose(paste0("ll = ",             fmtx(llval,  digits[["fit"]], flag=" "))), "  ")
      if (withS)
         cat(mstyle$verbose(paste0("sigma2 =", paste(fmtx(sigma2, digits[["var"]], flag=" "), collapse=" "), "  ")))
      if (withG) {
         cat(mstyle$verbose(paste0("tau2 =",   paste(fmtx(tau2,   digits[["var"]], flag=" "), collapse=" "), "  ")))
         cat(mstyle$verbose(paste0("rho =",    paste(fmtx(rho,    digits[["var"]], flag=" "), collapse=" "), "  ")))
      }
      if (withH) {
         cat(mstyle$verbose(paste0("gamma2 =", paste(fmtx(gamma2, digits[["var"]], flag=" "), collapse=" "), "  ")))
         cat(mstyle$verbose(paste0("phi =",    paste(fmtx(phi,    digits[["var"]], flag=" "), collapse=" "), "  ")))
      }
      cat("\n")

   }

   return(-1 * c(llval))

}

############################################################################

.cooks.distance.rma.mv <- function(i, obj, parallel, svb, cluster, ids, reestimate, btt) {

   if (parallel == "snow")
      library(metafor)

   incl <- cluster %in% ids[i]

   ### elements that need to be returned

   outlist <- "coef.na=coef.na, beta=beta"

   ### note: not.na=FALSE only when there are missings in data, not when model below cannot be fitted or results in dropped coefficients

   if (reestimate) {

      ### set initial values to estimates from full model

      control             <- obj$control
      control$sigma2.init <- obj$sigma2
      control$tau2.init   <- obj$tau2
      control$rho.init    <- obj$rho
      control$gamma2.init <- obj$gamma2
      control$phi.init    <- obj$phi

      ### fit model without data from ith cluster

      args <- list(yi=obj$yi, V=obj$V, W=obj$W, mods=obj$X, random=obj$random, struct=obj$struct, intercept=FALSE, data=obj$mf.r, method=obj$method, test=obj$test, dfs=obj$dfs, level=obj$level, R=obj$R, Rscale=obj$Rscale, sigma2=ifelse(obj$vc.fix$sigma2, obj$sigma2, NA), tau2=ifelse(obj$vc.fix$tau2, obj$tau2, NA), rho=ifelse(obj$vc.fix$rho, obj$rho, NA), gamma2=ifelse(obj$vc.fix$gamma2, obj$gamma2, NA), phi=ifelse(obj$vc.fix$phi, obj$phi, NA), sparse=obj$sparse, dist=obj$dist, control=control, subset=!incl, outlist=outlist)
      res <- try(suppressWarnings(.do.call(rma.mv, args)), silent=TRUE)

   } else {

      ### set values of variance/correlation components to those from the 'full' model

      args <- list(yi=obj$yi, V=obj$V, W=obj$W, mods=obj$X, random=obj$random, struct=obj$struct, intercept=FALSE, data=obj$mf.r, method=obj$method, test=obj$test, dfs=obj$dfs, level=obj$level, R=obj$R, Rscale=obj$Rscale, sigma2=obj$sigma2, tau2=obj$tau2, rho=obj$rho, gamma2=obj$gamma2, phi=obj$phi, sparse=obj$sparse, dist=obj$dist, control=obj$control, subset=!incl, outlist=outlist)
      res <- try(suppressWarnings(.do.call(rma.mv, args)), silent=TRUE)

   }

   if (inherits(res, "try-error"))
      return(list(cook.d = NA))

   ### removing a cluster could lead to a model coefficient becoming inestimable

   if (any(res$coef.na))
      return(list(cook.d = NA))

   ### compute dfbeta value(s) (including coefficients as specified via btt)

   dfb <- obj$beta[btt] - res$beta[btt]

   ### compute Cook's distance

   return(list(cook.d = crossprod(dfb,svb) %*% dfb))

}

.rstudent.rma.mv <- function(i, obj, parallel, cluster, ids, reestimate) {

   if (parallel == "snow")
      library(metafor)

   incl <- cluster %in% ids[i]

   k.id <- sum(incl)

   ### elements that need to be returned

   outlist <- "coef.na=coef.na, sigma2=sigma2, tau2=tau2, rho=rho, gamma2=gamma2, phi=phi, beta=beta, vb=vb"

   if (reestimate) {

      ### set initial values to estimates from full model

      control             <- obj$control
      control$sigma2.init <- obj$sigma2
      control$tau2.init   <- obj$tau2
      control$rho.init    <- obj$rho
      control$gamma2.init <- obj$gamma2
      control$phi.init    <- obj$phi

      ### fit model without data from ith cluster

      args <- list(yi=obj$yi, V=obj$V, W=obj$W, mods=obj$X, random=obj$random, struct=obj$struct, intercept=FALSE, data=obj$mf.r, method=obj$method, test=obj$test, dfs=obj$dfs, level=obj$level, R=obj$R, Rscale=obj$Rscale, sigma2=ifelse(obj$vc.fix$sigma2, obj$sigma2, NA), tau2=ifelse(obj$vc.fix$tau2, obj$tau2, NA), rho=ifelse(obj$vc.fix$rho, obj$rho, NA), gamma2=ifelse(obj$vc.fix$gamma2, obj$gamma2, NA), phi=ifelse(obj$vc.fix$phi, obj$phi, NA), sparse=obj$sparse, dist=obj$dist, control=control, subset=!incl, outlist=outlist)
      res <- try(suppressWarnings(.do.call(rma.mv, args)), silent=TRUE)

   } else {

      ### set values of variance/correlation components to those from the 'full' model

      args <- list(yi=obj$yi, V=obj$V, W=obj$W, mods=obj$X, random=obj$random, struct=obj$struct, intercept=FALSE, data=obj$mf.r, method=obj$method, test=obj$test, dfs=obj$dfs, level=obj$level, R=obj$R, Rscale=obj$Rscale, sigma2=obj$sigma2, tau2=obj$tau2, rho=obj$rho, gamma2=obj$gamma2, phi=obj$phi, sparse=obj$sparse, dist=obj$dist, control=obj$control, subset=!incl, outlist=outlist)
      res <- try(suppressWarnings(.do.call(rma.mv, args)), silent=TRUE)

   }

   if (inherits(res, "try-error"))
      return(list(delresid = rep(NA, k.id), sedelresid = rep(NA, k.id), X2 = NA, k.id = NA, pos = which(incl)))

   ### removing a cluster could lead to a model coefficient becoming inestimable

   if (any(res$coef.na))
      return(list(delresid = rep(NA, k.id), sedelresid = rep(NA, k.id), X2 = NA, k.id = NA, pos = which(incl)))

   ### elements that need to be returned

   outlist <- "M=M"

   ### fit model based on all data but with var/cor components fixed to those from res

   args <- list(yi=obj$yi, V=obj$V, W=obj$W, mods=obj$X, random=obj$random, struct=obj$struct, intercept=FALSE, data=obj$mf.r, method=obj$method, test=obj$test, dfs=obj$dfs, level=obj$level, R=obj$R, Rscale=obj$Rscale, sigma2=res$sigma2, tau2=res$tau2, rho=res$rho, gamma2=res$gamma2, phi=res$phi, sparse=obj$sparse, dist=obj$dist, control=obj$control, outlist=outlist)
   tmp <- try(suppressWarnings(.do.call(rma.mv, args)), silent=TRUE)

   Xi <- obj$X[incl,,drop=FALSE]
   delpred  <- Xi %*% res$beta
   vdelpred <- Xi %*% res$vb %*% t(Xi)
   delresid <- c(obj$yi[incl] - delpred)
   sedelresid <- c(sqrt(diag(tmp$M[incl,incl,drop=FALSE] + vdelpred)))

   sve <- try(chol2inv(chol(tmp$M[incl,incl,drop=FALSE] + vdelpred)), silent=TRUE)
   #sve <- try(solve(tmp$M[incl,incl,drop=FALSE] + vdelpred), silent=TRUE)

   if (inherits(sve, "try-error"))
      return(list(delresid = delresid, sedelresid = sedelresid, X2 = NA, k.id = k.id, pos = which(incl)))

   X2 <- c(rbind(delresid) %*% sve %*% cbind(delresid))

   if (is.list(X2)) # when sparse=TRUE, this is a list with a one-element matrix
      X2 <- X2[[1]][1]

   return(list(delresid = delresid, sedelresid = sedelresid, X2 = X2, k.id = k.id, pos = which(incl)))

}

.dfbetas.rma.mv <- function(i, obj, parallel, cluster, ids, reestimate) {

   if (parallel == "snow")
      library(metafor)

   incl <- cluster %in% ids[i]

   ### elements that need to be returned

   outlist <- "coef.na=coef.na, sigma2=sigma2, tau2=tau2, rho=rho, gamma2=gamma2, phi=phi, beta=beta"

   if (reestimate) {

      ### set initial values to estimates from full model

      control             <- obj$control
      control$sigma2.init <- obj$sigma2
      control$tau2.init   <- obj$tau2
      control$rho.init    <- obj$rho
      control$gamma2.init <- obj$gamma2
      control$phi.init    <- obj$phi

      ### fit model without data from ith cluster

      args <- list(yi=obj$yi, V=obj$V, W=obj$W, mods=obj$X, random=obj$random, struct=obj$struct, intercept=FALSE, data=obj$mf.r, method=obj$method, test=obj$test, dfs=obj$dfs, level=obj$level, R=obj$R, Rscale=obj$Rscale, sigma2=ifelse(obj$vc.fix$sigma2, obj$sigma2, NA), tau2=ifelse(obj$vc.fix$tau2, obj$tau2, NA), rho=ifelse(obj$vc.fix$rho, obj$rho, NA), gamma2=ifelse(obj$vc.fix$gamma2, obj$gamma2, NA), phi=ifelse(obj$vc.fix$phi, obj$phi, NA), sparse=obj$sparse, dist=obj$dist, control=control, subset=!incl, outlist=outlist)
      res <- try(suppressWarnings(.do.call(rma.mv, args)), silent=TRUE)

   } else {

      ### set values of variance/correlation components to those from the 'full' model

      args <- list(yi=obj$yi, V=obj$V, W=obj$W, mods=obj$X, random=obj$random, struct=obj$struct, intercept=FALSE, data=obj$mf.r, method=obj$method, test=obj$test, dfs=obj$dfs, level=obj$level, R=obj$R, Rscale=obj$Rscale, sigma2=obj$sigma2, tau2=obj$tau2, rho=obj$rho, gamma2=obj$gamma2, phi=obj$phi, sparse=obj$sparse, dist=obj$dist, control=obj$control, subset=!incl, outlist=outlist)
      res <- try(suppressWarnings(.do.call(rma.mv, args)), silent=TRUE)

   }

   if (inherits(res, "try-error"))
      return(list(dfbs = NA))

   ### removing a cluster could lead to a model coefficient becoming inestimable

   if (any(res$coef.na))
      return(list(dfbs = NA))

   ### elements that need to be returned

   outlist <- "vb=vb"

   ### fit model based on all data but with var/cor components fixed to those from res

   args <- list(yi=obj$yi, V=obj$V, W=obj$W, mods=obj$X, random=obj$random, struct=obj$struct, intercept=FALSE, data=obj$mf.r, method=obj$method, test=obj$test, dfs=obj$dfs, level=obj$level, R=obj$R, Rscale=obj$Rscale, sigma2=res$sigma2, tau2=res$tau2, rho=res$rho, gamma2=res$gamma2, phi=res$phi, sparse=obj$sparse, dist=obj$dist, control=obj$control, outlist=outlist)
   tmp <- try(suppressWarnings(.do.call(rma.mv, args)), silent=TRUE)

   ### compute dfbeta value(s)

   dfb <- obj$beta - res$beta

   ### compute dfbetas

   dfbs <- c(dfb / sqrt(diag(tmp$vb)))

   return(list(dfbs = dfbs))

}

############################################################################

.ddf.calc <- function(dfs, X, k, p, mf.s=NULL, mf.g=NULL, mf.h=NULL, beta=TRUE) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (beta) {

      if (is.numeric(dfs)) {
         ddf <- dfs
         if (length(ddf) == 1L)
            ddf <- rep(ddf, p)
         if (length(ddf) != p)
            stop(mstyle$stop(paste0("Length of 'dfs' argument (", length(dfs), ") does not match the number of model coefficient (", p, ").")), call.=FALSE)
      }

      if (is.character(dfs) && dfs == "residual")
         ddf <- rep(k-p, p)

      if (is.character(dfs) && dfs == "contain") {

         if (!is.null(mf.g))
            mf.g <- cbind(inner=apply(mf.g, 1, paste, collapse=" + "), outer=mf.g[ncol(mf.g)])
         if (!is.null(mf.h))
            mf.h <- cbind(inner=apply(mf.h, 1, paste, collapse=" + "), outer=mf.h[ncol(mf.h)])

         s.nlevels <- sapply(mf.s, function(x) length(unique(x))) # list() if no S
         g.nlevels <- c(length(unique(mf.g[[1]])), length(unique(mf.g[[2]]))) # c(0,0) if no G
         h.nlevels <- c(length(unique(mf.h[[1]])), length(unique(mf.h[[2]]))) # c(0,0) if no H

         #print(list(s.nlevels, g.nlevels, h.nlevels))

         s.ddf <- rep(k, p)
         g.ddf <- rep(k, p)
         h.ddf <- rep(k, p)

         for (j in seq_len(p)) {
            if (!is.null(mf.s)) {
               s.lvl <- sapply(seq_along(mf.s), function(i) all(apply(table(X[,j], mf.s[[i]]) > 0, 2, sum) == 1))
               if (any(s.lvl))
                  s.ddf[j] <- min(s.nlevels[s.lvl])
            }
            if (!is.null(mf.g)) {
               g.lvl <- sapply(seq_along(mf.g), function(i) all(apply(table(X[,j], mf.g[[i]]) > 0, 2, sum) == 1))
               if (any(g.lvl))
                  g.ddf[j] <- min(g.nlevels[g.lvl])
            }
            if (!is.null(mf.h)) {
               h.lvl <- sapply(seq_along(mf.h), function(i) all(apply(table(X[,j], mf.h[[i]]) > 0, 2, sum) == 1))
               if (any(h.lvl))
                  h.ddf[j] <- min(h.nlevels[h.lvl])
            }
         }

         #return(list(s.ddf, g.ddf, h.ddf))
         ddf <- pmin(s.ddf, g.ddf, h.ddf)
         ddf <- ddf - p

      }

      names(ddf) <- colnames(X)

   } else {

      if (is.numeric(dfs))
         dfs <- "contain"

      if (dfs == "residual")
         ddf <- k-p

      if (dfs == "contain") {

         if (!is.null(mf.s))
            ddf <- length(unique(mf.s))
         if (!is.null(mf.g))
            ddf <- length(unique(mf.g))
         if (!is.null(mf.h))
            ddf <- length(unique(mf.h))

         ddf <- ddf - p

      }

   }

   ddf[ddf < 1] <- 1

   return(ddf)

}

############################################################################
