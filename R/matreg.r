matreg <- function(y, x, R, n, V, cov=FALSE, means, ztor=FALSE, nearpd=FALSE, level=95, digits, ...) {

   mstyle <- .get.mstyle()

   if (missing(digits))
      digits <- 4

   level <- .level(level)

   ### check/process R argument

   if (missing(R))
      stop(mstyle$stop("Must specify the 'R' argument."))

   R <- as.matrix(R)

   if (nrow(R) != ncol(R))
      stop(mstyle$stop("Argument 'R' must be a square matrix."))

   if (is.null(rownames(R)))
      rownames(R) <- colnames(R)
   if (is.null(colnames(R)))
      colnames(R) <- rownames(R)

   p <- nrow(R)

   if (p <= 1L)
      stop(mstyle$stop("The 'R' matrix must be at least of size 2x2."))

   ### check/process y argument

   if (length(y) != 1L)
      stop(mstyle$stop("Argument 'y' should specify a single variable."))

   if (is.character(y)) {

      if (is.null(rownames(R)))
         stop(mstyle$stop("'R' must have dimension names when specifying a variable name for 'y'."))

      if (anyDuplicated(rownames(R)))
         stop(mstyle$stop("Dimension names of 'R' must be unique."))

      y.pos <- pmatch(y, rownames(R)) # NA if no match or there are duplicates

      if (is.na(y.pos))
         stop(mstyle$stop(paste0("Could not find variable '", y, "' in the 'R' matrix.")))

      y <- y.pos

   }

   y <- round(y)

   if (y < 1 || y > p)
      stop(mstyle$stop(paste0("Index 'y' must be >= 1 or <= ", p, ".")))

   ### check/process x argument

   if (missing(x)) # if not specified, use all other variables in R as predictors
      x <- seq_len(p)[-y]

   if (is.character(x)) {

      if (is.null(rownames(R)))
         stop(mstyle$stop("'R' must have dimension names when specifying variable names for 'x'."))

      if (anyDuplicated(rownames(R)))
         stop(mstyle$stop("Dimension names of 'R' must be unique."))

      x.pos <- pmatch(x, rownames(R)) # NA if no match or there are duplicates

      if (anyNA(x.pos))
         stop(mstyle$stop(paste0("Could not find variable", ifelse(sum(is.na(x.pos)) > 1L, "s", ""), " '", paste(x[is.na(x.pos)], collapse=", "), "' in the 'R' matrix.")))

      x <- x.pos

   }

   x <- round(x)

   if (anyDuplicated(x))
      stop(mstyle$stop("Argument 'x' should not contain duplicated elements."))

   if (any(x < 1 | x > p))
      stop(mstyle$stop(paste0("Indices in 'x' must be >= 1 or <= ", p, ".")))

   if (y %in% x)
      stop(mstyle$stop("Variable 'y' should not be an element of 'x'."))

   ### check/process V/n arguments

   if (missing(V))
      V <- NULL

   if (is.null(V) && missing(n))
      stop(mstyle$stop("Either 'V' or 'n' must be specified."))

   if (!is.null(V) && !missing(n))
      stop(mstyle$stop("Either 'V' or 'n' must be specified, not both."))

   if (cov && ztor)
      stop(mstyle$stop("Cannot use a covariance matrix as input when 'ztor=TRUE'."))

   ### get ... argument and check for extra/superfluous arguments

   ddd <- list(...)

   .chkdots(ddd, c("nearPD"))

   if (.isTRUE(ddd$nearPD))
      nearpd <- TRUE

   ############################################################################

   m <- length(x)

   R[upper.tri(R)] <- t(R)[upper.tri(R)]

   if (!is.null(V)) {

      V <- as.matrix(V)

      if (nrow(V) != ncol(V))
      stop(mstyle$stop("Argument 'V' must be a square matrix."))

      V[upper.tri(V)] <- t(V)[upper.tri(V)]

      if (cov) {
         s <- p*(p+1)/2
      } else {
         s <- p*(p-1)/2
      }

      if (nrow(V) != s)
         stop(mstyle$stop(paste0("Dimensions of 'V' (", nrow(V), "x", ncol(V), ") do not match the number of elements in 'R' (", s, ").")))

   }

   ############################################################################

   if (ztor) {

      if (!is.null(V)) {
         zij <- R[lower.tri(R)]
         Dmat <- diag(2 / (cosh(2*zij) + 1), nrow=length(zij), ncol=length(zij), names=FALSE)
         V <- Dmat %*% V %*% Dmat
      }

      R <- tanh(R)
      diag(R) <- 1

   }

   if (cov) {

      S <- R
      R <- cov2cor(R)
      sdy <- sqrt(diag(S)[y])
      sdx <- sqrt(diag(S)[x])

   } else {

      if (any(abs(R) > 1, na.rm=TRUE))
         stop(mstyle$stop("Argument 'R' must be a correlation matrix, but contains values outside [-1,1]."))

      diag(R) <- 1
      sdy <- 1

   }

   ############################################################################

   Rxy <- R[x, y, drop=FALSE]
   Rxx <- R[x, x, drop=FALSE]

   #invRxx <- solve(Rxx)
   invRxx <- try(chol2inv(chol(Rxx)), silent=TRUE)

   if (inherits(invRxx, "try-error")) {
      if (nearpd) {
         message(mstyle$message("Cannot invert R[x,x] matrix. Using nearPD(). Treat results with caution."))
         Rxx <- as.matrix(nearPD(Rxx, corr=TRUE)$mat)
      } else {
         stop(mstyle$stop("Cannot invert R[x,x] matrix."))
      }
      invRxx <- try(chol2inv(chol(Rxx)), silent=TRUE)
      if (inherits(invRxx, "try-error"))
         stop(mstyle$stop("Still cannot invert R[x,x] matrix."))
   }

   b <- invRxx %*% Rxy

   if (!is.null(rownames(Rxx))) {
      rownames(b) <- rownames(Rxx)
   } else {
      rownames(b) <- paste0("x", x)
   }

   colnames(b) <- NULL

   ############################################################################

   if (cov) {

      if (missing(means)) {

         means <- rep(0,p)
         has.means <- FALSE

      } else {

         if (length(means) != p)
            stop(mstyle$stop(paste0("Length of 'means' (", length(means), ") does not match the dimensions of 'R' (", p, "x", p, ").")))

         has.means <- TRUE

      }

   }

   ############################################################################

   if (is.null(V)) {

      # when no V matrix is specified

      if (length(n) != 1L)
         stop(mstyle$stop("Argument 'n' should be a single number."))

      df <- n - m - ifelse(cov, 1, 0)

      if (df <= 0)
         stop(mstyle$stop("Cannot fit model when 'n' is equal to or less than the number of regression coefficients."))

      sse <- 1 - c(t(b) %*% Rxy)
      mse <- sse / df
      vb  <- mse * invRxx

      R2    <- 1 - sse
      R2adj <- 1 - (1 - R2) * ((n-ifelse(cov, 1, 0)) / df)

      F  <- c(value = (R2 / m) / mse, df1=m, df2=df)
      Fp <- pf(F[[1]], df1=m, df2=df, lower.tail=FALSE)

      mse <- sdy^2 * (n-1) * (1 - R2) / df

      if (cov) {

         b <- b * sdy / sdx
         b <- rbind(means[y] - means[x] %*% b, b)
         rownames(b)[1] <- "intrcpt"

         XtX    <- (n-1) * bldiag(0,S[x,x]) + n * tcrossprod(c(1,means[x]))
         invXtX <- try(suppressWarnings(chol2inv(chol(XtX))), silent=TRUE)

         if (inherits(invXtX, "try-error")) {
            vb <- matrix(NA_real_, nrow=(m+1), ncol=(m+1))
            warning(mstyle$warning("Cannot obtain var-cov matrix of the regression coefficients."), call.=FALSE)
         } else {
            vb <- mse * invXtX
         }

         if (!has.means) {
            b[1,]  <- NA_real_
            vb[1,] <- NA_real_
            vb[,1] <- NA_real_
         }

      }

      rownames(vb) <- colnames(vb) <- rownames(b)

      se    <- sqrt(diag(vb))
      tval  <- c(b / se)
      pval  <- 2*pt(abs(tval), df=df, lower.tail=FALSE)
      crit  <- qt(level/2, df=df, lower.tail=FALSE)
      ci.lb <- c(b - crit * se)
      ci.ub <- c(b + crit * se)

      res <- list(tab = data.frame(beta=b, se=se, tval=tval, df=df, pval=pval, ci.lb=ci.lb, ci.ub=ci.ub), vb=vb, R2=R2, R2adj=R2adj, F=F, Fdf=c(m,df), Fp=Fp, mse=mse, digits=digits, test="t")

   } else {

      # when a V matrix is specified

      R2 <- c(t(b) %*% Rxy) # as in Becker & Aloe (2019); assume that this also applies for Cov matrices

      if (cov) {

         b <- b * sdy / sdx

         Rxy <- S[x, y, drop=FALSE]
         invRxx <- diag(1/sdx, nrow=m, ncol=m) %*% invRxx %*% diag(1/sdx, nrow=m, ncol=m)

         Udiag <- TRUE

      } else {

         Udiag <- FALSE

      }

      U <- matrix(NA_integer_, nrow=p, ncol=p)
      U[lower.tri(U, diag=Udiag)] <- seq_len(s)
      U[upper.tri(U, diag=Udiag)] <- t(U)[upper.tri(U, diag=Udiag)]

      Uxx <- U[x, x, drop=FALSE]
      Uxy <- U[x, y, drop=FALSE]
      uxx <- unique(c(na.omit(c(Uxx))))
      uxy <- c(Uxy)

      A <- matrix(0, nrow=m, ncol=s)

      for (a in 1:ncol(A)) {
         if (a %in% uxx) {
            pos <- c(which(a == Uxx, arr.ind=TRUE))
            J <- matrix(0, nrow=m, ncol=m)
            J[pos[1],pos[2]] <- J[pos[2],pos[1]] <- 1
            A[,a] <- - invRxx %*% J %*% invRxx %*% Rxy
         }
         if (a %in% uxy) {
            pos <- c(which(a == Uxy, arr.ind=TRUE))
            A[,a] <- invRxx[,pos[1]]
         }
      }

      vb <- A %*% V %*% t(A)

      if (cov) {

         b <- rbind(means[y] - means[x] %*% b, b)
         rownames(b)[1] <- "intrcpt"
         X <- rbind(means[x], diag(m))
         vb <- X %*% vb %*% t(X)

         if (!has.means) {
            b[1,]  <- NA_real_
            vb[1,] <- NA_real_
            vb[,1] <- NA_real_
         }

      }

      se    <- sqrt(diag(vb))
      zval  <- c(b / se)
      pval  <- 2*pnorm(abs(zval), lower.tail=FALSE)
      crit  <- qnorm(level/2, lower.tail=FALSE)
      ci.lb <- c(b - crit * se)
      ci.ub <- c(b + crit * se)

      if (cov) {
         QM <- try(as.vector(t(b[-1,,drop=FALSE]) %*% chol2inv(chol(vb[-1,-1,drop=FALSE])) %*% b[-1,,drop=FALSE]), silent=TRUE)
      } else {
         QM <- try(as.vector(t(b) %*% chol2inv(chol(vb)) %*% b), silent=TRUE)
      }

      if (inherits(QM, "try-error"))
         QM <- NA_real_

      QMp <- pchisq(QM, df=m, lower.tail=FALSE)

      rownames(vb) <- colnames(vb) <- rownames(b)

      res <- list(tab = data.frame(beta=b, se=se, zval=zval, pval=pval, ci.lb=ci.lb, ci.ub=ci.ub), vb=vb, R2=R2, QM=QM, QMdf=c(m,NA_integer_), QMp=QMp, digits=digits, test="z")

   }

   class(res) <- c("matreg")
   return(res)

}
