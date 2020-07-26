matreg <- function(y, x, R, n, V, nearPD=FALSE, level=95, digits) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (missing(digits))
      digits <- 4

   level <- ifelse(level == 0, 1, ifelse(level >= 1, (100-level)/100, ifelse(level > .5, 1-level, level)))

   p <- nrow(R)

   #y <- round(y)
   #x <- round(x)

   if (length(y) != 1L)
      stop(mstyle$stop("Argument 'y' should be a single index."))

   if (y < 1 || y > p)
      stop(mstyle$stop("Index 'y' must be >= 1 or <= ", p, "."))

   if (missing(x))
      x <- seq_len(p)[-y]

   if (anyDuplicated(x))
      stop(mstyle$stop("Argument 'x' should not contain duplicated elements."))

   if (any(x < 1 | x > p))
      stop(mstyle$stop("Indices in 'x' must be >= 1 or <= ", p, "."))

   if (y %in% x)
      stop(mstyle$stop("Index 'y' should not be an element of 'x'."))

   if (nrow(R) != ncol(R))
      stop(mstyle$stop("Argument 'R' must be a square matrix."))

   R[upper.tri(R)] <- t(R)[upper.tri(R)]

   diag(R) <- 1

   m <- length(x)

   ############################################################################

   R01 <- R[x, y, drop=FALSE]
   R11 <- R[x, x, drop=FALSE]

   #invR11 <- solve(R11)
   invR11 <- try(chol2inv(chol(R11)), silent=TRUE)

   if (inherits(invR11, "try-error")) {
      if (nearPD) {
         message(mstyle$message("Cannot invert R[x,x] matrix. Using nearPD(). Treat results with caution."))
         R11 <- as.matrix(nearPD(R11, corr=TRUE)$mat)
      } else {
         stop(mstyle$stop("Cannot invert R[x,x] matrix."))
      }
      invR11 <- try(chol2inv(chol(R11)), silent=TRUE)
      if (inherits(invR11, "try-error"))
         stop(mstyle$stop("Still cannot invert R[x,x] matrix."))
   }

   b <- invR11 %*% R01

   if (!is.null(rownames(R11))) {
      rownames(b) <- rownames(R11)
   } else {
      rownames(b) <- x
   }

   ############################################################################

   if (missing(V) && missing(n))
      stop(mstyle$stop("Either 'V' or 'n' must be specified."))

   if (!missing(V) && !missing(n))
      stop(mstyle$stop("Either 'V' or 'n' must be specified, not both."))

   if (!missing(V)) {

      if (nrow(V) != ncol(V))
         stop(mstyle$stop("Argument 'V' must be a square matrix."))

      s <- p*(p-1)/2

      if (nrow(V) != s)
         stop(mstyle$stop("Dimensions of 'V' do not match the number of elements in 'R'."))

      U <- matrix(NA, nrow=p, ncol=p)
      U[lower.tri(U)] <- seq_len(s)
      U[upper.tri(U)] <- t(U)[upper.tri(U)]
      U11 <- U[x, x, drop=FALSE]
      U01 <- U[x, y, drop=FALSE]
      u11 <- unique(c(na.omit(c(U11))))
      u01 <- c(U01)

      A <- matrix(0, nrow=m, ncol=s)

      for (a in 1:ncol(A)) {
         if (a %in% u11) {
            pos <- c(which(a == U11, arr.ind=TRUE))
            J <- matrix(0, nrow=m, ncol=m)
            J[pos[1],pos[2]] <- J[pos[2],pos[1]] <- 1
            A[,a] <- - invR11 %*% J %*% invR11 %*% R01
         }
         if (a %in% u01) {
            pos <- c(which(a == U01, arr.ind=TRUE))
            A[,a] <- invR11[,pos[1]]
         }
      }

      vb <- A %*% V %*% t(A)
      se <- sqrt(diag(vb))
      zval <- c(b / se)
      pval <- 2*pnorm(abs(zval), lower.tail=FALSE)
      crit <- qnorm(level/2, lower.tail=FALSE)
      ci.lb <- c(b - crit * se)
      ci.ub <- c(b + crit * se)

      QM <- try(t(b) %*% chol2inv(chol(vb)) %*% b, silent=TRUE)

      if (inherits(QM, "try-error"))
         QM <- NA

      QMp <- pchisq(QM, df=m, lower.tail=FALSE)

      R2 <- c(t(b) %*% R01)

      res <- list(tab = data.frame(b=b, se=se, zval=zval, pval=pval, ci.lb=ci.lb, ci.ub=ci.ub), vb=vb, R2=R2, digits=digits)

   } else {

      if (length(n) != 1L)
         stop(mstyle$stop("Argument 'n' should be a single number."))

      if (n <= m)
         stop(mstyle$stop("Cannot fit model when 'n' is equal to or less than the number of predictors."))

      df <- n - m
      sse <- 1 - c(t(b) %*% R01)
      mse <- sse / df
      vb <- mse * invR11
      se <- sqrt(diag(vb))
      tval <- c(b / se)
      pval <- 2*pt(abs(tval), df=df, lower.tail=FALSE)
      crit <- qt(level/2, df=df, lower.tail=FALSE)
      ci.lb <- c(b - crit * se)
      ci.ub <- c(b + crit * se)

      R2 <- 1 - sse
      F <- c(value = (R2 / m) / mse, df1=m, df2=df)
      Fp <- pf(F[[1]], df1=m, df2=df, lower.tail=FALSE)

      res <- list(tab = data.frame(b=b, se=se, tval=tval, pval=pval, ci.lb=ci.lb, ci.ub=ci.ub), vb=vb, R2=R2, F=F, Fp=Fp, digits=digits)

   }

   class(res) <- c("matreg")
   return(res)

}
