############################################################################

# function to generate all possible permutations

# .genperms <- function(k) {
#
#    v <- seq_len(k)
#
#    sub <- function(k, v) {
#       if (k==1L) {
#          matrix(v,1,k)
#       } else {
#          X  <-  NULL
#          for(i in seq_len(k)) {
#             X <- rbind(X, cbind(v[i], Recall(k-1, v[-i])))
#          }
#       X
#       }
#    }
#
#    return(sub(k, v[seq_len(k)]))
#
# }

# function to generate all possible unique permutations

.genuperms <- function(x) {

   z <- NULL

   sub <- function(x, y) {
      len.x <- length(x)
      if (len.x == 0L) {
         return(y)
      } else {
         prev.num <- 0
         for (i in seq_len(len.x)) {
            num <- x[i]
            if (num > prev.num) {
               prev.num <- num
               z <- rbind(z, Recall(x[-i], c(y,num)))
            }
         }
         return(z)
      }
   }

   return(sub(x, y=NULL))

}

.permci <- function(val, obj, j, exact, iter, progbar, level, digits, control) {

   mstyle <- .get.mstyle()

   # fit model with shifted outcome
   args <- list(yi=obj$yi - c(val*obj$X[,j]), vi=obj$vi, weights=obj$weights, mods=obj$X, intercept=FALSE, method=obj$method, weighted=obj$weighted,
                test=obj$test, tau2=ifelse(obj$tau2.fix, obj$tau2, NA), control=obj$control, skipr2=TRUE)
   res <- try(suppressWarnings(.do.call(rma.uni, args)), silent=TRUE)

   if (inherits(res, "try-error"))
      stop()

   # p-value based on permutation test
   pval <- permutest(res, exact=exact, iter=iter, progbar=FALSE, control=control)$pval[j]

   # get difference between p-value and level
   diff <- pval - level / ifelse(control$alternative == "two.sided", 1, 2)

   # show progress
   if (progbar)
      cat(mstyle$verbose(paste("pval =", fmtx(pval, digits[["pval"]]), " diff =", fmtx(diff, digits[["pval"]], flag=" "), " val =", fmtx(val, digits[["est"]], flag=" "), "\n")))

   # penalize negative differences, which should force the CI bound to correspond to a p-value of *at least* level
   diff <- ifelse(diff < 0, diff*10, diff)

   return(diff)

}

############################################################################
