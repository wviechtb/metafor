nobs.rma <- function(object, all=FALSE, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(object), must="rma")

   if (all) {
      n.obs <- c(studies     = object$k,
                 data        = object$k.all,
                 subset      = sum(object$subset),
                 not.na      = sum(object$not.na),
                 effective   = object$k.eff,
                 df.residual = object$k.eff - object$p.eff)
   } else {
      #n.obs <- object$k.eff - ifelse(object$method == "REML", 1, 0) * object$p.eff
      n.obs <- object$k
   }

   return(n.obs)

}
