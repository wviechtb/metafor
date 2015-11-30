confint.rma.mh <- function(object, parm, level, digits, transf, targs, ...) {

   if (!is.element("rma.mh", class(object)))
      stop("Argument 'object' must be an object of class \"rma.mh\".")

   x <- object

   if (missing(level))
      level <- x$level

   if (missing(digits))
      digits <- x$digits

   if (missing(transf))
      transf <- FALSE

   if (missing(targs))
      targs <- NULL

   #########################################################################

   alpha <- ifelse(level > 1, (100-level)/100, 1-level)
   crit  <- qnorm(alpha/2, lower.tail=FALSE)

   b <- x$b
   ci.lb <- x$b - crit * x$se
   ci.ub <- x$b + crit * x$se

   ### if requested, apply transformation function

   if (is.logical(transf) && transf && is.element(x$measure, c("OR","RR","IRR"))) ### if transf=TRUE, apply exp transformation to ORs, RRs, and IRRs
      transf <- exp

   if (is.function(transf)) {
      if (is.null(targs)) {
         b     <- sapply(b, transf)
         ci.lb <- sapply(ci.lb, transf)
         ci.ub <- sapply(ci.ub, transf)
      } else {
         b     <- sapply(b, transf, targs)
         ci.lb <- sapply(ci.lb, transf, targs)
         ci.ub <- sapply(ci.ub, transf, targs)
      }
   }

   #########################################################################

   res <- cbind(estimate=b, ci.lb, ci.ub)
   res <- list(fixed=res)
   rownames(res$fixed) <- ""

   res$digits <- digits

   class(res) <- "confint.rma"
   return(res)

}
