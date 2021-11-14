### based on stats:::update.default but with some adjustments

update.rma <- function(object, formula., ..., evaluate=TRUE) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(object), must="rma", notav="robust.rma")

   if (is.null(call <- getCall(object)))
      stop(mstyle$stop("Need an object with call component."))

   extras <- match.call(expand.dots = FALSE)$...

   if (!missing(formula.)) {

      if (inherits(object, c("rma.uni","rma.mv"))) {
         if (inherits(object$call$yi, "call")) {
            call$yi <- update.formula(object$call$yi, formula.)
         } else {
            if (is.null(object$call$mods)) {
               object$call$mods <- ~ 1
               call$mods <- update.formula(object$call$mods, formula.)
            } else {
               if (!any(grepl("~", object$call$mods))) {
                  stop(mstyle$stop("The 'mods' argument in 'object' must be a formula for updating to work."))
               } else {
                  call$mods <- update.formula(object$call$mods, formula.)
               }
            }
         }
      }

      if (inherits(object, "rma.glmm"))
         call$mods <- update.formula(object$call$mods, formula.)

   }

   if (length(extras)) {
      existing <- !is.na(match(names(extras), names(call)))
      for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
      if (any(!existing)) {
         call <- c(as.list(call), extras[!existing])
         call <- as.call(call)
      }
   }

   if (evaluate)
      eval(call, parent.frame())
   else call

}
