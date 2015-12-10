### based on stats:::update.default but with some adjustments

update.rma <- function (object, formula., ..., evaluate=TRUE) {

   if (!is.element("rma", class(object)))
      stop("Argument 'object' must be an object of class \"rma\".")

   if (is.element("robust.rma", class(object)))
      stop("Function not applicable to objects of class \"robust.rma\".")

   if (is.null(call <- getCall(object)))
      stop("Need an object with call component.")

   extras <- match.call(expand.dots = FALSE)$...

   if (!missing(formula.)) {

      if (any(is.element(c("rma.uni","rma.mv"), class(object)))) {
         if (class(object$call$yi) == "call") {
            call$yi <- update.formula(object$call$yi, formula.)
         } else {
            if (is.null(object$call$mods)) {
               object$call$mods <- ~ 1
               call$mods <- update.formula(object$call$mods, formula.)
            } else {
               if (!any(grepl("~", object$call$mods))) {
                  stop("The 'mods' argument in 'object' must be a formula for updating to work.")
               } else {
                  call$mods <- update.formula(object$call$mods, formula.)
               }
            }
         }
      }

      if (is.element("rma.glmm", class(object)))
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
