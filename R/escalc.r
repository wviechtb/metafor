escalc <- function(measure, formula, ...) {

   if (missing(measure))
      stop("Must specify an effect size or outcome measure via the 'measure' argument.")

   ### catch cases where the user specifies a formula for the 'measure' argument
   ### (i.e., incorrectly uses the first argument to specify the formula)

   if (inherits(measure, "formula"))
      stop("Must specify the formula via the second ('formula') argument.")

   ### if formula argument is used, dispatch accordingly

   if (!missing(formula))
      UseMethod("escalc", formula)

   ### otherwise dispatch to escalc.default

   UseMethod("escalc", NULL)

}
