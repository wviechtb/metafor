escalc <- function(measure, formula, ...) {

   if (missing(measure) || class(measure) == "formula")
      stop("Must specify an effect size or outcome measure.")

   if (missing(formula))
      formula <- NULL

   ### this will fall back on escalc.default when formula is NULL

   UseMethod("escalc", formula)

}
