formula.rma <- function(x, type="mods", ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(x, "rma"))
      stop(mstyle$stop("Argument 'x' must be an object of class \"rma\"."))

   type <- match.arg(type, c("mods", "yi", "scale"))

   if (type == "scale" && x$model != "rma.ls")
      stop(mstyle$stop("Can only use type='scale' for location-scale models."))

   if (type == "mods")
      return(x$formula.mods)

   if (type == "yi")
      return(x$formula.yi)

   if (type == "scale")
      return(x$formula.scale)

}
