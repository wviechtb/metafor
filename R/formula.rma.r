formula.rma <- function(x, type="mods", ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(x), must="rma")

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
