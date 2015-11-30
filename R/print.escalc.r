print.escalc <- function(x, digits, ...) {

   if (!is.element("escalc", class(x)))
      stop("Argument 'x' must be an object of class \"escalc\".")

   attr(x, "class") <- NULL

   if (missing(digits))
      digits <- attr(x, "digits")

   if (is.null(digits))
      digits <- 4

   ### get all positions of the variable names in the object

   yi.pos    <- na.omit(match(attr(x, "yi.names"),    names(x))) ### if the object no longer contains that variable, get NA, so use na.omit()
   vi.pos    <- na.omit(match(attr(x, "vi.names"),    names(x))) ### if the object no longer contains that variable, get NA, so use na.omit()
   sei.pos   <- na.omit(match(attr(x, "sei.names"),   names(x))) ### if the object no longer contains that variable, get NA, so use na.omit()
   zi.pos    <- na.omit(match(attr(x, "zi.names"),    names(x))) ### if the object no longer contains that variable, get NA, so use na.omit()
   ci.lb.pos <- na.omit(match(attr(x, "ci.lb.names"), names(x))) ### if the object no longer contains that variable, get NA, so use na.omit()
   ci.ub.pos <- na.omit(match(attr(x, "ci.ub.names"), names(x))) ### if the object no longer contains that variable, get NA, so use na.omit()

   x <- data.frame(x)

   if (length(yi.pos) > 0)
      x[yi.pos] <- apply(x[yi.pos], 2, formatC, digits=digits, format="f")

   if (length(vi.pos) > 0)
      x[vi.pos] <- apply(x[vi.pos], 2, formatC, digits=digits, format="f")

   if (length(sei.pos) > 0)
      x[sei.pos] <- apply(x[sei.pos], 2, formatC, digits=digits, format="f")

   if (length(zi.pos) > 0)
      x[zi.pos] <- apply(x[zi.pos], 2, formatC, digits=digits, format="f")

   if (length(ci.lb.pos) > 0)
      x[ci.lb.pos] <- apply(x[ci.lb.pos], 2, formatC, digits=digits, format="f")

   if (length(ci.ub.pos) > 0)
      x[ci.ub.pos] <- apply(x[ci.ub.pos], 2, formatC, digits=digits, format="f")

   print(x, ...)

}
