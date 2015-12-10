print.escalc <- function(x, digits, ...) {

   if (!is.element("escalc", class(x)))
      stop("Argument 'x' must be an object of class \"escalc\".")

   attr(x, "class") <- NULL

   if (missing(digits))
      digits <- attr(x, "digits")

   if (is.null(digits))
      digits <- 4

   ### get positions of the variable names in the object
   ### note: if the object no longer contains a particular variable, match() returns NA;
   ### use na.omit(), so that length() is then zero (as needed for if() statements below)

   yi.pos    <- na.omit(match(attr(x, "yi.names"),    names(x)))
   vi.pos    <- na.omit(match(attr(x, "vi.names"),    names(x)))
   sei.pos   <- na.omit(match(attr(x, "sei.names"),   names(x)))
   zi.pos    <- na.omit(match(attr(x, "zi.names"),    names(x)))
   ci.lb.pos <- na.omit(match(attr(x, "ci.lb.names"), names(x)))
   ci.ub.pos <- na.omit(match(attr(x, "ci.ub.names"), names(x)))

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
