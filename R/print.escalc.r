print.escalc <- function(x, digits=attr(x,"digits"), ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(x), must="escalc")

   attr(x, "class") <- NULL

   digits <- .get.digits(digits=digits, xdigits=attr(x, "digits"), dmiss=FALSE)

   ### get positions of the variable names in the object
   ### note: if the object no longer contains a particular variable, match() returns NA;
   ### use na.omit(), so that length() is then zero (as needed for if statements below)

   yi.pos    <- na.omit(match(attr(x, "yi.names"),    names(x)))
   vi.pos    <- na.omit(match(attr(x, "vi.names"),    names(x)))
   sei.pos   <- na.omit(match(attr(x, "sei.names"),   names(x)))
   zi.pos    <- na.omit(match(attr(x, "zi.names"),    names(x)))
   pval.pos  <- na.omit(match(attr(x, "pval.names"),  names(x)))
   ci.lb.pos <- na.omit(match(attr(x, "ci.lb.names"), names(x)))
   ci.ub.pos <- na.omit(match(attr(x, "ci.ub.names"), names(x)))

   ### get rownames attribute so we can back-assign it

   rnames <- attr(x, "row.names")

   ### for printing, turn expressions into strings

   is.expr <- sapply(x, is.expression)
   x[is.expr] <- lapply(x[is.expr], as.character)

   ### turn x into a regular data frame

   x <- data.frame(x)
   rownames(x) <- rnames

   ### round variables according to the digits argument

   if (length(yi.pos) > 0L)
      x[yi.pos] <- lapply(x[yi.pos], fmtx, digits[["est"]])

   if (length(vi.pos) > 0L)
      x[vi.pos] <- lapply(x[vi.pos], fmtx, digits[["var"]])

   if (length(sei.pos) > 0L)
      x[sei.pos] <- lapply(x[sei.pos], fmtx, digits[["se"]])

   if (length(zi.pos) > 0L)
      x[zi.pos] <- lapply(x[zi.pos], fmtx, digits[["test"]])

   if (length(pval.pos) > 0L)
      x[pval.pos] <- lapply(x[pval.pos], fmtp, digits[["pval"]]) # note: using fmtp here

   if (length(ci.lb.pos) > 0L)
      x[ci.lb.pos] <- lapply(x[ci.lb.pos], fmtx, digits[["ci"]])

   if (length(ci.ub.pos) > 0L)
      x[ci.ub.pos] <- lapply(x[ci.ub.pos], fmtx, digits[["ci"]])

   ### print data frame with styling

   .space()

   tmp <- capture.output(print(x, ...))
   .print.table(tmp, mstyle)

   .space()

}
