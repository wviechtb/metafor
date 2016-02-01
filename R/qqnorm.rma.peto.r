qqnorm.rma.peto <- function(y, type="rstandard", pch=19, label=FALSE, offset=0.3, ...) {

   if (!inherits(y, "rma.peto"))
      stop("Argument 'y' must be an object of class \"rma.peto\".")

   x <- y

   type <- match.arg(type, c("rstandard", "rstudent"))

   if (x$k == 1)
      stop("Stopped because k = 1.")

   if (length(label) != 1)
      stop("Argument 'label' should be of length 1.")

   #########################################################################

   if (type == "rstandard") {
      res    <- rstandard(x)
      not.na <- !is.na(res$z)
      zi     <- res$z[not.na]
      slab   <- res$slab[not.na]
      pos    <- order(zi)
      slab   <- slab[pos]
   } else {
      res    <- rstudent(x)
      not.na <- !is.na(res$z)
      zi     <- res$z[not.na]
      slab   <- res$slab[not.na]
      pos    <- order(zi)
      slab   <- slab[pos]
   }

   sav <- qqnorm(zi, pch=pch, bty="l", ...)
   abline(a=0, b=1, lty="solid", ...)
   #qqline(zi, ...)
   #abline(h=0, lty="dotted", ...)
   #abline(v=0, lty="dotted", ...)

   #########################################################################

   ### labeling of points

   if (is.numeric(label)) {

      label <- round(label)

      if (label < 1 | label > x$k)
         stop("Out of range value for 'label' argument.")

      pos.x <- sav$x[pos]
      pos.y <- sav$y[pos]

      dev <- abs(pos.x - pos.y)

      for (i in seq_len(x$k)) {

         if (sum(dev > dev[i]) < label)
            text(pos.x[i], pos.y[i], slab[i], pos=ifelse(pos.x[i] >= 0, 2, 4), offset=offset, ...)

      }

   } else {

      if (is.logical(label))
         label <- ifelse(label, "all", "none")

      pos.x <- sav$x[pos]
      pos.y <- sav$y[pos]

      if (label != "none") {

         for (i in seq_len(x$k)) {

            text(pos.x[i], pos.y[i], slab[i], pos=ifelse(pos.x[i] >= 0, 2, 4), offset=offset, ...)

         }

      }

   }

   #########################################################################

   invisible(sav)

}
