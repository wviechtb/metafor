qqnorm.rma.mh <- function(y, type="rstandard", pch=19, label=FALSE, offset=0.3, pos=13, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(y, "rma.mh"))
      stop(mstyle$stop("Argument 'y' must be an object of class \"rma.mh\"."))

   x <- y

   type <- match.arg(type, c("rstandard", "rstudent"))

   if (x$k == 1)
      stop(mstyle$stop("Stopped because k = 1."))

   if (length(label) != 1L)
      stop(mstyle$stop("Argument 'label' should be of length 1."))

   #########################################################################

   if (type == "rstandard") {
      res    <- rstandard(x)
      not.na <- !is.na(res$z)
      zi     <- res$z[not.na]
      slab   <- res$slab[not.na]
      ord    <- order(zi)
      slab   <- slab[ord]
   } else {
      res    <- rstudent(x)
      not.na <- !is.na(res$z)
      zi     <- res$z[not.na]
      slab   <- res$slab[not.na]
      ord    <- order(zi)
      slab   <- slab[ord]
   }

   sav <- qqnorm(zi, pch=pch, bty="l", ...)
   abline(a=0, b=1, lty="solid", ...)
   #qqline(zi, ...)
   #abline(h=0, lty="dotted", ...)
   #abline(v=0, lty="dotted", ...)

   #########################################################################

   ### labeling of points

   if ((is.character(label) && label=="none") || .isFALSE(label))
      return(invisible(sav))

   if ((is.character(label) && label=="all") || .isTRUE(label))
      label <- x$k

   if (is.numeric(label)) {

      label <- round(label)

      if (label < 1 | label > x$k)
         stop(mstyle$stop("Out of range value for 'label' argument."))

      pos.x <- sav$x[ord]
      pos.y <- sav$y[ord]

      dev <- abs(pos.x - pos.y)

      for (i in seq_len(x$k)) {

         if (sum(dev > dev[i]) < label) {
            if (pos <= 4)
               text(pos.x[i], pos.y[i], slab[i], pos=pos, offset=offset, ...)
            if (pos == 13)
               text(pos.x[i], pos.y[i], slab[i], pos=ifelse(pos.x[i]-pos.y[i] >= 0, 1, 3), offset=offset, ...)
            if (pos == 24)
               text(pos.x[i], pos.y[i], slab[i], pos=ifelse(pos.x[i]-pos.y[i] <= 0, 2, 4), offset=offset, ...)
               #text(pos.x[i], pos.y[i], slab[i], pos=ifelse(pos.x[i] >= 0, 2, 4), offset=offset, ...)
         }

      }

   }

   #########################################################################

   invisible(sav)

}
