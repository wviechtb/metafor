qqnorm.rma.uni <- function(y, type="rstandard", pch=19, envelope=TRUE,
level=y$level, bonferroni=FALSE, reps=1000, smooth=TRUE, bass=0,
label=FALSE, offset=0.3, ...) {

   if (!is.element("rma.uni", class(y)))
      stop("Argument 'y' must be an y of class \"rma.uni\".")

   na.act <- getOption("na.action")

   x <- y

   type <- match.arg(type, c("rstandard", "rstudent"))

   if (x$k == 1)
      stop("Stopped because k = 1.")

   draw.envelope <- envelope

   if (label == "out" & !envelope) {
      envelope      <- TRUE
      draw.envelope <- FALSE
   }

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

   ### construct simulation based pseudo confidence envelope

   if (envelope) {

      alpha <- ifelse(level > 1, (100-level)/100, 1-level)

      dat <- matrix(rnorm(x$k*reps), nrow=x$k, ncol=reps)

      options(na.action="na.omit")
      H <- hatvalues(x, type="matrix")
      options(na.action = na.act)

      ImH <- diag(x$k) - H
      ei  <- ImH %*% dat
      ei  <- apply(ei, 2, sort)
      if (bonferroni) {
         lb  <- apply(ei, 1, quantile,   (alpha/2)/x$k) ### consider using rowQuantiles() from matrixStats package
         ub  <- apply(ei, 1, quantile, 1-(alpha/2)/x$k) ### consider using rowQuantiles() from matrixStats package
      } else {
         lb  <- apply(ei, 1, quantile,   (alpha/2)) ### consider using rowQuantiles() from matrixStats package
         ub  <- apply(ei, 1, quantile, 1-(alpha/2)) ### consider using rowQuantiles() from matrixStats package
      }

      temp.lb <- qqnorm(lb, plot.it=FALSE)
      if (smooth)
         temp.lb <- supsmu(temp.lb$x, temp.lb$y, bass=bass)
      if (draw.envelope)
         lines(temp.lb$x, temp.lb$y, lty="dotted", ...)
         #lines(temp.lb$x, temp.lb$y, lty="12", lwd=1.5, ...)
      temp.ub <- qqnorm(ub, plot.it=FALSE)
      if (smooth)
         temp.ub <- supsmu(temp.ub$x, temp.ub$y, bass=bass)
      if (draw.envelope)
         lines(temp.ub$x, temp.ub$y, lty="dotted", ...)
         #lines(temp.ub$x, temp.ub$y, lty="12", lwd=1.5, , ...)

   }

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
         label <- ifelse(label, "out", "none")

      pos.x <- sav$x[pos]
      pos.y <- sav$y[pos]

      if (label != "none") {

         for (i in seq_len(x$k)) {

            if (label == "all") {
               text(pos.x[i], pos.y[i], slab[i], pos=ifelse(pos.x[i] >= 0, 2, 4), offset=offset, ...)
            } else {
               if (pos.y[i] < temp.lb$y[i] || pos.y[i] > temp.ub$y[i])
                  text(pos.x[i], pos.y[i], slab[i], pos=ifelse(pos.x[i] >= 0, 2, 4), offset=offset, ...)
            }

         }

      }

   }

   #########################################################################

   #if (envelope) {
   #   invisible(list(pts=sav, ci.lb=temp.lb, ci.ub=temp.ub))
   #} else {
   #   invisible(sav)
   #}

   invisible(sav)

}
