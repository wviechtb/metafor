qqnorm.rma.uni <- function(y, type="rstandard", pch=19, envelope=TRUE,
level=y$level, bonferroni=FALSE, reps=1000, smooth=TRUE, bass=0,
label=FALSE, offset=0.3, pos=13, lty, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(y), must="rma.uni", notav=c("rma.nn", "rma.uni.selmodel"))

   na.act <- getOption("na.action")
   on.exit(options(na.action=na.act), add=TRUE)

   x <- y

   type <- match.arg(type, c("rstandard", "rstudent"))

   if (x$k == 1)
      stop(mstyle$stop("Stopped because k = 1."))

   draw.envelope <- envelope

   if (label == "out" & !envelope) {
      envelope <- TRUE
      draw.envelope <- FALSE
   }

   if (length(label) != 1L)
      stop(mstyle$stop("Argument 'label' should be of length 1."))

   if (missing(lty)) {
      lty <- c("solid", "dotted") ### 1st value = diagonal line, 2nd value = pseudo confidence envelope
   } else {
      if (length(lty) == 1L)
         lty <- c(lty, lty)
   }

   ddd <- list(...)

   lqqnorm <- function(..., seed) qqnorm(...)
   labline <- function(..., seed) abline(...)
   llines  <- function(..., seed) lines(...)
   ltext   <- function(..., seed) text(...)

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

   sav <- lqqnorm(zi, pch=pch, bty="l", ...)
   labline(a=0, b=1, lty=lty[1], ...)
   #qqline(zi, ...)
   #abline(h=0, lty="dotted", ...)
   #abline(v=0, lty="dotted", ...)

   #########################################################################

   ### construct simulation based pseudo confidence envelope

   if (envelope) {

      level <- .level(level)

      if (!is.null(ddd$seed))
         set.seed(ddd$seed)

      dat <- matrix(rnorm(x$k*reps), nrow=x$k, ncol=reps)

      options(na.action="na.omit")
      H <- hatvalues(x, type="matrix")
      options(na.action = na.act)

      ImH <- diag(x$k) - H
      ei  <- ImH %*% dat
      ei  <- apply(ei, 2, sort)
      if (bonferroni) {
         lb <- apply(ei, 1, quantile,   (level/2)/x$k) ### consider using rowQuantiles() from matrixStats package
         ub <- apply(ei, 1, quantile, 1-(level/2)/x$k) ### consider using rowQuantiles() from matrixStats package
      } else {
         lb <- apply(ei, 1, quantile,   (level/2)) ### consider using rowQuantiles() from matrixStats package
         ub <- apply(ei, 1, quantile, 1-(level/2)) ### consider using rowQuantiles() from matrixStats package
      }

      temp.lb <- qqnorm(lb, plot.it=FALSE)
      if (smooth)
         temp.lb <- supsmu(temp.lb$x, temp.lb$y, bass=bass)
      if (draw.envelope)
         llines(temp.lb$x, temp.lb$y, lty=lty[2], ...)
         #llines(temp.lb$x, temp.lb$y, lty="12", lwd=1.5, ...)
      temp.ub <- qqnorm(ub, plot.it=FALSE)
      if (smooth)
         temp.ub <- supsmu(temp.ub$x, temp.ub$y, bass=bass)
      if (draw.envelope)
         llines(temp.ub$x, temp.ub$y, lty=lty[2], ...)
         #llines(temp.ub$x, temp.ub$y, lty="12", lwd=1.5, , ...)

   }

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
               ltext(pos.x[i], pos.y[i], slab[i], pos=pos, offset=offset, ...)
            if (pos == 13)
               ltext(pos.x[i], pos.y[i], slab[i], pos=ifelse(pos.x[i]-pos.y[i] >= 0, 1, 3), offset=offset, ...)
            if (pos == 24)
               ltext(pos.x[i], pos.y[i], slab[i], pos=ifelse(pos.x[i]-pos.y[i] <= 0, 2, 4), offset=offset, ...)
               #ltext(pos.x[i], pos.y[i], slab[i], pos=ifelse(pos.x[i] >= 0, 2, 4), offset=offset, ...)
         }

      }

   } else {

      pos.x <- sav$x[ord]
      pos.y <- sav$y[ord]

      for (i in seq_len(x$k)) {

         if (pos.y[i] < temp.lb$y[i] || pos.y[i] > temp.ub$y[i]) {
            if (pos <= 4)
               ltext(pos.x[i], pos.y[i], slab[i], pos=pos, offset=offset, ...)
            if (pos == 13)
               ltext(pos.x[i], pos.y[i], slab[i], pos=ifelse(pos.x[i]-pos.y[i] >= 0, 1, 3), offset=offset, ...)
            if (pos == 24)
               ltext(pos.x[i], pos.y[i], slab[i], pos=ifelse(pos.x[i]-pos.y[i] <= 0, 2, 4), offset=offset, ...)
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
