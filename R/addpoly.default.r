# Note: If x and vi (or sei) are specified, the CI bounds for the polygon are
# calculated based on a normal distribution. But the Knapp and Hartung method
# may have been used to obtain vi (or sei), in which case we would want to use
# a t-distribution. Adding a corresponding argument would be a bit awkward,
# since the user would then have to specify the degrees of freedom. Instead,
# the user can just pass the CI (and PI) bounds (that were calculated with
# test="knha") directly to the function via the ci.lb and ci.ub (and pi.lb and
# pi.ub) arguments.

addpoly.default <- function(x, vi, sei, ci.lb, ci.ub, pi.lb, pi.ub,
rows=-1, level, annotate, digits, width, mlab, transf,
atransf, targs, efac, col, border, lty, fonts, cex, annosym, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   if (missing(level))
      level <- .getfromenv("forest", "level", default=95)

   if (missing(annotate))
      annotate <- .getfromenv("forest", "annotate", default=TRUE)

   if (missing(digits))
      digits <- .getfromenv("forest", "digits", default=2)

   if (missing(width))
      width <- .getfromenv("forest", "width", default=NULL)

   if (missing(transf))
      transf <- .getfromenv("forest", "transf", default=FALSE)

   if (missing(atransf))
      atransf <- .getfromenv("forest", "atransf", default=FALSE)

   if (is.function(transf) && is.function(atransf))
      stop(mstyle$stop("Use either 'transf' or 'atransf' to specify a transformation (not both)."))

   if (missing(targs))
      targs <- .getfromenv("forest", "targs", default=NULL)

   if (missing(efac))
      efac <- .getfromenv("forest", "efac", default=1)

   if (missing(fonts))
      fonts <- .getfromenv("forest", "fonts", default=NULL)

   if (missing(mlab))
      mlab <- NULL

   if (missing(col))
      col <- "black"

   if (missing(border))
      border <- "black"

   if (missing(lty))
      lty <- "dotted"

   if (missing(cex))
      cex <- .getfromenv("forest", "cex", default=NULL)

   if (missing(annosym))
      annosym <- .getfromenv("forest", "annosym", default=NULL)

   ddd <- list(...)

   if (!is.null(ddd$cr.lb))
      pi.lb <- ddd$cr.lb
   if (!is.null(ddd$cr.ub))
      pi.ub <- ddd$cr.ub

   if (length(lty) == 1L)
      lty <- c(lty, "solid")

   ### annotation symbols vector

   if (is.null(annosym))
      annosym <- c(" [", ", ", "]", "-") # 4th element for minus sign symbol
   if (length(annosym) == 3L)
      annosym <- c(annosym, "-")
   if (length(annosym) != 4L)
      stop(mstyle$stop("Argument 'annosym' must be a vector of length 3."))

   if (!is.null(ddd$lcol)) {
      lcol <- ddd$lcol
   } else {
      lcol <- "gray50"
   }

   lsegments <- function(..., cr.lb, cr.ub, addcred, pi.type, lcol) segments(...)
   ltext     <- function(..., cr.lb, cr.ub, addcred, pi.type, lcol) text(...)
   lpolygon  <- function(..., cr.lb, cr.ub, addcred, pi.type, lcol) polygon(...)

   ### set/get fonts (1st for labels, 2nd for annotations)
   ### when passing a named vector, the names are for 'family' and the values are for 'font'

   if (is.null(fonts)) {
      fonts <- rep(par("family"), 2L)
   } else {
      if (length(fonts) == 1L)
         fonts <- rep(fonts, 2L)
   }

   if (is.null(names(fonts)))
      fonts <- setNames(c(1L,1L), nm=fonts)

   par(family=names(fonts)[1], font=fonts[1])

   #########################################################################

   level <- .level(level)

   yi <- x

   if (!missing(vi) && is.function(vi)) # if vi is utils::vi()
      stop(mstyle$stop("Cannot find variable specified for 'vi' argument."), call.=FALSE)

   if (hasArg(ci.lb) && hasArg(ci.ub)) {

      ### CI bounds are specified by user

      if (length(ci.lb) != length(ci.ub))
         stop(mstyle$stop("Length of 'ci.lb' and 'ci.ub' is not the same."))

      if (missing(vi) && missing(sei)) {

         ### vi/sei not specified, so calculate vi based on CI bounds
         ### note: assumes that the CI is a symmetric Wald-type CI
         ###       computed based on a standard normal distribution

         vi <- ((ci.ub - ci.lb) / (2*qnorm(level/2, lower.tail=FALSE)))^2

      } else {

         ### vi not specified, but sei is, so set vi = sei^2

         if (missing(vi))
            vi <- sei^2

      }

      if (length(ci.lb) != length(vi))
         stop(mstyle$stop("Length of 'vi' (or 'sei') does not match length of ('ci.lb', 'ci.ub') pairs."))

   } else {

      ### CI bounds are not specified by user

      if (missing(vi)) {
         if (missing(sei)) {
            stop(mstyle$stop("Must specify either 'vi', 'sei', or ('ci.lb', 'ci.ub') pairs."))
         } else {
            vi <- sei^2
         }
      }

      if (length(yi) != length(vi))
         stop(mstyle$stop("Length of 'vi' (or 'sei') does not match length of 'x'."))

      ci.lb <- yi - qnorm(level/2, lower.tail=FALSE) * sqrt(vi)
      ci.ub <- yi + qnorm(level/2, lower.tail=FALSE) * sqrt(vi)

   }

   if (hasArg(pi.lb) && hasArg(pi.ub))  {

      if (length(pi.lb) != length(pi.ub))
         stop(mstyle$stop("Length of 'pi.lb' and 'pi.ub' is not the same."))

      if (length(pi.lb) != length(yi))
         stop(mstyle$stop("Length of ('pi.lb', 'pi.ub') does not match length of 'x'."))

   } else {

      pi.lb <- rep(NA, length(yi))
      pi.ub <- rep(NA, length(yi))

   }

   k <- length(yi)

   ### set rows value

   if (is.null(rows)) {
      rows <- -1:(-k)
   } else {
      if (length(rows) == 1L)
         rows <- rows:(rows-k+1)
   }

   if (length(rows) != length(yi))
      stop(mstyle$stop(paste0("Length of the 'rows' argument (", length(rows), ") does not correspond to the number of polygons to be plotted (", length(yi), ").")))

   ### check for NAs in yi/vi and act accordingly

   yivi.na <- is.na(yi) | is.na(vi)

   if (any(yivi.na)) {

      not.na <- !yivi.na

      if (na.act == "na.omit") {
         yi    <- yi[not.na]
         vi    <- vi[not.na]
         ci.lb <- ci.lb[not.na]
         ci.ub <- ci.ub[not.na]
         pi.lb <- pi.lb[not.na]
         pi.ub <- pi.ub[not.na]
         mlab  <- mlab[not.na]

         ### rearrange rows due to NAs being omitted

         rows.new <- rows
         rows.na  <- rows[!not.na]
         for (j in seq_len(length(rows.na))) {
            rows.new[rows <= rows.na[j]] <- rows.new[rows <= rows.na[j]] + 1
         }
         rows <- rows.new[not.na]

      }

      if (na.act == "na.fail")
         stop(mstyle$stop("Missing values in results."))

   }

   k <- length(yi)

   ### if requested, apply transformation to yi's and CI bounds

   if (is.function(transf)) {
      if (is.null(targs)) {
         yi    <- sapply(yi, transf)
         ci.lb <- sapply(ci.lb, transf)
         ci.ub <- sapply(ci.ub, transf)
         pi.lb <- sapply(pi.lb, transf)
         pi.ub <- sapply(pi.ub, transf)
      } else {
         yi    <- sapply(yi, transf, targs)
         ci.lb <- sapply(ci.lb, transf, targs)
         ci.ub <- sapply(ci.ub, transf, targs)
         pi.lb <- sapply(pi.lb, transf, targs)
         pi.ub <- sapply(pi.ub, transf, targs)
      }
   }

   ### make sure order of intervals is always increasing

   tmp <- .psort(ci.lb, ci.ub)
   ci.lb <- tmp[,1]
   ci.ub <- tmp[,2]
   tmp <- .psort(pi.lb, pi.ub)
   pi.lb <- tmp[,1]
   pi.ub <- tmp[,2]

   ### determine height of plot and set cex accordingly (if not specified)

   par.usr <- par("usr")
   height  <- par.usr[4]-par.usr[3]
   ### cannot use this since the value of k used in creating the plot is unknown
   #lheight <- strheight("O")
   #cex.adj <- ifelse(k * lheight > height * 0.8, height/(1.25 * k * lheight), 1)
   cex.adj <- min(1,20/height)
   xlim    <- par.usr[1:2]

   if (is.null(cex))
      cex <- par("cex") * cex.adj

   ### add annotations

   if (annotate) {

      if (is.function(atransf)) {
         if (is.null(targs)) {
            annotext <- cbind(sapply(yi, atransf), sapply(ci.lb, atransf), sapply(ci.ub, atransf))
         } else {
            annotext <- cbind(sapply(yi, atransf, targs), sapply(ci.lb, atransf, targs), sapply(ci.ub, atransf, targs))
         }

         ### make sure order of intervals is always increasing

         tmp <- .psort(annotext[,2:3])
         annotext[,2:3] <- tmp

      } else {

         annotext <- cbind(yi, ci.lb, ci.ub)

      }

      annotext <- .fcf(annotext, digits)
      annotext <- sub("-", annosym[4], annotext, fixed=TRUE)

      if (is.null(width)) {
         width <- apply(annotext, 2, function(x) max(nchar(x)))
      } else {
         if (length(width) == 1L)
            width <- rep(width, ncol(annotext))
      }

      for (j in seq_len(ncol(annotext))) {
         annotext[,j] <- formatC(annotext[,j], width=width[j])
      }

      annotext <- cbind(annotext[,1], annosym[1], annotext[,2], annosym[2], annotext[,3], annosym[3])
      annotext <- apply(annotext, 1, paste, collapse="")
      par(family=names(fonts)[2], font=fonts[2])
      ltext(x=xlim[2], rows, labels=annotext, pos=2, cex=cex, ...)
      par(family=names(fonts)[1], font=fonts[1])

   }

   if (length(col) == 1L)
      col <- rep(col, k)

   if (length(border) == 1L)
      border <- rep(border, k)

   if (length(lcol) == 1L)
      lcol <- rep(lcol, k)

   ### add polygon(s)

   for (i in seq_len(k)) {

      lsegments(pi.lb[i], rows[i], pi.ub[i], rows[i], lty=lty[1], col=lcol[i], ...)
      lsegments(pi.lb[i], rows[i]-(height/150)*cex*efac, pi.lb[i], rows[i]+(height/150)*cex*efac, col=lcol[i], lty=lty[2], ...)
      lsegments(pi.ub[i], rows[i]-(height/150)*cex*efac, pi.ub[i], rows[i]+(height/150)*cex*efac, col=lcol[i], lty=lty[2], ...)

      lpolygon(x=c(ci.lb[i], yi[i], ci.ub[i], yi[i]), y=c(rows[i], rows[i]+(height/100)*cex*efac, rows[i], rows[i]-(height/100)*cex*efac), col=col[i], border=border[i], ...)

      if (!is.null(mlab)) {
         if (is.list(mlab)) {
            ltext(xlim[1], rows[i], mlab[[i]], pos=4, cex=cex, ...)
         } else {
            ltext(xlim[1], rows[i], mlab[i], pos=4, cex=cex, ...)
         }
      }

   }

}
