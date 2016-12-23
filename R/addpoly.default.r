# Note: If 'x' and 'vi' (or 'sei') are specified, the CI bounds for the
# polygon are calculated based on a standard normal distribution. But the
# Knapp and Hartung method may have been used to obtain 'vi' (or 'sei'),
# in which case we would want to use a t-distribution. But adding a
# corresponding argument would be a bit awkward, since the user would then
# also have to specify the degrees of freedom. Instead, the user can just
# pass the CI bounds (that were calculated with 'test="knha"') directly to
# the function via the 'ci.lb' and 'ci.ub' argument.

addpoly.default <- function(x, vi, sei, ci.lb, ci.ub, rows=-1, level=95,
annotate=TRUE, digits=2, width, mlab, transf, atransf, targs,
efac=1, col, border, cex, ...) {

   #########################################################################

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop("Unknown 'na.action' specified under options().")

   if (missing(transf))
      transf <- FALSE

   if (missing(atransf))
      atransf <- FALSE

   if (is.function(transf) && is.function(atransf))
      stop("Use either 'transf' or 'atransf' to specify a transformation (not both).")

   if (missing(targs))
      targs <- NULL

   if (missing(mlab))
      mlab <- NULL

   if (missing(cex))
      cex <- NULL

   if (missing(col))
      col <- "black"

   if (missing(border))
      border <- "black"

   #########################################################################

   alpha <- ifelse(level > 1, (100-level)/100, 1-level)

   yi <- x

   if (hasArg(ci.lb) && hasArg(ci.ub)) {

      ### CI bounds are specified by user

      if (length(ci.lb) != length(ci.ub))
         stop("Length of 'ci.lb' and 'ci.ub' do not match.")

      if (missing(vi) && missing(sei)) {

         ### vi/sei not specified, so calculate vi based on CI bounds
         ### note: technically assumes that the CI is a symmetric Wald-type
         ###       CI computed based on a standard normal distribution

         vi <- ((ci.ub - ci.lb) / (2*qnorm(alpha/2, lower.tail=FALSE)))^2

      } else {

         ### vi not specified, but sei is, so set vi = sei^2

         if (missing(vi))
            vi <- sei^2

      }

      if (length(ci.lb) != length(vi))
         stop("Length of 'vi' (or 'sei') does not match length of ('ci.lb', 'ci.ub') pairs.")

   } else {

      ### CI bounds are not specified by user

      if (missing(vi)) {
         if (missing(sei)) {
            stop("Must specify either 'vi', 'sei', or ('ci.lb', 'ci.ub') pairs.")
         } else {
            vi <- sei^2
            ci.lb <- yi - qnorm(alpha/2, lower.tail=FALSE) * sei
            ci.ub <- yi + qnorm(alpha/2, lower.tail=FALSE) * sei
         }
      } else {
         ci.lb <- yi - qnorm(alpha/2, lower.tail=FALSE) * sqrt(vi)
         ci.ub <- yi + qnorm(alpha/2, lower.tail=FALSE) * sqrt(vi)
      }

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
      stop("Number of outcomes does not correspond to the length of the 'rows' argument.")

   ### check for NAs in yi/vi and act accordingly

   yivi.na <- is.na(yi) | is.na(vi)

   if (any(yivi.na)) {

      not.na <- !yivi.na

      if (na.act == "na.omit") {
         yi    <- yi[not.na]
         vi    <- vi[not.na]
         ci.lb <- ci.lb[not.na]
         ci.ub <- ci.ub[not.na]
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
         stop("Missing values in results.")

   }

   k <- length(yi)

   ### if requested, apply transformation to yi's and CI bounds

   if (is.function(transf)) {
      if (is.null(targs)) {
         yi    <- sapply(yi, transf)
         ci.lb <- sapply(ci.lb, transf)
         ci.ub <- sapply(ci.ub, transf)
      } else {
         yi    <- sapply(yi, transf, targs)
         ci.lb <- sapply(ci.lb, transf, targs)
         ci.ub <- sapply(ci.ub, transf, targs)
      }
   }

   ### make sure order of intervals is always increasing

   tmp <- .psort(ci.lb, ci.ub)
   ci.lb <- tmp[,1]
   ci.ub <- tmp[,2]

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

      annotext <- formatC(annotext, format="f", digits=digits)

      if (missing(width) || is.null(width)) {
         width <- apply(annotext, 2, function(x) max(nchar(x)))
      } else {
         if (length(width) == 1L)
            width <- rep(width, ncol(annotext))
      }

      for (j in seq_len(ncol(annotext))) {
         annotext[,j] <- formatC(annotext[,j], width=width[j])
      }

      annotext <- cbind(annotext[,1], " [", annotext[,2], ", ", annotext[,3], "]")
      annotext <- apply(annotext, 1, paste, collapse="")
      text(x=xlim[2], rows, labels=annotext, pos=2, cex=cex, ...)

   }

   ### add polygon(s)

   for (i in seq_len(k)) {

      polygon(x=c(ci.lb[i], yi[i], ci.ub[i], yi[i]), y=c(rows[i], rows[i]+(height/100)*cex*efac, rows[i], rows[i]-(height/100)*cex*efac), col=col, border=border, ...)

      if (!is.null(mlab))
         text(xlim[1], rows[i], mlab[i], pos=4, cex=cex, ...)

   }

}
