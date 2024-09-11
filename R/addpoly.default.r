addpoly.default     <- function(x, vi, sei, ci.lb, ci.ub, pi.lb, pi.ub,
rows=-1, level,         annotate,                predstyle, predlim, digits, width, mlab,
transf, atransf, targs, efac, col, border, lty, fonts, cex, constarea=FALSE, ...) {

   #########################################################################

   mstyle <- .get.mstyle()

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   if (missing(x))
      stop(mstyle$stop("Must specify 'x' argument."))

   k <- length(x)

   ddd <- list(...)

   if (!is.null(ddd$cr.lb))
      pi.lb <- ddd$cr.lb
   if (!is.null(ddd$cr.ub))
      pi.ub <- ddd$cr.ub

   if (missing(level))
      level <- .getfromenv("forest", "level", default=95)

   level <- .level(level)

   if (hasArg(pi.lb) && !is.null((pi.lb))) {

      pi.level <- attributes(pi.lb)$level

      if (is.null(pi.level))
         pi.level <- level

      pi.dist <- attributes(pi.lb)$dist

      if (is.null(pi.dist))
         pi.dist <- "norm"

      pi.ddf <- attributes(pi.lb)$ddf

      if (is.null(pi.ddf))
         pi.ddf <- Inf

      pi.se <- attributes(pi.lb)$se

   } else {

      pi.level <- level

   }

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

   transf.char <- deparse(transf)

   if (is.function(transf) && is.function(atransf))
      stop(mstyle$stop("Use either 'transf' or 'atransf' to specify a transformation (not both)."))

   if (missing(targs))
      targs <- .getfromenv("forest", "targs", default=NULL)

   if (missing(predstyle))
      predstyle <- "line"

   predstyle <- match.arg(predstyle, c("line", "bar", "shade", "dist"))

   if (missing(predlim))
      predlim <- NULL

   if (missing(efac))
      efac <- .getfromenv("forest", "efac", default=1)

   ### vertical expansion factor: 1st = polygon(s), 2nd = PI end lines

   ### note: forest.rma() puts 'efac' into .metafor in the order:
   ### 1st = CI/PI end lines, 2nd = arrows, 3rd = polygons, 4th = bar/shade/dist height
   ### so need to pick out the 3rd and 1st/4th in that order
   ### (so 1st = polygon, 2nd = PI end lines or bar/shade/dist height)

   if (predstyle == "line") {

      if (length(efac) == 4L)
         efac <- efac[c(3,1)]

   } else {

      if (length(efac) == 3L)
         efac <- efac[c(3,4)]

   }

   efac <- .expand1(efac, 2L)

   ### annotation symbols vector

   annosym <- .chkddd(ddd$annosym, .getfromenv("forest", "annosym", default=NULL))

   if (is.null(annosym))
      annosym <- c(" [", ", ", "]", "-", " ") # 4th element for minus sign symbol; 5th for space (in place of numbers and +)
   if (length(annosym) == 3L)
      annosym <- c(annosym, "-", " ")
   if (length(annosym) == 4L)
      annosym <- c(annosym, " ")
   if (length(annosym) != 5)
      stop(mstyle$stop("Argument 'annosym' must be a vector of length 3 (or 4 or 5)."))

   if (missing(fonts))
      fonts <- .getfromenv("forest", "fonts", default=NULL)

   if (missing(mlab))
      mlab <- NULL

   if (k == 1L) {

      if (predstyle=="dist") {
         col2 <- .coladj(par("bg","fg"), dark=0.60, light=-0.60)
      } else {
         col2 <- par("fg")
      }

      if (predstyle=="shade") {
         col3 <- .coladj(par("bg","fg"), dark=0.05, light=-0.05)
      } else {
         col3 <- .coladj(par("bg","fg"), dark=0.20, light=-0.20)
      }

      if (missing(col)) { # 1st = summary polygon, 2nd = PI line/bar / shade center / tails, 3rd = shade end / ><0 region, 4th = <>0 region
         col <- c(par("fg"), col2, col3, NA)
      } else {
         if (length(col) == 1L)
            col <- c(col, col2, col3, NA)
         if (length(col) == 2L)
            col <- c(col, col3, NA)
         if (length(col) == 3L)
            col <- c(col, NA)
      }

      if (missing(border)) {
         border <- c(par("fg"), par("fg")) # 1st = summary polygon, 2nd = bar for predstyle="bar" and distribution for predstyle="dist"
      } else {
         if (length(border) == 1L)
            border <- c(border, par("fg")) # if user only specified one value, assume it is for the summary polygon
      }

      if (missing(border)) {
         border <- c(par("fg"), par("fg")) # 1st = summary polygon, 2nd = bar for predstyle="bar"
      } else {
         if (length(border) == 1L)
            border <- c(border, par("fg"))
      }

   } else {

      if (predstyle != "line")
         stop(mstyle$stop(paste0("Can only use predstyle='", predstyle, "' when plotting a single polygon.")))

      if (missing(col))
         col <- par("fg") # color of the polygons (can be a vector)

      if (missing(border))
         border <- par("fg") # border color of the polygons (can be a vector)

   }

   lcol <- .chkddd(ddd$lcol, par("fg")) # color of PI lines (can be a vector)

   if (missing(lty))
      lty <- "dotted"

   if (length(lty) == 1L)
      lty <- c(lty, "solid") # 1st for PI line, 2nd for PI end

   if (missing(cex))
      cex <- .getfromenv("forest", "cex", default=NULL)

   if (is.null(mlab)) {
      if (predstyle == "line") {
         mlab <- rep("", k)
      } else {
         if (predstyle %in% c("bar","shade"))
            mlab <- c("", paste0("Prediction Interval", annosym[1], round(100*(1-pi.level),digits[[1]]), "% PI", annosym[3]))
         if (predstyle == "dist")
            mlab <- c("", paste0("Predictive Distribution", annosym[1], round(100*(1-pi.level),digits[[1]]), "% PI", annosym[3]))
            # note: this assumes that the PI actually is a 100*(1-pi.level) PI, which may not be true
      }
   } else {
      if (predstyle == "line") {
         mlab <- .expand1(mlab, k)
         if (length(mlab) != k)
            stop(mstyle$stop(paste0("Length of the 'mlab' argument (", length(mlab), ") does not correspond to the number of polygons to be plotted (", k, ").")))
      } else {
         if (length(mlab) == 1L && predstyle %in% c("bar","shade"))
            mlab <- c(mlab, paste0("Prediction Interval", annosym[1], round(100*(1-pi.level),digits[[1]]), "% PI", annosym[3]))
         if (length(mlab) == 1L && predstyle == "dist")
            mlab <- c(mlab, paste0("Predictive Distribution", annosym[1], round(100*(1-pi.level),digits[[1]]), "% PI", annosym[3]))
      }
   }

   lsegments <- function(..., cr.lb, cr.ub, addcred, pi.type, lcol, annosym, textpos) segments(...)
   ltext     <- function(..., cr.lb, cr.ub, addcred, pi.type, lcol, annosym, textpos) text(...)
   lpolygon  <- function(..., cr.lb, cr.ub, addcred, pi.type, lcol, annosym, textpos) polygon(...)
   lrect     <- function(..., cr.lb, cr.ub, addcred, pi.type, lcol, annosym, textpos) rect(...)
   llines    <- function(..., cr.lb, cr.ub, addcred, pi.type, lcol, annosym, textpos) lines(...)

   ### set/get fonts (1st for labels, 2nd for annotations)
   ### when passing a named vector, the names are for 'family' and the values are for 'font'

   if (is.null(fonts)) {
      fonts <- rep(par("family"), 2L)
   } else {
      fonts <- .expand1(fonts, 2L)
   }

   if (is.null(names(fonts)))
      fonts <- setNames(c(1L,1L), nm=fonts)

   par(family=names(fonts)[1], font=fonts[1])

   #########################################################################

   yi <- x

   if (!missing(vi) && is.function(vi)) # if vi is utils::vi()
      stop(mstyle$stop("Cannot find variable specified for the 'vi' argument."))

   if (hasArg(ci.lb) && hasArg(ci.ub) && !is.null(ci.lb) && !is.null(ci.ub)) {

      ### CI bounds are specified by user

      if (length(ci.lb) != length(ci.ub))
         stop(mstyle$stop("Length of 'ci.lb' and 'ci.ub' is not the same."))

      if (length(ci.lb) != k)
         stop(mstyle$stop("Length of ('ci.lb','ci.ub') does not match length of 'x'."))

      vi <- ifelse(is.na(ci.lb) | is.na(ci.ub), NA_real_, 1) # need this below for checking for NAs

   } else {

      ### CI bounds are not specified by user

      if (missing(vi)) {
         if (missing(sei)) {
            stop(mstyle$stop("Must specify either 'vi', 'sei', or ('ci.lb','ci.ub')."))
         } else {
            vi <- sei^2
         }
      }

      if (length(vi) != k)
         stop(mstyle$stop("Length of 'vi' (or 'sei') does not match length of 'x'."))

      # note: the CI bounds are calculated based on a normal distribution, but
      # the Knapp and Hartung method may have been used to obtain vi (or sei),
      # in which case we would want to use a t-distribution; instead, the user
      # should pass the CI/PI bounds (calculated with test="knha") directly to
      # the function via the ci.lb/ci.ub and pi.lb/pi.ub arguments

      ci.lb <- yi - qnorm(level/2, lower.tail=FALSE) * sqrt(vi)
      ci.ub <- yi + qnorm(level/2, lower.tail=FALSE) * sqrt(vi)

   }

   if (hasArg(pi.lb) && hasArg(pi.ub) && !is.null(pi.lb) && !is.null(pi.ub)) {

      if (length(pi.lb) != length(pi.ub))
         stop(mstyle$stop("Length of 'pi.lb' and 'pi.ub' is not the same."))

      if (length(pi.lb) != k)
         stop(mstyle$stop("Length of ('pi.lb', 'pi.ub') does not match length of 'x'."))

   } else {

      if (predstyle != "line")
         stop(mstyle$stop("Cannot draw prediction interval if 'pi.lb' and 'pi.ub' are unspecified."))

      pi.lb <- rep(NA_real_, k)
      pi.ub <- rep(NA_real_, k)

   }

   ### set rows value

   if (is.null(rows)) {
      rows <- -1:(-k)
   } else {
      if (length(rows) == 1L)
         rows <- rows:(rows-k+1)
   }

   if (predstyle == "line") {

      if (length(rows) != k)
         stop(mstyle$stop(paste0("Length of the 'rows' argument (", length(rows), ") does not correspond to the number of polygons to be plotted (", k, ").")))

   } else {

      if (length(rows) == 1L)
         rows <- c(rows, rows-1)

   }

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
         if (predstyle == "line")
            mlab  <- mlab[not.na]

         ### rearrange rows due to NAs being omitted

         if (predstyle == "line") {

            rows.new <- rows
            rows.na  <- rows[!not.na]
            for (j in seq_along(rows.na)) {
               rows.new[rows <= rows.na[j]] <- rows.new[rows <= rows.na[j]] + 1
            }
            rows <- rows.new[not.na]

         }

      }

      if (na.act == "na.fail")
         stop(mstyle$stop("Missing values in results."))

   }

   k <- length(yi)

   if (k == 0L)
      stop(mstyle$stop("Processing terminated since k = 0."))

   ### if requested, apply transformation to yi's and CI bounds

   yi.utransf <- yi

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

   ### allow adjustment of position of study labels and annotations via textpos argument

   textpos <- .chkddd(ddd$textpos, .getfromenv("forest", "textpos", default=xlim))

   if (length(textpos) != 2L)
      stop(mstyle$stop("Argument 'textpos' must be of length 2."))

   if (is.na(textpos[1]))
      textpos[1] <- xlim[1]

   if (is.na(textpos[2]))
      textpos[2] <- xlim[2]

   ### add annotations

   if (annotate) {

      if (is.function(atransf)) {
         if (is.null(targs)) {
            if (predstyle %in% c("bar","shade","dist")) {
               annotext <- cbind(sapply(c(yi, NA_real_), atransf), sapply(c(ci.lb, pi.lb), atransf), sapply(c(ci.ub, pi.ub), atransf))
            } else {
               annotext <- cbind(sapply(yi, atransf), sapply(ci.lb, atransf), sapply(ci.ub, atransf))
            }
         } else {
            if (predstyle %in% c("bar","shade","dist")) {
               annotext <- cbind(sapply(c(yi, NA_real_), atransf, targs), sapply(c(ci.lb, pi.lb), atransf, targs), sapply(c(ci.ub, pi.ub), atransf, targs))
            } else {
               annotext <- cbind(sapply(yi, atransf, targs), sapply(ci.lb, atransf, targs), sapply(ci.ub, atransf, targs))
            }
         }

         ### make sure order of intervals is always increasing

         tmp <- .psort(annotext[,2:3])
         annotext[,2:3] <- tmp

      } else {

         if (predstyle %in% c("bar","shade","dist")) {
            annotext <- cbind(c(yi, NA_real_), c(ci.lb, pi.lb), c(ci.ub, pi.ub))
         } else {
            annotext <- cbind(yi, ci.lb, ci.ub)
         }

      }

      annotext <- fmtx(annotext, digits[[1]])

      if (is.null(width)) {
         width <- apply(annotext, 2, function(x) max(nchar(x)))
      } else {
         width <- .expand1(width, ncol(annotext))
      }

      for (j in seq_len(ncol(annotext))) {
         annotext[,j] <- formatC(annotext[,j], width=width[j])
      }

      annotext <- cbind(annotext[,1], annosym[1], annotext[,2], annosym[2], annotext[,3], annosym[3])

      annotext <- apply(annotext, 1, paste, collapse="")
      if (predstyle %in% c("bar","shade","dist"))
         annotext[2] <- gsub("NA", "", annotext[2], fixed=TRUE)
      annotext <- gsub("-", annosym[4], annotext, fixed=TRUE)
      annotext <- gsub(" ", annosym[5], annotext, fixed=TRUE)

      par(family=names(fonts)[2], font=fonts[2])
      if (predstyle %in% c("bar","shade","dist")) {
         ltext(x=textpos[2], c(rows[1],rows[2]), labels=annotext, pos=2, cex=cex, ...)
      } else {
         ltext(x=textpos[2], rows, labels=annotext, pos=2, cex=cex, ...)
      }
      par(family=names(fonts)[1], font=fonts[1])

   }

   col    <- .expand1(col, k)
   border <- .expand1(border, k)
   lcol   <- .expand1(lcol, k)

   if (isTRUE(constarea)) {
      area <- (ci.ub - ci.lb) * (height/100)*cex*efac[1]
      area <- area / min(area, na.rm=TRUE)
      invarea <- 1 / area
      polheight <- (height/100)*cex*efac[1]*invarea
   } else {
      polheight <- rep((height/100)*cex*efac[1], k)
   }

   piendheight <- height / 150 * cex * efac[2]
   barheight   <- min(0.25, height / 150 * cex * efac[2])

   for (i in seq_len(k)) {

      ### add prediction interval(s)
      ### note: in contrast to forest.rma(), these do not respect 'alim' (could in principle
      ### store 'alim' in .metafor environment and extract these limits from there, but it
      ### is also nice to have the option to draw PIs without being bounded by 'alim')

      if (predstyle == "line") {

         lsegments(pi.lb[i], rows[i], pi.ub[i], rows[i], lty=lty[1], col=lcol[i], ...)
         lsegments(pi.lb[i], rows[i]-piendheight, pi.lb[i], rows[i]+piendheight, col=lcol[i], lty=lty[2], ...)
         lsegments(pi.ub[i], rows[i]-piendheight, pi.ub[i], rows[i]+piendheight, col=lcol[i], lty=lty[2], ...)

      }

      if (predstyle == "bar") {

         lrect(pi.lb[i], rows[2]-barheight, yi[i], rows[2]+barheight, col=col[2], border=border[2], ...)
         lrect(pi.ub[i], rows[2]-barheight, yi[i], rows[2]+barheight, col=col[2], border=border[2], ...)

      }

      if (predstyle %in% c("shade","dist")) {

         if (is.null(pi.se))
            stop(mstyle$stop("Cannot extract SE of the prediction interval."))

         if (is.function(transf)) {
            funlist <- lapply(list("1"=exp, "2"=transf.ztor, "3"=tanh, "4"=transf.ilogit, "5"=plogis, "6"=transf.iarcsin), deparse)
            funmatch <- sapply(funlist, identical, transf.char)
            if (!any(funmatch))
               stop(mstyle$stop("Chosen transformation not (currently) possible with this 'predstyle'."))
         }

         if (pi.dist != "norm" && pi.ddf <= 1L)
            stop(mstyle$stop("Cannot shade/draw prediction distribution when df <= 1."))

         if (predstyle == "shade") {
            x.len <- 100
            q.lo <- pi.level/2
            q.hi <- 1-pi.level/2
         } else {
            x.len <- 10000
            q.lo <- 0.0001
            q.hi <- 0.9999
         }

         if (is.null(predlim) || predstyle == "shade") {
            if (pi.dist == "norm") {
               crits <- qnorm(c(q.lo,q.hi), mean=yi.utransf[i], sd=pi.se)
               xs <- seq(crits[1], crits[2], length.out=x.len)
               ys <- dnorm(xs, mean=yi.utransf[i], sd=pi.se)
            } else {
               crits <- qt(c(q.lo,q.hi), df=pi.ddf) * pi.se + yi.utransf[i]
               xs <- seq(crits[1], crits[2], length.out=x.len)
               ys <- dt((xs - yi.utransf[i]) / pi.se, df=pi.ddf) / pi.se
            }
         } else {
            if (length(predlim) != 2L)
               stop(mstyle$stop("Argument 'predlim' must be of length 2."))
            xs <- seq(predlim[1], predlim[2], length.out=x.len)
            if (is.function(transf)) {
               if (funmatch[1])
                  xs <- suppressWarnings(log(xs))
               if (any(funmatch[2:3]))
                  xs <- suppressWarnings(atanh(xs))
               if (any(funmatch[4:5]))
                  xs <- suppressWarnings(qlogis(xs))
               if (funmatch[6])
                  xs <- suppressWarnings(transf.arcsin(xs))
               sel <- is.finite(xs) # FALSE for +-Inf and NA/NaN
               x.len <- sum(sel)
               xs <- xs[sel]
            }
            if (pi.dist == "norm") {
               ys <- dnorm(xs, mean=yi.utransf[i], sd=pi.se)
            } else {
               ys <- dt((xs - yi.utransf[i]) / pi.se, df=pi.ddf) / pi.se
            }
         }

         sel.l0 <- xs < 0
         sel.g0 <- xs > 0

         if (is.function(transf)) {
            xs <- sapply(xs, transf)
            if (funmatch[1]) {
               ys <- ys / xs
               x.lo <- 0.01
               x.hi <- Inf
            }
            if (any(funmatch[2:3])) {
               ys <- ys / (1-xs^2)
               x.lo <- -0.99
               x.hi <-  0.99
            }
            if (any(funmatch[4:5])) {
               ys <- ys / (xs*(1-xs))
               x.lo <- 0.01
               x.hi <- 0.99
            }
            if (funmatch[6]) {
               ys <- ys / (2*sqrt(xs*(1-xs)))
               x.lo <- 0.01
               x.hi <- 0.99
            }
            if (is.null(predlim)) {
               sel <- xs > x.lo & xs < x.hi
               sel.l0 <- sel.l0[sel]
               sel.g0 <- sel.g0[sel]
               ys <- ys[sel]
               xs <- xs[sel]
            }
         }

      }

      if (predstyle == "shade") {

         intensity <- 1 - (ys - min(ys)) / (max(ys) - min(ys))

         colfun <- colorRamp(c(col[2], col[3]))
         rectcol <- colfun(intensity)
         rectcol <- apply(rectcol, 1, function(x) if (anyNA(x)) NA else rgb(x[1], x[2], x[3], maxColorValue=255))

         lrect(xs[-1], rows[2]-barheight, xs[-length(xs)], rows[2]+barheight, col=rectcol, border=rectcol, ...)

      }

      if (predstyle == "dist") {

         ys <- ys / max(ys) * efac[2]

         if (is.null(predlim)) {
            sel <- ys > 0.005
         } else {
            sel <- rep(TRUE, length(ys))
         }

         xs.sel.l0 <- xs[sel.l0 & sel]
         xs.sel.g0 <- xs[sel.g0 & sel]
         ys.sel.l0 <- ys[sel.l0 & sel]
         ys.sel.g0 <- ys[sel.g0 & sel]

         xs <- xs[sel]
         ys <- ys[sel]

         drow <- rows[2] - 0.5
         ys <- ys + drow
         ys.sel.l0 <- ys.sel.l0 + drow
         ys.sel.g0 <- ys.sel.g0 + drow

         ### shade regions above/below 0

         if (yi.utransf[i] > 0) {
            lpolygon(c(xs.sel.g0,rev(xs.sel.g0)), c(ys.sel.g0,rep(drow,length(ys.sel.g0))), col=col[4], border=ifelse(is.na(col[4]),NA,border[2]), ...)
            lpolygon(c(xs.sel.l0,rev(xs.sel.l0)), c(ys.sel.l0,rep(drow,length(ys.sel.l0))), col=col[3], border=ifelse(is.na(col[3]),NA,border[2]), ...)
         } else {
            lpolygon(c(xs.sel.g0,rev(xs.sel.g0)), c(ys.sel.g0,rep(drow,length(ys.sel.g0))), col=col[3], border=ifelse(is.na(col[3]),NA,border[2]), ...)
            lpolygon(c(xs.sel.l0,rev(xs.sel.l0)), c(ys.sel.l0,rep(drow,length(ys.sel.l0))), col=col[4], border=ifelse(is.na(col[4]),NA,border[2]), ...)
         }

         ### shade tail areas

         sel <- xs <= pi.lb
         xs.sel <- xs[sel]
         ys.sel <- ys[sel]
         lpolygon(c(xs.sel,rev(xs.sel)), c(ys.sel,rep(drow,length(ys.sel))), col=col[2], border=border[2], ...)

         sel <- xs >= pi.ub
         xs.sel <- xs[sel]
         ys.sel <- ys[sel]
         lpolygon(c(xs.sel,rev(xs.sel)), c(ys.sel,rep(drow,length(ys.sel))), col=col[2], border=border[2], ...)

         ### add horizontal and distribution lines

         llines(xs, rep(drow,length(ys)), col=border[2], ...)
         llines(xs, ys, col=border[2], ...)

      }

      ### add polygon(s)

      lpolygon(x=c(ci.lb[i], yi[i], ci.ub[i], yi[i]),
               y=c(rows[i], rows[i]+polheight[i], rows[i], rows[i]-polheight[i]),
               col=col[i], border=border[i], ...)

      ### add label(s)

      if (!is.null(mlab)) {

         ltext(x=textpos[1], rows[i], mlab[[i]], pos=4, cex=cex, ...)

         if (predstyle %in% c("bar","shade","dist"))
            ltext(textpos[1], rows[2], mlab[[2]], pos=4, cex=cex, ...)

      }

   }

}
