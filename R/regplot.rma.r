regplot.rma <- function(x, mod, pred=TRUE, ci=TRUE, pi=FALSE, shade=TRUE,
xlim, ylim, predlim, olim, xlab, ylab, at, digits=2L,
transf, atransf, targs, level=x$level,
pch, psize, plim=c(0.5,3), col, bg, slab,
grid=FALSE, refline, label=FALSE, offset=c(1,1), labsize=1,
lcol, lwd, lty, legend=FALSE, xvals, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="rma", notav=c("rma.mh","rma.peto"))

   if (x$int.only)
      stop(mstyle$stop("Plot not applicable to intercept-only models."))

   na.act <- getOption("na.action")
   on.exit(options(na.action=na.act), add=TRUE)

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   if (missing(transf))
      transf <- FALSE

   if (missing(atransf))
      atransf <- FALSE

   transf.char  <- deparse(transf)
   atransf.char <- deparse(atransf)

   if (is.function(transf) && is.function(atransf))
      stop(mstyle$stop("Use either 'transf' or 'atransf' to specify a transformation (not both)."))

   if (exists(".darkplots") && .isTRUE(.darkplots))
      par(fg="gray95", bg="gray10", col="gray95", col.axis="gray95", col.lab="gray95", col.main="gray95", col.sub="gray95")

   mf <- match.call()

   if (missing(pch)) {
      pch <- 21
   } else {
      pch <- .getx("pch", mf=mf, data=x$data)
   }

   if (missing(psize)) {
      psize <- NULL
   } else {
      psize <- .getx("psize", mf=mf, data=x$data)
   }

   if (missing(col)) {
      col <- par("fg")
   } else {
      col <- .getx("col", mf=mf, data=x$data)
   }

   if (missing(bg)) {
      if (.is.dark(par("bg"))) {
         bg <- "gray40"
      } else {
         bg <- "darkgray"
      }
   } else {
      bg <- .getx("bg", mf=mf, data=x$data)
   }

   if (missing(slab)) {
      slab <- x$slab
   } else {
      slab <- .getx("slab", mf=mf, data=x$data)
      if (length(slab) != x$k.all)
         stop(mstyle$stop(paste0("Length of the 'slab' argument (", length(slab), ") does not correspond to the size of the original dataset (", x$k.all, ").")))
      slab <- .getsubset(slab, x$subset)
   }

   if (missing(label)) {
      label <- NULL
   } else {
      label <- .getx("label", mf=mf, data=x$data)
   }

   if (missing(targs))
      targs <- NULL

   if (missing(ylab))
      ylab <- .setlab(x$measure, transf.char, atransf.char, gentype=1, short=FALSE)

   if (missing(at))
      at <- NULL

   ### grid argument can either be a logical or a color

   if (is.logical(grid)) {
      if (.is.dark(par("bg"))) {
         gridcol <- "gray30"
      } else {
         gridcol <- "gray70"
      }
   }
   if (is.character(grid)) {
      gridcol <- grid
      grid <- TRUE
   }

   ### shade argument can either be a logical or a color (first for ci, second for pi)

   if (is.logical(shade)) {
      if (.is.dark(par("bg"))) {
         shadecol <- c("gray25", "gray15")
      } else {
         shadecol <- c("gray85", "gray95")
      }
   }
   if (is.character(shade)) {
      if (length(shade) == 1L)
         shade <- c(shade, shade)
      shadecol <- shade
      shade <- TRUE
   }

   ### copy pred to addpred (since using pred below for predicted values)

   if (inherits(pred, "list.rma")) {
      addpred <- TRUE
      if (missing(xvals))
         stop(mstyle$stop("Must specify 'xvals' argument when passing an object from predict() to 'pred'."))
      if (length(xvals) != length(pred$pred))
         stop(mstyle$stop(paste0("Length of the 'xvals' argument (", length(xvals), ") does not correspond to the number of predicted values (", length(pred$pred), ").")))
   } else {
      addpred <- pred
   }

   ### set refline to NA if it is not specified

   if (missing(refline))
      refline <- NA_real_

   ### set lcol, lty, and lwd (1 = reg line, 2 = ci bounds, 3 = pi bounds, 4 = refline)

   if (missing(lcol)) {
      if (.is.dark(par("bg"))) {
         lcol <- c(rep("gray80", 3), "gray50")
      } else {
         lcol <- c(rep(par("fg"), 3), "gray50")
      }
   } else {
      if (length(lcol) == 1L)
         lcol <- rep(lcol, 4L)
      if (length(lcol) == 2L)
         lcol <- c(lcol[c(1,2,2)], "gray50")
      if (length(lcol) == 3L)
         lcol <- c(lcol, "gray50")
   }

   if (missing(lty)) {
      lty <- c("solid", "dashed", "dotted", "solid")
   } else {
      if (length(lty) == 1L)
         lty <- rep(lty, 4L)
      if (length(lty) == 2L)
         lty <- c(lty[c(1,2,2)], "solid")
      if (length(lty) == 3L)
         lty <- c(lty, "solid")
   }

   if (missing(lwd)) {
      lwd <- c(3,1,1,2)
   } else {
      if (length(lwd) == 1L)
         lwd <- rep(lwd, 4L)
      if (length(lwd) == 2L)
         lwd <- c(lwd[c(1,2,2)], 2)
      if (length(lwd) == 3L)
         lwd <- c(lwd, 2)
   }

   level <- .level(level)

   ddd <- list(...)

   lplot    <- function(..., grep, fixed) plot(...)
   laxis    <- function(..., grep, fixed) axis(...)
   lpolygon <- function(..., grep, fixed) polygon(...)
   llines   <- function(..., grep, fixed) lines(...)
   lpoints  <- function(..., grep, fixed) points(...)
   labline  <- function(..., grep, fixed) abline(...)
   lbox     <- function(..., grep, fixed) box(...)
   ltext    <- function(..., grep, fixed) text(...)

   if (is.null(ddd$fixed)) {
      fixed <- FALSE
   } else {
      fixed <- .isTRUE(ddd$fixed)
   }

   if (is.null(ddd$grep)) {
      grep <- FALSE
   } else {
      grep <- .isTRUE(ddd$grep)
   }

   ############################################################################

   ### checks on mod argument

   if (missing(mod)) {
      if (x$p == 2L && x$int.incl) {
         mod <- 2
      } else {
         if (x$p == 1L) {
            mod <- 1
         } else {
            stop(mstyle$stop("Must specify 'mod' argument for models with multiple predictors."))
         }
      }
   }

   if (length(mod) != 1L)
      stop(mstyle$stop("Can only specify a single variable via argument 'mod'."))

   if (!(is.character(mod) || is.numeric(mod)))
      stop(mstyle$stop("Argument 'mod' must either be a character string or a scalar."))

   if (is.character(mod)) {

      if (grep) {

         mod.pos <- grep(mod, colnames(x$X), fixed=fixed)

         if (length(mod.pos) != 1L)
            stop(mstyle$stop("Could not find or uniquely identify the moderator variable specified via the 'mod' argument."))

      } else {

         mod.pos <- charmatch(mod, colnames(x$X))

         if (is.na(mod.pos) || mod.pos == 0L)
            stop(mstyle$stop("Could not find or uniquely identify the moderator variable specified via the 'mod' argument."))

      }

   } else {

      mod.pos <- round(mod)

      if (mod.pos < 1 | mod.pos > x$p)
         stop(mstyle$stop("Specified position of 'mod' variable does not exist in the model."))

   }

   ### extract the observed outcomes, corresponding sampling variances, model matrix, and ids

   yi   <- c(x$yi.f)
   vi   <- x$vi.f
   X    <- x$X.f
   ids  <- x$ids

   ### get weights

   options(na.action = "na.pass") # using na.pass to get the entire vector (length of yi.f)

   weights <- try(weights(x), silent=TRUE) # does not work for rma.glmm and rma.uni.selmodel objects

   if (inherits(weights, "try-error"))
      weights <- rep(1, x$k.f)

   options(na.action = na.act)

   ### note: pch (if vector), psize (if vector), col (if vector), bg (if vector)
   ###       must have the same length as the original dataset so we have to
   ###       apply the same subsetting (if necessary) and removing of NAs as was
   ###       done during the model fitting (note: NAs are removed further below)

   if (length(pch) == 1L)
      pch <- rep(pch, x$k.all)

   if (length(pch) != x$k.all)
      stop(mstyle$stop(paste0("Length of the 'pch' argument (", length(pch), ") does not correspond to the size of the original dataset (", x$k.all, ").")))

   pch <- .getsubset(pch, x$subset)

   psize.char <- FALSE

   if (!is.null(psize)) {

      if (is.character(psize)) {

         psize <- match.arg(psize, c("seinv", "vinv"))

         if (psize  == "seinv") {
            psize <- 1 / sqrt(vi)
         } else {
            psize <- 1 / vi
         }

         psize.char <- TRUE

      } else {

         if (length(psize) == 1L)
            psize <- rep(psize, x$k.all)

         if (length(psize) != x$k.all)
            stop(mstyle$stop(paste0("Length of the 'psize' argument (", length(psize), ") does not correspond to the size of the original dataset (", x$k.all, ").")))

         psize <- .getsubset(psize, x$subset)

      }

   }

   if (length(col) == 1L)
      col <- rep(col, x$k.all)

   if (length(col) != x$k.all)
      stop(mstyle$stop(paste0("Length of the 'col' argument (", length(col), ") does not correspond to the size of the original dataset (", x$k.all, ").")))

   col <- .getsubset(col, x$subset)

   if (length(bg) == 1L)
      bg <- rep(bg, x$k.all)

   if (length(bg) != x$k.all)
      stop(mstyle$stop(paste0("Length of the 'bg' argument (", length(bg), ") does not correspond to the size of the original dataset (", x$k.all, ").")))

   bg <- .getsubset(bg, x$subset)

   if (!is.null(label)) {

      if (is.character(label)) {

         label <- match.arg(label, c("all", "ciout", "piout"))

         if (label == "all") {

            label <- rep(TRUE, x$k.all)

            label <- .getsubset(label, x$subset)

         }

      } else if (is.logical(label)) {

         #if (!is.logical(label))
         #   stop(mstyle$stop("Argument 'label' must be a logical vector (or a single character string)."))

         if (length(label) == 1L)
            label <- rep(label, x$k.all)

         if (length(label) != x$k.all)
            stop(mstyle$stop(paste0("Length of the 'label' argument (", length(label), ") does not correspond to the size of the original dataset (", x$k.all, ").")))

         label <- .getsubset(label, x$subset)

      } else if (is.numeric(label)) {

         label <- round(label)
         label <- seq(x$k.all) %in% label

         label <- .getsubset(label, x$subset)

      }

   }

   ############################################################################

   has.na <- is.na(yi) | is.na(vi) | apply(is.na(X), 1, any)
   not.na <- !has.na

   if (any(has.na)) {

      yi      <- yi[not.na]
      vi      <- vi[not.na]
      X       <- X[not.na,,drop=FALSE]
      slab    <- slab[not.na]
      ids     <- ids[not.na]
      weights <- weights[not.na]
      pch     <- pch[not.na]
      psize   <- psize[not.na] # if NULL, remains NULL
      col     <- col[not.na]
      bg      <- bg[not.na]
      if (!is.character(label))
         label <- label[not.na]

   }

   k <- length(yi)

   ############################################################################

   ### extract values for moderator of interest

   xi <- X[,mod.pos]

   if (inherits(pred, "list.rma")) {

      xs <- xvals

      ci.lb <- pred$ci.lb
      ci.ub <- pred$ci.ub
      if (is.null(pred$pi.lb) || anyNA(pred$pi.lb)) {
         pi.lb <- pred$ci.lb
         pi.ub <- pred$ci.ub
         if (pi)
            warning(mstyle$warning("Object passed to 'pred' argument does not contain prediction interval information."), call.=FALSE)
         pi <- FALSE
      } else {
         pi.lb <- pred$pi.lb
         pi.ub <- pred$pi.ub
      }

      pred <- pred$pred

      if (!is.null(label) && is.character(label) && label %in% c("ciout", "piout")) {
         warning(mstyle$stop("Cannot label points based on the confidence/prediction interval when passing an object to 'pred'."))
         label <- NULL
      }

      yi.pred  <- NULL
      yi.ci.lb <- NULL
      yi.ci.ub <- NULL
      yi.pi.lb <- NULL
      yi.pi.ub <- NULL

   } else {

      ### get predicted values

      if (!missing(xvals)) {
         xs  <- xvals
         len <- length(xs)
         predlim <- range(xs)
      } else {
         len <- 1000
         if (missing(predlim)) {
            range.xi <- max(xi) - min(xi)
            predlim  <- c(min(xi) - .04*range.xi, max(xi) + .04*range.xi)
            xs <- seq(predlim[1], predlim[2], length=len)
         } else {
            if (length(predlim) != 2L)
               stop(mstyle$stop("Argument 'predlim' must be of length 2."))
            xs <- seq(predlim[1], predlim[2], length=len)
         }
      }

      Xnew <- rbind(colMeans(X))[rep(1,len),,drop=FALSE]
      Xnew[,mod.pos] <- xs

      if (x$int.incl)
         Xnew <- Xnew[,-1,drop=FALSE]

      tmp <- predict(x, newmods=Xnew, level=level)

      pred  <- tmp$pred
      ci.lb <- tmp$ci.lb
      ci.ub <- tmp$ci.ub
      if (is.null(tmp$pi.lb) || anyNA(tmp$pi.lb)) {
         pi.lb <- ci.lb
         pi.ub <- ci.ub
         if (pi)
            warning(mstyle$warning("Cannot draw prediction interval for the given model."), call.=FALSE)
         pi <- FALSE
      } else {
         pi.lb <- tmp$pi.lb
         pi.ub <- tmp$pi.ub
      }

      Xnew <- rbind(colMeans(X))[rep(1,k),,drop=FALSE]
      Xnew[,mod.pos] <- xi

      if (x$int.incl)
         Xnew <- Xnew[,-1,drop=FALSE]

      tmp <- predict(x, newmods=Xnew, level=level)

      yi.pred  <- tmp$pred
      yi.ci.lb <- tmp$ci.lb
      yi.ci.ub <- tmp$ci.ub
      if (is.null(tmp$pi.lb) || anyNA(tmp$pi.lb)) {
         yi.pi.lb <- yi.ci.lb
         yi.pi.ub <- yi.ci.ub
         if (!is.null(label) && is.character(label) && label == "piout") {
            warning(mstyle$warning("Cannot label points based on the prediction interval for the given model."), call.=FALSE)
            label <- NULL
         }
      } else {
         yi.pi.lb <- tmp$pi.lb
         yi.pi.ub <- tmp$pi.ub
      }

   }

   ############################################################################

   ### if requested, apply transformation to yi's and CI/PI bounds

   if (is.function(transf)) {
      if (is.null(targs)) {
         yi       <- sapply(yi, transf)
         pred     <- sapply(pred, transf)
         ci.lb    <- sapply(ci.lb, transf)
         ci.ub    <- sapply(ci.ub, transf)
         pi.lb    <- sapply(pi.lb, transf)
         pi.ub    <- sapply(pi.ub, transf)
         yi.pred  <- sapply(yi.pred, transf)
         yi.ci.lb <- sapply(yi.ci.lb, transf)
         yi.ci.ub <- sapply(yi.ci.ub, transf)
         yi.pi.lb <- sapply(yi.pi.lb, transf)
         yi.pi.ub <- sapply(yi.pi.ub, transf)
      } else {
         yi       <- sapply(yi, transf, targs)
         pred     <- sapply(pred, transf, targs)
         ci.lb    <- sapply(ci.lb, transf, targs)
         ci.ub    <- sapply(ci.ub, transf, targs)
         pi.lb    <- sapply(pi.lb, transf, targs)
         pi.ub    <- sapply(pi.ub, transf, targs)
         yi.pred  <- sapply(yi.pred, transf, targs)
         yi.ci.lb <- sapply(yi.ci.lb, transf, targs)
         yi.ci.ub <- sapply(yi.ci.ub, transf, targs)
         yi.pi.lb <- sapply(yi.pi.lb, transf, targs)
         yi.pi.ub <- sapply(yi.pi.ub, transf, targs)
      }
   }

   ### make sure order of intervals is always increasing

   tmp <- .psort(ci.lb, ci.ub)
   ci.lb <- tmp[,1]
   ci.ub <- tmp[,2]

   tmp <- .psort(pi.lb, pi.ub)
   pi.lb <- tmp[,1]
   pi.ub <- tmp[,2]

   ### apply observation/outcome limits if specified

   if (!missing(olim)) {
      if (length(olim) != 2L)
         stop(mstyle$stop("Argument 'olim' must be of length 2."))
      olim <- sort(olim)
      yi[yi < olim[1]] <- olim[1]
      yi[yi > olim[2]] <- olim[2]
      pred[pred < olim[1]] <- olim[1]
      pred[pred > olim[2]] <- olim[2]
      ci.lb[ci.lb < olim[1]] <- olim[1]
      ci.ub[ci.ub > olim[2]] <- olim[2]
      pi.lb[pi.lb < olim[1]] <- olim[1]
      pi.ub[pi.ub > olim[2]] <- olim[2]
   }

   ### set default point sizes (if not specified by user)

   if (is.null(psize) || psize.char) {
      if (length(plim) < 2L)
         stop(mstyle$stop("Argument 'plim' must be of length 2 or 3."))
      if (psize.char) {
         wi <- psize
      } else {
         wi <- sqrt(weights)
      }
      if (!is.na(plim[1]) && !is.na(plim[2])) {
         rng <- max(wi, na.rm=TRUE) - min(wi, na.rm=TRUE)
         if (rng <= .Machine$double.eps^0.5) {
            psize <- rep(1, k)
         } else {
            psize <- (wi - min(wi, na.rm=TRUE)) / rng
            psize <- (psize * (plim[2] - plim[1])) + plim[1]
         }
      }
      if (is.na(plim[1]) && !is.na(plim[2])) {
         psize <- wi / max(wi, na.rm=TRUE) * plim[2]
         if (length(plim) == 3L)
            psize[psize <= plim[3]] <- plim[3]
      }
      if (!is.na(plim[1]) && is.na(plim[2])) {
         psize <- wi / min(wi, na.rm=TRUE) * plim[1]
         if (length(plim) == 3L)
            psize[psize >= plim[3]] <- plim[3]
      }
      if (all(is.na(psize)))
         psize <- rep(1, k)
   }

   ############################################################################

   if (missing(xlab))
      xlab <- colnames(X)[mod.pos]

   if (!is.expression(xlab) && xlab == "")
      xlab <- "Moderator"

   if (missing(xlim)) {
      xlim <- range(xi)
   } else {
      if (length(xlim) != 2L)
         stop(mstyle$stop("Argument 'xlim' must be of length 2."))
   }

   if (missing(ylim)) {
      if (pi) {
         ylim <- range(c(yi, pi.lb, pi.ub))
      } else if (ci) {
         ylim <- range(c(yi, ci.lb, ci.ub))
      } else {
         ylim <- range(yi)
      }
   } else {
      if (length(ylim) != 2L)
         stop(mstyle$stop("Argument 'ylim' must be of length 2."))
   }

   ### if user has specified 'at' argument, make sure ylim actually contains the min and max 'at' values

   if (!is.null(at)) {
      ylim[1] <- min(c(ylim[1], at), na.rm=TRUE)
      ylim[2] <- max(c(ylim[2], at), na.rm=TRUE)
   }

   ############################################################################

   ### set up plot

   lplot(NA, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, yaxt="n", ...)

   ### generate y-axis positions if none are specified

   if (is.null(at)) {
      at <- axTicks(side=2)
   } else {
      at <- at[at > par("usr")[3]]
      at <- at[at < par("usr")[4]]
   }

   ### y-axis labels (apply transformation to axis labels if requested)

   at.lab <- at

   if (is.function(atransf)) {
      if (is.null(targs)) {
         at.lab <- fmtx(sapply(at.lab, atransf), digits[[1]], drop0ifint=TRUE)
      } else {
         at.lab <- fmtx(sapply(at.lab, atransf, targs), digits[[1]], drop0ifint=TRUE)
      }
   } else {
      at.lab <- fmtx(at.lab, digits[[1]], drop0ifint=TRUE)
   }

   ### add y-axis

   laxis(side=2, at=at, labels=at.lab, ...)

   ### add predicted values / CI bounds

   if (shade) {
      if (pi)
         lpolygon(c(xs, rev(xs)), c(pi.lb, rev(pi.ub)), border=NA, col=shadecol[2], ...)
      if (ci)
         lpolygon(c(xs, rev(xs)), c(ci.lb, rev(ci.ub)), border=NA, col=shadecol[1], ...)
   }

   if (ci) {
      llines(xs, ci.lb, col=lcol[2], lty=lty[2], lwd=lwd[2], ...)
      llines(xs, ci.ub, col=lcol[2], lty=lty[2], lwd=lwd[2], ...)
   }

   if (pi) {
      llines(xs, pi.lb, col=lcol[3], lty=lty[3], lwd=lwd[3], ...)
      llines(xs, pi.ub, col=lcol[3], lty=lty[3], lwd=lwd[3], ...)
   }

   ### add grid

   if (.isTRUE(grid))
      grid(col=gridcol) # grid needs to be at x and y tick positions also if using y-axis transformation

   ### add refline

   labline(h=refline, col=lcol[4], lty=lty[4], lwd=lwd[4], ...)

   if (addpred)
      llines(xs, pred, col=lcol[1], lty=lty[1], lwd=lwd[1], ...)

   ### redraw box

   lbox(...)

   ### order points by psize for plotting

   order.vec <- order(psize, decreasing=TRUE)

   xi.o    <- xi[order.vec]
   yi.o    <- yi[order.vec]
   pch.o   <- pch[order.vec]
   psize.o <- psize[order.vec]
   col.o   <- col[order.vec]
   bg.o    <- bg[order.vec]

   ### add points

   lpoints(x=xi.o, y=yi.o, pch=pch.o, col=col.o, bg=bg.o, cex=psize.o, ...)

   ### labeling of points

   if (!is.null(label)) {

      if (!is.null(label) && is.character(label) && label %in% c("ciout", "piout")) {

         if (label == "ciout") {
            label <- yi < yi.ci.lb | yi > yi.ci.ub
            label[xi < predlim[1] | xi > predlim[2]] <- FALSE
         } else {
            label <- yi < yi.pi.lb | yi > yi.pi.ub
            label[xi < predlim[1] | xi > predlim[2]] <- FALSE
         }

      }

      yrange <- ylim[2] - ylim[1]

      if (length(offset) == 2L)
         offset <- c(offset[1]/100 * yrange, offset[2]/100 * yrange, 1)
      if (length(offset) == 1L)
         offset <- c(0, offset/100 * yrange, 1)

      for (i in which(label)) {

         if (isTRUE(yi[i] > yi.pred[i])) { # yi.pred might be NULL, so use isTRUE()
            ltext(xi[i], yi[i] + offset[1] + offset[2]*psize[i]^offset[3], slab[i], cex=labsize, ...)
         } else {
            ltext(xi[i], yi[i] - offset[1] - offset[2]*psize[i]^offset[3], slab[i], cex=labsize, ...)
         }

      }

   } else {

      label <- rep(FALSE, k)

   }

   ### add legend (if requested)

   if (is.logical(legend) && isTRUE(legend))
      lpos <- "topright"

   if (is.character(legend)) {
      lpos <- legend
      legend <- TRUE
   }

   if (legend) {

      pch.l  <- NULL
      col.l  <- NULL
      bg.l   <- NULL
      lty.l  <- NULL
      lwd.l  <- NULL
      tcol.l <- NULL
      ltxt   <- NULL

      if (length(unique(pch)) == 1L && length(unique(col)) == 1L && length(unique(bg)) == 1L) {
         pch.l  <- NA
         col.l  <- NA
         bg.l   <- NA
         lty.l  <- "blank"
         lwd.l  <- NA
         tcol.l <- "transparent"
         ltxt   <- "Studies"
      }

      if (addpred) {
         pch.l  <- c(pch.l, NA)
         col.l  <- c(col.l, NA)
         bg.l   <- c(bg.l,  NA)
         lty.l  <- c(lty.l, NA)
         lwd.l  <- c(lwd.l, NA)
         tcol.l <- c(tcol.l, "transparent")
         ltxt   <- c(ltxt, "Regression Line")
      }

      if (ci) {
         pch.l  <- c(pch.l, 22)
         col.l  <- c(col.l, lcol[2])
         bg.l   <- c(bg.l,  shadecol[1])
         lty.l  <- c(lty.l, NA)
         lwd.l  <- c(lwd.l, 1)
         tcol.l <- c(tcol.l, "transparent")
         ltxt   <- c(ltxt, paste0(round(100*(1-level), digits[[1]]), "% Confidence Interval"))
      }

      if (pi) {
         pch.l  <- c(pch.l, 22)
         col.l  <- c(col.l, lcol[3])
         bg.l   <- c(bg.l,  shadecol[2])
         lty.l  <- c(lty.l, NA)
         lwd.l  <- c(lwd.l, 1)
         tcol.l <- c(tcol.l, "transparent")
         ltxt   <- c(ltxt, paste0(round(100*(1-level), digits[[1]]), "% Prediction Interval"))
      }

      if (length(ltxt) >= 1L) {
         if (.is.dark(par("bg"))) {
            legend(lpos, inset=.01, bg="gray10", pch=pch.l, col=col.l, pt.bg=bg.l, lty=lty.l, lwd=lwd.l, text.col=tcol.l, pt.cex=1.5, seg.len=3, legend=ltxt)
         } else {
            legend(lpos, inset=.01, bg="white", pch=pch.l, col=col.l, pt.bg=bg.l, lty=lty.l, lwd=lwd.l, text.col=tcol.l, pt.cex=1.5, seg.len=3, legend=ltxt)
         }
      }

      pch.l  <- NULL
      col.l  <- NULL
      bg.l   <- NULL
      lty.l  <- NULL
      lwd.l  <- NULL
      tcol.l <- NULL
      ltxt   <- NULL

      if (length(unique(pch)) == 1L && length(unique(col)) == 1L && length(unique(bg)) == 1L) {
         pch.l  <- pch[1]
         col.l  <- col[1]
         bg.l   <- bg[1]
         lty.l  <- "blank"
         lwd.l  <- 1
         tcol.l <- par("fg")
         ltxt   <- "Studies"
      }

      if (addpred) {
         pch.l  <- c(pch.l, NA)
         col.l  <- c(col.l, lcol[1])
         bg.l   <- c(bg.l,  NA)
         lty.l  <- c(lty.l, lty[1])
         lwd.l  <- c(lwd.l, lwd[1])
         tcol.l <- c(tcol.l, par("fg"))
         ltxt   <- c(ltxt, "Regression Line")
      }

      if (ci) {
         pch.l  <- c(pch.l, NA)
         col.l  <- c(col.l, lcol[2])
         bg.l   <- c(bg.l,  NA)
         lty.l  <- c(lty.l, lty[2])
         lwd.l  <- c(lwd.l, lwd[2])
         tcol.l <- c(tcol.l, par("fg"))
         ltxt   <- c(ltxt, paste0(round(100*(1-level), digits[[1]]), "% Confidence Interval"))
      }

      if (pi) {
         pch.l  <- c(pch.l, NA)
         col.l  <- c(col.l, lcol[3])
         bg.l   <- c(bg.l,  NA)
         lty.l  <- c(lty.l, lty[3])
         lwd.l  <- c(lwd.l, lwd[3])
         tcol.l <- c(tcol.l, par("fg"))
         ltxt   <- c(ltxt, paste0(round(100*(1-level), digits[[1]]), "% Prediction Interval"))
      }

      if (length(ltxt) >= 1L)
         legend(lpos, inset=.01, bg=NA, pch=pch.l, col=col.l, pt.bg=bg.l, lty=lty.l, lwd=lwd.l, text.col=tcol.l, pt.cex=1.5, seg.len=3, legend=ltxt)

   }

   ############################################################################

   sav <- data.frame(slab, ids, xi, yi, pch, psize, col, bg, label, order=order.vec)
   if (length(yi.pred) != 0L) # yi.pred might be NULL or list()
      sav$pred <- yi.pred
   attr(sav, "offset")  <- offset
   attr(sav, "labsize") <- labsize
   class(sav) <- "regplot"
   invisible(sav)

}
