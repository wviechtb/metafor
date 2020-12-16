funnel.rma <- function(x, yaxis="sei", xlim, ylim, xlab, ylab,
steps=5, at, atransf, targs, digits, level=x$level, addtau2=FALSE,
type="rstandard", back="lightgray", shade="white", hlines="white",
refline, lty=3, pch=19, pch.fill=21, col, bg,
label=FALSE, offset=0.4, legend=FALSE, ci.res=1000, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(x, "rma"))
      stop(mstyle$stop("Argument 'x' must be an object of class \"rma\"."))

   na.act <- getOption("na.action")
   on.exit(options(na.action=na.act))

   yaxis <- match.arg(yaxis, c("sei", "vi", "seinv", "vinv", "ni", "ninv", "sqrtni", "sqrtninv", "lni", "wi"))
   type  <- match.arg(type,  c("rstandard", "rstudent"))

   if (missing(atransf))
      atransf <- FALSE

   atransf.char <- deparse(substitute(atransf))

   ### check if sample size information is available if plotting (some function of) the sample sizes on the y-axis

   if (is.element(yaxis, c("ni", "ninv", "sqrtni", "sqrtninv", "lni"))) {
      if (is.null(x$ni))
         stop(mstyle$stop("No sample size information stored in model object."))
      if (anyNA(x$ni))
         warning(mstyle$warning("Sample size information stored in model object \n  contains NAs. Not all studies will be plotted."), call.=FALSE)
   }

   ### set y-axis label if not specified

   if (missing(ylab)) {
      if (yaxis == "sei")
         ylab <- "Standard Error"
      if (yaxis == "vi")
         ylab <- "Variance"
      if (yaxis == "seinv")
         ylab <- "Inverse Standard Error"
      if (yaxis == "vinv")
         ylab <- "Inverse Variance"
      if (yaxis == "ni")
         ylab <- "Sample Size"
      if (yaxis == "ninv")
         ylab <- "Inverse Sample Size"
      if (yaxis == "sqrtni")
         ylab <- "Square Root Sample Size"
      if (yaxis == "sqrtninv")
         ylab <- "Inverse Square Root Sample Size"
      if (yaxis == "lni")
         ylab <- "Log Sample Size"
      if (yaxis == "wi")
         ylab <- "Weight (in %)"
   }

   if (missing(at))
      at <- NULL

   if (missing(targs))
      targs <- NULL

   ### default number of digits (if not specified)

   if (missing(digits)) {
      if (yaxis == "sei")
         digits <- c(2L,3L)
      if (yaxis == "vi")
         digits <- c(2L,3L)
      if (yaxis == "seinv")
         digits <- c(2L,3L)
      if (yaxis == "vinv")
         digits <- c(2L,3L)
      if (yaxis == "ni")
         digits <- c(2L,0L)
      if (yaxis == "ninv")
         digits <- c(2L,3L)
      if (yaxis == "sqrtni")
         digits <- c(2L,3L)
      if (yaxis == "sqrtninv")
         digits <- c(2L,3L)
      if (yaxis == "lni")
         digits <- c(2L,3L)
      if (yaxis == "wi")
         digits <- c(2L,2L)
   } else {
      if (length(digits) == 1L)     ### digits[1] for x-axis labels
         digits <- c(digits,digits) ### digits[2] for y-axis labels
   }

   ### note: digits can also be a list (e.g., digits=list(2L,3))

   if (length(lty) == 1L)
      lty <- rep(lty, 2)

   ### note: 'pch', 'col', and 'bg' are assumed to be of the same length as the original data passed to rma()
   ###       so we have to apply the same subsetting (if necessary) and removing of NAs as done during the
   ###       model fitting (note: NAs are removed further below)

   if (length(pch) == 1L) {
      pch.vec <- FALSE
      pch <- rep(pch, x$k.all)
   } else {
      pch.vec <- TRUE
   }
   if (length(pch) != x$k.all)
      stop(mstyle$stop(paste0("Number of outcomes (", x$k.all, ") does not correspond to the length of the 'pch' argument (", length(pch), ").")))

   if (!is.null(x$subset))
      pch <- pch[x$subset]

   if (!inherits(x, "rma.uni.trimfill")) {

      if (missing(col))
         col <- "black"
      if (length(col) == 1L) {
         col.vec <- FALSE
         col <- rep(col, x$k.all)
      } else {
         col.vec <- TRUE
      }
      if (length(col) != x$k.all)
         stop(mstyle$stop(paste0("Number of outcomes (", x$k.all, ") does not correspond to the length of the 'col' argument (", length(col), ").")))

      if (!is.null(x$subset))
         col <- col[x$subset]

      if (missing(bg))
         bg <- "white"
      if (length(bg) == 1L) {
         bg.vec <- FALSE
         bg <- rep(bg, x$k.all)
      } else {
         bg.vec <- TRUE
      }
      if (length(bg) != x$k.all)
         stop(mstyle$stop(paste0("Number of outcomes (", x$k.all, ") does not correspond to the length of the 'bg' argument (", length(bg), ").")))

      if (!is.null(x$subset))
         bg <- bg[x$subset]

   } else {

      ### for trimfill objects, 'col' and 'bg' are used to specify the colors of the observed and imputed data

      if (missing(col))
         col <- c("black", "black")
      if (length(col) == 1L)
         col <- c(col, "black")
      col.vec <- FALSE

      if (missing(bg))
         bg <- c("white", "white")
      if (length(bg) == 1L)
         bg <- c(bg, "white")
      bg.vec <- FALSE

   }

   if (length(label) != 1L)
      stop(mstyle$stop("Argument 'label' should be of length 1."))

   ddd <- list(...)

   lplot     <- function(..., refline2, level2, lty2) plot(...)
   labline   <- function(..., refline2, level2, lty2) abline(...)
   lsegments <- function(..., refline2, level2, lty2) segments(...)
   laxis     <- function(..., refline2, level2, lty2) axis(...)
   lpolygon  <- function(..., refline2, level2, lty2) polygon(...)
   llines    <- function(..., refline2, level2, lty2) lines(...)
   lpoints   <- function(..., refline2, level2, lty2) points(...)
   lrect     <- function(..., refline2, level2, lty2) rect(...)
   ltext     <- function(..., refline2, level2, lty2) text(...)

   if (!is.null(ddd$refline2)) {
      refline2 <- ddd$refline2
   } else {
      refline2 <- NULL
   }

   if (!is.null(ddd$level2)) {
      level2 <- ddd$level2
   } else {
      level2 <- x$level
   }

   if (!is.null(ddd$lty2)) {
      lty2 <- ddd$lty2
   } else {
      lty2 <- 3
   }

   #########################################################################

   ### get values for the x-axis (and corresponding vi, sei, and ni values)
   ### if int.only, get the observed values; otherwise, get the (deleted) residuals

   if (x$int.only) {

      if (missing(refline))
         refline <- c(x$beta)

      if (inherits(x, "rma.mv") && addtau2) {
         warning(mstyle$warning("Argument 'addtau2' ignored for 'rma.mv' models."), call.=FALSE)
         addtau2 <- FALSE
      }

      yi   <- x$yi             ### yi/vi/ni is already subsetted and NAs are removed
      vi   <- x$vi
      ni   <- x$ni             ### ni can be NULL (and there may be 'additional' NAs)
      sei  <- sqrt(vi)
      if (!is.null(x$not.na.yivi))
         x$not.na <- x$not.na.yivi
      slab <- x$slab[x$not.na] ### slab is subsetted but NAs are not removed, so still need to do this here
      pch  <- pch[x$not.na]    ### same for pch

      if (!inherits(x, "rma.uni.trimfill")) {
         col <- col[x$not.na]
         bg  <- bg[x$not.na]
      } else {
         fill <- x$fill[x$not.na]
      }

      if (missing(xlab))
         xlab <- .setlab(x$measure, transf.char="FALSE", atransf.char, gentype=1)

   } else {

      if (missing(refline))
         refline <- 0

      if (addtau2) {
         warning(mstyle$warning("Argument 'addtau2' ignored for models that contain moderators."), call.=FALSE)
         addtau2 <- FALSE
      }

      options(na.action = "na.pass") ### note: subsetted but include the NAs (there may be more
                                     ###       NAs than the ones in x$not.na (rstudent() can fail),
      if (type == "rstandard") {     ###       so we don't use x$not.na below
         res <- rstandard(x)
      } else {
         res <- rstudent(x)
      }

      options(na.action = na.act)

      ### need to check for missings here

      not.na <- !is.na(res$resid) ### vector of residuals is of size k.f and can includes NAs
      yi     <- res$resid[not.na]
      sei    <- res$se[not.na]
      ni     <- x$ni.f[not.na]    ### ni can be NULL and can still include NAs
      vi     <- sei^2
      slab   <- x$slab[not.na]
      pch    <- pch[not.na]
      col    <- col[not.na]
      bg     <- bg[not.na]

      if (missing(xlab))
         xlab <- "Residual Value"

   }

   if (inherits(x, "rma.ls") && addtau2) {
      warning(mstyle$warning("Argument 'addtau2' ignored for 'rma.ls' models."), call.=FALSE)
      addtau2 <- FALSE
   }

   tau2 <- ifelse(addtau2, x$tau2, 0)

   ### get weights (omit any NAs)

   if (yaxis == "wi") {
      options(na.action = "na.omit")
      weights <- weights(x)
      options(na.action = na.act)
   }

   #########################################################################

   ### set y-axis limits

   if (missing(ylim)) {

      ### 1st ylim value is always the lowest precision  (should be at the bottom of the plot)
      ### 2nd ylim value is always the highest precision (should be at the top of the plot)

      if (yaxis == "sei")
         ylim <- c(max(sei), 0)
      if (yaxis == "vi")
         ylim <- c(max(vi), 0)
      if (yaxis == "seinv")
         ylim <- c(min(1/sei), max(1/sei))
      if (yaxis == "vinv")
         ylim <- c(min(1/vi), max(1/vi))
      if (yaxis == "ni")
         ylim <- c(min(ni, na.rm=TRUE), max(ni, na.rm=TRUE))
      if (yaxis == "ninv")
         ylim <- c(max(1/ni, na.rm=TRUE), min(1/ni, na.rm=TRUE))
      if (yaxis == "sqrtni")
         ylim <- c(min(sqrt(ni), na.rm=TRUE), max(sqrt(ni), na.rm=TRUE))
      if (yaxis == "sqrtninv")
         ylim <- c(max(1/sqrt(ni), na.rm=TRUE), min(1/sqrt(ni), na.rm=TRUE))
      if (yaxis == "lni")
         ylim <- c(min(log(ni), na.rm=TRUE), max(log(ni), na.rm=TRUE))
      if (yaxis == "wi")
         ylim <- c(min(weights), max(weights))

      ### infinite y-axis limits can happen with "seinv" and "vinv" when one or more sampling variances are 0

      if (any(is.infinite(ylim)))
         stop(mstyle$stop("Setting 'ylim' automatically not possible (must set y-axis limits manually)."))

   } else {

      ### make sure that user supplied limits are in the right order

      if (is.element(yaxis, c("sei", "vi", "ninv", "sqrtninv")))
         ylim <- c(max(ylim), min(ylim))

      if (is.element(yaxis, c("seinv", "vinv", "ni", "sqrtni", "lni", "wi")))
         ylim <- c(min(ylim), max(ylim))

      ### make sure that user supplied limits are in the appropriate range

      if (is.element(yaxis, c("sei", "vi", "ni", "ninv", "sqrtni", "sqrtninv", "lni"))) {
         if (ylim[1] < 0 || ylim[2] < 0)
            stop(mstyle$stop("Both limits for the y axis must be >= 0."))
      }

      if (is.element(yaxis, c("seinv", "vinv"))) {
         if (ylim[1] <= 0 || ylim[2] <= 0)
            stop(mstyle$stop("Both limits for the y axis must be > 0."))
      }

      if (is.element(yaxis, c("wi"))) {
         if (ylim[1] < 0 || ylim[2] < 0)
            stop(mstyle$stop("Both limits for the y axis must be >= 0."))
      }

   }

   #########################################################################

   ### set x-axis limits

   if (is.element(yaxis, c("sei", "vi", "seinv", "vinv"))) {

      level     <- ifelse(level == 0, 1, ifelse(level >= 1, (100-level)/100, ifelse(level > .5, 1-level, level)))
      level2    <- ifelse(level2 == 0, 1, ifelse(level2 >= 1, (100-level2)/100, ifelse(level2 > .5, 1-level2, level2)))
      #level    <- ifelse(level >= 1, (100-level)/100, ifelse(level > .5, 1-level, level)) ### note: there may be multiple level values
      level.min <- min(level)                                                              ### note: smallest level is the widest CI
      lvals     <- length(level)

      ### calculate the CI bounds at the bottom of the figure (for the widest CI if there are multiple)

      if (yaxis == "sei") {
         x.lb.bot <- refline - qnorm(level.min/2, lower.tail=FALSE) * sqrt(ylim[1]^2 + tau2)
         x.ub.bot <- refline + qnorm(level.min/2, lower.tail=FALSE) * sqrt(ylim[1]^2 + tau2)
      }
      if (yaxis == "vi") {
         x.lb.bot <- refline - qnorm(level.min/2, lower.tail=FALSE) * sqrt(ylim[1] + tau2)
         x.ub.bot <- refline + qnorm(level.min/2, lower.tail=FALSE) * sqrt(ylim[1] + tau2)
      }
      if (yaxis == "seinv") {
         x.lb.bot <- refline - qnorm(level.min/2, lower.tail=FALSE) * sqrt(1/ylim[1]^2 + tau2)
         x.ub.bot <- refline + qnorm(level.min/2, lower.tail=FALSE) * sqrt(1/ylim[1]^2 + tau2)
      }
      if (yaxis == "vinv") {
         x.lb.bot <- refline - qnorm(level.min/2, lower.tail=FALSE) * sqrt(1/ylim[1] + tau2)
         x.ub.bot <- refline + qnorm(level.min/2, lower.tail=FALSE) * sqrt(1/ylim[1] + tau2)
      }

      if (missing(xlim)) {
         xlim    <- c(min(x.lb.bot,min(yi)), max(x.ub.bot,max(yi))) ### make sure x-axis not only includes widest CI, but also all yi values
         rxlim   <- xlim[2] - xlim[1]        ### calculate range of the x axis limits
         xlim[1] <- xlim[1] - (rxlim * 0.10) ### subtract 10% of range from lower x-axis bound
         xlim[2] <- xlim[2] + (rxlim * 0.10) ### add      10% of range to   upper x-axis bound
      } else {
         xlim <- sort(xlim) ### just in case the user supplies the limits in the wrong order
      }

   }

   if (is.element(yaxis, c("ni", "ninv", "sqrtni", "sqrtninv", "lni", "wi"))) {

      if (missing(xlim)) {
         xlim    <- c(min(yi), max(yi))
         rxlim   <- xlim[2] - xlim[1]        ### calculate range of the x axis limits
         xlim[1] <- xlim[1] - (rxlim * 0.10) ### subtract 10% of range from lower x-axis bound
         xlim[2] <- xlim[2] + (rxlim * 0.10) ### add      10% of range to   upper x-axis bound
      } else {
         xlim <- sort(xlim) ### just in case the user supplies the limits in the wrong order
      }

   }

   ### if user has specified 'at' argument, make sure xlim actually contains the min and max 'at' values

   if (!is.null(at)) {
      xlim[1] <- min(c(xlim[1], at), na.rm=TRUE)
      xlim[2] <- max(c(xlim[2], at), na.rm=TRUE)
   }

   #########################################################################

   ### set up plot

   lplot(NA, NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, xaxt="n", yaxt="n", bty="n", ...)

   ### add background shading

   par.usr <- par("usr")
   lrect(par.usr[1], par.usr[3], par.usr[2], par.usr[4], col=back, border=NA, ...)

   ### add y-axis

   laxis(side=2, at=seq(from=ylim[1], to=ylim[2], length.out=steps), labels=formatC(seq(from=ylim[1], to=ylim[2], length.out=steps), digits=digits[[2]], format="f", drop0trailing=ifelse(class(digits[[2]]) == "integer", TRUE, FALSE)), ...)

   ### add horizontal lines

   labline(h=seq(from=ylim[1], to=ylim[2], length.out=steps), col=hlines, ...)

   #########################################################################

   ### add CI region(s)

   if (is.element(yaxis, c("sei", "vi", "seinv", "vinv"))) {

      ### add a bit to the top/bottom ylim so that the CI region(s) fill out the entire figure

      if (yaxis == "sei") {
         rylim   <- ylim[1] - ylim[2]
         ylim[1] <- ylim[1] + (rylim * 0.10)
         ylim[2] <- max(0, ylim[2] - (rylim * 0.10))
      }

      if (yaxis == "vi") {
         rylim   <- ylim[1] - ylim[2]
         ylim[1] <- ylim[1] + (rylim * 0.10)
         ylim[2] <- max(0, ylim[2] - (rylim * 0.10))
      }

      if (yaxis == "seinv") {
         rylim   <- ylim[2] - ylim[1]
         #ylim[1] <- max(.0001, ylim[1] - (rylim * 0.10)) ### not clear how much to add to bottom
         ylim[2] <- ylim[2] + (rylim * 0.10)
      }

      if (yaxis == "vinv") {
         rylim   <- ylim[2] - ylim[1]
         #ylim[1] <- max(.0001, ylim[1] - (rylim * 0.10)) ### not clear how much to add to bottom
         ylim[2] <- ylim[2] + (rylim * 0.10)
      }

      yi.vals <- seq(from=ylim[1], to=ylim[2], length.out=ci.res)

      if (yaxis == "sei")
         vi.vals  <- yi.vals^2

      if (yaxis == "vi")
         vi.vals  <- yi.vals

      if (yaxis == "seinv")
         vi.vals  <- 1/yi.vals^2

      if (yaxis == "vinv")
         vi.vals  <- 1/yi.vals

      for (m in lvals:1) {

         ci.left  <- refline - qnorm(level[m]/2, lower.tail=FALSE) * sqrt(vi.vals + tau2)
         ci.right <- refline + qnorm(level[m]/2, lower.tail=FALSE) * sqrt(vi.vals + tau2)

         lpolygon(c(ci.left,ci.right[ci.res:1]), c(yi.vals,yi.vals[ci.res:1]), border=NA, col=shade[m], ...)
         llines(ci.left,  yi.vals, lty=lty[1], ...)
         llines(ci.right, yi.vals, lty=lty[1], ...)

      }

      if (!is.null(refline2)) {

         ci.left  <- refline2 - qnorm(level2/2, lower.tail=FALSE) * sqrt(vi.vals + tau2)
         ci.right <- refline2 + qnorm(level2/2, lower.tail=FALSE) * sqrt(vi.vals + tau2)
         llines(ci.left,  yi.vals, lty=lty2, ...)
         llines(ci.right, yi.vals, lty=lty2, ...)

      }

   }

   ### add vertical reference line
   ### use segments so that line does not extent beyond tip of CI region

   if (is.element(yaxis, c("sei", "vi", "seinv", "vinv")))
      lsegments(refline, ylim[1], refline, ylim[2], lty=lty[2], ...)

   if (is.element(yaxis, c("ni", "ninv", "sqrtni", "sqrtninv", "lni", "wi")))
      labline(v=refline, lty=lty[2], ...)

   #########################################################################

   ### add points

   xaxis.vals <- yi

   if (yaxis == "sei")
      yaxis.vals <- sei

   if (yaxis == "vi")
      yaxis.vals <- vi

   if (yaxis == "seinv")
      yaxis.vals <- 1/sei

   if (yaxis == "vinv")
      yaxis.vals <- 1/vi

   if (yaxis == "ni")
      yaxis.vals <- ni

   if (yaxis == "ninv")
      yaxis.vals <- 1/ni

   if (yaxis == "sqrtni")
      yaxis.vals <- sqrt(ni)

   if (yaxis == "sqrtninv")
      yaxis.vals <- 1/sqrt(ni)

   if (yaxis == "lni")
      yaxis.vals <- log(ni)

   if (yaxis == "wi")
      yaxis.vals <- weights

   if (!inherits(x, "rma.uni.trimfill")) {

      lpoints(xaxis.vals, yaxis.vals, pch=pch, col=col, bg=bg, ...)

   } else {

      lpoints(xaxis.vals[!fill], yaxis.vals[!fill], pch=pch,      col=col[1], bg=bg[1], ...)
      lpoints(xaxis.vals[fill],  yaxis.vals[fill],  pch=pch.fill, col=col[2], bg=bg[2], ...)

   }

   #########################################################################

   ### add L-shaped box around plot

   box(bty="l")

   ### generate x axis positions if none are specified

   if (is.null(at)) {
      at <- axTicks(side=1)
      #at <- pretty(x=c(alim[1], alim[2]), n=steps-1)
      #at <- pretty(x=c(min(ci.lb), max(ci.ub)), n=steps-1)
   } else {
      at <- at[at > par("usr")[1]]
      at <- at[at < par("usr")[2]]
   }

   at.lab <- at

   if (is.function(atransf)) {
      if (is.null(targs)) {
         at.lab <- formatC(sapply(at.lab, atransf), digits=digits[[1]], format="f", drop0trailing=ifelse(class(digits[[1]]) == "integer", TRUE, FALSE))
      } else {
         at.lab <- formatC(sapply(at.lab, atransf, targs), digits=digits[[1]], format="f", drop0trailing=ifelse(class(digits[[1]]) == "integer", TRUE, FALSE))
      }
   } else {
      at.lab <- formatC(at.lab, digits=digits[[1]], format="f", drop0trailing=ifelse(class(digits[[1]]) == "integer", TRUE, FALSE))
   }

   ### add x-axis

   laxis(side=1, at=at, labels=at.lab, ...)

   ############################################################################

   ### labeling of points

   k <- length(yi)

   if (is.numeric(label) || is.character(label) || .isTRUE(label)) {

      if (is.numeric(label)) {
         label <- round(label)
         if (label < 0)
            label <- 0
         if (label > k)
            label <- k
         label <- order(abs(yi - refline), decreasing=TRUE)[seq_len(label)]
      } else if ((is.character(label) && label == "all") || .isTRUE(label)) {
         label <- seq_len(k)
      } else if ((is.character(label) && label == "out")) {
         if (!is.element(yaxis, c("sei", "vi", "seinv", "vinv"))) {
            label <- seq_len(k)
         } else {
            label <- which(abs(yi - refline) / sqrt(vi + tau2) >= qnorm(level.min/2, lower.tail=FALSE))
         }
      } else {
         label <- NULL
      }

      for (i in label)
            ltext(yi[i], yaxis.vals[i], slab[i], pos=ifelse(yi[i]-refline >= 0, 4, 2), offset=offset, ...)

   }

   #########################################################################

   ### add legend (if requested)

   if (is.logical(legend) && isTRUE(legend))
      lpos <- "topright"

   if (is.character(legend)) {
      lpos <- legend
      legend <- TRUE
   }

   if (legend && !is.element(yaxis, c("sei", "vi", "seinv", "vinv"))) {
      legend <- FALSE
      warning(mstyle$warning("Argument 'legend' only applicable if 'yaxis' is 'sei', 'vi', 'seinv', or 'vinv'."), call.=FALSE)
   }

   if (legend) {

      level <- c(level, 0)
      lvals <- length(level)

      add.studies <- !pch.vec && !col.vec && !bg.vec # only add 'Studies' to legend if pch, col, and bg were not vectors to begin with

      scipen <- options(scipen=100)
      lchars <- max(nchar(level))-2
      options(scipen=scipen$scipen)

      pval1 <- NULL
      pval2 <- NULL
      phantom <- NULL

      ltxt <- sapply(1:lvals, function(i) {
         if (i == 1)
            return(as.expression(bquote(paste(.(pval1) < p, phantom() <= .(pval2)), list(pval1=.fcf(level[i], lchars), pval2=.fcf(1, lchars)))))
            #return(as.expression(bquote(p > .(pval), list(pval=.fcf(level[i], lchars)))))
         if (i > 1 && i < lvals)
            return(as.expression(bquote(paste(.(pval1) < p, phantom() <= .(pval2)), list(pval1=.fcf(level[i], lchars), pval2=.fcf(level[i-1], lchars)))))
         if (i == lvals)
            return(as.expression(bquote(paste(.(pval1) < p, phantom() <= .(pval2)), list(pval1=.fcf(0, lchars), pval2=.fcf(level[i-1], lchars)))))
      })

      pch.l  <- rep(22, lvals)
      col.l  <- rep("black", lvals)
      pt.cex <- rep(2, lvals)
      pt.bg  <- c(shade, back)

      if (add.studies) {
         ltxt   <- c(ltxt, expression(plain(Studies)))
         pch.l  <- c(pch.l, pch[1])
         col.l  <- c(col.l, col[1])
         pt.cex <- c(pt.cex, 1)
         pt.bg  <- c(pt.bg, bg[1])
      }

      if (inherits(x, "rma.uni.trimfill")) {
         ltxt   <- c(ltxt, expression(plain(Filled~Studies)))
         pch.l  <- c(pch.l, pch.fill[1])
         col.l  <- c(col.l, col[2])
         pt.cex <- c(pt.cex, 1)
         pt.bg  <- c(pt.bg, bg[2])
      }

      legend(lpos, inset=.01, bg="white", pch=pch.l, col=col.l, pt.cex=pt.cex, pt.bg=pt.bg, legend=ltxt)

   }

   ############################################################################

   ### prepare data frame to return

   sav <- data.frame(x=xaxis.vals, y=yaxis.vals, slab=slab, stringsAsFactors=FALSE)

   if (inherits(x, "rma.uni.trimfill"))
      sav$fill <- fill

   invisible(sav)

}
