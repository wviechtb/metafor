funnel.default <- function(x, vi, sei, ni, subset, yaxis="sei", xlim, ylim, xlab, ylab,
steps=5, at, atransf, targs, digits, level=95,
back="lightgray", shade="white", hlines="white",
refline=0, pch=19, pch.fill=21, ci.res=1000, ...) {

   #########################################################################

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop("Unknown 'na.action' specified under options().")

   yaxis <- match.arg(yaxis, c("sei", "vi", "seinv", "vinv", "ni", "ninv", "sqrtni", "sqrtninv", "lni", "wi"))

   if (missing(atransf))
      atransf <- FALSE

   atransf.char <- deparse(substitute(atransf))

   yi <- x

   ### check if sample size information is available if plotting (some function of) sample sizes

   if (missing(ni))
      ni <- NULL

   if (is.element(yaxis, c("ni", "ninv", "sqrtni", "sqrtninv", "lni"))) {
      if (is.null(ni) && !is.null(attr(yi, "ni")))
         ni <- attr(yi, "ni")
      if (!is.null(ni) && length(ni) != length(yi))
         ni <- NULL
      if (is.null(ni))
         stop("No sample size information available.")
   }

   ### check if sampling variances and/or standard errors are available

   if (missing(vi))
      vi <- NULL

   if (missing(sei))
      sei <- NULL

   if (is.null(vi)) {
      if (!is.null(sei))
         vi <- sei^2
   }
   if (is.null(sei)) {
      if (!is.null(vi))
         sei <- sqrt(vi)
   }

   if (is.element(yaxis, c("sei", "vi", "seinv", "vinv", "wi")) && is.null(vi))
      stop("Need to specify 'vi' or 'sei' argument.")

   ### set negative variances and/or standard errors to 0

   if (!is.null(vi))
      vi[vi < 0] <- 0
   if (!is.null(sei))
      sei[sei < 0] <- 0

   ### get slab from attributes of yi; if not available or it doesn't have the right length, set slab <- 1:k

   slab <- attr(yi, "slab")

   if (is.null(slab) || length(slab) != length(yi))
      slab <- seq_along(yi)

   ### set y-axis label if not specified

   if (missing(ylab)) {
      if(yaxis == "sei")
         ylab <- "Standard Error"
      if(yaxis == "vi")
         ylab <- "Variance"
      if(yaxis == "seinv")
         ylab <- "Inverse Standard Error"
      if(yaxis == "vinv")
         ylab <- "Inverse Variance"
      if(yaxis == "ni")
         ylab <- "Sample Size"
      if(yaxis == "ninv")
         ylab <- "Inverse Sample Size"
      if(yaxis == "sqrtni")
         ylab <- "Square Root Sample Size"
      if(yaxis == "sqrtninv")
         ylab <- "Inverse Square Root Sample Size"
      if(yaxis == "lni")
         ylab <- "Log Sample Size"
      if(yaxis == "wi")
         ylab <- "Weight (in %)"
   }

   if (missing(at))
      at <- NULL

   if (missing(targs))
      targs <- NULL

   ### default number of digits (if not specified)

   if (missing(digits)) {
      if(yaxis == "sei")
         digits <- c(2L,3L)
      if(yaxis == "vi")
         digits <- c(2L,3L)
      if(yaxis == "seinv")
         digits <- c(2L,3L)
      if(yaxis == "vinv")
         digits <- c(2L,3L)
      if(yaxis == "ni")
         digits <- c(2L,0L)
      if(yaxis == "ninv")
         digits <- c(2L,3L)
      if(yaxis == "sqrtni")
         digits <- c(2L,3L)
      if(yaxis == "sqrtninv")
         digits <- c(2L,3L)
      if(yaxis == "lni")
         digits <- c(2L,3L)
      if(yaxis == "wi")
         digits <- c(2L,2L)
   } else {
      if (length(digits) == 1L)     ### digits[1] for x-axis labels
         digits <- c(digits,digits) ### digits[2] for y-axis labels
   }

   ### note: digits can also be a list (e.g., digits=list(2L,3))

   #########################################################################

   ### if a subset of studies is specified

   if (!is.null(subset)) {
      yi   <- yi[subset]
      vi   <- vi[subset]
      sei  <- sei[subset]
      ni   <- ni[subset]
      slab <- slab[subset]
   }

   ### check for NAs and act accordingly

   has.na <- is.na(yi) | if (is.element(yaxis, c("vi", "vinv"))) is.na(vi) else FALSE | if (is.element(yaxis, c("sei", "seinv"))) is.na(vi) else FALSE | if (is.element(yaxis, c("ni", "ninv", "sqrtni", "sqrtninv", "lni"))) is.na(ni) else FALSE

   if (any(has.na)) {

      not.na <- !has.na

      if (na.act == "na.omit" || na.act == "na.exclude" || na.act == "na.pass") {

         yi   <- yi[not.na]
         vi   <- vi[not.na]
         sei  <- sei[not.na]
         ni   <- ni[not.na]
         slab <- slab[not.na]

      }

      if (na.act == "na.fail")
         stop("Missing values in data.")

   }

   if (missing(xlab))
      xlab <- .setlab(attr(yi, "measure"), transf.char="FALSE", atransf.char, gentype=1)

   ### at least two studies left?

   if (length(yi) < 2)
      stop("Plotting terminated since k < 2.")

   ### get weights

   if (yaxis == "wi") {
      weights <- 1/vi
      weights <- weights / sum(weights) * 100
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
         ylim <- c(min(ni), max(ni))
      if (yaxis == "ninv")
         ylim <- c(max(1/ni), min(1/ni))
      if (yaxis == "sqrtni")
         ylim <- c(min(sqrt(ni)), max(sqrt(ni)))
      if (yaxis == "sqrtninv")
         ylim <- c(max(1/sqrt(ni)), min(1/sqrt(ni)))
      if (yaxis == "lni")
         ylim <- c(min(log(ni)), max(log(ni)))
      if (yaxis == "wi")
         ylim <- c(min(weights), max(weights))

      ### infinite y-axis limits can happen with "seinv" and "vinv" when one or more sampling variances are 0

      if (any(is.infinite(ylim)))
         stop("Setting 'ylim' automatically not possible (must set y-axis limits manually).")

   } else {

      ### make sure that user supplied limits are in the right order

      if (is.element(yaxis, c("sei", "vi", "ninv", "sqrtninv")))
         ylim <- c(max(ylim), min(ylim))

      if (is.element(yaxis, c("seinv", "vinv", "ni", "sqrtni", "lni", "wi")))
         ylim <- c(min(ylim), max(ylim))

      ### make sure that user supplied limits are in the appropriate range

      if (is.element(yaxis, c("sei", "vi", "ni", "ninv", "sqrtni", "sqrtninv", "lni"))) {
         if (ylim[1] < 0 || ylim[2] < 0)
            stop("Both limits for the y axis must be >= 0.")
      }

      if (is.element(yaxis, c("seinv", "vinv"))) {
         if (ylim[1] <= 0 || ylim[2] <= 0)
            stop("Both limits for the y axis must be > 0.")
      }

      if (is.element(yaxis, c("wi"))) {
         if (ylim[1] < 0 || ylim[2] < 0)
            stop("Both limits for the y axis must be >= 0.")
      }

   }

   #########################################################################

   ### set x-axis limits

   if (is.element(yaxis, c("sei", "vi", "seinv", "vinv"))) {

      level     <- ifelse(level > 1, (100-level)/100, ifelse(level > .5, 1-level, level)) ### note: there may be multiple level values
      level.min <- min(level)                                                             ### note: smallest level is the widest CI
      lvals     <- length(level)

      ### calculate the CI bounds at the bottom of the figure (for the widest CI if there are multiple)

      if (yaxis == "sei") {
         x.lb.bot <- refline - qnorm(level.min/2, lower.tail=FALSE) * sqrt(ylim[1]^2)
         x.ub.bot <- refline + qnorm(level.min/2, lower.tail=FALSE) * sqrt(ylim[1]^2)
      }
      if (yaxis == "vi") {
         x.lb.bot <- refline - qnorm(level.min/2, lower.tail=FALSE) * sqrt(ylim[1])
         x.ub.bot <- refline + qnorm(level.min/2, lower.tail=FALSE) * sqrt(ylim[1])
      }
      if (yaxis == "seinv") {
         x.lb.bot <- refline - qnorm(level.min/2, lower.tail=FALSE) * sqrt(1/ylim[1]^2)
         x.ub.bot <- refline + qnorm(level.min/2, lower.tail=FALSE) * sqrt(1/ylim[1]^2)
      }
      if (yaxis == "vinv") {
         x.lb.bot <- refline - qnorm(level.min/2, lower.tail=FALSE) * sqrt(1/ylim[1])
         x.ub.bot <- refline + qnorm(level.min/2, lower.tail=FALSE) * sqrt(1/ylim[1])
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

   plot(NA, NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, xaxt="n", yaxt="n", bty="n", ...)

   ### add background shading

   par.usr <- par("usr")
   rect(par.usr[1], par.usr[3], par.usr[2], par.usr[4], col=back, border=NA, ...)

   ### add y-axis

   axis(side=2, at=seq(from=ylim[1], to=ylim[2], length.out=steps), labels=formatC(seq(from=ylim[1], to=ylim[2], length.out=steps), digits=digits[[2]], format="f", drop0trailing=ifelse(class(digits[[2]]) == "integer", TRUE, FALSE)), ...)

   ### add horizontal lines

   abline(h=seq(from=ylim[1], to=ylim[2], length.out=steps), col=hlines, ...)

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

         ci.left  <- refline - qnorm(level[m]/2, lower.tail=FALSE) * sqrt(vi.vals)
         ci.right <- refline + qnorm(level[m]/2, lower.tail=FALSE) * sqrt(vi.vals)

         polygon(c(ci.left,ci.right[ci.res:1]), c(yi.vals,yi.vals[ci.res:1]), border=NA, col=shade[m], ...)
         lines(ci.left,  yi.vals, lty="dotted", ...)
         lines(ci.right, yi.vals, lty="dotted", ...)

      }

   }

   ### add vertical reference line
   ### use segments so that line does not extent beyond tip of CI region

   if (is.element(yaxis, c("sei", "vi", "seinv", "vinv")))
      segments(refline, ylim[1], refline, ylim[2], ...)

   if (is.element(yaxis, c("ni", "ninv", "sqrtni", "sqrtninv", "lni", "wi")))
      abline(v=refline, ...)

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

   points(xaxis.vals, yaxis.vals, pch=pch, ...)

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

   axis(side=1, at=at, labels=at.lab, ...)

   ### prepare data frame to return

   sav <- data.frame(x=xaxis.vals, y=yaxis.vals, slab=slab)

   invisible(sav)

}
