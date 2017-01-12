radial.rma <- function(x, center=FALSE, xlim=NULL, zlim,
xlab, zlab, atz, aty, steps=7, level=x$level, digits=2,
back="lightgray", transf, targs, pch=19, arc.res=100, cex, ...) {

   #########################################################################

   if (!inherits(x, "rma"))
      stop("Argument 'x' must be an object of class \"rma\".")

   if (inherits(x, "robust.rma"))
      stop("Function not applicable to objects of class \"robust.rma\".")

   if (missing(transf))
      transf <- FALSE

   if (missing(targs))
      targs <- NULL

   if (missing(atz))
      atz <- NULL

   if (missing(aty))
      aty <- NULL

   #########################################################################

   ### radial plots only for intercept-only models

   if (x$int.only) {

      yi   <- x$yi
      yi.c <- yi
      vi   <- x$vi
      b    <- c(x$b)
      ci.lb <- x$ci.lb
      ci.ub <- x$ci.ub
      tau2 <- 1/mean(1/x$tau2) ### geometric mean of tau^2 values (hackish solution for models with multiple tau^2 values)
                               ### note: this works for 1/mean(1/0) = 0; TODO: consider something more sophisticated here
      if (is.null(aty)) {
         atyis <- range(yi)
      } else {
         atyis <- range(aty)
         aty.c <- aty
      }

   } else {
      stop("Radial plots only applicable for models without moderators.")
   }

   if (center) {
      yi    <- yi - x$b
      b     <- 0
      ci.lb <- ci.lb - x$b
      ci.ub <- ci.ub - x$b
      atyis <- atyis - x$b
      if (!is.null(aty))
         aty <- aty - x$b
   }

   #########################################################################

   level <- ifelse(level > 1, (100-level)/100, ifelse(level > .5, 1-level, level))
   zcrit <- qnorm(level/2, lower.tail=FALSE)

   zi <- yi / sqrt(vi+tau2)
   xi <-  1 / sqrt(vi+tau2)

   ### if vi=0 and tau2=0, then zi and xi will be Inf

   if (any(is.infinite(c(xi,zi))))
      stop("Setting 'xlim' and 'zlim' automatically not possible (must set axis limits manually).")

   ### set x axis limits if none are specified

   if (missing(xlim)) {
      xlims <- c(0, (1.30*max(xi)))                   ### add 30% to upper bound
   } else {
      xlims <- sort(xlim)
   }

   ### x axis position of the confidence interval

   ci.xpos <- xlims[2] + 0.12*(xlims[2]-xlims[1])     ### add 12% of range to upper bound

   ### x axis position of the y axis on the right

   ya.xpos <- xlims[2] + 0.14*(xlims[2]-xlims[1])     ### add 14% of range to upper bound

   xaxismax <- xlims[2]

   ### set z axis limits if none are specified (these are the actual y axis limits of the plot)

   if (missing(zlim)) {
      zlims <- c(min(-5, 1.10*min(zi), 1.10*ci.lb*ci.xpos, 1.10*min(atyis)*ya.xpos, 1.10*min(yi)*ya.xpos, -1.10*zcrit+xaxismax*b), max(5, 1.10*max(zi), 1.10*ci.ub*ci.xpos, 1.10*max(atyis)*ya.xpos, 1.10*max(yi)*ya.xpos, 1.10*zcrit+xaxismax*b))
   } else {
      zlims <- sort(zlim)
   }

   ### adjust margins

   par.mar <- par("mar")
   par.mar.adj <- par.mar - c(0,-3,0,-5)
   par.mar.adj[par.mar.adj < 1] <- 1
   par(mar = par.mar.adj)
   on.exit(par(mar = par.mar))

   ### label for the x axis

   if (missing(xlab)) {
      if (x$method == "FE") {
         xlab <- expression(x[i]==1/sqrt(v[i]), ...)
      } else {
         xlab <- expression(x[i]==1/sqrt(v[i]+tau^2), ...)
      }
   }

   par.pty <- par("pty")
   par(pty="s")
   on.exit(par(pty = par.pty), add=TRUE)

   plot(NA, NA, ylim=zlims, xlim=xlims, bty="n", xaxt="n", yaxt="n", xlab=xlab, ylab="", xaxs="i", yaxs="i", ...)

   if (missing(cex))
      cex <- par("cex")

   ### add polygon and +-zcrit lines

   polygon(c(0,xaxismax,xaxismax,0), c(zcrit, zcrit+xaxismax*b, -zcrit+xaxismax*b, -zcrit), border=NA, col=back, ...)
   segments(0, 0, xaxismax, xaxismax*b, lty="solid", ...)
   segments(0, -zcrit, xaxismax, -zcrit+xaxismax*b, lty="dotted", ...)
   segments(0,  zcrit, xaxismax,  zcrit+xaxismax*b, lty="dotted", ...)

   ### add x axis

   axis(side=1, ...)

   ### add z axis

   if (is.null(atz)) {
      axis(side=2, at=seq(-4,4, length=9), labels=NA, las=1, tcl=par("tcl")/2, ...)
      axis(side=2, at=seq(-2,2, length=3), las=1, ...)
   } else {
      axis(side=2, at=atz, labels=atz, las=1, ...)
   }

   ### add label for the z axis

   if (missing(zlab)) {
      if (center) {
         if (x$method == "FE") {
            mtext(expression(z[i]==frac(y[i]-hat(theta),sqrt(v[i]))), side=2, line=par.mar.adj[2]-1, at=0, adj=0, las=1, cex=cex, ...)
         } else {
            mtext(expression(z[i]==frac(y[i]-hat(mu),sqrt(v[i]+tau^2))), side=2, line=par.mar.adj[2]-1, adj=0, at=0, las=1, cex=cex, ...)
         }
      } else {
         if (x$method == "FE") {
            mtext(expression(z[i]==frac(y[i],sqrt(v[i]))), side=2, line=par.mar.adj[2]-2, at=0, adj=0, las=1, cex=cex, ...)
         } else {
            mtext(expression(z[i]==frac(y[i],sqrt(v[i]+tau^2))), side=2, line=par.mar.adj[2]-1, at=0, adj=0, las=1, cex=cex, ...)
         }

      }
   } else {
      mtext(zlab, side=2, line=par.mar.adj[2]-4, at=0, cex=cex, ...)
   }


   #########################################################################

   ### add y axis arc and CI arc on the right

   par.xpd <- par("xpd")
   par(xpd=TRUE)

   par.usr <- par("usr")
   asp.rat <- (par.usr[4]-par.usr[3])/(par.usr[2]-par.usr[1])

   if (length(arc.res) == 1L)
      arc.res <- c(arc.res, arc.res/4)

   ### add y axis arc

   if (is.null(aty)) {
      atyis <- seq(min(yi), max(yi), length=arc.res[1])
   } else {
      atyis <- seq(min(aty), max(aty), length=arc.res[1])
   }

   len <- ya.xpos
   xis <- rep(NA_real_,length(atyis))
   zis <- rep(NA_real_,length(atyis))
   for (i in seq_len(length(atyis))) {
      xis[i] <- sqrt(len^2/(1+(atyis[i]/asp.rat)^2))
      zis[i] <- xis[i]*atyis[i]
   }

   valid <- zis > zlims[1] & zis < zlims[2]
   lines(xis[valid], zis[valid], ...)

   ### add y axis tick marks

   if (is.null(aty)) {
      atyis <- seq(min(yi), max(yi), length=steps)
   } else {
      atyis <- aty
   }

   len.l <- ya.xpos
   len.u <- ya.xpos + .015*(xlims[2]-xlims[1])
   xis.l <- rep(NA_real_,length(atyis))
   zis.l <- rep(NA_real_,length(atyis))
   xis.u <- rep(NA_real_,length(atyis))
   zis.u <- rep(NA_real_,length(atyis))
   for (i in seq_len(length(atyis))) {
      xis.l[i] <- sqrt(len.l^2/(1+(atyis[i]/asp.rat)^2))
      zis.l[i] <- xis.l[i]*atyis[i]
      xis.u[i] <- sqrt(len.u^2/(1+(atyis[i]/asp.rat)^2))
      zis.u[i] <- xis.u[i]*atyis[i]
   }

   valid <- zis.l > zlims[1] & zis.u > zlims[1] & zis.l < zlims[2] & zis.u < zlims[2]

   if (any(valid))
      segments(xis.l[valid], zis.l[valid], xis.u[valid], (xis.u*atyis)[valid], ...)

   ### add y axis labels

   if (is.null(aty)) {
      atyis     <- seq(min(yi),   max(yi),   length=steps)
      atyis.lab <- seq(min(yi.c), max(yi.c), length=steps)
   } else {
      atyis     <- aty
      atyis.lab <- aty.c
   }

   len <- ya.xpos+.02*(xlims[2]-xlims[1])
   xis <- rep(NA_real_,length(atyis))
   zis <- rep(NA_real_,length(atyis))
   for (i in seq_len(length(atyis))) {
      xis[i] <- sqrt(len^2/(1+(atyis[i]/asp.rat)^2))
      zis[i] <- xis[i]*atyis[i]
   }

   if (is.function(transf)) {
      if (is.null(targs)) {
         atyis.lab <- sapply(atyis.lab, transf)
      } else {
         atyis.lab <- sapply(atyis.lab, transf, targs)
      }
   }

   valid <- zis > zlims[1] & zis < zlims[2]

   if (any(valid))
      text(xis[valid], zis[valid], formatC(atyis.lab[valid], digits=digits, format="f"), pos=4, cex=cex, ...)

   ### add CI arc

   atyis <- seq(ci.lb, ci.ub, length=arc.res[2])
   len <- ci.xpos
   xis <- rep(NA_real_,length(atyis))
   zis <- rep(NA_real_,length(atyis))
   for (i in seq_len(length(atyis))) {
      xis[i] <- sqrt(len^2/(1+(atyis[i]/asp.rat)^2))
      zis[i] <- xis[i]*atyis[i]
   }

   valid <- zis > zlims[1] & zis < zlims[2]

   if (any(valid))
      lines(xis[valid], zis[valid], ...)

   ### add CI tick marks

   atyis <- c(ci.lb, b, ci.ub)
   len.l <- ci.xpos-.007*(xlims[2]-xlims[1])
   len.u <- ci.xpos+.007*(xlims[2]-xlims[1])
   xis.l <- rep(NA_real_,3)
   zis.l <- rep(NA_real_,3)
   xis.u <- rep(NA_real_,3)
   zis.u <- rep(NA_real_,3)
   for (i in seq_len(length(atyis))) {
      xis.l[i] <- sqrt(len.l^2/(1+(atyis[i]/asp.rat)^2))
      zis.l[i] <- xis.l[i]*atyis[i]
      xis.u[i] <- sqrt(len.u^2/(1+(atyis[i]/asp.rat)^2))
      zis.u[i] <- xis.u[i]*atyis[i]
   }

   valid <- zis.l > zlims[1] & zis.u > zlims[1] & zis.l < zlims[2] & zis.u < zlims[2]

   if (any(valid))
      segments(xis.l[valid], zis.l[valid], xis.u[valid], (xis.u*atyis)[valid], ...)

   par(xpd=par.xpd)

   #########################################################################

   ### add points to the plot

   points(xi, zi, pch=pch, cex=cex, ...)

   invisible(data.frame(x=xi, y=zi, slab=x$slab[x$not.na]))

}
