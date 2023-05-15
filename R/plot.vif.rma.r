plot.vif.rma <- function(x,
   breaks="Scott", freq=FALSE, col, border, col.out, col.density,
   trim=0, adjust=1, lwd=c(2,0), layout, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="vif.rma")

   if (missing(col)) {
      if (.is.dark(par("bg"))) {
         col <- "gray50"
      } else {
         col <- "gray"
      }
   }

   if (missing(border)) {
      if (.is.dark(par("bg"))) {
         border <- par("bg")
      } else {
         border <- "white"
      }
   }

   if (missing(col.out)) {
      if (.is.dark(par("bg"))) {
         col.out <- rgb(1,0,0,0.4)
      } else {
         col.out <- rgb(1,0,0,0.5)
      }
   }

   if (missing(col.density)) {
      if (.is.dark(par("bg"))) {
         col.density <- "dodgerblue"
      } else {
         col.density <- "blue"
      }
   }

   par.mfrow <- par("mfrow")

   if (!is.null(x$alpha)) {

      if (is.null(x[[2]]$sim)) {
         plot(x[[1]], breaks=breaks, freq=freq, col=col, border=border, trim=trim, col.out=col.out,
              col.density=col.density, adjust=adjust, lwd=lwd, layout=layout, mainadd="Location ", ...)
         return(invisible())
      }

      if (is.null(x[[1]]$sim)) {
         plot(x[[2]], breaks=breaks, freq=freq, col=col, border=border, trim=trim, col.out=col.out,
              col.density=col.density, adjust=adjust, lwd=lwd, layout=layout, mainadd="Scale ", ...)
         return(invisible())
      }

      np <- length(x[[1]]$vifs) + length(x[[2]]$vifs)

      ### set/check layout argument

      if (missing(layout)) {
         layout <- n2mfrow(np)
      } else {
         layout <- layout[layout >= 1]
         layout <- round(layout)
         if (length(layout) != 2L)
            stop(mstyle$stop("Incorrect specification of 'layout' argument."))
      }

      plot(x[[1]], breaks=breaks, freq=freq, col=col, border=border, trim=trim, col.out=col.out,
           col.density=col.density, adjust=adjust, lwd=lwd, layout=layout, mainadd="Location ", ...)
      plot(x[[2]], breaks=breaks, freq=freq, col=col, border=border, trim=trim, col.out=col.out,
           col.density=col.density, adjust=adjust, lwd=lwd, mainadd="Scale ", new=FALSE, par.mfrow=par.mfrow, ...)

      return(invisible())

   }

   ddd <- list(...)

   if (is.null(ddd$tail)) {
      tail <- "upper"
   } else {
      tail <- match.arg(ddd$tail, c("lower", "upper"))
   }

   if (is.null(ddd$new)) {
      new <- TRUE
   } else {
      new <- FALSE
   }

   if (is.null(ddd$mainadd)) {
      mainadd <- ""
   } else {
      mainadd <- ddd$mainadd
   }

   ### check if 'sim' was actually used

   if (is.null(x$sim))
      stop(mstyle$stop("Can only plot 'vif.rma' objects when 'sim=TRUE' was used."))

   ### number of plots

   np <- length(x$vifs)

   ### set/check layout argument

   if (missing(layout)) {
      layout <- n2mfrow(np)
   } else {
      layout <- layout[layout >= 1]
      layout <- round(layout)
      if (length(layout) != 2L)
         stop(mstyle$stop("Incorrect specification of 'layout' argument."))
   }

   ### 1st: obs stat, 2nd: density

   if (length(lwd) == 1L)
      lwd <- lwd[c(1,1)]

   ### cannot plot density when freq=TRUE

   if (freq)
      lwd[2] <- 0

   ### check trim

   if (trim >= 0.5)
      stop(mstyle$stop("The value of 'trim' must be < 0.5."))

   ### local plotting functions

   lhist     <- function(..., tail, new, par.mfrow, mainadd) hist(...)
   labline   <- function(..., tail, new, par.mfrow, mainadd) abline(...)
   lsegments <- function(..., tail, new, par.mfrow, mainadd) segments(...)
   llines    <- function(..., tail, new, par.mfrow, mainadd) lines(...)

   ############################################################################

   if (new) {
      par(mfrow=layout)
   } else {
      on.exit(par(mfrow = ddd$par.mfrow), add=TRUE)
   }

   for (i in seq_len(np)) {

      pvif <- x$sim[,i]
      pvif <- pvif[is.finite(pvif)]

      den <- density(pvif, adjust=adjust)

      if (trim > 0) {
         bound <- quantile(pvif, probs=1-trim)
         pvif <- pvif[pvif <= bound]
      }

      tmp <- lhist(pvif, breaks=breaks, plot=FALSE)

      ylim <- c(0, max(ifelse(lwd[2] == 0, 0, max(den$y)), max(tmp$density)))

      tmp <- lhist(pvif, breaks=breaks, col=col, border=border,
                   main=paste0(mainadd, "Coefficient", ifelse(x$vif[[i]]$m > 1, "s", ""), ": ", names(x$vifs)[i]),
                   xlab="Value of VIF",
                   freq=freq, ylim=ylim, xaxt="n", ...)

      xat <- axTicks(side=1)
      xlabels <- xat

      axis(side=1, at=xat, labels=xlabels)

      .coltail(tmp, val=x$vifs[i], col=col.out, border=border, freq=freq, ...)

      usr <- par()$usr

      if (x$vifs[i] > usr[2] && lwd[1] > 0) {
         ya <- mean(par()$yaxp[1:2])
         arrows(usr[2] - .08*(usr[2]-usr[1]), ya, usr[2] - .01*(usr[2]-usr[1]), ya,
                length = .02*(grconvertY(usr[4], from="user", to="inches")-
                             (grconvertY(usr[3], from="user", to="inches"))))
      }

      x$vifs[i] <- min(x$vifs[i], usr[2])

      par(xpd = TRUE)
      lsegments(x$vifs[i], usr[3], x$vifs[i], usr[4], lwd=lwd[1], lty="dashed", ...)
      par(xpd = FALSE)

      #den$y <- den$y[den$x <= par()$xaxp[2]]
      #den$x <- den$x[den$x <= par()$xaxp[2]]
      llines(den, lwd=lwd[2], col=col.density, ...)

   }

   ############################################################################

   invisible()

}
