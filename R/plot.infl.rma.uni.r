plot.infl.rma.uni <- function(x, plotinf=TRUE, plotdfbs=FALSE, dfbsnew=FALSE, logcov=TRUE,
slab.style=1, las=0, pch=21, bg, bg.infl, col.na, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(x), must="infl.rma.uni")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   ddd <- list(...)

   if (!is.null(ddd$layout))
      warning(mstyle$warning("Argument 'layout' has been deprecated."), call.=FALSE)

   .start.plot()

   if (missing(bg))
      bg <- .coladj(par("bg","fg"), dark=0.35, light=-0.35)

   if (missing(bg.infl))
      bg.infl <- "red"

   if (missing(col.na))
      col.na <- .coladj(par("bg","fg"), dark=0.2, light=-0.2)

   #########################################################################

   ### check for NAs and stop if there are any when na.act == "na.fail"

   any.na <- is.na(as.data.frame(x$inf))

   if (any(any.na) && na.act == "na.fail")
      stop(mstyle$stop("Missing values in results."))

   #########################################################################

   ### process plotinf argument

   if (is.logical(plotinf)) {
      if (plotinf)
         which.inf <- seq_len(8)
   } else {
      which.inf <- plotinf
      which.inf <- which.inf[(which.inf >= 1) & (which.inf <= 8)]
      which.inf <- unique(round(which.inf))
      if (length(which.inf) == 0L)
         stop(mstyle$stop("Incorrect specification of the 'plotinf' argument."))
      plotinf <- TRUE
   }

   ### process plotdfbs argument

   if (is.logical(plotdfbs)) {
      if (plotdfbs)
         which.dfbs <- seq_len(x$p)
   } else {
      which.dfbs <- plotdfbs
      which.dfbs <- which.dfbs[(which.dfbs >= 1) & (which.dfbs <= x$p)]
      which.dfbs <- unique(round(which.dfbs))
      if (length(which.dfbs) == 0L)
         stop(mstyle$stop("Incorrect specification of the 'plotdfbs' argument."))
      plotdfbs <- TRUE
   }

   #########################################################################

   if (!plotinf & !plotdfbs)
      stop(mstyle$stop("At least one of the arguments 'plotinf' or 'plotdfbs' must be TRUE."))

   if (!plotinf & dfbsnew)
      dfbsnew <- FALSE

   par.mar <- par("mar")
   par.mar.adj <- par.mar - c(2,1,2,0)
   par.mar.adj[par.mar.adj < 1] <- 1
   par(mar=par.mar.adj)
   on.exit(par(mar=par.mar), add=TRUE)

   #########################################################################

   ### filter out potential arguments to abbreviate() (which cause problems with the various plot functions)

   lplot   <- function(..., minlength, strict, layout) plot(...)
   lpoints <- function(..., minlength, strict, layout) points(...)
   llines  <- function(..., minlength, strict, layout) lines(...)
   laxis   <- function(..., minlength, strict, layout) axis(...)
   labline <- function(..., minlength, strict, layout) abline(...)

   #########################################################################

   ids <- switch(slab.style, "1" = x$ids, "2" = x$inf$slab, "3" = abbreviate(x$inf$slab, ...))
   #print(ids)

   #########################################################################

   ### plot inf values if requested

   if (plotinf) {

      np.inf <- length(which.inf)

      if (np.inf > 1L) {

         # if no plotting device is open or mfrow is too small, set mfrow appropriately

         if (dev.cur() == 1L || prod(par("mfrow")) < np.inf) {

            #par(mfrow=n2mfrow(np.inf)) # this behaves slightly differently (see below)

            if (np.inf == 2L)
               par(mfrow=c(2,1))
            if (np.inf == 3L)
               par(mfrow=c(3,1))
            if (np.inf == 4L)
               par(mfrow=c(2,2))
            if (np.inf == 5L)
               par(mfrow=c(5,1)) # n2mfrow(5) yields c(3,2)
            if (np.inf == 6L)
               par(mfrow=c(3,2))
            if (np.inf == 7L)
               par(mfrow=c(7,1)) # n2mfrow(7) yields c(3,3)
            if (np.inf == 8L)
               par(mfrow=c(4,2)) # n2mfrow(8) yields c(3,3)

         }

         on.exit(par(mfrow=c(1L,1L)), add=TRUE)

      }

      ######################################################################

      for (i in seq_along(which.inf)) {

         if (which.inf[i] == 1) {

            zi     <- x$inf$rstudent
            not.na <- !is.na(zi)

            if (na.act == "na.omit") {
               zi       <- zi[not.na]
               len.ids  <- length(x$ids)-sum(!not.na)
               ids.infl <- x$is.infl[not.na]
               lab.ids  <- ids[not.na]
            }
            if (na.act == "na.exclude" || na.act == "na.pass") {
               len.ids  <- length(x$ids)
               ids.infl <- x$is.infl
               lab.ids  <- ids
            }
            if (any(!is.na(zi))) {
               zi.min <- min(min(zi,-2,na.rm=TRUE), qnorm(0.025))*1.05
               zi.max <- max(max(zi, 2,na.rm=TRUE), qnorm(0.975))*1.05
               lplot(NA, NA, xlim=c(1,len.ids), ylim=c(zi.min,zi.max), xaxt="n", main="rstudent", xlab="", ylab="", las=las, ...)
               laxis(side=1, at=seq_len(len.ids), labels=lab.ids, xlab="", las=las, ...)
               labline(h=0, lty="dashed", ...)
               labline(h=c(qnorm(0.025),qnorm(0.975)), lty="dotted", ...)
               if (na.act == "na.exclude" || na.act == "na.pass")
                  llines(seq_len(len.ids)[not.na], zi[not.na], col=col.na, ...)
               llines(seq_len(len.ids), zi, ...)
               lpoints(x=seq_len(len.ids),           y=zi,           bg=bg,      pch=pch, ...)
               lpoints(x=seq_len(len.ids)[ids.infl], y=zi[ids.infl], bg=bg.infl, pch=pch, ...)
               #if (num.infl)
               #   text(seq_len(len.ids)[ids.infl], zi[ids.infl], seq_len(len.ids)[ids.infl], pos=ifelse(zi[ids.infl] > 0, 3, 1), ...)
            } else {
               lplot(NA, NA, xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n", main="rstudent", xlab="", ylab="", ...)
            }

         }

         if (which.inf[i] == 2) {

            zi     <- x$inf$dffits
            not.na <- !is.na(zi)

            if (na.act == "na.omit") {
               zi       <- zi[not.na]
               len.ids  <- length(x$ids)-sum(!not.na)
               ids.infl <- x$is.infl[not.na]
               lab.ids  <- ids[not.na]
            }
            if (na.act == "na.exclude" || na.act == "na.pass") {
               len.ids  <- length(x$ids)
               ids.infl <- x$is.infl
               lab.ids  <- ids
            }
            if (any(!is.na(zi))) {
               zi.min <- min(min(zi,na.rm=TRUE), -3*sqrt(x$p/(x$k-x$p)))*1.05
               zi.max <- max(max(zi,na.rm=TRUE),  3*sqrt(x$p/(x$k-x$p)))*1.05
               lplot(NA, NA, xlim=c(1,len.ids), ylim=c(zi.min,zi.max), xaxt="n", main="dffits", xlab="", ylab="", las=las, ...)
               laxis(side=1, at=seq_len(len.ids), labels=lab.ids, xlab="", las=las, ...)
               labline(h= 0, lty="dashed", ...)
               labline(h= 3*sqrt(x$p/(x$k-x$p)), lty="dotted", ...)
               labline(h=-3*sqrt(x$p/(x$k-x$p)), lty="dotted", ...)
               if (na.act == "na.exclude" || na.act == "na.pass")
                  llines(seq_len(len.ids)[not.na], zi[not.na], col=col.na, ...)
               llines(seq_len(len.ids), zi, ...)
               lpoints(x=seq_len(len.ids),           y=zi,           bg=bg,      pch=pch, ...)
               lpoints(x=seq_len(len.ids)[ids.infl], y=zi[ids.infl], bg=bg.infl, pch=pch, ...)
               #if (num.infl)
               #   text(seq_len(len.ids)[ids.infl], zi[ids.infl], seq_len(len.ids)[ids.infl], pos=ifelse(zi[ids.infl] > 0, 3, 1), ...)
            } else {
               lplot(NA, NA, xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n", main="dffits", xlab="", ylab="", ...)
            }

         }

         if (which.inf[i] == 3) {

            zi     <- x$inf$cook.d
            not.na <- !is.na(zi)

            if (na.act == "na.omit") {
               zi       <- zi[not.na]
               len.ids  <- length(x$ids)-sum(!not.na)
               ids.infl <- x$is.infl[not.na]
               lab.ids  <- ids[not.na]
            }
            if (na.act == "na.exclude" || na.act == "na.pass") {
               len.ids  <- length(x$ids)
               ids.infl <- x$is.infl
               lab.ids  <- ids
            }
            if (any(!is.na(zi))) {
               zi.min <- 0
               zi.max <- max(zi,na.rm=TRUE)*1.05
               lplot(NA, NA, xlim=c(1,len.ids), ylim=c(zi.min,zi.max), xaxt="n", main="cook.d", xlab="", ylab="", las=las, ...)
               laxis(side=1, at=seq_len(len.ids), labels=lab.ids, xlab="", las=las, ...)
               labline(h=qchisq(0.5, df=x$m), lty="dotted", ...)
               if (na.act == "na.exclude" || na.act == "na.pass")
                  llines(seq_len(len.ids)[not.na], zi[not.na], col=col.na, ...)
               llines(seq_len(len.ids), zi, ...)
               lpoints(x=seq_len(len.ids),           y=zi,           bg=bg,      pch=pch, ...)
               lpoints(x=seq_len(len.ids)[ids.infl], y=zi[ids.infl], bg=bg.infl, pch=pch, ...)
               #if (num.infl)
               #   text(seq_len(len.ids)[ids.infl], zi[ids.infl], seq_len(len.ids)[ids.infl], pos=3, ...)
            } else {
               lplot(NA, NA, xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n", main="cook.d", xlab="", ylab="", ...)
            }

         }

         if (which.inf[i] == 4) {

            zi     <- x$inf$cov.r
            not.na <- !is.na(zi)

            if (na.act == "na.omit") {
               zi       <- zi[not.na]
               len.ids  <- length(x$ids)-sum(!not.na)
               ids.infl <- x$is.infl[not.na]
               lab.ids  <- ids[not.na]
            }
            if (na.act == "na.exclude" || na.act == "na.pass") {
               len.ids  <- length(x$ids)
               ids.infl <- x$is.infl
               lab.ids  <- ids
            }
            if (any(!is.na(zi))) {
               zi.min <- min(zi,na.rm=TRUE)
               zi.max <- max(zi,na.rm=TRUE)
               if (logcov) {
                  lplot(NA, NA, xlim=c(1,len.ids), ylim=c(zi.min,zi.max), xaxt="n", main="cov.r", xlab="", ylab="", las=las, log="y", ...)
               } else {
                  lplot(NA, NA, xlim=c(1,len.ids), ylim=c(zi.min,zi.max), xaxt="n", main="cov.r", xlab="", ylab="", las=las, ...)
               }
               laxis(side=1, at=seq_len(len.ids), labels=lab.ids, xlab="", las=las, ...)
               labline(h=1, lty="dashed", ...)
               #labline(h=1+3*x$m/(x$k-x$m), lty="dotted", ...)
               #labline(h=1-3*x$m/(x$k-x$m), lty="dotted", ...)
               if (na.act == "na.exclude" || na.act == "na.pass")
                  llines(seq_len(len.ids)[not.na], zi[not.na], col=col.na, ...)
               llines(seq_len(len.ids), zi, ...)
               lpoints(x=seq_len(len.ids),           y=zi,           bg=bg,      pch=pch, ...)
               lpoints(x=seq_len(len.ids)[ids.infl], y=zi[ids.infl], bg=bg.infl, pch=pch, ...)
               #if (num.infl)
               #   text(seq_len(len.ids)[ids.infl], zi[ids.infl], seq_len(len.ids)[ids.infl], pos=ifelse(zi[ids.infl] > 1, 3, 1), ...)
            } else {
               lplot(NA, NA, xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n", main="cov.r", xlab="", ylab="", ...)
            }

         }

         if (which.inf[i] == 5) {

            zi     <- x$inf$tau2.del
            not.na <- !is.na(zi)

            if (na.act == "na.omit") {
               zi       <- zi[not.na]
               len.ids  <- length(x$ids)-sum(!not.na)
               ids.infl <- x$is.infl[not.na]
               lab.ids  <- ids[not.na]
            }
            if (na.act == "na.exclude" || na.act == "na.pass") {
               len.ids  <- length(x$ids)
               ids.infl <- x$is.infl
               lab.ids  <- ids
            }
            if (any(!is.na(zi))) {
               zi.min <- min(zi,na.rm=TRUE)
               zi.max <- max(zi,na.rm=TRUE)
               lplot(NA, NA, xlim=c(1,len.ids), ylim=c(zi.min,zi.max), xaxt="n", main="tau2.del", xlab="", ylab="", las=las, ...)
               laxis(side=1, at=seq_len(len.ids), labels=lab.ids, xlab="", las=las, ...)
               labline(h=x$tau2, lty="dashed", ...)
               if (na.act == "na.exclude" || na.act == "na.pass")
                  llines(seq_len(len.ids)[not.na], zi[not.na], col=col.na, ...)
               llines(seq_len(len.ids), zi, ...)
               lpoints(x=seq_len(len.ids),           y=zi,           bg=bg,      pch=pch, ...)
               lpoints(x=seq_len(len.ids)[ids.infl], y=zi[ids.infl], bg=bg.infl, pch=pch, ...)
            } else {
               lplot(NA, NA, xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n", main="tau2.del", xlab="", ylab="", ...)
            }

         }

         if (which.inf[i] == 6) {

            zi     <- x$inf$QE.del
            not.na <- !is.na(zi)

            if (na.act == "na.omit") {
               zi       <- zi[not.na]
               len.ids  <- length(x$ids)-sum(!not.na)
               ids.infl <- x$is.infl[not.na]
               lab.ids  <- ids[not.na]
            }
            if (na.act == "na.exclude" || na.act == "na.pass") {
               len.ids  <- length(x$ids)
               ids.infl <- x$is.infl
               lab.ids  <- ids
            }
            if (any(!is.na(zi))) {
               zi.min <- min(zi,na.rm=TRUE)
               zi.max <- max(zi,na.rm=TRUE)
               lplot(NA, NA, xlim=c(1,len.ids), ylim=c(zi.min,zi.max), xaxt="n", main="QE.del", xlab="", ylab="", las=las, ...)
               laxis(side=1, at=seq_len(len.ids), labels=lab.ids, xlab="", las=las, ...)
               labline(h=x$QE, lty="dashed", ...)
               #labline(h=qchisq(0.95, df=x$k-x$p), lty="dotted", ...)
               labline(h=x$k-x$p, lty="dotted", ...)
               if (na.act == "na.exclude" || na.act == "na.pass")
                  llines(seq_len(len.ids)[not.na], zi[not.na], col=col.na, ...)
               llines(seq_len(len.ids), zi, ...)
               lpoints(x=seq_len(len.ids),           y=zi,           bg=bg,      pch=pch, ...)
               lpoints(x=seq_len(len.ids)[ids.infl], y=zi[ids.infl], bg=bg.infl, pch=pch, ...)
            } else {
               lplot(NA, NA, xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n", main="QE.del", xlab="", ylab="", ...)
            }

         }

         if (which.inf[i] == 7) {

            zi     <- x$inf$hat
            not.na <- !is.na(zi)

            if (na.act == "na.omit") {
               zi       <- zi[not.na]
               len.ids  <- length(x$ids)-sum(!not.na)
               ids.infl <- x$is.infl[not.na]
               lab.ids  <- ids[not.na]
            }
            if (na.act == "na.exclude" || na.act == "na.pass") {
               len.ids  <- length(x$ids)
               ids.infl <- x$is.infl
               lab.ids  <- ids
            }
            if (any(!is.na(zi))) {
               zi.min <- 0
               zi.max <- max(max(zi,na.rm=TRUE), 3*x$p/x$k)*1.05
               lplot(NA, NA, xlim=c(1,len.ids), ylim=c(zi.min,zi.max), xaxt="n", main="hat", xlab="", ylab="", las=las, ...)
               laxis(side=1, at=seq_len(len.ids), labels=lab.ids, xlab="", las=las, ...)
               labline(h=x$p/x$k, lty="dashed", ...)
               labline(h=3*x$p/x$k, lty="dotted", ...)
               if (na.act == "na.exclude" || na.act == "na.pass")
                  llines(seq_len(len.ids)[not.na], zi[not.na], col=col.na, ...)
               llines(seq_len(len.ids), zi, ...)
               lpoints(x=seq_len(len.ids),           y=zi,           bg=bg,      pch=pch, ...)
               lpoints(x=seq_len(len.ids)[ids.infl], y=zi[ids.infl], bg=bg.infl, pch=pch, ...)
            } else {
               lplot(NA, NA, xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n", main="hat", xlab="", ylab="", ...)
            }

         }

         if (which.inf[i] == 8) {

            zi     <- x$inf$weight
            not.na <- !is.na(zi)

            if (na.act == "na.omit") {
               zi       <- zi[not.na]
               len.ids  <- length(x$ids)-sum(!not.na)
               ids.infl <- x$is.infl[not.na]
               lab.ids  <- ids[not.na]
            }
            if (na.act == "na.exclude" || na.act == "na.pass") {
               len.ids  <- length(x$ids)
               ids.infl <- x$is.infl
               lab.ids  <- ids
            }
            if (any(!is.na(zi))) {
               zi.min <- 0
               zi.max <- max(zi,na.rm=TRUE)*1.05
               lplot(NA, NA, xlim=c(1,len.ids), ylim=c(zi.min,zi.max), xaxt="n", main="weight", xlab="", ylab="", las=las, ...)
               laxis(side=1, at=seq_len(len.ids), labels=lab.ids, xlab="", las=las, ...)
               labline(h=100/x$k, lty="dashed", ...)
               if (na.act == "na.exclude" || na.act == "na.pass")
                  llines(seq_len(len.ids)[not.na], zi[not.na], col=col.na, ...)
               llines(seq_len(len.ids), zi, ...)
               lpoints(x=seq_len(len.ids),           y=zi,           bg=bg,      pch=pch, ...)
               lpoints(x=seq_len(len.ids)[ids.infl], y=zi[ids.infl], bg=bg.infl, pch=pch, ...)
            } else {
               lplot(NA, NA, xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n", main="weight", xlab="", ylab="", ...)
            }

         }

      }

   }

   #########################################################################

   ### plot dfbs values if requested

   if (plotdfbs) {

      np.dfbs <- length(which.dfbs)

      if (plotinf && (np.inf + np.dfbs <= prod(par("mfrow")))) {

         # if np.inf + np.dfbs is small enough to fit on the same multi-panel
         # plot, then do so, but reset mfrow to c(1L,1L) for consistency

         on.exit(par(mfrow=c(1L,1L)), add=TRUE)

      } else {

         if (dfbsnew) { # this is always FALSE when plotinf=FALSE
            dev.new()
            .start.plot()
            par(mar=par.mar.adj)
         } else {
            if (plotinf) {
               caps <- dev.capabilities()$events
               if (any(is.element(c("MouseDown","Keybd"), caps))) {
                  message(mstyle$message("Press any key or click on the plot to show the DFBETAS values ...."), appendLF=FALSE)
                  getGraphicsEvent(prompt="", onMouseDown=function(button,x,y) return(1), onKeybd=function(key) return(1))
               } else {
                  par.ask <- par("ask")
                  par(ask=TRUE)
                  on.exit(par(ask=par.ask), add=TRUE)
               }
            }
         }

         # if no plotting device is open or mfrow is too small, set mfrow appropriately
         if (plotinf || dev.cur() == 1L || prod(par("mfrow")) < np.dfbs)
            par(mfrow=n2mfrow(np.dfbs))
         on.exit(par(mfrow=c(1L,1L)), add=TRUE)

      }

      for (i in seq_along(which.dfbs)) {

         zi     <- x$dfbs[[which.dfbs[i]]]
         not.na <- !is.na(zi)

         if (na.act == "na.omit") {
            zi       <- zi[not.na]
            len.ids  <- length(x$ids)-sum(!not.na)
            ids.infl <- x$is.infl[not.na]
            lab.ids  <- ids[not.na]
         }
         if (na.act == "na.exclude" || na.act == "na.pass") {
            len.ids  <- length(x$ids)
            ids.infl <- x$is.infl
            lab.ids  <- ids
         }
         zi.min <- min(zi,na.rm=TRUE)*1.05
         zi.max <- max(zi,na.rm=TRUE)*1.05
         lplot(NA, NA, xlim=c(1,len.ids), ylim=c(zi.min,zi.max), xaxt="n", main=paste("dfbs: ", names(x$dfbs)[which.dfbs[i]]), xlab="", ylab="", las=las, ...)
         laxis(side=1, at=seq_len(len.ids), labels=lab.ids, xlab="", las=las, ...)
         labline(h= 0, lty="dashed", ...)
         labline(h= 1, lty="dotted", ...)
         labline(h=-1, lty="dotted", ...)
         if (na.act == "na.exclude" || na.act == "na.pass")
            llines(seq_len(len.ids)[not.na], zi[not.na], col=col.na, ...)
         llines(seq_len(len.ids), zi, ...)
         lpoints(x=seq_len(len.ids),           y=zi,           bg=bg,      pch=pch, ...)
         lpoints(x=seq_len(len.ids)[ids.infl], y=zi[ids.infl], bg=bg.infl, pch=pch, ...)

      }

   }

   #########################################################################

   invisible()

}
