plot.permutest.rma.uni <- function(x, beta, alpha, QM=FALSE, QS=FALSE,
   breaks="Scott", freq=FALSE, col, border, col.out, col.ref, col.density,
   trim=0, adjust=1, lwd=c(2,0,0,4), layout, legend=FALSE, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="permutest.rma.uni")

   .start.plot()

   if (missing(col))
      col <- .coladj(par("bg","fg"), dark=0.3, light=-0.3)

   if (missing(border))
      border <- .coladj(par("bg"), dark=0.1, light=-0.1)

   if (missing(col.out))
      col.out <- rgb(1,0,0,0.5)

   if (missing(col.ref))
      col.ref <- .coladj(par("fg"), dark=-0.3, light=0.3)

   if (missing(col.density))
      col.density <- ifelse(.is.dark(), "dodgerblue", "blue")

   ddd <- list(...)

   if (is.null(ddd$alternative)) {
      alternative <- x$alternative
   } else {
      alternative <- match.arg(ddd$alternative, c("two.sided", "less", "greater"))
   }

   if (is.null(ddd$p2defn)) {
      p2defn <- x$p2defn
   } else {
      p2defn <- match.arg(ddd$p2defn, c("abs", "px2"))
   }

   if (is.null(ddd$stat)) {
      stat <- x$stat
   } else {
      stat <- match.arg(ddd$stat, c("test", "coef"))
   }

   ### check trim

   if (trim >= 0.5)
      stop(mstyle$stop("The value of 'trim' must be < 0.5."))

   # 1st: obs stat, 2nd: ref dist, 3rd: density, 4th: refline

   if (length(lwd) == 1L)
      lwd <- c(lwd[c(1,1,1)], 4)
   if (length(lwd) == 2L)
      lwd <- c(lwd[c(1,2,2)], 4)
   if (length(lwd) == 3L)
      lwd <- c(lwd[c(1,2,2,3)])

   # cannot plot ref dist and density when freq=TRUE

   if (freq)
      lwd[c(2,3)] <- 0

   lhist   <- function(..., alternative, p2defn, stat) hist(...)
   labline <- function(..., alternative, p2defn, stat) abline(...)
   llines  <- function(..., alternative, p2defn, stat) lines(...)

   ############################################################################

   if (x$skip.beta) {
      beta <- NULL
   } else {
      if (missing(beta)) {
         if (x$int.only) {
            beta <- 1
         } else {
            if (x$int.incl) {
               beta <- 2:x$p
            } else {
               beta <- 1:x$p
            }
         }
      } else {
         if (all(is.na(beta))) { # set beta=NA to not plot any location coefficients
            beta <- NULL
         } else {
            beta <- .set.btt(beta, x$p, x$int.incl, names(x$zval.perm))
         }
      }
   }

   if (stat == "test") {
      perm1 <- x$zval.perm[beta]
      obs1  <- x$zval[beta]
   } else {
      perm1 <- x$beta.perm[beta]
      obs1  <- x$beta[beta,1]
   }

   if (x$int.only || x$skip.beta) {
      QM.perm <- NULL
   } else {
      if (QM) {
         QM.perm <- x$QM.perm
      } else {
         QM.perm <- NULL
      }
   }

   if (inherits(x, "permutest.rma.ls") && !x$skip.alpha) {

      if (missing(alpha)) {
         if (x$Z.int.only) {
            alpha <- 1
         } else {
            if (x$int.incl) {
               alpha <- 2:x$q
            } else {
               alpha <- 1:x$q
            }
         }
      } else {
         if (all(is.na(alpha))) { # set alpha=NA to not plot any scale coefficients
            alpha <- NULL
         } else {
            alpha <- .set.btt(alpha, x$q, x$Z.int.incl, names(x$zval.perm.alpha))
         }
      }

      if (stat == "test") {
         perm2 <- x$zval.alpha.perm[alpha]
         obs2  <- x$zval.alpha[alpha]
      } else {
         perm2 <- x$alpha.perm[alpha]
         obs2  <- x$alpha[alpha,1]
      }

      if (QS) {
         QS.perm <- x$QS.perm
      } else {
         QS.perm <- NULL
      }

   } else {

      alpha <- NULL
      QS.perm <- NULL

   }

   ############################################################################

   ### function to add legend

   addlegend <- function(legend) {

      if (is.logical(legend) && isTRUE(legend))
         lpos <- "topright"

      if (is.character(legend)) {
         lpos <- legend
         legend <- TRUE
      }

      if (legend && any(lwd[2:3] > 0)) {

         ltxt  <- c("Kernel Density Estimate of\nthe Permutation Distribution", "Theoretical Null Distribution")
         lwds  <- lwd[2:3]
         lcols <- c(col.density, col.ref)
         ltys  <- c("solid", "solid")
         #pchs  <- c("","","\u2506") # \u250a
         ltxt  <- ltxt[lwds > 0]
         lcols <- lcols[lwds > 0]
         ltys  <- ltys[lwds > 0]
         #pchs  <- pchs[lwds > 0]
         lwds  <- lwds[lwds > 0]
         legend(lpos, inset=.01, bg=.coladj(par("bg"), dark=0, light=0), lwd=lwds, col=lcols, lty=ltys, legend=ltxt)

      }

      return(FALSE)

   }

   ############################################################################

   # number of plots

   np <- length(beta) + length(alpha) + ifelse(is.null(QM.perm), 0L, 1L) + ifelse(is.null(QS.perm), 0L, 1L)

   if (np == 0L)
      stop(mstyle$stop("Must select at least one elements to plot."))

   # set/check layout argument

   if (missing(layout)) {
      layout <- n2mfrow(np)
   } else {
      layout <- layout[layout >= 1]
      layout <- round(layout)
      if (length(layout) != 2L)
         stop(mstyle$stop("Incorrect specification of 'layout' argument."))
   }

   #print(list(np, layout))

   ############################################################################

   par.mfrow <- par("mfrow")
   par(mfrow=layout)
   on.exit(par(mfrow = par.mfrow), add=TRUE)

   if (!is.null(QM.perm)) {

      pdist <- QM.perm

      if (is.na(x$ddf)) {
         xs <- seq(0, max(qchisq(.995, df=length(x$btt)), max(pdist, na.rm=TRUE)), length=1000)
         ys <- dchisq(xs, df=length(x$btt))
      } else {
         xs <- seq(0, max(qf(.995, df1=length(x$btt), df2=x$ddf), max(pdist, na.rm=TRUE)), length=1000)
         ys <- df(xs, df1=length(x$btt), df2=x$ddf)
      }

      den <- density(pdist, adjust=adjust, na.rm=TRUE, n=8192)

      if (trim > 0) {
         bound <- quantile(pdist, probs=1-trim, na.rm=TRUE)
         pdist <- pdist[pdist <= bound]
      }

      if (lwd[2] == 0 && lwd[3] == 0) {

         tmp <- lhist(pdist, breaks=breaks, col=col, border=border,
                      main=ifelse(inherits(x, "permutest.rma.ls"), "Omnibus Test of Location Coefficients", "Omnibus Test of Coefficients"),
                      xlab="Value of Test Statistic",
                      freq=freq, ...)

      } else {

         tmp <- lhist(pdist, breaks=breaks, plot=FALSE)

         ylim <- c(0, max(ifelse(lwd[2] == 0, 0, max(ys)), ifelse(lwd[3] == 0, 0, max(den$y)), max(tmp$density)))

         tmp <- lhist(pdist, breaks=breaks, col=col, border=border,
                      main=ifelse(inherits(x, "permutest.rma.ls"), "Omnibus Test of Location Coefficients", "Omnibus Test of Coefficients"),
                      xlab="Value of Test Statistic",
                      freq=freq, ylim=ylim, ...)

      }

      .coltail(tmp, val=x$QM, tail="upper", col=col.out, border=border, freq=freq, ...)

      if (lwd[1] > 0)
         labline(v=x$QM, lwd=lwd[1], lty="dashed", ...)
      if (lwd[2] > 0)
         llines(xs, ys, lwd=lwd[2], col=col.ref, ...)
      if (lwd[3] > 0)
         llines(den, lwd=lwd[3], col=col.density, ...)

      legend <- addlegend(legend)

   }

   for (i in seq_len(ncol(perm1))) {

      pdist <- perm1[[i]]

      if (is.na(x$ddf)) {
         xs <- seq(min(-qnorm(.995), min(pdist, na.rm=TRUE)), max(qnorm(.995), max(pdist, na.rm=TRUE)), length=1000)
         ys <- dnorm(xs)
      } else {
         xs <- seq(min(-qt(.995, df=x$ddf), min(pdist, na.rm=TRUE)), max(qt(.995, df=x$ddf), max(pdist, na.rm=TRUE)), length=1000)
         ys <- dt(xs, df=x$ddf)
      }

      den <- density(pdist, adjust=adjust, na.rm=TRUE, n=8192)

      if (trim > 0) {
         bounds <- quantile(pdist, probs=c(trim/2, 1-trim/2), na.rm=TRUE)
         pdist <- pdist[pdist >= bounds[1] & pdist <= bounds[2]]
      }

      if (lwd[2] == 0 && lwd[3] == 0) {

         tmp <- lhist(pdist, breaks=breaks, col=col, border=border,
                      main=ifelse(x$int.only, "", paste0(ifelse(inherits(x, "permutest.rma.ls"), "Location Coefficient: ", "Coefficient: "), names(perm1)[i])),
                      xlab=ifelse(stat == "test", "Value of Test Statistic", "Value of Coefficient"),
                      freq=freq, ...)

      } else {

         tmp <- lhist(pdist, breaks=breaks, plot=FALSE)

         ylim <- c(0, max(ifelse(lwd[2] == 0, 0, max(ys)), ifelse(lwd[3] == 0, 0, max(den$y)), max(tmp$density)))

         tmp <- lhist(pdist, breaks=breaks, col=col, border=border,
                      main=ifelse(x$int.only, "", paste0(ifelse(inherits(x, "permutest.rma.ls"), "Location Coefficient: ", "Coefficient: "), names(perm1)[i])),
                      xlab=ifelse(stat == "test", "Value of Test Statistic", "Value of Coefficient"),
                      freq=freq, ylim=ylim, ...)

      }

      if (alternative == "two.sided") {

         if (p2defn == "abs") {

            .coltail(tmp, val=-abs(obs1[i]), tail="lower", col=col.out, border=border, freq=freq, ...)
            .coltail(tmp, val= abs(obs1[i]), tail="upper", col=col.out, border=border, freq=freq, ...)

            if (lwd[1] > 0)
               labline(v=c(-obs1[i],obs1[i]), lwd=lwd[1], lty="dashed", ...)

         } else {

            if (obs1[i] > median(pdist, na.rm=TRUE)) {

               .coltail(tmp, val= abs(obs1[i]), tail="upper", mult=2, col=col.out, border=border, freq=freq, ...)

               if (lwd[1] > 0)
                  labline(v=obs1[i], lwd=lwd[1], lty="dashed", ...)

            } else {

               .coltail(tmp, val=-abs(obs1[i]), tail="lower", mult=2, col=col.out, border=border, freq=freq, ...)

               if (lwd[1] > 0)
                  labline(v=-abs(obs1[i]), lwd=lwd[1], lty="dashed", ...)

            }

         }
      }

      if (alternative == "less") {

         .coltail(tmp, val=obs1[i], tail="lower", col=col.out, border=border, freq=freq, ...)

         if (lwd[1] > 0)
            labline(v=obs1[i], lwd=lwd[1], lty="dashed", ...)

      }

      if (alternative == "greater") {

         .coltail(tmp, val=obs1[i], tail="upper", col=col.out, border=border, freq=freq, ...)

         if (lwd[1] > 0)
            labline(v=obs1[i], lwd=lwd[1], lty="dashed", ...)

      }

      if (lwd[2] > 0)
         llines(xs, ys, lwd=lwd[2], col=col.ref, ...)
      if (lwd[3] > 0)
         llines(den, lwd=lwd[3], col=col.density, ...)

      if (lwd[4] > 0)
         labline(v=0, lwd=lwd[4], ...)

      legend <- addlegend(legend)

   }

   if (inherits(x, "permutest.rma.ls")) {

      if (!is.null(QS.perm)) {

         pdist <- QS.perm

         if (is.na(x$ddf.alpha)) {
            xs <- seq(0, max(qchisq(.995, df=length(x$att)), max(pdist, na.rm=TRUE)), length=1000)
            ys <- dchisq(xs, df=length(x$att))
         } else {
            xs <- seq(0, max(qf(.995, df1=length(x$att), df2=x$ddf.alpha), max(pdist, na.rm=TRUE)), length=1000)
            ys <- df(xs, df1=length(x$att), df2=x$ddf.alpha)
         }

         den <- density(pdist, adjust=adjust, na.rm=TRUE, n=8192)

         if (trim > 0) {
            bound <- quantile(pdist, probs=1-trim, na.rm=TRUE)
            pdist <- pdist[pdist <= bound]
         }

         if (lwd[2] == 0 && lwd[3] == 0) {

            tmp <- lhist(pdist, breaks=breaks, col=col, border=border,
                         main="Omnibus Test of Scale Coefficients",
                         xlab="Value of Test Statistic",
                         freq=freq, ...)

         } else {

            tmp <- lhist(pdist, breaks=breaks, plot=FALSE)

            ylim <- c(0, max(ifelse(lwd[2] == 0, 0, max(ys)), ifelse(lwd[3] == 0, 0, max(den$y)), max(tmp$density)))

            tmp <- lhist(pdist, breaks=breaks, col=col, border=border,
                         main="Omnibus Test of Scale Coefficients",
                         xlab="Value of Test Statistic",
                         freq=freq, ylim=ylim, ...)

         }

         .coltail(tmp, val=x$QS, tail="upper", col=col.out, border=border, freq=freq, ...)

         if (lwd[1] > 0)
            labline(v=x$QS, lwd=lwd[1], lty="dashed", ...)
         if (lwd[2] > 0)
            llines(xs, ys, lwd=lwd[2], col=col.ref, ...)
         if (lwd[3] > 0)
            llines(den, lwd=lwd[3], col=col.density, ...)

         legend <- addlegend(legend)

      }

      for (i in seq_len(ncol(perm2))) {

         pdist <- perm2[[i]]

         if (is.na(x$ddf.alpha)) {
            xs <- seq(min(-qnorm(.995), min(pdist, na.rm=TRUE)), max(qnorm(.995), max(pdist, na.rm=TRUE)), length=1000)
            ys <- dnorm(xs)
         } else {
            xs <- seq(min(-qt(.995, df=x$ddf.alpha), min(pdist, na.rm=TRUE)), max(qt(.995, df=x$ddf.alpha), max(pdist, na.rm=TRUE)), length=1000)
            ys <- dt(xs, df=x$ddf.alpha)
         }

         den <- density(pdist, adjust=adjust, na.rm=TRUE, n=8192)

         if (trim > 0) {
            bounds <- quantile(pdist, probs=c(trim/2, 1-trim/2), na.rm=TRUE)
            pdist <- pdist[pdist >= bounds[1] & pdist <= bounds[2]]
         }

         if (lwd[2] == 0 && lwd[3] == 0) {

            tmp <- lhist(pdist, breaks=breaks, col=col, border=border,
                         main=ifelse(x$Z.int.only, "", paste0("Scale Coefficient: ", names(perm2)[i])),
                         xlab=ifelse(stat == "test", "Value of Test Statistic", "Value of Coefficient"),
                         freq=freq, ...)

         } else {

            tmp <- lhist(pdist, breaks=breaks, plot=FALSE)

            ylim <- c(0, max(ifelse(lwd[2] == 0, 0, max(ys)), ifelse(lwd[3] == 0, 0, max(den$y)), max(tmp$density)))

            tmp <- lhist(pdist, breaks=breaks, col=col, border=border,
                         main=ifelse(x$Z.int.only, "", paste0("Scale Coefficient: ", names(perm2)[i])),
                         xlab=ifelse(stat == "test", "Value of Test Statistic", "Value of Coefficient"),
                         freq=freq, ylim=ylim, ...)

         }

         if (alternative == "two.sided") {

            if (p2defn == "abs") {

               .coltail(tmp, val=-abs(obs2[i]), tail="lower", col=col.out, border=border, freq=freq, ...)
               .coltail(tmp, val= abs(obs2[i]), tail="upper", col=col.out, border=border, freq=freq, ...)

               if (lwd[1] > 0)
                  labline(v=c(-obs2[i],obs2[i]), lwd=lwd[1], lty="dashed", ...)

            } else {

               if (obs2[i] > median(pdist, na.rm=TRUE)) {

                  .coltail(tmp, val= abs(obs2[i]), tail="upper", mult=2, col=col.out, border=border, freq=freq, ...)

                  if (lwd[1] > 0)
                     labline(v=obs2[i], lwd=lwd[1], lty="dashed", ...)

               } else {

                  .coltail(tmp, val=-abs(obs2[i]), tail="lower", mult=2, col=col.out, border=border, freq=freq, ...)

                  if (lwd[1] > 0)
                     labline(v=-abs(obs2[i]), lwd=lwd[1], lty="dashed", ...)

               }

            }
         }

         if (alternative == "less") {

            .coltail(tmp, val=obs2[i], tail="lower", col=col.out, border=border, freq=freq, ...)

            if (lwd[1] > 0)
               labline(v=obs2[i], lwd=lwd[1], lty="dashed", ...)

         }

         if (alternative == "greater") {

            .coltail(tmp, val=obs2[i], tail="upper", col=col.out, border=border, freq=freq, ...)

            if (lwd[1] > 0)
               labline(v=obs2[i], lwd=lwd[1], lty="dashed", ...)

         }

         if (lwd[2] > 0)
            llines(xs, ys, lwd=lwd[2], col=col.ref, ...)
         if (lwd[3] > 0)
            llines(den, lwd=lwd[3], col=col.density, ...)
         if (lwd[4] > 0)
            labline(v=0, lwd=lwd[4], ...)

         legend <- addlegend(legend)

      }

   }

   ############################################################################

   invisible()

}
