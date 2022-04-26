plot.permutest.rma.uni <- function(x, beta, alpha, QM=FALSE, QS=FALSE,
   breaks="Scott", freq=FALSE, col="lightgray", border=NULL,
   col.out=rgb(1,0,0,0.5), col.ref="black", col.density="blue", adjust=1,
   lwd=c(2,0,0,4), layout, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="permutest.rma.uni")

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

   # 1st: obs stat, 2nd: ref dist, 3rd: density, 4th: refline

   if (length(lwd) == 1L)
      lwd <- c(lwd[c(1,1,1)], 4)
   if (length(lwd) == 2L)
      lwd <- c(lwd[c(1,2,2)], 4)
   if (length(lwd) == 3L)
      lwd <- c(lwd[c(1,2,2,3)])

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

   coltail <- function(h, val, tail="lower", mult=1, col, border, freq, ...) {

      h$counts  <- h$counts  * mult
      h$density <- h$density * mult

      if (tail == "lower") {

         above <- which(h$breaks > val)
         if (length(above) > 0L) {
            pos <- above[1]
            h$breaks[pos] <- val
         }
         sel <- h$breaks <= val
         if (sum(sel) >= 2L) {
            h$breaks  <- h$breaks[sel]
            h$counts  <- h$counts[sel[-1]]
            h$density <- h$density[sel[-1]]
            h$mids    <- h$mids[sel[-1]]
            lines(h, col=col, border=border, freq=freq, ...)
         }

      } else {

         below <- which(h$breaks < val)
         if (length(below) > 0L) {
            pos <- below[length(below)]
            h$breaks[pos] <- val
         }
         sel <- h$breaks >= val
         if (sum(sel) >= 2L) {
            len <- length(below)
            h$breaks  <- h$breaks[sel]
            h$counts  <- h$counts[sel[-len]]
            h$density <- h$density[sel[-len]]
            h$mids    <- h$mids[sel[-len]]
            lines(h, col=col, border=border, freq=freq, ...)
         }

      }

   }

   par.mfrow <- par("mfrow")
   par(mfrow=layout)
   on.exit(par(mfrow = par.mfrow), add=TRUE)

   if (!is.null(QM.perm)) {

      tmp <- lhist(QM.perm, breaks=breaks, col=col, border=border,
                   main=ifelse(inherits(x, "permutest.rma.ls"), "Omnibus Test of Location Coefficients", "Omnibus Test of Coefficients"),
                   xlab="Value of Test Statistic",
                   freq=freq, ...)

      coltail(tmp, val=x$QM, tail="upper", col=col.out, border=border, freq=freq, ...)
      labline(v=x$QM, lwd=lwd[1], lty="dashed", ...)

      if (is.na(x$ddf)) {
         xs <- seq(0, max(qchisq(.995, df=length(x$btt)), max(QM.perm, na.rm=TRUE)), length=1000)
         llines(xs, dchisq(xs, df=length(x$btt)), lwd=lwd[2], col=col.ref, ...)
      } else {
         xs <- seq(0, max(qf(.995, df1=length(x$btt), df2=x$ddf), max(QM.perm, na.rm=TRUE)), length=1000)
         llines(xs, df(xs, df1=length(x$btt), df2=x$ddf), lwd=lwd[2], col=col.ref, ...)
      }

      llines(density(QM.perm, adjust=adjust, na.rm=TRUE), lwd=lwd[3], col=col.density, ...)

   }

   for (i in seq_len(ncol(perm1))) {

      tmp <- lhist(perm1[[i]], breaks=breaks, col=col, border=border,
                   main=ifelse(x$int.only, "", paste0(ifelse(inherits(x, "permutest.rma.ls"), "Location Coefficient: ", "Coefficient: "), names(perm1)[i])),
                   xlab=ifelse(stat == "test", "Value of Test Statistic", "Value of Coefficient"),
                   freq=freq, ...)

      if (alternative == "two.sided") {

         if (p2defn == "abs") {

            coltail(tmp, val=-abs(obs1[i]), tail="lower", col=col.out, border=border, freq=freq, ...)
            coltail(tmp, val= abs(obs1[i]), tail="upper", col=col.out, border=border, freq=freq, ...)
            labline(v=c(-obs1[i],obs1[i]), lwd=lwd[1], lty="dashed", ...)

         } else {

            if (obs1[i] > median(perm1[[i]], na.rm=TRUE)) {

               coltail(tmp, val= abs(obs1[i]), tail="upper", mult=2, col=col.out, border=border, freq=freq, ...)
               labline(v=obs1[i], lwd=lwd[1], lty="dashed", ...)

            } else {

               coltail(tmp, val=-abs(obs1[i]), tail="lower", mult=2, col=col.out, border=border, freq=freq, ...)
               labline(v=-abs(obs1[i]), lwd=lwd[1], lty="dashed", ...)

            }

         }
      }

      if (alternative == "less") {

         coltail(tmp, val=obs1[i], tail="lower", col=col.out, border=border, freq=freq, ...)
         labline(v=obs1[i], lwd=lwd[1], lty="dashed", ...)

      }

      if (alternative == "greater") {

         coltail(tmp, val=obs1[i], tail="upper", col=col.out, border=border, freq=freq, ...)
         labline(v=obs1[i], lwd=lwd[1], lty="dashed", ...)

      }

      if (is.na(x$ddf)) {
         xs <- seq(min(-qnorm(.995), min(perm1[[i]], na.rm=TRUE)), max(qnorm(.995), max(perm1[[i]], na.rm=TRUE)), length=1000)
         llines(xs, dnorm(xs), lwd=lwd[2], col=col.ref, ...)
      } else {
         xs <- seq(min(-qt(.995, df=x$ddf), min(perm1[[i]], na.rm=TRUE)), max(qt(.995, df=x$ddf), max(perm1[[i]], na.rm=TRUE)), length=1000)
         llines(xs, dt(xs, df=x$ddf), lwd=lwd[2], col=col.ref, ...)
      }

      llines(density(perm1[[i]], adjust=adjust, na.rm=FALSE), lwd=lwd[3], col=col.density, ...)

      labline(v=0, lwd=lwd[4], ...)

   }

   if (inherits(x, "permutest.rma.ls")) {

      if (!is.null(QS.perm)) {

         tmp <- lhist(QS.perm, breaks=breaks, col=col, border=border,
                      main="Omnibus Test of Scale Coefficients",
                      xlab="Value of Test Statistic",
                      freq=freq, ...)

         coltail(tmp, val=x$QS, tail="upper", col=col.out, border=border, freq=freq, ...)
         labline(v=x$QS, lwd=lwd[1], lty="dashed", ...)

         if (is.na(x$ddf.alpha)) {
            xs <- seq(0, max(qchisq(.995, df=length(x$att)), max(QS.perm, na.rm=TRUE)), length=1000)
            llines(xs, dchisq(xs, df=length(x$att)), lwd=lwd[2], col=col.ref, ...)
         } else {
            xs <- seq(0, max(qf(.995, df1=length(x$att), df2=x$ddf.alpha), max(QS.perm, na.rm=TRUE)), length=1000)
            llines(xs, df(xs, df1=length(x$att), df2=x$ddf.alpha), lwd=lwd[2], col=col.ref, ...)
         }

         llines(density(QS.perm, adjust=adjust, na.rm=TRUE), lwd=lwd[3], col=col.density, ...)

      }

      for (i in seq_len(ncol(perm2))) {

         tmp <- lhist(perm2[[i]], breaks=breaks, col=col, border=border,
                      main=ifelse(x$Z.int.only, "", paste0("Scale Coefficient: ", names(perm2)[i])),
                      xlab=ifelse(stat == "test", "Value of Test Statistic", "Value of Coefficient"),
                      freq=freq, ...)

         if (alternative == "two.sided") {

            if (p2defn == "abs") {

               coltail(tmp, val=-abs(obs2[i]), tail="lower", col=col.out, border=border, freq=freq, ...)
               coltail(tmp, val= abs(obs2[i]), tail="upper", col=col.out, border=border, freq=freq, ...)
               labline(v=c(-obs2[i],obs2[i]), lwd=lwd[1], lty="dashed", ...)

            } else {

               if (obs2[i] > median(perm2[[i]], na.rm=TRUE)) {

                  coltail(tmp, val= abs(obs2[i]), tail="upper", mult=2, col=col.out, border=border, freq=freq, ...)
                  labline(v=obs2[i], lwd=lwd[1], lty="dashed", ...)

               } else {

                  coltail(tmp, val=-abs(obs2[i]), tail="lower", mult=2, col=col.out, border=border, freq=freq, ...)
                  labline(v=-abs(obs2[i]), lwd=lwd[1], lty="dashed", ...)

               }

            }
         }

         if (alternative == "less") {

            coltail(tmp, val=obs2[i], tail="lower", col=col.out, border=border, freq=freq, ...)
            labline(v=obs2[i], lwd=lwd[1], lty="dashed", ...)

         }

         if (alternative == "greater") {

            coltail(tmp, val=obs2[i], tail="upper", col=col.out, border=border, freq=freq, ...)
            labline(v=obs2[i], lwd=lwd[1], lty="dashed", ...)

         }

         if (is.na(x$ddf.alpha)) {
            xs <- seq(min(-qnorm(.995), min(perm2[[i]], na.rm=TRUE)), max(qnorm(.995), max(perm2[[i]], na.rm=TRUE)), length=1000)
            llines(xs, dnorm(xs), lwd=lwd[2], col=col.ref, ...)
         } else {
            xs <- seq(min(-qt(.995, df=x$ddf.alpha), min(perm2[[i]], na.rm=TRUE)), max(qt(.995, df=x$ddf.alpha), max(perm2[[i]], na.rm=TRUE)), length=1000)
            llines(xs, dt(xs, df=x$ddf.alpha), lwd=lwd[2], col=col.ref, ...)
         }

         llines(density(perm2[[i]], adjust=adjust, na.rm=TRUE), lwd=lwd[3], col=col.density, ...)

         labline(v=0, lwd=lwd[4], ...)

      }

   }

   ############################################################################

   invisible()

}
