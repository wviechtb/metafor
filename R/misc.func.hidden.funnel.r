############################################################################

.funnel.legend <- function(legend, level, shade, back, yaxis, trimfill, pch, col, bg, pch.fill, pch.vec, col.vec, bg.vec, colci) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   lopts <- list(x         = "topright",
                 y         = NULL,
                 inset     = 0.01,
                 bty       = "o",
                 bg        = .coladj(par("bg","fg"), dark=c(0,-0.9), light=c(0,0.9)),
                 studies   = TRUE,
                 show      = "pvals",
                 cex       = c(1,2,1),
                 x.intersp = 1,
                 y.intersp = 1)

   if (is.list(legend)) {

      ### replace defaults with any user-defined values
      lopts.pos <- pmatch(names(legend), names(lopts))
      lopts[c(na.omit(lopts.pos))] <- legend[!is.na(lopts.pos)]
      legend <- TRUE

      if (length(lopts$cex) == 1L)
         lopts$cex <- c(lopts$cex, 2*lopts$cex, lopts$cex)
      if (length(lopts$cex) == 2L)
         lopts$cex <- c(lopts$cex[1], lopts$cex[2], lopts$cex[1])

   } else {

      if (is.character(legend)) {
         lopts$x <- legend
         legend <- TRUE
      } else {
         if (!is.logical(legend))
            stop(mstyle$stop("Argument 'legend' must either be logical, a string, or a list."), call.=FALSE)
      }

   }

   if (!is.na(lopts$show) && !is.element(lopts$show, c("pvals","cis")))
      stop(mstyle$stop("Valid options for 'show' are 'pvals, 'cis', or NA."), call.=FALSE)

   ### can only add p-values / CI regions if 'yaxis' is 'sei', 'vi', 'seinv', or 'vinv'

   if (legend && !is.element(yaxis, c("sei", "vi", "seinv", "vinv")))
      lopts$show <- NA

   ### only add 'Studies' to legend if pch, col, and bg are not vectors

   if (pch.vec || col.vec || bg.vec)
      lopts$studies <- FALSE

   ### if neither studies nor p-values / CI regions are shown, then omit the legend

   if (!lopts$studies && is.na(lopts$show))
      legend <- FALSE

   if (legend) {

      ltxt   <- NULL
      pch.l  <- NULL
      col.l  <- NULL
      pt.cex <- NULL
      pt.bg  <- NULL

      if (isTRUE(lopts$show == "pvals")) {

         level <- c(level, 0)
         lvals <- length(level)

         scipen <- options(scipen=100)
         lchars <- max(nchar(level))-2L
         options(scipen=scipen$scipen)

         ltxt <- sapply(seq_len(lvals), function(i) {
            if (i == 1)
               return(as.expression(bquote(paste(.(pval1) < p, phantom() <= .(pval2)), list(pval1=fmtx(level[i], lchars), pval2=fmtx(1, lchars)))))
            if (i > 1 && i < lvals)
               return(as.expression(bquote(paste(.(pval1) < p, phantom() <= .(pval2)), list(pval1=fmtx(level[i], lchars), pval2=fmtx(level[i-1], lchars)))))
            if (i == lvals)
               return(as.expression(bquote(paste(.(pval1) < p, phantom() <= .(pval2)), list(pval1=fmtx(0, lchars), pval2=fmtx(level[i-1], lchars)))))
         })

         pch.l  <- rep(22, lvals)
         col.l  <- rep(colci, lvals)
         pt.cex <- rep(lopts$cex[2], lvals)
         pt.bg  <- c(shade, back)

      }

      if (isTRUE(lopts$show == "cis")) {

         level <- 100-100*level
         lvals <- length(level)

         scipen <- options(scipen=100)
         lchars <- max(nchar(level))-2L
         options(scipen=scipen$scipen)

         ltxt <- sapply(seq_len(lvals), function(i) as.expression(bquote(paste(.(ci)*"% CI Region"), list(ci=fmtx(level[i], lchars)))))

         pch.l  <- rep(22, lvals)
         col.l  <- rep(colci, lvals)
         pt.cex <- rep(lopts$cex[2], lvals)
         pt.bg  <- c(shade)

      }

      if (isTRUE(lopts$studies)) {

         if (trimfill) {
            ltxt <- c(ltxt, expression(plain(Observed~Studies)))
         } else {
            ltxt <- c(ltxt, expression(plain(Studies)))
         }
         pch.l  <- c(pch.l, pch[1])
         col.l  <- c(col.l, col[1])
         pt.cex <- c(pt.cex, lopts$cex[3])
         pt.bg  <- c(pt.bg, bg[1])

         if (trimfill) {
            ltxt   <- c(ltxt, expression(plain(Imputed~Studies)))
            pch.l  <- c(pch.l, pch.fill[1])
            col.l  <- c(col.l, col[2])
            pt.cex <- c(pt.cex, lopts$cex[3])
            pt.bg  <- c(pt.bg, bg[2])
         }

      }

      legend(x=lopts$x, y=lopts$y, inset=lopts$inset, bty=lopts$bty, bg=lopts$bg,
             cex=lopts$cex[1], x.intersp=lopts$x.intersp, y.intersp=lopts$y.intersp,
             pch=pch.l, col=col.l, pt.cex=pt.cex, pt.bg=pt.bg, legend=ltxt)

   }

}

############################################################################
