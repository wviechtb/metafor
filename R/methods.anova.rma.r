############################################################################

as.data.frame.anova.rma <- function(x, ...) {

   .chkclass(class(x), must="anova.rma")

   if (x$type == "Wald.btt") {

      tab <- data.frame(coefs = .format.btt(x$btt),
                        QM    = x$QM,
                        df    = round(x$QMdf[1], 2),
                        pval  = x$QMp)

      if (is.element(x$test, c("knha","adhoc","t"))) {
         names(tab)[2:3] <- c("Fval", "df1")
         tab <- cbind(tab[1:3], df2 = round(x$QMdf[2], 2), tab[4])
      }

   }

   if (x$type == "Wald.att") {

      tab <- data.frame(coefs = .format.btt(x$att),
                        QS    = x$QS,
                        df    = round(x$QSdf[1], 2),
                        pval  = x$QSp)

      if (is.element(x$test, c("knha","adhoc","t"))) {
         names(tab)[2:3] <- c("Fval", "df1")
         tab <- cbind(tab[1:3], df2 = round(x$QSdf[2], 2), tab[4])
      }

   }

   if (x$type == "Wald.Xb") {

      if (is.element(x$test, c("knha","adhoc","t"))) {
         tab <- data.frame(hyp=x$hyp[[1]], estimate=c(x$Xb), se=x$se, tval=x$zval, df=round(x$ddf,2), pval=x$pval)
      } else {
         tab <- data.frame(hyp=x$hyp[[1]], estimate=c(x$Xb), se=x$se, zval=x$zval, pval=x$pval)
      }
      rownames(tab) <- paste0(seq_len(x$m), ":")

      return(tab)

   }

   if (x$type == "Wald.Za") {

      if (is.element(x$test, c("knha","adhoc","t"))) {
         tab <- data.frame(hyp=x$hyp[[1]], estimate=c(x$Za), se=x$se, tval=x$zval, df=round(x$ddf,2), pval=x$pval)
      } else {
         tab <- data.frame(hyp=x$hyp[[1]], estimate=c(x$Za), se=x$se, zval=x$zval, pval=x$pval)
      }
      rownames(tab) <- paste0(seq_len(x$m), ":")

      return(tab)

   }

   if (x$type == "LRT") {

      tab <- data.frame(c(x$parms.f, x$parms.r),
                        c(x$fit.stats.f["AIC"], x$fit.stats.r["AIC"]),
                        c(x$fit.stats.f["BIC"], x$fit.stats.r["BIC"]),
                        c(x$fit.stats.f["AICc"], x$fit.stats.r["AICc"]),
                        c(x$fit.stats.f["ll"], x$fit.stats.r["ll"]),
                        c(NA_real_, x$LRT),
                        c(NA_real_, x$pval),
                        c(x$QE.f, x$QE.r),
                        c(x$tau2.f, x$tau2.r),
                        c(NA_real_, NA_real_))

      colnames(tab) <- c("df", "AIC", "BIC", "AICc", "logLik", "LRT", "pval", "QE", "tau^2", "R^2")
      rownames(tab) <- c("Full", "Reduced")

      tab["Full",c("LRT","pval")] <- NA_real_
      tab["Full","R^2"] <- NA_real_
      tab["Reduced","R^2"] <- x$R2

      ### remove tau^2 column if full model is a FE/EE/CE model or tau2.f/tau2.r is NA

      if (is.element(x$method, c("FE","EE","CE")) || (is.na(x$tau2.f) || is.na(x$tau2.r)))
         tab <- tab[-which(names(tab) == "tau^2")]

      ### remove R^2 column if full model is a rma.mv or rma.ls model

      if (is.element("rma.mv", x$class.f) || is.element("rma.ls", x$class.f))
         tab <- tab[-which(names(tab) == "R^2")]

   }

   return(tab)

}

as.data.frame.list.anova.rma <- function(x, ...) {

   .chkclass(class(x), must="list.anova.rma")

   if (x[[1]]$type == "Wald.btt") {

      tab <- data.frame(spec  = names(x),
                        coefs = sapply(x, function(x) .format.btt(x$btt)),
                        QM    = sapply(x, function(x) x$QM),
                        df    = sapply(x, function(x) round(x$QMdf[1], 2)),
                        pval  = sapply(x, function(x) x$QMp))

   }

   if (x[[1]]$type == "Wald.att") {

      tab <- data.frame(spec  = names(x),
                        coefs = sapply(x, function(x) .format.btt(x$att)),
                        QS    = sapply(x, function(x) x$QS),
                        df    = sapply(x, function(x) round(x$QSdf[1], 2)),
                        pval  = sapply(x, function(x) x$QSp))

   }

   if (is.element(x[[1]]$test, c("knha","adhoc","t"))) {
      names(tab)[3:4] <- c("Fval", "df1")
      if (x[[1]]$type == "Wald.btt")
         tab <- cbind(tab[1:4], df2 = sapply(x, function(x) round(x$QMdf[2], 2)), tab[5])
      if (x[[1]]$type == "Wald.att")
         tab <- cbind(tab[1:4], df2 = sapply(x, function(x) round(x$QSdf[2], 2)), tab[5])
   }

   # if all btt/att specifications are numeric, remove the 'spec' column
   if (all(substr(tab$spec, 1, 1) %in% as.character(1:9)))
       tab$spec <- NULL

   # just use numbers for row names
   rownames(tab) <- NULL

   return(tab)

}

############################################################################
