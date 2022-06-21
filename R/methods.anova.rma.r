############################################################################

as.data.frame.anova.rma <- function(x, ...) {

   .chkclass(class(x), must="anova.rma")

   if (x$type == "Wald.btt") {

      df <- data.frame(coefs = .format.btt(x$btt),
                       QM    = x$QM,
                       df    = round(x$QMdf[1], 2),
                       pval  = x$QMp)

      if (is.element(x$test, c("knha","adhoc","t"))) {
         names(df)[2:3] <- c("Fval", "df1")
         df <- cbind(df[1:3], df2 = round(x$QMdf[2], 2), df[4])
      }

   }

   if (x$type == "Wald.att") {

      df <- data.frame(coefs = .format.btt(x$att),
                       QS    = x$QS,
                       df    = round(x$QSdf[1], 2),
                       pval  = x$QSp)

      if (is.element(x$test, c("knha","adhoc","t"))) {
         names(df)[2:3] <- c("Fval", "df1")
         df <- cbind(df[1:3], df2 = round(x$QSdf[2], 2), df[4])
      }

   }

   if (x$type == "Wald.Xb") {

      if (is.element(x$test, c("knha","adhoc","t"))) {
         df <- data.frame(hyp=x$hyp[[1]], estimate=c(x$Xb), se=x$se, tval=x$zval, df=round(x$ddf,2), pval=x$pval)
      } else {
         df <- data.frame(hyp=x$hyp[[1]], estimate=c(x$Xb), se=x$se, zval=x$zval, pval=x$pval)
      }
      rownames(df) <- paste0(seq_len(x$m), ":")

      return(df)

   }

   if (x$type == "Wald.Za") {

      if (is.element(x$test, c("knha","adhoc","t"))) {
         df <- data.frame(hyp=x$hyp[[1]], estimate=c(x$Za), se=x$se, tval=x$zval, df=round(x$ddf,2), pval=x$pval)
      } else {
         df <- data.frame(hyp=x$hyp[[1]], estimate=c(x$Za), se=x$se, zval=x$zval, pval=x$pval)
      }
      rownames(df) <- paste0(seq_len(x$m), ":")

      return(df)

   }

   if (x$type == "LRT") {

      df <- data.frame(c(x$parms.f, x$parms.r),
                       c(x$fit.stats.f["AIC"], x$fit.stats.r["AIC"]),
                       c(x$fit.stats.f["BIC"], x$fit.stats.r["BIC"]),
                       c(x$fit.stats.f["AICc"], x$fit.stats.r["AICc"]),
                       c(x$fit.stats.f["ll"], x$fit.stats.r["ll"]),
                       c(NA, x$LRT),
                       c(NA, x$pval),
                       c(x$QE.f, x$QE.r),
                       c(x$tau2.f, x$tau2.r),
                       c(NA, NA))

      colnames(df) <- c("df", "AIC", "BIC", "AICc", "logLik", "LRT", "pval", "QE", "tau^2", "R^2")
      rownames(df) <- c("Full", "Reduced")

      df["Full",c("LRT","pval")] <- NA
      df["Full","R^2"] <- NA
      df["Reduced","R^2"] <- x$R2

      ### remove tau^2 column if full model is a FE/EE/CE model or tau2.f/tau2.r is NA

      if (is.element(x$method, c("FE","EE","CE")) || (is.na(x$tau2.f) || is.na(x$tau2.r)))
         df <- df[-which(names(df) == "tau^2")]

      ### remove R^2 column if full model is a rma.mv or rma.ls model

      if (is.element("rma.mv", x$class.f) || is.element("rma.ls", x$class.f))
         df <- df[-which(names(df) == "R^2")]

   }

   return(df)

}

as.data.frame.list.anova.rma <- function(x, ...) {

   .chkclass(class(x), must="list.anova.rma")

   df <- data.frame(spec  = names(x),
                    coefs = sapply(x, function(x) .format.btt(x$btt)),
                    QM    = sapply(x, function(x) x$QM),
                    df    = sapply(x, function(x) round(x$QMdf[1], 2)),
                    pval  = sapply(x, function(x) x$QMp))

   if (is.element(x[[1]]$test, c("knha","adhoc","t"))) {
      names(df)[3:4] <- c("Fval", "df1")
      df <- cbind(df[1:4], df2 = round(x[[1]]$QMdf[2], 2), df[5])
   }

   # if all btt specifications were numeric, remove the 'spec' column
   if (all(substr(df$spec, 1, 1) %in% as.character(1:9)))
       df$spec <- NULL

   # just use numbers for row names
   rownames(df) <- NULL

   return(df)

}

############################################################################
