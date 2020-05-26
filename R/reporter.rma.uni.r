reporter.rma.uni <- function(x, dir, filename, format="html_document", open=TRUE, digits, forest, funnel, footnotes=FALSE, verbose=TRUE, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(x, "rma.uni"))
      stop(mstyle$stop("Argument 'x' must be an object of class \"rma.uni\"."))

   if (!suppressMessages(suppressWarnings(requireNamespace("rmarkdown", quietly=TRUE))))
      stop(mstyle$stop("Please install the 'rmarkdown' package to use the reporter function."))

   if (!is.element(x$test, c("z", "knha")))
      stop(mstyle$stop("Cannot only use reporter function when test='z' or test='knha'."))

   if (x$model == "rma.ls")
      stop(mstyle$stop("Cannot use reporter function for location-scale models."))

   if (!x$weighted)
      stop(mstyle$stop("Cannot use reporter function when 'weighted=FALSE'."))

   if (!is.null(x$weights))
      stop(mstyle$stop("Cannot use reporter function for models with custom weights."))

   if (is.null(x$tau2.fix))
      stop(mstyle$stop("Cannot use reporter function for models with a fixed tau^2 value."))

   if (!x$int.only)
      stop(mstyle$stop("Cannot currently use reporter function for models with moderators. This will be implemented shortly."))

   if (x$k == 1)
      stop(mstyle$stop("Cannot use reporter function when k = 1."))

   if (inherits(x, "robust.rma"))
      stop(mstyle$stop("Cannot use reporter function for objects of class \"robust.rma\"."))

   if (missing(digits)) {
      digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
   }

   format <- match.arg(format, c("html_document", "pdf_document", "word_document")) # , "bookdown::pdf_document2"))

   if (format == "pdf_document" && (Sys.which("pdflatex") == ""))
      warning(mstyle$warning("Cannot detect pdflatex executable. Rendering the pdf is likely to fail."), immediate.=TRUE)

   ### set/get directory for generating the report

   if (missing(dir)) {
      dir <- normalizePath(tempdir(), winslash="/")
      success <- file.exists(dir)
      if (!success)
         stop(mstyle$stop("No temporary directory available for creating the report."))
   } else {
      if (!is.character(dir))
         stop(mstyle$stop("Argument 'dir' must be a character string."))
      success <- file.exists(dir)
      if (!success)
         stop(mstyle$stop("Specified directory does not exist."))
   }

   if (verbose)
      message(mstyle$message(paste0("\nDirectory for generating the report is: ", dir, "\n")))

   ### copy references.bib and apa.csl files to directory for generating the report

   if (verbose)
      message(mstyle$message("Copying references.bib and apa.csl to report directory ..."))

   success <- file.copy(system.file("reporter", "references.bib", package = "metafor"), dir, overwrite=TRUE)

   if (!success)
      stop(mstyle$stop("Could not copy 'references.bib' file to report directory."))

   success <- file.copy(system.file("reporter", "apa.csl", package = "metafor"), dir, overwrite=TRUE)

   if (!success)
      stop(mstyle$stop("Could not copy 'apa.csl' file to report directory."))

   ### set default filenames

   object.name <- deparse(substitute(x))
   has.object.name <- TRUE

   if (grepl("rma(", object.name, fixed=TRUE) || grepl("rma.uni(", object.name, fixed=TRUE)) { # check for 'reporter(rma(yi, vi))' usage
      has.object.name <- FALSE
      object.name <- "res"
   }

   if (missing(filename)) {
      file.rmd <- paste0("report_", object.name, ".rmd")
      file.obj <- paste0("report_", object.name, ".rdata")
      file.tex <- paste0("report_", object.name, ".tex")
   } else {
      if (!is.character(filename))
         stop(mstyle$stop("Argument 'filename' must be a character string."))
      file.rmd <- paste0(filename, ".rmd")
      file.obj <- paste0(filename, ".rdata")
      file.tex <- paste0(filename, ".tex")
   }

   ### process forest argument

   if (missing(forest)) {
      args.forest <- ""
   } else {
      if (!is.character(args.forest))
         stop(mstyle$stop("Argument 'args.forest' must be a character string."))
      args.forest <- paste0(", ", forest)
   }

   ### process funnel argument

   if (missing(funnel)) {
      args.funnel <- ""
   } else {
      if (!is.character(args.funnel))
         stop(mstyle$stop("Argument 'args.funnel' must be a character string."))
      args.funnel <- paste0(", ", funnel)
   }

   ### save model object

   if (verbose)
      message(mstyle$message(paste0("Saving model object to ", file.obj, " ...")))

   success <- try(save(x, file=file.path(dir, file.obj)))

   if (inherits(success, "try-error"))
      stop(mstyle$stop("Could not save model object to report directory."))

   ### open rmd file connection

   if (verbose)
      message(mstyle$message(paste0("Creating ", file.rmd, " file ...")))

   con <- try(file(file.path(dir, file.rmd), "w"))

   if (inherits(con, "try-error"))
      stop(mstyle$stop("Could not create .rmd file in report directory."))

   ### get measure name

   measure <- tolower(.setlab(x$measure, transf.char=FALSE, atransf.char=FALSE, gentype=1))
   measure <- sub("observed outcome", "outcome", measure)
   measure <- sub("fisher's z", "Fisher r-to-z", measure)
   measure <- sub("yule", "Yule", measure)
   measure <- sub("freeman", "Freeman", measure)
   measure <- sub("tukey", "Tukey", measure)
   measure <- sub("log ratio of means", "response ratio", measure)

   ### model type

   model <- ifelse(x$method == "FE", ifelse(x$int.only, "FE", "MR"), ifelse(x$int.only, "RE", "ME"))
   model.name <- c(FE = "fixed-effects", MR = "(fixed-effects) meta-regression", RE = "random-effects", ME = "(mixed-effects) meta-regression")[model]

   ### get tau^2 estimator name and set reference

   tau2.method <- c(FE = "", HS = "Hunter-Schmidt", HE = "Hedges'", DL = "DerSimonian-Laird", GENQ = "generalized Q-statistic", GENQM = "(median-unbiased) generalized Q-statistic", SJ = "Sidik-Jonkman", ML = "maximum-likelihood", REML = "restricted maximum-likelihood", EB = "empirical Bayes", PM = "Paule-Mandel", PMM = "(median-unbiased) Paule-Mandel")[x$method]

   if (x$method == "HS" && model == "RE")
      tau2.ref <- "[@hunter1990; @viechtbauer2005]"
   if (x$method == "HS" && model == "ME")
      tau2.ref <- "[@hunter1990; @viechtbauer2015]"

   if (x$method == "HE" && model == "RE")
      tau2.ref <- "[@hedges1985]"
   if (x$method == "HE" && model == "ME")
      tau2.ref <- "[@hedges1992]"

   if (x$method == "DL" && model == "RE")
      tau2.ref <- "[@dersimonian1986]"
   if (x$method == "DL" && model == "ME")
      tau2.ref <- "[@raudenbush2009]"

   if (x$method == "GENQ" && model == "RE")
      tau2.ref <- "[@dersimonian2007]"
   if (x$method == "GENQ" && model == "ME")
      tau2.ref <- "[@jackson2014]"

   if (x$method == "GENQM" && model == "RE")
      tau2.ref <- "[@dersimonian2007]"
   if (x$method == "GENQM" && model == "ME")
      tau2.ref <- "[@jackson2014]"

   if (x$method == "SJ")
      tau2.ref <- "[@sidik2005]"

   if (x$method == "ML" && model == "RE")
      tau2.ref <- "[@hardy1996]"
   if (x$method == "ML" && model == "ME")
      tau2.ref <- "[@raudenbush2009]"

   if (x$method == "REML" && model == "RE")
      tau2.ref <- "[@viechtbauer2005]"
   if (x$method == "REML" && model == "ME")
      tau2.ref <- "[@raudenbush2009]"

   if (x$method == "EB" && model == "RE")
      tau2.ref <- "[@morris1983]"
   if (x$method == "EB" && model == "ME")
      tau2.ref <- "[@berkey1995]"

   if (x$method == "PM" && model == "RE")
      tau2.ref <- "[@paule1982]"
   if (x$method == "PM" && model == "ME")
      tau2.ref <- "[@viechtbauer2015]"

   if (x$method == "PMM" && model == "RE")
      tau2.ref <- "[@paule1982]"
   if (x$method == "PMM" && model == "ME")
      tau2.ref <- "[@viechtbauer2015]"

   ### Q-test reference

   if (is.element(model, c("FE", "RE"))) {
      qtest.ref <- "[@cochran1954]"
   } else {
      qtest.ref <- "[@hedges1983]"
   }

   ### CI level

   level <- 100 * (1-x$level)

   ### Bonferroni-corrected critical value for studentized residuals

   crit <- qnorm(x$level/(2*x$k), lower.tail=FALSE)

   ### get influence results

   infres <- influence(x)

   ### formating function for p-values

   fpval <- function(p, pdigits=digits[["pval"]])
      paste0("$p ", ifelse(p < 10^(-pdigits), paste0("< ", .fcf(10^(-pdigits), pdigits)), paste0("= ", .fcf(p, pdigits))), "$")
   # consider giving only 2 digits for p-value if p > .05 or p > .10

   #########################################################################

   ### yaml header

   header <- paste0("---\n")
   header <- paste0(header, "output: ", format, "\n")
   header <- paste0(header, "title: Analysis Report\n")
   header <- paste0(header, "author: Generated with the reporter() Function of the metafor Package\n")
   header <- paste0(header, "bibliography: references.bib\n")
   header <- paste0(header, "csl: apa.csl\n")
   header <- paste0(header, "date: \"`r format(Sys.time(), '%d %B, %Y')`\"\n")
   header <- paste0(header, "---\n")

   #########################################################################

   ### rsetup

   rsetup <- paste0("```{r, setup, include=FALSE}\n")
   rsetup <- paste0(rsetup, "library(metafor)\n")
   rsetup <- paste0(rsetup, "load('", file.path(dir, file.obj), "')\n")
   rsetup <- paste0(rsetup, "```")

   #########################################################################

   ### methods section

   methods <- "\n## Methods\n\n"

   if (x$measure != "GEN")
      methods <- paste0(methods, "The analysis was carried out using the ", measure, " as the outcome measure. ")

   methods <- paste0(methods, "A ", model.name, " model was fitted to the data. ")

   if (is.element(model, c("RE", "ME")))
      methods <- paste0(methods, "The amount of ", ifelse(x$int.only, "", "residual "), "heterogeneity (i.e., $\\tau^2$), was estimated using the ", tau2.method, " estimator ", tau2.ref, ". ")

   if (model == "FE")
      methods <- paste0(methods, "The $Q$-test for heterogeneity ", qtest.ref, " and the $I^2$ statistic [@higgins2002] are reported. ")

   if (model == "MR")
      methods <- paste0(methods, "The $Q$-test for residual heterogeneity ", qtest.ref, " is reported. ")

   if (model == "RE")
      methods <- paste0(methods, "In addition to the estimate of $\\tau^2$, the $Q$-test for heterogeneity ", qtest.ref, " and the $I^2$ statistic [@higgins2002] are reported. ")

   if (model == "ME")
      methods <- paste0(methods, "In addition to the estimate of $\\tau^2$, the $Q$-test for residual heterogeneity ", qtest.ref, " is reported. ")

   if (model == "RE")
      methods <- paste0(methods, "In case any amount of heterogeneity is detected (i.e., $\\hat{\\tau}^2 > 0$, regardless of the results of the $Q$-test), a credibility/prediction interval for the true outcomes is also provided [@riley2011]. ")

   if (x$test == "knha")
      methods <- paste0(methods, "Tests and confidence intervals were computed using the Knapp and Hartung method [@knapp2003]. ")

   methods <- paste0(methods, "Studentized residuals and Cook's distances are used to examine whether studies may be outliers and/or influential in the context of the model [@viechtbauer2010b]. ")

   #methods <- paste0(methods, "Studies with a studentized residual larger than $\\pm 1.96$ are considered potential outliers. ")
   methods <- paste0(methods, "Studies with a studentized residual larger than the $100 \\times (1 - ", x$level, "/(2 \\times k))$th percentile of a standard normal distribution are considered potential outliers (i.e., using a Bonferroni correction with two-sided $\\alpha = ", x$level, "$ for $k$ studies included in the meta-analysis). ") # $\\pm ", .fcf(crit, digits[["test"]]), "$ (

   #methods <- paste0(methods, "Studies with a Cook's distance larger than ", .fcf(qchisq(0.5, df=infres$m), digits[["test"]]), " (the 50th percentile of a $\\chi^2$-distribution with ", infres$m, " degree", ifelse(infres$m > 1, "s", ""), " of freedom) are considered to be influential. ")
   methods <- paste0(methods, "Studies with a Cook's distance larger than the median plus six times the interquartile range of the Cook's distances are considered to be influential.")

   methods <- if (footnotes) paste0(methods, "[^cook] ") else paste0(methods, " ")

   if (is.element(model, c("FE", "RE")))
      methods <- paste0(methods, "The rank correlation test [@begg1994] and the regression test [@sterne2005], using the standard error of the observed outcomes as predictor, are used to check for funnel plot asymmetry. ")

   if (is.element(model, c("MR", "ME")))
      methods <- paste0(methods, "The regression test [@sterne2005], using the standard error of the observed outcomes as predictor (in addition to the moderators already included in the model), is used to check for funnel plot asymmetry. ")

   methods <- paste0(methods, "The analysis was carried out using R (version ", getRversion(), ") [@rcore2018] and the **metafor** package (version ", x$version, ") [@viechtbauer2010a]. ")

   #########################################################################

   ### results section

   results <- "\n## Results\n\n"

   ### number of studies
   results <- paste0(results, "A total of $k=", x$k, "$ studies were included in the analysis. ")

   ### range of observed outcomes
   results <- paste0(results, "The observed ", measure, "s ranged from $", .fcf(min(x$yi), digits[["est"]]), "$ to $", .fcf(max(x$yi), digits[["est"]]), "$, ")

   ### percent positive/negative
   results <- paste0(results, "with the majority of estimates being ", ifelse(mean(x$yi > 0) > .50, "positive", "negative"), " (", ifelse(mean(x$yi > 0) > .50, round(100*mean(x$yi > 0)), round(100*mean(x$yi < 0))), "%). ")

   if (is.element(model, c("FE", "RE"))) {

      ### estimated average outcome with CI
      results <- paste0(results, "The estimated average ", measure, " based on the ", model.name, " model was ", ifelse(model == "FE", "$\\hat{\\theta} = ", "$\\hat{\\mu} = "), .fcf(c(x$beta), digits[["est"]]), "$ ")
      results <- paste0(results, "(", level, "% CI: $", .fcf(x$ci.lb, digits[["ci"]]), "$ to $", .fcf(x$ci.ub, digits[["ci"]]), "$). ")

      ### note: for some outcome measures (e.g., proportions), the test H0: mu/theta = 0 is not really relevant; maybe check for this
      results <- paste0(results, "Therefore, the average outcome ", ifelse(x$pval > 0.05, "did not differ", "differed"), " significantly from zero ($", ifelse(x$test == "z", "z", paste0("t(", x$k-1, ")")), " = ", .fcf(x$zval, digits[["test"]]), "$, ", fpval(x$pval), "). ")

      ### forest plot
      results <- paste0(results, "A forest plot showing the observed outcomes and the estimate based on the ", model.name, " model is shown in Figure 1.\n\n")
      if (is.element(format, c("pdf_document", "bookdown::pdf_document2")))
         results <- paste0(results, "```{r, forestplot, echo=FALSE, fig.align=\"center\", fig.cap=\"Forest plot showing the observed outcomes and the estimate of the ", model.name, " model\"")
      if (format == "html_document")
         results <- paste0(results, "```{r, forestplot, echo=FALSE, fig.align=\"center\", fig.cap=\"Figure 1: Forest plot showing the observed outcomes and the estimate of the ", model.name, " model\"")
      if (format == "word_document")
         results <- paste0(results, "```{r, forestplot, echo=FALSE, fig.cap=\"Figure 1: Forest plot showing the observed outcomes and the estimate of the ", model.name, " model\"")
      results <- paste0(results, ", dev.args=list(pointsize=9)}\npar(family=\"mono\")\ntmp <- metafor::forest(x, addcred=TRUE", args.forest, ")\ntext(tmp$xlim[1], x$k+2, \"Study\", pos=4, font=2, cex=tmp$cex)\ntext(tmp$xlim[2], x$k+2, \"Outcome [", level, "% CI]\", pos=2, font=2, cex=tmp$cex)\n```")

      results <- paste0(results, "\n\n")

      ### test for heterogeneity
      if (x$QEp > 0.10)
         results <- paste0(results, "According to the $Q$-test, there was no significant amount of heterogeneity in the true outcomes ")
      if (x$QEp > 0.05 && x$QEp <= 0.10)
         results <- paste0(results, "The $Q$-test for heterogeneity was not significant, but some heterogeneity may still be present in the true outcomes ")
      if (x$QEp <= 0.05)
         results <- paste0(results, "According to the $Q$-test, the true outcomes appear to be heterogeneous ")
      results <- paste0(results, "($Q(", x$k-1, ") = ", .fcf(x$QE, digits[["test"]]), "$, ", fpval(x$QEp))

      ### tau^2 estimate (only for RE models)
      if (model == "RE")
         results <- paste0(results, ", $\\hat{\\tau}^2 = ", .fcf(x$tau2, digits[["var"]]), "$")

      ### I^2 statistic
      results <- paste0(results, ", $I^2 = ", .fcf(x$I2, digits[["het"]]), "$%). ")

      ### for the RE model, when any amount of heterogeneity is detected, provide credibility/prediction interval and note whether the directionality of effects is consistent or not
      if (model == "RE" && x$tau2 > 0) {
         pred <- predict(x)
         results <- paste0(results, "A ", level, "% credibility/prediction interval for the true outcomes is given by $", .fcf(pred$cr.lb, digits[["ci"]]), "$ to $", .fcf(pred$cr.ub, digits[["ci"]]), "$. ")
         if (c(x$beta) > 0 && pred$cr.lb < 0)
            results <- paste0(results, "Hence, although the average outcome is estimated to be positive, in some studies the true outcome may in fact be negative.")
         if (c(x$beta) < 0 && pred$cr.ub > 0)
            results <- paste0(results, "Hence, although the average outcome is estimated to be negative, in some studies the true outcome may in fact be positive.")
         if ((c(x$beta) > 0 && pred$cr.lb > 0) || (c(x$beta) < 0 && pred$cr.ub < 0))
            results <- paste0(results, "Hence, even though there may be some heterogeneity, the true outcomes of the studies are generally in the same direction as the estimated average outcome.")
      }

      results <- paste0(results, "\n\n")

      ### check if some studies have very large weights relatively speaking
      largeweight <- weights(x)/100 >= 3 / x$k
      if (any(largeweight)) {
         if (sum(largeweight) == 1)
            results <- paste0(results, "One study (", names(largeweight)[largeweight], ") had a relatively large weight ")
         if (sum(largeweight) == 2)
            results <- paste0(results, "Two studies (", paste(names(largeweight)[largeweight], collapse="; "), ") had relatively large weights ")
         if (sum(largeweight) >= 3)
            results <- paste0(results, "Several studies (", paste(names(largeweight)[largeweight], collapse="; "), ") had relatively large weights ")
         results <- paste0(results, "compared to the rest of the studies (i.e., $\\mbox{weight} \\ge 3/k$, so a weight at least 3 times as large as having equal weights across studies). ")
      }

      ### check for outliers
      zi <- infres$inf$rstudent
      abszi <- abs(zi)
      results <- paste0(results, "An examination of the studentized residuals revealed that ")
      if (all(abszi < crit, na.rm=TRUE))
         results <- paste0(results, "none of the studies had a value larger than $\\pm ", .fcf(crit, digits[["test"]]), "$ and hence there was no indication of outliers ")
      if (sum(abszi >= crit, na.rm=TRUE) == 1)
         results <- paste0(results, "one study (", infres$inf$slab[abszi >= crit & !is.na(abszi)], ") had a value larger than $\\pm ", .fcf(crit, digits[["test"]]), "$ and may be a potential outlier ")
      if (sum(abszi >= crit, na.rm=TRUE) == 2)
         results <- paste0(results, "two studies (", paste(infres$inf$slab[abszi >= crit & !is.na(abszi)], collapse="; "), ") had values larger than $\\pm ", .fcf(crit, digits[["test"]]), "$ and may be potential outliers ")
      if (sum(abszi >= crit, na.rm=TRUE) >= 3)
         results <- paste0(results, "several studies (", paste(infres$inf$slab[abszi >= crit & !is.na(abszi)], collapse="; "), ") had values larger than $\\pm ", .fcf(crit, digits[["test"]]), "$ and may be potential outliers ")
      results <- paste0(results, "in the context of this model. ")

      ### check for influential cases
      #is.infl <- pchisq(infres$inf$cook.d, df=1) > .50
      is.infl <- infres$inf$cook.d > median(infres$inf$cook.d, na.rm=TRUE) + 6 * IQR(infres$inf$cook.d, na.rm=TRUE)
      results <- paste0(results, "According to the Cook's distances, ")
      if (all(!is.infl, na.rm=TRUE))
         results <- paste0(results, "none of the studies ")
      if (sum(is.infl, na.rm=TRUE) == 1)
         results <- paste0(results, "one study (", infres$inf$slab[is.infl & !is.na(abszi)], ") ")
      if (sum(is.infl, na.rm=TRUE) == 2)
         results <- paste0(results, "two studies (", paste(infres$inf$slab[is.infl & !is.na(abszi)], collapse="; "), ") ")
      if (sum(is.infl, na.rm=TRUE) >= 3)
         results <- paste0(results, "several studies (", paste(infres$inf$slab[is.infl & !is.na(abszi)], collapse="; "), ") ")
      results <- paste0(results, "could be considered to be overly influential.")

      results <- paste0(results, "\n\n")

      ### publication bias
      ranktest <- suppressWarnings(ranktest(x))
      regtest  <- regtest(x)
      results <- paste0(results, "A funnel plot of the estimates is shown in Figure 2. ")
      if (ranktest$pval > .05 && regtest$pval > .05) {
         results <- paste0(results, "Neither the rank correlation nor the regression test indicated any funnel plot asymmetry ")
         results <- paste0(results, "(", fpval(ranktest$pval), " and ", fpval(regtest$pval), ", respectively). ")
      }
      if (ranktest$pval <= .05 && regtest$pval <= .05) {
         results <- paste0(results, "Both the rank correlation and the regression test indicated potential funnel plot asymmetry ")
         results <- paste0(results, "(", fpval(ranktest$pval), " and ", fpval(regtest$pval), ", respectively). ")
      }
      if (ranktest$pval > .05 && regtest$pval <= .05)
         results <- paste0(results, "The regression test indicated funnel plot asymmetry (", fpval(regtest$pval), ") but not the rank correlation test (", fpval(ranktest$pval), "). ")
      if (ranktest$pval <= .05 && regtest$pval > .05)
         results <- paste0(results, "The rank correlation test indicated funnel plot asymmetry ($", fpval(ranktest$pval), ") but not the regression test (", fpval(regtest$pval), "). ")

      ### funnel plot
      if (is.element(format, c("pdf_document", "bookdown::pdf_document2")))
         results <- paste0(results, "\n\n```{r, funnelplot, echo=FALSE, fig.align=\"center\", fig.cap=\"Funnel plot\", dev.args=list(pointsize=9)}\nmetafor::funnel(x", args.funnel, ")\n```")
      if (format == "html_document")
         results <- paste0(results, "\n\n```{r, funnelplot, echo=FALSE, fig.align=\"center\", fig.cap=\"Figure 2: Funnel plot\", dev.args=list(pointsize=9)}\nmetafor::funnel(x", args.funnel, ")\n```")
      if (format == "word_document")
         results <- paste0(results, "\n\n```{r, funnelplot, echo=FALSE, fig.cap=\"Figure 2: Funnel plot\", dev.args=list(pointsize=9)}\nmetafor::funnel(x", args.funnel, ")\n```")

   }

   if (is.element(model, c("MR", "ME"))) {

      if (x$int.incl) {
         mods <- colnames(x$X)[-1]
         p <- x$p - 1
      } else {
         mods <- colnames(x$X)
         p <- x$p
      }

      results <- paste0(results, "The meta-regression model included ", p, " predictor", ifelse(p > 1, "s ", " "))
      if (p == 1)
         results <- paste0(results, "(i.e., '", mods, "').")
      if (p == 2)
         results <- paste0(results, "(i.e., '", mods[1], "' and '", mods[2], "').")
      if (p >= 3)
         results <- paste0(results, "(i.e., ", paste0("'", mods[-p], "'", collapse=", "), " and ", mods[p], ").")

   }

   # 95% CI for tau^2 and I^2
   # table for meta-regression model
   # links to help pages for functions used

   #########################################################################

   ### notes section

   notes <- "\n## Notes\n\n"

   notes <- paste0(notes, "This analysis report was dynamically generated ", ifelse(has.object.name, paste0("for model object '`", object.name, "`'"), ""), " with the `reporter()` function of the **metafor** package. ")

   call <- capture.output(x$call)
   call <- trimws(call, which="left")
   call <- paste(call, collapse="")

   notes <- paste0(notes, "The model call that was used to fit the model was '`", call, "`'. ")
   notes <- paste0(notes, "This report provides an illustration of how the results of the model can be reported, but is not a substitute for a careful examination of the results.")

   #########################################################################

   ### references section

   references <- "\n## References\n"

   #########################################################################

   if (footnotes) {
      fnotes <- ""
      fnotes <- paste0(fnotes, "[^cook]: This is a somewhat arbitrary rule, but tends to detect 'spikes' in a plot of the Cook's distances fairly accurately. A better rule may be implemented in the future.")
   }

   #########################################################################

   ### write sections to rmd file

   writeLines(header, con)
   writeLines(rsetup, con)
   writeLines(methods, con)
   writeLines(results, con)
   writeLines(notes, con)
   writeLines(references, con)

   if (footnotes)
      writeLines(fnotes, con)

   ### close rmd file connection

   close(con)

   ### render rmd file

   if (verbose)
      message(mstyle$message(paste0("Rendering ", file.rmd, " file ...")))

   if (verbose >= 2) {
      file.out <- rmarkdown::render(file.path(dir, file.rmd), output_format=format, quiet=ifelse(verbose <= 1, TRUE, FALSE))
   } else {
      file.out <- suppressWarnings(rmarkdown::render(file.path(dir, file.rmd), output_format=format, quiet=ifelse(verbose <= 1, TRUE, FALSE)))
   }

   if (verbose)
      message(mstyle$message(paste0("Generated ", file.out, " ...")))

   ### render() sometimes fails to delete the intermediate tex file, so in case this happens clean up
   ### see also: https://github.com/rstudio/rmarkdown/issues/1308

   if (file.exists(file.path(dir, file.tex)))
      unlink(file.path(dir, file.tex))

   ### try to open output file

   if (open) {

      if (verbose)
         message(mstyle$message(paste0("Opening report ...\n")))

      if (.Platform$OS.type == "windows") {
         shell.exec(file.out)
      } else {
         optb <- getOption("browser")
         if (is.function(optb)) {
            invisible(optb(file.out))
         } else {
            system(paste0(optb, " '", file.out, "'"))
         }
      }

   }

   invisible(file.out)

}
