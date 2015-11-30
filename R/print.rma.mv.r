print.rma.mv <- function(x, digits, showfit=FALSE, signif.stars=getOption("show.signif.stars"), signif.legend=signif.stars, ...) {

   if (!is.element("rma.mv", class(x)))
      stop("Argument 'x' must be an object of class \"rma.mv\".")

   if (missing(digits))
      digits <- x$digits

   cat("\n")

   cat("Multivariate Meta-Analysis Model (k = ", x$k, "; ", sep="")
   cat("method: ", x$method, ")", sep="")

   if (showfit) {
      cat("\n")
      if (x$method == "REML") {
         fs <- c(formatC(round(x$fit.stats$REML, digits=digits), digits=digits, format="f"))
      } else {
         fs <- c(formatC(round(x$fit.stats$ML, digits=digits), digits=digits, format="f"))
      }
      names(fs) <- c("logLik", "Deviance", "AIC", "BIC", "AICc")
      cat("\n")
      print(fs, quote=FALSE, print.gap=2)
      cat("\n")
   } else {
      cat("\n\n")
   }

   sigma2 <- formatC(x$sigma2, digits=digits, format="f")
   tau2   <- formatC(x$tau2,   digits=digits, format="f")
   rho    <- formatC(x$rho,    digits=digits, format="f")
   gamma2 <- formatC(x$gamma2, digits=digits, format="f")
   phi    <- formatC(x$phi,    digits=digits, format="f")
   sigma  <- formatC(sqrt(x$sigma2), digits=digits, format="f")
   tau    <- formatC(sqrt(x$tau2),   digits=digits, format="f")
   gamma  <- formatC(sqrt(x$gamma2), digits=digits, format="f")

   cat("Variance Components: ")

   right <- TRUE

   if (!x$withS && !x$withG && !x$withH) {
      cat("none\n\n")
   } else {
      cat("\n\n")

      if (x$withS) {

         vc <- cbind(estim=sigma2, sqrt=sigma, nlvls=x$s.nlevels, fixed=ifelse(x$vc.fix$sigma2, "yes", "no"), factor=x$s.names, R=ifelse(x$Rfix, "yes", "no"))
         colnames(vc) <- c("estim", "sqrt", "nlvls", "fixed", "factor", "R")
         if (!x$withR)
            vc <- vc[,-6,drop=FALSE]
         if (length(x$sigma2) == 1) {
            rownames(vc) <- "sigma^2  "
         } else {
            rownames(vc) <- paste("sigma^2.", 1:length(x$sigma2), sep="")
         }
         print(vc, quote=FALSE, right=right, print.gap=2)
         cat("\n")

      }

      if (x$withG) {

         ### note: use g.nlevels.f[1] since the number of arms is based on all data (i.e., including NAs), but use
         ### g.nlevels[2] since the number of studies is based on what is actually available (i.e., excluding NAs)

         mng <- max(nchar(x$g.names))

         cat("outer factor: ", paste0(x$g.names[2], paste(rep(" ", max(0,mng-nchar(x$g.names[2]))), collapse=""), collapse=""), " (nlvls = ", x$g.nlevels[2], ")\n", sep="")
         cat("inner factor: ", paste0(x$g.names[1], paste(rep(" ", max(0,mng-nchar(x$g.names[1]))), collapse=""), collapse=""), " (nlvls = ", x$g.nlevels.f[1], ")\n", sep="")

         cat("\n")

         if (is.element(x$struct[1], c("CS","AR","ID"))) {

            vc <- cbind(tau2, tau, ifelse(x$vc.fix$tau2, "yes", "no"))
            vc <- rbind(vc, c(rho, "", ifelse(x$vc.fix$rho, "yes", "no")))
            colnames(vc) <- c("estim", "sqrt", "fixed")
            rownames(vc) <- c("tau^2    ", "rho")
            if (x$struct[1] == "ID")
               vc <- vc[1,,drop=FALSE]
            print(vc, quote=FALSE, right=right, print.gap=2)

         }

         if (is.element(x$struct[1], c("HCS","HAR","DIAG"))) {

            vc <- cbind(tau2, tau, x$g.levels.k, ifelse(x$vc.fix$tau2, "yes", "no"), x$g.levels.f[[1]])
            vc <- rbind(vc, c(rho, "", "", ifelse(x$vc.fix$rho, "yes", "no"), ""))
            colnames(vc) <- c("estim", "sqrt", "k.lvl", "fixed", "level")
            if (length(x$tau2) == 1) {
               rownames(vc) <- c("tau^2   ", "rho")
            } else {
               rownames(vc) <- c(paste("tau^2.", 1:length(x$tau2), "  ", sep=""), "rho")
            }
            if (x$struct[1] == "DIAG")
               vc <- vc[1:length(tau2),,drop=FALSE]
            print(vc, quote=FALSE, right=right, print.gap=2)

         }

         if (is.element(x$struct[1], c("UN","UNHO"))) {

            if (x$struct[1] == "UN") {
               vc <- cbind(tau2, tau, x$g.levels.k, ifelse(x$vc.fix$tau2, "yes", "no"), x$g.levels.f[[1]])
            } else {
               vc <- cbind(rep(tau2, length(x$g.levels.k)), rep(tau, length(x$g.levels.k)), x$g.levels.k, ifelse(rep(x$vc.fix$tau2,length(x$g.levels.k)), "yes", "no"), x$g.levels.f[[1]])
            }
            colnames(vc) <- c("estim", "sqrt", "k.lvl", "fixed", "level")
            if (length(x$g.levels.k) == 1) {
               rownames(vc) <- c("tau^2")
            } else {
               rownames(vc) <- paste("tau^2.", 1:length(x$g.levels.k), "  ", sep="")
            }
            print(vc, quote=FALSE, right=right, print.gap=2)
            cat("\n")

            if (length(x$rho) == 1) {
               G <- matrix(NA_real_, nrow=2, ncol=2)
            } else {
               G <- matrix(NA_real_, nrow=x$g.nlevels.f[1], ncol=x$g.nlevels.f[1])
            }
            G[upper.tri(G)] <- rho
            G[lower.tri(G)] <- t(G)[lower.tri(G)]
            diag(G) <- 1
            #G[upper.tri(G)] <- ""

            if (length(x$rho) == 1) {
               G.info <- matrix(NA_real_, nrow=2, ncol=2)
            } else {
               G.info <- matrix(NA_real_, nrow=x$g.nlevels.f[1], ncol=x$g.nlevels.f[1])
            }
            G.info[upper.tri(G.info)] <- x$g.levels.comb.k
            G.info[lower.tri(G.info)] <- t(G.info)[lower.tri(G.info)]
            G.info[upper.tri(G.info)] <- ifelse(x$vc.fix$rho, "yes", "no")
            diag(G.info) <- "-"

            vc <- cbind(G, "", G.info)
            colnames(vc) <- c(paste("rho.", abbreviate(x$g.levels.f[[1]]), sep=""), "", abbreviate(x$g.levels.f[[1]]))
            rownames(vc) <- x$g.levels.f[[1]]
            print(vc, quote=FALSE, right=right, print.gap=2)

         }

         cat("\n")

      }

      if (x$withH) {

         ### note: use h.nlevels.f[1] since the number of arms is based on all data (i.e., including NAs), but use
         ### h.nlevels[2] since the number of studies is based on what is actually available (i.e., excluding NAs)

         mng <- max(nchar(x$h.names))

         cat("outer factor: ", paste0(x$h.names[2], paste(rep(" ", max(0,mng-nchar(x$h.names[2]))), collapse=""), collapse=""), " (nlvls = ", x$h.nlevels[2], ")\n", sep="")
         cat("inner factor: ", paste0(x$h.names[1], paste(rep(" ", max(0,mng-nchar(x$h.names[1]))), collapse=""), collapse=""), " (nlvls = ", x$h.nlevels.f[1], ")\n", sep="")

         cat("\n")

         if (is.element(x$struct[2], c("CS","AR"))) {

            vc <- cbind(gamma2, gamma, ifelse(x$vc.fix$gamma2, "yes", "no"))
            vc <- rbind(vc, c(phi, "", ifelse(x$vc.fix$phi, "yes", "no")))
            colnames(vc) <- c("estim", "sqrt", "fixed")
            rownames(vc) <- c("gamma^2  ", "phi")
            if (x$struct[2] == "ID")
               vc <- vc[1,,drop=FALSE]
            print(vc, quote=FALSE, right=right, print.gap=2)

         }

         if (is.element(x$struct[2], c("HCS","HAR"))) {

            vc <- cbind(gamma2, gamma, x$h.levels.k, ifelse(x$vc.fix$gamma2, "yes", "no"), x$h.levels.f[[1]])
            vc <- rbind(vc, c(phi, "", "", ifelse(x$vc.fix$phi, "yes", "no"), ""))
            colnames(vc) <- c("estim", "sqrt", "k.lvl", "fixed", "level")
            if (length(x$gamma2) == 1) {
               rownames(vc) <- c("gamma^2 ", "rho")
            } else {
               rownames(vc) <- c(paste("gamma^2.", 1:length(x$gamma2), "  ", sep=""), "phi")
            }
            if (x$struct[2] == "DIAG")
               vc <- vc[1:length(gamma2),,drop=FALSE]
            print(vc, quote=FALSE, right=right, print.gap=2)

         }

         if (is.element(x$struct[2], c("UN","UNHO"))) {

            if (x$struct[2] == "UN") {
               vc <- cbind(gamma2, gamma, x$h.levels.k, ifelse(x$vc.fix$gamma2, "yes", "no"), x$h.levels.f[[1]])
            } else {
               vc <- cbind(rep(gamma2, length(x$h.levels.k)), rep(gamma, length(x$h.levels.k)), x$h.levels.k, ifelse(rep(x$vc.fix$gamma2,length(x$h.levels.k)), "yes", "no"), x$h.levels.f[[1]])
            }
            colnames(vc) <- c("estim", "sqrt", "k.lvl", "fixed", "level")
            if (length(x$h.levels.k) == 1) {
               rownames(vc) <- c("gamma^2")
            } else {
               rownames(vc) <- paste("gamma^2.", 1:length(x$h.levels.k), "  ", sep="")
            }
            print(vc, quote=FALSE, right=right, print.gap=2)
            cat("\n")

            if (length(x$phi) == 1) {
               H <- matrix(NA_real_, nrow=2, ncol=2)
            } else {
               H <- matrix(NA_real_, nrow=x$h.nlevels.f[1], ncol=x$h.nlevels.f[1])
            }
            H[upper.tri(H)] <- phi
            H[lower.tri(H)] <- t(H)[lower.tri(H)]
            diag(H) <- 1
            #H[upper.tri(H)] <- ""

            if (length(x$phi) == 1) {
               H.info <- matrix(NA_real_, nrow=2, ncol=2)
            } else {
               H.info <- matrix(NA_real_, nrow=x$h.nlevels.f[1], ncol=x$h.nlevels.f[1])
            }
            H.info[upper.tri(H.info)] <- x$h.levels.comb.k
            H.info[lower.tri(H.info)] <- t(H.info)[lower.tri(H.info)]
            H.info[upper.tri(H.info)] <- ifelse(x$vc.fix$phi, "yes", "no")
            diag(H.info) <- "-"

            vc <- cbind(H, "", H.info)
            colnames(vc) <- c(paste("phi.", abbreviate(x$h.levels.f[[1]]), sep=""), "", abbreviate(x$h.levels.f[[1]]))
            rownames(vc) <- x$h.levels.f[[1]]
            print(vc, quote=FALSE, right=right, print.gap=2)

         }

         cat("\n")

      }

   }

   if (!is.na(x$QE)) {
      if (x$int.only) {
         cat("Test for Heterogeneity: \n")
         cat("Q(df = ", x$k-x$p, ") = ", formatC(x$QE, digits=digits, format="f"), ", p-val ", .pval(x$QEp, digits=digits, showeq=TRUE, sep=" "), "\n\n", sep="")
      } else {
         cat("Test for Residual Heterogeneity: \n")
         cat("QE(df = ", x$k-x$p, ") = ", formatC(x$QE, digits=digits, format="f"), ", p-val ", .pval(x$QEp, digits=digits, showeq=TRUE, sep=" "), "\n\n", sep="")
      }
   }

   if (x$p > 1) {
      cat("Test of Moderators (coefficient(s) ", paste(x$btt, collapse=","),"): \n", sep="")
      if (x$knha) {
         cat("F(df1 = ", x$m, ", df2 = ", x$dfs, ") = ", formatC(x$QM, digits=digits, format="f"), ", p-val ", .pval(x$QMp, digits=digits, showeq=TRUE, sep=" "), "\n\n", sep="")
      } else {
         cat("QM(df = ", x$m, ") = ", formatC(x$QM, digits=digits, format="f"), ", p-val ", .pval(x$QMp, digits=digits, showeq=TRUE, sep=" "), "\n\n", sep="")
      }
   }

   res.table <- cbind(estimate=c(x$b), se=x$se, zval=x$zval, pval=x$pval, ci.lb=x$ci.lb, ci.ub=x$ci.ub)
   rownames(res.table) <- rownames(x$b)
   if (x$knha)
      colnames(res.table)[3] <- "tval"
   signif <- symnum(x$pval, corr=FALSE, na=FALSE, cutpoints=c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))
   if (signif.stars) {
      res.table <- cbind(formatC(res.table, digits=digits, format="f"), signif)
      colnames(res.table)[7] <- ""
   } else {
      res.table <- formatC(res.table, digits=digits, format="f")
   }
   res.table[,4] <- .pval(x$pval, digits=digits)

   if (x$int.only)
      res.table <- res.table[1,]

   cat("Model Results:\n\n")
   if (x$int.only) {
      print(res.table, quote=FALSE, right=TRUE)
   } else {
      print(res.table, quote=FALSE, right=TRUE, print.gap=2)
   }
   cat("\n")
   if (signif.legend)
      cat("---\nSignif. codes: ", attr(signif, "legend"), "\n\n")

   invisible()

}
