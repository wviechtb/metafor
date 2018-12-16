print.rma.mv <- function(x, digits, showfit=FALSE, signif.stars=getOption("show.signif.stars"), signif.legend=signif.stars, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(x, "rma.mv"))
      stop(mstyle$stop("Argument 'x' must be an object of class \"rma.mv\"."))

   if (missing(digits))
      digits <- x$digits

   if (!exists(".rmspace"))
      cat("\n")

   cat(mstyle$section("Multivariate Meta-Analysis Model"))
   cat(mstyle$section(paste0(" (k = ", x$k, "; ")))
   cat(mstyle$section(paste0("method: ", x$method, ")")))

   if (showfit) {
      cat("\n")
      if (x$method == "REML") {
         fs <- c(formatC(round(x$fit.stats$REML, digits=digits), digits=digits, format="f"))
      } else {
         fs <- c(formatC(round(x$fit.stats$ML, digits=digits), digits=digits, format="f"))
      }
      names(fs) <- c("logLik", "Deviance", "AIC", "BIC", "AICc")
      cat("\n")
      tmp <- capture.output(print(fs, quote=FALSE, print.gap=2))
      .print.table(tmp, mstyle)
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

   cat(mstyle$section("Variance Components:"))

   right <- TRUE

   if (!x$withS && !x$withG && !x$withH) {
      cat(mstyle$text(" none"))
      cat("\n\n")
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
            rownames(vc) <- paste("sigma^2.", seq_along(x$sigma2), sep="")
         }
         tmp <- capture.output(print(vc, quote=FALSE, right=right, print.gap=2))
         .print.table(tmp, mstyle)
         cat("\n")

      }

      if (x$withG) {

         ### note: use g.nlevels.f[1] since the number of arms is based on all data (i.e., including NAs), but use
         ### g.nlevels[2] since the number of studies is based on what is actually available (i.e., excluding NAs)

         if (is.element(x$struct[1], c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","GEN"))) {
            inner <- trimws(paste0(strsplit(paste0(x$formulas[[1]], collapse=""), "|", fixed=TRUE)[[1]][1], collapse=""))
            if (nchar(inner) > 15)
               inner <- paste0(substr(inner, 1, 15), "[...]", collapse="")
         } else {
            inner <- x$g.names[1]
         }
         outer <- tail(x$g.names, 1)

         mng <- max(nchar(c(inner, outer)))

         cat(mstyle$text(paste0("outer factor: ", paste0(outer, paste(rep(" ", max(0,mng-nchar(outer))), collapse=""), collapse=""), " (nlvls = ", x$g.nlevels[2], ")")))
         cat("\n")
         if (is.element(x$struct[1], c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","GEN"))) {
            cat(mstyle$text(paste0("inner term:   ", paste0(inner, paste(rep(" ", max(0,mng-nchar(inner))), collapse=""), collapse=""), " (nlvls = ", x$g.nlevels.f[1], ")")))
         } else {
            cat(mstyle$text(paste0("inner factor: ", paste0(inner, paste(rep(" ", max(0,mng-nchar(inner))), collapse=""), collapse=""), " (nlvls = ", x$g.nlevels.f[1], ")")))
         }

         cat("\n\n")

         if (is.element(x$struct[1], c("CS","AR","CAR","ID","SPEXP","SPGAU","SPLIN","SPRAT","SPSPH"))) {

            vc <- cbind(tau2, tau, ifelse(x$vc.fix$tau2, "yes", "no"))
            vc <- rbind(vc, c(rho, "", ifelse(x$vc.fix$rho, "yes", "no")))
            colnames(vc) <- c("estim", "sqrt", "fixed")
            rownames(vc) <- c("tau^2    ", "rho")
            if (x$struct[1] == "ID")
               vc <- vc[1,,drop=FALSE]
            tmp <- capture.output(print(vc, quote=FALSE, right=right, print.gap=2))
            .print.table(tmp, mstyle)

         }

         if (is.element(x$struct[1], c("HCS","HAR","DIAG"))) {

            vc <- cbind(tau2, tau, x$g.levels.k, ifelse(x$vc.fix$tau2, "yes", "no"), x$g.levels.f[[1]])
            vc <- rbind(vc, c(rho, "", "", ifelse(x$vc.fix$rho, "yes", "no"), ""))
            colnames(vc) <- c("estim", "sqrt", "k.lvl", "fixed", "level")
            if (length(x$tau2) == 1) {
               rownames(vc) <- c("tau^2   ", "rho")
            } else {
               rownames(vc) <- c(paste("tau^2.", seq_along(x$tau2), "  ", sep=""), "rho")
            }
            if (x$struct[1] == "DIAG")
               vc <- vc[seq_along(tau2),,drop=FALSE]
            tmp <- capture.output(print(vc, quote=FALSE, right=right, print.gap=2))
            .print.table(tmp, mstyle)

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
               rownames(vc) <- paste("tau^2.", seq_along(x$g.levels.k), "  ", sep="")
            }
            tmp <- capture.output(print(vc, quote=FALSE, right=right, print.gap=2))
            .print.table(tmp, mstyle)
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
            colnames(vc) <- c(paste("rho.", abbreviate(x$g.levels.f[[1]]), sep=""), "", abbreviate(x$g.levels.f[[1]])) ### FIXME: x$g.levels.f[[1]] may be numeric, in which case a wrapping 'header' is not recognized
            rownames(vc) <- x$g.levels.f[[1]]
            tmp <- capture.output(print(vc, quote=FALSE, right=right, print.gap=2))
            .print.table(tmp, mstyle)

         }

         if (is.element(x$struct[1], c("GEN"))) {

            vc <- cbind(tau2, tau, ifelse(x$vc.fix$tau2, "yes", "no"), "")
            colnames(vc) <- c("estim", "sqrt", "fixed", "rho:")
            rownames(vc) <- x$g.names[-length(x$g.names)]

            G.info <- formatC(cov2cor(x$G), digits=digits, format="f")
            diag(G.info) <- "-"
            G.info[upper.tri(G.info)] <- ifelse(x$vc.fix$rho, "yes", "no")
            colnames(G.info) <- abbreviate(x$g.names[-length(x$g.names)])
            vc <- cbind(vc, G.info)
            tmp <- capture.output(print(vc, quote=FALSE, right=right, print.gap=2))
            .print.table(tmp, mstyle)

         }

         cat("\n")

      }

      if (x$withH) {

         ### note: use h.nlevels.f[1] since the number of arms is based on all data (i.e., including NAs), but use
         ### h.nlevels[2] since the number of studies is based on what is actually available (i.e., excluding NAs)

         if (is.element(x$struct[2], c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","GEN"))) {
            inner <- trimws(paste0(strsplit(paste0(x$formulas[[2]], collapse=""), "|", fixed=TRUE)[[1]][1], collapse=""))
            if (nchar(inner) > 15)
               inner <- paste0(substr(inner, 1, 15), "[...]", collapse="")
         } else {
            inner <- x$h.names[1]
         }
         outer <- tail(x$h.names, 1)

         mng <- max(nchar(c(inner, outer)))

         cat(mstyle$text(paste0("outer factor: ", paste0(outer, paste(rep(" ", max(0,mng-nchar(outer))), collapse=""), collapse=""), " (nlvls = ", x$h.nlevels[2], ")")))
         cat("\n")
         if (is.element(x$struct[2], c("SPEXP","SPGAU","SPLIN","SPRAT","SPSPH","GEN"))) {
            cat(mstyle$text(paste0("inner term:   ", paste0(inner, paste(rep(" ", max(0,mng-nchar(inner))), collapse=""), collapse=""), " (nlvls = ", x$h.nlevels.f[1], ")")))
         } else {
            cat(mstyle$text(paste0("inner factor: ", paste0(inner, paste(rep(" ", max(0,mng-nchar(inner))), collapse=""), collapse=""), " (nlvls = ", x$h.nlevels.f[1], ")")))
         }

         cat("\n\n")

         if (is.element(x$struct[2], c("CS","AR","CAR","ID","SPEXP","SPGAU","SPLIN","SPRAT","SPSPH"))) {

            vc <- cbind(gamma2, gamma, ifelse(x$vc.fix$gamma2, "yes", "no"))
            vc <- rbind(vc, c(phi, "", ifelse(x$vc.fix$phi, "yes", "no")))
            colnames(vc) <- c("estim", "sqrt", "fixed")
            rownames(vc) <- c("gamma^2  ", "phi")
            if (x$struct[2] == "ID")
               vc <- vc[1,,drop=FALSE]
            tmp <- capture.output(print(vc, quote=FALSE, right=right, print.gap=2))
            .print.table(tmp, mstyle)

         }

         if (is.element(x$struct[2], c("HCS","HAR","DIAG"))) {

            vc <- cbind(gamma2, gamma, x$h.levels.k, ifelse(x$vc.fix$gamma2, "yes", "no"), x$h.levels.f[[1]])
            vc <- rbind(vc, c(phi, "", "", ifelse(x$vc.fix$phi, "yes", "no"), ""))
            colnames(vc) <- c("estim", "sqrt", "k.lvl", "fixed", "level")
            if (length(x$gamma2) == 1) {
               rownames(vc) <- c("gamma^2 ", "phi")
            } else {
               rownames(vc) <- c(paste("gamma^2.", seq_along(x$gamma2), "  ", sep=""), "phi")
            }
            if (x$struct[2] == "DIAG")
               vc <- vc[seq_along(gamma2),,drop=FALSE]
            tmp <- capture.output(print(vc, quote=FALSE, right=right, print.gap=2))
            .print.table(tmp, mstyle)

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
               rownames(vc) <- paste("gamma^2.", seq_along(x$h.levels.k), "  ", sep="")
            }
            tmp <- capture.output(print(vc, quote=FALSE, right=right, print.gap=2))
            .print.table(tmp, mstyle)
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
            colnames(vc) <- c(paste("phi.", abbreviate(x$h.levels.f[[1]]), sep=""), "", abbreviate(x$h.levels.f[[1]])) ### FIXME: x$h.levels.f[[1]] may be numeric, in which case a wrapping 'header' is not recognized
            rownames(vc) <- x$h.levels.f[[1]]
            tmp <- capture.output(print(vc, quote=FALSE, right=right, print.gap=2))
            .print.table(tmp, mstyle)

         }

         if (is.element(x$struct[2], c("GEN"))) {

            vc <- cbind(gamma2, gamma, ifelse(x$vc.fix$gamma2, "yes", "no"), "")
            colnames(vc) <- c("estim", "sqrt", "fixed", "phi:")
            rownames(vc) <- x$h.names[-length(x$h.names)]

            H.info <- formatC(cov2cor(x$H), digits=digits, format="f")
            diag(H.info) <- "-"
            H.info[upper.tri(H.info)] <- ifelse(x$vc.fix$phi, "yes", "no")
            colnames(H.info) <- abbreviate(x$h.names[-length(x$h.names)])
            vc <- cbind(vc, H.info)
            tmp <- capture.output(print(vc, quote=FALSE, right=right, print.gap=2))
            .print.table(tmp, mstyle)

         }

         cat("\n")

      }

   }

   if (!is.na(x$QE)) {
      if (x$int.only) {
         cat(mstyle$section("Test for Heterogeneity:"))
         cat("\n")
         cat(mstyle$result(paste0("Q(df = ", x$k-x$p, ") = ", formatC(x$QE, digits=digits, format="f"), ", p-val ", .pval(x$QEp, digits=digits, showeq=TRUE, sep=" "))))
      } else {
         cat(mstyle$section("Test for Residual Heterogeneity:"))
         cat("\n")
         cat(mstyle$result(paste0("QE(df = ", x$k-x$p, ") = ", formatC(x$QE, digits=digits, format="f"), ", p-val ", .pval(x$QEp, digits=digits, showeq=TRUE, sep=" "))))
      }
      cat("\n\n")
   }

   if (x$p > 1 && !is.na(x$QM)) {
      cat(mstyle$section(paste0("Test of Moderators (coefficient", ifelse(x$m == 1, " ", "s "), .format.btt(x$btt),"):")))
      cat("\n")
      if (is.element(x$test, c("t"))) {
         cat(mstyle$result(paste0("F(df1 = ", x$m, ", df2 = ", x$dfs, ") = ", formatC(x$QM, digits=digits, format="f"), ", p-val ", .pval(x$QMp, digits=digits, showeq=TRUE, sep=" "))))
      } else {
         cat(mstyle$result(paste0("QM(df = ", x$m, ") = ", formatC(x$QM, digits=digits, format="f"), ", p-val ", .pval(x$QMp, digits=digits, showeq=TRUE, sep=" "))))
      }
      cat("\n\n")
   }

   res.table <- cbind(estimate=c(x$beta), se=x$se, zval=x$zval, pval=x$pval, ci.lb=x$ci.lb, ci.ub=x$ci.ub)
   rownames(res.table) <- rownames(x$beta)
   if (is.element(x$test, c("t")))
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

   cat(mstyle$section("Model Results:"))
   cat("\n\n")
   if (x$int.only) {
      tmp <- capture.output(.print.vector(res.table))
   } else {
      tmp <- capture.output(print(res.table, quote=FALSE, right=TRUE, print.gap=2))
   }
   .print.table(tmp, mstyle)

   if (signif.legend) {
      cat("\n")
      cat(mstyle$legend("---\nSignif. codes: "), mstyle$legend(attr(signif, "legend")))
      cat("\n")
   }

   if (!exists(".rmspace"))
      cat("\n")

   invisible()

}
