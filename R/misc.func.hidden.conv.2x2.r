############################################################################

.rec2x2diag <- function(sens, spec, ppv, ni, round=TRUE) {

   #eps <- .Machine$double.eps
   #eps <- 0

   #ai <- ifelse(sens <= 0 + eps | ppv <= 0 + eps, 0, ni * sens * ppv * (1 - spec) / (sens * (1 - ppv) + ppv * (1 - spec)))
   #bi <- ifelse(spec >= 1 - eps | ppv >= 1 - eps, 0, ni * sens * (1 - ppv) * (1 - spec) / (sens * (1 - ppv) + ppv * (1 - spec)))
   #ci <- ifelse(sens >= 1 - eps | npv >= 1 - eps, 0, ni * (1 - sens) * ppv * (1 - spec) / (sens * (1 - ppv) + ppv * (1 - spec)))
   #di <- ifelse(spec <= 0 + eps | npv <= 0 + eps, 0, ni * sens * spec * (1 - ppv) / (sens * (1 - ppv) + ppv * (1 - spec)))

   ai <- ni * sens * ppv * (1 - spec) / (sens * (1 - ppv) + ppv * (1 - spec))
   bi <- ni * sens * (1 - ppv) * (1 - spec) / (sens * (1 - ppv) + ppv * (1 - spec))
   ci <- ni * (1 - sens) * ppv * (1 - spec) / (sens * (1 - ppv) + ppv * (1 - spec))
   di <- ni * sens * spec * (1 - ppv) / (sens * (1 - ppv) + ppv * (1 - spec))

   out <- c(ai, bi, ci, di)

   if (round)
      out <- round(out)

   return(out)

}

.funconv2x2diag <- function(par, obs, ni, round=FALSE) {

   if (par[1] * (1 - par[3]) + par[3] * (1 - par[2]) <= .Machine$double.eps)
      return(10)

   rec <- .rec2x2diag(sens=par["sens"], spec=par["spec"], ppv=par["ppv"], ni=ni, round=FALSE)

   if (round)
      rec <- round(rec)

   ni.rec <- sum(rec)
   sens.rec <- rec[1] / (rec[1] + rec[3])
   spec.rec <- rec[4] / (rec[2] + rec[4])
   ppv.rec  <- rec[1] / (rec[1] + rec[2])
   npv.rec  <- rec[4] / (rec[3] + rec[4])

   #loss <- (obs["sens"] - sens.rec)^2 + (obs["spec"] - spec.rec)^2 + (obs["ppv"] - ppv.rec)^2 + (obs["npv"] - npv.rec)^2
   loss <- (obs["sens"] - sens.rec)^2 + (obs["spec"] - spec.rec)^2 + (obs["ppv"] - ppv.rec)^2 + (obs["npv"] - npv.rec)^2 + (ni - ni.rec)^2
   #loss <- (qlogis(obs["sens"]) - qlogis(sens.rec))^2 + (qlogis(obs["spec"]) - qlogis(spec.rec))^2 + (qlogis(obs["ppv"]) - qlogis(ppv.rec))^2 + (qlogis(obs["npv"]) - qlogis(npv.rec))^2 + (ni - ni.rec)^2

   return(loss)

}

.largestremaindermethod <- function(x, n) {

   xfloor <- floor(x)
   remainder <- n - sum(xfloor)
   x.int <- xfloor
   if (isTRUE(remainder > 0)) {
      frac.part <- x - xfloor
      idx <- order(-frac.part)
      x.int[idx[seq_len(remainder)]] <- x.int[idx[seq_len(remainder)]] + 1
   }
   return(x.int)

}

############################################################################

.rec2x2diagold <- function(sens, spec, ppv, npv, ni, round=TRUE) {

   # using sens, spec, and ppv (without npv) (this is also used if all four are available)
   if (!is.na(sens) && !is.na(spec) && !is.na(ppv)) {
      ai <- ni * sens * ppv * (1 - spec) / (sens * (1 - ppv) + ppv * (1 - spec))
      bi <- ni * sens * (1 - ppv) * (1 - spec) / (sens * (1 - ppv) + ppv * (1 - spec))
      ci <- ni * (1 - sens) * ppv * (1 - spec) / (sens * (1 - ppv) + ppv * (1 - spec))
      di <- ni * sens * spec * (1 - ppv) / (sens * (1 - ppv) + ppv * (1 - spec))
   }

   # using sens, spec, and npv (without ppv)
   if (!is.na(sens) && !is.na(spec) && is.na(ppv) && !is.na(npv)) {
      ai <- ni * sens * spec * (1 - npv) / (spec * (1 - npv) + (1 - sens) * npv)
      bi <- ni * (1 - sens) * (1 - spec) * npv / (spec * (1 - npv) + (1 - sens) * npv)
      ci <- ni * (1 - sens) * spec * (1 - npv) / (spec * (1 - npv) + (1 - sens) * npv)
      di <- ni * (1 - sens) * spec * npv / (spec * (1 - npv) + (1 - sens) * npv)
   }

   # using sens, ppv, and npv (without spec)
   if (!is.na(sens) && is.na(spec) && !is.na(ppv) && !is.na(npv)) {
      ai <- ni * sens * ppv * (1 - npv) / (sens * (1 - npv) + ppv * (1 - sens))
      bi <- ni * sens * (1 - ppv) * (1 - npv) / (sens * (1 - npv) + ppv * (1 - sens))
      ci <- ni * (1 - sens) * ppv * (1 - npv) / (sens * (1 - npv) + ppv * (1 - sens))
      di <- ni * (1 - sens) * ppv * npv / (sens * (1 - npv) + ppv * (1 - sens))
   }

   # using spec, ppv, and npv (without sens)
   if (is.na(sens) && !is.na(spec) && !is.na(ppv) && !is.na(npv)) {
      ai <- ni * (1 - spec) * ppv * npv / (spec * (1 - ppv) + (1 - spec) * npv)
      bi <- ni * (1 - spec) * (1 - ppv) * npv / (spec * (1 - ppv) + (1 - spec) * npv)
      ci <- ni * spec * (1 - ppv) * (1 - npv) / (spec * (1 - ppv) + (1 - spec) * npv)
      di <- ni * spec * (1 - ppv) * npv / (spec * (1 - ppv) + (1 - spec) * npv)
   }

   out <- c(ai, bi, ci, di)

   if (round)
      out <- round(out)

   return(out)

}

.funconv2x2diagold <- function(par, obs, ni, round=FALSE) {

   sens <- obs["sens"]
   spec <- obs["spec"]
   ppv  <- obs["ppv"]
   npv  <- obs["npv"]

   if (!is.na(sens) && !is.na(spec) && !is.na(ppv) && !is.na(npv))
      rec <- .rec2x2diagold(sens=par["sens"], spec=par["spec"], ppv=par["ppv"], npv=par["npv"], ni=ni, round=FALSE)
   if (!is.na(sens) && !is.na(spec) && !is.na(ppv) && is.na(npv))
      rec <- .rec2x2diagold(sens=par["sens"], spec=par["spec"], ppv=par["ppv"], npv=NA, ni=ni, round=FALSE)
   if (!is.na(sens) && !is.na(spec) && is.na(ppv) && !is.na(npv))
      rec <- .rec2x2diagold(sens=par["sens"], spec=par["spec"], ppv=NA, npv=par["npv"], ni=ni, round=FALSE)
   if (!is.na(sens) && is.na(spec) && !is.na(ppv) && !is.na(npv))
      rec <- .rec2x2diagold(sens=par["sens"], spec=NA, ppv=par["ppv"], npv=par["npv"], ni=ni, round=FALSE)
   if (is.na(sens) && !is.na(spec) && !is.na(ppv) && !is.na(npv))
      rec <- .rec2x2diagold(sens=NA, spec=par["spec"], ppv=par["ppv"], npv=par["npv"], ni=ni, round=FALSE)

   if (round)
      rec <- round(rec)

   ni.rec <- sum(rec)
   sens.rec <- rec[1] / (rec[1] + rec[3])
   spec.rec <- rec[4] / (rec[2] + rec[4])
   ppv.rec  <- rec[1] / (rec[1] + rec[2])
   npv.rec  <- rec[4] / (rec[3] + rec[4])

   loss.sens <- (sens - sens.rec)^2
   loss.spec <- (spec - spec.rec)^2
   loss.ppv  <- (ppv  - ppv.rec)^2
   loss.npv  <- (npv  - npv.rec)^2
   loss.n    <- (ni   - ni.rec)^2

   loss <- sum(loss.sens, loss.spec, loss.ppv, loss.npv, loss.n, na.rm=TRUE)

   return(loss)

}

############################################################################
