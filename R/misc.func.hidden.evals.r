############################################################################

### to register getfit method for 'rma.uni' and 'rma.mv' objects: eval(metafor:::.glmulti)

.glmulti <- str2expression("

if (!(\"glmulti\" %in% .packages()))
   stop(\"Must load the 'glmulti' package first to use this code.\")

setOldClass(\"rma.uni\")

setMethod(\"getfit\", \"rma.uni\", function(object, ...) {
   if (object$test==\"z\") {
      cbind(estimate=coef(object), se=sqrt(diag(vcov(object))), df=Inf)
   } else {
      cbind(estimate=coef(object), se=sqrt(diag(vcov(object))), df=object$k-object$p)
   }
})

setOldClass(\"rma.mv\")

setMethod(\"getfit\", \"rma.mv\", function(object, ...) {
   if (object$test==\"z\") {
      cbind(estimate=coef(object), se=sqrt(diag(vcov(object))), df=Inf)
   } else {
      cbind(estimate=coef(object), se=sqrt(diag(vcov(object))), df=object$k-object$p)
   }
})

setOldClass(\"rma.glmm\")

setMethod(\"getfit\", \"rma.glmm\", function(object, ...) {
   if (object$test==\"z\") {
      cbind(estimate=coef(object), se=sqrt(diag(vcov(object))), df=Inf)
   } else {
      cbind(estimate=coef(object), se=sqrt(diag(vcov(object))), df=object$k-object$p)
   }
})

")

### helper functions to make MuMIn work together with metafor: eval(metafor:::.MuMIn)

.MuMIn <- str2expression("

makeArgs.rma <- function (obj, termNames, comb, opt, ...) {
   ret <- MuMIn:::makeArgs.default(obj, termNames, comb, opt)
   names(ret)[1L] <- \"mods\"
   ret
}

coefTable.rma <- function (model, ...) {
  MuMIn:::.makeCoefTable(model$b, model$se, coefNames = rownames(model$b))
}

")

### helper functions to make mice work together with metafor (note: no longer
### needed, as there are glance and tidy methods for rma objects in broom now)

#.mice <- str2expression("
#
#glance.rma <- function (x, ...)
#   data.frame(df.residual=df.residual(x))
#
#tidy.rma <- function (x, ...) {
#   ret <- coef(summary(x))
#   colnames(ret)[2] <- \"std.error\"
#   ret$term <- rownames(ret)
#   return(ret)
#}
#
#")

############################################################################
