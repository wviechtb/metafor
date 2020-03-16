vif.rma <- function(x, intercept=FALSE, table=FALSE, digits, ...) {
   #########################################################################
   ### CJ: I have edited this function to compute the GVIF, which is valid for
   ### CJ: factor variables. This is important, because nearly all meta-analyses
   ### CJ: contain factor variables for between-study differences. The columns of
   ### CJ: the output object are labeled:
   ### CJ: GVIF:    This is the generalized VIF, as introduced by Fox & Monette, 1992
   ### CJ:          doi.org/10.1080/01621459.1992.10475190
   ### CJ: p:       Number of parameters used to compute the GVIF
   ### CJ: cGVIF:   The comparable GVIF, which can be compared across GVIFs with  
   ### CJ:          different numbers of parameters. Computed as GVIF^(1/2p)
   ### CJ: cGVIF^2: Squared comparable GVIF, which can be interpreted along the
   ### CJ:          same rules of thumb as the familiar VIF. This statistic is
   ### CJ:          typically the one that should be reported to screen for
   ### CJ:          multicolinearity.
   
   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(x, "rma"))
      stop(mstyle$stop("Argument 'x' must be an object of class \"rma\"."))

   if (inherits(x, "robust.rma"))
      stop(mstyle$stop("Function not applicable to objects of class \"robust.rma\"."))

   if (x$int.only)
      stop(mstyle$stop("VIF not applicable for intercept-only models."))
   
   ### CJ: Identify the number of predictors in the model
   mods <- labels(terms(x$formula.mods))
   if (length(mods) < 2) 
      stop("Cannot compute VIF with less than two moderators.")

   if (missing(digits)) {
      digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
   }

   #########################################################################

   vb <- vcov(x)

   if (inherits(x, "rma.ls")){
      vb <- vb$vb
   }
   
   ### CJ: Identify factor variable assignment
   var_ass <- attr(model.matrix(x$formula.mods), "assign")
   
   ### remove intercept row/colum from vb if model includes one and intercept=FALSE
   ### CJ: I am concerned that the VIF might not make sense for an incercept.
   ### CJ: although I guess the default argument intercept = FALSE addresses this.
   if(x$intercept){
     if(!intercept){
        vb <- vb[-1,-1,drop=FALSE]
        var_ass <- var_ass[-1] 
     } else {
        mods <- c("intrcpt", mods)
     }
   }

   ### rescale vb to correlation matrix
   rb <- cov2cor(vb)
   determinant <- det(rb)
   # Create output object
   out <- data.frame(
      "GVIF" = sapply(min(var_ass):max(var_ass), function(this_mod){
         i <- which(var_ass == this_mod)
         det(as.matrix(rb[i, i])) * det(as.matrix(rb[-i, -i]))/determinant
      }),
      "p" = as.vector(table(var_ass))
   )
   # Compute derived columns
   out$cGVIF <- out$GVIF^(1/(2 * out$p))
   out[["cGVIF^2"]] <- out$cGVIF^2
   rownames(out) <- mods
   
   ### add NA for intercept if model includes one and intercept=FALSE
   ### CJ: This part makes less sense if the default output of vif.rma() is always
   ### CJ: a table, which would be the case when returning the GVIF. I have tried to 
   ### CJ: make a draft of a way the previous functionality of the function can
   ### CJ: be preserved as much as possible.
   if (x$intercept && !intercept && table){
      out <- rbind(NA, out)
      rownames(out)[1] <- "intrcpt"
      # Repeat each row based on the number of factor levels, so that the cbind()
      # call below will not throw an error about different numbers of rows.
      out <- out[rep(1:nrow(out), c(1, as.vector(table(var_ass)))), ]
   }
   if (table) {
      out <- cbind(coef(summary(x)), out)
   }
   out <- round(out, digits=digits[["est"]])

   return(out)

}
