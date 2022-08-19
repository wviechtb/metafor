############################################################################

as.data.frame.vif.rma <- function(x, ...) {

   .chkclass(class(x), must="vif.rma")

   if (!is.null(x$alpha)) {

      tab <- list(beta = as.data.frame(x[[1]], ...), alpha = as.data.frame(x[[2]], ...))

   } else {

      tab <- data.frame(spec  = sapply(x$vif, function(x) x$spec),
                        coefs = sapply(x$vif, function(x) x$coefs),
                        m     = sapply(x$vif, function(x) x$m),
                        vif   = sapply(x$vif, function(x) x$vif),
                        sif   = sapply(x$vif, function(x) x$sif))

      # add proportions if they are available

      if (!is.null(x$prop))
         tab$prop <- x$prop

      #names(tab)[2] <- "coef(s)"
      #names(tab)[4] <- "(g)vif"
      #names(tab)[5] <- "(g)sif"

      # if all btt/att specifications are numeric, remove the 'spec' column

      if (all(substr(tab$spec, 1, 1) %in% as.character(1:9)))
          tab$spec <- NULL

      # just use numbers for row names when btt was specified

      if (isTRUE(x$bttspec) || isTRUE(x$attspec))
         rownames(tab) <- NULL

   }

   return(tab)

}

############################################################################
