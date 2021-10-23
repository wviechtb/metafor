points.regplot <- function(x, ...) {

   .chkclass(class(x), must="regplot")

   ### redraw points

   points(x=x$xi[x$order], y=x$yi[x$order], pch=x$pch[x$order], cex=x$psize[x$order], col=x$col[x$order], bg=x$bg[x$order], ...)

   ### redraw labels

   if (any(x$label)) {

      offset  <- attr(x, "offset")
      labsize <- attr(x, "labsize")

      for (i in which(x$label)) {

         if (isTRUE(x$yi[i] > x$pred[i])) { # x$pred might be NULL, so use isTRUE()
            text(x$xi[i], x$yi[i] + offset[1] + offset[2]*x$psize[i]^offset[3], x$slab[i], cex=labsize, ...)
         } else {
            text(x$xi[i], x$yi[i] - offset[1] - offset[2]*x$psize[i]^offset[3], x$slab[i], cex=labsize, ...)
         }

      }

   }

   invisible()

}
