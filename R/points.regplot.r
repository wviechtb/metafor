points.regplot <- function(x, ...) {

   .chkclass(class(x), must="regplot")

   points(x=x$xi[x$order], y=x$yi[x$order], pch=x$pch[x$order], cex=x$psize[x$order], col=x$col[x$order], bg=x$bg[x$order])

   invisible()

}
