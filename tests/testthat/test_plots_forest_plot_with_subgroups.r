### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: http://www.metafor-project.org/doku.php/plots:forest_plot_with_subgroups

context("Checking plots example: forest plot with subgroups")

test_that("plot can be drawn.", {

   skip_on_cran()

   opar <- par(no.readonly=TRUE)

   ### decrease margins so the full space is used
   par(mar=c(4,4,1,2))

   ### load BCG vaccine data
   data(dat.bcg, package="metafor")

   ### fit random-effects model (use slab argument to define study labels)
   res <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, measure="RR",
              slab=paste(author, year, sep=", "), method="REML")

   ### set up forest plot (with 2x2 table counts added; rows argument is used
   ### to specify exactly in which rows the outcomes will be plotted)
   forest(res, xlim=c(-16, 6), at=log(c(.05, .25, 1, 4)), atransf=exp,
          ilab=cbind(dat.bcg$tpos, dat.bcg$tneg, dat.bcg$cpos, dat.bcg$cneg),
          ilab.xpos=c(-9.5,-8,-6,-4.5), cex=.75, ylim=c(-1, 27),
          order=order(dat.bcg$alloc), rows=c(3:4,9:15,20:23),
          xlab="Risk Ratio", mlab="", psize=1)

   ### add text with Q-value, dfs, p-value, and I^2 statistic
   text(-16, -1, pos=4, cex=0.75, bquote(paste("RE Model for All Studies (Q = ",
        .(formatC(res$QE, digits=2, format="f")), ", df = ", .(res$k - res$p),
        ", p = ", .(formatC(res$QEp, digits=2, format="f")), "; ", I^2, " = ",
        .(formatC(res$I2, digits=1, format="f")), "%)")))

   ### set font expansion factor (as in forest() above) and use bold italic
   ### font and save original settings in object 'op'
   op <- par(cex=.75, font=4)

   ### add text for the subgroups
   text(-16, c(24,16,5), pos=4, c("Systematic Allocation",
                                  "Random Allocation",
                                  "Alternate Allocation"))

   ### switch to bold font
   par(font=2)

   ### add column headings to the plot
   text(c(-9.5,-8,-6,-4.5), 26, c("TB+", "TB-", "TB+", "TB-"))
   text(c(-8.75,-5.25),     27, c("Vaccinated", "Control"))
   text(-16,                26, "Author(s) and Year",  pos=4)
   text(6,                  26, "Risk Ratio [95% CI]", pos=2)

   ### set par back to the original settings
   par(op)

   ### fit random-effects model in the three subgroups
   res.s <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, measure="RR",
                subset=(alloc=="systematic"), method="REML")
   res.r <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, measure="RR",
                subset=(alloc=="random"), method="REML")
   res.a <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, measure="RR",
                subset=(alloc=="alternate"), method="REML")

   ### add summary polygons for the three subgroups
   addpoly(res.s, row=18.5, cex=0.75, atransf=exp, mlab="")
   addpoly(res.r, row= 7.5, cex=0.75, atransf=exp, mlab="")
   addpoly(res.a, row= 1.5, cex=0.75, atransf=exp, mlab="")

   ### add text with Q-value, dfs, p-value, and I^2 statistic for subgroups
   text(-16, 18.5, pos=4, cex=0.75, bquote(paste("RE Model for Subgroup (Q = ",
        .(formatC(res.s$QE, digits=2, format="f")), ", df = ", .(res.s$k - res.s$p),
        ", p = ", .(formatC(res.s$QEp, digits=2, format="f")), "; ", I^2, " = ",
        .(formatC(res.s$I2, digits=1, format="f")), "%)")))
   text(-16, 7.5, pos=4, cex=0.75, bquote(paste("RE Model for Subgroup (Q = ",
        .(formatC(res.r$QE, digits=2, format="f")), ", df = ", .(res.r$k - res.r$p),
        ", p = ", .(formatC(res.r$QEp, digits=2, format="f")), "; ", I^2, " = ",
        .(formatC(res.r$I2, digits=1, format="f")), "%)")))
   text(-16, 1.5, pos=4, cex=0.75, bquote(paste("RE Model for Subgroup (Q = ",
        .(formatC(res.a$QE, digits=2, format="f")), ", df = ", .(res.a$k - res.a$p),
        ", p = ", .(formatC(res.a$QEp, digits=2, format="f")), "; ", I^2, " = ",
        .(formatC(res.a$I2, digits=1, format="f")), "%)")))

   par(opar)

})
