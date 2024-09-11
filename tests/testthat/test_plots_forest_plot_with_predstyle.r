### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true")

source("settings.r")

context("Checking plots example: forest plot with adjusted predstyle")

dat <- dat.bcg
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat,
              slab=paste(author, year, sep=", "))
res <- rma(yi, vi, data=dat)
pred <- predict(res)

test_that("plot can be drawn with predstyle='l'.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   png("images/test_plots_forest_plot_with_predstyle_l_test.png", res=240, width=1800, height=1800, type="cairo")

   forest(res, atransf=exp, at=log(c(0.05, 0.25, 1, 4)), xlim=c(-16,6),
          ilab=cbind(tpos, tneg, cpos, cneg), ilab.lab=c("TB+","TB-","TB+","TB-"),
          ilab.xpos=c(-9.5,-8,-6,-4.5), cex=0.75, header="Author(s) and Year",
          psize=1, shade=TRUE, addpred=TRUE, predstyle="l")
   text(c(-8.75,-5.25), res$k+2.8, c("Vaccinated", "Control"), cex=0.75, font=2)

   dev.off()

   expect_true(.vistest("images/test_plots_forest_plot_with_predstyle_l_test.png", "images/test_plots_forest_plot_with_predstyle_l.png"))

   png("images/test_plots_forest_plot_with_predstyle_l_test.png", res=240, width=1800, height=1800, type="cairo")

   with(dat, forest(yi, vi, atransf=exp, at=log(c(0.05, 0.25, 1, 4)), xlim=c(-16,6),
          ilab=cbind(tpos, tneg, cpos, cneg), ilab.lab=c("TB+","TB-","TB+","TB-"),
          ilab.xpos=c(-9.5,-8,-6,-4.5), cex=0.75, header="Author(s) and Year",
          psize=1, shade=TRUE, ylim=c(-2,16)))
   text(c(-8.75,-5.25), res$k+2.8, c("Vaccinated", "Control"), cex=0.75, font=2)
   abline(h=0)
   with(pred, addpoly(pred, ci.lb=ci.lb, ci.ub=ci.ub, pi.lb=pi.lb, pi.ub=pi.ub, mlab="Random-Effects Model"))

   dev.off()

   expect_true(.vistest("images/test_plots_forest_plot_with_predstyle_l_test.png", "images/test_plots_forest_plot_with_predstyle_l.png"))

   png("images/test_plots_forest_plot_with_predstyle_l_test.png", res=240, width=1800, height=1800, type="cairo")

   with(dat, forest(yi, vi, atransf=exp, at=log(c(0.05, 0.25, 1, 4)), xlim=c(-16,6),
          ilab=cbind(tpos, tneg, cpos, cneg), ilab.lab=c("TB+","TB-","TB+","TB-"),
          ilab.xpos=c(-9.5,-8,-6,-4.5), cex=0.75, header="Author(s) and Year",
          psize=1, shade=TRUE, ylim=c(-2,16)))
   text(c(-8.75,-5.25), res$k+2.8, c("Vaccinated", "Control"), cex=0.75, font=2)
   abline(h=0)
   addpoly(res, row=-1, addpred=TRUE)

   dev.off()

   expect_true(.vistest("images/test_plots_forest_plot_with_predstyle_l_test.png", "images/test_plots_forest_plot_with_predstyle_l.png"))

   png("images/test_plots_forest_plot_with_predstyle_l_test.png", res=240, width=1800, height=1800, type="cairo")

   with(dat, forest(yi, vi, atransf=exp, at=log(c(0.05, 0.25, 1, 4)), xlim=c(-16,6),
          ilab=cbind(tpos, tneg, cpos, cneg), ilab.lab=c("TB+","TB-","TB+","TB-"),
          ilab.xpos=c(-9.5,-8,-6,-4.5), cex=0.75, header="Author(s) and Year",
          psize=1, shade=TRUE, ylim=c(-2,16)))
   text(c(-8.75,-5.25), res$k+2.8, c("Vaccinated", "Control"), cex=0.75, font=2)
   abline(h=0)
   addpoly(pred, rows=-1, addpred=TRUE, mlab="Random-Effects Model")

   dev.off()

   expect_true(.vistest("images/test_plots_forest_plot_with_predstyle_l_test.png", "images/test_plots_forest_plot_with_predstyle_l.png"))

})

test_that("plot can be drawn with predstyle='b'.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   png("images/test_plots_forest_plot_with_predstyle_b_test.png", res=240, width=1800, height=1800, type="cairo")

   forest(res, atransf=exp, at=log(c(0.05, 0.25, 1, 4)), xlim=c(-16,6),
          ilab=cbind(tpos, tneg, cpos, cneg), ilab.lab=c("TB+","TB-","TB+","TB-"),
          ilab.xpos=c(-9.5,-8,-6,-4.5), cex=0.75, header="Author(s) and Year",
          psize=1, shade=TRUE, addpred=TRUE, predstyle="b")
   text(c(-8.75,-5.25), res$k+2.8, c("Vaccinated", "Control"), cex=0.75, font=2)

   dev.off()

   expect_true(.vistest("images/test_plots_forest_plot_with_predstyle_b_test.png", "images/test_plots_forest_plot_with_predstyle_b.png"))

   png("images/test_plots_forest_plot_with_predstyle_b_test.png", res=240, width=1800, height=1800, type="cairo")

   with(dat, forest(yi, vi, atransf=exp, at=log(c(0.05, 0.25, 1, 4)), xlim=c(-16,6),
          ilab=cbind(tpos, tneg, cpos, cneg), ilab.lab=c("TB+","TB-","TB+","TB-"),
          ilab.xpos=c(-9.5,-8,-6,-4.5), cex=0.75, header="Author(s) and Year",
          psize=1, shade=TRUE, ylim=c(-3,16)))
   text(c(-8.75,-5.25), res$k+2.8, c("Vaccinated", "Control"), cex=0.75, font=2)
   abline(h=0)
   with(pred, addpoly(pred, ci.lb=ci.lb, ci.ub=ci.ub, pi.lb=pi.lb, pi.ub=pi.ub, mlab="Random-Effects Model", predstyle="b"))

   dev.off()

   expect_true(.vistest("images/test_plots_forest_plot_with_predstyle_b_test.png", "images/test_plots_forest_plot_with_predstyle_b.png"))

   png("images/test_plots_forest_plot_with_predstyle_b_test.png", res=240, width=1800, height=1800, type="cairo")

   with(dat, forest(yi, vi, atransf=exp, at=log(c(0.05, 0.25, 1, 4)), xlim=c(-16,6),
          ilab=cbind(tpos, tneg, cpos, cneg), ilab.lab=c("TB+","TB-","TB+","TB-"),
          ilab.xpos=c(-9.5,-8,-6,-4.5), cex=0.75, header="Author(s) and Year",
          psize=1, shade=TRUE, ylim=c(-3,16)))
   text(c(-8.75,-5.25), res$k+2.8, c("Vaccinated", "Control"), cex=0.75, font=2)
   abline(h=0)
   addpoly(res, row=-1, predstyle="b")

   dev.off()

   expect_true(.vistest("images/test_plots_forest_plot_with_predstyle_b_test.png", "images/test_plots_forest_plot_with_predstyle_b.png"))

   png("images/test_plots_forest_plot_with_predstyle_b_test.png", res=240, width=1800, height=1800, type="cairo")

   with(dat, forest(yi, vi, atransf=exp, at=log(c(0.05, 0.25, 1, 4)), xlim=c(-16,6),
          ilab=cbind(tpos, tneg, cpos, cneg), ilab.lab=c("TB+","TB-","TB+","TB-"),
          ilab.xpos=c(-9.5,-8,-6,-4.5), cex=0.75, header="Author(s) and Year",
          psize=1, shade=TRUE, ylim=c(-3,16)))
   text(c(-8.75,-5.25), res$k+2.8, c("Vaccinated", "Control"), cex=0.75, font=2)
   abline(h=0)
   addpoly(pred, rows=-1, predstyle="b", mlab="Random-Effects Model")

   dev.off()

   expect_true(.vistest("images/test_plots_forest_plot_with_predstyle_b_test.png", "images/test_plots_forest_plot_with_predstyle_b.png"))

})

test_that("plot can be drawn with predstyle='s'.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   png("images/test_plots_forest_plot_with_predstyle_s_test.png", res=240, width=1800, height=1800, type="cairo")

   forest(res, atransf=exp, at=log(c(0.05, 0.25, 1, 4)), xlim=c(-16,6),
          ilab=cbind(tpos, tneg, cpos, cneg), ilab.lab=c("TB+","TB-","TB+","TB-"),
          ilab.xpos=c(-9.5,-8,-6,-4.5), cex=0.75, header="Author(s) and Year",
          psize=1, shade=TRUE, addpred=TRUE, predstyle="s")
   text(c(-8.75,-5.25), res$k+2.8, c("Vaccinated", "Control"), cex=0.75, font=2)

   dev.off()

   expect_true(.vistest("images/test_plots_forest_plot_with_predstyle_s_test.png", "images/test_plots_forest_plot_with_predstyle_s.png"))

   png("images/test_plots_forest_plot_with_predstyle_s_test.png", res=240, width=1800, height=1800, type="cairo")

   with(dat, forest(yi, vi, atransf=exp, at=log(c(0.05, 0.25, 1, 4)), xlim=c(-16,6),
          ilab=cbind(tpos, tneg, cpos, cneg), ilab.lab=c("TB+","TB-","TB+","TB-"),
          ilab.xpos=c(-9.5,-8,-6,-4.5), cex=0.75, header="Author(s) and Year",
          psize=1, shade=TRUE, ylim=c(-3,16)))
   text(c(-8.75,-5.25), res$k+2.8, c("Vaccinated", "Control"), cex=0.75, font=2)
   abline(h=0)
   with(pred, addpoly(pred, ci.lb=ci.lb, ci.ub=ci.ub, pi.lb=pi.lb, pi.ub=pi.ub, mlab="Random-Effects Model", predstyle="s"))

   dev.off()

   expect_true(.vistest("images/test_plots_forest_plot_with_predstyle_s_test.png", "images/test_plots_forest_plot_with_predstyle_s.png"))

   png("images/test_plots_forest_plot_with_predstyle_s_test.png", res=240, width=1800, height=1800, type="cairo")

   with(dat, forest(yi, vi, atransf=exp, at=log(c(0.05, 0.25, 1, 4)), xlim=c(-16,6),
          ilab=cbind(tpos, tneg, cpos, cneg), ilab.lab=c("TB+","TB-","TB+","TB-"),
          ilab.xpos=c(-9.5,-8,-6,-4.5), cex=0.75, header="Author(s) and Year",
          psize=1, shade=TRUE, ylim=c(-3,16)))
   text(c(-8.75,-5.25), res$k+2.8, c("Vaccinated", "Control"), cex=0.75, font=2)
   abline(h=0)
   addpoly(res, row=-1, predstyle="s")

   dev.off()

   expect_true(.vistest("images/test_plots_forest_plot_with_predstyle_s_test.png", "images/test_plots_forest_plot_with_predstyle_s.png"))

   png("images/test_plots_forest_plot_with_predstyle_s_test.png", res=240, width=1800, height=1800, type="cairo")

   with(dat, forest(yi, vi, atransf=exp, at=log(c(0.05, 0.25, 1, 4)), xlim=c(-16,6),
          ilab=cbind(tpos, tneg, cpos, cneg), ilab.lab=c("TB+","TB-","TB+","TB-"),
          ilab.xpos=c(-9.5,-8,-6,-4.5), cex=0.75, header="Author(s) and Year",
          psize=1, shade=TRUE, ylim=c(-3,16)))
   text(c(-8.75,-5.25), res$k+2.8, c("Vaccinated", "Control"), cex=0.75, font=2)
   abline(h=0)
   addpoly(pred, rows=-1, predstyle="s", mlab="Random-Effects Model")

   dev.off()

   expect_true(.vistest("images/test_plots_forest_plot_with_predstyle_s_test.png", "images/test_plots_forest_plot_with_predstyle_s.png"))

})

test_that("plot can be drawn with predstyle='d'.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   png("images/test_plots_forest_plot_with_predstyle_d_test.png", res=240, width=1800, height=1800, type="cairo")

   forest(res, atransf=exp, at=log(c(0.05, 0.25, 1, 4)), xlim=c(-16,6),
          ilab=cbind(tpos, tneg, cpos, cneg), ilab.lab=c("TB+","TB-","TB+","TB-"),
          ilab.xpos=c(-9.5,-8,-6,-4.5), cex=0.75, header="Author(s) and Year",
          psize=1, shade=TRUE, addpred=TRUE, predstyle="d")
   text(c(-8.75,-5.25), res$k+2.8, c("Vaccinated", "Control"), cex=0.75, font=2)

   dev.off()

   expect_true(.vistest("images/test_plots_forest_plot_with_predstyle_d_test.png", "images/test_plots_forest_plot_with_predstyle_d.png"))

   png("images/test_plots_forest_plot_with_predstyle_d_test.png", res=240, width=1800, height=1800, type="cairo")

   with(dat, forest(yi, vi, atransf=exp, at=log(c(0.05, 0.25, 1, 4)), xlim=c(-16,6),
          ilab=cbind(tpos, tneg, cpos, cneg), ilab.lab=c("TB+","TB-","TB+","TB-"),
          ilab.xpos=c(-9.5,-8,-6,-4.5), cex=0.75, header="Author(s) and Year",
          psize=1, shade=TRUE, ylim=c(-3,16)))
   text(c(-8.75,-5.25), res$k+2.8, c("Vaccinated", "Control"), cex=0.75, font=2)
   abline(h=0)
   with(pred, addpoly(pred, ci.lb=ci.lb, ci.ub=ci.ub, pi.lb=pi.lb, pi.ub=pi.ub, mlab="Random-Effects Model", predstyle="d"))

   dev.off()

   expect_true(.vistest("images/test_plots_forest_plot_with_predstyle_d_test.png", "images/test_plots_forest_plot_with_predstyle_d.png"))

   png("images/test_plots_forest_plot_with_predstyle_d_test.png", res=240, width=1800, height=1800, type="cairo")

   with(dat, forest(yi, vi, atransf=exp, at=log(c(0.05, 0.25, 1, 4)), xlim=c(-16,6),
          ilab=cbind(tpos, tneg, cpos, cneg), ilab.lab=c("TB+","TB-","TB+","TB-"),
          ilab.xpos=c(-9.5,-8,-6,-4.5), cex=0.75, header="Author(s) and Year",
          psize=1, shade=TRUE, ylim=c(-3,16)))
   text(c(-8.75,-5.25), res$k+2.8, c("Vaccinated", "Control"), cex=0.75, font=2)
   abline(h=0)
   addpoly(res, row=-1, predstyle="d")

   dev.off()

   expect_true(.vistest("images/test_plots_forest_plot_with_predstyle_d_test.png", "images/test_plots_forest_plot_with_predstyle_d.png"))

   png("images/test_plots_forest_plot_with_predstyle_d_test.png", res=240, width=1800, height=1800, type="cairo")

   with(dat, forest(yi, vi, atransf=exp, at=log(c(0.05, 0.25, 1, 4)), xlim=c(-16,6),
          ilab=cbind(tpos, tneg, cpos, cneg), ilab.lab=c("TB+","TB-","TB+","TB-"),
          ilab.xpos=c(-9.5,-8,-6,-4.5), cex=0.75, header="Author(s) and Year",
          psize=1, shade=TRUE, ylim=c(-3,16)))
   text(c(-8.75,-5.25), res$k+2.8, c("Vaccinated", "Control"), cex=0.75, font=2)
   abline(h=0)
   addpoly(pred, rows=-1, predstyle="d", mlab="Random-Effects Model")

   dev.off()

   expect_true(.vistest("images/test_plots_forest_plot_with_predstyle_d_test.png", "images/test_plots_forest_plot_with_predstyle_d.png"))

})

test_that("plot can be drawn with predstyle='d' and transf=exp.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   png("images/test_plots_forest_plot_with_predstyle_d_transf_test.png", res=240, width=1800, height=1800, type="cairo")

   forest(res, transf=exp, alim=c(0,3), steps=4, xlim=c(-2,4.2), cex=0.75,
          header="Author(s) and Year", psize=1, refline=1, shade=TRUE,
          addpred=TRUE, predstyle="d")

   dev.off()

   expect_true(.vistest("images/test_plots_forest_plot_with_predstyle_d_transf_test.png", "images/test_plots_forest_plot_with_predstyle_d_transf.png"))

   png("images/test_plots_forest_plot_with_predstyle_d_transf_test.png", res=240, width=1800, height=1800, type="cairo")

   with(dat, forest(yi, vi, transf=exp, alim=c(0,3), steps=4, xlim=c(-2,4.2), cex=0.75,
          header="Author(s) and Year", psize=1, refline=1, shade=TRUE,
          ylim=c(-3,16)))
   abline(h=0)
   with(pred, addpoly(pred, ci.lb=ci.lb, ci.ub=ci.ub, pi.lb=pi.lb, pi.ub=pi.ub, mlab="Random-Effects Model", predstyle="d"))

   dev.off()

   expect_true(.vistest("images/test_plots_forest_plot_with_predstyle_d_transf_test.png", "images/test_plots_forest_plot_with_predstyle_d_transf.png"))

   png("images/test_plots_forest_plot_with_predstyle_d_transf_test.png", res=240, width=1800, height=1800, type="cairo")

   with(dat, forest(yi, vi, transf=exp, alim=c(0,3), steps=4, xlim=c(-2,4.2), cex=0.75,
          header="Author(s) and Year", psize=1, refline=1, shade=TRUE,
          ylim=c(-3,16)))
   abline(h=0)
   addpoly(res, row=-1, predstyle="d")

   dev.off()

   expect_true(.vistest("images/test_plots_forest_plot_with_predstyle_d_transf_test.png", "images/test_plots_forest_plot_with_predstyle_d_transf.png"))

   png("images/test_plots_forest_plot_with_predstyle_d_transf_test.png", res=240, width=1800, height=1800, type="cairo")

   with(dat, forest(yi, vi, transf=exp, alim=c(0,3), steps=4, xlim=c(-2,4.2), cex=0.75,
          header="Author(s) and Year", psize=1, refline=1, shade=TRUE,
          ylim=c(-3,16)))
   abline(h=0)
   addpoly(pred, rows=-1, predstyle="d", mlab="Random-Effects Model")

   dev.off()

   expect_true(.vistest("images/test_plots_forest_plot_with_predstyle_d_transf_test.png", "images/test_plots_forest_plot_with_predstyle_d_transf.png"))

})

rm(list=ls())
