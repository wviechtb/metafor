# metafor 4.6-0 (2024-03-28)

- the `steps` argument in the various `profile()` functions can now also be a numeric vector to specify for which parameter values the likelihood should be evaluated

- a few minor fixes to the dynamic theming of plots based on the foreground and background colors of the plotting device

- slightly improved flexibility for setting package options

- new measures added to `escalc()`: `"SMN"` for the single-group standardized mean / single-group standardized mean difference, `"SMCRP"` for the standardized mean change using raw score standardization with pooled standard deviations, and `"SMCRPH"` for the standardized mean change using raw score standardization with pooled standard deviations and heteroscedastic population variances at the two measurement occasions

- calculation of the sampling variances for measures `"SMDH"`, `"SMD1H"`, and `"SMCRH"` was slightly adjusted for consistency

- in `plot.gosh.rma()`, can also set `het="tau"` (to plot the square root of tau^2 as the measure of heterogeneity)

- in the various `forest()` functions, argument `ylim` can now only be a single value to specify the lower bound (while the upper bound is still set automatically)

- in `forest()` and `regplot()`, observation limits set via `olim` are now properly applied to all elements

- various internal improvements to `selmodel()`

- `selmodel()` no longer stops with an error when one or more intervals defined by the `steps` argument do not contain any observed p-values (instead a warning is issued and model fitting proceeds, but may fail)

- added `decreasing` argument to `selmodel()` for enforcing that the delta estimates must be a monotonically decreasing function of the p-values in the step function model

- added the undocumented argument `pval` to `selmodel()` for passing p-values directly to the function (doing this is highly experimental)

- some internal refactoring of the code

- improved the documentation a bit

# metafor 4.4-0 (2023-09-27)

- added `getmfopt()` and `setmfopt()` functions for getting and setting package options and made some of the options more flexible

- removed argument `weighted` from `fsn()` (whether weighted or unweighted averages are used in Orwin's method is now simply determined by whether sampling variances are specified or not); added `type="General"` to `fsn()` as a generalization of the Orwin and Rosenberg methods (that allows for a fail-safe N calculation based on a random-effects model); can now pass an `rma` object to the `fsn()` function

- further improved the theming of all plots based on the foreground and background colors; within RStudio, plot colors can also be automatically chosen based on the theme (with `setmfopt(theme="auto")`)

- added additional/optional argument `tabfig` to the various `forest()` functions, for easily setting the `annosym` argument to an appropriate vector for exactly aligning numbers (when using a matching font)

- added (for now undocumented) `vccon` argument to `rma.mv()` for setting equality constraints on variance/correlation components

- `replace` argument in `conv.2x2()`, `conv.delta()`, `conv.fivenum()`, and `conv.wald()` can now also be a logical

- added `summary.matreg()` and `print.summary.matreg()` methods for including additional statistics in the output (R^2 and the omnibus test) and added `coef.matreg()` and `vcov.matreg()` extractor functions

- formatting functions `fmtp()`, `fmtx()`, and `fmtt()` gain a `quote` argument, which is set to `FALSE` by default

- for measures `"PCOR"`, `"ZPCOR"`, `"SPCOR"`, and `"ZSPCOR"`, argument `mi` in `escalc()` now refers to the total number of predictors in the regression models (i.e., also counting the focal predictor of interest)

- added measures `"R2"` and "`ZR2"` to `escalc()`

- `addpoly.default()` and `addpoly.rma.predict()` gain a `constarea` argument (for the option to draw the polygons with a constant area)

- `plot.rma.uni.selmodel()` gains a `shade` argument (for shading the confidence interval region)

- `plot.permutest.rma.uni()` gains a `legend` argument

- `vcalc()` gains a `sparse` argument

- `aggregate.escalc` gains `var.names` argument

- made the `legend` argument more flexible in `funnel()`

- made the `append` argument more flexible in `to.long()`

- added a few more transformation functions

- small bug fixes

- added automated visual comparison tests of plots

- improved the documentation a bit

# metafor 4.2-0 (2023-05-08)

- improved the various plotting functions so they respect `par("fg")`; as a result, one can now create plots with a dark background and light plotting colors

- also allow two or three values for `xlab` in the various `forest()` functions (for adding labels at the ends of the x-axis limits)

- better default choices for `xlim` in the various `forest()` functions; also, argument `ilab.xpos` is now optional when using the `ilab` argument

- added `shade` and `colshade` arguments to the various `forest()` functions

- the various `forest()` functions no longer enforce that `xlim` must be at least as wide as `alim`

- added `link` argument to `rma.glmm()`

- `rma.glmm()` with `measure="OR", model="CM.EL", method="ML"` now treats tau^2 values below 1e-04 effectively as zero before computing the standard errors of the fixed effects; this helps to avoid numerical problems in approximating the Hessian; similarly, `selmodel()` now treats tau^2 values below 1e-04 or min(vi/10) effectively as zero before computing the standard errors

- for measure `SMCC`, can now specify d-values, t-test statistics, and p-values via arguments `di`, `ti`, and `pi`

- functions that issue a warning when omitting studies due to NAs now indicate how many were omitted

- properly documented the `level` argument

- added a few more transformation functions

- small bug fixes

- improved the documentation a bit

# metafor 4.0-0 (2023-03-19)

- added `conv.2x2()` function for reconstructing the cell frequencies in 2x2 tables based on other summary statistics

- added `conv.wald()` function for converting Wald-type confidence intervals and test statistics to sampling variances

- added `conv.fivenum()` function for estimating means and standard deviations from five-number summary values

- added `conv.delta()` function for transforming observed effect sizes or outcomes and their sampling variances using the delta method

- added `emmprep()` function to create a reference grid for use with the `emmeans()` function from the package of the same name

- exposed formatter functions `fmtp()`, `fmtx()`, and `fmtt()`

- package `numDeriv` moved from `Suggests` to `Depends`

- `model.matrix.rma()` gains `asdf` argument

- corrected bug in `vcalc()` (values for `obs` and `type` were taken directly as indices instead of using them as identifiers)

- improved efficiency of `vif()` when `sim=TRUE` by reshuffling only the data needed in the model matrix; due to some edge cases, the simulation approach cannot be used when some redundant predictors were dropped from the original model; and when redundancies occur after reshuffling the data, the simulated (G)VIF value(s) are now set to `Inf` instead of `NA`

- `selmodel()` gains `type='trunc'` and `type='truncest'` models (the latter should be considered experimental)

- added `exact="i"` option in `permutest()` (to just return the number of iterations required for an exact permutation test)

- `escalc()` now provides more informative error messages when not specifying all required arguments to compute a particular measure

- added measures `"ZPHI"`, `"ZTET"`, `"ZPB"`, `"ZBIS"`, and `"ZSPCOR"` to `escalc()` (but note that Fisher's r-to-z transformation is not a variance-stabilizing transformation for these measures)

- the variance of measure `ZPCOR` is now calculated with `1/(ni-mi-3)` (instead of `1/(ni-mi-1)`), which provides a better approximation in small samples (and analogous to how the variance of `ZCOR` is calculated with `1/(ni-3)`)

- as with `measure="SMD"`, one can now also use arguments `di` and `ti` to specify d-values and t-test statistics for measures `RPB`, `RBIS`, `D2ORN`, and `D2ORL` in `escalc()`

- for measures `COR`, `UCOR`, and `ZCOR`, can now use argument `ti` to specify t-test statistics in `escalc()`

- can also specify (two-sided) p-values (of the respective t-tests) for these measures (and for measures `PCOR`, `ZPCOR`, `SPCOR`, and `ZSPCOR`) via argument `pi` (the sign of the p-value is taken to be the sign of the measure)

- can also specify (semi-)partial correlations directly via argument `ri` for measures `PCOR`, `ZPCOR`, `SPCOR`, and `ZSPCOR`

- when passing a correlation marix to `rcalc()`, it now orders the elements (columnwise) based on the lower triangular part of the matrix, not the upper one (which is more consistent with what `matreg()` expects as input when using the `V` argument)

- optimizers `Rcgmin` and `Rvmmin` are now available in `rma.uni()`, `rma.mv()`, `rma.glmm()`, and `selmodel()`

- improved the documentation a bit

# metafor 3.8-1 (2022-08-26)

- `funnel.default()`, `funnel.rma()`, and `regplot.rma()` gain `slab` argument

- `vif()` was completely refactored and gains `reestimate`, `sim`, and `parallel` arguments; added `as.data.frame.vif.rma()` and `plot.vif.rma()` methods

- `plot.permutest.rma.uni()` function sets the y-axis limits automatically and in a smarter way when also drawing the reference/null distribution and the density estimate

- added possibility to specify a list for `btt` in `anova.rma()`; added `print.list.anova.rma()` to print the resulting object

- added `as.data.frame.anova.rma()` and `as.data.frame.list.anova.rma()` methods

- documented the possibility to use an identity link (with `link="identity"`) in `rma.uni()` when fitting location-scale models (although this will often lead to estimation problems); added `solnp()` as an additional optimizer for this case

- optimizers `nloptr` and `constrOptim.nl` (the latter from the `alabama` package) are now available in `rma.uni()` for location-scale models when using an identity link

- added measure `SMD1H` to `escalc()`

- for `measure="SMD"`, `escalc()` now also allows the user to specify d-values and t-test statistics via arguments `di` and `ti`, respectively

- `aggregate.escalc()` gains `addk` argument

- added (experimental!) support for measures `"RR"`, `"RD"`, `"PLN"`, and `"PR"` to `rma.glmm()` (but using these measures will often lead to estimation problems)

- `replmiss()` gains `data` argument

- `cumul()` functions also store data, so that arguments `ilab`, `col`, `pch`, and `psize` in the `forest.cumul.rma()` function can look for variables therein

- fixed issue with rendering Rmarkdown documents with `metafor` output due to the use of a zero-width space

# metafor 3.4-0 (2022-04-21)

- added `misc-models`, `misc-recs`, and `misc-options` help pages

- added `as.data.frame.confint.rma()` and `as.data.frame.list.confint.rma` methods

- `permutest()` can now also do permutation tests for location-scale models; it also always returns the permutation distributions; hence, argument `retpermdist` was removed

- added `plot.permutest.rma.uni()` function to plot the permutation distributions

- simplified `regtest()`, `ranktest()`, and `tes()` to single functions instead of using generics and methods; this way, a `data` argument could be added

- added `vcalc()` and `blsplit()` functions

- `robust()` gains `clubSandwich` argument; if set to `TRUE`, the methods from the `clubSandwich` package (https://cran.r-project.org/package=clubSandwich) are used to obtain the cluster-robust results; `anova.rma()` and `predict.rma()` updated to work appropriately in this case

- results from `robust()` are no longer printed with `print.robust.rma()` but with the print methods `print.rma.uni()` and `print.rma.mv()`

- `anova.rma()` now gives a warning when running LRTs not based on ML/REML estimation and gains `rhs` argument; it also now has a `refit` argument (to refit REML fits with ML in case the fixed effects of the models differ)

- setting `dfs="contain"` in `rma.mv()` automatically sets `test="t"` for convenience

- elements of `rho` and `phi` in `rma.mv()` are now based on the lower triangular part of the respective correlation matrix (instead of the upper triangular part) for consistency with other functions; note that this is in principle a backwards incompatible change, although this should only be a concern in very special circumstances

- `rma.mv()` gains `cvvc` argument (for calculating the var-cov matrix of the variance/correlation/covariance components)

- added measure `"MPORM"` to `escalc()` for computing marginal log odds ratios based on marginal 2x2 tables directly (which requires specification of the correlation coefficients in the paired tables for the calculation of the sampling variances via the `ri` argument)

- added measure `"REH"` to `escalc()` for computing the (log transformed) relative excess heterozygosity (to assess deviations from the Hardy-Weinberg equilibrium)

- `aggregate.escalc()` gains `checkpd` argument and `struct="CS+CAR"`

- `rma.glmm()` now has entire array of optimizers available for `model="CM.EL"` and `measure="OR"`; switched the default from `optim()` with method `BFGS` to `nlminb()` for consistency with `rma.mv()`, `rma.uni()`, and `selmodel.rma.uni()`

- `rma.glmm()` gains `coding` and `cor` arguments and hence more flexibility how the group variable should be coded in the random effects structure and whether the random study effects should be allowed to be correlated with the random group effects

- `rma.uni()` now also provides R^2 for fixed-effects models

- `matreg()` can now also analyze a covariance matrix with a corresponding `V` matrix; can also specify variable names (instead of indices) for arguments `x` and `y`

- renamed argument `nearPD` to `nearpd` in `matreg()` (but `nearPD` continues to work)

- `plot.profile.rma()` gains `refline` argument

- added `addpoly.rma.predict()` method

- `addpoly.default()` and `addpoly.rma()` gain `lty` and `annosym` arguments; if unspecified, arguments `annotate`, `digits`, `width`, `transf`, `atransf`, `targs`, `efac`, `fonts`, `cex`, and `annosym` are now automatically set equal to the same values that were used when creating the forest plot

- documented `textpos` and `rowadj` arguments for the various `forest` functions and moved the `top` and `annosym` arguments to 'additional arguments'

- fixed that `level` argument in `addpoly.rma()` did not affect the CI width

- `points.regplot()` function now also redraws the labels (if there were any to begin with)

- added `lbfgsb3c`, `subplex`, and `BBoptim` as possible optimizer in `rma.mv()`, `rma.glmm()`, `rma.uni()`, and `selmodel.rma.uni()`

- the object returned by model fitting functions now includes the data frame specified via the `data` argument; various method functions now automatically look for specified variables within this data frame first

- datasets moved to the `metadat` package (https://cran.r-project.org/package=metadat)

- improved the documentation a bit

# metafor 3.0-2 (2021-06-09)

- the `metafor` package now makes use of the `mathjaxr` package to nicely render equations shown in the HTML help pages

- `rma()` can now also fit location-scale models

- added `selmodel()` for fitting a wide variety of selection models (and added the corresponding `plot.rma.uni.selmodel()` function for drawing the estimated selection function)

- `rma.mv()` gains `dfs` argument and now provides an often better way for calculating the (denominator) degrees of freedom for approximate t- and F-tests when `dfs="contain"`

- added `tes()` function for the test of excess significance

- added `regplot()` function for drawing scatter plots / bubble plots based on meta-regression models

- added `rcalc()` for calculating the variance-covariance matrix of correlation coefficients and `matreg()` for fitting regression models based on correlation/covariance matrices

- added convenience functions `dfround()` and `vec2mat()`

- added `aggregate.escalc()` function to aggregate multiple effect sizes or outcomes within studies/clusters

- `regtest()` now shows the 'limit estimate' of the (average) true effect when using `sei`, `vi`, `ninv`, or `sqrtninv` as predictors (and the model does not contain any other moderators)

- `vif()` gains `btt` argument and can now also compute generalized variance inflation factors; a proper `print.vif.rma()` function was also added

- `anova.rma()` argument `L` renamed to `X` (the former still works, but is no longer documented)

- argument `order` in `cumul()` should now just be a variable, not the order of the variable, to be used for ordering the studies and must be of the same length as the original dataset that was used in the model fitting

- similarly, vector arguments in various plotting functions such as `forest.rma()` must now be of the same length as the original dataset that was used in the model fitting (any subsetting and removal of `NA`s is automatically applied)

- the various `leave1out()` and `cumul()` functions now provide `I^2` and `H^2` also for fixed-effects models; accordingly, `plot.cumul.rma()` now also works with such models

- fixed `level` not getting passed down to the various `cumul()` functions

- `plot.cumul.rma()` argument `addgrid` renamed to `grid` (the former still works, but is no longer documented)

- `forest.default()`, `forest.rma()`, and `labbe()` gain `plim` argument and now provide more flexibility in terms of the scaling of the points

- `forest.rma()` gains `colout` argument (to adjust the color of the observed effect sizes or outcomes)

- in the various `forest()` functions, the right header is now suppressed when `annotate=FALSE` and `header=TRUE`

- `funnel.default()` and `funnel.rma()` gain `label` and `offset` arguments

- `funnel.default()` and `funnel.rma()` gain `lty` argument; the reference line is now drawn by default as a dotted line (like the line for the pseudo confidence region)

- the `forest` and `funnel` arguments of `reporter.rma.uni()` can now also be logicals to suppress the drawing of these plots

- added `weighted` argument to `fsn()` (for Orwin's method)

- added some more transformation functions

- `bldiag()` now properly handles ?x0 or 0x? matrices

- p-values are still given to 2 digits even when `digits=1`

- `summary.escalc()` also provides the p-values (of the Wald-type tests); but when using the `transf` argument, the sampling variances, standard errors, test statistics, and p-values are no longer shown

- `rma.uni()` no longer constrains a fixed tau^2 value to 0 when k=1

- slight speedup in functions that repeatedly fit `rma.uni()` models by skipping the computation of the pseudo R^2 statistic

- started using the `pbapply` package for showing progress bars, also when using parallel processing

- to avoid potential confusion, all references to 'credibility intervals' have been removed from the documentation; these intervals are now exclusively referred to as 'prediction intervals'; in the output, the bounds are therefore indicated now as `pi.lb` and `pi.ub` (instead of `cr.lb` and `cr.ub`); the corresponding argument names were changed in `addpoly.default()`; argument `addcred` was changed to `addpred` in `addpoly.rma()` and `forest.rma()`; however, code using the old arguments names should continue to work

- one can now use `weights(...,  type="rowsum")` for intercept-only `rma.mv` models (to obtain 'row-sum weights')

- `simulate.rma()` gains `olim` argument; renamed the `clim` argument in `summary.escalc()` and the various `forest()` functions to `olim` for consistency (the old `clim` argument should continue to work)

- show nicer network graphs for `dat.hasselblad1998` and `dat.senn2013` in the help files

- added 24 datasets (`dat.anand1999`, `dat.assink2016`, `dat.baskerville2012`, `dat.bornmann2007`, `dat.cannon2006`, `dat.cohen1981`, `dat.craft2003`, `dat.crede2010`, `dat.dagostino1998`, `dat.damico2009`, `dat.dorn2007`, `dat.hahn2001`, `dat.kalaian1996`, `dat.kearon1998`, `dat.knapp2017`, `dat.landenberger2005`, `dat.lau1992`, `dat.lim2014`, `dat.lopez2019`, `dat.maire2019, `, `dat.moura2021` `dat.obrien2003`, `dat.vanhowe1999`, `dat.viechtbauer2021`)

- the package now runs a version check on startup in interactive sessions; setting the environment variable `METAFOR_VERSION_CHECK` to `FALSE` disables this

- refactored various functions (for cleaner/simpler code)

- improved the documentation a bit

# metafor 2.4-0 (2020-03-19)

- version jump to 2.4-0 for CRAN release (from now on, even minor numbers for CRAN releases, odd numbers for development versions)

- the various `forest()` functions gain `header` argument

- `escalc()` gains `include` argument

- setting `verbose=3` in model fitting functions sets `options(warn=1)`

- `forest.rma()` and `forest.default()` now throw informative errors when misusing `order` and `subset` arguments

- fixed failing tests due to the `stringsAsFactors=FALSE` change in the upcoming version of R

- `print.infl.rma.uni()` gains `infonly` argument, to only show the influential studies

- removed `MASS` from `Suggests` (no longer needed)

- argument `btt` can now also take a string to grep for

- added `optimParallel` as possible optimizer in `rma.mv()`

- added (for now undocumented) option to fit models in `rma.glmm()` via the `GLMMadaptive` package (instead of `lme4`); to try this, use: `control=list(package="GLMMadaptive")`

- started to use numbering scheme for devel version (the number after the dash indicates the devel version)

- added `contrmat()` function (for creating a matrix that indicates which groups have been compared against each other in each row of a dataset)

- added `to.wide()` function (for restructuring long format datasets into the wide format needed for contrast-based analyses)

- `I^2` and `H^2` are also shown in output for fixed-effects models

- argument `grid` in `baujat()` can now also be a color name

- added (for now undocumented) `time` argument to more functions that are computationally expensive

- added (for now undocumented) `textpos` argument to the various forest functions

- added a new dataset (`dat.graves2010`)

- added more tests

# metafor 2.1-0 (2019-05-13)

- added `formula()` method for objects of class `rma`

- `llplot()` now also allows for `measure="GEN"`; also, the documentation and y-axis label have been corrected to indicate that the function plots likelihoods (not log likelihoods)

- `confint.rma.mv()` now returns an object of class `list.confint.rma` when obtaining CIs for all variance and correlation components of the model; added corresponding `print.list.confint.rma()` function

- moved `tol` argument in `permutest()` to `control` and renamed to `comptol`

- added `PMM` and `GENQM` estimators in `rma.uni()`

- added `vif()` function to get variance inflation factors

- added `.glmulti` object for making the interaction with `glmulti` easier

- added `reporter()` and `reporter.rma.uni()` for dynamically generating analysis reports for objects of class `rma.uni`

- output is now styled/colored when `crayon` package is loaded (this only works on a 'proper' terminal with color support; also works in RStudio)

- overhauled `plot.gosh.rma()`; when `out` is specified, it now shows two distributions, one for the values when the outlier is included and one for the values when for outlier is excluded; dropped the `hcol` argument and added `border` argument

- refactored `influence.rma.uni()` to be more consistent internally with other functions; `print.infl.rma.uni()` and `plot.infl.rma.uni()` adjusted accordingly; functions `cooks.distance.rma.uni()`, `dfbetas.rma.uni()`, and `rstudent.rma.uni()` now call `influence.rma.uni()` for the computations

- `rstudent.rma.uni()` now computes the SE of the deleted residuals in such a way that it will yield identical results to a mean shift outlier model even when that model is fitted with `test="knha"`

- `rstandard.rma.uni()` gains `type` argument, and can now also compute conditional residuals (it still computes marginal residuals by default)

- `cooks.distance.rma.mv()` gains `cluster` argument, so that the Cook's distances can be computed for groups of estimates

- `cooks.distance.rma.mv()` gains `parallel`, `ncpus`, and `cl` arguments and can now make use of parallel processing

- `cooks.distance.rma.mv()` should be faster by using the estimates from the full model as starting values when fitting the models with the ith study/cluster deleted from the dataset

- `cooks.distance.rma.mv()` gains `reestimate` argument; when set to `FALSE`, variance/correlation components are not reestimated

- `rstandard.rma.mv()` gains `cluster` argument for computing cluster-level multivariate standardized residuals

- added `rstudent.rma.mv()` and `dfbetas.rma.mv()`

- smarter matching of elements in `newmods` (when using a named vector) in `predict()` that also works for models with interactions (thanks to Nicole Erler for pointing out the problem)

- `rma.uni()` and `rma.mv()` no longer issue (obvious) warnings when user constrains `vi` or `V` to 0 (i.e., `vi=0` or `V=0`, respectively)

- `rma.mv()` does more intelligent filtering based on `NA`s in `V` matrix

- `rma.mv()` now ensures strict symmetry of any (var-cov or correlation) matrices specified via the `R` argument

- fixed `rma.mv()` so checks on `R` argument run as intended; also fixed an issue when multiple formulas with slashes are specified via `random` (thanks to Andrew Loignon for pointing out the problem)

- suppressed showing calls on some warnings/errors in `rma.mv()`

- `rma.mv()` now allows for a continuous-time autoregressive random effects structure (`struct="CAR"`) and various spatial correlation structures (`struct="SPEXP"`, `"SPGAU"`, `"SPLIN"`, `"SPRAT"`, and `"SPSPH"`)

- `rma.mv()` now allows for `struct="GEN"` which models correlated random effects for any number of predictors, including continuous ones (i.e., this allows for 'random slopes')

- in the various `forest()` functions, when `options(na.action="na.pass")` or `options(na.action="na.exclude")` and an annotation contains `NA`, this is now shown as a blank (instead of `NA [NA, NA]`)

- the various `forest()` and `addpoly()` functions gain a `fonts` argument

- the various `forest()` functions gain a `top` argument

- the various `forest()` functions now show correct point sizes when the weights of the studies are exactly the same

- `forest.cumul.rma()` gains a `col` argument

- `funnel.default()` and `funnel.rma()` can now take vectors as input for the `col` and `bg` arguments (and also for `pch`); both functions also gain a `legend` argument

- `addpoly()` functions can now also show prediction interval bounds

- removed 'formula interface' from `escalc()`; until this actually adds some kind of extra functionality, this just makes `escalc()` more confusing to use

- `escalc()` can now compute the coefficient of variation ratio and the variability ratio for pre-post or matched designs (`"CVRC"`, `"VRC"`)

- `escalc()` does a bit more housekeeping

- added (currently undocumented) arguments `onlyo1`, `addyi`, and `addvi` to `escalc()` that allow for more flexibility when computing certain bias corrections and when computing sampling variances for measures that make use of the `add` and `to` arguments

- `escalc()` now sets `add=0` for measures where the use of such a bias correction makes little sense; this applies to the following measures: `"AS"`, `"PHI"`, `"RTET"`, `"IRSD"`, `"PAS"`, `"PFT"`, `"IRS"`, and `"IRFT"`; one can still force the use of the bias correction by explicitly setting the `add` argument to some non-zero value

- added `clim` argument to `summary.escalc()`

- added `ilim` argument to `trimfill()`

- `labbe()` gains `lty` argument

- `labbe()` now (invisibly) returns a data frame with the coordinates of the points that were drawn (which may be useful for manual labeling of points in the plot)

- added a print method for `profile.rma` objects

- `profile.rma.mv()` now check whether any of the profiled log-likelihood values is larger than the log-likelihood of the fitted model (using numerical tolerance given by `lltol`) and issues a warning if so

- `profile.rma.uni()`, `profile.rma.mv()`, and `plot.profile.rma()` gain `cline` argument; `plot.profile.rma()` gains `xlim`, `ylab`, and `main` arguments

- fixed an issue with `robust.rma.mv()` when the model was fitted with `sparse=TRUE` (thanks to Roger Martineau for noting the problem)

- various method functions (`fitted()`, `resid()`, `predict()`, etc.) behave in a more consistent manner when model omitted studies with missings

- `predict.rma()` gains `vcov` argument; when set to `TRUE`, the variance-covariance matrix of the predicted values is also returned

- `vcov.rma()` can now also return the variance-covariance matrix of the fitted values (`type="fitted"`) and the residuals (`type="resid"`)

- added `$<-` and `as.matrix()` methods for `list.rma` objects

- fixed error in `simulate.rma()` that would generate too many samples for `rma.mv` models

- added undocumented argument `time` to all model fitting functions; if set to `TRUE`, the model fitting time is printed

- added more tests (also for parallel operations); also, all tests updated to use proper tolerances instead of rounding

- reorganized the documentation a bit

# metafor 2.0-0 (2017-06-22)

- added `simulate()` method for `rma` objects; added `MASS` to `Suggests` (since simulating for `rma.mv` objects requires `mvrnorm()` from `MASS`)

- `cooks.distance.rma.mv()` now works properly even when there are missing values in the data

- `residuals()` gains `type` argument and can compute Pearson residuals

- the `newmods` argument in `predict()` can now be a named vector or a matrix/data frame with column names that get properly matched up with the variables in the model

- added `ranef.rma.mv()` for extracting the BLUPs of the random effects for `rma.mv` models

- all functions that repeatedly refit models now have the option to show a progress bar

- added `ranktest.default()`, so user can now pass the outcomes and corresponding sampling variances directly to the function

- added `regtest.default()`, so user can now pass the outcomes and corresponding sampling variances directly to the function

- `funnel.default()` gains `subset` argument

- `funnel.default()` and `funnel.rma()` gain `col` and `bg` arguments

- `plot.profile.rma()` gains `ylab` argument

- more consistent handling of `robust.rma` objects

- added a print method for `rma.gosh` objects

- the (log) relative risk is now called the (log) risk ratio in all help files, plots, code, and comments

- `escalc()` can now compute outcome measures based on paired binary data (`"MPRR"`, `"MPOR"`, `"MPRD"`, `"MPORC"`, and `"MPPETO"`)

- `escalc()` can now compute (semi-)partial correlation coefficients (`"PCOR"`, `"ZPCOR"`, `"SPCOR"`)

- `escalc()` can now compute measures of variability for single groups (`"CVLN"`, `"SDLN"`) and for the difference in variability between two groups (`"CVR"`, `"VR"`); also the log transformed mean (`"MNLN"`) has been added for consistency

- `escalc()` can now compute the sampling variance for `measure="PHI"` for studies using stratified sampling (`vtpye="ST"`)

- the `[` method for `escalc` objects now properly handles the `ni` and `slab` attributes and does a better job of cleaning out superfluous variable name information

- added `rbind()` method for `escalc` objects

- added `as.data.frame()` method for `list.rma` objects

- added a new dataset (`dat.pagliaro1992`) for another illustration of a network meta-analysis

- added a new dataset (`dat.laopaiboon2015`) on the effectiveness of azithromycin for treating lower respiratory tract infections

- `rma.uni()` and `rma.mv()` now check if the ratio of the largest to smallest sampling variance is very large; results may not be stable then (and very large ratios typically indicate wrongly coded data)

- model fitting functions now check if extra/superfluous arguments are specified via `...` and issues are warning if so

- instead of defining own generic `ranef()`, import `ranef()` from `nlme`

- improved output formatting

- added more tests (but disabled a few tests on CRAN to avoid some issues when R is compiled with `--disable-long-double`)

- some general code cleanup

- renamed `diagram_metafor.pdf` vignette to just `diagram.pdf`

- minor updates in the documentation

# metafor 1.9-9 (2016-09-25)

- started to use git as version control system, GitHub to host the repository (https://github.com/wviechtb/metafor) for the development version of the package, Travis CI as continuous integration service (https://travis-ci.org/wviechtb/metafor), and Codecov for automated code coverage reporting (https://app.codecov.io/gh/wviechtb/metafor)

- argument `knha` in `rma.uni()` and argument `tdist` in `rma.glmm()` and `rma.mv()` are now superseded by argument `test` in all three functions; for backwards compatibility, the `knha` and `tdist` arguments still work, but are no longer documented

- `rma(yi, vi, weights=1, test="knha")` now yields the same results as `rma(yi, vi, weighted=FALSE, test="knha")` (but use of the Knapp and Hartung method in the context of an unweighted analysis remains an experimental feature)

- one can now pass an `escalc` object directly to `rma.uni()`, which then tries to automatically determine the `yi` and `vi` variables in the data frame (thanks to Christian Roever for the suggestion)

- `escalc()` can now also be used to convert a regular data frame to an `escalc` object

- for `measure="UCOR"`, the exact bias-correction is now used (instead of the approximation); when `vtype="UB"`, the exact equation is now used to compute the unbiased estimate of the variance of the bias-corrected correlation coefficient; hence `gsl` is now a suggested package (needed to compute the hypergeometric function) and is loaded when required

- `cooks.distance()` now also works with `rma.mv` objects; and since model fitting can take some time, an option to show a progress bar has been added

- fixed an issue with `robust.rma.mv()` throwing errors when the model was fitted with `sparse=TRUE`

- fixed an error with `robust.rma.mv()` when the model was fitted with user-defined weights (or a user-defined weight matrix)

- added `ranef()` for extracting the BLUPs of the random effects (only for `rma.uni` objects at the moment)

- reverted back to the pre-1.1-0 way of computing p-values for individual coefficients in `permutest.rma.uni()`, that is, the p-value is computed with `mean(abs(z_perm) >= abs(z_obs) - tol)` (where `tol` is a numerical tolerance)

- `permutest.rma.uni()` gains `permci` argument, which can be used to obtain permutation-based CIs of the model coefficients (note that this is computationally very demanding and may take a long time to complete)

- `rma.glmm()` continues to work even when the saturated model cannot be fitted (although the tests for heterogeneity are not available then)

- `rma.glmm()` now allows control over the arguments used for `method.args` (via `control=list(hessianCtrl=list(...))`) passed to `hessian()` (from the `numDeriv` package) when using `model="CM.EL"` and `measure="OR"`

- in `rma.glmm()`, default `method.args` value for `r` passed to `hessian()` has been increased to 16 (while this slows things down a bit, this appears to improve the accuracy of the numerical approximation to the Hessian, especially when tau^2 is close to 0)

- the various `forest()` and `addpoly()` functions now have a new argument called `width`, which provides manual control over the width of the annotation columns; this is useful when creating complex forest plots with a monospaced font and we want to ensure that all annotations are properly lined up at the decimal point

- the annotations created by the various `forest()` and `addpoly()` functions are now a bit more compact by default

- more flexible `efac` argument in the various `forest()` functions

- trailing zeros in the axis labels are now dropped in forest and funnel plots by default; but trailing zeros can be retained by specifying a numeric (and not an integer) value for the `digits` argument

- added `funnel.default()`, which directly takes as input a vector with the observed effect sizes or outcomes and the corresponding sampling variances, standard errors, and/or sample sizes

- added `plot.profile.rma()`, a plot method for objects returned by the `profile.rma.uni()` and `profile.rma.mv()` functions

- simplified `baujat.rma.uni()`, `baujat.rma.mh()`, and `baujat.rma.peto()` to `baujat.rma()`, which now handles objects of class `rma.uni`, `rma.mh`, and `rma.peto`

- `baujat.rma()` gains argument `symbol` for more control over the plotting symbol

- `labbe()` gains a `grid` argument

- more logical placement of labels in `qqnorm.rma.uni()`, `qqnorm.rma.mh()`, and `qqnorm.rma.peto()` functions (and more control thereof)

- `qqnorm.rma.uni()` gains `lty` argument

- added `gosh.rma()` and `plot.gosh.rma()` for creating GOSH (i.e., graphical display of study heterogeneity) plots based on Olkin et al. (2012)

- in the (rare) case where all observed outcomes are exactly equal to each other, `test="knha"` (i.e., `knha=TRUE`) in `rma()` now leads to more appropriate results

- updated datasets so those containing precomputed effect size estimates or observed outcomes are already declared to be `escalc` objects

- added new datasets (`dat.egger2001` and `dat.li2007`) on the effectiveness of intravenous magnesium in acute myocardial infarction

- `methods` package is now under `Depends` (in addition to `Matrix`), so that `rma.mv(..., sparse=TRUE)` always works, even under Rscript

- some general code cleanup

- added more tests (and used a more consistent naming scheme for tests)

# metafor 1.9-8 (2015-09-28)

- due to more stringent package testing, it is increasingly difficult to ensure that the package passes all checks on older versions of R; from now on, the package will therefore require, and be checked under, only the current (and the development) version of R

- added `graphics`, `grDevices`, and `methods` to `Imports` (due to recent change in how CRAN checks packages)

- the `struct` argument for `rma.mv()` now also allows for `"ID"` and `"DIAG"`, which are identical to the `"CS"` and `"HCS"` structures, but with the correlation parameter fixed to 0

- added `robust()` for (cluster) robust tests and confidence intervals for `rma.uni` and `rma.mv` models (this uses a robust sandwich-type estimator of the variance-covariance matrix of the fixed effects along the lines of the Eicker-Huber-White method)

- `confint()` now works for models fitted with the `rma.mv()` function; for variance and correlation parameters, the function provides profile likelihood confidence intervals; the output generated by the `confint()` function has been adjusted in general to make the formatting more consistent across the different model types

- for objects of class `rma.mv`, `profile()` now provides profile plots for all (non-fixed) variance and correlation components of the model when no component is specified by the user (via the `sigma2`, `tau2`, `rho`, `gamma2`, or `phi` arguments)

- for `measure="MD"` and `measure="ROM"`, one can now choose between `vtype="LS"` (the default) and `vtype="HO"`; the former computes the sampling variances without assuming homoscedasticity, while the latter assumes homoscedasticity

- multiple model objects can now be passed to the `fitstats()`, `AIC()`, and `BIC()` functions

- check for duplicates in the `slab` argument is now done *after* any subsetting is done (as suggested by Michael Dewey)

- `rma.glmm()` now again works when using `add=0`, in which case some of the observed outcomes (e.g., log odds or log odds ratios) may be `NA`

- when using `rma.glmm()` with `model="CM.EL"`, the saturated model (used to compute the Wald-type and likelihood ratio tests for the presence of (residual) heterogeneity) often fails to converge; the function now continues to run (instead of stopping with an error) and simply omits the test results from the output

- when using `rma.glmm()` with `model="CM.EL"` and inversion of the Hessian fails via the Choleski factorization, the function now makes another attempt via the QR decomposition (even when this works, a warning is issued)

- for `rma.glmm()`, BIC and AICc values were switched around; corrected

- more use of `suppressWarnings()` is made when functions repeatedly need to fit the same model, such as `cumul()`, `influence()`, and `profile()`; that way, one does not get inundated with the same warning(s)

- some (overdue) updates to the documentation

# metafor 1.9-7 (2015-05-22)

- default optimizer for `rma.mv()` changed to `nlminb()` (instead of `optim()` with `"Nelder-Mead"`); extensive testing indicated that `nlminb()` (and also `optim()` with `"BFGS"`) is typically quicker and more robust; note that this is in principle a non-backwards compatible change, but really a necessary one; and you can always revert to the old behavior with `control=list(optimizer="optim", optmethod="Nelder-Mead")`

- all tests have been updated in accordance with the recommended syntax of the `testthat` package; for example, `expect_equivalent(x,y)` is used instead of `test_that(x, is_equivalent_to(y))`

- changed a few `is_identical_to()` comparisons to `expect_equivalent()` ones (that failed on Sparc Solaris)

# metafor 1.9-6 (2015-05-07)

- `funnel()` now works again for `rma.glmm` objects (note to self: quit breaking things that work!)

- `rma.glmm()` will now only issue a warning (and not an error) when the Hessian for the saturated model cannot be inverted (which is needed to compute the Wald-type test for heterogeneity, so the test statistic is then simply set to `NA`)

- `rma.mv()` now allows for two terms of the form `~ inner | outer`; the variance components corresponding to such a structure are called `gamma2` and correlations are called `phi`; other functions that work with objects of class `rma.mv` have been updated accordingly

- `rma.mv()` now provides (even) more optimizer choices: `nlm()` from the `stats` package, `hjk()` and `nmk()` from the `dfoptim` package, and `ucminf()` from the `ucminf` package; choose the desired optimizer via the control argument (e.g., `control=list(optimizer="nlm")`)

- `profile.rma.uni()` and `profile.rma.mv()` now can do parallel processing (which is especially relevant for `rma.mv` objects, where profiling is crucial and model fitting can be slow)

- the various `confint()` functions now have a `transf` argument (to apply some kind of transformation to the model coefficients and confidence interval bounds); coefficients and bounds for objects of class `rma.mh` and `rma.peto` are no longer automatically transformed

- the various `forest()` functions no longer enforce that the actual x-axis limits (`alim`) encompass the observed outcomes to be plotted; also, outcomes below or above the actual x-axis limits are no longer shown

- the various `forest()` functions now provide control over the horizontal lines (at the top/bottom) that are automatically added to the plot via the `lty` argument (this also allows for removing them); also, the vertical reference line is now placed *behind* the points/CIs

- `forest.default()` now has argument `col` which can be used to specify the color(s) to be used for drawing the study labels, points, CIs, and annotations

- the `efac` argument for `forest.rma()` now also allows two values, the first for the arrows and CI limits, the second for summary estimates

- corrected some axis labels in various plots when `measure="PLO"`

- axes in `labbe()` plots now have `"(Group 1)"` and `"(Group 2)"` added by default

- `anova.rma()` gains argument `L` for specifying linear combinations of the coefficients in the model that should be tested to be zero

- in case removal of a row of data would lead to one or more inestimable model coefficients, `baujat()`, `cooks.distance()`, `dfbetas()`, `influence()`, and `rstudent()` could fail for `rma.uni` objects; such cases are now handled properly

- for models with moderators, the `predict()` function now shows the study labels when they have been specified by the user (and `newmods` is not used)

- if there is only one fixed effect (model coefficient) in the model, the `print.infl.rma.uni()` function now shows the DFBETAS values with the other case diagnostics in a single table (for easier inspection); if there is more than one fixed effect, a separate table is still used for the DFBETAS values (with one column for each coefficient)

- added `measure="SMCRH"` to the `escalc()` function for the standardized mean change using raw score standardization with heteroscedastic population variances at the two measurement occasions

- added `measure="ROMC"` to the `escalc()` function for the (log transformed) ratio of means (response ratio) when the means reflect two measurement occasions (e.g., for a single group of people) and hence are correlated

- added own function for computing/estimating the tetrachoric correlation coefficient (for `measure="RTET"`); package therefore no longer suggests `polycor` but now suggest `mvtnorm` (which is loaded as needed)

- element `fill` returned by `trimfill.rma.uni()` is now a logical vector (instead of a 0/1 dummy variable)

- `print.list.rma()` now also returns the printed results invisibly as a data frame

- added a new dataset (`dat.senn2013`) as another illustration of a network meta-analysis

- `metafor` now depends on at least version 3.1.0 of R

# metafor 1.9-5 (2014-11-24)

- moved the `stats` and `Matrix` packages from `Depends` to `Imports`; as a result, had to add `utils` to `Imports`; moved the `Formula` package from `Depends` to `Suggests`

- added `update.rma()` function (for updating/refitting a model); model objects also now store and keep the call

- the `vcov()` function now also extracts the marginal variance-covariance matrix of the observed effect sizes or outcomes from a fitted model (of class `rma.uni` or `rma.mv`)

- `rma.mv()` now makes use of the Cholesky decomposition when there is a `random = ~ inner | outer` formula and `struct="UN"`; this is numerically more stable than the old approach that avoided non-positive definite solutions by forcing the log-likelihood to be -Inf in those cases; the old behavior can be restored with `control = list(cholesky=FALSE)`

- `rma.mv()` now requires the `inner` variable in an `~ inner | outer` formula to be a factor or character variable (except when `struct` is `"AR"` or `"HAR"`); use `~ factor(inner) | outer` in case it isn't

- `anova.rma.uni()` function changed to `anova.rma()` that works now for both `rma.uni` and `rma.mv` objects

- the `profile.rma.mv()` function now omits the number of the variance or correlation component from the plot title and x-axis label when the model only includes one of the respective parameters

- `profile()` functions now pass on the `...` argument also to the `title()` function used to create the figure titles (esp. relevant when using the `cex.main` argument)

- the `drop00` argument of the `rma.mh()` and `rma.peto()` functions now also accepts a vector with two logicals, the first applies when calculating the observed outcomes, the second when applying the Mantel-Haenszel or Peto's method

- `weights.rma.uni()` now shows the correct weights when `weighted=FALSE`

- argument `showweight` renamed to `showweights` in the `forest.default()` and `forest.rma()` functions (more consistent with the naming of the various `weights()` functions)

- added `model.matrix.rma()` function (to extract the model matrix from objects of class `rma`)

- `funnel()` and `radial()` now (invisibly) return data frames with the coordinates of the points that were drawn (may be useful for manual labeling of points in the plots)

- `permutest.rma.uni()` function now uses a numerical tolerance when making comparisons (>= or <=) between an observed test statistic and the test statistic under the permuted data; when using random permutations, the function now ensures that the very first permutation correspond to the original data

- corrected some missing/redundant row/column labels in some output

- most `require()` calls replaced with `requireNamespace()` to avoid altering the search path (hopefully this won't break stuff ...)

- some non-visible changes including more use of some (non-exported) helper functions for common tasks

- dataset `dat.collins91985a` updated (including all reported outcomes and some more information about the various trials)

- oh, and guess what? I updated the documentation ...

# metafor 1.9-4 (2014-07-30)

- added `method="GENQ"` to `rma.uni()` for the generalized Q-statistic estimator of tau^2, which allows for used-defined weights (note: the DL and HE estimators are just special cases of this method)

- when the model was fitted with `method="GENQ"`, then `confint()` will now use the generalized Q-statistic method to construct the corresponding confidence interval for tau^2 (thanks to Dan Jackson for the code); the iterative method used to obtain the CI makes use of Farebrother's algorithm as implemented in the `CompQuadForm` package

- slight improvements in how the `rma.uni()` function handles non-positive sampling variances

- `rma.uni()`, `rma.mv()`, and `rma.glmm()` now try to detect and remove any redundant predictors before the model fitting; therefore, if there are exact linear relationships among the predictor variables (i.e., perfect multicollinearity), terms are removed to obtain a set of predictors that is no longer perfectly multicollinear (a warning is issued when this happens); note that the order of how the variables are specified in the model formula can influence which terms are removed

- the last update introduced an error in how hat values were computed when the model was fitted with the `rma()` function using the Knapp & Hartung method (i.e., when `knha=TRUE`); this has been fixed

- `regtest()` no longer works (for now) with `rma.mv` objects (it wasn't meant to in the first place); if you want to run something along the same lines, just consider adding some measure of the precision of the observed outcomes (e.g., their standard errors) as a predictor to the model

- added `"sqrtni"` and `"sqrtninv"` as possible options for the `predictor` argument of `regtest()`

- more optimizers are now available for the `rma.mv()` function via the `nloptr` package by setting `control = list(optimizer="nloptr")`; when using this optimizer, the default is to use the BOBYQA implementation from that package with a relative convergence criterion of 1e-8 on the function value (see documentation on how to change these defaults)

- `predict.rma()` function now works for `rma.mv` objects with multiple tau^2 values even if the user specifies the `newmods` argument but not the `tau2.levels` argument (but a warning is issued and the prediction intervals are not computed)

- argument `var.names` now works properly in `escalc()` when the user has not made use of the `data` argument (thanks to Jarrett Byrnes for bringing this to my attention)

- added `plot()` function for cumulative random-effects models results as obtained with the `cumul.rma.uni()` function; the plot shows the model estimate on the x-axis and the corresponding tau^2 estimate on the y-axis in the cumulative order of the results

- fixed the omitted offset term in the underlying model fitted by the `rma.glmm()` function when `method="ML"`, `measure="IRR"`, and `model="UM.FS"`, that is, when fitting a mixed-effects Poisson regression model with fixed study effects to two-group event count data (thanks to Peter Konings for pointing out this error)

- added two new datasets (`dat.bourassa1996`, `dat.riley2003`)

- added function `replmiss()` (just a useful helper function)

- package now uses `LazyData: TRUE`

- some improvements to the documentation (do I still need to mention this every time?)

# metafor 1.9-3 (2014-05-05)

- some minor tweaks to `rma.uni()` that should be user transparent

- `rma.uni()` now has a `weights` argument, allowing the user to specify arbitrary user-defined weights; all functions affected by this have been updated accordingly

- better handling of mismatched length of `yi` and `ni` vectors in `rma.uni()` and `rma.mv()` functions

- subsetting is now handled as early as possible within functions with subsetting capabilities; this avoids some (rare) cases where studies ultimately excluded by the subsetting could still affect the results

- some general tweaks to `rma.mv()` that should make it a bit faster

- argument `V` of `rma.mv()` now also accepts a list of var-cov matrices for the observed effects or outcomes; from the list elements, the full (block diagonal) var-cov matrix `V` is then automatically constructed

- `rma.mv()` now has a new argument `W` allowing the user to specify arbitrary user-defined weights or an arbitrary weight matrix

- `rma.mv()` now has a new argument `sparse`; by setting this to `TRUE`, the function uses sparse matrix objects to the extent possible; this can speed up model fitting substantially for certain models (hence, the `metafor` package now depends on the `Matrix` package)

- `rma.mv()` now allows for `struct="AR"` and `struct="HAR"`, to fit models with (heteroscedastic) autoregressive (AR1) structures among the true effects (useful for meta-analyses of studies reporting outcomes at multiple time points)

- `rma.mv()` now has a new argument `Rscale` which can be used to control how matrices specified via the `R` argument are scaled (see docs for more details)

- `rma.mv()` now only checks for missing values in the rows of the lower triangular part of the `V` matrix (including the diagonal); this way, if `Vi = matrix(c(.5,NA,NA,NA), nrow=2, ncol=2)` is the var-cov matrix of the sampling errors for a particular study with two outcomes, then only the second row/column needs to be removed before the model fitting (and not the entire study)

- added five new datasets (`dat.begg1989`, `dat.ishak2007`, `dat.fine1993`, `dat.konstantopoulos2011`, and `dat.hasselblad1998`) to provide further illustrations of the use of the `rma.mv()` function (for meta-analyses combining controlled and uncontrolled studies, for meta-analyses of longitudinal studies, for multilevel meta-analyses, and for network meta-analyses / mixed treatment comparison meta-analyses)

- added `rstandard.rma.mv()` function to compute standardized residuals for models fitted with the `rma.mv()` function (`rstudent.rma.mv()` to be added at a later point); also added `hatvalues.rma.mv()` for computing the hat values and `weights.rma.uni()` for computing the weights (i.e., the diagonal elements of the weight matrix)

- the various `weights()` functions now have a new argument `type` to indicate whether only the diagonal elements of the weight matrix (default) or the entire weight matrix should be returned

- the various `hatvalues()` functions now have a new argument `type` to indicate whether only the diagonal elements of the hat matrix (default) or the entire hat matrix should be returned

- `predict.rma()` function now works properly for `rma.mv` objects (also has a new argument `tau2.levels` to specify, where applicable, the levels of the inner factor when computing prediction intervals)

- `forest.rma()` function now provides a bit more control over the color of the summary polygon and is now compatible with `rma.mv` objects; also, has a new argument `lty`, which provides more control over the line type for the individual CIs and the prediction interval

- `addpoly.default()` and `addpoly.rma()` now have a `border` argument (for consistency with the `forest.rma()` function); `addpoly.rma()` now yields the correct CI bounds when the model was fitted with `knha=TRUE`

- `forest.cumul.rma()` now provides the correct CI bounds when the models were fitted with the Knapp & Hartung method (i.e., when `knha=TRUE` in the original `rma()` function call)

- the various `forest()` functions now return information about the chosen values for arguments `xlim`, `alim`, `at`, `ylim`, `rows`, `cex`, `cex.lab`, and `cex.axis` invisibly (useful for tweaking the default values); thanks to Michael Dewey for the suggestion

- the various `forest()` functions now have a new argument, `clim`, to set limits for the confidence/prediction interval bounds

- `cumul.mh()` and `cumul.peto()` now get the order of the studies right when there are missing values in the data

- the `transf` argument of `leave1out.rma.mh()`, `leave1out.rma.peto()`, `cumul.rma.mh()`, and `cumul.rma.peto()` should now be used to specify the actual function for the transformation (the former behavior of setting this argument to `TRUE` to exponentiate log RRs, log ORs, or log IRRs still works for back-compatibility); this is more consistent with how the `cumul.rma.uni()` and `leave1out.rma.uni()` functions work and is also more flexible

- added `bldiag()` function to construct a block diagonal matrix from (a list of) matrices (may be needed to construct the `V` matrix when using the `rma.mv()` function); `bdiag()` function from the `Matrix` package does the same thing, but creates sparse matrix objects

- `profile.rma.mv()` now has a `startmethod` argument; by setting this to `"prev"`, successive model fits are started at the parameter estimates from the previous model fit; this may speed things up a bit; also, the method for automatically choosing the `xlim` values has been changed

- slight improvement to `profile.rma.mv()` function, which would throw an error if the last model fit did not converge

- added a new dataset (`dat.linde2005`) for replication of the analyses in Viechtbauer (2007)

- added a new dataset (`dat.molloy2014`) for illustrating the meta-analysis of (r-to-z transformed) correlation coefficients

- added a new dataset (`dat.gibson2002`) to illustrate the combined analysis of standardized mean differences and probit transformed risk differences

- computations in `weights.mh()` slightly changed to prevent integer overflows for large counts

- unnecessary warnings in `transf.ipft.hm()` are now suppressed (cases that raised those warnings were already handled correctly)

- in `predict()`, `blup()`, `cumul()`, and `leave1out()`, when using the `transf` argument, the standard errors (which are `NA`) are no longer shown in the output

- argument `slab` in various functions will now also accept non-unique study labels; `make.unique()` is used as needed to make them unique

- `vignettes("metafor")` and `vignettes("metafor_diagram")` work again (yes, I know they are not true vignettes in the strict sense, but I think they should show up on the CRAN website for the package and using a minimal valid Sweave document that is recognized by the R build system makes that happen)

- `escalc()` and its `summary()` method now keep better track when the data frame contains multiple columns with outcome or effect size values (and corresponding sampling variances) for print formatting; also simplified the class structure a bit (and hence, `print.summary.escalc()` removed)

- `summary.escalc()` has a new argument `H0` to specify the value of the outcome under the null hypothesis for computing the test statistics

- added measures `"OR2DN"` and `"D2ORN"` to `escalc()` for transforming log odds ratios to standardized mean differences and vice-versa, based on the method of Cox & Snell (1989), which assumes normally distributed response variables within the two groups before the dichotomization

- `permutest.rma.uni()` function now catches an error when the number of permutations requested is too large (for R to even create the objects to store the results in) and produces a proper error message

- `funnel.rma()` function now allows the `yaxis` argument to be set to `"wi"` so that the actual weights (in %) are placed on the y-axis (useful when arbitrary user-defined have been specified)

- for `rma.glmm()`, the control argument `optCtrl` is now used for passing control arguments to all of the optimizers (hence, control arguments `nlminbCtrl` and `minqaCtrl` are now defunct)

- `rma.glmm()` should not throw an error anymore when including only a single moderator/predictor in the model

- `predict.rma()` now returns an object of class `list.rma` (therefore, function `print.predict.rma()` has been removed)

- for `rma.list` objects, added `[`, `head()`, and `tail()` methods

- automated testing using the `testthat` package (still many more tests to add, but finally made a start on this)

- encoding changed to UTF-8 (to use 'foreign characters' in the docs and to make the HTML help files look a bit nicer)

- guess what? some improvements to the documentation! (also combined some of the help files to reduce the size of the manual a bit; and yes, it's still way too big)

# metafor 1.9-2 (2013-10-07)

- added function `rma.mv()` to fit multivariate/multilevel meta-analytic models via appropriate linear (mixed-effects) models; this function allows for modeling of non-independent sampling errors and/or true effects and can be used for network meta-analyses, meta-analyses accounting for phylogenetic relatedness, and other complicated meta-analytic data structures

- added the AICc to the information criteria computed by the various model fitting functions

- if the value of tau^2 is fixed by the user via the corresponding argument in `rma.uni()`, then tau^2 is no longer counted as an additional parameter for the computation of the information criteria (i.e., AIC, BIC, and AICc)

- `rma.uni()`, `rma.glmm()`, and `rma.mv()` now use a more stringent check whether the model matrix is of full rank

- added `profile()` method functions for objects of class `rma.uni` and `rma.mv` (can be used to obtain a plot of the profiled log-likelihood as a function of a specific variance component or correlation parameter of the model)

- `predict.rma()` function now has an `intercept` argument that allows the user to decide whether the intercept term should be included when calculating the predicted values (rare that this should be changed from the default)

- for `rma.uni()`, `rma.glmm()`, and `rma.mv()`, the `control` argument can now also accept an integer value; values > 1 generate more verbose output about the progress inside of the function

- `rma.glmm()` has been updated to work with `lme4` 1.0.x for fitting various models; as a result, `model="UM.RS"` can only use `nAGQ=1` at the moment (hopefully this will change in the future)

- the `control` argument of `rma.glmm()` can now be used to pass all desired control arguments to the various functions and optimizers used for the model fitting (admittedly the use of lists within this argument is a bit unwieldy, but much more flexible)

- `rma.mh()` and `rma.peto()` also now have a `verbose` argument (not really needed, but added for sake of consistency across functions)

- fixed (silly) error that would prevent `rma.glmm()` from running for measures `"IRR"`, `"PLO"`, and `"IRLN"` when there are missing values in the data (lesson: add some missing values to datasets for the unit tests!)

- a bit of code reorganization (should be user transparent)

- vignettes (`"metafor"` and `"metafor_diagram"`) are now just 'other files' in the doc directory (as these were not true vignettes to begin with)

- some improvements to the documentation (as always)

# metafor 1.9-1 (2013-07-20)

- `rma.mh()` now also implements the Mantel-Haenszel method for incidence rate differences (`measure="IRD"`)

- when analyzing incidence rate ratios (`measure="IRR"`) with the `rma.mh()` function, the Mantel-Haenszel test for person-time data is now also provided

- `rma.mh()` has a new argument `correct` (default is `TRUE`) to indicate whether the continuity correction should be applied when computing the (Cochran-)Mantel-Haenszel test statistic

- renamed elements `CMH` and `CMHp` (for the Cochran-Mantel-Haenszel test statistic and corresponding p-value) to `MH` and `MHp`

- added function `baujat()` to create Baujat plots

- added a new dataset (`dat.pignon2000`) to illustrate the use of the `baujat()` function

- added function `to.table()` to convert data from vector format into the corresponding table format

- added function `to.long()` to convert data from vector format into the corresponding long format

- `rma.glmm()` now even runs when k=1 (yielding trivial results)

- for models with an intercept and moderators, `rma.glmm()` now internally rescales (non-dummy) variables to z-scores during the model fitting (this improves the stability of the model fitting, especially when `model="CM.EL"`); results are given after back-scaling, so this should be transparent to the user

- in `rma.glmm()`, default number of quadrature points (`nAGQ`) is now 7 (setting this to 100 was a bit overkill)

- a few more error checks here and there for misspecified arguments

- some improvements to the documentation

# metafor 1.9-0 (2013-06-21)

- vignette renamed to `metafor` so `vignette("metafor")` works now

- added a diagram to the documentation, showing the various functions in the `metafor` package (and how they relate to each other); can be loaded with `vignette("metafor_diagram")`

- `anova.rma.uni()` function can now also be used to test (sub)sets of model coefficients with a Wald-type test when a single model is passed to the function

- the pseudo R^2 statistic is now automatically calculated by the `rma.uni()` function and supplied in the output (only for mixed-effects models and when the model includes an intercept, so that the random- effects model is clearly nested within the mixed-effects model)

- component `VAF` is now called `R2` in `anova.rma.uni()` function

- added function `hc()` that carries out a random-effects model analysis using the method by Henmi and Copas (2010); thanks to Michael Dewey for the suggestion and providing the code

- added new dataset (`dat.lee2004`), which was used in the article by Henmi and Copas (2010) to illustrate their method

- fixed missing x-axis labels in the `forest()` functions

- `rma.glmm()` now computes Hessian matrices via the `numDeriv` package when `model="CM.EL"` and `measure="OR"` (i.e., for the conditional logistic model with exact likelihood); so `numDeriv` is now a suggested package and is loaded within `rma.glmm()` when required

- `trimfill.rma.uni()` now also implements the `"Q0"` estimator (although the `"L0"` and `"R0"` estimators are generally to be preferred)

- `trimfill.rma.uni()` now also calculates the SE of the estimated number of missing studies and, for estimator `"R0"`, provides a formal test of the null hypothesis that the number of missing studies on a given side is zero

- added new dataset (`dat.bangertdrowns2004`)

- the `level` argument in various functions now either accepts a value representing a percentage or a proportion (values greater than 1 are assumed to be a percentage)

- `summary.escalc()` now computes confidence intervals correctly when using the `transf` argument

- computation of Cochran-Mantel-Haenszel statistic in `rma.mh()` changed slightly to avoid integer overflow with very big counts

- some internal improvements with respect to object attributes that were getting discarded when subsetting

- some general code cleanup

- some improvements to the documentation

# metafor 1.8-0 (2013-04-11)

- added additional clarifications about the change score outcome measures (`"MC"`, `"SMCC"`, and `"SMCR"`) to the help file for the `escalc()` function and changed the code so that `"SMCR"` no longer expects argument `sd2i` to be specified (which is not needed anyways) (thanks to Markus Kösters for bringing this to my attention)

- sampling variance for the biserial correlation coefficient (`"RBIS"`) is now calculated in a slightly more accurate way

- `llplot()` now properly scales the log-likelihoods

- argument `which` in the `plot.infl.rma.uni()` function has been replaced with argument `plotinf` which can now also be set to `FALSE` to suppress plotting of the various case diagnostics altogether

- labeling of the axes in `labbe()` plots is now correct for odds ratios (and transformations thereof)

- added two new datasets (`dat.nielweise2007` and `dat.nielweise2008`) to illustrate some methods/models from the `rma.glmm()` function

- added a new dataset (`dat.yusuf1985`) to illustrate the use of `rma.peto()`

- test for heterogeneity is now conducted by the `rma.peto()` function exactly as described by Yusuf et al. (1985)

- in `rma.glmm()`, default number of quadrature points (`nAGQ`) is now 100 (which is quite a bit slower, but should provide more than sufficient accuracy in most cases)

- the standard errors of the HS and DL estimators of tau^2 are now correctly computed when tau^2 is prespecified by the user in the `rma()` function; in addition, the standard error of the SJ estimator is also now provided when tau^2 is prespecified

- `rma.uni()` and `rma.glmm()` now use a better method to check whether the model matrix is of full rank

- I^2 and H^2 statistics are now also calculated for mixed-effects models by the `rma.uni()` and `rma.glmm()` function; `confint.rma.uni()` provides the corresponding confidence intervals for `rma.uni` models

- various `print()` methods now have a new argument called `signif.stars`, which defaults to `getOption("show.signif.stars")` (which by default is `TRUE`) to determine whether the infamous 'significance stars' should be printed

- slight changes in wording in the output produced by the `print.rma.uni()` and `print.rma.glmm()` functions

- some improvements to the documentation

# metafor 1.7-0 (2013-02-06)

- added `rma.glmm()` function for fitting of appropriate generalized linear (mixed-effects) models when analyzing odds ratios, incidence rate ratios, proportions, or rates; the function makes use of the `lme4` and `BiasedUrn` packages; these are now suggested packages and loaded within `rma.glmm()` only when required (this makes for faster loading of the `metafor` package)

- added several method functions for objects of class `rma.glmm` (not all methods yet implemented; to be completed in the future)

- `rma.uni()` now allows the user to specify a formula for the `yi` argument, so instead of rma(yi, vi, mods=~mod1+mod2), one can specify the same model with rma(yi~mod1+mod2, vi)

- `rma.uni()` now has a `weights` argument to specify the inverse of the sampling variances (instead of using the `vi` or `sei` arguments); for now, this is all this argument should be used for (in the future, this argument may potentially be used to allow the user to define alternative weights)

- `rma.uni()` now checks whether the model matrix is not of full rank and issues an error accordingly (instead of the rather cryptic error that was issued before)

- `rma.uni()` now has a `verbose` argument

- `coef.rma()` now returns only the model coefficients (this change was necessary to make the package compatible with the `multcomp` package; see `help(rma)` for an example); use `coef(summary())` to obtain the full table of results

- the `escalc()` function now does some more extensive error checking for misspecified data and some unusual cases

- `append` argument is now `TRUE` by default in the `escalc()` function

- objects generated by the `escalc()` function now have their own class

- added `print()` and `summary()` methods for objects of class `escalc`

- added `[` and `cbind()` methods for objects of class `escalc`

- added a few additional arguments to the `escalc()` function (i.e., `slab`, `subset`, `var.names`, `replace`, `digits`)

- added `drop00` argument to the `escalc()`, `rma.uni()`, `rma.mh()`, and `rma.peto()` functions

- added `"MN"`, `"MC"`, `"SMCC"`, and `"SMCR"` measures to the `escalc()` and `rma.uni()` functions for the raw mean, the raw mean change, and the standardized mean change (with change score or raw score standardization) as possible outcome measures

- the `"IRFT"` measure in the `escalc()` and `rma.uni()` functions is now computed with `1/2*(sqrt(xi/ti) + sqrt(xi/ti+1/ti))` which is more consistent with the definition of the Freeman-Tukey transformation for proportions

- added `"RTET"` measure to the `escalc()` and `rma.uni()` functions to compute the tetrachoric correlation coefficient based on 2x2 table data (the `polycor` package is therefore now a suggested package, which is loaded within `escalc()` only when required)

- added `"RPB"` and `"RBIS"` measures to the `escalc()` and `rma.uni()` functions to compute the point-biserial and biserial correlation coefficient based on means and standard deviations

- added `"PBIT"` and `"OR2D"` measures to the `escalc()` and `rma.uni()` functions to compute the standardized mean difference based on 2x2 table data

- added the `"D2OR"` measure to the `escalc()` and `rma.uni()` functions to compute the log odds ratio based on the standardized mean difference

- added `"SMDH"` measure to the `escalc()` and `rma.uni()` functions to compute the standardized mean difference without assuming equal population variances

- added `"ARAW"`, `"AHW"`, and `"ABT"` measures to the `escalc()` and `rma.uni()` functions for the raw value of Cronbach's alpha, the transformation suggested by Hakstian & Whalen (1976), and the transformation suggested by Bonett (2002) for the meta-analysis of reliability coefficients (see `help(escalc)` for details)

- corrected a small mistake in the equation used to compute the sampling variance of the phi coefficient (`measure="PHI"`) in the `escalc()` function

- the `permutest.rma.uni()` function now uses an algorithm to find only the unique permutations of the model matrix (which may be much smaller than the total number of permutations), making the exact permutation test feasible in a larger set of circumstances (thanks to John Hodgson for making me aware of this issue and to Hans-Jörg Viechtbauer for coming up with a recursive algorithm for finding the unique permutations)

- prediction interval in `forest.rma()` is now indicated with a dotted (instead of a dashed) line; ends of the interval are now marked with vertical bars

- completely rewrote the `funnel.rma()` function which now supports many more options for the values to put on the y-axis; `trimfill.rma.uni()` function was adapted accordingly

- removed the `ni` argument from the `regtest.rma()` function; instead, sample sizes can now be explicitly specified via the `ni` argument when using the `rma.uni()` function (i.e., when `measure="GEN"`); the `escalc()` function also now adds information on the `ni` values to the resulting data frame (as an attribute of the `yi` variable), so, if possible, this information is passed on to `regtest.rma()`

- added switch so that `regtest()` can also provide the full results from the fitted model (thanks to Michael Dewey for the suggestion)

- `weights.rma.mh()` now shows the weights in % as intended (thanks to Gavin Stewart for pointing out this error)

- more flexible handling of the `digits` argument in the various forest functions

- forest functions now use `pretty()` by default to set the x-axis tick locations (`alim` and `at` arguments can still be used for complete control)

- studies that are considered to be 'influential' are now marked with an asterisk when printing the results returned by the `influence.rma.uni()` function (see the documentation of this function for details on how such studies are identified)

- added additional extractor functions for some of the influence measures (i.e., `cooks.distance()`, `dfbetas()`); unfortunately, the `covratio()` and `dffits()` functions in the `stats` package are not generic; so, to avoid masking, there are currently no extractor functions for these measures

- better handling of missing values in some unusual situations

- corrected small bug in `fsn()` that would not allow the user to specify the standard errors instead of the sampling variances (thanks to Bernd Weiss for pointing this out)

- `plot.infl.rma.uni()` function now allows the user to specify which plots to draw (and the layout) and adds the option to show study labels on the x-axis

- added proper `print()` method for objects generated by the `confint.rma.uni()`, `confint.rma.mh()`, and `confint.rma.peto()` functions

- when `transf` or `atransf` argument was a monotonically *decreasing* function, then confidence and prediction interval bounds were in reversed order; various functions now check for this and order the bounds correctly

- `trimfill.rma.uni()` now only prints information about the number of imputed studies when actually printing the model object

- `qqnorm.rma.uni()`, `qqnorm.rma.mh()`, and `qqnorm.rma.peto()` functions now have a new argument called `label`, which allows for labeling of points; the functions also now return (invisibly) the x and y coordinates of the points drawn

- `rma.mh()` with `measure="RD"` now computes the standard error of the estimated risk difference based on Sato, Greenland, & Robins (1989), which provides a consistent estimate under both large-stratum and sparse-data limiting models

- the restricted maximum likelihood (REML) is now calculated using the full likelihood equation (without leaving out additive constants)

- the model deviance is now calculated as -2 times the difference between the model log-likelihood and the log-likelihood under the saturated model (this is a more appropriate definition of the deviance than just taking -2 times the model log-likelihood)

- naming scheme of illustrative datasets bundled with the package has been changed; now datasets are called `<dat.authoryear>`; therefore, the datasets are now called (`old name -> new name`):
  - `dat.bcg      -> dat.colditz1994`
  - `dat.warfarin -> dat.hart1999`
  - `dat.los      -> dat.normand1999`
  - `dat.co2      -> dat.curtis1998`
  - `dat.empint   -> dat.mcdaniel1994`

- but `dat.bcg` has been kept as an alias for `dat.colditz1994`, as it has been referenced under that name in some publications

- added new dataset (`dat.pritz1997`) to illustrate the meta-analysis of proportions (raw values and transformations thereof)

- added new dataset (`dat.bonett2010`) to illustrate the meta-analysis of Cronbach's alpha values (raw values and transformations thereof)

- added new datasets (`dat.hackshaw1998`, `dat.raudenbush1985`)

- (approximate) standard error of the tau^2 estimate is now computed and shown for most of the (residual) heterogeneity estimators

- added `nobs()` and `df.residual()` methods for objects of class `rma`

- `metafor.news()` is now simply a wrapper for `news(package="metafor")`

- the package code is now byte-compiled, which yields some modest increases in execution speed

- some general code cleanup

- the `metafor` package no longer depends on the `nlme` package

- some improvements to the documentation

# metafor 1.6-0 (2011-04-13)

- `trimfill.rma.uni()` now returns a proper object even when the number of missing studies is estimated to be zero

- added the (log transformed) ratio of means as a possible outcome measure to the `escalc()` and `rma.uni()` functions (`measure="ROM"`)

- added new dataset (`dat.co2`) to illustrate the use of the ratio of means outcome measure

- some additional error checking in the various forest functions (especially when using the `ilab` argument)

- in `labbe.rma()`, the solid and dashed lines are now drawn behind (and not on top of) the points

- slight change to `transf.ipft.hm()` so that missing values in `targs$ni` are ignored

- some improvements to the documentation

# metafor 1.5-0 (2010-12-16)

- the `metafor` package now has its own project website at: https://www.metafor-project.org

- added `labbe()` function to create L'Abbe plots

- the `forest.default()` and `addpoly.default()` functions now allow the user to directly specify the lower and upper confidence interval bounds (this can be useful when the CI bounds have been calculated with other methods/functions)

- added the incidence rate for a single group and for two groups (and transformations thereof) as possible outcome measures to the `escalc()` and `rma.uni()` functions (`measure="IRR"`, `"IRD"`, `"IRSD"`, `"IR"`, `"IRLN"`, `"IRS"`, and `"IRFT"`)

- added the incidence rate ratio as a possible outcome measure to the `rma.mh()` function

- added transformation functions related to incidence rates

- added the Freeman-Tukey double arcsine transformation and its inverse to the transformation functions

- added some additional error checking for out-of-range p-values in the `permutest.rma.uni()` function

- added some additional checking for out-of-range values in several transformation functions

- added `confint()` methods for `rma.mh` and `rma.peto` objects (only for completeness sake; print already provides CIs)

- added new datasets (`dat.warfarin`, `dat.los`, `dat.empint`)

- some improvements to the documentation

# metafor 1.4-0 (2010-07-30)

- a paper about the package has now been published in the Journal of Statistical Software (https://www.jstatsoft.org/v36/i03/)

- added citation info; see: `citation("metafor")`

- the `metafor` package now depends on the `nlme` package

- added extractor functions for the AIC, BIC, and deviance

- some updates to the documentation

# metafor 1.3-0 (2010-06-25)

- the `metafor` package now depends on the `Formula` package

- made `escalc()` generic and implemented a default and a formula interface

- added the (inverse) arcsine transformation to the set of transformation functions

# metafor 1.2-0 (2010-05-18)

- cases where k is very small (e.g., k equal to 1 or 2) are now handled more gracefully

- added sanity check for cases where all observed outcomes are equal to each other (this led to division by zero when using the Knapp & Hartung method)

- the "smarter way to set the number of iterations for permutation tests" (see notes for previous version below) now actually works like it is supposed to

- the `permutest.rma.uni()` function now provides more sensible results when k is very small; the documentation for the function has also been updated with some notes about the use of permutation tests under those circumstances

- made some general improvements to the various forest plot functions making them more flexible in particular when creating more complex displays; most importantly, added a `rows` argument and removed the `addrows` argument

- some additional examples have been added to the help files for the forest and addpoly functions to demonstrate how to create more complex displays with these functions

- added `showweight` argument to the `forest.default()` and `forest.rma()` functions

- `cumul()` functions not showing all of the output columns when using fixed-effects models has been corrected

- `weights.rma.uni()` function now handles `NA`s appropriately

- `weights.rma.mh()` and `weights.rma.peto()` functions added

- `logLik.rma()` function now behaves more like other `logLik()` functions (such as `logLik.lm()` and `logLik.lme()`)

# metafor 1.1-0 (2010-04-28)

- `cint()` generic removed and replaced with `confint()` method for objects of class `rma.uni`

- slightly improved the code to set the x-axis title in the `forest()` and `funnel()` functions

- added `coef()` method for `permutest.rma.uni` objects

- added `append` argument to `escalc()` function

- implemented a smarter way to set the number of iterations for permutation tests (i.e., the `permutest.rma.uni()` function will now switch to an exact test if the number of iterations required for an exact test is actually smaller than the requested number of iterations for an approximate test)

- changed the way how p-values for individual coefficients are calculated in `permutest.rma.uni()` to 'two times the one-tailed area under the permutation distribution' (more consistent with the way we typically define two-tailed p-values)

- added `retpermdist` argument to `permutest.rma.uni()` to return the permutation distributions of the test statistics

- slight improvements to the various transformation functions to cope better with some extreme cases

- p-values are now calculated in such a way that very small p-values stored in fitted model objects are no longer truncated to 0 (the printed results are still truncated depending on the number of digits specified)

- changed the default number of iterations for the ML, REML, and EB estimators from 50 to 100

# metafor 1.0-1 (2010-02-02)

- version jump in conjunction with the upcoming publication of a paper in the Journal of Statistical Software describing the `metafor` package

- instead of specifying a model matrix, the user can now specify a model formula for the `mods` argument in the `rma()` function (e.g., like in the `lm()` function)

- `permutest()` function now allows exact permutation tests (but this is only feasible when k is not too large)

- `forest()` function now uses the `level` argument properly to adjust the CI level of the summary estimate for models without moderators (i.e., for fixed- and random-effets models)

- `forest()` function can now also show the prediction interval as a dashed line for a random-effects model

- information about the measure used is now passed on to the `forest()` and `funnel()` functions, which try to set an appropriate x-axis title accordingly

- `funnel()` function now has more arguments (e.g., `atransf`, `at`) providing more control over the display of the x-axis

- `predict()` function now has its own `print()` method and has a new argument called `addx`, which adds the values of the moderator variables to the returned object (when `addx=TRUE`)

- functions now properly handle the `na.action` `"na.pass"` (treated essentially like `"na.exclude"`)

- added method for `weights()` to extract the weights used when fitting models with `rma.uni()`

- some small improvements to the documentation

# metafor 0.5-7 (2009-12-06)

- added `permutest()` function for permutation tests

- added `metafor.news()` function to display the `NEWS` file of the `metafor` package within R (based on same idea in the `animate` package by Yihui Xie)

- added some checks for values below machine precision

- a bit of code reorganization (nothing that affects how the functions work)

# metafor 0.5-6 (2009-10-19)

- small changes to the computation of the DFFITS and DFBETAS values in the `influence()` function, so that these statistics are more in line with their definitions in regular linear regression models

- added option to the plot function for objects returned by `influence()` to allow plotting the covariance ratios on a log scale (now the default)

- slight adjustments to various `print()` functions (to catch some errors when certain values were `NA`)

- added a control option to `rma()` to adjust the step length of the Fisher scoring algorithm by a constant factor (this may be useful when the algorithm does not converge)

# metafor 0.5-5 (2009-10-08)

- added the phi coefficient (`measure="PHI"`), Yule's Q (`"YUQ"`), and Yule's Y (`"YUY"`) as additional measures to the `escalc()` function for 2x2 table data

- forest plots now order the studies so that the first study is at the top of the plot and the last study at the bottom (the order can still be set with the `order` or `subset` argument)

- added `cumul()` function for cumulative meta-analyses (with a corresponding `forest()` method to plot the cumulative results)

- added `leave1out()` function for leave-one-out diagnostics

- added option to `qqnorm.rma.uni()` so that the user can choose whether to apply the Bonferroni correction to the bounds of the pseudo confidence envelope

- some internal changes to the class and methods names

- some small corrections to the documentation

# metafor 0.5-4 (2009-09-18)

- corrected the `trimfill()` function

- improvements to various print functions

- added a `regtest()` function for various regression tests of funnel plot asymmetry (e.g., Egger's regression test)

- made `ranktest()` generic and added a method for objects of class `rma` so that the test can be carried out after fitting

- added `anova()` function for full vs reduced model comparisons via fit statistics and likelihood ratio tests

- added the Orwin and Rosenberg approaches to `fsn()`

- added H^2 measure to the output for random-effects models

- in `escalc()`, `measure="COR"` is now used for the (usual) raw correlation coefficient and `measure="UCOR"` for the bias corrected correlation coefficients

- some small corrections to the documentation

# metafor 0.5-3 (2009-07-31)

- small changes to some of the examples

- added the log transformed proportion (`measure="PLN"`) as another measure to the `escalc()` function; changed `"PL"` to `"PLO"` for the logit (i.e., log odds) transformation for proportions

# metafor 0.5-2 (2009-07-06)

- added an option in `plot.infl.rma.uni()` to open a new device for plotting the DFBETAS values

- thanks to Jim Lemon, added a much better method for adjusting the size of the labels, annotations, and symbols in the `forest()` function when the number of studies is large

# metafor 0.5-1 (2009-06-14)

- made some small changes to the documentation (some typos corrected, some confusing points clarified)

# metafor 0.5-0 (2009-06-05)

- first version released on CRAN
