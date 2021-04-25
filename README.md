metafor: A Meta-Analysis Package for R
======================================

[![License: GPL (>=2)](https://img.shields.io/badge/license-GPL-blue)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
[![R build status](https://github.com/wviechtb/metafor/workflows/R-CMD-check/badge.svg)](https://github.com/wviechtb/metafor/actions)
[![Code Coverage](https://codecov.io/gh/wviechtb/metafor/branch/master/graph/badge.svg)](https://codecov.io/gh/wviechtb/metafor)
[![CRAN Version](https://www.r-pkg.org/badges/version/metafor)](https://cran.r-project.org/package=metafor)
[![devel Version](https://img.shields.io/badge/devel-2.5--88-brightgreen.svg)](https://www.metafor-project.org/doku.php/installation#development_version)
[![Monthly Downloads](https://cranlogs.r-pkg.org/badges/metafor)](https://cranlogs.r-pkg.org/badges/metafor)
[![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/metafor)](https://cranlogs.r-pkg.org/badges/grand-total/metafor)

## Description

The `metafor` package is a comprehensive collection of functions for conducting meta-analyses in R. The package includes functions to calculate various effect sizes or outcome measures, fit fixed-, random-, and mixed-effects models to such data, carry out moderator and meta-regression analyses, and create various types of meta-analytical plots (e.g., forest, funnel, radial, L'Abbé, Baujat, GOSH plots). For meta-analyses of binomial and person-time data, the package also provides functions that implement specialized methods, including the Mantel-Haenszel method, Peto's method, and a variety of suitable generalized linear (mixed-effects) models (i.e., mixed-effects logistic and Poisson regression models). Finally, the package provides functionality for fitting meta-analytic multivariate/multilevel models that account for non-independent sampling errors and/or true effects (e.g., due to the inclusion of multiple treatment studies, multiple endpoints, or other forms of clustering). Network meta-analyses and meta-analyses accounting for known correlation structures (e.g., due to phylogenetic relatedness) can also be conducted.

## Package Website

The `metafor` package website can be found at [https://www.metafor-project.org](https://www.metafor-project.org). On the website, you can find:

* some [news](https://www.metafor-project.org/doku.php/news:news) concerning the package and/or its development,
* a more detailed description of the [package features](https://www.metafor-project.org/doku.php/features),
* a log of the [package updates](https://www.metafor-project.org/doku.php/updates) that have been made over the years,
* a [to-do list](https://www.metafor-project.org/doku.php/todo) and a description of planned features to be implemented in the future,
* information on how to [download and install](https://www.metafor-project.org/doku.php/installation) the package,
* information on how to obtain [documentation and help](https://www.metafor-project.org/doku.php/help) with using the package,
* some [analysis examples](https://www.metafor-project.org/doku.php/analyses) that illustrate various models, methods, and techniques,
* a little showcase of [plots and figures](https://www.metafor-project.org/doku.php/plots) that can be created with the package,
* some [tips and notes](https://www.metafor-project.org/doku.php/tips) that may be useful when working with the package,
* a list of people that have in some shape or form [contributed](https://www.metafor-project.org/doku.php/contributors) to the development of the package,
* a [frequently asked questions](https://www.metafor-project.org/doku.php/faq) section, and
* some [links](https://www.metafor-project.org/doku.php/links) to other websites related to software for meta-analysis.

## Documentation

A good starting place for those interested in using the `metafor` package is the following paper:

Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. *Journal of Statistical Software, 36*(3), 1-48. [https://www.jstatsoft.org/v36/i03/](https://www.jstatsoft.org/v36/i03/).

In addition to reading the paper, carefully read the [package intro](https://wviechtb.github.io/metafor/reference/metafor-package.html) and then the help pages for the [`escalc`](https://wviechtb.github.io/metafor/reference/escalc.html) and the [`rma.uni`](https://wviechtb.github.io/metafor/reference/rma.uni.html) functions (or the [`rma.mh`](https://wviechtb.github.io/metafor/reference/rma.mh.html), [`rma.peto`](https://wviechtb.github.io/metafor/reference/rma.peto.html), [`rma.glmm`](https://wviechtb.github.io/metafor/reference/rma.glmm.html), [`rma.mv`](https://wviechtb.github.io/metafor/reference/rma.mv.html) functions if you intend to use these methods). The help pages for these functions provide links to many additional functions, which can be used after fitting a model. You can also read the entire documentation online at [https://wviechtb.github.io/metafor/](https://wviechtb.github.io/metafor/) (where it is nicely formatted, equations are shown correctly, and the output from all examples is provided).

## Installation

The current official (i.e., [CRAN](https://cran.r-project.org/package=metafor)) release can be installed within R with:
```r
install.packages("metafor")
```

The development version of the package can be installed with:
```r
install.packages("remotes")
remotes::install_github("wviechtb/metafor")
```
This builds the package from source based on the current version on [GitHub](https://github.com/wviechtb/metafor).

## Meta

The metafor package was written by [Wolfgang Viechtbauer](https://www.wvbauer.com/). It is licensed under the [GNU General Public License](https://www.gnu.org/licenses/old-licenses/gpl-2.0.txt). For citation info, type `citation(package='metafor')` in R. To report any issues or bugs or to suggest enhancements to the package, please go [here](https://github.com/wviechtb/metafor/issues).
