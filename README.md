metafor: A Meta-Analysis Package for R
======================================

[![Build Status](https://travis-ci.org/wviechtb/metafor.svg?branch=master)](https://travis-ci.org/wviechtb/metafor)
[![Code Coverage](https://img.shields.io/codecov/c/github/wviechtb/metafor.svg)](https://codecov.io/github/wviechtb/metafor?branch=master)
[![CRAN Version](http://www.r-pkg.org/badges/version/metafor)](http://cran.rstudio.com/web/packages/metafor)
[![Monthly Downloads](http://cranlogs.r-pkg.org/badges/metafor)](http://cranlogs.r-pkg.org/badges/metafor)
[![Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/metafor)](http://cranlogs.r-pkg.org/badges/grand-total/metafor)

## Description

The `metafor` package is a comprehensive collection of functions for conducting meta-analyses in R. The package includes functions to calculate various effect sizes or outcome measures, fit fixed-, random-, and mixed-effects models to such data, carry out moderator and meta-regression analyses, and create various types of meta-analytical plots (e.g., forest, funnel, radial, L'Abbé, Baujat, GOSH plots). For meta-analyses of binomial and person-time data, the package also provides functions that implement specialized methods, including the Mantel-Haenszel method, Peto's method, and a variety of suitable generalized linear (mixed-effects) models (i.e., mixed-effects logistic and Poisson regression models). Finally, the package provides functionality for fitting meta-analytic multivariate/multilevel models that account for non-independent sampling errors and/or true effects (e.g., due to the inclusion of multiple treatment studies, multiple endpoints, or other forms of clustering). Network meta-analyses and meta-analyses accounting for known correlation structures (e.g., due to phylogenetic relatedness) can also be conducted.

## Package Website

The `metafor` package website can be found at [http://www.metafor-project.org](http://www.metafor-project.org). On the website, you can find:

* some [news](http://www.metafor-project.org/doku.php/news) concerning the package and/or its development,
* a more detailed description of the [package features](http://www.metafor-project.org/doku.php/features),
* a log of the [package updates](http://www.metafor-project.org/doku.php/updates) that have been made over the years,
* a [to-do list](http://www.metafor-project.org/doku.php/todo) and a description of planned features to be implemented in the future,
* information on how to [download and install](http://www.metafor-project.org/doku.php/installation) the package,
* information on how to obtain [documentation and help](http://www.metafor-project.org/doku.php/help) with using the package,
* some [analysis examples](http://www.metafor-project.org/doku.php/analyses) that illustrate various models, methods, and techniques,
* a little showcase of [plots and figures](http://www.metafor-project.org/doku.php/plots) that can be created with the package,
* some [tips and notes](http://www.metafor-project.org/doku.php/tips) that may be useful when working with the package,
* a list of people that have in some shape or form [contributed](http://www.metafor-project.org/doku.php/contributors) to the development of the package,
* an (incomplete) [list of articles](http://www.metafor-project.org/doku.php/articles) that have used the package as part of the analyses,
* a [frequently asked questions](http://www.metafor-project.org/doku.php/faq) section, and
* some [links](http://www.metafor-project.org/doku.php/links) to other websites related to software for meta-analysis.

## Installation

The current official (i.e., [CRAN](https://cran.r-project.org/package=metafor)) release can be installed directly within R with:
```r
install.packages("metafor")
```

After installing the [devtools](https://cran.r-project.org/package=devtools) package with ```install.packages("devtools")```, the development version of the metafor package can be installed with:
```r
library("devtools")
install_github("wviechtb/metafor")
```
This approach builds the package from source based on the development branch on GitHub. On Windows, [Rtools](http://cran.r-project.org/bin/windows/Rtools/) must be installed for this to work.

## Meta

The metafor package was written by [Wolfgang Viechtbauer](http://www.wvbauer.com/). It is licensed under the [GNU General Public License Version 2](http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt). For citation info, type `citation(package='metafor')` in R. To report any issues or bugs, please go [here](https://github.com/wviechtb/metafor/issues).
