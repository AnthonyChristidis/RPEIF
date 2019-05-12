
[![Build Status](https://travis-ci.org/AnthonyChristidis/IFs.svg?branch=master)](https://travis-ci.com/AnthonyChristidis/IFs) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/IFs)](https://cran.r-project.org/package=IFs) [![Downloads](http://cranlogs.r-pkg.org/badges/IFs)](https://cran.r-project.org/package=IFs)

IFs
===

This package provides functions for computing the influence functions of risk and performance measures.

------------------------------------------------------------------------

### Installation

You can install the **stable** version on [R CRAN](https://cran.r-project.org/package=IFs).

``` r
install.packages("IFs", dependencies = TRUE)
```

You can install the **development** version from [GitHub](https://github.com/AnthonyChristidis/IFs)

``` r
library(devtools)
devtools::install_github("AnthonyChristidis/IFs")
```

### Usage

``` r
# A small example
library(IFs)
# Computing the IF of the returns (with outlier cleaning and prewhitening)
# Loading the data
data(edhec, package="PerformanceAnalytics")
colnames(edhec) = c("CA", "CTAG", "DIS", "EM","EMN", "ED", "FIA",
                    "GM", "LS", "MA", "RV", "SS", "FoF")
outIF <- IF(risk="mean",
            returns=edhec[,"CA"], evalShape=FALSE, retVals=NULL, nuisance.par=NULL,
            IFplot=TRUE, IFprint=TRUE,
            compile=TRUE, prewhiten=FALSE,
            cleanOutliers=TRUE, cleanMethod=c("locScaleRob", "Boudt")[1], eff=0.99, alpha.robust=0.05)
```

### License

This package is free and open source software, licensed under GPL (&gt;= 2).
