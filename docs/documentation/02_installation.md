---
layout: page
title: Installation
description: ~
---

## Dependencies

* R version >= 3.6.0
* R packages: Rcpp, RcppArmadillo, RcppDist, dplyr, magrittr, readr

### 1. Install devtools if necessary

```{r eval=FALSE}
install.packages('devtools')
```

### 2. Install MAPLE

```{r eval=FALSE}
devtools::install_github('yuanzhongshang/MAPLE')
```

### 3. Load package

```{r eval=FALSE}
library(MAPLE)
```

### 4. Some possible issues when installing the package, especially on the MacOS system

(1)Cannot find tools necessary when using R in the MacOS system.

```{r eval=FALSE}
Error: Failed to install 'MAPLE' from GitHub:
  Could not find tools necessary to compile a package
Call `pkgbuild::check_build_tools(debug = TRUE)` to diagnose the problem.
```

possible solution: in R, type the code  options(buildtools.check = function(action) TRUE ), see the discussion about this error, [link](https://stackoverflow.com/questions/37776377/error-when-installing-an-r-package-from-github-could-not-find-build-tools-neces)

(2)library not found for -lgfortran when using R in the MacOS system.

```{r eval=FALSE}
ld: library not found for -lgfortran
```

It seems the gfortran is not well installed on the MacOS system. Please check this [link](https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/) for the gfortran installation to see if it helps.