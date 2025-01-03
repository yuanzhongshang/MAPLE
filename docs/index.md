---
layout: full
homepage: true
disable_anchors: true
description: Mendelian randomization method with self-Adaptive determination of samPle structure and multiple pLEiotropic effects
---
# MAPLE Overview

![](MAPLE.png)

MAPLE is an R package for efficient statistical inference of mendelian randomization analysis. MAPLE utilizes a set of correlated SNPs, self-adaptively accounts for the sample structure and the uncertainty that these correlated SNPs may exhibit multiple pleiotropic effects, as well as explicitly models both uncorrelated and correlated horizontal pleiotropy. MAPLE is implemented as an open-source R package, freely available at <https://github.com/yuanzhongshang/MAPLE>.

# Installation
You can install the released version of MAPLE from Github with the following code. This package is supported for Windows 10/11, and Linux. The package has been tested on the following systems:
* Windows 10, 11
* Linux: Ubuntu (22.04.4)

Dependencies

* R version >= 3.6.0
* R packages: Rcpp, RcppArmadillo, RcppDist, dplyr, magrittr, readr

## 1. Install devtools if necessary

```
install.packages('devtools')
```

## 2. Install MAPLE

```
devtools::install_github('yuanzhongshang/MAPLE')
```

We tested the install time on Windows 11:
```
system.time({devtools::install_github('yuanzhongshang/MAPLE')})
user  system elapsed
  2.66   0.78   90.55
```

## 3. Load package

```
library(MAPLE)
```

## Issues

All feedback, bug reports and suggestions are warmly welcomed! Please make sure to raise issues with a detailed and reproducible example and also please provide the output of your sessionInfo() in R!

# How to use MAPLE
The MAPLE User Manual: [here](MAPLE_user_manual.pdf)
