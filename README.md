# MAPLE

MAPLE (a novel Mendelian randomization method with self-Adaptive determination of samPle structure and multiple pLEiotropic effects), is an R package for efficient statistical inference of mendelian randomization analysis. MAPLE utilizes a set of correlated SNPs, self-adaptively accounts for the sample structure and the uncertainty that these correlated SNPs may exhibit multiple pleiotropic effects, as well as explicitly models both uncorrelated and correlated horizontal pleiotropy. The term ‘self-adaptive’ represents MAPLE is able to automatically infer the sample structure and the probability that a SNP has one specific pleiotropy effect from the data at hand. In particular, MAPLE first acquires the accurate estimate of the nuisance error parameter using the genome-wide summary statistics, then places the inference of the causal effect into a likelihood-framework and relies on a scalable sampling-based algorithm to obtain calibrated $p$-values.

# Installation

It is easy to install the development version of MAPLE package using the 'devtools' package.

```{r eval=FALSE}
#install.packages("devtools")
library(devtools)
install_github("yuanzhongshang/MAPLE")
```

# Usage

The main function in the package is MAPLE, you can find the instructions by `'?MAPLE'`.

```{r eval=FALSE}
library(MAPLE)

?MAPLE
```

# Quick Start

See [Tutorial](https://yuanzhongshang.github.io/MAPLE/) for detailed documentation and examples.

# Development

This R package is developed by Liye Zhang, Jiadong Ji, Lu Liu and Zhongshang Yuan.
