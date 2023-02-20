# causalKNN: causal KNN and TEP treatment effect estimation

## Overview

This package implements the causal KNN and treatment effect projection (TEP) heterogeneous treatment effect estimation algorithms introduced in Hitsch, Misra, and Zhang (2023).


## Installation

To install this package in R, run the following commands in R:

```R
install.packages("devtools")
devtools::install_github("walterwzhang/causalKNN-dev", build_vignettes = TRUE)
```

Alternatively, if you are installing the package from source, run the following commands in R:

```R
install.packages(path_to_file, repos = NULL, type="source")
```

where `path_to_file` is the local file path to the repository.

## Usage

For how to use this package, please see the [Intro Vignette](https://walterwzhang.github.io/causalKNN-dev/articles/causalKNN-Intro.html) or access it from R:

```R
library(causalKNN)
vignette("causalKNN-Intro")
```

To bootstrap aggregate the estimators, please follow the [Bootstrap Vignette](https://walterwzhang.github.io/causalKNN-dev/articles/causalKNN-Bootstrap.html):

```R
vignette("causalKNN-Bootstrap")
```

## Notes

- The changelog can be found in [NEWS.md](NEWS.md)
- The K nearest neighbor index matrices are computed using the `knn.index.dist` function from the `KernelKnn` package
- The Elastic-Net (using `glmnet`) is offered as the treatment effect projection's regression algorithm
