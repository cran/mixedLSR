
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mixedLSR <img src="man/figures/logo.png" align="right" width="25%"/>

<!-- badges: start -->

[![R-CMD-check](https://github.com/alexanderjwhite/mixedLSR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/alexanderjwhite/mixedLSR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Mixed, low-rank, and sparse multivariate regression (mixedLSR) provides
tools for performing mixture regression when the coefficient matrix is
low-rank and sparse. mixedLSR allows subgroup identification by
alternating optimization with simulated annealing to encourage global
optimum convergence. This method is data-adaptive, automatically
performing parameter selection to identify low-rank substructures in the
coefficient matrix.

## Installation

You can install the development version of mixedLSR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("alexanderjwhite/mixedLSR")
```
