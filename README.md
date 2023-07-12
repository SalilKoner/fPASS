
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- # pkgdown <img src="man/figures/logo.png" align="right" alt="" width="120" /> -->
<!-- <!-- badges: start -->
<!-- [![CRAN Status](https://www.r-pkg.org/badges/version/pkgdown)](https://cran.r-project.org/package=pkgdown){.pkgdown-release} -->
<!-- [![R-CMD-check](https://github.com/r-lib/pkgdown/workflows/R-CMD-check/badge.svg)](https://github.com/r-lib/pkgdown/actions){.pkgdown-devel} -->
<!-- [![Codecov test coverage](https://codecov.io/gh/r-lib/pkgdown/branch/main/graph/badge.svg)](https://app.codecov.io/gh/r-lib/pkgdown?branch=main) -->
<!-- badges: end -->

# fPASS: An R package for Power and Sample Size analysis (PASS) for Projection-based Two-Sample test for functional data.

[Salil Koner](https://biostat.duke.edu/profile/salil-koner)

For further details about the testing procedure, please see [Wang
(2021)](https://doi.org/10.1214/21-EJS1802) for the univariate case, and
[Koner and Luo (2023)](https://arxiv.org/abs/2302.05612) for the
multivariate case.

fPASS is designed to make it quick and easy software for randomized
clinical trial simulation tool for determining treatment efficacy where
the response collected under a longitudinal or functional design. The
current development version of the package can be installed by running
the following.

## Installation

<!-- ::: .pkgdown-release -->
<!-- ```{r, eval = FALSE} -->
<!-- # Install released version from CRAN -->
<!-- install.packages("pkgdown") -->
<!-- ``` -->
<!-- ::: -->

<div class=".pkgdown-devel">

``` r
# Install development version from GitHub
devtools::install_github("SalilKoner/fPASS", build_vignettes = TRUE) # Vignettes takes about 20 minutes to run. 
```

</div>
