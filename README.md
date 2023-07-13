
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- # pkgdown <img src="man/figures/logo.png" align="right" alt="" width="120" /> -->
<!-- badges: start -->

[![R-CMD-check](https://github.com/SalilKoner/fPASS/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/SalilKoner/fPASS/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# fPASS: Power and Sample size for projection-based testing under a repeatedly measured design

[Salil Koner](https://biostat.duke.edu/profile/salil-koner)

The details of the power and sample size formula and the relevant
computational details are documented in the
[manuscript](https://salilkoner.github.io/assets/PASS_manuscript.pdf).
The users are encourage to see [Wang
(2021)](https://doi.org/10.1214/21-EJS1802) and [Koner and Luo
(2023)](https://arxiv.org/abs/2302.05612) for further details about the
testing procedure.

**fPASS** is designed to make it quick and easy software for randomized
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
devtools::install_github("SalilKoner/fPASS") # Vignettes takes about 20 minutes to run. 
```

</div>

## Vignettes

If you want to install the package with the vignettes to be built, then
run

``` r
# Install development version from GitHub
devtools::install_github("SalilKoner/fPASS", build_vignettes = TRUE) # Vignettes takes about 20 minutes to run. 
```

followed by `browseVignettes("fPASS")` to see the application of the
package in real life case studies.
