.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Power and Sample size (PASS) Calculation
                         for Projection-based Testing of
                         Sparsely Observed (Univariate and Multivariate)
                         Functional Data: Version ",
                         utils::packageVersion("fPASS"))
}
