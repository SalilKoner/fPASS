.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Power and Sample size (PASS) calculation
                         for two sample projection-based testing of group difference
                         under a repeatedly measured longitudinal or sparsely functional
                         design : Version ",
                         utils::packageVersion("fPASS"))
}
