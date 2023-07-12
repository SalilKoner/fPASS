#' \code{fPASS} package
#'
#' See the README on \href{https://github.com/SalilKoner/fPASS}{GitHub}
#' @keywords internal
"_PACKAGE"


#' @details R package \pkg{fPASS} is designed to perform the power and sample size analysis
#' for functional under a dense and sparse (random) design and longitudinal data. The
#' function can handle data from wide variety of covariance structure, can be parametric,
#' or non-parametric. The user have a lot of flexibility into tweaking the arguments of the function
#' to assess the power function of the test under different sampling design and covariance process
#' of the response trajectory, and for any arbitrary mean difference function. Overall, the
#' functionality of the module is quite comprehensive and includes all the different cases
#' considered in the NCSS PASS software. We believe that this software can be an effective
#' clinical trial design tools when considering the projection-based test as the primary
#' decision making method.
#' @seealso See the primary function [fPASS::PASS_Proj_Test_ufDA()] to compute
#' the power and sample size of the test.
## usethis namespace: start
#' @docType package
#' @importFrom utils globalVariables
#' @name fPASS-package
## usethis namespace: end
NULL

if(getRversion() >= "2.15.1")  utils::globalVariables(c(".id.", "argvals", "subj", "y"))
