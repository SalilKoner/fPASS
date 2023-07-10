#' Power of the Two-sample Projection-based test for functional data with known (or estimated)
#' eigencomponents.
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' @description
#' The function `Power_Proj_Test_ufDA()` computes the power of
#' of the two-sample projection-based test for functional response data
#' setting, when the group difference, the eigenfunctions of the covariance
#' of the data are specified at dense grid of time points, along with the
#' (estimated) covariance of the `shrinkage` scores.
#'
#' @details The projection-based test first extracts K eigenfunctions from the data, and then
#' project the mean difference function onto each of the eigenfunctions to obtain a K-dimensional
#' projection vector that reflects the group difference. Wang (2021) pointed that under the null
#' hypothesis the covariance of K-dimensional functional principal component analysis (fPCA) scores
#' are the same, and thus a Hotelling \eqn{T^2} test with assuming equal variance of the shrinkage scores
#' is a valid test. However, Koner and Luo (2023) pointed out that under the alternate hypothesis,
#' when the difference is mean is significant, the covariance of the shrinkage scores also differ
#' between the groups. Therefore, while computing the power of test, we must have to derive the
#' distribution of the Hotelling \eqn{T^2} statistic under the assumption of unequal variance. The
#' alogrithm for the power of multivariate Hotelling \eqn{T^2} under unequal variance
#' is coded in [pHotellingT()] function. This particular function is a wrapper around that
#' function, which inputs the mean difference as a function, and the eigenfunctions and
#' the scores, and subsequently call the [fPASS::pHotellingT()] function to compute the power
#' under unequal variance. See Koner and Luo (2023) for more details on the
#' formula of the non-null distribution.
#'
#' @author Salil Koner \cr Maintainer: Salil Koner
#' \email{salil.koner@@duke.edu}
#' @import Matrix
#' @import rlang
#' @import lifecycle
#' @import dplyr
#' @importFrom testthat expect_true expect_gte expect_no_error
#' @importFrom stats cov rnorm
#' @param total_sample_size Total sample size combing the two groups, must be a positive integer.
#' @param argvals The working grid of timepoints to evaluate the eigenfunctions and the mean functions.
#'                It is preferred to take the working grid as dense grid so that
#'                \eqn{\int [\mu_1(t) - \mu_2(t)]\phi_k(t) \,dt} can be calculated with a required precision.
#' @param mean_vector The difference in the mean function evaluated at argvals, must be a numeric vector of length same
#'                    as that that of argavls.
#' @param eigen_matrix The matrix of eigenfunctions evaluated at argvals,
#'                     must be a length(argvals) by K matrix, where K is the number of eigenfunctions.
#' @param scores_var1 The true (or estimate) of covariance matrix of the shrinkage scores for the first group.
#' Must be symmetric (\code{is.symmetric(scores_var1) == TRUE}) and positive definite
#' (\code{chol(scores_var1)} without an error!).
#' @param scores_var2 The true (or estimate) of covariance matrix of the shrinkage scores for the second group.
#' Must be symmetric (\code{is.symmetric(scores_var2) == TRUE}) and positive definite
#' (\code{chol(scores_var2)} without an error!).
#' @param weights The weights to put to compute the projection \eqn{\int [\mu_1(t) - \mu_2(t)]\phi_k(t) \,dt},
#' for each \eqn{k=1,\dots, K}. The integral is numerically approximated as
#' \code{sum(mean_diff(argvals)*eigen_matrix[,k]*weights)}.
#' @param sig.level Significance level of the test, default set at 0.05, must be less than 0.2.
#'                  This is used to compute the critical value of the test.
#' @param alloc.ratio Allocation of total sample size into the two groups. Must set as a vector of two
#'                    positive numbers. For equal allocation it should be put as c(1,1), for non-equal
#'                    allocation one can put c(2,1) or c(3,1) etc. Default set as c(1,1).
#' @param npc_to_pick Number of eigenfunction to be used to compute the power. Typically this is
#'                    becomes handy when the user want to discard few of the last eigenfunctions,
#'                    typically with a very small eigenvalues.
#' @return Power of the projection-based test for specified difference in the mean function
#' and the eigencomponents of the covariance of the functional data.
#' @seealso See [fPASS::pHotellingT()] and [fPASS::Sim_HotellingT_unequal_var()] for samples
#' from Hotelling T distribution.
#' @references Wang, Qiyao (2021)
#' \emph{Two-sample inference for sparse functional data,  Electronic Journal of Statistics,
#' Vol. 15, 1395-1423} \cr
#' \doi{https://doi.org/10.1214/21-EJS1802}.
#' @export Power_Proj_Test_ufDA
#' @examples
#' \dontrun{
#' ngrid          <- 101
#' interval       <- c(-1,1)
#' gauss.quad.pts <- gss::gauss.quad(ngrid,interval) # evaluation points
#' working.grid   <- gauss.quad.pts$pt
#' mean_fn        <- function(t) {0.4*sin(2*pi*t)}
#' mean_vector    <- mean_fn(working.grid)
#' eigen_fn       <- function(t, k){ sqrt(2)*{(k==2)*sin(2*pi*t) + (k==1)*cos(2*pi*t)} }
#' eigen_matrix   <- cbind(eigen_fn(working.grid,1), eigen_fn(working.grid,2))
#' mean_proj      <- sapply(1:2, function(r) integrate(function(x)
#' eigen_fn(x,r)*mean_fn(x), interval[1], interval[2])$value)
#' sig1           <- diag(2)
#' sig2           <- 2*diag(2)
#' alp            <- 0.05
#' n              <- 100
#' k              <- ncol(eigen_matrix)
#' cutoff         <- {(n - 2)*k/(n - k -1)}*qf(1-alp, k, n-k-1)
#' func_power     <- Power_Proj_Test_ufDA(total_sample_size=n,
#' argvals=working.grid,
#' mean_vector = mean_vector, eigen_matrix = eigen_matrix,
#' scores_var1 = sig1, scores_var2= sig2, weights = gauss.quad.pts$wt,
#' sig.level=alp, alloc.ratio = c(1,1), npc_to_pick=ncol(eigen_matrix))
#' }
Power_Proj_Test_ufDA <- function(total_sample_size, argvals,
                                 mean_vector, eigen_matrix,
                                 scores_var1, scores_var2, weights,
                                 sig.level=0.05, alloc.ratio = c(1,1),
                                 npc_to_pick=ncol(eigen_matrix)){

  #*********************************************************************************************%
  #*                                  Missing argument checking                                 %
  #*********************************************************************************************%
  call.f       <- rlang::call_match(defaults = FALSE)
  call.f.all   <- rlang::call_match(defaults = TRUE)
  default_args <- setdiff(rlang::call_args_names(call.f.all), rlang::call_args_names(call.f))
  if (length(default_args) > 0 ){
    rlang::inform(paste0("User did not specify values to the following arguments: ",
                         paste0(default_args, collapse = ","), " ---- default values are used"))
  }
  # Checking null values of argument
  non_null_arg_vals <- !sapply(rlang::call_args(call.f), is.null)
  testthat::expect_true(all(non_null_arg_vals),
                        info = paste0("NULL values are specified for the arguments with no default values: ",
                                      paste0(rlang::call_args_names(call.f)[!non_null_arg_vals], collapse = ",")) )

  #*********************************************************************************************%
  #*                         Compatibility checking of arguments                                %
  #*********************************************************************************************%
  # argument checking: total_sample_size
  testthat::expect_true(rlang::is_double(total_sample_size, n=1, finite = TRUE) & (total_sample_size > 0),
                        info = "total_sample_size must be a positive integer")
  # argument checking: sig.level
  testthat::expect_true(rlang::is_double(sig.level, n=1, finite=TRUE) &
                          (sig.level > 0) & (sig.level <= 0.2),
                        info = "level of significance sig.level must be number between 0 and 0.2");
  # argument checking: alloc.ratio
  testthat::expect_true(rlang::is_double(alloc.ratio, n=2, finite=TRUE) & all(alloc.ratio > 0) ,
                        info = "alloc.ratio must be positive numeric vector of length 2");
  # argument checking : argvals
  testthat::expect_true(rlang::is_double(argvals, finite=TRUE) & !is.unsorted(argvals),
                        info = "argvals must be a strictly increasing numeric vector")
  # argument checking: mean_vector
  testthat::expect_true(rlang::is_double(mean_vector, n=length(argvals), finite=TRUE),
                        info = "mean_vector must a numeric vector evaluated at argvals")
  # argument checking: eigen_matrix
  testthat::expect_true(rlang::is_double(eigen_matrix, finite=TRUE) & is.matrix(eigen_matrix) &
                        (nrow(eigen_matrix) == length(argvals)),
                        info ="eigen_matrix has to be a length(argvals) by K numeric matrix
                               where the K eigenfunctions are evaluated at the argvals")
  # argument checking: scores_var1
  testthat::expect_true(rlang::is_double(scores_var1, n=(ncol(eigen_matrix))^2, finite=TRUE) &
                        isSymmetric(scores_var1),
                        info = "scores_var1 must be a symmetric matrix with dimension same as
                                the number of columns of eigen_matrix.");
  # argument checking: scores_var2
  testthat::expect_true(rlang::is_double(scores_var2, n=length(scores_var1), finite=TRUE) &
                          isSymmetric(scores_var2),
                        info = "scores_var2 must be a symmetric matrix of
                                dimension same as that of scores_var1");
  # argument checking: weights
  testthat::expect_true(rlang::is_double(weights, n=length(argvals), finite=TRUE),
                        info = "weights must be numeric vector")
  # argument checking: npc_to_pick
  testthat::expect_true(rlang::is_integerish(npc_to_pick,n=1, finite=TRUE) &
                          (npc_to_pick > 0) & (npc_to_pick <= ncol(eigen_matrix)),
                        info = "npc_to_pick must be a positive number less than the
                                number of columns of eigen_matrix")
  # argument checking : scores_var1, scores_var2
  testthat::expect_no_error(chol(scores_var1), message = "scores_var1 is not positive definite")
  testthat::expect_no_error(chol(scores_var2), message = "scores_var2 is not positive definite")

  #*********************************************************************************************%
  #*                                          Main body                                         %
  #*********************************************************************************************%

  projection     <- colSums(sweep(eigen_matrix[, 1:npc_to_pick, drop=FALSE], 1,
                                  mean_vector*weights, FUN = "*"))
  sig1           <- scores_var1[1:npc_to_pick,  1:npc_to_pick, drop=FALSE]
  sig2           <- scores_var2[1:npc_to_pick,  1:npc_to_pick, drop=FALSE]
  critical.value <- {{(total_sample_size - 2)*npc_to_pick}/(total_sample_size - npc_to_pick -1)}*
                       qf(1-sig.level, npc_to_pick, total_sample_size-npc_to_pick-1)

  pHotellingT(q=critical.value, total_sample_size=total_sample_size, mean_diff=projection,
              sig1=sig1, sig2=sig2, alloc.ratio=alloc.ratio, lower.tail=FALSE)
}

