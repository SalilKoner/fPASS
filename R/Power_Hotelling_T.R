#' The approximate degrees of freedom formula for sum of Wishart.
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' @description
#' The approximate degrees of freedom formula for sum of two
#' independent Wishart random variable
#' with parameter sig1 and sig2, and degrees of freedom n1-1 and n2-1
#' where n1 + n2 is equal to the `total_sample_size`.
#'
#' See Koner and Luo (2023) for more details on the formula for degrees of freedom.
#' @import Matrix
#' @import lifecycle
#' @import rlang
#' @importFrom testthat expect_true expect_gte expect_no_error
#' @importFrom expm sqrtm
#' @param total_sample_size Target sample size, must be a positive integer.
#' @param sig1 The true (or estimate) of covariance matrix for the first group. Must be symmetric
#' (\code{is.symmetric(sig1) == TRUE}) and positive definite (\code{chol(sig1)} without an error!).
#' @param sig2 The true (or estimate) of covariance matrix for the second group. Must be symmetric
#' (\code{is.symmetric(sig2) == TRUE}) and positive definite (\code{chol(sig2)} without an error!).
#' @param alloc.ratio Allocation of total sample size into the two groups. Must set as a vector of two
#'                    positive numbers. For equal allocation it should be put as c(1,1), for non-equal
#'                    allocation one can put c(2,1) or c(3,1) etc.
#' @return The approximate degrees of freedom.
#' @export Sum_of_Wishart_df
#' @author Salil Koner \cr Maintainer: Salil Koner
#' \email{salil.koner@@duke.edu}
#' @seealso [fPASS::Sim_HotellingT_unequal_var()] and [fPASS::pHotellingT()].
#' @examples
#' \dontrun{
#' k <- 8
#' mu1  <- rep(0,k); del  <- 0.4; mu2 <- mu1 + rep(del, k);
#' sig1 <- diag(k); sig2 <- sig1 + del*toeplitz(c(1,rep(0.5, k-1)))
#' alt.dist.samples <- Sum_of_Wishart_df(total_sample_size=150, k=k,
#' sig1=sig1, sig2=sig2, alloc.ratio=c(2,1))
#' }
Sum_of_Wishart_df <- function(total_sample_size, alloc.ratio, sig1, sig2){


  #*********************************************************************************************%
  #*                                  Missing argument checking                                 %
  #*********************************************************************************************%
  # set.seed(12345)
  call.f       <- rlang::call_match(defaults = FALSE)
  call.f.all   <- rlang::call_match(defaults = TRUE)
  default_args <- setdiff(rlang::call_args_names(call.f.all), rlang::call_args_names(call.f))
  if (length(default_args) > 0 ){
    rlang::inform(paste0("User did not specify values to the following arguments: ",
                         paste0(default_args, collapse = ","), " ---- default values are used"))
  }
  non_null_arg_vals <- !sapply(rlang::call_args(call.f), is.null)
  testthat::expect_true(all(non_null_arg_vals),
                        info = paste0("NULL values are specified for the arguments with no default values: ",
                                      paste0(rlang::call_args_names(call.f)[!non_null_arg_vals], collapse = ",")) )

  #*********************************************************************************************%
  #*                         Compatibility checking of arguments                                %
  #*********************************************************************************************%
  # argument checking : total_sample_size
  testthat::expect_true(rlang::is_double(total_sample_size, n=1, finite=TRUE) & (total_sample_size > 0),
                        info = "total_sample_size must be positive")
  # argument checking : sig1
  testthat::expect_true(rlang::is_double(sig1, finite=TRUE) & isSymmetric(sig1),
                        info = "sig1 must be a symmetric matrix");
  # argument checking : sig2
  testthat::expect_true(rlang::is_double(sig2, n=length(sig1), finite=TRUE) & isSymmetric(sig2),
                        info = "sig2 must be a symmetric matrix of
                                dimension same as that of sig1");
  # argument checking : alloc.ratio
  testthat::expect_true(rlang::is_double(alloc.ratio, n=2, finite=TRUE) & all(alloc.ratio > 0),
                        info = "alloc.ratio must be positive numeric vector of length 2");

  testthat::expect_no_error(chol(sig1), message = "sig1 is not positive definite")
  testthat::expect_no_error(chol(sig2), message = "sig2 is not positive definite")

  k               <- ncol(sig1)
  n1              <- total_sample_size * {alloc.ratio[1]/sum(alloc.ratio)}
  n1_by_n2        <- alloc.ratio[1]/alloc.ratio[2]
  n2              <- total_sample_size - n1
  sig.til         <- sig1 + n1_by_n2 * sig2
  sqinv.sig.til   <- expm::sqrtm(solve(sig.til))
  delta           <- (sqinv.sig.til %*% sig1) %*% sqinv.sig.til
  I_minus_delta   <- diag(k)-delta

  delta.str       <- {n1_by_n2*(n1_by_n2 - 1/n2)}*delta + (1-1/n2)*I_minus_delta
  denom           <- {(n1_by_n2 - 1/n2)*(n1_by_n2^2)} *
                     {sum(diag(delta %*% delta)) + {sum(diag(delta))}^2} + (1-1/n2) *
                     {sum(diag(I_minus_delta %*% I_minus_delta)) + {sum(diag(I_minus_delta))}^2}
  numer           <- sum(diag(delta.str %*% delta.str)) + {sum(diag(delta.str))}^2
  denom.df        <- n2 * numer/denom

  denom.df
}


#' Samples from the non-null distribution of the Hotelling-\eqn{T^2} statistic under unequal covariance.
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' @description
#' The function \code{Sim_HotellingT_unequal_var()} generates samples from the
#' (non-null) distribution of the two-sample Hotelling-\eqn{T^2} statistic
#' under the assuming of unequal covariance of the multivariate response
#' between the two groups. This function is used to compute the power function
#' of Two-Sample (TS) Projection-based test (Wang 2021, EJS.)
#' for sparsely observed univariate functional data.
#'
#' @details
#' Under the assumption of the equal variance, we know that the alternative
#' distribution of the Hotelling-\eqn{T^2} statistic has an F distribution with the
#' non-centrality depending on the difference between the true mean vectors and the
#' (common) covariance of the response. However, when the true covariance of the true groups
#'  of responses differ, the alternate distribution becomes non-trivial. Koner and Luo (2023)
#'  proved that the alternate distribution of the test-statistic approximately follows
#'  a ratio of the linear combination of the K (dimension of the response) non-central
#'  chi-squared random variables (where the non-centrality parameter depends on the mean difference)
#'  and a chi-squared distribution whose degrees of freedom depends on a complicated functions of
#'  sample size in the two groups.
#'  See Koner and Luo (2023) for more details on the formula of the non-null distribution.
#'
#' @author Salil Koner \cr Maintainer: Salil Koner
#' \email{salil.koner@@duke.edu}
#'
#' @import Matrix
#' @import lifecycle
#' @import rlang
#' @importFrom testthat expect_true expect_gte expect_no_error
#' @importFrom stats cov rnorm complete.cases predict rchisq weighted.mean pf qf
#' @importFrom expm sqrtm
#' @param total_sample_size Target sample size, must be a positive integer.
#' @param mean_diff The difference in the mean vector between the two groups, must be a vector.
#' @param sig1 The true (or estimate) of covariance matrix for the first group. Must be symmetric
#' (\code{is.symmetric(sig1) == TRUE}) and positive definite (\code{chol(sig1)} without an error!).
#' @param sig2 The true (or estimate) of covariance matrix for the second group. Must be symmetric
#' (\code{is.symmetric(sig2) == TRUE}) and positive definite (\code{chol(sig2)} without an error!).
#' @param alloc.ratio Allocation of total sample size into the two groups. Must set as a vector of two
#'                    positive numbers. For equal allocation it should be put as c(1,1), for non-equal
#'                    allocation one can put c(2,1) or c(3,1) etc.
#' @param nsim The number of samples to be generated from the alternate distribution.
#' @return A named list with two elements.
#' \enumerate{
#' \item \code{samples} - a vector of length \code{nsim}, containing
#' The samples from the distribution of the Hotelling T statistic
#' under unequal variance.
#' \item \code{denom.df} - The denominator degrees of freedom of the chi-square statistic
#' obtained by approximation of the sum of two Wishart distribution under unequal variance.
#' }
#' @references Wang, Qiyao (2021)
#' \emph{Two-sample inference for sparse functional data,  Electronic Journal of Statistics,
#' Vol. 15, 1395-1423} \cr
#' \doi{https://doi.org/10.1214/21-EJS1802}.
#' @seealso [Hotelling::hotelling.test()], [Hotelling::hotelling.stat()] to generate empirical samples
#' from the Hotelling T-statistic from empirical data.
#' @export
#' @examples
#' \dontrun{
#' # In both the samples we assume that the covariances differ as the
#' # mean vector between the two groups differ.
#' # Case 1: Null hypothesis is true. True mean difference is zero, and the true
#' # covariance of the two groups are same.
#' k <- 5
#' mu1  <- rep(0,k); del  <- 0; mu2 <- mu1 + rep(del, k);
#' sig1 <- diag(k); sig2 <- sig1 + del*toeplitz(c(1,rep(0.5, k-1))); n <- 200;
#' null.dist.samples <- Sim_HotellingT_unequal_var(total_sample_size=n, mean_diff=mu1-mu2,
#'                      sig1=sig1, sig2=sig2, alloc.ratio=c(1,1), nsim=1e3)
#' # The following Kolmogorov Smirnov test confirms that under null hypothesis
#' # and when the covariances are same, the distribution is exactly a
#' # central F distribution with \eqn{k} and \eqn{n-k}  degrees of freedom.
#' ks.test(null.dist.samples$samples, {{(n - 2) * k}/(n - k -1)} * {rf(n=1e3, k, n-k-1)} )
#' }
#' \dontrun{
#' # Case 2: Alternate hypothesis is true. The mean difference is non-zero,
#' # and the covariances of the two groups are same:
#' k <- 6
#' mu1  <- rep(0,k); del  <- 0.15; mu2 <- mu1 + rep(del, k);
#' sig1 <- diag(k); sig2 <- sig1;
#' n1 <- 100; n2 <- 100;
#' alt.dist.samples <- Sim_HotellingT_unequal_var(total_sample_size=n1+n2, mean_diff=mu1-mu2,
#'                                                sig1=sig1, sig2=sig2, alloc.ratio=c(1,1), nsim=1e3)
#' ks.test(alt.dist.samples$samples,
#'         {(n1+n2 - 2) * k /(n1+n2 - k -1)}*rf(n=B, k, n1+n2-k-1,
#'           ncp = {(n1*n2)/(n1+n2)}*as.vector(crossprod(mu1-mu2, solve(sig1, mu1-mu2))) ) )
#' }
#' \dontrun{
#' # Case 3: Alternate hypothesis is true. The mean difference is non-zero,
#' # and the covariances of the two groups are different
#' k <- 5
#' mu1  <- rep(0,k); del  <- 0.25; mu2 <- mu1 + rep(del, k);
#' sig1 <- diag(k); sig2 <- sig1 + del*toeplitz(c(1,rep(0.5, k-1)))
#' alt.dist.samples <- Sim_HotellingT_unequal_var(total_sample_size=200, mean_diff=mu1-mu2,
#' sig1=sig1, sig2=sig2, alloc.ratio=c(1,1), nsim=1e3)
#'}
#'\dontrun{
#' # Generate samples with unequal allocation ratio:
#' k <- 8
#' mu1  <- rep(0,k); del  <- 0.4; mu2 <- mu1 + rep(del, k);
#' sig1 <- diag(k); sig2 <- sig1 + del*toeplitz(c(1,rep(0.5, k-1)))
#' alt.dist.samples <- Sim_HotellingT_unequal_var(total_sample_size=150, mean_diff=mu1-mu2,
#' sig1=sig1, sig2=sig2, alloc.ratio=c(2,1), nsim=1e3)
#'}
Sim_HotellingT_unequal_var <- function(total_sample_size, mean_diff,
                                         sig1, sig2, alloc.ratio=c(1,1), nsim=1e4){


  #*********************************************************************************************%
  #*                                  Missing argument checking                                 %
  #*********************************************************************************************%
  # set.seed(12345)
  call.f       <- rlang::call_match(defaults = FALSE)
  call.f.all   <- rlang::call_match(defaults = TRUE)
  default_args <- setdiff(rlang::call_args_names(call.f.all), rlang::call_args_names(call.f))
  if (length(default_args) > 0 ){
    rlang::inform(paste0("User did not specify values to the following arguments: ",
                         paste0(default_args, collapse = ","), " ---- default values are used"))
  }
  non_null_arg_vals <- !sapply(rlang::call_args(call.f), is.null)
  testthat::expect_true(all(non_null_arg_vals),
                        info = paste0("NULL values are specified for the arguments with no default values: ",
                                      paste0(rlang::call_args_names(call.f)[!non_null_arg_vals], collapse = ",")) )

  #*********************************************************************************************%
  #*                         Compatibility checking of arguments                                %
  #*********************************************************************************************%
  # argument checking : total_sample_size
  testthat::expect_true(rlang::is_double(total_sample_size, n=1, finite=TRUE) & (total_sample_size > 0),
                        info = "total_sample_size must be a positive integer")
  # argument checking : sig1
  testthat::expect_true(rlang::is_double(sig1, finite=TRUE) & isSymmetric(sig1),
                        info = "sig1 must be a symmetric matrix");
  # argument checking : sig2
  testthat::expect_true(rlang::is_double(sig2, n=length(sig1), finite=TRUE) & isSymmetric(sig2),
                        info = "sig2 must be a symmetric matrix of
                                dimension same as that of sig1");
  # argument checking : mean_diff
  testthat::expect_true(rlang::is_double(mean_diff, n=nrow(sig1), finite=TRUE),
                        info = "mean_diff must be numeric vector of
                                length same as dimension of sig1")
  # argument checking : alloc.ratio
  testthat::expect_true(rlang::is_double(alloc.ratio, n=2, finite=TRUE) & all(alloc.ratio > 0),
                        info = "alloc.ratio must be positive numeric vector of length 2");
  # argument checking : nsim
  testthat::expect_true(rlang::is_integerish(nsim, n=1, finite = TRUE) & (nsim > 0),
                         info = "nsim must be a positive integer")
  # argument checking : sig1, sig2
  testthat::expect_no_error(chol(sig1), message = "sig1 is not positive definite")
  testthat::expect_no_error(chol(sig2), message = "sig2 is not positive definite")

  #*********************************************************************************************%
  #*                                         Main body                                          %
  #*********************************************************************************************%
  mean_diff       <- as.vector(mean_diff)
  alloc.ratio     <- as.vector(alloc.ratio)
  B               <- nsim # Simulation replication
  k               <- length(mean_diff)
  n1              <- total_sample_size * {alloc.ratio[1]/sum(alloc.ratio)}
  n1_by_n2        <- alloc.ratio[1]/alloc.ratio[2]
  n2              <- total_sample_size - n1
  sig.til         <- sig1 + n1_by_n2 * sig2
  sqinv.sig.til   <- expm::sqrtm(solve(sig.til))
  delta           <- (sqinv.sig.til %*% sig1) %*% sqinv.sig.til
  I_minus_delta   <- diag(k)-delta

  delta.str       <- {n1_by_n2*(n1_by_n2 - 1/n2)}*delta + (1-1/n2)*I_minus_delta
  denom           <- {(n1_by_n2 - 1/n2)*(n1_by_n2^2)} *
                     {sum(diag(delta %*% delta)) + {sum(diag(delta))}^2} +
                     (1-1/n2) *
                     {sum(diag(I_minus_delta %*% I_minus_delta)) + {sum(diag(I_minus_delta))}^2}
  numer           <- sum(diag(delta.str %*% delta.str)) + {sum(diag(delta.str))}^2
  denom.df        <- n2 * numer/denom
  norma.vars      <- matrix(rnorm(B*k), B, k)
  ncp.vec         <- as.vector(sqrt(n1) * {sqinv.sig.til %*% mean_diff})
  norma.vars.e    <- sweep(norma.vars, 2, ncp.vec, FUN="+")
  chi_ncp         <- apply(norma.vars.e, 1, function(x) sum(x*as.vector(solve(delta.str, x))) )
  alt.dist.den    <- rchisq(n=B, df=denom.df-k+1)/denom.df
  samples         <- {chi_ncp/alt.dist.den}*{(n1_by_n2 + 1 - 2/n2)/(1 + 1/n1_by_n2)}

  ret.obj         <- c("samples", "denom.df")
  ret             <- lapply(1:length(ret.obj), function(u) get(ret.obj[u]))
  names(ret)      <- ret.obj
  return(ret)
}



#' CDF of Hotelling-\eqn{T^2} statistic.
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' @description
#' The function \code{pHotellingT()} computes the cumulative distribution function (CDF)
#' of the two-sample Hotelling-\eqn{T^2} statistic (\eqn{P(T > q)}) in the multivariate response
#' setting. This function is used to compute the power function
#' of Two-Sample (TS) Projection-based test (Wang 2021, EJS.)
#' for sparsely observed univariate functional data.
#'
#' @details Under the assumption of the equal variance, we know that the alternative
#' distribution of the Hotelling-\eqn{T^2} statistic (\eqn{(n-k-1)T/(n-2)*K}) has an
#' F distribution with the
#' non-centrality depending on the difference between the true mean vectors and the
#' (common) covariance of the response. However, when the true covariance of the true groups
#'  of responses differ, the alternate distribution becomes non-trivial. Koner and Luo (2023)
#'  proved that the alternate distribution of the test-statistic approximately follows
#'  a ratio of the linear combination of the K (dimension of the response) non-central
#'  chi-squared random variables (where the non-centrality parameter depends on the mean difference)
#'  and a chi-squared distribution whose degrees of freedom depends on a complicated functions of
#'  sample size in the two groups. This function initially calls the
#'  [fPASS::Sim_HotellingT_unequal_var] function to obtain the samples from the non-null distribution
#'  and computes the CDF numerically with high precision based on a large number of samples.
#'  See Koner and Luo (2023) for more details on the formula of the non-null distribution.
#' @import Matrix
#' @import lifecycle
#' @import rlang
#' @import dplyr
#' @importFrom testthat expect_true expect_gte expect_no_error
#' @importFrom MASS mvrnorm
#' @importFrom stats cov rnorm uniroot
#' @importFrom expm sqrtm
#' @author Salil Koner \cr Maintainer: Salil Koner
#' \email{salil.koner@@duke.edu}
#' @param q The point at which the CDF needs to be evaluated
#' @param total_sample_size Target sample size, must be a positive integer.
#' @param mean_diff The difference in the mean vector between the two groups, must be a vector.
#' @param sig1 The true (or estimate) of covariance matrix for the first group. Must be symmetric
#' (\code{is.symmetric(sig1) == TRUE}) and positive definite (\code{chol(sig1)} without an error!).
#' @param sig2 The true (or estimate) of covariance matrix for the second group. Must be symmetric
#' (\code{is.symmetric(sig2) == TRUE}) and positive definite (\code{chol(sig2)} without an error!).
#' @param alloc.ratio Allocation of total sample size into the two groups. Must set as a vector of two
#'                    positive numbers. For equal allocation it should be put as c(1,1), for non-equal
#'                    allocation one can put c(2,1) or c(3,1) etc.
#' @param lower.tail if TRUE, the CDF is returned, otherwise right tail probability is returned.
#' @return The CDF of the Hotelling T statistic, if \code{lower.tail == TRUE},
#' otherwise the right tail probability is returned.
#' @seealso [Hotelling::hotelling.test()], [Hotelling::hotelling.stat()] to generate empirical samples
#' from the Hotelling T-statistic from empirical data.
#' @export
#' @examples
#' \dontrun{
#' B           <- 10000
#' k           <- 4
#' n2          <- 60
#' n1_by_n2    <- 2
#' n1          <- n1_by_n2 * n2
#' mu1         <- rep(0,k)
#' del         <- 0.4
#' mu2         <- mu1 + rep(del, k) # rep(0.19,k)  # 0.23 (0.9), 0.18 (0.7) 0.20 (0.8)
#' sig1        <- diag(k)
#' sig2        <- sig1
#' cutoff      <- seq(0,30, length.out=20)
#' the_cdf     <- round(pHotellingT(cutoff, n1+n2, mu1 - mu2,
#'                                  sig1, sig2, alloc.ratio=c(2,1),
#'                                  lower.tail=FALSE),3)
#' }
pHotellingT <- function(q, total_sample_size, mean_diff,
                            sig1, sig2, alloc.ratio=c(1,1),
                            lower.tail=TRUE){

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
  non_null_arg_vals <- !sapply(rlang::call_args(call.f), is.null)
  testthat::expect_true(all(non_null_arg_vals),
                        info = paste0("NULL values are specified for the arguments with no default values: ",
                                      paste0(rlang::call_args_names(call.f)[!non_null_arg_vals], collapse = ",")) )

  #*********************************************************************************************%
  #*                              compatibility check of arguments                              %
  #*********************************************************************************************%
  # argument checking : q
  testthat::expect_true(rlang::is_double(q, finite=TRUE),
                        info = "q must be a numeric vector");
  # argument checking : total_sample_size
  testthat::expect_true(rlang::is_double(total_sample_size, n=1, finite=TRUE) & (total_sample_size > 0),
                        info = "total_sample_size must be a positive integer")
  # argument checking : sig1
  testthat::expect_true(rlang::is_double(sig1, finite=TRUE) & isSymmetric(sig1),
                        info = "sig1 must be a symmetric matrix");
  # argument checking : sig2
  testthat::expect_true(rlang::is_double(sig2, n=length(sig1), finite=TRUE) & isSymmetric(sig2),
                        info = "sig2 must be a symmetric matrix of
                                dimension same as that of sig1");
  # argument checking : mean_diff
  testthat::expect_true(rlang::is_double(mean_diff, n=nrow(sig1), finite=TRUE),
                        info = "mean_diff must be numeric vector of
                                length same as dimension of sig1")
  # argument checking : alloc.ratio
  testthat::expect_true(rlang::is_double(alloc.ratio, n=2, finite=TRUE) & all(alloc.ratio > 0),
                        info = "alloc.ratio must be positive numeric vector of length 2");
  # argument checking : lower.tail
  testthat::expect_true(rlang::is_logical(lower.tail),
                        info = "Argument lower.tail must a logical with length 1")
  # argument checking : sig1, sig2
  testthat::expect_no_error(chol(sig1), message = "sig1 is not positive definite")
  testthat::expect_no_error(chol(sig2), message = "sig2 is not positive definite")

  #*********************************************************************************************%
  #*                                         Main body                                          %
  #*********************************************************************************************%
  if (norm(sig1-sig2) < 1e-12){
    sig.common   <- (sig1 + sig2)/2
    k            <- length(mean_diff)
    lamda        <- as.numeric(crossprod(mean_diff, solve(sig.common, as.vector(mean_diff))))
    prob         <- pf({(total_sample_size-k-1)*as.vector(q)}/{(total_sample_size-2)*k},
                       df1=k, df2=total_sample_size-k-1,
                       ncp = {total_sample_size/{sum(alloc.ratio) * sum(1/alloc.ratio)}}*lamda,
                       lower.tail = lower.tail) # power for n
  } else{
    alt.samples  <- Sim_HotellingT_unequal_var(total_sample_size=total_sample_size,
                                               mean_diff=mean_diff,
                                               sig1=sig1, sig2=sig2,
                                               alloc.ratio=alloc.ratio, nsim=1e4)
    cdf          <- sapply(as.vector(q), function(cc) mean(alt.samples$samples <= cc, na.rm=TRUE))
    if(lower.tail) prob <- cdf
    else prob <- 1 - cdf

  }
  testthat::expect_true(rlang::is_double(prob, n=length(q), finite=TRUE) &
                        all(prob >= 0) & all(prob <= 1),
                        info = "prob must a numeric vector of length same as q lying between 0 and 1")
  return(prob)
}
