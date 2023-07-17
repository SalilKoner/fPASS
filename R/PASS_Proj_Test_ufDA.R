#' Power and Sample size (PASS) calculation of
#' Two-Sample Projection-based test for sparsely observed univariate functional data.
#' @description
#'
#' `r lifecycle::badge("stable")`
#'
#' @description
#'
#' The function `PASS_Proj_Test_ufDA()` computes the power and sample size (PASS) required to conduct
#' the projection-based test of mean function between two groups of longitudinal data
#' or sparsely observed functional data under a random irregular design, under
#' common covariance structure between the groups. See Wang (2021) for more details
#' of the testing procedure.
#'
#' @details The function is designed to perform the power and sample size analysis
#' for functional under a dense and sparse (random) design and longitudinal data. The
#' function can handle data from wide variety of covariance structure, can be parametric,
#' or non-parametric. Additional with traditional stationary structures assumed for longitudinal
#' data (see [nlme::corClasses]), the user can specify any other non-stationary covariance
#' function in the form of either a covariance function or in terms of eigenfunctions and
#' eigenvalues. The user have a lot of flexibility into tweaking the arguments of the function
#' to assess the power function of the test under different sampling design and covariance process
#' of the response trajectory, and for any arbitrary mean difference function. Overall, the
#' functionality of the module is quite comprehensive and includes all the different cases
#' considered in the 'NCSS PASS (2023)' software. We believe that this software can be an effective
#' clinical trial design tools when considering the projection-based test as the primary
#' decision making method.
#'
#' @inheritSection Extract_Eigencomp_fDA Specification of key arguments
#' @author Salil Koner \cr Maintainer: Salil Koner
#' \email{salil.koner@@duke.edu}
#' @import lifecycle
#' @importFrom testthat expect_true expect_gte expect_no_error expect_named
#' @param sample_size Total sample size combining both the groups, must be a positive integer.
#' @param target.power Target power to achieve, must be a number between 0 and 1.
#'                     Only one of \code{sample_size} and \code{target.power} should be non-null. The
#'                     function will return sample size if \code{sample_size} is NULL, and
#'                     return power if \code{target.power} is NULL.
#' @param sig.level Significance level of the test, default set at 0.05, must be less than 0.2.
#' @inheritParams Extract_Eigencomp_fDA
#' @param nWgrid The length of the working grid based in the domain of the function on which
#'        the eigenfunctions will be estimated. The actual working grid will be calculated using
#'        the [gss::gauss.quad()] function (so that it facilitates the numerical integration of
#'        the eigenfunction with the mean function using gaussian quadrature rule)
#' @param npc_to_use Number of eigenfunctions to use to compute the power. Default is NULL, in
#' which case all the eigenfunctions estimated from the data will be used.
#' @param return.eigencomp Indicates whether to return the eigencomponents obtained from the fPCA
#' on the large data with sample size equal to \code{eval_SS}. Default is \code{FALSE}.
#' @param nsim The number of samples to be generated from the alternate distribution of
#'             Hotelling T statistic. Default value is 10000.
#' @return A list with following elements, \code{power_value} if \code{is.null(target.power)}
#' then returns the power of the test when n equal to \code{sample_size}, otherwise \code{required_SS},
#' the sample size required to achieve the power of the test at \code{target.power}.
#' If \code{return.eigencomp == TRUE} then \code{est_eigencomp} is also returned, containing
#'  the entire output obtained from internal call of [fPASS::Extract_Eigencomp_fDA()].
#' @seealso See [fPASS::Power_Proj_Test_ufDA()] and [fPASS::Extract_Eigencomp_fDA()].
#' @references Wang, Qiyao (2021)
#' \emph{Two-sample inference for sparse functional data,  Electronic Journal of Statistics,
#' Vol. 15, 1395-1423}
#' \doi{https://doi.org/10.1214/21-EJS1802}. \cr \cr
#' PASS 2023 Power Analysis and Sample Size Software (2023). NCSS, LLC. Kaysville, Utah, USA, ncss.com/software/pass.
#' @export PASS_Proj_Test_ufDA
#' @examples
#'
#' # Example 1: Power analysis for stationary exponential covariance.
#' # Should return a power same as the size because
#' # the true mean difference is zero.
#'
#' set.seed(12345)
#' mean.diff <- function(t) {0*t};
#' obs.design = list("design" = "longitudinal",
#'                   "visit.schedule" = seq(0.1, 0.9, length.out=7),
#'                    "visit.window" = 0.05)
#' cor.str <- nlme::corExp(1, form = ~ time | Subject);
#' sigma2 <- 1; sigma2.e <- 0.25; nobs_per_subj <- 8;
#' missing_type <- "constant"; missing_percent <- 0.01;
#' # Please increase `eval_SS` argument from 1000 to 5000 to get
#' # accurate precision on the estimated eigenfunctions.
#' pow  <- PASS_Proj_Test_ufDA(sample_size = 100, target.power = NULL, sig.level = 0.05,
#'                             obs.design = obs.design,
#'                             mean_diff_fnm = "mean.diff", cov.type = "ST",
#'                             cov.par = list("var" = sigma2, "cor" = cor.str),
#'                             sigma2.e = sigma2.e, nobs_per_subj = nobs_per_subj,
#'                             missing_type = missing_type,
#'                             missing_percent = missing_percent, eval_SS = 1000,
#'                             alloc.ratio = c(1,1), nWgrid = 201,
#'                             fpca_method = "fpca.sc",
#'                             mean_diff_add_args = list(), fpca_optns = list("pve" = 0.95),
#'                             nsim = 1e3)
#'
#' print(pow$power_value)
#'
#'# Example 2: Sample size calculation for a non-stationary covariance:
#'
#' alloc.ratio  <- c(1,1)
#' mean.diff    <- function(t) {3 * (t^3)};
#' eig.fun <- function(t, k) {
#'   if (k==1) ef <- sqrt(2)*sin(2*pi*t)
#'   else if (k==2) ef <- sqrt(2)*cos(2*pi*t)
#'   return(ef)}
#' eig.fun.vec  <- function(t){cbind(eig.fun(t, 1),eig.fun(t, 2))}
#' eigen.comp   <- list("eig.val" = c(1, 0.5), "eig.obj" = eig.fun.vec)
#' obs.design   <- list(design = "functional", fun.domain = c(0,1))
#' cov.par      <- list("cov.obj" = NULL, "eigen.comp" = eigen.comp)
#' sigma2.e     <- 0.001; nobs_per_subj <- 4:7;
#' missing_type <- "nomiss"; missing_percent <- 0;
#' fpca_method  <- "fpca.sc"
#' # Please increase `eval_SS` argument from 1000 to 5000 to get
#' # accurate precision on the estimated eigenfunctions.
#' pow  <- PASS_Proj_Test_ufDA(sample_size = NULL, target.power = 0.8,
#'                             sig.level = 0.05, obs.design = obs.design,
#'                             mean_diff_fnm = "mean.diff", cov.type = "NS",
#'                             cov.par = cov.par, sigma2.e = sigma2.e,
#'                             nobs_per_subj = nobs_per_subj, missing_type = missing_type,
#'                             missing_percent = missing_percent, eval_SS = 1000,
#'                             alloc.ratio = alloc.ratio, fpca_method = "fpca.sc",
#'                             mean_diff_add_args = list(), fpca_optns = list(pve = 0.95),
#'                             nsim = 1e3, nWgrid = 201)
#'
#' print(pow$required_SS)
#'
PASS_Proj_Test_ufDA  <- function(sample_size, target.power, sig.level = 0.05,
                                 nobs_per_subj, obs.design, mean_diff_fnm,
                                 cov.type = c("ST", "NS"), cov.par, sigma2.e,
                                 missing_type = c("nomiss", "constant"),
                                 missing_percent = 0,
                                 eval_SS = 5000, alloc.ratio = c(1,1),
                                 fpca_method = c("fpca.sc", "face"),
                                 mean_diff_add_args=list(),
                                 fpca_optns = list(pve = 0.95),
                                 nWgrid = 201,
                                 npc_to_use = NULL, return.eigencomp = FALSE,
                                 nsim = 1e4){

  testthat::expect_equal(xor(!is.null(sample_size), !is.null(target.power)), TRUE,
                         info = "Only one of sample_size or target.power should be NULL")
  if(!is.null(sample_size)){
    message("Computing power of Projection-based test for total sample size = ", sample_size, "\n") # Added by SK on Jan 17
    # argument checking: total_sample_size
    testthat::expect_true(rlang::is_integerish(sample_size, n=1, finite = TRUE) & (sample_size > 0),
                          info = "sample_size must be a positive integer with value greater than 10")
  }
  if(!is.null(target.power)){
    message("Computing sample size required to achieve power = ", target.power, "\n") # Added by SK on Jan 17
    testthat::expect_true(rlang::is_double(target.power, n=1, finite = TRUE) & (target.power <= 1) &
                          (target.power > 0),
                          info = "target.power must be a positive number between 0 and 1.")
  }
  testthat::expect_true(rlang::is_scalar_logical(return.eigencomp), info = "return.eigencomp must be logical")
  assign(mean_diff_fnm, match.fun(mean_diff_fnm))

  testthat::expect_true(rlang::is_list(obs.design))
  testthat::expect_true("design" %in% names(obs.design),
                        info = "obs.design must be a named list with name 'design'")
  design <- match.arg(obs.design$design, c("longitudinal", "functional"))
  if (design == "functional") {
    testthat::expect_true(rlang::is_integerish(nobs_per_subj, finite=TRUE),
                          info = "nobs_per_subj must be a positive integer scalar/vector
                          (for varying number of obs)")
    testthat::expect_named(obs.design[!names(obs.design) %in% "design"], "fun.domain", ignore.order = TRUE,
                           info = "Names of obs.design must be 'design' and 'fun.domain'.")
    # argument checking: interval
    interval <- obs.design[["fun.domain"]]
    testthat::expect_true(rlang::is_double(interval, n=2, finite=TRUE) & !is.unsorted(interval),
                          info = "obs.design$interval of the observation points must be
                          two-length numeric vector")
  } else{
    interval <- c(0,1)
  }
  testthat::expect_true(rlang::is_integerish(nWgrid, n=1, finite=TRUE) & (nWgrid > 1),
                        info = "nWgrid must be a positive integer greater than 1.")
  if (nWgrid < 20){
    warning("nWgrid is less than 20, consider increasing it to at least 20.")
  }
  gauss.quad.pts <- gss::gauss.quad(nWgrid, interval)
  work.grid      <- gauss.quad.pts$pt

  est_eigencomp <- Extract_Eigencomp_fDA(obs.design = obs.design, mean_diff_fnm = mean_diff_fnm, cov.type = cov.type,
                                         cov.par = cov.par, sigma2.e = sigma2.e, nobs_per_subj = nobs_per_subj,
                                         missing_type = missing_type, missing_percent = missing_percent,
                                         eval_SS = eval_SS, alloc.ratio = alloc.ratio, fpca_method = fpca_method,
                                         data.driven.scores = FALSE, mean_diff_add_args = mean_diff_add_args,
                                         fpca_optns = fpca_optns, work.grid = work.grid, nWgrid = nWgrid)
  if(is.null(npc_to_use)) npc_to_use <- ncol(est_eigencomp$est_eigenfun) else{
    npc_to_use  <- min(npc_to_use, ncol(est_eigencomp$est_eigenfun))
  }
  min_search <- uniroot(function(n){
    Sum_of_Wishart_df(total_sample_size = n, alloc.ratio = alloc.ratio,
                             sig1 = est_eigencomp$score_var1[1:npc_to_use, 1:npc_to_use, drop=FALSE],
                             sig2 = est_eigencomp$score_var2[1:npc_to_use, 1:npc_to_use, drop=FALSE]) -
      npc_to_use + 1}, c(3, 5000), tol = .Machine$double.eps^0.75, extendInt = "upX")$root
  # cat("min_search = ", min_search, "\n")
  if (is.null(sample_size)){
    required_SS <- uniroot(function(n){
      Power_Proj_Test_ufDA(total_sample_size = n, argvals = est_eigencomp$working.grid,
                           mean_vector = est_eigencomp$mean_diff_vec, eigen_matrix = est_eigencomp$est_eigenfun,
                           scores_var1 = est_eigencomp$score_var1, scores_var2 = est_eigencomp$score_var2,
                           weights = gauss.quad.pts$wt, sig.level=sig.level, alloc.ratio = alloc.ratio,
                           npc_to_pick = npc_to_use, nsim = nsim) - target.power
    }, c(max(3,min_search + 1, npc_to_use+2), 1000), tol = .Machine$double.eps^0.25, extendInt = "upX")$root
  } else if (is.null(target.power)){
    power_value  <- Power_Proj_Test_ufDA(total_sample_size = sample_size, argvals = est_eigencomp$working.grid,
                                         mean_vector = est_eigencomp$mean_diff_vec, eigen_matrix = est_eigencomp$est_eigenfun,
                                         scores_var1 = est_eigencomp$score_var1, scores_var2 = est_eigencomp$score_var2,
                                         weights = gauss.quad.pts$wt, sig.level=sig.level, alloc.ratio = alloc.ratio,
                                         npc_to_pick = npc_to_use, nsim = nsim)
  }
  if (is.null(sample_size)) ret.objects <- "required_SS" else ret.objects <- "power_value"
  if (return.eigencomp){
    ret.objects <- c(ret.objects, "est_eigencomp")
  }
  ret.val        <- lapply(ret.objects, function(obj) get(obj))
  names(ret.val) <- ret.objects
  ret.val
}
