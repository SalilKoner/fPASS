#' Power and Sample size (PASS) calculation of
#' Two-Sample (TS) Projection-based test (ProjTest) for sparsely observed univariate functional data.
#' @description
#'
#' `r lifecycle::badge("experimental")`
#'
#' @description
#'
#' This function `PASS_Proj_Test_ufDA()` computes the power and sample size required to conduct
#' the projection-based test of mean function between two groups of longitudinal data
#' or sparsely observed functional data under a random irregular design, under
#' common covariance structure between the groups.
#' @details The projection-based test assumes the two groups of the data share the identical
#' covariance structure. We allow the covariance structure of the data as both parametric
#' and non-parametric. If the covariance structure is non-parametric, then, the user must specify the
#' the covariance function in the form of a function.
#' The projection-based test represents the sparsely observed
#' functions parsimoniously in terms of Karhunen-Loeve (KL) expansion and use the
#' functional principal component analysis scores to test the difference in the mean
#' function between two groups, by taking advantage of multivariate Hotelling-\eqn{T^2} test.
#' See Wang (2021) for more details of the testing procedure. We use the `fpca.sc()` which
#' is a just a copy of `refund::fpca.sc()` function to account for the shrinkage correctly,
#' or `face::face.sparse()` function to conduct the functional principal component analysis (fPCA).
#'
#' ## Details on the specification of arguments.
#'
#' If `obs.design$design == "functional"` then a dense grid of length, specified by ngrid (typically 101/201) is internally created, and
#' the observation points will be randomly chosen from them. The time points could also randomly chosen between
#' any number between the interval, but then for large number of subject, [fpca.sc()] function will take huge
#' time to estimate the eigenfunction. For dense design, the user must set a large value of the argument
#' `nobs_per_subj` and for sparse (random) design, `nobs_per_subj` should be set small (and varying).
#' On the other hand, typical to longitudinal data, if the measurements are taken at fixed time points (from baseline)
#' for each subject, then the user must set `obs.design$design == "longitudinal"` and the time points must be accordingly specified
#' in the argument `obs.design$visit.schedule`. The length of `obs.design$visit.schedule` must match `length(nobs_per_subj)-1`.
#  Internally for any case, when `design == "longitudinal"`, we will scale the visit times
#' so that it lies between \eqn{[0,1]}, so user should not specify any "fun.domain" in the
#' list for obs.design$design == "longitudinal". Make sure that the mean_function and the covariance function specified
#' in the `cov.par` and `mean_diff_fnm` parameter also scaled to take argument between \eqn{[0,1]}.
#' Also, it is imperative to say that `nobs_per_subj` must be of a scalar positive inter for `design == "longitudinal"`.
#'
#' @author Salil Koner \cr Maintainer: Salil Koner
#' \email{salil.koner@@duke.edu}
#' @import lifecycle
#' @importFrom testthat expect_true expect_gte expect_no_error expect_named
#' @param sample_size Total sample size combining both the groups, must be a positive integer more than 3.
#' @param target.power Target power, must be a number between 0 and 1.
#'                     Only one of sample_size and target.power should be non-null. The
#'                     function will automatically solve for sample size if sample_size is NULL, and
#'                     solve for power if target.power is NULL.
#' @param sig.level Significance level of the test, default set at 0.05, must be less than 0.2.
#'                  This is used to compute the critical value of the test.
#' @param nobs_per_subj The number of observations per subject. Each element of it greater than 3.
#' Could be a vector to specify that the number of observation for each is randomly varying
#' between the elements of the vector, or a scalar to ensure that the number of observations are same
#' for each subject.
#' @param obs.design The sampling design of the observations. Must be provided as a list with the following elements.
#' For functional design it must be provided as a named list with elements "design" and "fun.domain", where for example,
#' design = "functional" and `fun.domain = c(0,1)` specifying the domain of the functions. For longitudinal design,
#' it must be a named list with elements "design", "visit.schedule" and "visit.window", where the last two elements
#' specifying schedule of visits (in months or days or any scale), other than the baseline visit
#' and the maximum time window for every visit, which is also for every visit other than the baseline. See examples.
#' @param mean_diff_fnm The name of the function representing the group difference. Must be supplied as character.
#' @param cov.type  The type of the covariance structure of the data, must be either of "ST" (stationary) or
#'                  "NS" (non-stationary)
#' @param cov.par The covariance parameters that needs to be specified.
#'                If `cov.type='ST'` then, cov.par
#'                must be a named list of two elements, "var" and "cor", where "var" denotes
#'                the common variance of the observations, must be a positive number; and "cor" denotes
#'                The correlation structure between the observations.
#'                It must be provided in the form specified in the
#'                available in the documentation of R package \pkg{nlme}. Check the package documentation for more details.
#'                The argument of this function is passed onto the [nlme::corMatrix()] to extract the subject-specific
#'                covariance matrix.
#'                If `cov.type='NS'` then, cov.par
#'                must be a named list of two elements, "cov.obj" and "eigen.comp".
#'                Only one of the "cov.obj" or "eigen.comp" must be non-null. If the "cov.obj" is
#'                specified, then it must be either a bivariate function. Alternatively,
#'                if the true eigenfunctions (and number of eigenfunctions, say K) are known,
#'                then the user can specify that by specifying eigen.comp.
#'                In this case, the eigen.comp must be a named list with two elements,
#'                "eig.obj" and "eig.val".
#'                The object "eig.obj" must be a function so that the code `eigen.comp$eig.obj(timepoints)`
#'                returns a matrix of dimension r by K, where r is the length of `timepoints` and K
#'                is the number of eigenfunctions.
#' @param sigma2.e Measurement error variance, should be set as zero if there is no measurement error.
#' @param missing_type The type of missing in the number of observations of the subjects. Can be one of
#'                     \code{missing_type = "nomiss"} for no missing observations
#'                     or \code{missing_type = "constant"} for constant
#'                     missing. Only supports \code{missing_type = "constant"} now. Note that, if
#'                     nobs_per_subj is specified as vector,
#'                     `missing_type` is forced to set as "nomiss" and \code{missing_percent = 0}, because
#'                     the \code{missing_type = "constant"} has no meaning if the number of observations are
#'                     varying between the subject at the first, typically considered in
#'                     the case of sparse random functional design.
#' @param missing_percent The percentage of missing at each observation points for each subject.
#' @param eval_SS The sample size based on which the eigenfunctions will be estimated from data.
#' To compute the theoretical power of the test we must make sure that we use a large enough sample size
#' to generate the data such that the estimated eigenfunctions are very close to the true eigenfunctions
#' and that the sampling design will not have much effect on the loss of precision. Default value 5000.
#' @param alloc.ratio The allocation ratio of samples in the each group. Note that the eigenfunctions will still
#' be estimated based on the total sample_size, however, the variance of the shrinkage scores will be
#' estimated based on the allocation of the samples in each group. Must be given as vector of
#' length 2. Default value is set at c(1,1), indicating equal sample size.
#' @param fpca_method The method by which the FPCA is computed. Must be one of c("fpca.sc", "face").
#' @param mean_diff_add_args Additional arguments to be passed to the mean_diff
#'                           function specified by the name `mean_diff_fnm`.
#' @param fpca_optns Additional options to be passed onto either of `refund::fpca.sc()`
#'                  or `face::face.sparse()` function in order
#'                   to estimate the eigencomponents. It must be a named list with elements
#'                   to be passed onto the respective function, depending on the `fpca_method`.
#' @param npc_to_use Number of eigenfunctions to use to compute the power.
#' @return Returns the estimated eigenfunctions and the covariance of the shrinkage scores for the
#'         for the two groups, as well as the combined covariances assuming that they have the same o
#'         covariances.
#' @seealso See [fPASS::Power_Proj_Test_ufDA()] and [fPASS::pHotellingT()].
#' @references Wang, Qiyao (2021)
#' \emph{Two-sample inference for sparse functional data,  Electronic Journal of Statistics,
#' Vol. 15, 1395-1423} \cr
#' \doi{https://doi.org/10.1214/21-EJS1802}.
#' @export PASS_Proj_Test_ufDA
#' @examples
#'
#' \dontrun{
#' set.seed(12345)
#' mean.diff <- function(t) {t};
#' obs.design = list("design" = "functional",
#' "visit.schedule" = seq(0.1, 0.9, length.out=7), "visit.window" = 0.05)
#' cor.str <- nlme::corExp(1, form = ~ time | Subject);
#' sigma2 <- 1; sigma2.e <- 0.25; nobs_per_subj <- 6;
#' missing_type <- "constant"; missing_percent <- 0.01;
#' eigencomp  <- PASS_Proj_Test_ufDA(obs.design = obs.design,
#'  mean_diff_fnm = "mean.diff", cov.type = "ST",
#'    cov.par = list("var" = sigma2, "cor" = cor.str),
#'    sigma2.e = sigma2.e, nobs_per_subj = nobs_per_subj,
#'    missing_type = missing_type,
#'    missing_percent = missing_percent, eval_SS = 5000,
#'   alloc.ratio = c(1,1),
#'   fpca_method = "fpca.sc",
#'   mean_diff_add_args = list(), fpca_optns = list("pve" = 0.95))
#'
#'}
#'\dontrun{
#' alloc.ratio  <- c(1,1)
#' mean.diff    <- function(t) {1 * (t^3)};
#' eig.fun <- function(t, k) {
#'   if (k==1) ef <- sqrt(2)*sin(2*pi*t)
#'   else if (k==2) ef <- sqrt(2)*cos(2*pi*t)
#'   return(ef)}
#' eig.fun.vec  <- function(t){cbind(eig.fun(t, 1),eig.fun(t, 2))}
#' eigen.comp   <- list("eig.val" = c(1, 0.5), "eig.obj" = eig.fun)
#' obs.design   <- list(design = "functional", fun.domain = c(0,1))
#' cov.par      <- list("cov.obj" = NULL, "eigen.comp" = eigen.comp)
#' sigma2.e     <- 0.001; nobs_per_subj <- 4:7;
#' missing_type <- "nomiss"; missing_percent <- 0;
#' fpca_method  <- "fpca.sc"
#' eigencomp  <- PASS_Proj_Test_ufDA(obs.design = obs.design,
#' mean_diff_fnm = "mean.diff", cov.type = "NS",
#' cov.par = cov.par, sigma2.e = sigma2.e,
#' nobs_per_subj = nobs_per_subj, missing_type = missing_type,
#' missing_percent = missing_percent, eval_SS = 5000,
#' alloc.ratio = alloc.ratio, fpca_method = "fpca.sc",
#' mean_diff_add_args = list(), fpca_optns = list(pve = 0.95))
#'}
PASS_Proj_Test_ufDA  <- function(sample_size, target.power, sig.level = 0.05,
                                 nobs_per_subj, obs.design, mean_diff_fnm,
                                 cov.type = c("ST", "NS"), cov.par, sigma2.e,
                                 missing_type = c("nomiss", "constant"),
                                 missing_percent = 0,
                                 eval_SS = 5000, alloc.ratio = c(1,1),
                                 fpca_method = c("fpca.sc", "face"),
                                 mean_diff_add_args=list(),
                                 fpca_optns = list(pve = 0.95),
                                 npc_to_use = 2){

  testthat::expect_equal(xor(!is.null(sample_size), !is.null(target.power)), TRUE,
                         info = "Only one of sample_size or target.power should be NULL")
  if(!is.null(sample_size)){
    cat("Computing power of Projection-based test for total sample size = ", sample_size, "\n") # Added by SK on Jan 17
    # argument checking: total_sample_size
    testthat::expect_true(rlang::is_integerish(sample_size, n=1, finite = TRUE) & (sample_size > 10),
                          info = "sample_size must be a positive integer with value greater than 10")
  }
  if(!is.null(target.power)){
    cat("Computing sample size required to achieve power = ", target.power, "\n") # Added by SK on Jan 17
    testthat::expect_true(rlang::is_double(target.power, n=1, finite = TRUE) & (target.power <= 1) &
                          (target.power > 0),
                          info = "target.power must be a positive number between 0 and 1.")
  }
  assign(mean_diff_fnm, match.fun(mean_diff_fnm))
  est_eigencomp <- Extract_Eigencomp_fDA(obs.design = obs.design, mean_diff_fnm = mean_diff_fnm, cov.type = cov.type,
                                         cov.par = cov.par, sigma2.e = sigma2.e, nobs_per_subj = nobs_per_subj,
                                         missing_type = missing_type, missing_percent = missing_percent,
                                         eval_SS = eval_SS, alloc.ratio = alloc.ratio, fpca_method = fpca_method,
                                         data.driven.scores = FALSE, mean_diff_add_args = mean_diff_add_args,
                                         fpca_optns = fpca_optns)
  if (is.null(sample_size)){
    required_SS <- uniroot(function(n){
      Power_Proj_Test_ufDA(total_sample_size = n, argvals = est_eigencomp$working.grid,
                           mean_vector = est_eigencomp$mean_diff_vec, eigen_matrix = est_eigencomp$est_eigenfun,
                           scores_var1 = est_eigencomp$score_var1, scores_var2 = est_eigencomp$score_var2,
                           weights = est_eigencomp$weights, sig.level=sig.level, alloc.ratio = alloc.ratio,
                           npc_to_pick = npc_to_use) - target.power
    }, c(npc_to_use+1, 1e+07), tol = .Machine$double.eps^0.25, extendInt = "upX")$root
  } else if (is.null(target.power)){
    power_value  <- Power_Proj_Test_ufDA(total_sample_size = sample_size, argvals = est_eigencomp$working.grid,
                                         mean_vector = est_eigencomp$mean_diff_vec, eigen_matrix = est_eigencomp$est_eigenfun,
                                         scores_var1 = est_eigencomp$score_var1, scores_var2 = est_eigencomp$score_var2,
                                         weights = est_eigencomp$weights, sig.level=sig.level, alloc.ratio = alloc.ratio,
                                         npc_to_pick = npc_to_use)
  }
  ret.objects    <- "est_eigencomp"
  if (is.null(sample_size)) ret.objects <- c(ret.objects, "required_SS") else ret.objects <- c(ret.objects, "power_value")
  ret.val        <- lapply(ret.objects, function(obj) get(obj))
  names(ret.val) <- ret.objects

  ret.val
}
