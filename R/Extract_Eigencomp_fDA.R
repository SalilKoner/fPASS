
#' Extract/estimate eigenfunction from a sparse functional or longitudinal design
#' by simulating from a large number of subjects.
#' @description
#'
#' `r lifecycle::badge("experimental")`
#'
#' @description
#'
#' The function `Extract_Eigencomp_fDA()` computes the eigenfunctions and the
#' covariance of the shrinkage scores required to conduct
#' the projection-based test of mean function between two groups of longitudinal data
#' or sparsely observed functional data under a random irregular design, as developed by Wang (2021).
#'
#' @details The function can handle data from wide variety of covariance structure, can be parametric,
#' or non-parametric. Additional with traditional stationary structures assumed for longitudinal
#' data (see [nlme::corClasses]), the user can specify any other non-stationary covariance
#' function in the form of either a covariance function or in terms of eigenfunctions and
#' eigenvalues. The user have a lot of flexibility into tweaking the arguments \code{nobs_per_subject},
#' \code{obs.design}, and \code{cov.par} to compute the eigencomponents
#' under different sampling design and covariance process of the response trajectory, and
#' for any arbitrary mean difference function. Internally, using the sampling
#' design and the covariance structure specified, we generate a large data with
#' large number of subjects, and estimate the eigenfunctions and the covariance of the estimated
#' `shrinkage` scores by means of functional principal component analysis (fPCA). We put the option of using
#' two most commonly used softwares for fPCA in the functional data literature, [refund::fpca.sc()]
#' and [face::face.sparse()]. However, since the [refund::fpca.sc()] do not compute the `shrinkage`
#' scores correctly, especially when the measurement error variance is estimated to be zero,
#' we made a duplicate version of that version where we write out the scoring part on our own.
#' The new function is named as [fPASS::fpca.sc()], please check it out.
#'
#' ## Details on the specification of arguments.
#'
#' If \code{obs.design$design == 'functional'} then a dense grid of length,
#' specified by ngrid (typically 101/201) is internally created, and
#' the observation points will be randomly chosen from them.
#' The time points could also randomly chosen between
#' any number between the interval, but then for large number of subject,
#' [fPASS::fpca.sc()] function will take huge
#' time to estimate the eigenfunction. For dense design, the user must set
#' a large value of the argument \code{nobs_per_subj} and for sparse (random) design,
#' \code{nobs_per_subj} should be set small (and varying).
#' On the other hand, typical to longitudinal data, if the measurements are
#' taken at fixed time points (from baseline)
#' for each subject, then the user must set \code{obs.design$design == 'longitudinal'} and
#' the time points must be accordingly specified
#' in the argument \code{obs.design$visit.schedule}. The length of \code{obs.design$visit.schedule}
#' must match \code{length(nobs_per_subj)-1}. Internally, when
#' \code{obs.design$design == 'longitudinal'}, the function scale the visit times
#' so that it lies between \eqn{[0,1]}, so the user should not
#' specify any element named \code{fun.domain} in the
#' list for \code{obs.design$design == 'longitudinal'}. Make sure that
#'  the mean function and the covariance function specified
#' in the \code{cov.par} and \code{mean_diff_fnm} parameter also scaled to
#' take argument between \eqn{[0,1]}. Also, it is imperative to say that `nobs_per_subj` must
#' be of a scalar positive integer for \code{design == 'longitudinal'}.
#'
#' @author Salil Koner \cr Maintainer: Salil Koner
#' \email{salil.koner@@duke.edu}
#' @import Matrix
#' @import rlang
#' @import lifecycle
#' @import dplyr
#' @importFrom face face.sparse
#' @import mgcv
#' @importFrom testthat expect_true expect_gte expect_no_error expect_named
#' @importFrom nlme corMatrix Initialize getCovariateFormula getGroupsFormula
#' @importFrom MASS mvrnorm
#' @importFrom stats cov rnorm complete.cases predict rchisq weighted.mean pf qf
#' @importFrom stringr str_split_1
#' @importFrom gss gauss.quad
#' @importFrom purrr map2
#' @param nobs_per_subj The number of observations per subject. Each element of it must be greater than 3.
#' It could also be a vector to indicate that the number of observation for each is randomly varying
#' between the elements of the vector, or a scalar to ensure that the number of observations are same
#' for each subject. See examples.
#' @param obs.design The sampling design of the observations. Must be provided as
#' a list with the following elements. If the design is longitudinal (e.g. a clinical trial
#' where there is pre-specified schedule of visit for the participants) it must be
#' a named list with elements \code{design}, \code{visit.schedule} and \code{visit.window}, where
#' \code{obs.design$design} must be specified as `'longitudinal'`, \code{visit.schedule}
#' specifying schedule of visits (in months or days or any unit of time), other than the baseline visit
#' and \code{visit.window} denoting the maximum time window for every visit.
#' For functional design (where the observation points are either densely observed within a
#' compact interval or under a sparse random design), the argument must be provided
#' as a named list with elements\code{design} and \code{fun.domain}, where
#' \code{obs.design$design} must be specified as `'functional'` and \code{obs.design$fun.domain}
#' must be specified as a two length vector indicating the domain of the function.
#' See Details on the specification of arguments section below more details.
#' @param mean_diff_fnm The name of the function that output of the difference of the mean between the
#' two groups at any given time. It must be supplied as character, so that \code{match.fun(mean_diff_fnm)}
#' returns a valid function, that takes a vector input, and returns a vector of the same length of the input.
#' @param cov.type  The type of the covariance structure of the data, must be either of 'ST' (stationary) or
#'                  'NS' (non-stationary). This argument along with the \code{cov.par} argument must be
#'                   specified compatibly to ensure that the function does not return an error. See the details
#'                   of \code{cov.par} argument.
#' @param cov.par The covariance structure of the latent response trajectory.
#'                If \code{cov.type == 'ST'} then, \code{cov.par}
#'                must be specified a named list of two elements, \code{var} and \code{cor},
#'                where \code{var} is the common variance of the observations, which must be a
#'                positive number; and \code{cor} specifies the correlation structure between
#'                the observations. \code{cov.par$cor} must be specified in the form of the
#'                [nlme::corClasses] specified in R package \pkg{nlme}.
#'                Check the package documentation for more details for each of the correlation classes.
#'                The \code{cov.par$cor} must be a \code{corStruct} class so it can be
#'                passed onto the [nlme::corMatrix()] to extract the subject-specific covariance matrix.
#'                If `cov.type='NS'` then, \code{cov.par}
#'                must be a named list of two elements, \code{cov.obj} and \code{eigen.comp},
#'                where only one of the \code{cov.par$cov.obj} or \code{cov.par$eigen.comp}
#'                must be non-null. This is to specify that the covariance structure of the
#'                latent trajectory can be either provided in the form of covariance function or
#'                in the form of eigenfunction and eigenvalues (Spectral decomposition).
#'                If the \code{cov.par$cov.obj} is specified, then it must be a bivariate function,
#'                with two arguments. Alternatively, if the true eigenfunctions  are known,
#'                then the user can specify that by specifying \code{cov.par$eigen.comp}.
#'                In this case, the \code{cov.par$eigen.comp} must be a named list with two elements,
#'                \code{eig.obj} and \code{eig.val}, where \code{cov.par$eigen.comp$eig.val}
#'                must be positive vector and \code{cov.par$eigen.comp$eig.obj}
#'                must be a vectorized function so that its evaluation at a vector of time points
#'                returns a matrix of dimension r by \code{length(cov.par$eigen.comp$eig.val)},
#'                with r being the length of time points.
#' @param sigma2.e Measurement error variance, should be set as zero or a very small number
#'                 if the measurement error is not significant.
#' @param missing_type The type of missing in the number of observations of the subjects. Can be one of
#'                     \code{'nomiss'} for no missing observations
#'                     or \code{'constant'} for constant
#'                     missing percentage at every time point. The current version of package only supports
#'                     \code{missing_type = 'constant'}.
#' @param missing_percent The percentage of missing at each observation points for each subject.
#' If \code{nobs_per_subj} is supplied as vector, then \code{missing_type}
#'                      is forced to set as `'nomiss'` and \code{missing_percent = 0}, because
#'                     the \code{missing_type = 'constant'} has no meaning if the number of observations are
#'                     varying between the subject at the first, typically considered in
#'                     the case of sparse random functional design.
#' @param eval_SS The sample size based on which the eigencomponents will be estimated from data.
#' To compute the theoretical power of the test we must make sure that we use a large enough sample size
#' to generate the data such that the estimated eigenfunctions are very close to the true eigenfunctions
#' and that the sampling design will not have much effect on the loss of precision. Default value 5000.
#' @param alloc.ratio The allocation ratio of samples in the each group. Note that the eigenfunctions
#'  will still be estimated based on the total sample_size, however, the variance
#'  of the `shrinkage` scores (which is required to compute the power function) will be
#' estimated based on the allocation of the samples in each group. Must be given as vector of
#' length 2. Default value is set at \code{c(1,1)}, indicating equal sample size.
#' @param fpca_method The method by which the FPCA is computed. Must be one of
#' `'fpca.sc'` and `'face'`. If \code{fpca_method == 'fpca.sc'} then the eigencomponents
#' are estimated using the function [refund::fpca.sc()]. However, since the [refund::fpca.sc()]
#' function fails to estimate the correct `shrinkage` scores, and throws \code{NA} values
#' when the measurement errors is estimated to be zero, we wrote out a similar function
#' where we corrected those error in current version of [refund::fpca.sc()]. Check out
#' the [fPASS::fpca.sc()] function for details. If \code{fpca_method == 'face'}, then
#' the eigencomponents are estimated using [face::face.sparse()] function.
#' @param data.driven.scores Indicates whether the scores are estimated from the full data, WITHOUT
#'        assuming the mean function is unknown, rather the mean function is estimated using
#'        [mgcv::gam()] function.
#' @param mean_diff_add_args Additional arguments to be passed to group difference
#'                           function specified in the argument \code{mean_diff_fnm}.
#' @param fpca_optns Additional options to be passed onto either of [fPASS::fpca.sc()]
#'                  or [face::face.sparse()] function in order
#'                   to estimate the eigencomponents. It must be a named list with elements
#'                   to be passed onto the respective function, depending on the \code{fpca_method}.
#' @return A list with the elements listed below.
#' \enumerate{
#' \item \code{mean_diff_vec} - The evaluation of the mean function at the working grid.
#' \item \code{est_eigenfun} - The evaluation of the estimated eigenfunctions at the working grid.
#' \item \code{est_eigval} - Estimated eigen values.
#' \item \code{working.grid} - The grid points at which \code{mean_diff_vec} and
#' \code{est_eigenfun} are evaluated.
#' \item \code{weights} - The Gaussian quadrature weights obtained from [gss::gauss.quad()]
#' to compute the projection \eqn{\int [\mu_1(t) - \mu_2(t)]\phi_k(t) \,dt},
#' for each \eqn{k=1,\dots, K}.
#' \item \code{fpcCall} - The exact call of either of the [fPASS::fpca.sc()] or [face::face.sparse()]
#' used to compute the eigencomponents.
#' \item \code{scores_var1} - Estimated covariance of the `shrinkage` scores for the treatment group.
#' \item \code{scores_var2} - Estimated covariance of the `shrinkage` scores for the placebo group.
#' \item \code{pooled_var}  - Pooled covariance of the scores combining both the groups. This is required
#' if the user wants to compute the power of Hotelling T statistic under equal variance assumption.
#' }
#' If \code{data.driven.scores ==  TRUE} additional components are returned
#' \enumerate{
#' \item \code{scores_1} - Estimated `shrinkage` scores for all the subjects in treatment group.
#' \item \code{scores_2} - Estimated `shrinkage` scores for all the subjects in placebo group.
#' }
#' The output of this function is designed such a way the
#' user can directly input the output obtained from this function into the arguments of
#' [fPASS::Power_Proj_Test_ufDA()] function to obtain the power and the sample size right away. The function
#' [fPASS::PASS_Proj_Test_ufDA] does the same, it is essentially a wrapper of[fPASS::Extract_Eigencomp_fDA()]
#' and [fPASS::Power_Proj_Test_ufDA()] together.
#' @seealso See [fPASS::Power_Proj_Test_ufDA()], [refund::fpca.sc()] and [face::face.sparse()].
#' @export Extract_Eigencomp_fDA
#' @references Wang, Qiyao (2021)
#' \emph{Two-sample inference for sparse functional data,  Electronic Journal of Statistics,
#' Vol. 15, 1395-1423} \cr
#' \doi{https://doi.org/10.1214/21-EJS1802}.
#' @examples
#' \dontrun{
#' set.seed(12345)
#' mean.diff <- function(t) {t};
#' obs.design <- list("design" = "longitudinal",
#' "visit.schedule" = seq(0.1, 0.9, length.out=7),
#' "visit.window" = 0.05)
#' cor.str <- nlme::corExp(1, form = ~ time | Subject);
#' sigma2 <- 1; sigma2.e <- 0.25; nobs_per_subj <- 8;
#' missing_type <- "constant"; missing_percent <- 0.01;
#' eigencomp  <- Extract_Eigencomp_fDA(obs.design = obs.design,
#' mean_diff_fnm = "mean.diff", cov.type = "ST",
#' cov.par = list("var" = sigma2, "cor" = cor.str),
#' sigma2.e = sigma2.e, nobs_per_subj = nobs_per_subj,
#' missing_type = missing_type,
#' missing_percent = missing_percent, eval_SS = 5000,
#' alloc.ratio = c(1,1),
#' fpca_method = "fpca.sc", data.driven.scores = FALSE,
#' mean_diff_add_args = list(), fpca_optns = list(pve = 0.95))
#' }
#'\dontrun{
#' alloc.ratio  <- c(1,1)
#' mean.diff    <- function(t) {1 * (t^3)};
#' eig.fun <- function(t, k) { if (k==1) {
#' ef <- sqrt(2)*sin(2*pi*t)
#' } else if (k==2) {ef <- sqrt(2)*cos(2*pi*t)}
#' return(ef)}
#' eig.fun.vec  <- function(t){cbind(eig.fun(t, 1),eig.fun(t, 2))}
#' eigen.comp   <- list("eig.val" = c(1, 0.5), "eig.obj" = eig.fun.vec)
#' obs.design   <- list(design = "functional", fun.domain = c(0,1))
#' cov.par      <- list("cov.obj" = NULL, "eigen.comp" = eigen.comp)
#' sigma2.e     <- 0.001; nobs_per_subj <- 4:7;
#' missing_type <- "nomiss"; missing_percent <- 0;
#' fpca_method  <- "fpca.sc"
#' eigencomp  <- Extract_Eigencomp_fDA(obs.design = obs.design,
#'  mean_diff_fnm = "mean.diff",
#'  cov.type = "NS", cov.par = cov.par,
#'  sigma2.e = sigma2.e, nobs_per_subj = nobs_per_subj,
#'  missing_type = missing_type,
#'  missing_percent = missing_percent, eval_SS = 5000,
#'  alloc.ratio = alloc.ratio,
#'  fpca_method = "fpca.sc", data.driven.scores = FALSE,
#'  mean_diff_add_args = list(), fpca_optns = list(pve = 0.95))
#' }
Extract_Eigencomp_fDA  <- function(nobs_per_subj, obs.design, mean_diff_fnm,
                                   cov.type = c("ST", "NS"), cov.par, sigma2.e,
                                   missing_type = c("nomiss", "constant"),
                                   missing_percent = 0,
                                   eval_SS = 5000, alloc.ratio = c(1,1),
                                   fpca_method = c("fpca.sc", "face"),
                                   data.driven.scores = FALSE,
                                   mean_diff_add_args=list(),
                                   fpca_optns = list()){

  #*********************************************************************************************%
  #*                                  Missing argument checking                                 %
  #*********************************************************************************************%
  # Functional argument checking : default value specification
  call.f       <- rlang::call_match(defaults = FALSE)
  call.f.all   <- rlang::call_match(defaults = TRUE)
  default_args <- setdiff(rlang::call_args_names(call.f.all), rlang::call_args_names(call.f))
  if (length(default_args) > 0 ){
    rlang::inform(paste0("User did not specify values to the following arguments: ",
                         paste0(default_args, collapse = ","), " ---- default values are used"))
  }
  # Non-null argument checking
  non_null_arg_vals <- !sapply(rlang::call_args(call.f), is.null)
  testthat::expect_true(all(non_null_arg_vals),
                        info = paste0("NULL values are specified for the arguments: ",
                                      paste0(rlang::call_args_names(call.f)[!non_null_arg_vals],
                                             collapse = ", "), "\n") )
  #*********************************************************************************************%
  #*                         Compatibility checking of arguments                                %
  #*********************************************************************************************%

  testthat::expect_true(rlang::is_list(obs.design))
  testthat::expect_true("design" %in% names(obs.design), info = "obs.design must be a named list with name 'design'")
  design <- match.arg(obs.design$design, c("longitudinal", "functional"))
  if (design == "functional") {
    testthat::expect_true(rlang::is_integerish(nobs_per_subj, finite=TRUE),
                          info = "nobs_per_subj must be a positive integer scalar/vector (for varying number of obs)")
    testthat::expect_named(obs.design[!names(obs.design) %in% "design"], "fun.domain", ignore.order = TRUE,
                           info = "Names of obs.design must be 'design' and 'fun.domain'.")
    # argument checking: interval
    interval <- obs.design[["fun.domain"]]
    testthat::expect_true(rlang::is_double(interval, n=2, finite=TRUE) & !is.unsorted(interval),
                          info = "obs.design$interval of the observation points must be two-length numeric vector")
  } else{
    testthat::expect_true(rlang::is_integerish(nobs_per_subj, n=1, finite=TRUE),
                          info = "nobs_per_subj must be a scalar integer for longitudinal design.")
    testthat::expect_named(obs.design[!names(obs.design) %in% "design"], c("visit.schedule", "visit.window"), ignore.order = TRUE,
                           info = "Names of obs.design must be 'design' and 'visit.schedule' and 'visit.window'.")
    visit.sched <- obs.design[["visit.schedule"]]; visit.sched.range <-  obs.design[["visit.window"]];
    testthat::expect_true(rlang::is_double(visit.sched, n=nobs_per_subj-1, finite=TRUE) & all(visit.sched > 0) &
                            !is.unsorted(visit.sched),
                          info = "obs.design$visit.schedule must be a strictly increasing sequence of positive numbers
                                  with same as the length of nobs_per_subj minus 1.")
    testthat::expect_true(rlang::is_double(visit.sched.range, n=1, finite=TRUE) |
                          rlang::is_double(visit.sched.range, n=length(visit.sched), finite=TRUE),
                          info = "obs.design$visit.window must be a positive scalar or a vector of length of
                                  visit.sched minus 1 (excluding for the baseline)")
    if (length(visit.sched) > 1)  visit.sched.range <-  rep(visit.sched.range, length(visit.sched))
    visit.window        <- purrr::map2(visit.sched - visit.sched.range, visit.sched + visit.sched.range, ~ c(.x, .y))
    testthat::expect_true(!is.unsorted(unlist(visit.window), strictly = TRUE),
                          info = "The range of visit schedule parameters specified in
                          visit.sched.range parameter is specified in such a way that
                          there is an overlap in the previous visit time.")
    visit.window.scaled <- lapply(visit.window, function(i) i/max(unlist(visit.window)))
    interval            <- c(0,1)
  }

  ngrid          <- 201
  gauss.quad.pts <- gss::gauss.quad(ngrid,interval) # evaluation points
  working.grid   <- gauss.quad.pts$pt

  #************************************************************%
  #*       Argument checking of mean function                  %
  #************************************************************%
  # argument checking: mean_diff_fnm
  mean_diff      <- match.fun(mean_diff_fnm)
  testthat::expect_true(rlang::is_function(mean_diff),
                        info = paste0("There does not exist any function with name ", mean_diff_fnm))
  # argument checking: mean_fun_args
  mean_diff_args <- rlang::fn_fmls_names(mean_diff)
  testthat::expect_true(rlang::is_list(mean_diff_add_args, n=length(mean_diff_args)-1),
                        info="mean_diff_add_args must be supplied as a named list
                              containing all the arguments of a function")
  if (length(mean_diff_args) > 1) {
    testthat::expect_setequal(names(mean_diff_add_args), mean_diff_args[-1])
  }
  mean_diff_vec  <- eval(rlang::call2("mean_diff", working.grid, !!!mean_diff_add_args))
  testthat::expect_true(rlang::is_double(mean_diff_vec, n=length(working.grid),finite=TRUE),
                        info = paste0(mean_diff_fnm, "function does not produce results as expected."))

  #************************************************************%
  #*    Argument checking of all covariance parameters         %
  #************************************************************%
  is.eig.given   <- FALSE;
  # argument checking: cov.type
  cov.type       <- match.arg(cov.type)
  # argument checking: cov.par
  testthat::expect_true(is.list(cov.par), info = "cov.par must be a list")
  if (cov.type == "ST"){
    # argument checking: cov.par
    testthat::expect_named(cov.par, c("var", "cor"), ignore.order = TRUE)
    sigma2       <- cov.par[["var"]]
    cor.str      <- cov.par[["cor"]]
    testthat::expect_true(rlang::is_double(sigma2, n=1, finite=TRUE) & (sigma2 > 0),
                          info = "variance sigma2 must be a positive number");
    testthat::expect_s3_class(cor.str, "corStruct")
    if (design == "functional"){
       testthat::expect_false(class(cor.str)[1] %in% c("corAR1", "corARMA", "corSymm"),
                             info = "cov.par[['cor']] can not be any one of c('corAR1', 'corARMA', 'corSymm')")
    }
  } else{
    # argument checking: cov.par
    testthat::expect_named(cov.par, c("cov.obj", "eigen.comp"), ignore.order = TRUE)
    cov.obj      <- cov.par[["cov.obj"]]
    eig.comp     <- cov.par[["eigen.comp"]]
    testthat::expect_true(xor(is.null(cov.obj), is.null(eig.comp)),
                          info = "Only one of the covariance parameters
                                   among covarinace function (or matrix) or
                                   eigencomponent (eigenvalues and eigen functions)
                                   must be specified")
    if (!is.null(cov.obj)){
      message("Covariance function has been provided")
      testthat::expect_true(is.function(cov.obj))
      cov.obj      <- match.fun(cov.obj)
      cov.mat      <- outer(working.grid, working.grid, cov.obj)
      testthat::expect_true(rlang::is_double(cov.mat, n=(length(working.grid))^2, finite=TRUE) &
                              isSymmetric(cov.mat),
                            info = "cov.obj function does not produce a proper covariance matrix")
    } else{
      message("Eigen function has been provided")
      testthat::expect_true(is.list(eig.comp), info = "eig.comp must be a list")
      testthat::expect_named(eig.comp, c("eig.obj", "eig.val"), ignore.order = TRUE)
      eig.val      <- eig.comp[["eig.val"]]
      testthat::expect_true(rlang::is_double(eig.val, finite=TRUE) & all(eig.val > 0),
                            info = "eigenvalues must be positive with length same as
                                    the number of eigenfunctions");
      eig.obj        <- eig.comp[["eig.obj"]]
      eig.obj        <- match.fun(eig.obj)
      testthat::expect_true(is.function(eig.obj), info="eig.obj is not a function.")
      eig.mat        <- eig.obj(working.grid)
      testthat::expect_true(rlang::is_double(eig.mat, finite = TRUE) &
                            is.matrix(eig.mat) &
                            all(dim(eig.mat) == c(length(working.grid), length(eig.val))),
                            info = "eig.obj function is not propertly vectorized to produce a matrix of
                                    eigenfunctions where the number of rows are the length of evaluation points
                                    and the columns representing the eigenfunctions.")
      is.eig.given <- TRUE
    }
  }
  # argument checking: data.driven.scores
  testthat::expect_true(rlang::is_scalar_logical(data.driven.scores));
  # argument checking: sigma2.e
  testthat::expect_true(rlang::is_double(sigma2.e, n=1, finite=TRUE) & (sigma2.e > 0),
                        info = "variance sigma2.e must be a positive number");
  #************************************************************%
  #*    Argument checking of additional parameters             %
  #*    involved in the data generation.                       %
  #************************************************************%
  # argument checking: eval_SS
  testthat::expect_true(rlang::is_integerish(eval_SS, n=1, finite = TRUE) & (eval_SS > 10),
                        info = "eval_SS must be a positive integer with value greater than 10")
  # argument checking: alloc.ratio
  testthat::expect_true(rlang::is_double(alloc.ratio, n=2, finite=TRUE) & all(alloc.ratio > 0) ,
                        info = "alloc.ratio must be positive numeric vector of length 2");
  # argument checking: missing_type
  missing_type   <- match.arg(missing_type)
  # argument checking: missing_percent
  testthat::expect_true(rlang::is_double(missing_percent, n=1, finite=TRUE) &
                        (missing_percent >= 0) & (missing_percent <= 0.8),
                        info = "missing_percent must be number between 0 and 0.8");
  # argument checking: fpca_method
  fpca_method    <- match.arg(fpca_method)
  # argument checking: fpca_options
  if (fpca_method == "fpca.sc"){
    fpca_mandatory_args_names <- list("ydata", "Y.pred", "center")
    fpca_args_all <- rlang::fn_fmls_names(fpca.sc)
  } else {
    fpca_mandatory_args_names <- c("data", "newdata", "center", "argvals.new", "calculate.scores")
    fpca_args_all  <- rlang::fn_fmls_names(face::face.sparse)
  }
  testthat::expect_equal(rlang::is_list(fpca_optns), TRUE, info = "fpca_optns must be a list");
  testthat::expect_true(dplyr::setequal(intersect(names(fpca_optns),  fpca_args_all), names(fpca_optns)),
                        info = "Unnecessary arguments to provided to pass to the fpca function")

  #*********************************************************************************************%
  #*                                  Main body of the function                                 %
  #*********************************************************************************************%
  message("Estimation of eigenfunction starts")
  if (length(nobs_per_subj) > 1){
    message("Number of observations per subject varies,
                internally setting missing_type must be `nomiss` and
                 missing percent to zero")
    missing_type    <- "nomiss"
    missing_percent <- 0
  }
  n.big             <- eval_SS # Edited by SK on Jan 17
  if (length(nobs_per_subj) > 1){
    cat("nobs_per_subj is a vector: number of measurements
         for the subjects will be taken randomly between",
        range(nobs_per_subj)[1], "and", range(nobs_per_subj)[2], "\n")
    m.orig         <- sample(nobs_per_subj, n.big, replace = TRUE)
  } else {
    m.orig         <- rep(nobs_per_subj, n.big)
  }
  visit.orig       <- lapply(m.orig, seq)
  if (design == "functional"){
    tind.orig      <- lapply(1:n.big, function(i) sort(sample(1:ngrid, m.orig[i], replace = FALSE)) )
    tvals.orig     <- lapply(1:n.big, function(i) working.grid[tind.orig[[i]]] )
  } else{
    long.grid      <- lapply(visit.window.scaled, function(il) {seq(il[1], il[2],
                                                                length.out=floor(ngrid/length(visit.window.scaled)))})
    tind.orig      <- lapply(seq_along(2:nobs_per_subj), function(i) sample(seq_along(long.grid[[i]]),
                                                                            n.big, replace = TRUE))
    tvals.orig     <- lapply(seq_along(2:nobs_per_subj), function(i) long.grid[[i]][tind.orig[[i]]])
    tind.orig      <- split(do.call(cbind, tind.orig), f=1:n.big)
    tvals.orig     <- lapply(split(do.call(cbind, tvals.orig), f=1:n.big), function(.x) c(0, .x))
  }
  group1_size      <- round(n.big * alloc.ratio[1]/sum(alloc.ratio)) # Added by SK on Jan 17
  group2_size      <- n.big - group1_size # Added by SK on Jan 17
  g                <- sample(c(rep(TRUE, group1_size), rep(FALSE, group2_size))) # Added by SK on Jan 17
  if ((missing_type == "constant") & (missing_percent > 0)){
    miss_pattern   <- sapply(1:nobs_per_subj, function(obs) { sample(x=c(TRUE, FALSE),
                                                                   size=n.big,
                                                                   prob = c(1-missing_percent, missing_percent),
                                                                   replace = TRUE) } )
    tind           <- lapply(1:n.big, function(i) tind.orig[[i]][miss_pattern[i,]])
    tvals          <- lapply(1:n.big, function(i) tvals.orig[[i]][miss_pattern[i,]])
    visit          <- lapply(1:n.big, function(i) visit.orig[[i]][miss_pattern[i,]])
    m              <- sapply(tvals, length)
  } else{
    tind           <- tind.orig
    tvals          <- tvals.orig
    m              <- m.orig
    visit          <- visit.orig
  }
  ids_with_1obs    <- m > 1
  tind             <- tind[ids_with_1obs]
  tvals            <- tvals[ids_with_1obs]
  visit            <- visit[ids_with_1obs]
  m                <- m[ids_with_1obs]
  n.new            <- sum(ids_with_1obs)
  g                <- g[ids_with_1obs]
  if (cov.type == "ST"){
    covar_form     <- stringr::str_split_1(deparse(nlme::getCovariateFormula(cor.str)), "~")[-1]
    group_form     <- stringr::str_split_1(deparse(nlme::getGroupsFormula(cor.str)), "~")[-1]
    covar_vars     <- stringr::str_split_1(covar_form, " \\+ ")
    groups_vars    <- stringr::str_split_1(group_form, " \\+ ")
    cat("grouping variable is", groups_vars, "\n")
    cat("covariate variable is", covar_vars, "\n")
    testthat::expect_true((length(groups_vars)  == 1) & (length(covar_vars)  == 1) &
                          all(c(covar_vars, groups_vars) !=  ".id.") & (covar_vars != groups_vars),
                          info = "more than one grouping vars or one covariate vars
                          present in the cor.str formula. The covariate var must be
                          indicating the visit times and the groups must be subjects, with
                          both the names anything other than .id.")
    if (class(cor.str)[1] %in% c("corAR1", "corARMA", "corSymm")){
      sim.dat      <- data.frame(".id." = rep(1:n.new, m)) %>%
                      mutate(!!covar_vars := unlist(visit)) %>% mutate(!!groups_vars := factor(.id.))
    } else{
      sim.dat      <- data.frame(".id." = rep(1:n.new, m)) %>%
                      mutate(!!covar_vars := unlist(tvals)) %>% mutate(!!groups_vars := factor(.id.))
    }
    cs1Exp         <- cor.str
    cs1Exp         <- nlme::Initialize(cs1Exp, sim.dat)
    cor.mat        <- nlme::corMatrix(cs1Exp)
    Xlist          <- lapply(1:n.new, function(i) MASS::mvrnorm(n=1, mu=rep(0, m[i]),
                                                              Sigma = sigma2*cor.mat[[i]]) )
  } else{
    if (!is.eig.given){
      Xlist        <- lapply(1:n.new, function(i) MASS::mvrnorm(n=1, mu=rep(0, m[i]),
                                                              Sigma = cov.obj[tvals[[i]], tvals[[i]] ] ) )
    } else{
      Xis          <- sapply(eig.val, function(lam) rnorm(n.new, mean = 0, sd = sqrt(lam)))
      Xlist        <- lapply(1:n.new, function(i) as.vector(eig.obj(tvals[[i]]) %*% Xis[i, ]) )
    }
  }
  Yc_list          <- lapply(1:n.new, function(i) Xlist[[i]] +
                             rnorm(m[i], mean = 0, sd = sqrt(sigma2.e)))
  mean_diff_list   <- lapply(1:n.new, function(i) eval(rlang::call2("mean_diff", tvals[[i]], !!!mean_diff_add_args)) )
  Ylist            <- lapply(1:n.new, function(i) mean_diff_list[[i]]*g[i] + Yc_list[[i]])
  if (data.driven.scores){
    gamDat         <- data.frame("subj"=rep(1:n.new, m), "y" = unlist(Ylist),
                                 "argvals" = unlist(tvals), "Group" = rep(g, m)) %>%
                      mutate(trt.ind = as.numeric(g))
    fit.m          <- mgcv::gam(y ~ s(argvals, k=12) + s(argvals, k=12, by=trt.ind), data=gamDat)
    y_mean         <- gamDat$y - as.vector(mgcv::predict.gam(fit.m, newdata=gamDat %>% mutate(trt.ind = 0)))
    fuldat         <- gamDat %>% dplyr::select(subj, argvals) %>% dplyr::mutate(y = y_mean)
    fuldat_c       <- gamDat %>% dplyr::select(subj, argvals) %>% dplyr::mutate(y=fit.m$residuals)
  } else{
    fuldat         <- data.frame("subj" = rep(1:n.new, m), "argvals" = unlist(tvals), "y" = unlist(Ylist))
    fuldat_c       <- data.frame("subj" = rep(1:n.new, m), "argvals" = unlist(tvals), "y" = unlist(Yc_list))
  }
  message("Eigenfunction is estimated based on a sample of size = ", length(Ylist), "\n") # Added by SK on Jan 17.
  if (fpca_method == "face"){
    fpca_mandatory_args <- list("data"= quote(fuldat_c), "newdata" = quote(fuldat),
                                "center"= FALSE, "argvals.new" = quote(working.grid),
                                "calculate.scores" = TRUE)
    fpca_args      <- c(fpca_mandatory_args[setdiff(names(fpca_mandatory_args), names(fpca_optns))], fpca_optns)
    fpcCall        <- rlang::call2("face.sparse", !!!fpca_args,  .ns = 'face')
    fpcObj         <- eval(fpcCall)
    est_eigenfun   <- fpcObj$eigenfunctions
    est_scores     <- fpcObj$rand_eff$scores
    est_eigval     <- fpcObj$eigenvalues
  } else{
    # renaming the dataset for for fpca.sc function
    fuldat_c       <- fuldat_c %>% dplyr::rename(.id = subj, .index = argvals, .value = y)
    fuldat         <- fuldat %>% dplyr::rename(.id = subj, .index = argvals, .value = y)
    fpca_mandatory_args <- list("ydata" = quote(fuldat_c), "Y.pred" = quote(irreg2mat(fuldat)), "center" = FALSE)
    fpca_args      <- c(fpca_mandatory_args[setdiff(names(fpca_mandatory_args), names(fpca_optns))], fpca_optns)
    fpcCall        <- rlang::call2("fpca.sc", !!!fpca_args)
    fpcObj         <- eval(fpcCall)
    est_eigenfun   <- matrix(NA, nrow = length(working.grid), ncol = ncol(fpcObj$efunctions))
    for (col in 1:ncol(fpcObj$efunctions)){
      effit        <- mgcv::gam(ef ~ s(argvals), data=data.frame("ef" = fpcObj$efunctions[, col], "argvals" = fpcObj$argvals))
      est_eigenfun[, col] <- as.vector(mgcv::predict.gam(effit, newdata = data.frame("argvals" = working.grid)) )
    }
    est_scores     <- fpcObj$scores
    est_eigval     <- fpcObj$evalues
  }
  weights          <- gauss.quad.pts$wt
  ret.objects      <- c("mean_diff_vec", "est_eigenfun", "est_eigval", "working.grid", "weights", "fpcCall")
  if(!is.null(est_scores)){
    scores_1       <- est_scores[g,  , drop=FALSE]
    scores_2       <- est_scores[!g, , drop=FALSE]
    score_var1     <- cov(scores_1)
    score_var2     <- cov(scores_2)
    pooled_var     <- cov(est_scores)
    ret.objects    <- c(ret.objects, c("score_var1", "score_var2", "pooled_var"))
    if (data.driven.scores){
      ret.objects  <- c(ret.objects, c("scores_1", "scores_2"))
    }
  } else{
    testthat::expect_false(data.driven.scores, info="You can't expect to get scores if
                           scores are not estimated from FPCA. Set data.driven.scores = FALSE instead.")
  }
  ret.val          <- lapply(ret.objects, function(obj) get(obj))
  names(ret.val)   <- ret.objects

  ret.val
}


