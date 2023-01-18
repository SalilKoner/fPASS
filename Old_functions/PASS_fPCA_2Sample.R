#' Power and Sample size (PASS) calculation of
#' Two-Sample (TS) Projection-based test (ProjTest) for sparsely observed univariate functional data.
#' @description This function `PASS_TS_ProjTest_sUfDA()` computes the power and sample size required to conduct
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
#' See Wang (2021) for more details of the testing procedure. We use the `FPCA()`
#' function in the R package `fdapace` to conduct the functional principal component analysis.
#' @author Salil Koner \cr Maintainer: Salil Koner
#' \email{salil.koner@@duke.edu}
#' @import Matrix
#' @import dplyr
#' @import fdapace
#' @import refund
#' @import face
#' @importFrom nlme corMatrix Initialize
#' @importFrom MASS mvrnorm
#' @importFrom stats cov rnorm uniroot
#' @param target.SS Target sample size, must be a positive interger more than 3.
#' @param target.power Target power, must be a number between 0 and 1.
#'                     Only one of target.SS and target.power should be non-null. The
#'                     function will automatically solve for sample size if target.SS is NULL, and
#'                     solve for power if target.power is NULL.
#' @param alpha.fix Level of significance of the test. Must be a number between 0 and 1, possibly less than 0.2.
#' @param ngrid The number of equidistant grid points between (0,1) at which the eigencomponents are estimated/calculated.
#'              The function internally creates a regular grid of points of length = ngrid between (0,1),
#'              i.e. `tgrid = seq(0,1,length.out = ngrid)` to conduct the test.
#' @param mean.diff The difference in the mean function between the two groups. Must be supplied as a function class.
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
#'                specified, then it must be either a bivariate function or a ngrid \times ngrid matrix,
#'                containing the covariance between two timepoints in tgrid. Alternatively,
#'                if the true eigenfunctions (and number of eigenfunctions, say K) are known,
#'                then the user can specify that by specifying eigen.comp.
#'                In this case, the eigen.comp must be a named list with two elements, "eig.obj" and "eig.val".
#'                The object "eig.obj" can either a function that returns a vector of length K for fixed timepoint
#'                or it can be a matrix of dimension ngrid by K where eig.mat(i,k) = phi_k(t_i) for i=1,2,..ngrid and k=1,..K.
#'                And the quantity "eig.val" must be a numeric vector of length K with positive elements.
#' @param sigma2.e Measurement error variance, should be left as NULL if there is no measurement error.
#' @param nobs_per_subj The number of observations per subject. Must be a positive integer greater than 3.
#' @param missing_type The type pf missing in the number of observations of the subjects.
#' Only supports \code{missing_type = "constant"} now.
#' @param missing_percent The percentage of missing at each observation points for each subject.
#' @param fpca_optns Additional options to be passed onto the [fdapace::FPCA()] function in order
#' to estimate the eigencomponent required for the sample size calculation.
#' @return If target.SS is null, it returns the optimal sample size required for specified target.power.
#'         Otherwise returns the empirical power for the sample size specified by target.SS.
#' @examples
#' set.seed(12345)
#' target.SS <- NULL; target.power <- 0.8; alpha.fix <- 0.05; mean.diff <- function(t) {t};
#' cor.str <- nlme::corExp(1, form = ~ time | Subject);
#' sigma2 <- 1; sigma2.e <- 0.25; nobs_per_subj <- 6;
#' missing_type <- "constant"; missing_percent <- 0.01;
#' fpca_optns <- list("FVEthreshold" = 0.95)
#' \dontrun{
#' Sample_size <- PASS_TS_ProjTest_sUfDA(target.SS = target.SS,
#' target.power = target.power, alpha.fix = alpha.fix,
#' ngrid = 101, mean.diff= mean.diff,
#' cov.par = list("var" = sigma2, "cor" = cor.str), sigma2.e = sigma2.e,
#' nobs_per_subj = nobs_per_subj, missing_type = missing_type,
#' missing_percent = missing_percent, fpca_optns =fpca_optns)
#' }
#'
#' @export
#' @references Wang, Qiyao (2021)
#' \emph{Two-sample inference for sparse functional data,  Electronic Journal of Statistics,
#' Vol. 15, 1395-1423} \cr
#' \doi{https://doi.org/10.1214/21-EJS1802}.
#'
PASS_TS_ProjTest_sUfDA <- function(target.SS, target.power, alpha.fix, ngrid = 101,
                       mean.diff, cov.type = c("ST", "NS"),
                       cov.par,
                       sigma2.e, nobs_per_subj,
                       missing_type = c("nomiss", "constant"),
                       missing_percent = NULL, fpca_optns = list()){

  stopifnot("Only one of target.SS or target.power should be NULL" =
            xor(!is.null(target.SS), !is.null(target.power)),
            "Mean function not provided" = (!is.null(mean.diff) && is.function(mean.diff) )
            )

  tgrid          <- seq(0,1,length.out = ngrid)
  sigma2.err     <- ifelse(is.null(sigma2.e), 0, sigma2.e)
  cov.type       <- match.arg(cov.type)
  missing_type   <- match.arg(missing_type)

  eig.compute    <- TRUE
  is.eig.given   <- FALSE

  if (cov.type == "ST"){

    stopifnot("Names of the covariance parameters must be a list with names
              var (common variance) and cor (correlation structure)" =
              is.list(cov.par) & all(!is.na(pmatch(names(cov.par), c("var", "cor"))))
    )

    sigma2       <- cov.par[["var"]]
    cor.str      <- cov.par[["cor"]]

    stopifnot("Either of the Covariance parameters are not specified" =
                (!is.null(sigma2) & !is.null(cor.str)))
    stopifnot("common variance of stationary covariance must be positive" = sigma2 > 0)

  } else{

    stopifnot("Names of the covariance parameters must be a list with names
              cov.obj (covariance object) and/or eigen.comp (eigencomponents)" =
              is.list(cov.par) & all(!is.na(pmatch(names(cov.par), c("cov.obj", "eigen.comp"))))
    )

    cov.obj      <- cov.par[["cov.obj"]]
    eig.comp     <- cov.par[["eigen.comp"]]

    stopifnot("Only one of the covariance parameters among covarinace function (or matrix) or
               eigencomponent (eigenvalues and eigen functions) must be specified" =
               xor(is.null(cov.obj), is.null(eig.comp)))

    if (!is.null(cov.obj)){

      if (is.function(cov.obj)){
        message("Covariance function has been provided")
        cov.mat      <- outer(tgrid, tgrid, cov.obj)
      } else if (is.matrix(cov.obj)) {
        message("Covariance matrix has been provided")
        stopifnot("Dimension of the covariance matrix does not match with the ngrid argument" =
                  all(dim(cov.obj) == ngrid))
        cov.mat      <- cov.obj
      } else{
        stop("The covariance object for non-stationary must be a function of matrix")
      }
    } else{
      message("Eigen function or matrix has been provided, no need to compute the eigenfunctions")
      stopifnot("Eigencomponents must be a list with two entries, eig.fun and eig.val" =
                  is.list(eig.comp) &
                  all(!is.na(pmatch(names(eig.comp), c("eig.obj", "eig.val", "eig.compute"))))
                )
      eig.obj      <- eig.comp[["eig.obj"]]
      eig.val      <- eig.comp[["eig.val"]]
      eig.compute  <- eig.comp[["eig.compute"]]

      stopifnot("The eigenfunctions must be a numeric vector with positive elements" = all(eig.val > 0))
      stopifnot("The eig.compute element must be logical" = is.logical(eig.compute))

      if (is.function(eig.obj)){
        message("Eigen function has been provided")
        eig.mat      <- t(sapply(tgrid, function(ti) eig.obj(ti)))
      } else if (is.matrix(eig.obj)) {
        message("Eigen matrix has been provided")
        stopifnot("Dimension of the covariance matrix
                  does not match with the ngrid argument" = (nrow(eig.obj) == ngrid) )
        eig.mat      <- eig.obj
      } else{
        stop("The eigen object must be a function of matrix")
      }

      if (ncol(eig.mat) != length(eig.val) ){
        stop("Number of eigen values and the eignfunctions does not match")
      }

      is.eig.given <- TRUE
    }

  }

  if (eig.compute){

    message("True eigenfunctions are unknown, it will be estimated from the data!")
    stopifnot("Number of measurements per subjects must be a positive integer" =
                all(nobs_per_subj == floor(nobs_per_subj)) )

    stopifnot("Missing percent must be a number between 0 and 100" = is.numeric(missing_percent))

    if (length(nobs_per_subj) > 1){
      stopifnot("Number of observations per subject varies,
                 so missing_type must be `nomiss` and
                 missing percent must be set to zero" =
                 (missing_type == "nomiss") & (missing_percent == 0))
    }

    if (missing_type != "nomiss"){
      cat("Working with missing type =", missing_type, "\n")
      stopifnot("missing percent must be between 0 to 90 when missing type i" =
                  ((missing_percent >= 0) & (missing_percent < 90)) )
    }

    n.big          <- 1e3
    if (length(nobs_per_subj) > 1){
      cat("nobs_per_subj is a vector: number of measurements
        for the subjects will be taken randomly between",
          range(nobs_per_subj)[1], "and", range(nobs_per_subj)[2], "\n")
      m.orig       <- sample(nobs_per_subj, n.big, replace = TRUE)
    } else {
      m.orig       <- rep(nobs_per_subj, n.big)
    }

    tind.orig      <- lapply(1:n.big, function(i) sort(sample(1:ngrid, m.orig[i], replace = FALSE)) )
    tvals.orig     <- lapply(1:n.big, function(i) tgrid[tind.orig[[i]]] )
    g              <- sample(1:2, n.big, replace = TRUE)

    if (!is.null(missing_percent) || (missing_percent > 0)){

      miss_pattern <- sapply(1:nobs_per_subj, function(obs)
        sample(x=c(TRUE, FALSE), size=n.big,
               prob = c(1-missing_percent, missing_percent),
               replace = TRUE))
      tind         <- lapply(1:n.big, function(i) tind.orig[[i]][miss_pattern[i,]])
      tvals        <- lapply(1:n.big, function(i) tvals.orig[[i]][miss_pattern[i,]])
      m            <- sapply(tvals, length)

    } else{
      tind         <- tind.orig
      tvals        <- tvals.orig
      m            <- m.orig
    }

    ids_with_1obs  <- m > 1
    tind           <- tind[ids_with_1obs]
    tvals          <- tvals[ids_with_1obs]
    m              <- m[ids_with_1obs]
    n.new          <- sum(ids_with_1obs)
    g              <- g[ids_with_1obs]

    if (cov.type == "ST"){
      sim.dat      <- data.frame("id" = rep(1:n.new, m),
                                 "time" = unlist(tvals)) %>% mutate(Subject = factor(id))
      cs1Exp       <- cor.str
      cs1Exp       <- nlme::Initialize(cs1Exp, sim.dat)
      cor.mat      <- nlme::corMatrix(cs1Exp)

      Xlist        <- lapply(1:n.new, function(i) MASS::mvrnorm(n=1, mu=rep(0, m[i]),
                                                      Sigma = sigma2*cor.mat[[i]]) )
    } else{

      if (!is.eig.given){
        Xlist      <- lapply(1:n.new, function(i) MASS::mvrnorm(n=1, mu=rep(0, m[i]),
                                                    Sigma = cov.mat[tind[[i]], tind[[i]] ]) )
      } else{
        Xis        <- sapply(eig.val, function(lam) rnorm(n.new, mean = 0, sd = sqrt(lam)))
        X.all      <- Xis %*% t(eig.mat)
        Xlist      <- lapply(1:n.new, function(i) X.all[i, tind[[i]]])
      }

    }

    Ylist        <- lapply(1:n.new, function(i) Xlist[[i]] +
                             rnorm(m[i], mean = 0, sd = sqrt(sigma2.err)))
    optns        <-  fdapace::SetOptions(Ylist, tvals, fpca_optns)
    CheckOptions(tvals, optns, length(Ylist))

    fpcObj       <- fdapace::FPCA(Ly = Ylist, Lt = tvals,
                                    optns = list(dataType = "Sparse",
                                                 userMu = list(t = tgrid, mu = rep(0, ngrid)),
                                                 nRegGrid = ngrid,
                                                 FVEthreshold = ifelse(is.null(fpca_optns$FVEthreshold),
                                                                       0.95, fpca_optns$FVEthreshold),
                                                 methodSelectK = ifelse(is.null(fpca_optns$methodSelectK),
                                                                        "FVE", fpca_optns$methodSelectK) ))

    eff.size     <- colMeans(sweep(fpcObj$phi, 1, mean.diff(tgrid), FUN = "*"))
    sigma.mat    <- cov(fpcObj$xiEst)
  } else{

    message("True eigenfunctions are known, will be directly used to obtain the projection!")
    eff.size     <- colMeans(sweep(eig.mat, 1, mean.diff(tgrid), FUN = "*"))
    sigma.mat    <- diag(eig.val)
  }

  p.body <- quote({

    k            <-length(eff.size)
    k1           <-1
    k2           <-1
    lamda        <-as.numeric(eff.size%*%solve(sigma.mat*((k1+k2)/(k1*k2)))%*%eff.size) # ncp
    cr           <-qf(1-alpha.fix,k,k1*n+k2*n-k-1) # critical region
    pwr          <-1-pf(cr,k,k1*n+k2*n-k-1,n*lamda) # power for n
    pwr

  })

  if (is.null(target.SS)){
    opt.SS <- uniroot(function(n) eval(p.body) - target.power, c(4, 1e+07),
                      tol = .Machine$double.eps^0.25, extendInt = "upX")$root
  } else if (is.null(target.power)){
    pwr.fn       <- function(n) eval(p.body)
    emp.power    <- pwr.fn(target.SS)
  }

  ret.objects    <- ifelse(is.null(target.SS), "opt.SS", "emp.power")
  ret.val        <- lapply(ret.objects, function(obj) get(obj))
  names(ret.val) <- ret.objects

  ret.val
}




