#' Twp-Sample Projection-based test for sparsely observed univariate functional data
#' @description The function `TS_ProjTest_sUfDA()` conducts
#' the Two-sample projection-based test of significant difference in mean function
#' for twp groups sparsely observed univariate functional data
#' or sparsely observed functional data under a random irregular design, under
#' the assumption of common covariance structure between the two groups.
#' @details The projection-based test assumes the two groups of the data share the identical
#' covariance structure. The projection-based test represents the sparsely observed
#' functions parsimoniously in terms of Karhunen-Loeve (KL) expansion and use the
#' functional principal component analysis scores to test the difference in the mean
#' function between two groups, by taking advantage of multivariate Hotelling-\eqn{T^2} test.
#' See Wang (2021) for more details of the testing procedure. We use the `face.sparse()`
#' function in the R package `face` to conduct the functional principal component analysis.
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
#' @param Samp.Size The sample size (combined of two groups)
#' at which the power to be calculated
#' @param mean.diff The difference in the mean function between the two groups. Must be supplied as a function class.
#' @param sigma2 The common variance of the observations. Only implementing common variance across all measurement points.
#' @param cor.str The correlation structure between the observations. It must be provided in the form specified in the
#'                available in the documentation of R package \pkg{nlme}. Check the package documentation for more details.
#'                The argument of this function is passed onto the [nlme::corMatrix()] to extract the subject-specific
#'                covariance matrix.
#' @param sigma2.e Measurement error variance, should be left as NULL if there is no measurement error.
#' @param nobs_per_subj The number of observations per subject. Must be a positive integer greater than 3.
#' @param missing_type The type pf missing in the number of observations of the subjects.
#' Only supports \code{missing_type = "constant"} now.
#' @param missing_percent The percentage of missing at each observation points for each subject.
#' @param fpca_optns Additional options to be passed onto the [face::face.sparse()] function in order
#' to estimate the eigencomponents required for the test.
#' @return Returns the p-value of the projection-based hypothesis test.
#' @examples
#' set.seed(12345)
#' Samp.Size <- 480; delta <- 0.78;
#' mean.diff <- function(t) {delta*(t^3)};
#' alloc.ration <- c(1,1); ngrid <- 101;
#' cor.str <- nlme::corCompSymm(0.5, form = ~ time | Subject);
#' sigma2 <- 1.5; sigma2.e <- 0.25; nobs_per_subj <- 4:7;
#' cov.par <- list("var" = sigma2, "cor" = cor.str);
#' missing_type <- "nomiss"; missing_percent <- 0;
#' fpca_optns  <- list("FVEthreshold" = 0.95)
#' \dontrun{
#' est.power   <- TS_ProjTest_sUfDA(Samp.Size, alloc.ratio = c(1,1),
#'                                  ngrid = ngrid, mean.diff = mean.diff,
#'                                  cov.type = "ST", cov.par = cov.par, sigma2.e = sigma2.e,
#'                                  nobs_per_subj = nobs_per_subj, missing_type = missing_type,
#'                                  missing_percent = missing_percent,
#'                                  fpca_method = fpca_method, fpca_optns =fpca_optns)
#' }
#' #' cor.str     <- nlme::corAR1(0.5, form = ~ time | Subject);
#' cov.par     <- list("var" = sigma2, "cor" = cor.str);
#' \dontrun{
#' est.power   <- TS_ProjTest_sUfDA(Samp.Size = Samp.Size,
#'                                  alloc.ratio = alloc.ratio,
#'                                  ngrid = ngrid, mean.diff = mean.diff,
#'                                  cov.type = "ST", cov.par = cov.par, sigma2.e = sigma2.e,
#'                                  nobs_per_subj = nobs_per_subj, missing_type = missing_type,
#'                                  missing_percent = missing_percent,
#'                                  fpca_method = fpca_method, fpca_optns =fpca_optns)
#' }
#' @export
#' @references Wang, Qiyao (2021)
#' \emph{Two-sample inference for sparse functional data,  Electronic Journal of Statistics,
#' Vol. 15, 1395-1423} \cr
#' \doi{https://doi.org/10.1214/21-EJS1802}.
#'
TS_ProjTest_sUfDA <- function(Samp.Size, alloc.ratio = c(1,1), ngrid = 101,
                              mean.diff, cov.type = c("ST", "NS"), cov.par,
                              sigma2.e, nobs_per_subj,
                              missing_type = c("nomiss", "constant"),
                              missing_percent = 0,
                              fpca_method = c("face", "fpca.sc"), fpca_optns = list()){

  #set.seed(12345)

  # Null condition checking
  if(!is.null(Samp.Size)){
    testthat::expect_length(Samp.Size, 1);
    testthat::expect_is(Samp.Size, "numeric", info = "Sample size must be numeric")
    testthat::expect_equal(Samp.Size == floor(Samp.Size), TRUE,
                           info = "Sample size must be a positive integer")
    testthat::expect_gte(Samp.Size, 10)
  }

  testthat::expect_equal(is.numeric(alloc.ratio), TRUE, info = "alloc.ratio must be numeric");
  testthat::expect_length(alloc.ratio, 2);
  testthat::expect_equal(all(alloc.ratio > 0), TRUE, info = "allocation numbers must be positive");

  testthat::expect_length(ngrid, 1);
  testthat::expect_is(ngrid, c("numeric", "integer"), info = "ngrid must be numeric")
  testthat::expect_equal(ngrid == floor(ngrid), TRUE, info = "ngrid must be a positive integer")
  testthat::expect_lte(ngrid, 501); testthat::expect_gte(ngrid, 51);

  testthat::expect_is(sigma2.e, "numeric", info = "sigma2.e must be numeric");
  testthat::expect_length(sigma2.e, 1);
  testthat::expect_gte(sigma2.e, 0);

  testthat::expect_equal(is.numeric(nobs_per_subj), TRUE, info = "nobs_per_subj must be numeric")
  testthat::expect_equal(all(nobs_per_subj == floor(nobs_per_subj)), TRUE,
                         info = "each element of nobs_per_subj must be a positive integer")

  testthat::expect_is(missing_percent, "numeric", info = "missing_percent must be numeric");
  testthat::expect_length(missing_percent, 1);
  testthat::expect_gte(missing_percent, 0); testthat::expect_lte(missing_percent, 80);

  testthat::expect_equal(is.function(mean.diff), TRUE, info = "mean.diff must be a function class")
  testthat::expect_equal(is.list(cov.par), TRUE, info = "cov.par must be a list")

  testthat::expect_equal(is.list(fpca_optns), TRUE, info = "fpca_optns must be a named list");
  # testthat::expect_named(fpca_optns, c("FVEthreshold", "methodSelectK"), ignore.order = TRUE);

  tgrid          <- seq(0,1,length.out = ngrid)
  cov.type       <- match.arg(cov.type)
  missing_type   <- match.arg(missing_type)
  fpca_method    <- match.arg(fpca_method)

  eig.compute    <- TRUE
  is.eig.given   <- FALSE

  if (cov.type == "ST"){

    testthat::expect_named(cov.par, c("var", "cor"), ignore.order = TRUE)
    sigma2       <- cov.par[["var"]]
    cor.str      <- cov.par[["cor"]]
    testthat::expect_type(sigma2, "double");
    testthat::expect_length(sigma2, 1);
    testthat::expect_gt(sigma2, 0);
    testthat::expect_s3_class(cor.str, "corStruct")

  } else{
    testthat::expect_named(cov.par, c("cov.obj", "eigen.comp"), ignore.order = TRUE)
    cov.obj      <- cov.par[["cov.obj"]]
    eig.comp     <- cov.par[["eigen.comp"]]

    testthat::expect_equal(xor(is.null(cov.obj), is.null(eig.comp)), TRUE,
                           info = "Only one of the covariance parameters
                                   among covarinace function (or matrix) or
                                   eigencomponent (eigenvalues and eigen functions)
                                   must be specified")

    if (!is.null(cov.obj)){
      testthat::expect_is(cov.obj, c("function", "matrix"))
      if (is.function(cov.obj)){
        message("Covariance function has been provided")
        cov.mat      <- outer(tgrid, tgrid, cov.obj)
      } else if (is.matrix(cov.obj)) {
        message("Covariance matrix has been provided")
        stopifnot("Dimension of the covariance matrix does not match with the ngrid argument" =
                    all(dim(cov.obj) == ngrid))
        cov.mat      <- cov.obj
      }
    } else{
      message("Eigen function or matrix has been provided")
      testthat::expect_equal(is.list(eig.comp), TRUE, info = "eig.comp must be a list")
      testthat::expect_named(eig.comp, c("eig.obj", "eig.val", "eig.compute"), ignore.order = TRUE)

      eig.obj      <- eig.comp[["eig.obj"]]
      eig.val      <- eig.comp[["eig.val"]]
      eig.compute  <- eig.comp[["eig.compute"]]

      testthat::expect_equal(is.numeric(eig.val), TRUE, info = "eigenvalues must be numeric vector");
      testthat::expect_equal(all(eig.val > 0), TRUE, info = "eigenvalues must be positive");
      testthat::expect_length(eig.compute, 1); testthat::expect_type(eig.compute, "logical");
      testthat::expect_is(eig.obj, c("function", "matrix"))

      if (is.function(eig.obj)){
        message("Eigen function has been provided")
        eig.mat      <- t(sapply(tgrid, function(ti) eig.obj(ti)))
      } else if (is.matrix(eig.obj)) {
        message("Eigen matrix has been provided")
        testthat::expect_equal(nrow(eig.obj), ngrid,
                               info = "number of rows of eigen matrix and
                                       ngrid do not match")
        eig.mat      <- eig.obj
      }

      testthat::expect_equal(ncol(eig.mat), length(eig.val),
                             info = "Number of eigen values and the
                                    eignfunctions do not match")
      is.eig.given <- TRUE
    }
  }

    message("True eigenfunctions are unknown, it will be estimated from the data!")
    if (length(nobs_per_subj) > 1){
      message("Number of observations per subject varies,
                internally setting missing_type must be `nomiss` and
                 missing percent to zero")
      missing_type    <- "nomiss"
      missing_percent <- 0
    }

    n.big          <- Samp.Size
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

    group1_size    <- round(n.big * alloc.ratio[1]/sum(alloc.ratio))
    group2_size    <- n.big - group1_size
    g              <- sample(c(rep(1, group1_size), rep(2, group2_size)))

    if ((missing_type == "constant") & (missing_percent > 0)){

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
                                                  Sigma = sigma2*cor.mat[[i]])
                                          + as.numeric(g[i]==2)*mean.diff(tvals[[i]]) )
    } else{

      if (!is.eig.given){
        Xlist      <- lapply(1:n.new, function(i) MASS::mvrnorm(n=1, mu=rep(0, m[i]),
                                                  Sigma = cov.mat[tind[[i]], tind[[i]] ])
                                            + as.numeric(g[i]==2)*mean.diff(tvals[[i]]) )
      } else{
        Xis        <- sapply(eig.val, function(lam) rnorm(n.new, mean = 0, sd = sqrt(lam)))
        X.all      <- Xis %*% t(eig.mat)
        Xlist      <- lapply(1:n.new, function(i) X.all[i, tind[[i]]] +
                               as.numeric(g[i]==2)*mean.diff(tvals[[i]]) )
      }

    }

  Ylist   <- lapply(1:n.new, function(i) Xlist[[i]] +
                             rnorm(m[i], mean = 0, sd = sqrt(sigma2.e)))

  gamDat    <- data.frame("subj"=rep(1:n.new, m), "y" = unlist(Ylist),
                          "argvals" = unlist(tvals), "Group" = rep(g, m)) %>%
               mutate(trt.ind = as.numeric(Group == 2))

  fit.m     <- mgcv::gam(y ~ s(argvals, k=12) + s(argvals, k=12, by=trt.ind), data=gamDat)
  y_mean    <- gamDat$y - as.vector(predict(fit.m, newdata=gamDat %>% mutate(trt.ind = 0)))
  fpcDt.pr  <- gamDat %>% dplyr::select(subj, argvals) %>% mutate(y = y_mean)
  fpcDt     <- gamDat %>% dplyr::select(subj, argvals) %>% mutate(y=fit.m$residuals)

  if (fpca_method == "face"){
    fpcObj   <- face::face.sparse(data=fpcDt, newdata=fpcDt.pr, center = F,
                                   knots=12, calculate.scores = TRUE,
                                   pve=ifelse(is.null(fpca_optns$FVEthreshold),
                                              0.95, fpca_optns$FVEthreshold))

    grp1_xi   <- fpcObj$rand_eff$scores[g==1, , drop=FALSE]
    grp2_xi   <- fpcObj$rand_eff$scores[g==2, , drop=FALSE]

  } else if (fpca_method == "fpca.sc"){
    fuldat       <- fpcDt %>% dplyr::select(subj, argvals, y) %>%
                              dplyr::rename(.id = subj, .index = argvals, .value = y)
    fuldat.pr    <- fpcDt.pr %>% dplyr::select(subj, argvals, y) %>%
                                 dplyr::rename(.id = subj, .index = argvals, .value = y)
    fpcObj       <- fpca.sc(ydata = fuldat, Y.pred = irreg2mat(fuldat.pr), nbasis = 12,
                            center = FALSE, var = TRUE,
                            pve = ifelse(is.null(fpca_optns$FVEthreshold), 0.95, fpca_optns$FVEthreshold))

    grp1_xi   <- fpcObj$scores[g==1, , drop=FALSE]
    grp2_xi   <- fpcObj$scores[g==2, , drop=FALSE]
  }



  HT.test   <- Hotelling::hotelling.test(x = grp2_xi, y = grp1_xi,
                                         shrinkage = FALSE, var.equal = TRUE,
                                         perm = FALSE, B = 10000, progBar = (perm && TRUE))

  pVal_HT   <- HT.test$pval

  pVal_HT

}

