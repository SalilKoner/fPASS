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
#' Samp.Size <- 480; alpha.fix <- 0.05; delta <- 0.78;
#' mean.diff <- function(t) {delta*(t^3)};
#' cor.str <- nlme::corCompSymm(0.5, form = ~ time | Subject);
#' sigma2 <- 1.5; sigma2.e <- 0.25; nobs_per_subj <- 4:7;
#' missing_type <- NULL; missing_percent <- 0;
#' fpca_optns  <- list()
#' \dontrun{
#' est.power   <- TS_ProjTest_sUfDA(Samp.Size, alpha.fix,
#'                mean.diff, sigma2, cor.str, sigma2.e, nobs_per_subj,
#'                missing_type, missing_percent, fpca_optns)
#' }
#' \dontrun{
#' cor.str     <- nlme::corAR1(0.5, form = ~ time | Subject);
#' est.power   <- TS_ProjTest_sUfDA(Samp.Size, alpha.fix,
#'                mean.diff, sigma2, cor.str, sigma2.e, nobs_per_subj,
#'                missing_type, missing_percent, fpca_optns)
#' }
#' @export
#' @references Wang, Qiyao (2021)
#' \emph{Two-sample inference for sparse functional data,  Electronic Journal of Statistics,
#' Vol. 15, 1395-1423} \cr
#' \doi{https://doi.org/10.1214/21-EJS1802}.
#'
TS_ProjTest_sUfDA <- function(Samp.Size, mean.diff, sigma2, cor.str, sigma2.e,
                       nobs_per_subj, missing_type = NULL, missing_percent = NULL,
                       fpca_optns = list("PVE" = 0.95)){

  if (!is.null(fpca_optns)){
    stopifnot("Argument fpca_options must be a list" = (class(fpca_optns) == "list") )
    stopifnot("Elements of fpca_options is not valid" = all(!is.na(pmatch(names(fpca_optns),
                                                                          c("PVE")))) )
  }

  stopifnot("Covariance structure must be specified" = !is.null(cor.str),
            "Mean function not provided" = (!is.null(mean.diff) && class(mean.diff) == "function"),
            "Number of measurements per subjects must be a positive integer" =
              all(nobs_per_subj == floor(nobs_per_subj)) )

  if (length(nobs_per_subj) > 1){
    stopifnot("Number of observations per subject varies,
              so missing percent must be null" =
                is.null(missing_percent) || (missing_percent == 0))
  }

  if (!is.null(missing_type)){
    cat("Working with missing type =", missing_type, "\n")
    stopifnot("missing percent must be greater than 0 and less than 90" =
                ((missing_percent > 0) & (missing_percent < 90)) )
  }

  sigma2.err     <- ifelse(is.null(sigma2.e), 0, sigma2.e)
  n.big          <- Samp.Size
  N              <- 101
  tgrid          <- seq(0,1,length.out = N)

  if (length(nobs_per_subj) > 1){
    cat("nobs_per_subj is a vector: number of measurements
        for the subjects will be taken randomly between",
        range(nobs_per_subj)[1], "and", range(nobs_per_subj)[2], "\n")
    m.orig       <- sample(nobs_per_subj, Samp.Size, replace = TRUE)
  } else {
    m.orig       <- rep(nobs_per_subj, Samp.Size)
  }

  tind           <- lapply(1:Samp.Size, function(i) sort(sample(1:N, m.orig[i], replace = FALSE)) )
  tvals.orig     <- lapply(1:Samp.Size, function(i) tgrid[tind[[i]]] )
  g              <- sample(1:2, Samp.Size, replace = TRUE)

  if (!is.null(missing_percent) & (missing_percent > 0)){

    miss_pattern <- sapply(1:nobs_per_subj, function(obs)
                           sample(x=c(TRUE, FALSE), size=Samp.Size,
                                  prob = c(1-missing_percent, missing_percent),
                                  replace = TRUE))
    tvals        <- lapply(1:Samp.Size, function(i) tvals.orig[[i]][miss_pattern[i,]])
    m            <- sapply(tvals, length)

  } else{
    tvals        <- tvals.orig
    m            <- m.orig
  }

  ids_with_1obs  <- m > 1
  tvals          <- tvals[ids_with_1obs]
  m              <- m[ids_with_1obs]
  n.new          <- sum(ids_with_1obs)
  g              <- g[ids_with_1obs]

  sim.dat        <- data.frame("id" = rep(1:n.new, m),
                               "time" = unlist(tvals)) %>%
                    mutate(Subject = factor(id))
  cs1Exp         <- cor.str
  cs1Exp         <- nlme::Initialize(cs1Exp, sim.dat)
  cor.mat        <- nlme::corMatrix(cs1Exp)

  Xlist          <- lapply(1:n.new, function(i) MASS::mvrnorm(n=1, mu=rep(0, m[i]),
                                                              Sigma = sigma2*cor.mat[[i]]) +
                             + as.numeric(g[i]==2)*mean.diff(tvals[[i]]))

  Ylist          <- lapply(1:n.new, function(i) Xlist[[i]] +
                             rnorm(m[i], mean = 0, sd = sqrt(sigma2.err)))

  gamDat    <- data.frame("subj"=rep(1:n.new, m), "y" = unlist(Ylist),
                          "argvals" = unlist(tvals), "Group" = rep(g, m)) %>%
               mutate(trt.ind = as.numeric(Group == 2))
  fit.m     <- mgcv::gam(y ~ s(argvals, k=12) + s(argvals, k=12, by=trt.ind), data=gamDat)
  y_mean    <- gamDat$y - as.vector(predict(fit.m, newdata=gamDat %>% mutate(trt.ind = 0)))
  fpcDt.pr  <- gamDat %>% dplyr::select(subj, argvals) %>% mutate(y = y_mean)
  fpcDt     <- gamDat %>% dplyr::select(subj, argvals) %>% mutate(y=fit.m$residuals)


  faceObj   <- face::face.sparse(data=fpcDt, newdata=fpcDt.pr, center = F,
                                 knots=12, calculate.scores = TRUE,
                                 pve=ifelse(is.null(fpca_optns$PVE),
                                            0.99, fpca_optns$PVE))

  grp1_xi   <- faceObj$rand_eff$scores[g==1, , drop=FALSE]
  grp2_xi   <- faceObj$rand_eff$scores[g==2, , drop=FALSE]

  HT.test   <- Hotelling::hotelling.test(x = grp2_xi, y = grp1_xi,
                                         shrinkage = FALSE, var.equal = TRUE,
                                         perm = FALSE, B = 10000, progBar = (perm && TRUE))

  pVal_HT   <- HT.test$pval

  pVal_HT

}

