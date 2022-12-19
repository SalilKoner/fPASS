#' Power and Sample size calculation of
#' Projection-based test for difference between mean function of
#' two groups.
#' @description This function `PASS_UfPCA()`` computes the power and sample size required to conduct
#' the projection-based test of mean function between two groups of longitudinal data
#' or sparsely observed functional data under a random irregular design, under
#' different covariance structure of the data.
#' @details The projection-based test assumes the two groups of the data share the identical
#' covariance structure. The projection-based test represents the sparsely observed
#' functions parsimoniously in terms of Kahrunen-Loeve (KL) expansion and use the
#' functional principal component analysis scores to test the difference in the mean
#' function between two groups, by taking advantage of multivariate Hotelling-$T^2$ test.
#' See the EJS paper by Qiao Wang for more details of the testing procedure. We use the `FPCA()`
#' function in the R package `fdapace` to conduct the functional principal component analysis.
#' @author Salil Koner \cr Maintainer: Salil Koner
#' \email{salil.koner@@duke.edu}
#' @import Matrix
#' @import dplyr
#' @import fdapace
#' @import refund
#' @import face
#' @import nlme
#' @import MASS
#' @param target.SS Target sample size, must be a positive interger more than 3.
#' @param target.power Target power, must be a number between 0 and 1.
#'                     Only one of target.SS and target.power should be non-null. The
#'                     function will automatically solve for sample size if target.SS is NULL, and
#'                     solve for power if target.power is NULL.
#' @param alpha.fix Level of significance of the test. Must be a number between 0 and 1, possibly less than 0.2.
#' @param mean.diff The difference in the mean function between the two groups. Must be supplied as a function class.
#' @param sigma2 The common variance of the observations. Only implementing common variance across all measurement points.
#' @param cor.str The correlation structure between the observations. It must be provided in the form specified in the
#'                available in the documentation of nlme R package. Check the package documentation for more details.
#'                The argument of this function is passed onto the [nlme::corMatrix()] to extract the subject-specific
#'                covariance matrix.
#' @param sigma2.e Measurement error variance, should be left as NULL if there is no measurement error.
#' @param nobs_per_subj The number of observations per subject. Must be a positive integer greater than 3.
#' @return If target.SS is null, it returns the optimal sample size required for specified target.power.
#'         Otherwise returns the empirical power for the sample size specified by target.SS.
#' @export
#' @examples
#' set.seed(12345)
#' target.SS <- NULL; target.power <- 0.8; alpha.fix <- 0.05; mean.diff <- function(t) {t};
#' cor.str <- corExp(1, form = ~ time | Subject); sigma2 <- 1; sigma2.e <- 0.25; nobs <- 8;
#' Sample_size <- PASS_UfPCA(target.SS, target.power, alpha.fix, mean.diff, sigma2, cor.str, sigma2.e, nobs)
#'
PASS_UfPCA <- function(target.SS, target.power, alpha.fix,
                       mean.diff, sigma2, cor.str, sigma2.e,
                       nobs_per_subj, missing_type = "constant",
                       missing_percent = NULL){

  require(tidyverse)
  require(fdapace)

  stopifnot("Only one of target.SS or target.power should be NULL" =
            xor(!is.null(target.SS), !is.null(target.power)),
            "Covariance structure must be specified" = !is.null(cor.str),
            "Mean function not provided" = (!is.null(mean.diff) && class(mean.diff) == "function"),
            "Number of measurements per subjects must be a positive integer" =
             all(nobs_per_subj == floor(nobs_per_subj)) )

  if (length(nobs_per_subj) > 1){
    stopifnot("Number of observations per subject varies, so missing percent must be null" =
              is.null(missing_percent))
  }

  if (missing_type == "constant"){
    stopifnot("missing percent must be between 0 to 90" =
                ((missing_percent >= 0) & (missing_percent < 90)) )
  }



  sigma2.err     <- ifelse(is.null(sigma2.e), 0, sigma2.e)
  n.big          <- 1e3
  N              <- 101
  tgrid          <- seq(0,1,length.out = N)

  if (length(nobs) > 1){
    cat("nobs is a vector: number of measurements for the subjects will be taken randomly between",
                            range(obs)[1], "and", range(nobs)[2], "\n")
    m            <- sample(nobs, n.big, replace = TRUE)
  } else {
    m            <- rep(nobs, n.big)
  }

  tind           <- lapply(1:n.big, function(i) sort(sample(1:N, m[i], replace = FALSE)) )
  tvals          <- lapply(1:n.big, function(i) tgrid[tind[[i]]] )
  g              <- sample(1:2, n.big, replace = TRUE)

  sim.dat        <- data.frame("id" = rep(1:n.big, m),
                               "time" = unlist(tvals)) %>%
                    mutate(Subject = factor(id))
  cs1Exp         <- cor.str
  cs1Exp         <- nlme::Initialize(cs1Exp, sim.dat)
  cor.mat        <- nlme::corMatrix(cs1Exp)

  Xlist          <- lapply(1:n.big, function(i) MASS::mvrnorm(n=1, mu=rep(0, m[i]),
                                                    Sigma = sigma2*cor.mat[[i]]) )

  Ylist          <- lapply(1:n.big, function(i) Xlist[[i]] +
                             rnorm(m[i], mean = 0, sd = sqrt(sigma2.err)))

  fpcObj         <- fdapace::FPCA(Ly = Ylist, Lt = tvals, list(dataType = "Sparse",
                          userMu = list(t = tgrid, mu = rep(0, length(tgrid))),
                                        nRegGrid = 101, FVEthreshold = 0.95,
                                        methodSelectK = "FVE"))

  eff.size       <- colMeans(sweep(fpcObj$phi, 1, mean.diff(tgrid), FUN = "*"))
  sigma.mat      <- cov(fpcObj$xiEst)

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


