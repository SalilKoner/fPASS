#' Compute quadrature weights
#'
#' Utility function for numerical integration.
#' @param argvals function arguments.
#' @param method quadrature method. Can be either \code{trapedoidal} or \code{midpoint}.
#' @return a vector of quadrature weights for the points supplied in \code{argvals}.
#' @author Clara Happ, with modifications by Philip Reiss
#' @keywords internal
#' @noRd

# TODO: check about 'midpoint'
# TODO: have this function called by lf, af, etc.
# TODO: harmonize with quadrature implemented elsewhere: 'simpson' and 'riemann' options?

quadWeights<- function(argvals, method = "trapezoidal")
{
  ret <- switch(method,
                trapezoidal = {D <- length(argvals)
                1/2*c(argvals[2] - argvals[1], argvals[3:D] -argvals[1:(D-2)], argvals[D] - argvals[D-1])},
                midpoint = c(0,diff(argvals)),  # why is this called 'midpoint'???
                stop("function quadWeights: choose either trapezoidal or midpoint quadrature rule"))

  return(ret)
}

#' Functional principal components analysis by smoothed covariance
#'
#' @description
#'
#' `r lifecycle::badge("stable")`
#'
#' @description
#' This function computes a FPC decomposition for a set of observed curves,
#' which may be sparsely observed and/or measured with error. A mixed model
#' framework is used to estimate curve-specific scores and variances.
#'
#' @details This function is emulated from the [refund::fpca.sc()] function
#' where the estimation of covariance surface and the eigenfunctions are
#' exactly as that of [refund::fpca.sc()], but it rectifies the computational
#' intricacies involved in the estimation of `shrinkage` scores, and fixes
#' the issue of NA values in the score estimation when the measurement error
#' variance is estimated to be zero. Moreover, since this function is written
#' purely for the purpose of using it in the [fPASS::Extract_Eigencomp_fDA()]
#' function, where we do not need the usage of the arguments `var` and `simul`
#' and `sim.alpha` at all, we have deleted those arguments in the
#' [fPASS::fpca_sc()] function.
#'
#' The functional data must be supplied as either \itemize{ \item an \eqn{n
#' \times d} matrix \code{Y}, each row of which is one functional observation,
#' with missing values allowed; or \item a data frame \code{ydata}, with
#' columns \code{'.id'} (which curve the point belongs to, say \eqn{i}),
#' \code{'.index'} (function argument such as time point \eqn{t}), and
#' \code{'.value'} (observed function value \eqn{Y_i(t)}).}
#' @inheritParams refund::fpca.sc
#' @return An object of class \code{fpca} containing:
#' \item{Yhat}{FPC approximation (projection onto leading components)
#' of \code{Y.pred} if specified, or else of \code{Y}.}
#' \item{Y}{the observed data}\item{scores}{\eqn{n
#' \times npc} matrix of estimated FPC scores.} \item{mu}{estimated mean
#' function (or a vector of zeroes if \code{center==FALSE}).} \item{efunctions
#' }{\eqn{d \times npc} matrix of estimated eigenfunctions of the functional
#' covariance, i.e., the FPC basis functions.} \item{evalues}{estimated
#' eigenvalues of the covariance operator, i.e., variances of FPC scores.}
#' \item{npc }{number of FPCs: either the supplied \code{npc}, or the minimum
#' number of basis functions needed to explain proportion \code{pve} of the
#' variance in the observed curves.} \item{argvals}{argument values of
#' eigenfunction evaluations} \item{sigma2}{estimated measurement error
#' variance.} \item{diag.var}{diagonal elements of the covariance matrices for
#' each estimated curve.} \item{VarMats}{a list containing the estimated
#' covariance matrices for each curve in \code{Y}.} \item{crit.val}{estimated
#' critical values for constructing simultaneous confidence intervals.}
#' @author Salil Koner \email{salil.koner@@duke.edu}
#' @references Di, C., Crainiceanu, C., Caffo, B., and Punjabi, N. (2009).
#' Multilevel functional principal component analysis. \emph{Annals of Applied
#' Statistics}, 3, 458--488.
#'
#' Goldsmith, J., Greven, S., and Crainiceanu, C. (2013). Corrected confidence
#' bands for functional data using principal components. \emph{Biometrics},
#' 69(1), 41--51.
#'
#' Staniswalis, J. G., and Lee, J. J. (1998). Nonparametric regression
#' analysis of longitudinal data. \emph{Journal of the American Statistical
#' Association}, 93, 1403--1418.
#'
#' Yao, F., Mueller, H.-G., and Wang, J.-L. (2005). Functional data analysis
#' for sparse longitudinal data. \emph{Journal of the American Statistical
#' Association}, 100, 577--590.
#' @examples
#'
#' if(rlang::is_installed("refund")){
#'   library(refund)
#'   data(cd4)
#'   Fit.MM = fpca_sc(refund::cd4, pve = 0.95)
#' }
#'
#' # input a data frame instead of a matrix
#' nid <- 20
#' nobs <- sample(10:20, nid, rep=TRUE)
#' ydata <- data.frame(
#'     .id = rep(1:nid, nobs),
#'     .index = round(runif(sum(nobs), 0, 1), 3))
#' ydata$.value <- unlist(tapply(ydata$.index,
#'                               ydata$.id,
#'                               function(x)
#'                                   runif(1, -.5, .5) +
#'                                   dbeta(x, runif(1, 6, 8), runif(1, 3, 5))
#'                               )
#'                        )
#'
#' Fit.MM = fpca_sc(ydata=ydata)
#'
#'
#' @keywords internal
#' @importFrom Matrix nearPD Matrix t as.matrix
#' @importFrom mgcv gam predict.gam
#' @importFrom gamm4 gamm4
#' @importFrom stats complete.cases predict rchisq weighted.mean pf qf
#' @rdname fpca
#' @export
fpca_sc <- function(Y = NULL, ydata = NULL, Y.pred = NULL, argvals = NULL, random.int = FALSE,
                    nbasis = 10, pve = 0.95, npc = NULL,
                    useSymm = FALSE, makePD = FALSE, center = TRUE,
                    cov.est.method = 2, integration = "trapezoidal") {

  stopifnot((!is.null(Y) && is.null(ydata)) || (is.null(Y) && !is.null(ydata)))

  # if data.frame version of ydata is provided
  sparseOrNongrid <- !is.null(ydata)
  if (sparseOrNongrid) {
    stopifnot(ncol(ydata) == 3)
    stopifnot(c(".id", ".index", ".value") == colnames(ydata))
    stopifnot(is.null(argvals))
    Y = irreg2mat(ydata)
    argvals = sort(unique(ydata$.index))
  }

  if (is.null(Y.pred))
    Y.pred = Y
  D = NCOL(Y)
  I = NROW(Y)
  I.pred = NROW(Y.pred)

  if (is.null(argvals))
    argvals = seq(0, 1, length = D)

  d.vec = rep(argvals, each = I)
  id = rep(1:I, rep(D, I))

  if (center) {
    if (random.int) {
      ri_data <- data.frame(y = as.vector(Y), d.vec = d.vec, id = factor(id))
      gam0 = gamm4::gamm4(y ~ s(d.vec, k = nbasis), random = ~(1 | id), data = ri_data)$gam
      rm(ri_data)
    } else gam0 = mgcv::gam(as.vector(Y) ~ s(d.vec, k = nbasis))
    mu = predict(gam0, newdata = data.frame(d.vec = argvals))
    Y.tilde = Y - matrix(mu, I, D, byrow = TRUE)
  } else {
    Y.tilde = Y
    mu = rep(0, D)
  }

  if (cov.est.method == 2) {
    # smooth raw covariance estimate
    cov.sum = cov.count = cov.mean = matrix(0, D, D)
    for (i in 1:I) {
      obs.points = which(!is.na(Y[i, ]))
      cov.count[obs.points, obs.points] = cov.count[obs.points, obs.points] +
        1
      cov.sum[obs.points, obs.points] = cov.sum[obs.points, obs.points] + tcrossprod(Y.tilde[i,
                                                                                             obs.points])
    }
    G.0 = ifelse(cov.count == 0, NA, cov.sum/cov.count)
    diag.G0 = diag(G.0)
    diag(G.0) = NA
    if (!useSymm) {
      row.vec = rep(argvals, each = D)
      col.vec = rep(argvals, D)
      npc.0 = matrix(predict(mgcv::gam(as.vector(G.0) ~ te(row.vec, col.vec, k = nbasis),
                                 weights = as.vector(cov.count)), newdata = data.frame(row.vec = row.vec,
                                                                                       col.vec = col.vec)), D, D)
      npc.0 = (npc.0 + t(npc.0))/2
    } else {
      use <- upper.tri(G.0, diag = TRUE)
      use[2, 1] <- use[ncol(G.0), ncol(G.0) - 1] <- TRUE
      usecov.count <- cov.count
      usecov.count[2, 1] <- usecov.count[ncol(G.0), ncol(G.0) - 1] <- 0
      usecov.count <- as.vector(usecov.count)[use]
      use <- as.vector(use)
      vG.0 <- as.vector(G.0)[use]
      row.vec <- rep(argvals, each = D)[use]
      col.vec <- rep(argvals, times = D)[use]
      mCov <- mgcv::gam(vG.0 ~ te(row.vec, col.vec, k = nbasis), weights = usecov.count)
      npc.0 <- matrix(NA, D, D)
      spred <- rep(argvals, each = D)[upper.tri(npc.0, diag = TRUE)]
      tpred <- rep(argvals, times = D)[upper.tri(npc.0, diag = TRUE)]
      smVCov <- predict(mCov, newdata = data.frame(row.vec = spred, col.vec = tpred))
      npc.0[upper.tri(npc.0, diag = TRUE)] <- smVCov
      npc.0[lower.tri(npc.0)] <- t(npc.0)[lower.tri(npc.0)]
    }
  } else if (cov.est.method == 1) {
    # smooth y(s1)y(s2) values to obtain covariance estimate
    row.vec = col.vec = G.0.vec = c()
    cov.sum = cov.count = cov.mean = matrix(0, D, D)
    for (i in 1:I) {
      obs.points = which(!is.na(Y[i, ]))
      temp = tcrossprod(Y.tilde[i, obs.points])
      diag(temp) = NA
      row.vec = c(row.vec, rep(argvals[obs.points], each = length(obs.points)))
      col.vec = c(col.vec, rep(argvals[obs.points], length(obs.points)))
      G.0.vec = c(G.0.vec, as.vector(temp))
      # still need G.O raw to calculate to get the raw to get the diagonal
      cov.count[obs.points, obs.points] = cov.count[obs.points, obs.points] +
        1
      cov.sum[obs.points, obs.points] = cov.sum[obs.points, obs.points] + tcrossprod(Y.tilde[i,
                                                                                             obs.points])
    }
    row.vec.pred = rep(argvals, each = D)
    col.vec.pred = rep(argvals, D)
    npc.0 = matrix(predict(mgcv::gam(G.0.vec ~ mgcv::te(row.vec, col.vec, k = nbasis)), newdata = data.frame(row.vec = row.vec.pred,
                                                                                                 col.vec = col.vec.pred)), D, D)
    npc.0 = (npc.0 + t(npc.0))/2
    G.0 = ifelse(cov.count == 0, NA, cov.sum/cov.count)
    diag.G0 = diag(G.0)
  }

  if (makePD) {
    npc.0 <- {
      tmp <- Matrix::nearPD(npc.0, corr = FALSE, keepDiag = FALSE, do2eigen = TRUE,
                            trace = TRUE)
      as.matrix(tmp$mat)
    }
  }
  ### numerical integration for calculation of eigenvalues (see Ramsay & Silverman,
  ### Chapter 8)
  w <- quadWeights(argvals, method = integration)
  Wsqrt <- diag(sqrt(w))
  Winvsqrt <- diag(1/(sqrt(w)))
  V <- Wsqrt %*% npc.0 %*% Wsqrt
  evalues = eigen(V, symmetric = TRUE, only.values = TRUE)$values
  ###
  evalues = replace(evalues, which(evalues <= 0), 0)
  npc = ifelse(is.null(npc), min(which(cumsum(evalues)/sum(evalues) > pve)), npc)
  efunctions = matrix(Winvsqrt %*% eigen(V, symmetric = TRUE)$vectors[, seq(len = npc)],
                      nrow = D, ncol = npc)
  evalues = eigen(V, symmetric = TRUE, only.values = TRUE)$values[1:npc]  # use correct matrix for eigenvalue problem
  cov.hat = efunctions %*% tcrossprod(diag(evalues, nrow = npc, ncol = npc), efunctions)
  ### numerical integration for estimation of sigma2
  T.len <- argvals[D] - argvals[1]  # total interval length
  T1.min <- min(which(argvals >= argvals[1] + 0.25 * T.len))  # left bound of narrower interval T1
  T1.max <- max(which(argvals <= argvals[D] - 0.25 * T.len))  # right bound of narrower interval T1
  DIAG = (diag.G0 - diag(cov.hat))[T1.min:T1.max]  # function values
  w2 <- quadWeights(argvals[T1.min:T1.max], method = integration)
  sigma2 <- max(weighted.mean(DIAG, w = w2, na.rm = TRUE), 0)

  ####
  D.inv = diag(1/evalues, nrow = npc, ncol = npc)
  Z = efunctions
  Y.tilde = Y.pred - matrix(mu, I.pred, D, byrow = TRUE)
  scores = matrix(NA, nrow = I.pred, ncol = npc)
  if (sigma2 == 0) {
    message("Measurement error variance is estimated to be zero, setting to 1e-6 \n")
    sigma2 <- 1e-6
  }
  for (i.subj in 1:I.pred) {
    obs.points = which(!is.na(Y.pred[i.subj, ]))
    if (sigma2 == 0 & length(obs.points) < npc)
      stop("Measurement error estimated to be zero and there are fewer observed points than PCs; scores cannot be estimated.")
    Zcur = matrix(Z[obs.points, ], nrow = length(obs.points), ncol = dim(Z)[2])
    ZtZ_sD.inv = solve(crossprod(Zcur) + sigma2 * D.inv, t(Zcur) %*% Y.tilde[i.subj, obs.points] )
    scores[i.subj, ] = diag(evalues/sigma2, nrow = npc, ncol = npc) %*%
                       (crossprod(Zcur, Y.tilde[i.subj, obs.points]) - (crossprod(Zcur) %*% ZtZ_sD.inv))
  }

  ret.objects = c("scores", "mu", "efunctions", "evalues", "npc",
                  "argvals", "sigma2")
  ret = lapply(1:length(ret.objects), function(u) get(ret.objects[u]))
  names(ret) = ret.objects
  class(ret) = "fpca"
  return(ret)
}
