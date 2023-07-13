
testthat::test_that("Extract_Eigencomp_fDA() estimates the eigenfunctions correctly
                    using fpca.sc and functional design",{
  mean.diff    <- function(t) {1 * (t^3)};
  eig.fun <- function(t, k) {
    if (k==1) ef <- sqrt(2)*sin(2*pi*t)
    else if (k==2) ef <- sqrt(2)*cos(2*pi*t)
    return(ef)
  }
  eig.fun.vec  <- function(t){cbind(eig.fun(t, 1),eig.fun(t, 2))}
    eigencomp <- fPASS::Extract_Eigencomp_fDA(obs.design = list(design = "functional", fun.domain = c(0,1)),
                          mean_diff_fnm = "mean.diff", cov.type = "NS",
                          cov.par = list("cov.obj" = NULL, "eigen.comp" = list("eig.val" = c(1, 0.5), "eig.obj" = eig.fun.vec)),
                          sigma2.e = 0.001, nobs_per_subj = 4:7,
                          missing_type = "nomiss",
                          missing_percent = 0, eval_SS = 5000,
                          alloc.ratio = c(1,1),
                          fpca_method = "fpca.sc", data.driven.scores = FALSE,
                          mean_diff_add_args = list(), fpca_optns = list("pve" = 0.95))
    testthat::expect_lte(max(colSums(sweep({(abs(eigencomp$est_eigenfun[,1:2]) -
                                           abs(eig.fun.vec(eigencomp$working.grid)))^2},1,eigencomp$weights, FUN="*")))
                         , 0.01)
    })

testthat::test_that("Extract_Eigencomp_fDA() estimates the eigenfunctions correctly
                    using fpca.sc and longitudinal design",{
                      mean.diff    <- function(t) {1 * (t^3)};
                      eig.fun <- function(t, k) {
                        if (k==1) ef <- sqrt(2)*sin(2*pi*t)
                        else if (k==2) ef <- sqrt(2)*cos(2*pi*t)
                        return(ef)
                      }
                      eig.fun.vec  <- function(t){cbind(eig.fun(t, 1),eig.fun(t, 2))}
                      eigencomp <- fPASS::Extract_Eigencomp_fDA(obs.design = obs.design <- list("design" = "longitudinal",
                                                                                         "visit.schedule" = seq(0.1, 0.9, length.out=7),
                                                                                         "visit.window" = 0.05),
                                                         mean_diff_fnm = "mean.diff", cov.type = "NS",
                                                         cov.par = list("cov.obj" = NULL, "eigen.comp" = list("eig.val" = c(1, 0.5), "eig.obj" = eig.fun.vec)),
                                                         sigma2.e = 0.001, nobs_per_subj = 8,
                                                         missing_type = "nomiss",
                                                         missing_percent = 0, eval_SS = 5000,
                                                         alloc.ratio = c(1,1),
                                                         fpca_method = "fpca.sc", data.driven.scores = FALSE,
                                                         mean_diff_add_args = list(), fpca_optns = list("pve" = 0.95))
                      testthat::expect_lte(max(colSums(sweep({(abs(eigencomp$est_eigenfun[,1:2]) -
                                                                 abs(eig.fun.vec(eigencomp$working.grid)))^2},1,eigencomp$weights, FUN="*")))
                                           , 0.01)
                    })


testthat::test_that("Extract_Eigencomp_fDA() estimates the eigenfunctions correctly
                    using face and functional design",{
                      mean.diff    <- function(t) {1 * (t^3)};
                      eig.fun <- function(t, k) {
                        if (k==1) ef <- sqrt(2)*sin(2*pi*t)
                        else if (k==2) ef <- sqrt(2)*cos(2*pi*t)
                        return(ef)
                      }
                      eig.fun.vec  <- function(t){cbind(eig.fun(t, 1),eig.fun(t, 2))}
                      eigencomp <- fPASS::Extract_Eigencomp_fDA(obs.design = list(design = "functional", fun.domain = c(0,1)),
                                                         mean_diff_fnm = "mean.diff", cov.type = "NS",
                                                         cov.par = list("cov.obj" = NULL, "eigen.comp" = list("eig.val" = c(1, 0.5), "eig.obj" = eig.fun.vec)),
                                                         sigma2.e = 0.001, nobs_per_subj = 4:7,
                                                         missing_type = "nomiss",
                                                         missing_percent = 0, eval_SS = 5000,
                                                         alloc.ratio = c(1,1),
                                                         fpca_method = "fpca.sc", data.driven.scores = FALSE,
                                                         mean_diff_add_args = list(), fpca_optns = list("pve" = 0.95
                                                                                                        #calculate.scores = FALSE,
                                                                                                        #center = TRUE, knots = 10
                                                                                                        ))
                      testthat::expect_lte(max(colSums(sweep({(abs(eigencomp$est_eigenfun[,1:2]) -
                                                                 abs(eig.fun.vec(eigencomp$working.grid)))^2},1,eigencomp$weights, FUN="*")))
                                           , 0.1)
                    })

testthat::test_that("Extract_Eigencomp_fDA() estimates the eigenfunctions correctly
                    using face and longitudinal design",{
                      mean.diff    <- function(t) {1 * (t^3)};
                      eig.fun <- function(t, k) {
                        if (k==1) ef <- sqrt(2)*sin(2*pi*t)
                        else if (k==2) ef <- sqrt(2)*cos(2*pi*t)
                        return(ef)
                      }
                      eig.fun.vec  <- function(t){cbind(eig.fun(t, 1),eig.fun(t, 2))}
                      eigencomp <- fPASS::Extract_Eigencomp_fDA(obs.design = obs.design <- list("design" = "longitudinal",
                                                                                         "visit.schedule" = seq(0.1, 0.9, length.out=7),
                                                                                         "visit.window" = 0.05),
                                                         mean_diff_fnm = "mean.diff", cov.type = "NS",
                                                         cov.par = list("cov.obj" = NULL, "eigen.comp" = list("eig.val" = c(1, 0.5), "eig.obj" = eig.fun.vec)),
                                                         sigma2.e = 0.001, nobs_per_subj = 8,
                                                         missing_type = "nomiss",
                                                         missing_percent = 0, eval_SS = 5000,
                                                         alloc.ratio = c(1,1),
                                                         fpca_method = "fpca.sc", data.driven.scores = FALSE,
                                                         mean_diff_add_args = list(), fpca_optns = list("pve" = 0.95
                                                                                                        #calculate.scores = FALSE,
                                                                                                        #knots=10
                                                                                                        ))
                      testthat::expect_lte(max(colSums(sweep({(abs(eigencomp$est_eigenfun[,1:2]) -
                                                                 abs(eig.fun.vec(eigencomp$working.grid)))^2},1,eigencomp$weights, FUN="*")))
                                           , 0.1)
                    })


# Error checking
testthat::test_that("Extract_Eigencomp_fDA() throws an error when obs.design is specified wrongly",{
  alloc.ratio  <- c(1,1)
  mean.diff    <- function(t) {1 * (t^3)};
  eig.fun <- function(t, k) {
    if (k==1) ef <- sqrt(2)*sin(2*pi*t)
    else if (k==2) ef <- sqrt(2)*cos(2*pi*t)
    return(ef)
  }
  eig.fun.vec  <- function(t){cbind(eig.fun(t, 1),eig.fun(t, 2))}
  ttestthat::expect_error(
    fPASS::Extract_Eigencomp_fDA(obs.design = list("design" = "functional",
                                            "visit.schedule" = seq(0.1, 0.9, length.out=7),
                                            "visit.window" = 0.05),
                          mean_diff_fnm = "mean.diff", cov.type = "NS",
                          cov.par = list("cov.obj" = NULL, "eigen.comp" = list("eig.val" = c(1, 0.5), "eig.obj" = eig.fun.vec)),
                          sigma2.e = 0.001, nobs_per_subj = 8,
                          missing_type = "nomiss",
                          missing_percent = 0, eval_SS = 100,
                          alloc.ratio = c(1,1),
                          fpca_method = "fpca.sc", data.driven.scores = FALSE,
                          mean_diff_add_args = list(), fpca_optns = list("pve" = 0.95))
  )
  testthat::expect_error(
    fPASS::Extract_Eigencomp_fDA(obs.design = list("design" = "functional",
                                            "visit.schedule" = seq(0.1, 0.9, length.out=7),
                                            "visit.window" = 0.05),
                          mean_diff_fnm = "mean.diff", cov.type = "NS",
                          cov.par = list("cov.obj" = NULL, "eigen.comp" = list("eig.val" = c(1, 0.5), "eig.obj" = eig.fun.vec)),
                          sigma2.e = 0.001, nobs_per_subj = 8,
                          missing_type = "nomiss",
                          missing_percent = 0, eval_SS = 100,
                          alloc.ratio = c(1,1),
                          fpca_method = "fpca.sc", data.driven.scores = FALSE,
                          mean_diff_add_args = list(), fpca_optns = list("pve" = 0.95))
  )
  })

