
testthat::test_that("PASS_Proj_Test_ufDA() gives power equal to
                    sig.level when null hypothesis is true", {
                      mean.diff    <- function(t, delta) {delta * (t^3)};
                      eig.fun <- function(t, k) {
                        if (k==1) ef <- sqrt(2)*sin(2*pi*t)
                        else if (k==2) ef <- sqrt(2)*cos(2*pi*t)
                        return(ef)
                      }
                      eig.fun.vec  <- function(t){cbind(eig.fun(t, 1),eig.fun(t, 2))}
                      f <- fPASS::PASS_Proj_Test_ufDA(sample_size = 100, target.power = NULL, sig.level = 0.05,
                                                      nobs_per_subj = 4:7, obs.design = list(design = "functional", fun.domain = c(0,1)),
                                                      mean_diff_fnm = "mean.diff",
                                                      cov.type = "NS", cov.par = list("cov.obj" = NULL,
                                                                                      "eigen.comp" = list("eig.val" = c(1, 0.5),
                                                                                                          "eig.obj" = eig.fun.vec)),
                                                      sigma2.e = 0.001,
                                                      missing_type = "nomiss",
                                                      missing_percent = 0,
                                                      eval_SS = 5000, alloc.ratio = c(1,1),
                                                      fpca_method = "fpca.sc", nWgrid = 201,
                                                      mean_diff_add_args=list(delta=0),
                                                      fpca_optns = list(pve = 0.95), npc_to_use = 2,
                                                      nsim =1e4)
                      print(f$power_value)
                      testthat::expect_lte(f$power_value, 0.05 + 4*sqrt(0.05*0.95/5000))
                    })
