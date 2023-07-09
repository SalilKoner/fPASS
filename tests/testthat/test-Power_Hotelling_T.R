
testthat::test_that("pHotellingT() provides the correct power
                    when H0 is not true
                    and the covariance are unequal",{
            testthat::skip_if_not_installed("Hotelling")
            B           <- 10000
            k           <- 4
            n2          <- 60
            n1_by_n2    <- 2
            n1          <- n1_by_n2 * n2
            mu1         <- rep(0,k)
            del         <- 0.4
            mu2         <- mu1 + rep(del, k) # rep(0.19,k)  # 0.23 (0.9), 0.18 (0.7) 0.20 (0.8)
            sig1        <- diag(k)
            sig2        <- sig1  + del*toeplitz(c(1,rep(0.9, k-1)))
            cutoff      <- seq(0,30, length.out=20)
            the_cdf     <- round(fPASS::pHotellingT(cutoff, n1+n2, mu1 - mu2,
                                         sig1, sig2, alloc.ratio=c(2,1),
                                         lower.tail=FALSE),3)
            emp_samples <- replicate(B, {y1      <- MASS::mvrnorm(n=n1, mu=mu1, Sigma=sig1)
            y2          <- MASS::mvrnorm(n=n2, mu=mu2, Sigma=sig2)
            ht.test     <- Hotelling::hotelling.test(y1, y2, var.equal = TRUE)
            ht.test$stats$statistic })
            emp_cdf     <- round(sapply(cutoff, function(cc) mean(emp_samples > cc, na.rm=TRUE) ),3)
            obs.error   <- quantile(round(abs(emp_cdf - the_cdf),2), prob=0.95, names = F)
            exp.err     <- round(2*0.5*sqrt(1/B),2)
            testthat::expect_lte(obs.error, exp.err)
          })

testthat::test_that("pHotellingT() provides the correct power
                    when H0 is true
                    and the covariance are equal",{
                      testthat::skip_if_not_installed("Hotelling")
                      B           <- 10000
                      k           <- 4
                      n2          <- 40
                      n1_by_n2    <- 2
                      n1          <- n1_by_n2 * n2
                      mu1         <- rep(0,k)
                      del         <- 0
                      mu2         <- mu1 + rep(del, k) # rep(0.19,k)  # 0.23 (0.9), 0.18 (0.7) 0.20 (0.8)
                      sig1        <- diag(k)
                      sig2        <- sig1  + del*toeplitz(c(1,rep(0.9, k-1)))
                      cutoff      <- seq(0,30, length.out=20)
                      the_cdf     <- round(fPASS::pHotellingT(cutoff, n1+n2, mu1 - mu2,
                                           sig1, sig2, alloc.ratio=c(2,1),
                                           lower.tail=FALSE),3)
                      emp_samples <- replicate(B, {y1      <- MASS::mvrnorm(n=n1, mu=mu1, Sigma=sig1)
                      y2          <- MASS::mvrnorm(n=n2, mu=mu2, Sigma=sig2)
                      ht.test     <- Hotelling::hotelling.test(y1, y2, var.equal = TRUE)
                      ht.test$stats$statistic })
                      emp_cdf     <- round(sapply(cutoff, function(cc) mean(emp_samples > cc, na.rm=TRUE) ),3)
                      obs.error   <- quantile(round(abs(emp_cdf - the_cdf),2), prob=0.95, names = F)
                      exp.err     <- round(2*0.5*sqrt(1/B),2)
                      testthat::expect_lte(obs.error, exp.err)
                    })

testthat::test_that("pHotellingT() provides the correct power
                    when H0 is true
                    and the covariance are unequal",{
                      testthat::skip_if_not_installed("Hotelling")
                      B           <- 10000
                      k           <- 4
                      n2          <- 60
                      n1_by_n2    <- 2
                      n1          <- n1_by_n2 * n2
                      mu1         <- rep(0,k)
                      del         <- 0
                      mu2         <- mu1 + rep(del, k) # rep(0.19,k)  # 0.23 (0.9), 0.18 (0.7) 0.20 (0.8)
                      sig1        <- diag(k)
                      sig2        <- toeplitz(c(1,rep(0.9, k-1)))
                      cutoff      <- seq(0,30, length.out=20)
                      the_cdf     <- round(fPASS::pHotellingT(cutoff, n1+n2, mu1 - mu2,
                                                   sig1, sig2, alloc.ratio=c(2,1),
                                                   lower.tail=FALSE),3)
                      emp_samples <- replicate(B, {y1      <- MASS::mvrnorm(n=n1, mu=mu1, Sigma=sig1)
                      y2          <- MASS::mvrnorm(n=n2, mu=mu2, Sigma=sig2)
                      ht.test     <- Hotelling::hotelling.test(y1, y2, var.equal = TRUE)
                      ht.test$stats$statistic })
                      emp_cdf     <- round(sapply(cutoff, function(cc) mean(emp_samples > cc, na.rm=TRUE) ),3)
                      obs.error   <- quantile(round(abs(emp_cdf - the_cdf),2), prob=0.95, names = F)
                      exp.err     <- round(2*0.5*sqrt(1/B),2)
                      testthat::expect_lte(obs.error, exp.err)
                    })

testthat::test_that("pHotellingT() provides the correct power
                    when H0 is not true
                    and the covariance are equal",{
                      testthat::skip_if_not_installed("Hotelling")
                      B           <- 10000
                      k           <- 4
                      n2          <- 60
                      n1_by_n2    <- 2
                      n1          <- n1_by_n2 * n2
                      mu1         <- rep(0,k)
                      del         <- 0.4
                      mu2         <- mu1 + rep(del, k) # rep(0.19,k)  # 0.23 (0.9), 0.18 (0.7) 0.20 (0.8)
                      sig1        <- diag(k)
                      sig2        <- sig1
                      cutoff      <- seq(0,30, length.out=20)
                      the_cdf     <- round(fPASS::pHotellingT(cutoff, n1+n2, mu1 - mu2,
                                               sig1, sig2, alloc.ratio=c(2,1),
                                               lower.tail=FALSE),3)
                      emp_samples <- replicate(B, {y1      <- MASS::mvrnorm(n=n1, mu=mu1, Sigma=sig1)
                      y2          <- MASS::mvrnorm(n=n2, mu=mu2, Sigma=sig2)
                      ht.test     <- Hotelling::hotelling.test(y1, y2, var.equal = TRUE)
                      ht.test$stats$statistic })
                      emp_cdf     <- round(sapply(cutoff, function(cc) mean(emp_samples > cc, na.rm=TRUE) ),3)
                      obs.error   <- quantile(round(abs(emp_cdf - the_cdf),2), prob=0.95, names = F)
                      exp.err     <- round(2*0.5*sqrt(1/B),2)
                      testthat::expect_lte(obs.error, exp.err)
                    })

# Errors

testthat::test_that("Sim_HotellingT_unequal_var() throws an error when
          dimension of the covariance do not match",{
            testthat::expect_error(fPASS::Sim_HotellingT_unequal_var(50, rep(0.5,2), diag(2), diag(3), c(1,1), 1e4))
            testthat::expect_error(fPASS::Sim_HotellingT_unequal_var(50, rep(0.5,2), diag(3), diag(2), c(1,1), 1e4))
          })

testthat::test_that("Sim_HotellingT_unequal_var() throws an error when
          length of the mean vector do not match with the dimension of covariance",{
            testthat::expect_error(fPASS::Sim_HotellingT_unequal_var(50, rep(0.5,3), diag(2), diag(2), c(1,1), 1e4))
          })

testthat::test_that("Sim_HotellingT_unequal_var() throws an error when
          length of the covariance function is symmetric",{
            testthat::expect_error(fPASS::Sim_HotellingT_unequal_var(50, rep(0.5,3), matrix(rnorm(9),3,3),
                                                              diag(3), c(1,1), 1e4))
            testthat::expect_error(fPASS::Sim_HotellingT_unequal_var(50, rep(0.5,3), diag(3),
                                                              matrix(rnorm(9),3,3), c(1,1), 1e4))
          })

testthat::test_that("Sim_HotellingT_unequal_var() throws an error when
          length of the covariance function is not positive definite",{
            testthat::expect_error(fPASS::Sim_HotellingT_unequal_var(50, rep(0.5,3), -diag(3),
                                                              diag(3), c(1,1), 1e4))
            testthat::expect_error(fPASS::Sim_HotellingT_unequal_var(50, rep(0.5,3), diag(3),
                                                              -diag(3), c(1,1), 1e4))
            testthat::expect_error(fPASS::Sim_HotellingT_unequal_var(50, rep(0.5,3), diag(3),
                                                              diag(0,3), c(1,1), 1e4))
          })

testthat::test_that("Sim_HotellingT_unequal_var() throws an error when
          length of the total_sample_size is not positive and greater than 10",{
            testthat::expect_error(fPASS::Sim_HotellingT_unequal_var(-10, rep(0.5,3), diag(3),
                                                              diag(3), c(1,1), 1e4))
            testthat::expect_error(fPASS::Sim_HotellingT_unequal_var(5, rep(0.5,3), diag(3),
                                                              diag(3), c(1,1), 1e4))
          })

testthat::test_that("Sim_HotellingT_unequal_var() throws an error when
          length of the nsim is not positive",{
            testthat::expect_error(fPASS::Sim_HotellingT_unequal_var(100, rep(0.5,3), diag(3),
                                                              diag(3), c(1,1), 0))
            testthat::expect_error(fPASS::Sim_HotellingT_unequal_var(100, rep(0.5,3), diag(3),
                                                              diag(3), c(1,1), -10))
          })

