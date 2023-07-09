
testthat::test_that("Power_Proj_Test_ufDA() provides the correct power : case 1",{
  ngrid          <- 101
  interval       <- c(-1,1)
  gauss.quad.pts <- gss::gauss.quad(ngrid,interval) # evaluation points
  working.grid   <- gauss.quad.pts$pt
  mean_fn        <- function(t) {0.4*sin(2*pi*t)}
  mean_vector    <- mean_fn(working.grid)
  eigen_fn       <- function(t, k){ sqrt(2)*{(k==2)*sin(2*pi*t) + (k==1)*cos(2*pi*t)} }
  eigen_matrix   <- cbind(eigen_fn(working.grid,1), eigen_fn(working.grid,2))
  mean_proj      <- sapply(1:2, function(r) integrate(function(x)
                       eigen_fn(x,r)*mean_fn(x), interval[1], interval[2])$value)
  sig1           <- diag(2)
  sig2           <- 2*diag(2)
  alp            <- 0.05
  n              <- 100
  k              <- ncol(eigen_matrix)
  cutoff         <- {(n - 2)*k/(n - k -1)}*qf(1-alp, k, n-k-1)
  set.seed(456)
  mult_power     <- fPASS::pHotellingT(cutoff, n, mean_proj,
                                sig1, sig2, alloc.ratio=c(1,1),
                                lower.tail=FALSE)
  set.seed(456)
  func_power     <- fPASS::Power_Proj_Test_ufDA(total_sample_size=n, argvals=working.grid,
                                         mean_vector = mean_vector, eigen_matrix = eigen_matrix,
                                         scores_var1 = sig1, scores_var2= sig2, weights = gauss.quad.pts$wt,
                                         sig.level=alp, alloc.ratio = c(1,1), npc_to_pick=ncol(eigen_matrix))
  testthat::expect_equal(mult_power, func_power)

})

testthat::test_that("Power_Proj_Test_ufDA() provides the correct power: case 2",{
  ngrid          <- 101
  interval       <- c(-5,4)
  gauss.quad.pts <- gss::gauss.quad(ngrid,interval) # evaluation points
  working.grid   <- gauss.quad.pts$pt
  mean_fn        <- function(t) {0.4*cos(2*pi*t)}
  mean_vector    <- mean_fn(working.grid)
  eigen_fn       <- function(t, k){ sqrt(2)*{(k==2)*sin(2*pi*t) + (k==1)*cos(2*pi*t)} }
  eigen_matrix   <- cbind(eigen_fn(working.grid,1), eigen_fn(working.grid,2))
  mean_proj      <- sapply(1:2, function(r) integrate(function(x)
    eigen_fn(x,r)*mean_fn(x), interval[1], interval[2])$value)
  sig1           <- diag(2)
  sig2           <- 4*diag(2)
  alp            <- 0.1
  n              <- 100
  k              <- ncol(eigen_matrix)
  cutoff         <- {(n - 2)*k/(n - k -1)}*qf(1-alp, k, n-k-1)
  set.seed(1245)
  mult_power     <- fPASS::pHotellingT(cutoff, n, mean_proj,
                                sig1, sig2, alloc.ratio=c(1,1),
                                lower.tail=FALSE)
  set.seed(1245)
  func_power     <- fPASS::Power_Proj_Test_ufDA(total_sample_size=n, argvals=working.grid,
                                         mean_vector = mean_vector, eigen_matrix = eigen_matrix,
                                         scores_var1 = sig1, scores_var2= sig2, weights = gauss.quad.pts$wt,
                                         sig.level=alp, alloc.ratio = c(1,1), npc_to_pick=ncol(eigen_matrix))
  testthat::expect_equal(mult_power, func_power)

})

testthat::test_that("Power_Proj_Test_ufDA() provides the correct power: case 3",{
  ngrid          <- 101
  interval       <- c(0,1)
  gauss.quad.pts <- gss::gauss.quad(ngrid,interval) # evaluation points
  working.grid   <- gauss.quad.pts$pt
  mean_fn        <- function(t) {0.009*exp(2*pi*t)}
  mean_vector    <- mean_fn(working.grid)
  eigen_fn       <- function(t, k){ sqrt(2)*{(k==2)*sin(2*pi*t) + (k==1)*cos(2*pi*t)} }
  eigen_matrix   <- cbind(eigen_fn(working.grid,1), eigen_fn(working.grid,2))
  mean_proj      <- sapply(1:2, function(r) integrate(function(x)
                      eigen_fn(x,r)*mean_fn(x), interval[1], interval[2])$value)
  sig1           <- diag(2)
  sig2           <- 4*diag(2)
  alp            <- 0.1
  n              <- 100
  k              <- ncol(eigen_matrix)
  cutoff         <- {(n - 2)*k/(n - k -1)}*qf(1-alp, k, n-k-1)
  set.seed(1245)
  mult_power     <- fPASS::pHotellingT(cutoff, n, mean_proj,
                                sig1, sig2, alloc.ratio=c(1,1),
                                lower.tail=FALSE)
  set.seed(1245)
  func_power     <- fPASS::Power_Proj_Test_ufDA(total_sample_size=n, argvals=working.grid,
                                         mean_vector = mean_vector, eigen_matrix = eigen_matrix,
                                         scores_var1 = sig1, scores_var2= sig2, weights = gauss.quad.pts$wt,
                                         sig.level=alp, alloc.ratio = c(1,1), npc_to_pick=ncol(eigen_matrix))
  testthat::expect_equal(mult_power, func_power)

})

testthat::test_that("Power_Proj_Test_ufDA() provides the correct power: case 4",{
  ngrid          <- 101
  interval       <- c(0,1)
  gauss.quad.pts <- gss::gauss.quad(ngrid,interval) # evaluation points
  working.grid   <- gauss.quad.pts$pt
  mean_fn        <- function(t) {0*exp(2*pi*t)}
  mean_vector    <- mean_fn(working.grid)
  eigen_fn       <- function(t, k){ sqrt(2)*{(k==2)*sin(2*pi*t) + (k==1)*cos(2*pi*t)} }
  eigen_matrix   <- cbind(eigen_fn(working.grid,1), eigen_fn(working.grid,2))
  mean_proj      <- sapply(1:2, function(r) integrate(function(x)
    eigen_fn(x,r)*mean_fn(x), interval[1], interval[2])$value)
  sig1           <- diag(2)
  sig2           <- diag(2)
  alp            <- 0.1
  n              <- 300
  k              <- ncol(eigen_matrix)
  cutoff         <- {(n - 2)*k/(n - k -1)}*qf(1-alp, k, n-k-1)
  set.seed(1245)
  mult_power     <- fPASS::pHotellingT(cutoff, n, mean_proj,
                                sig1, sig2, alloc.ratio=c(1,1),
                                lower.tail=FALSE)
  set.seed(1245)
  func_power     <- fPASS::Power_Proj_Test_ufDA(total_sample_size=n, argvals=working.grid,
                                         mean_vector = mean_vector, eigen_matrix = eigen_matrix,
                                         scores_var1 = sig1, scores_var2= sig2, weights = gauss.quad.pts$wt,
                                         sig.level=alp, alloc.ratio = c(1,1), npc_to_pick=ncol(eigen_matrix))
  testthat::expect_equal(mult_power, func_power)

})


