target.SS <- NULL; target.power <- 0.8; alpha.fix <- 0.05; mean.diff <- function(t) {t};
cor.str <- nlme::corExp(1, form = ~ time | Subject);
sigma2 <- 1; sigma2.e <- 0.25; nobs_per_subj <- 6;
missing_type <- "constant"; missing_percent <- 0.01;
fpca_optns <- list("FVEthreshold" = 0.95)


target.SS    <- NULL; target.power <- 0.8; alpha.fix <- 0.05;
alloc.ratio  <- c(1,1)
mean.diff    <- function(t) {1 * (t^3)};
eig.fun      <- function(t) { sqrt(2)*c(sin(2*pi*t), cos(2*pi*t))}
eigen.comp   <- list("eig.val" = c(1, 0.5), "eig.obj" = eig.fun, "eig.compute" = TRUE)
cov.par      <- list("cov.obj" = NULL, "eigen.comp" = eigen.comp)
sigma2.e     <- 0.001; nobs_per_subj <- 4:7;
missing_type <- "nomiss"; missing_percent <- 0;
fpca_method  <- "fpca.sc"
fpca_optns   <- list("FVEthreshold" = 0.95, "methodSelectK" = "FVE")
ngrid        <- 101

Sample_size  <- PASS_TS_ProjTest_sUfDA(target.SS = target.SS, target.power = target.power,
                                      alpha.fix = alpha.fix, alloc.ratio = alloc.ratio,
                                      ngrid = ngrid, mean.diff = mean.diff,
                                      cov.type = "NS", cov.par = cov.par, sigma2.e = sigma2.e,
                                      nobs_per_subj = nobs_per_subj, missing_type = missing_type,
                                      missing_percent = missing_percent,
                                      fpca_method = fpca_method, fpca_optns =fpca_optns)



est.power   <- TS_ProjTest_sUfDA(Samp.Size = 500, alpha.fix = alpha.fix,
                                 alloc.ratio = c(2,1),
                                 ngrid = ngrid, mean.diff = mean.diff,
                                 cov.type = "NS", cov.par = cov.par, sigma2.e = sigma2.e,
                                 nobs_per_subj = nobs_per_subj, missing_type = missing_type,
                                 missing_percent = missing_percent,
                                 fpca_method = fpca_method, fpca_optns =fpca_optns)
