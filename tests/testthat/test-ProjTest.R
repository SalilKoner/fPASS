target.SS    <- NULL; target.power <- 0.75; alpha.fix <- 0.05;
alloc.ratio  <- c(1,1)
mean.diff    <- function(t) {1.5 * (t^3)};
eig.fun      <- function(t) { sqrt(2)*c(sin(2*pi*t), cos(2*pi*t))}
eigen.comp   <- list("eig.val" = c(1, 0.5), "eig.obj" = eig.fun, "eig.compute" = T)
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


Samp.Size    <- sum(floor(Sample_size$opt.SS)+1); target.power <- 0.7; alpha.fix <- 0.05;
alloc.ratio  <- c(1,1)
mean.diff    <- function(t) {1.5 * (t^3)};
eig.fun      <- function(t) { sqrt(2)*c(sin(2*pi*t), cos(2*pi*t))}
eigen.comp   <- list("eig.val" = c(1, 0.5), "eig.obj" = eig.fun, "eig.compute" = TRUE)
cov.par      <- list("cov.obj" = NULL, "eigen.comp" = eigen.comp)
sigma2.e     <- 0.001; nobs_per_subj <- 4:7;
missing_type <- "nomiss"; missing_percent <- 0;
fpca_method  <- "fpca.sc"
fpca_optns   <- list("FVEthreshold" = 0.95, "methodSelectK" = "FVE")
ngrid        <- 101

decision     <- function(rep){

  pVal         <- fPASS::TS_ProjTest_sUfDA(Samp.Size = Samp.Size,
                                    alloc.ratio = alloc.ratio,
                                    ngrid = ngrid, mean.diff = mean.diff,
                                    cov.type = "NS", cov.par = cov.par, sigma2.e = sigma2.e,
                                    nobs_per_subj = nobs_per_subj, missing_type = missing_type,
                                    missing_percent = missing_percent,
                                    fpca_method = fpca_method, fpca_optns =fpca_optns)
  as.numeric(pVal < alpha.fix)
}



library(foreach); library(doSNOW) ; library(doParallel);
cl                  <- makeSOCKcluster(8)
registerDoSNOW(cl)
progress            <- function(nfin, tag) { cat(sprintf('tasks completed:
                                                         %d; tag: %d\n', nfin,
                                                         tag)) }
opts                <- list(progress=progress)
packages_req        <- c("mgcv", "tidyverse", "refund", "face", "Hotelling", "fdapace", "fPASS")
reject              <- foreach(itid=1:2000, #.options.snow=opts,
                               .packages =packages_req) %dopar% {
                                 try(decision(rep=itid))}

stopCluster(cl)
notError        <- sapply(reject, function(el) class(el) != "try-error")
rejections.prob <- mean(unlist(reject[notError]), na.rm = TRUE)


# Varying target power, delta, covariance parameter, both ST and NS:
# Take the minimum and maximum of the sample sizes,
# for NS, take both small and medium sample sizes
