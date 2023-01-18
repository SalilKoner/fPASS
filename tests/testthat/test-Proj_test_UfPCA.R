
# Example 1: CS covariance
Samp.Size <- 480; alpha.fix <- 0.05; delta <- 0.78; mean.diff <- function(t) {delta*(t^3)};
cor.str <- nlme::corCompSymm(0.5, form = ~ time | Subject);
sigma2 <- 1.5; sigma2.e <- 0.25; nobs_per_subj <- 4:7;
missing_type <- NULL; missing_percent <- 0;
fpca_optns <- list()
library(foreach); library(doSNOW) ; library(doParallel);
cl                  <- makeSOCKcluster(8)
registerDoSNOW(cl)
progress            <- function(nfin, tag) { cat(sprintf('tasks completed:
                                                         %d; tag: %d\n', nfin, 
                                                         tag)) }
opts                <- list(progress=progress)
packages_req        <- c("mgcv", "tidyverse", "refund", "face", "Hotelling", "fdapace", "fPASS")
reject              <- foreach(itid=1:500, #.options.snow=opts, 
                               .packages =packages_req) %dopar% {
                                Proj_test_UfPCA(Samp.Size, mean.diff, sigma2, 
                                                        cor.str, sigma2.e, nobs_per_subj,
                                                        missing_type, missing_percent, fpca_optns)
                               }

stopCluster(cl)
rejections.prob <- mean(unlist(reject) < alpha.fix)




# Example 1: CS covariance
Samp.Size <- 480; alpha.fix <- 0.05; delta <- 0.78; mean.diff <- function(t) {delta*(t^3)};
cor.str <- nlme::corCompSymm(0.5, form = ~ time | Subject);
sigma2 <- 1.5; sigma2.e <- 0.25; nobs_per_subj <- 6;
missing_type <- NULL; missing_percent <- 0;
fpca_optns <- list()
library(foreach); library(doSNOW) ; library(doParallel);
cl                  <- makeSOCKcluster(8)
registerDoSNOW(cl)
progress            <- function(nfin, tag) { cat(sprintf('tasks completed:
                                                         %d; tag: %d\n', nfin, 
                                                         tag)) }
opts                <- list(progress=progress)
packages_req        <- c("mgcv", "tidyverse", "refund", "face", "Hotelling", "fdapace", "fPASS")
reject              <- foreach(itid=1:500, #.options.snow=opts, 
                               .packages =packages_req) %dopar% {
                                 Proj_test_UfPCA(Samp.Size, mean.diff, sigma2, 
                                                 cor.str, sigma2.e, nobs_per_subj,
                                                 missing_type, missing_percent, fpca_optns)
                               }

stopCluster(cl)
rejections.prob <- mean(unlist(reject) < alpha.fix)



# Example 1: CS covariance
Samp.Size <- 480; alpha.fix <- 0.05; delta <- 0.78; mean.diff <- function(t) {delta*(t^3)};
cor.str <- nlme::corCompSymm(0.5, form = ~ time | Subject);
sigma2 <- 1.5; sigma2.e <- 0.25; nobs_per_subj <- 6;
missing_type <- "constant"; missing_percent <- 0.2;
fpca_optns <- list()
library(foreach); library(doSNOW) ; library(doParallel);
cl                  <- makeSOCKcluster(8)
registerDoSNOW(cl)
progress            <- function(nfin, tag) { cat(sprintf('tasks completed:
                                                         %d; tag: %d\n', nfin, 
                                                         tag)) }
opts                <- list(progress=progress)
packages_req        <- c("mgcv", "tidyverse", "refund", "face", "Hotelling", "fdapace", "fPASS")
reject              <- foreach(itid=1:500, #.options.snow=opts, 
                               .packages =packages_req) %dopar% {
                                 Proj_test_UfPCA(Samp.Size, mean.diff, sigma2, 
                                                 cor.str, sigma2.e, nobs_per_subj,
                                                 missing_type, missing_percent, fpca_optns)
                               }

stopCluster(cl)
rejections.prob <- mean(unlist(reject) < alpha.fix)



# Example 1: CS covariance
Samp.Size <- 480; alpha.fix <- 0.05; delta <- 0.78; mean.diff <- function(t) {delta*(t^3)};
cor.str <- nlme::corCompSymm(0.5, form = ~ time | Subject);
sigma2 <- 1.5; sigma2.e <- 0.25; nobs_per_subj <- 6;
missing_type <- "constant"; missing_percent <- 0.4;
fpca_optns <- list()
library(foreach); library(doSNOW) ; library(doParallel);
cl                  <- makeSOCKcluster(8)
registerDoSNOW(cl)
progress            <- function(nfin, tag) { cat(sprintf('tasks completed:
                                                         %d; tag: %d\n', nfin, 
                                                         tag)) }
opts                <- list(progress=progress)
packages_req        <- c("mgcv", "tidyverse", "refund", "face", "Hotelling", "fdapace", "fPASS")
reject              <- foreach(itid=1:500, #.options.snow=opts, 
                               .packages =packages_req) %dopar% {
                                 Proj_test_UfPCA(Samp.Size, mean.diff, sigma2, 
                                                 cor.str, sigma2.e, nobs_per_subj,
                                                 missing_type, missing_percent, fpca_optns)
                               }

stopCluster(cl)
rejections.prob <- mean(unlist(reject) < alpha.fix)



miss_pct <- c(0.1, 0.25, 0.4, 0.55, 0.7, 0.8)
rejections.prob <- rep(NA, length(miss_pct))
for (i in seq_along(miss_pct)){
  Samp.Size <- 480; alpha.fix <- 0.05; delta <- 0.78; mean.diff <- function(t) {delta*(t^3)};
  cor.str <- nlme::corCompSymm(0.5, form = ~ time | Subject);
  sigma2 <- 1.5; sigma2.e <- 0.25; nobs_per_subj <- 6;
  missing_type <- "constant"; missing_percent <- miss_pct[i];
  fpca_optns <- list()
  library(foreach); library(doSNOW) ; library(doParallel);
  cl                  <- makeSOCKcluster(8)
  registerDoSNOW(cl)
  progress            <- function(nfin, tag) { cat(sprintf('tasks completed:
                                                         %d; tag: %d\n', nfin, 
                                                           tag)) }
  opts                <- list(progress=progress)
  packages_req        <- c("mgcv", "tidyverse", "refund", "face", "Hotelling", "fdapace", "fPASS")
  reject              <- foreach(itid=1:500, #.options.snow=opts, 
                                 .packages =packages_req) %dopar% {
                                   Proj_test_UfPCA(Samp.Size, mean.diff, sigma2, 
                                                   cor.str, sigma2.e, nobs_per_subj,
                                                   missing_type, missing_percent, fpca_optns)
                                 }
  
  stopCluster(cl)
  rejections.prob[i] <- mean(unlist(reject) < alpha.fix)
  
}

library(ggplot2)
qplot(c(0, miss_pct[1:5]), c(0.78, rejections.prob[1:5]), geom="line") + scale
