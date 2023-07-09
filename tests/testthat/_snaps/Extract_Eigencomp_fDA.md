# Extract_Eigencomp_fDA() throws an error when obs.design is specified wrongly

    Code
      fPASS::Extract_Eigencomp_fDA(obs.design = list(design = "functional",
        visit.schedule = seq(0.1, 0.9, length.out = 7), visit.window = 0.05),
      mean_diff_fnm = "mean.diff", cov.type = "NS", cov.par = list(cov.obj = NULL,
        eigen.comp = list(eig.val = c(1, 0.5), eig.obj = eig.fun.vec)), sigma2.e = 0.001,
      nobs_per_subj = 8, missing_type = "nomiss", missing_percent = 0, eval_SS = 100,
      alloc.ratio = c(1, 1), fpca_method = "fpca.sc", data.driven.scores = FALSE,
      mean_diff_add_args = list(), fpca_optns = list(pve = 0.95))
    Error <expectation_failure>
      Names of obs.design[!names(obs.design) %in% "design"] ('visit.schedule', 'visit.window') don't match 'fun.domain'
      Names of obs.design must be 'design' and 'fun.domain'.

