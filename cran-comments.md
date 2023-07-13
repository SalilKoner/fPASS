## R CMD check results

0 errors | 0 warnings | 2 notes

* This is a new release.

* Examples with CPU (user + system) or elapsed time > 5s
                         user system elapsed
  PASS_Proj_Test_ufDA   10.38   0.56   10.95
  Extract_Eigencomp_fDA  5.25   0.47    5.73

* Vignettes takes about 6-10 minutes on a Macbook Pro 2019 model with 
i9 processor and 16 GB RAM using single core. This could not be avoided in 
the interest of presenting an elaborate understanding of the package. An
already built version of vignette (named as fPASS.html) is provided in
the inst/doc folder. Please allow that time to build the vignette
while running R CMD CHECK or R CMD BUILD. 
