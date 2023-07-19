## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval=TRUE
)
#options(rmarkdown.html_vignette.check_title = FALSE)

## ---- eval = FALSE------------------------------------------------------------
#  # Install development version from GitHub
#  devtools::install_github("SalilKoner/fPASS", build_vignettes = FALSE)

## ----setup, eval=TRUE, message=FALSE------------------------------------------
library(fPASS)
library(nlme)
library(face)

## ---- include=F, eval=TRUE----------------------------------------------------
st <- proc.time()

## ---- eval=TRUE---------------------------------------------------------------
obs.design <- list("design" = "longitudinal",
                  "visit.schedule" = c(3,6,12,18,24), # pre-determined visit times, other than baseline
                  "visit.window" = 0.5) # visit window
nobs_per_subj <- length(obs.design$visit.schedule) + 1 # number of observation per subject is 6.

## ---- eval=TRUE---------------------------------------------------------------
# The function will internally create a data set
# and a factor variable Subject denoting subject id.
# You can change time and Subject by any other name as well.
cor.str <- nlme::corCAR1(0.5, form = ~ time | Subject);
# The marginal sd at each time-point is assumed to be 5. 
cov.par  <- list("var" = 25, "cor" = cor.str) 
# We must set cov.type = "ST" for any covariance 
# structure belonging to nlme::corClasses. 
cov.type <- "ST"

## ---- eval=TRUE---------------------------------------------------------------
mean_diff_fn  <- function(t){2*t}
mean_diff_fnm <- "mean_diff_fn" # the name of mean difference function needs to
                                # be specified with the argument. 

## ---- eval=TRUE---------------------------------------------------------------
missing_type <- "nomiss"
missing_percent <- 0 # Has be to a number between 0 and 0.8

## ---- eval=TRUE---------------------------------------------------------------
alloc.ratio <- c(1,1)

## ---- eval=TRUE---------------------------------------------------------------
target.power <- 0.8
sig.level  <- 0.05

## ---- eval=TRUE---------------------------------------------------------------
eval_SS <- 5000

## ---- eval=TRUE---------------------------------------------------------------
fpca_method <- "fpca.sc"

## ---- eval=TRUE---------------------------------------------------------------
sigma2.e <- 0.001

## ---- eval=TRUE---------------------------------------------------------------
fpca_optns  <- list("pve" = 0.95) 

## ---- eval=TRUE---------------------------------------------------------------
mean_diff_add_args  <- list() 

## ---- eval=TRUE---------------------------------------------------------------
nWgrid <- 201

## ---- message=FALSE, results='hide'-------------------------------------------
library(foreach)
required_sample_size <- foreach(i=1:10, .combine='c', .packages = "fPASS") %do% {
        mean_diff_fn <- mean_diff_fn
        PASS_Proj_Test_ufDA(sample_size = NULL, target.power = target.power,
                            sig.level = sig.level, nobs_per_subj = nobs_per_subj,
                            obs.design = obs.design, mean_diff_fnm = "mean_diff_fn",
                            cov.type = cov.type, cov.par = cov.par,
                            sigma2.e = sigma2.e, missing_type = missing_type,
                            missing_percent = missing_percent, eval_SS = eval_SS,
                            alloc.ratio = alloc.ratio, fpca_method = fpca_method,
                            mean_diff_add_args = mean_diff_add_args,
                            fpca_optns = fpca_optns, npc_to_use = NULL,
                            nsim = 1e3)$required_SS
}

## -----------------------------------------------------------------------------
quantile(required_sample_size, probs = c(0.25,0.5,0.75), names = TRUE)

## ---- message = F, results='hide'---------------------------------------------
fpca_method <- "fpca.sc"

library(foreach)
eigenlist <- foreach(i=1:3, .packages = "fPASS") %do% {
  mean_diff_fn <- mean_diff_fn
eigencomp <- Extract_Eigencomp_fDA(nobs_per_subj = nobs_per_subj,
                      obs.design = obs.design, mean_diff_fnm = mean_diff_fnm,
                      cov.type = cov.type, cov.par = cov.par,
                      sigma2.e = sigma2.e, missing_type = missing_type,
                      missing_percent = missing_percent, eval_SS = eval_SS,
                      alloc.ratio = alloc.ratio, fpca_method = fpca_method,
                      mean_diff_add_args = mean_diff_add_args,
                      fpca_optns = fpca_optns,
                      data.driven.scores = FALSE) 
 eigencomp[c("working.grid", "est_eigenfun", "est_eigenval")]
}

## -----------------------------------------------------------------------------
print(lapply(eigenlist, function(x) x$est_eigenval))

## ---- message = F, results='hide'---------------------------------------------

fpca_method <- "face"
# Setting eval_SS to 500, specifically for 
# the face::face.sparse() function
# to ensure that the vignette builds within
# real time. The user must set a much larger value
# of eval_SS e.g. ranging between 2000-5000, to ensure
# enough precision accuracy for the eigenfunctions. 
eval_SS  <- 500 # increase it to 2000 for practical case

library(foreach)
# library(doParallel)
# cl <- makeCluster(detectCores()-1)
# registerDoParallel(cl)
eigenlist <- foreach(i=1:3, .packages = "fPASS") %do% {
  mean_diff_fn <- mean_diff_fn
eigencomp <- Extract_Eigencomp_fDA(nobs_per_subj = nobs_per_subj,
                      obs.design = obs.design, mean_diff_fnm = mean_diff_fnm,
                      cov.type = cov.type, cov.par = cov.par,
                      sigma2.e = sigma2.e, missing_type = missing_type,
                      missing_percent = missing_percent, eval_SS = eval_SS,
                      alloc.ratio = alloc.ratio, fpca_method = fpca_method,
                      mean_diff_add_args = mean_diff_add_args,
                      fpca_optns = fpca_optns,
                      data.driven.scores = FALSE) 
 eigencomp[c("working.grid", "est_eigenfun", "est_eigenval")]
}
# stopCluster(cl)
# matplot(eigencomp$working.grid, eigencomp$est_eigenfun, type="l", 
#         xlab = "timepoints (scaled)", ylab = "eigenfunctions")

## ---- message=F---------------------------------------------------------------
print(lapply(eigenlist, function(x) x$est_eigenval))

## -----------------------------------------------------------------------------
npc_to_use <- 3

## ---- message=F, results='hide'-----------------------------------------------
fpca_method <- "face"

library(foreach)
required_sample_size <- foreach(i=1:5, .combine='c', .packages = "fPASS") %do% {
        mean_diff_fn <- mean_diff_fn
        PASS_Proj_Test_ufDA(sample_size = NULL, target.power = target.power,
                            sig.level = sig.level, nobs_per_subj = nobs_per_subj,
                            obs.design = obs.design, mean_diff_fnm = mean_diff_fnm,
                            cov.type = cov.type, cov.par = cov.par,
                            sigma2.e = sigma2.e, missing_type = missing_type,
                            missing_percent = missing_percent, eval_SS = 500,
                            alloc.ratio = alloc.ratio, fpca_method = fpca_method,
                            mean_diff_add_args = mean_diff_add_args,
                            fpca_optns = fpca_optns, npc_to_use = npc_to_use,
                            nsim = 1e3)$required_SS
}

## -----------------------------------------------------------------------------
quantile(required_sample_size, probs = c(0.25,0.5,0.75), names = TRUE)

## -----------------------------------------------------------------------------
missing_type <- "constant"
missing_percent <- 0.25

## ---- message=FALSE, results='hide'-------------------------------------------
fpca_method <- "fpca.sc"
eval_SS  <- 5000

library(foreach)
required_sample_size <- foreach(i=1:10, .combine='c', .packages = "fPASS") %do% {
        mean_diff_fn <- mean_diff_fn
        PASS_Proj_Test_ufDA(sample_size = NULL, target.power = target.power,
                            sig.level = sig.level, nobs_per_subj = nobs_per_subj,
                            obs.design = obs.design, mean_diff_fnm = "mean_diff_fn",
                            cov.type = cov.type, cov.par = cov.par,
                            sigma2.e = sigma2.e, missing_type = missing_type,
                            missing_percent = missing_percent, eval_SS = eval_SS,
                            alloc.ratio = alloc.ratio, fpca_method = fpca_method,
                            mean_diff_add_args = mean_diff_add_args,
                            fpca_optns = fpca_optns, npc_to_use = 2, nsim = 1e3)$required_SS
}

## ---- message=FALSE-----------------------------------------------------------
quantile(required_sample_size, probs = c(0.25,0.5,0.75), names = TRUE)

