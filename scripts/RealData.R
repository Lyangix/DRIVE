library(DRIVE)
library(LongCART)
library(randomForestSRC)
library(rpart)
load("AnalysisData2.rda")

# The installed examples define the callback contracts used by DRIVE.ML.
source(system.file("examples", "drive_ml_learners.R", package = "DRIVE"))
ml_fitting_rfsrc <- ml_survival_forest
ml_fitting_propensity <- ml_propensity_tree

ind_D = apply(D_status, 1, function(d) {
  if(length(unique(d)) <=1) return(max(stime_new)+1)
  else {
    return(stime_new[which(d != d[1])[1]])
  }
})

Covariates = dat[, c("Female", "Race", "AGE_AT_FIRSTMSICD", "FOLLOWUP_DURA", "DISEASE_DURA",
                     "noteall", "note3m", "MSICDall", "MSICD3m",  "MSCUI3m", "Soluall", "Solu3m",
                     "MRIall", "MRI6m", "HDall", "ERall", "PRIORDMT_DURA",
                     "PRIOR_RELAPSE_12MONS", "PRIOR_RELAPSE_24MONS")]
Covariates = scale(Covariates)

dat_DRIVE = list(T_D_c = round(dat$time, 1),
                event = dat$RELAPSE,
                stime = stime_new,
                W = ind_D,
                Z = Z,
                Covariates = as.matrix(Covariates),
                D_status = D_status)
results = TRTSWE(dat_DRIVE, max_t = max(stime_new), methods = c("DRIVE.joint",
                                                               "DRIVE.ML",
                                                               "ITT", "recensor", "remove",
                                                               "TimeVar"),
                 ml_fitting_propensity = ml_fitting_propensity,
                 ml_fitting_surv = ml_fitting_rfsrc,
                 Control = list(init_parameters = rep(0, 20),
                                seed = rpois(1, lambda = 10*abs(rnorm(1))),
                                nfolds = 796,
                                B = 100)
                 )
results


Covariates = dat[, c("Female", "Race", "AGE_AT_FIRSTMSICD", "FOLLOWUP_DURA", "DISEASE_DURA",
                     "Solu3m",
                     "HDall", "ERall",
                     "PRIOR_RELAPSE_12MONS", "PRIOR_RELAPSE_24MONS")]
Covariates = scale(Covariates)
Covariates2 = dat[, c("Female", "Race", "FOLLOWUP_DURA", "DISEASE_DURA",
                      "note3m",  "MSICD3m",  "MSCUI3m", "Soluall", "Solu3m",
                      "MRI6m", "HDall", "ERall", "PRIORDMT_DURA",
                      "PRIOR_RELAPSE_12MONS", "PRIOR_RELAPSE_24MONS")]
Covariates2 = scale(Covariates2)

dat_DRIVE = list(T_D_c = round(dat$time, 1),
                event = dat$RELAPSE,
                stime = stime_new,
                W = ind_D,
                Z = Z,
                Covariates = as.matrix(Covariates),
                Covariates2 = as.matrix(Covariates2),
                D_status = D_status)
results = TRTSWE(dat_DRIVE, max_t = max(stime_new), methods = c("DRIVE.joint",
                                                               "DRIVE.ML",
                                                               "ITT", "recensor", "remove",
                                                               "TimeVar"),
                 ml_fitting_propensity = ml_fitting_propensity,
                 ml_fitting_surv = ml_fitting_rfsrc,
                 Control = list(init_parameters = rep(0, 11),
                                seed = rpois(1, lambda = 10*abs(rnorm(1))),
                                nfolds = 796,
                                B = 100))
results

mean(apply(D_status[D_status[, 1] == 1, ], 1, function(d) {
  if(any(d!=d[1])) return(TRUE)
  else FALSE
}))



