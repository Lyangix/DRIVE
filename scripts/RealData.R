library(DRIVE)
library(LongCART)
library(randomForestSRC)
library(rpart)
load("AnalysisData2.rda")

ml_fitting_rfsrc = function(data, predictx) {
  ind_control = apply(data$D_status, 1, function(d) return(all(d == 0)))
  p = ncol(data$Covariates)
  tmp_data = data.frame(time = data$time[ind_control], event = data$event[ind_control],
                        data$Covariates[ind_control, , drop = FALSE])
  colnames(tmp_data) = c("time", "event", paste0("X", 1:p))
  colnames(predictx) = paste0("X", 1:p)
  mod = rfsrc(Surv(time, event) ~ ., data = tmp_data, ntime = NULL)
  pre_mod = predict(mod, newdata = predictx)
  pre_mod_stime = pre_mod$time.interest
  pre_mod_chf = pre_mod$chf
  pre_mod_chf = cbind(pre_mod_chf[, 1], t(apply(pre_mod_chf, 1, diff)))
  tmp_chf = matrix(0, nrow = nrow(predictx), ncol = length(data$stime))
  for(j in 1:length(pre_mod_stime)) {
    tmp_chf[, data$stime >= pre_mod_stime[j]] = pre_mod_chf[, j]
  }
  return(tmp_chf)
}


ml_fitting_propensity = function(data, predictx) {
  # p = ncol(data$Covariates)
  # colnames(data) = c("IV", paste0("X", 1:p))
  # colnames(predictx) = c("IV",  paste0("X", 1:p))
  mod = rpart(IV~., data = data, method = "class")
  return(predict(mod, newdata = predictx)[, 2])
}

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

dat_DRIV = list(T_D_c = round(dat$time, 1),
                event = dat$RELAPSE,
                stime = stime_new,
                W = ind_D,
                Z = Z,
                Covariates = as.matrix(Covariates),
                D_status = D_status)
results = TRTSWE(dat_DRIV, max_t = max(stime_new), methods = c("DRIV.s",
                                                               "DRIV.cf.hz.ml.est",
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

dat_DRIV = list(T_D_c = round(dat$time, 1),
                event = dat$RELAPSE,
                stime = stime_new,
                W = ind_D,
                Z = Z,
                Covariates = as.matrix(Covariates),
                Covariates2 = as.matrix(Covariates2),
                D_status = D_status)
results = TRTSWE(dat_DRIV, max_t = max(stime_new), methods = c("DRIV.s",
                                                               "DRIV.cf.hz.ml.est",
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



