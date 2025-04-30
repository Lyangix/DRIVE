require(timereg)



arg_filter <- function(args, fn) {
  args[names(formals(fn))]
}


easy_call = function(fn, args) {
  out = try(do.call(fn, arg_filter(args, fn)), FALSE)
  if (inherits(out, "try-error")) {
    para = names(formals(fn))
    stop("These arguments are needed: ", para[!(para %in% names(args))], "\n")
  } else{
    return(out)
  }
}


expit = function(d) return(1/(1 + exp(-d)))


DataGenerating = function(SimuArg) {
  UseMethod("DataGenerating")
}


DataGenerating.SimuArg.exogenous = function(SimuArg) {
  args = rlang::dots_list(!!!SimuArg$parameters, !!!SimuArg$initials)
  Z_proportion = NULL
  switching_rate_overall = NULL
  switching_rate_from_0 = NULL
  switching_rate_from_1 = NULL
  censoring_rate = NULL
  adcensoring_rate = NULL
  for (kk in 1:SimuArg$initials$nrep) {
    Covariates = easy_call(SimuArg$InitCovariates, rlang::dots_list(!!!args, unmeasured_Confounding = SimuArg$unmeasured_Confounding))
    Z = easy_call(SimuArg$InitAssignment, rlang::dots_list(!!!args, Covariates = Covariates))
    W = easy_call(SimuArg$SwitchingTime, rlang::dots_list(!!!args, Covariates = Covariates, Z = Z))
    T = easy_call(SimuArg$SurvTime, rlang::dots_list(!!!args, Covariates = Covariates, W = W, Z = Z))
    T_D = T$T_D
    T_0 = T$T_0
    C = easy_call(SimuArg$CensoringTime, rlang::dots_list(!!!args, Covariates = Covariates))
    T_D_c = ifelse(T_D <= C, T_D, C)
    T_D_c = ifelse(T_D_c <= args$max_t, T_D_c, args$max_t)
    event = T_D <= C & T_D <= args$max_t
    W_copy = W
    W[W < 0] = Inf
    if(any(T_D <= 0)) warning("Negative time-to-event outcome occurs")
    data = list(Covariates = Covariates,
                Z = Z,
                W = W,
                T_D_c = T_D_c,
                T_D = T_D,
                T_0 = T_0,
                C = C,
                event = event)
    if(SimuArg$Control$json_save){
      if(!dir.exists(paste0(SimuArg$Control$save_path, "DataGenerated"))){
        dir.create(paste0(SimuArg$Control$save_path, "DataGenerated"))
      }
      if(!dir.exists(paste0(SimuArg$Control$save_path, "DataGenerated/", 
                            SimuArg$initials$N, SimuArg$Control$Annotation))){
        dir.create(paste0(SimuArg$Control$save_path, "DataGenerated/", 
                          SimuArg$initials$N, SimuArg$Control$Annotation))
      }
      jsonlite::write_json(data, paste0(SimuArg$Control$save_path, "DataGenerated/", 
                                        SimuArg$initials$N, SimuArg$Control$Annotation,
                                        "/", kk, ".json"))
    }
    T_D_c = ifelse(T_D <= C, T_D, C)
    T_D_c = ifelse(T_D_c <= args$max_t, T_D_c, args$max_t)
    event = T_D <= args$max_t & T_D <= C
    Z_proportion = c(Z_proportion, mean(Z))
    switching_rate_overall = c(switching_rate_overall, mean(T_D_c > W))
    switching_rate_from_0 = c(switching_rate_from_0, mean((T_D_c>W)[Z == 0]))
    switching_rate_from_1 = c(switching_rate_from_1, mean((T_D_c>W)[Z == 1]))
    censoring_rate = c(censoring_rate, mean(1-event))
    adcensoring_rate = c(adcensoring_rate, mean(T_D > args$max_t))
  }
  list(Z_proportion = mean(Z_proportion),
       switching_rate_overall = mean(switching_rate_overall),
       switching_rate_from_0 = mean(switching_rate_from_0),
       switching_rate_from_1 = mean(switching_rate_from_1),
       censoring_rate = mean(censoring_rate),
       adcensoring_rate = mean(adcensoring_rate),
       anyNegetive = sum(W_copy < 0))
}



DataGenerating.SimuArg.endogenous = function(SimuArg) {
  args = rlang::dots_list(!!!SimuArg$parameters, !!!SimuArg$initials)
  Z_proportion = NULL
  switching_rate_overall = NULL
  switching_rate_from_0 = NULL
  switching_rate_from_1 = NULL
  censoring_rate = NULL
  adcensoring_rate = NULL
  for (kk in 1:SimuArg$initials$nrep) {
    Covariates = easy_call(SimuArg$InitCovariates, rlang::dots_list(!!!args, unmeasured_Confounding = SimuArg$unmeasured_Confounding))
    Z = easy_call(SimuArg$InitAssignment, rlang::dots_list(!!!args, Covariates = Covariates))
    T = rexp(args$N)
    W = easy_call(SimuArg$SwitchingTime, rlang::dots_list(!!!args, Covariates = Covariates, Z = Z, T = T))
    T = easy_call(SimuArg$SurvTime, rlang::dots_list(!!!args, Covariates = Covariates, W = W, Z = Z, T = T))
    T_D = T$T_D
    T_0 = T$T_0
    C = easy_call(SimuArg$CensoringTime, rlang::dots_list(!!!args, Covariates = Covariates))
    T_D_c = ifelse(T_D <= C, T_D, C)
    T_D_c = ifelse(T_D_c <= args$max_t, T_D_c, args$max_t)
    event = T_D <= C & T_D <= args$max_t
    W_copy = W
    W[W < 0] = Inf
    if(any(T_D <= 0)) warning("Negative time-to-event outcome occurs")
    data = list(Covariates = Covariates,
                Z = Z,
                W = W,
                T_D_c = T_D_c,
                T_D = T_D,
                T_0 = T_0,
                C = C,
                event = event)
    if(SimuArg$Control$json_save){
      if(!dir.exists(paste0(SimuArg$Control$save_path, "DataGenerated"))){
        dir.create(paste0(SimuArg$Control$save_path, "DataGenerated"))
      }
      if(!dir.exists(paste0(SimuArg$Control$save_path, "DataGenerated/", 
                            SimuArg$initials$N, SimuArg$Control$Annotation))){
        dir.create(paste0(SimuArg$Control$save_path, "DataGenerated/", 
                          SimuArg$initials$N, SimuArg$Control$Annotation))
      }
      jsonlite::write_json(data, paste0(SimuArg$Control$save_path, "DataGenerated/", 
                                        SimuArg$initials$N, SimuArg$Control$Annotation,
                                        "/", kk, ".json"))
    }
    T_D_c = ifelse(T_D <= C, T_D, C)
    T_D_c = ifelse(T_D_c <= args$max_t, T_D_c, args$max_t)
    event = T_D <= args$max_t & T_D <= C
    Z_proportion = c(Z_proportion, mean(Z))
    switching_rate_overall = c(switching_rate_overall, mean(T_D_c > W))
    switching_rate_from_0 = c(switching_rate_from_0, mean((T_D_c>W)[Z == 0]))
    switching_rate_from_1 = c(switching_rate_from_1, mean((T_D_c>W)[Z == 1]))
    censoring_rate = c(censoring_rate, mean(1-event))
    adcensoring_rate = c(adcensoring_rate, mean(T_D > args$max_t))
  }
  list(Z_proportion = mean(Z_proportion),
       switching_rate_overall = mean(switching_rate_overall),
       switching_rate_from_0 = mean(switching_rate_from_0),
       switching_rate_from_1 = mean(switching_rate_from_1),
       censoring_rate = mean(censoring_rate),
       adcensoring_rate = mean(adcensoring_rate),
       anyNegetive = sum(W_copy < 0))
}




DataFitting = function(ModelPar){
  UseMethod("DataFitting")
}


DataFitting.ModelPar.ITT = function(ModelPar){
  # event = ModelPar$dat$T_D < ModelPar$dat$C & ModelPar$dat$T_D < ModelPar$max_t
  event = ModelPar$dat$event
  surv = survival::Surv(ModelPar$dat$T_D_c + runif(ModelPar$N, 0, 0.01),
                        event = event,
                        type = "right")
  mod = ahaz::ahaz(surv, cbind(ModelPar$dat$Z, 
                               ModelPar$dat$Covariates[, ]))
  return(list(Coef = coef(mod),
              Var = diag(vcov(mod))))
}


DataFitting.ModelPar.remove = function(ModelPar){
  # event = ModelPar$dat$T_D < ModelPar$dat$C & ModelPar$dat$T_D < ModelPar$max_t
  event = ModelPar$dat$event
  event_remove = ModelPar$dat$T_D_c < ModelPar$dat$W
  surv = survival::Surv(ModelPar$dat$T_D_c[event_remove] + runif(sum(event_remove), 0, 0.01), 
                        event = event[event_remove], type = "right")
  mod = ahaz::ahaz(surv, cbind(ModelPar$dat$Z[event_remove], 
                               ModelPar$dat$Covariates[event_remove, ]))
  return(list(Coef = coef(mod),
              Var = diag(vcov(mod))))
}



DataFitting.ModelPar.recensor = function(ModelPar){
  # event = ModelPar$dat$T_D < ModelPar$dat$C & ModelPar$dat$T_D < ModelPar$max_t
  event = ModelPar$dat$event
  event_w = event & ModelPar$dat$T_D_c < ModelPar$dat$W
  T_D_c_w = ifelse(ModelPar$dat$T_D_c < ModelPar$dat$W, ModelPar$dat$T_D_c, ModelPar$dat$W)
  surv = survival::Surv(T_D_c_w + runif(ModelPar$N, 0, 0.01), event = event_w, type = "right")
  mod = ahaz::ahaz(surv, cbind(ModelPar$dat$Z, 
                               ModelPar$dat$Covariates))
  return(list(Coef = coef(mod),
              Var = diag(vcov(mod))))
}


DataFitting.ModelPar.TimeVar = function(ModelPar){
  # event = ModelPar$dat$T_D < ModelPar$dat$C & ModelPar$dat$T_D < ModelPar$max_t
  event = ModelPar$dat$event
  event_w = ModelPar$dat$T_D_c > ModelPar$dat$W
  p = ncol(ModelPar$dat$Covariates)
  colnames(ModelPar$dat$Covariates) = paste0("X", 1:p)
  tvdat = data.frame(id = 1:ModelPar$N, treatment = ModelPar$dat$Z, 
                     ModelPar$dat$Covariates, event = event, 
                     start_time = 0, end_time = ModelPar$dat$T_D_c)
  swdat = tvdat[event_w, ]
  tvdat$end_time[event_w] = ModelPar$dat$W[event_w]
  tvdat$event[event_w] = FALSE
  swdat$start_time = ModelPar$dat$W[event_w]
  swdat$treatment = 1 - swdat$treatment
  adat = rbind(tvdat, swdat)
  adat = adat[order(adat$id, adat$start_time), ]
  if(any(adat$end_time == 0)){
    nn = length(adat$end_time[adat$end_time == 0])
    adat$end_time[adat$end_time == 0] = runif(nn, 0, 0.01)
  }
  str_formula_aalen = stringr::str_c(paste0("const(X", 1:(ModelPar$p), ")"), collapse = "+")
  str_formula_aalen = paste0("Surv(start_time, end_time, event) ~ const(treatment) + ",
                             str_formula_aalen)
  mod = timereg::aalen(formula(str_formula_aalen), data = adat, 
                       max.time = ModelPar$max_t, id = adat$id)
  
  return(list(Coef = coef(mod)[, 1],
              Var = coef(mod)[, 2]^2))
}


DataFitting.ModelPar.DRIV = function(ModelPar){
  # event = ModelPar$dat$T_D < ModelPar$dat$C & ModelPar$dat$T_D < ModelPar$max_t
  event = ModelPar$dat$event
  T_D_c = ceiling(ModelPar$dat$T_D_c * ModelPar$Control$grid)/ModelPar$Control$grid
  stime = sort(T_D_c)
  stime = unique(stime)
  k = length(stime)
  D_status = matrix(nrow = ModelPar$N, ncol = k)
  for (i in 1:ModelPar$N) {
    if(T_D_c[i] > ModelPar$dat$W[i]){
      D_status[i, which(stime <= ModelPar$dat$W[i])] = ModelPar$dat$Z[i]
      D_status[i, which(stime > ModelPar$dat$W[i])] = 1 - ModelPar$dat$Z[i]
    } else {
      D_status[i, ] = ModelPar$dat$Z[i]
    }
  }
  if(!("Covariates2" %in% names(ModelPar))){
    args = rlang::dots_list(!!!ModelPar$Control, time = T_D_c, event = event,
                            IV = ModelPar$dat$Z, Covariates = ModelPar$dat$Covariates[, 1:ModelPar$p], 
                            Covariates2 = ModelPar$dat$Covariates[, 1:ModelPar$p], 
                            D_status = D_status, stime = stime)
  }else{
    args = rlang::dots_list(!!!ModelPar$Control, time = T_D_c, event = event,
                            IV = ModelPar$dat$Z, Covariates = ModelPar$dat$Covariates[, 1:ModelPar$p], 
                            Covariates2 = ModelPar$Covariates2, 
                            D_status = D_status, stime = stime)
  }
  mod = easy_call(integral_est_cpp, args)
  return(list(Coef = mod$x,
              Var = mod$var,
              Convergence = mod$Convergence))
}


DataFitting.ModelPar.DRIV.s = function(ModelPar){
  # event = ModelPar$dat$T_D < ModelPar$dat$C & ModelPar$dat$T_D < ModelPar$max_t
  ModelPar$Control = rlang::dots_list(!!!ModelPar$Control, 
                                      init_parameters = rep(0, ncol(ModelPar$dat$Covariates)+1),
                                      .homonyms = "first")
  event = ModelPar$dat$event
  # T_D_c = ceiling(ModelPar$dat$T_D_c * ModelPar$Control$grid)/ModelPar$Control$grid
  T_D_c = ModelPar$dat$T_D_c
  if(is.null(ModelPar$dat$stime)){
    stime = sort(T_D_c)
    stime = unique(stime)
  } else {
    stime = ModelPar$dat$stime
  }
  k = length(stime)
  if(is.null(ModelPar$dat$D_status)){
    D_status = matrix(nrow = ModelPar$N, ncol = k)
    for (i in 1:ModelPar$N) {
      if(T_D_c[i] > ModelPar$dat$W[i]){
        D_status[i, which(stime <= ModelPar$dat$W[i])] = ModelPar$dat$Z[i]
        D_status[i, which(stime > ModelPar$dat$W[i])] = 1 - ModelPar$dat$Z[i]
      } else {
        D_status[i, ] = ModelPar$dat$Z[i]
      }
    }
  } else {
    D_status = ModelPar$dat$D_status
  }
  
  if(!("Covariates2" %in% names(ModelPar))){
    args = rlang::dots_list(!!!ModelPar$Control, time = T_D_c, event = event,
                            IV = ModelPar$dat$Z, Covariates = ModelPar$dat$Covariates, 
                            Covariates2 = ModelPar$dat$Covariates, 
                            D_status = D_status, stime = stime)
  }else{
    args = rlang::dots_list(!!!ModelPar$Control, time = T_D_c, event = event,
                            IV = ModelPar$dat$Z, Covariates = ModelPar$dat$Covariates, 
                            Covariates2 = ModelPar$Covariates2, 
                            D_status = D_status, stime = stime)
  }
  mod = easy_call(integral2_est_cpp, args)
  return(list(Coef = mod$x,
              Var = mod$var,
              Convergence = mod$Convergence,
              dLam = mod$dLam))
}

DataFitting.ModelPar.DRIV.s.nleqslv = function(ModelPar){
  # event = ModelPar$dat$T_D < ModelPar$dat$C & ModelPar$dat$T_D < ModelPar$max_t
  event = ModelPar$dat$event
  # T_D_c = ceiling(ModelPar$dat$T_D_c * ModelPar$Control$grid)/ModelPar$Control$grid
  T_D_c = ModelPar$dat$T_D_c
  if(is.null(ModelPar$dat$stime)){
    stime = sort(T_D_c)
    stime = unique(stime)
  } else {
    stime = ModelPar$dat$stime
  }
  k = length(stime)
  if(is.null(ModelPar$dat$D_status)){
    D_status = matrix(nrow = ModelPar$N, ncol = k)
    for (i in 1:ModelPar$N) {
      if(T_D_c[i] > ModelPar$dat$W[i]){
        D_status[i, which(stime <= ModelPar$dat$W[i])] = ModelPar$dat$Z[i]
        D_status[i, which(stime > ModelPar$dat$W[i])] = 1 - ModelPar$dat$Z[i]
      } else {
        D_status[i, ] = ModelPar$dat$Z[i]
      }
    }
  } else {
    D_status = ModelPar$dat$D_status
  }
  if(!("Covariates2" %in% names(ModelPar))){
    args = rlang::dots_list(!!!ModelPar$Control, time = T_D_c, event = event,
                            IV = ModelPar$dat$Z, Covariates = ModelPar$dat$Covariates, 
                            Covariates2 = ModelPar$dat$Covariates, 
                            D_status = D_status, stime = stime)
  }else{
    args = rlang::dots_list(!!!ModelPar$Control, time = T_D_c, event = event,
                            IV = ModelPar$dat$Z, Covariates = ModelPar$dat$Covariates, 
                            Covariates2 = ModelPar$Covariates2, 
                            D_status = D_status, stime = stime)
  }
  mod = easy_call(integral2_cpp, args)
  print(mod$iter)
  print(mod$x)
  print(mod$fvec)
  print(mod$termcd)
  print(mod$jac)
  return(list(Coef = mod$x,
              Var = mod$var,
              Convergence = mod$Convergence))
}



DataFitting.ModelPar.DRIV.cf.hz = function(ModelPar){
  # event = ModelPar$dat$T_D < ModelPar$dat$C & ModelPar$dat$T_D < ModelPar$max_t
  event = ModelPar$dat$event
  # T_D_c = ceiling(ModelPar$dat$T_D_c * ModelPar$Control$grid)/ModelPar$Control$grid
  T_D_c = ModelPar$dat$T_D_c
  stime = sort(T_D_c)
  stime = unique(stime)
  k = length(stime)
  D_status = matrix(nrow = ModelPar$N, ncol = k)
  for (i in 1:ModelPar$N) {
    if(T_D_c[i] > ModelPar$dat$W[i]){
      D_status[i, which(stime <= ModelPar$dat$W[i])] = ModelPar$dat$Z[i]
      D_status[i, which(stime > ModelPar$dat$W[i])] = 1 - ModelPar$dat$Z[i]
    } else {
      D_status[i, ] = ModelPar$dat$Z[i]
    }
  }
  
  if(is.null(ModelPar$nfolds)) ModelPar$nfolds = 10
  if(is.null(ModelPar$seed)) ModelPar$seed = 5884419
  
  if(!("Covariates2" %in% names(ModelPar))){
    args = rlang::dots_list(!!!ModelPar$Control, time = T_D_c, event = event,
                            IV = ModelPar$dat$Z, Covariates = ModelPar$dat$Covariates, 
                            Covariates2 = ModelPar$dat$Covariates, 
                            D_status = D_status, stime = stime, nfolds = ModelPar$nfolds, seed = ModelPar$seed)
  }else{
    args = rlang::dots_list(!!!ModelPar$Control, time = T_D_c, event = event,
                            IV = ModelPar$dat$Z, Covariates = ModelPar$dat$Covariates, 
                            Covariates2 = ModelPar$Covariates2, 
                            D_status = D_status, stime = stime, nfolds = ModelPar$nfolds, seed = ModelPar$seed)
  }
  mod = easy_call(cfhaz_integral_cpp, args)
  return(list(Coef = mod$x,
              Var = mod$var,
              Convergence = mod$Convergence))
}



DataFitting.ModelPar.DRIV.cf.hz.est = function(ModelPar){
  # event = ModelPar$dat$T_D < ModelPar$dat$C & ModelPar$dat$T_D < ModelPar$max_t
  event = ModelPar$dat$event
  # T_D_c = ceiling(ModelPar$dat$T_D_c * ModelPar$Control$grid)/ModelPar$Control$grid
  T_D_c = ModelPar$dat$T_D_c
  if(is.null(ModelPar$dat$stime)){
    stime = sort(T_D_c)
    stime = unique(stime)
  } else {
    stime = ModelPar$dat$stime
  }
  k = length(stime)
  if(is.null(ModelPar$dat$D_status)){
    D_status = matrix(nrow = ModelPar$N, ncol = k)
    for (i in 1:ModelPar$N) {
      if(T_D_c[i] > ModelPar$dat$W[i]){
        D_status[i, which(stime <= ModelPar$dat$W[i])] = ModelPar$dat$Z[i]
        D_status[i, which(stime > ModelPar$dat$W[i])] = 1 - ModelPar$dat$Z[i]
      } else {
        D_status[i, ] = ModelPar$dat$Z[i]
      }
    }
  } else {
    D_status = ModelPar$dat$D_status
  }
  
  if(is.null(ModelPar$nfolds)) ModelPar$nfolds = 10
  if(is.null(ModelPar$seed)) ModelPar$seed = 5884419
  
  if(!("Covariates2" %in% names(ModelPar))){
    args = rlang::dots_list(!!!ModelPar$Control, time = T_D_c, event = event,
                            IV = ModelPar$dat$Z, Covariates = ModelPar$dat$Covariates, 
                            Covariates2 = ModelPar$dat$Covariates, 
                            D_status = D_status, stime = stime, nfolds = ModelPar$nfolds, seed = ModelPar$seed, .homonyms = "first")
  }else{
    args = rlang::dots_list(!!!ModelPar$Control, time = T_D_c, event = event,
                            IV = ModelPar$dat$Z, Covariates = ModelPar$dat$Covariates, 
                            Covariates2 = ModelPar$Covariates2, 
                            D_status = D_status, stime = stime, nfolds = ModelPar$nfolds, seed = ModelPar$seed, .homonyms = "first")
  }
  mod = easy_call(cfhaz_integral_est_cpp, args)
  return(list(Coef = mod$x,
              Var = mod$var,
              Convergence = mod$Convergence))
}

DataFitting.ModelPar.DRIV.cf.hz.est.ave = function(ModelPar){
  # event = ModelPar$dat$T_D < ModelPar$dat$C & ModelPar$dat$T_D < ModelPar$max_t
  event = ModelPar$dat$event
  # T_D_c = ceiling(ModelPar$dat$T_D_c * ModelPar$Control$grid)/ModelPar$Control$grid
  T_D_c = ModelPar$dat$T_D_c
  if(is.null(ModelPar$dat$stime)){
    stime = sort(T_D_c)
    stime = unique(stime)
  } else {
    stime = ModelPar$dat$stime
  }
  k = length(stime)
  if(is.null(ModelPar$dat$D_status)){
    D_status = matrix(nrow = ModelPar$N, ncol = k)
    for (i in 1:ModelPar$N) {
      if(T_D_c[i] > ModelPar$dat$W[i]){
        D_status[i, which(stime <= ModelPar$dat$W[i])] = ModelPar$dat$Z[i]
        D_status[i, which(stime > ModelPar$dat$W[i])] = 1 - ModelPar$dat$Z[i]
      } else {
        D_status[i, ] = ModelPar$dat$Z[i]
      }
    }
  } else {
    D_status = ModelPar$dat$D_status
  }
  
  if(is.null(ModelPar$nfolds)) ModelPar$nfolds = 10
  if(is.null(ModelPar$seed)) ModelPar$seed = 5884419
  
  if(!("Covariates2" %in% names(ModelPar))){
    args = rlang::dots_list(seed = NULL, !!!ModelPar$Control, time = T_D_c, event = event,
                            IV = ModelPar$dat$Z, Covariates = ModelPar$dat$Covariates, 
                            Covariates2 = ModelPar$dat$Covariates, 
                            D_status = D_status, stime = stime, nfolds = ModelPar$nfolds,  .homonyms = "first")
  }else{
    args = rlang::dots_list(seed = NULL, !!!ModelPar$Control, time = T_D_c, event = event,
                            IV = ModelPar$dat$Z, Covariates = ModelPar$dat$Covariates, 
                            Covariates2 = ModelPar$Covariates2, 
                            D_status = D_status, stime = stime, nfolds = ModelPar$nfolds, .homonyms = "first")
  }
  Coef = rep(0, ModelPar$Control$B)
  Var = rep(0, ModelPar$Control$B)
  Convergence = rep(FALSE, ModelPar$Control$B)
  for (i in 1:ModelPar$Control$B) {
    mod = easy_call(cfhaz_integral_est_cpp, args)
    Coef[i] = mod$x
    Var[i] = mod$var
    Convergence[i] = mod$Convergence
    cat("rep ", i, "\t", mod$x, "\t", mod$var, "\t", mod$Convergence, "\n")
  }
  return(list(Coef = mean(Coef),
              Var = mean(Var),
              Convergence = all(Convergence)))
}


DataFitting.ModelPar.DRIV.cf.hz.ml.est = function(ModelPar){
  # event = ModelPar$dat$T_D < ModelPar$dat$C & ModelPar$dat$T_D < ModelPar$max_t
  event = ModelPar$dat$event
  # T_D_c = ceiling(ModelPar$dat$T_D_c * ModelPar$Control$grid)/ModelPar$Control$grid
  T_D_c = ModelPar$dat$T_D_c
  if(is.null(ModelPar$dat$stime)){
    stime = sort(T_D_c)
    stime = unique(stime)
  } else {
    stime = ModelPar$dat$stime
  }
  k = length(stime)
  if(is.null(ModelPar$dat$D_status)){
    D_status = matrix(nrow = ModelPar$N, ncol = k)
    for (i in 1:ModelPar$N) {
      if(T_D_c[i] > ModelPar$dat$W[i]){
        D_status[i, which(stime <= ModelPar$dat$W[i])] = ModelPar$dat$Z[i]
        D_status[i, which(stime > ModelPar$dat$W[i])] = 1 - ModelPar$dat$Z[i]
      } else {
        D_status[i, ] = ModelPar$dat$Z[i]
      }
    }
  } else {
    D_status = ModelPar$dat$D_status
  }
  
  
  if(is.null(ModelPar$ml_fitting_surv)) stop("ml_fitting_surv is not specified")
  if(is.null(ModelPar$ml_fitting_propensity)) stop("ml_fitting_propensity is not specified")
  if(is.null(ModelPar$nfolds)) ModelPar$nfolds = 10
  if(is.null(ModelPar$seed)) ModelPar$seed = 5884419
  
  if(!("Covariates2" %in% names(ModelPar))){
    args = rlang::dots_list(!!!ModelPar$Control, time = T_D_c, event = event,
                            IV = ModelPar$dat$Z, Covariates = ModelPar$dat$Covariates,
                            ml_fitting_surv = ModelPar$ml_fitting_surv,
                            ml_fitting_propensity = ModelPar$ml_fitting_propensity,
                            Covariates2 = ModelPar$dat$Covariates, 
                            D_status = D_status, stime = stime, nfolds = ModelPar$nfolds, seed = ModelPar$seed, .homonyms = "first")
  }else{
    args = rlang::dots_list(!!!ModelPar$Control, time = T_D_c, event = event,
                            IV = ModelPar$dat$Z, Covariates = ModelPar$dat$Covariates, 
                            ml_fitting_surv = ModelPar$ml_fitting_surv,
                            ml_fitting_propensity = ModelPar$ml_fitting_propensity,
                            Covariates2 = ModelPar$Covariates2, 
                            D_status = D_status, stime = stime, nfolds = ModelPar$nfolds, seed = ModelPar$seed, .homonyms = "first")
  }
  mod = easy_call(cfhaz_ml_integral_est_cpp, args)
  return(list(Coef = mod$x,
              Var = mod$var,
              Convergence = mod$Convergence))
}

DataFitting.ModelPar.DRIV.cf.hz.ml.est.ave = function(ModelPar){
  # event = ModelPar$dat$T_D < ModelPar$dat$C & ModelPar$dat$T_D < ModelPar$max_t
  event = ModelPar$dat$event
  # T_D_c = ceiling(ModelPar$dat$T_D_c * ModelPar$Control$grid)/ModelPar$Control$grid
  T_D_c = ModelPar$dat$T_D_c
  if(is.null(ModelPar$dat$stime)){
    stime = sort(T_D_c)
    stime = unique(stime)
  } else {
    stime = ModelPar$dat$stime
  }
  k = length(stime)
  if(is.null(ModelPar$dat$D_status)){
    D_status = matrix(nrow = ModelPar$N, ncol = k)
    for (i in 1:ModelPar$N) {
      if(T_D_c[i] > ModelPar$dat$W[i]){
        D_status[i, which(stime <= ModelPar$dat$W[i])] = ModelPar$dat$Z[i]
        D_status[i, which(stime > ModelPar$dat$W[i])] = 1 - ModelPar$dat$Z[i]
      } else {
        D_status[i, ] = ModelPar$dat$Z[i]
      }
    }
  } else {
    D_status = ModelPar$dat$D_status
  }
  
  
  if(is.null(ModelPar$ml_fitting_surv)) stop("ml_fitting_surv is not specified")
  if(is.null(ModelPar$ml_fitting_propensity)) stop("ml_fitting_propensity is not specified")
  if(is.null(ModelPar$nfolds)) ModelPar$nfolds = 10
  # if(is.null(ModelPar$seed)) ModelPar$seed = 5884419
  
  if(!("Covariates2" %in% names(ModelPar))){
    args = rlang::dots_list(seed = NULL, !!!ModelPar$Control, time = T_D_c, event = event,
                            IV = ModelPar$dat$Z, Covariates = ModelPar$dat$Covariates,
                            ml_fitting_surv = ModelPar$ml_fitting_surv,
                            ml_fitting_propensity = ModelPar$ml_fitting_propensity,
                            Covariates2 = ModelPar$dat$Covariates, 
                            D_status = D_status, stime = stime, nfolds = ModelPar$nfolds,  .homonyms = "first")
  }else{
    args = rlang::dots_list(seed = NULL, !!!ModelPar$Control, time = T_D_c, event = event,
                            IV = ModelPar$dat$Z, Covariates = ModelPar$dat$Covariates, 
                            ml_fitting_surv = ModelPar$ml_fitting_surv,
                            ml_fitting_propensity = ModelPar$ml_fitting_propensity,
                            Covariates2 = ModelPar$Covariates2, 
                            D_status = D_status, stime = stime, nfolds = ModelPar$nfolds, .homonyms = "first")
  }
  Coef = rep(0, ModelPar$Control$B)
  Var = rep(0, ModelPar$Control$B)
  Convergence = rep(FALSE, ModelPar$Control$B)
  for (i in 1:ModelPar$Control$B) {
    mod = easy_call(cfhaz_ml_integral_est_cpp_strcf, args)
    Coef[i] = mod$x
    Var[i] = mod$var
    Convergence[i] = mod$Convergence
    cat("rep ", i, "\t", mod$x, "\t", mod$var, "\t", mod$Convergence, "\n")
  }
  return(list(Coef = mean(Coef),
              Var = mean(Var),
              Convergence = all(Convergence)))
}

DataFitting.ModelPar.DRIV.cf.hz.ml.est.rateCal = function(ModelPar){
  # event = ModelPar$dat$T_D < ModelPar$dat$C & ModelPar$dat$T_D < ModelPar$max_t
  event = ModelPar$dat$event
  # T_D_c = ceiling(ModelPar$dat$T_D_c * ModelPar$Control$grid)/ModelPar$Control$grid
  T_D_c = ModelPar$dat$T_D_c
  if(is.null(ModelPar$dat$stime)){
    stime = sort(T_D_c)
    stime = unique(stime)
  } else {
    stime = ModelPar$dat$stime
  }
  k = length(stime)
  if(is.null(ModelPar$dat$D_status)){
    D_status = matrix(nrow = ModelPar$N, ncol = k)
    for (i in 1:ModelPar$N) {
      if(T_D_c[i] > ModelPar$dat$W[i]){
        D_status[i, which(stime <= ModelPar$dat$W[i])] = ModelPar$dat$Z[i]
        D_status[i, which(stime > ModelPar$dat$W[i])] = 1 - ModelPar$dat$Z[i]
      } else {
        D_status[i, ] = ModelPar$dat$Z[i]
      }
    }
  } else {
    D_status = ModelPar$dat$D_status
  }
  
  if(is.null(ModelPar$ml_fitting_surv)) stop("ml_fitting_surv is not specified")
  if(is.null(ModelPar$ml_fitting_propensity)) stop("ml_fitting_propensity is not specified")
  if(is.null(ModelPar$nfolds)) ModelPar$nfolds = 10
  if(is.null(ModelPar$seed)) ModelPar$seed = 5884419
  
  if(!("Covariates2" %in% names(ModelPar))){
    args = rlang::dots_list(!!!ModelPar$Control, time = T_D_c, event = event,
                            IV = ModelPar$dat$Z, Covariates = ModelPar$dat$Covariates,
                            ml_fitting_surv = ModelPar$ml_fitting_surv,
                            ml_fitting_propensity = ModelPar$ml_fitting_propensity,
                            ml_fitting_surv_true = ModelPar$ml_fitting_surv_true,
                            ml_fitting_propensity_true = ModelPar$ml_fitting_propensity_true,
                            rate_propensity = ModelPar$rate_propensity , 
                            rate_hazard = ModelPar$rate_hazard, 
                            target_biases_hazard = ModelPar$target_biases_hazard,
                            target_biases_prop = ModelPar$target_biases_prop, learning_rate = ModelPar$learning_rate,
                            true_theta = ModelPar$true_theta,
                            Covariates2 = ModelPar$dat$Covariates, 
                            D_status = D_status, stime = stime, nfolds = ModelPar$nfolds, seed = ModelPar$seed)
  }else{
    args = rlang::dots_list(!!!ModelPar$Control, time = T_D_c, event = event,
                            IV = ModelPar$dat$Z, Covariates = ModelPar$dat$Covariates, 
                            ml_fitting_surv = ModelPar$ml_fitting_surv,
                            ml_fitting_propensity = ModelPar$ml_fitting_propensity,
                            ml_fitting_surv_true = ModelPar$ml_fitting_surv_true,
                            ml_fitting_propensity_true = ModelPar$ml_fitting_propensity_true,
                            rate_propensity = ModelPar$rate_propensity , 
                            rate_hazard = ModelPar$rate_hazard, 
                            target_biases_hazard = ModelPar$target_biases_hazard,
                            target_biases_prop = ModelPar$target_biases_prop, learning_rate = ModelPar$learning_rate,
                            true_theta = ModelPar$true_theta,
                            Covariates2 = ModelPar$Covariates2, 
                            D_status = D_status, stime = stime, nfolds = ModelPar$nfolds, seed = ModelPar$seed)
  }
  # print(names(args))
  # print(names(formals(cfhaz_ml_integral_est_cpp_rateCal)))
  mod = easy_call(cfhaz_ml_integral_est_cpp_rateCal, args)
  # Coef = NULL
  # Var = NULL
  # Convergence = NULL
  # mse_propensity = NULL
  # mse_surv = NULL
  # for(i in 1:(ncol(ModelPar$rate_hazard)*nrow(ModelPar$rate_propensity))) {
  #   Coef = c(Coef, mod[[i]]$x)
  #   Var = c(Var, mod[[i]]$var)
  #   Convergence = c(Convergence, mod[[i]]$Convergence)
  #   mse_propensity = c(mse_propensity, mod[[i]]$mse_propensity)
  #   mse_surv = c(mse_surv, mod[[i]]$mse_surv)
  # }
  Coef = sapply(mod, function(x) x$x)
  Var = sapply(mod, function(x) x$var)
  Convergence = sapply(mod, function(x) x$Convergence)
  mse_propensity = sapply(mod, function(x) x$mse_propensity)
  mse_surv = sapply(mod, function(x) x$mse_surv)
  rate_propensity = sapply(mod, function(x) x$rate_propensity)
  rate_hazard = sapply(mod, function(x) x$rate_hazard)
  
  return(list(Coef = Coef,
              Var = Var,
              Convergence = Convergence,
              mse_propensity = mse_propensity,
              mse_surv = mse_surv,
              rate_propensity= rate_propensity,
              rate_hazard = rate_hazard))
}


DataFitting.ModelPar.remove.control = function(ModelPar){
  # event = ModelPar$dat$T_D < ModelPar$dat$C & ModelPar$dat$T_D < ModelPar$max_t
  event = ModelPar$dat$event
  event_remove = ModelPar$dat$T_D_c < ModelPar$dat$W
  event_remove = event_remove & (ModelPar$dat$Z == 0)
  surv = survival::Surv(ModelPar$dat$T_D_c[event_remove] + runif(sum(event_remove), 0, 0.01), 
                        event = event[event_remove], type = "right")
  mod = ahaz::ahaz(surv, ModelPar$dat$Covariates[event_remove, ])
  return(list(Coef = c(0, coef(mod)),
              Var = c(0,diag(vcov(mod)))))
}


DataFitting.ModelPar.counterfactual = function(ModelPar){
  # event = ModelPar$dat$T_D < ModelPar$dat$C & ModelPar$dat$T_D < ModelPar$max_t
  event = rep(1, length(ModelPar$dat$event))
  surv = survival::Surv(ModelPar$dat$T_0 + runif(ModelPar$N, 0, 0.0001), 
                        event = event, type = "right")
  mod = ahaz::ahaz(surv, ModelPar$dat$Covariates[, ])
  return(list(Coef = c(0, coef(mod)),
              Var = c(0,diag(vcov(mod)))))
}


DataFitting.ModelPar.customized = function(ModelPar) {
  # let k = length(ModelPar$stime)
  # ModelPar$PredictedCovariates should be n $\times$ k
  # ModelPar$PredictedBaselineHaz should be n $\times$ k
  # ModelPar$PredictedHaz should be n $\times$ k
  # ModelPar$PredictedPropensityScore should be n $\times$ 1
  
  ModelPar$Control = rlang::dots_list(!!!ModelPar$Control, 
                                      init_parameters = 0.1,
                                      .homonyms = "first")
  event = ModelPar$dat$event
  T_D_c = ModelPar$dat$T_D_c
  stime = sort(T_D_c)
  stime = unique(stime)
  k = length(stime)
  D_status = matrix(nrow = ModelPar$N, ncol = k)
  for (i in 1:ModelPar$N) {
    if(T_D_c[i] > ModelPar$dat$W[i]){
      D_status[i, which(stime <= ModelPar$dat$W[i])] = ModelPar$dat$Z[i]
      D_status[i, which(stime > ModelPar$dat$W[i])] = 1 - ModelPar$dat$Z[i]
    } else {
      D_status[i, ] = ModelPar$dat$Z[i]
    }
  }
  
  if(is.null(ModelPar$PredictedCovariates)){
    if(is.null(ModelPar$PredictedHaz)) {
      stop("Either 'PredictedCovariates' or 'PredictedHaz' is needed")
    } else {
      if(ncol(ModelPar$PredictedHaz) != k){
        stop("'PredictedHaz' should be calculated at all time points")
      } else {
        ConfoundingPart = ModelPar$PredictedHaz
      }
    } 
  } else {
    if (ncol(ModelPar$PredictedCovariates) != k) {
      stop("'PredictedCovariates' should be calculated at all time points")
    } else {
      if(!is.null(ModelPar$PredictedHaz)){
        if(ncol(ModelPar$PredictedHaz) != k){
          stop("'PredictedHaz' should be calculated at all time points")
        } else {
          ConfoundingPart = ModelPar$PredictedHaz
        }
      } else {
        if(!is.null(ModelPar$PredictedBaselineHaz)) {
          if(ncol(ModelPar$PredictedBaselineHaz) != k) {
            ConfoundingPart = ModelPar$PredictedCovariates
            warning("'PredictedBaselineHaz' should be calculated at all time points")
          } else {
            ConfoundingPart = ModelPar$PredictedCovariates + ModelPar$PredictedBaselineHaz
          }
        } else {
          ConfoundingPart = ModelPar$PredictedCovariates
        }
      }
    }
  }
  
  if(is.null(ModelPar$PredictedPropensityScore)) {
    if(is.null(ModelPar$dat$Covariates)) {
      stop("Either 'PredictedPropensityScore' or 'Covariates' attr in 'dat' is needed")
    } else {
      mod = glm(ModelPar$dat$Z ~ ModelPar$dat$Covariates, family = binomial(link = "logit"))
      ModelPar$PredictedPropensityScore = expit(predict(mod))
    }
  }
  
  args = rlang::dots_list(!!!ModelPar$Control, time = T_D_c, event = event,
                          IV = ModelPar$dat$Z,
                          IV_c = ModelPar$dat$Z - ModelPar$PredictedPropensityScore, 
                          ConfoundingPart = ConfoundingPart, 
                          D_status = D_status, stime = stime)
  mod = easy_call(integral_customized_est_cpp, args)
  return(list(Coef = mod$x,
              Var = mod$var,
              Convergence = mod$Convergence))
}


print.SimuResults = function(SimuResults, ...){
  results = rlang::dots_list(!!!SimuResults, !!!list(...), 
                             Comp_parameters = rep(0, SimuResults$initials$p+1),
                             .homonyms = "first")
  cat("Simulation Results for ")
  cat(results$methods, ":\n")
  tb = NULL
  tb2 = NULL
  tb3 = NULL
  for (j in results$methods) {
    
    tb = rbind(tb, apply(results$SimuResults[[j]]$Coef, 1, mean) - results$Comp_parameters)
    tb2 = rbind(tb2, apply(results$SimuResults[[j]]$Coef, 1, sd))
    # browser()
    if(j %in% c("DRIV.s", "DRIV.cf.hz.ml.est", "DRIV.cf.hz.est")){
      tb3 = rbind(tb3, c(mean(sqrt(results$SimuResults[[j]]$Var)), rep(0, nrow(results$SimuResults[[j]]$Coef) - 1)))
    } else {
      tb3 = rbind(tb3, apply(sqrt(results$SimuResults[[j]]$Var), 1, mean))
    }
  }
  rownames(tb) = results$methods
  colnames(tb) = c("theta", paste0("alpha", 1:(SimuResults$initials$p)))
  rownames(tb2) = results$methods
  colnames(tb2) = c("theta", paste0("alpha", 1:(SimuResults$initials$p)))
  rownames(tb3) = results$methods
  colnames(tb3) = c("theta", paste0("alpha", 1:(SimuResults$initials$p)))
  cat("\t Mean bias or sampling mean: ", "\n")
  print.default(round(tb, 4), print.gap = 2L)
  cat("\n")
  cat("\t Standard deviation: ", "\n")
  print.default(round(tb2, 4), print.gap = 2L)
  cat("\n")
  cat("\t Mean standard error: ", "\n")
  print.default(round(tb3, 4), print.gap = 2L)
  cat("\n")
}

print.SimuResults = function(SimuResults, ...){
  results = rlang::dots_list(!!!SimuResults, !!!list(...), 
                             Comp_parameters = rep(0, SimuResults$initials$p+1),
                             .homonyms = "first")
  cat("Simulation Results for ")
  cat(results$methods, ":\n")
  tb = NULL
  tb2 = NULL
  tb3 = NULL
  for (j in results$methods) {
    
    tb = rbind(tb, apply(results$SimuResults[[j]]$Coef, 1, mean) - results$Comp_parameters)
    tb2 = rbind(tb2, apply(results$SimuResults[[j]]$Coef, 1, sd))
    # browser()
    if(j %in% c("DRIV.s", "DRIV.cf.hz.ml.est", "DRIV.cf.hz.est")){
      tb3 = rbind(tb3, c(mean(sqrt(results$SimuResults[[j]]$Var)), rep(0, nrow(results$SimuResults[[j]]$Coef) - 1)))
    } else {
      tb3 = rbind(tb3, apply(sqrt(results$SimuResults[[j]]$Var), 1, mean))
    }
  }
  rownames(tb) = results$methods
  colnames(tb) = c("theta", paste0("alpha", 1:(SimuResults$initials$p)))
  rownames(tb2) = results$methods
  colnames(tb2) = c("theta", paste0("alpha", 1:(SimuResults$initials$p)))
  rownames(tb3) = results$methods
  colnames(tb3) = c("theta", paste0("alpha", 1:(SimuResults$initials$p)))
  cat("\t Mean bias or sampling mean: ", "\n")
  print.default(round(tb, 4), print.gap = 2L)
  cat("\n")
  cat("\t Standard deviation: ", "\n")
  print.default(round(tb2, 4), print.gap = 2L)
  cat("\n")
  cat("\t Mean standard error: ", "\n")
  print.default(round(tb3, 4), print.gap = 2L)
  cat("\n")
}

print.TRTSWE = function(Results, all = FALSE, ...){
  results = rlang::dots_list(!!!Results, !!!list(...),
                             .homonyms = "first")
  
  cat("Results for ")
  cat(names(Results), ":\n")
  tb = NULL
  # tb2 = NULL
  tb3 = NULL
  tb4 = NULL
  p = 0
  for(j in names(Results)){
    if(length(results[[j]]$Estim$Coef)>p) p= length(results[[j]]$Estim$Coef)
  }
  
  
  for (j in names(Results)) {
    
    
    # browser()
    if(j %in% c("DRIV.s", "DRIV.cf.hz.ml.est", "DRIV.cf.hz.est")){
      if(j %in% "DRIV.s") {
        tb = rbind(tb, as.vector(results[[j]]$Estim$Coef))
      } else {
        tb = rbind(tb, c(results[[j]]$Estim$Coef, rep(NA, p-1)))
      }
      # tb2 = rbind(tb2, c(sd(results[[j]]$Estim$Coef), rep(NA, p-1)))
      tb3 = rbind(tb3, c(sqrt(results[[j]]$Estim$Var), rep(NA, p-1)))
      
    } else {
      tb = rbind(tb, results[[j]]$Estim$Coef)
      # tb2 = rbind(tb2, results[[j]]$Estim$Coef)
      tb3 = rbind(tb3, sqrt(results[[j]]$Estim$Var))
    }
  }
  rownames(tb) = names(Results)
  colnames(tb) = c("theta", paste0("alpha", 1:(p-1)))
  # rownames(tb2) = results$methods
  # colnames(tb2) = c("theta", paste0("alpha", 1:(SimuResults$initials$p)))
  rownames(tb3) = names(Results)
  colnames(tb3) = c("theta", paste0("alpha", 1:(p-1)))
  
  if(all){
    cat("\t Coef: ", "\n")
    print.default(round(tb, 4), print.gap = 2L)
    cat("\n")
    # cat("\t Standard deviation: ", "\n")
    # print.default(round(tb2, 4), print.gap = 2L)
    # cat("\n")
    cat("\t Standard error: ", "\n")
    print.default(round(tb3, 4), print.gap = 2L)
    cat("\n")
    tb4 = 2*(1- pnorm(abs(tb/tb3)))
    cat("\t P-value: ", "\n")
    print.default(round(tb4, 4), print.gap = 2L)
  }else {
    tb4 = 2*(1- pnorm(abs(tb/tb3)))
    tb = cbind(tb[, 1], tb3[, 1], tb4[, 1])
    colnames(tb) = c("theta", "SE", "Pval")
    print.default(round(tb, 4), print.gap = 2L)
  }
  
}


