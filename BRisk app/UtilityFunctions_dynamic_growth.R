# Define dBaranyi model 
dBaranyi <- function(time, state, pars, env_func, sec_models) {
  pars <- as.list(pars)
  state <- as.list(state)
  alpha <- state$Q/(1 + state$Q)
  beta <- 1 - state$N/pars$Nmax
  gamma <- calculate_gammas(time, env_func, sec_models)
  mu <- pars$mu_opt*prod(gamma)
  mu <- mu*log(10)  # Convert to ln CFU/[t]
  dN <- alpha * mu * beta * state$N 
  dQ <- mu*state$Q
  list(c(dQ = dQ,
         dN = dN))
}

# Define secondary models 
CPM_model <- function(x, xmin, xopt, xmax, n) {
  num <- (x-xmax)*(x-xmin)^n
  den <- (xopt-xmin)^(n-1)*( (xopt-xmin)*(x-xopt) - (xopt-xmax)*((n-1)*xopt + xmin-n*x) )
  gamma <- num/den
  gamma[x < xmin] <- 0
  gamma[x > xmax] <- 0
  return(gamma)
}

zwietering_gamma <- function(x, xmin, xopt, n) {
  gamma <- ((x-xmin)/(xopt-xmin))^n
  gamma[x < xmin] <- 0
  gamma[x > xopt] <- 0
  return(gamma)
}

full_Ratkowski <- function(x, xmin, xmax, c) {
  b <- 1 # Does not affect predictions (see supp. material)
  xopt <- (lambertW0(exp(-xmin*c + xmax*c + 1)) + c*xmin - 1)/c
  mu_opt <- b*(xopt - xmin)*(1 - exp(c*(xopt - xmax)))
  gamma <- b*(x - xmin)*(1 - exp(c*(x - xmax)))
  gamma <- gamma/mu_opt
  gamma <- gamma^2
  gamma[x < xmin] <- 0
  gamma[x > xmax] <- 0
  return(gamma)
}

reduced_Ratkowski <- function(x, xmin, b, clade){
  mu_opt <- b * (xopt_func(clade) - xmin)
  gamma <- b * (x - xmin)
  gamma <- gamma/mu_opt
  gamma <- gamma^2
  gamma[x < xmin] <- 0
  return(gamma)
}

# Define calculate_gammas
calculate_gammas <- function (this_t, env_func, sec_models) 
{
  out <- lapply(names(sec_models), function(this_condition) {
    this_x <- env_func[[this_condition]](this_t)
    this_sec <- sec_models[[this_condition]]
    this_gamma <- switch(this_sec$model, 
                         fullRatkowsky = full_Ratkowski(this_x, this_sec$xmin, this_sec$xmax, this_sec$c), 
                         CPM = CPM_model(this_x, this_sec$xmin, this_sec$xopt, this_sec$xmax, this_sec$n), 
                         Zwietering = zwietering_gamma(this_x, this_sec$xmin, this_sec$xopt, this_sec$n), 
                         reducedRatkowsky = reduced_Ratkowski(this_x, this_sec$xmin, this_sec$b, this_sec$clade), 
                         stop(paste("Model",this_sec$model, "not known.")))
    this_gamma
  })
  out <- unlist(out)
  names(out) <- names(sec_models)
  out
}

# Secondary model data
secondary_model_data <- function (model_name = NULL) {
  model_data <- list(CPM = list(identifier = "CPM", 
                                name = "Cardinal Parameter Model", 
                                pars = c("xmin", "xopt", "xmax", "n"), 
                                model = CPM_model, 
                                ref = paste("Rosso, L., Lobry, J. R., Bajard, S., and Flandrois, J. P. (1995).", 
                                            "Convenient Model To Describe the Combined Effects of Temperature and pH on", 
                                            "Microbial Growth. Applied and Environmental Microbiology, 61(2), 610-616.")), 
                     
                     Zwietering = list(identifier = "Zwietering", 
                                       name = "Zwietering gamma function", 
                                       pars = c("xmin", "xopt", "n"), 
                                       model = zwietering_gamma, 
                                       ref = paste("Zwietering, Marcel H., Wijtzes, T., De Wit, J. C., and Riet,", 
                                                   "K. V. (1992). A Decision Support System for Prediction of the Microbial", 
                                                   "Spoilage in Foods. Journal of Food Protection, 55(12), 973-979.", 
                                                   "https://doi.org/10.4315/0362-028X-55.12.973")), 
                     
                     fullRatkowsky = list(identifier = "fullRatkowsky", 
                                          name = "(Adapted) Full Ratkowsky model", 
                                          pars = c("xmin", "xmax", "c"), 
                                          model = full_Ratkowski, 
                                          ref = paste("Ratkowsky, D. A., Lowry, R. K., McMeekin, T. A.,", 
                                                      "Stokes, A. N., and Chandler, R. E. (1983). Model for", 
                                                      "bacterial culture growth rate throughout the entire", 
                                                      "biokinetic temperature range. Journal of Bacteriology,", 
                                                      "154(3), 1222-1226.")),
                     
                     reducedRatkowsky = list(identifier = "reducedRatkowsky", 
                                             name = "Reduced Ratkowsky model", 
                                             pars = c("xmin", "b", "clade"),
                                             model = reduced_Ratkowski, 
                                             ref = paste("Ratkowsky, D. A., Lowry, R. K., McMeekin, T. A.,", 
                                                         "Stokes, A. N., and Chandler, R. E. (1983). Model for", 
                                                         "bacterial culture growth rate throughout the entire", 
                                                         "biokinetic temperature range. Journal of Bacteriology,", 
                                                         "154(3), 1222-1226.")))
  
  
  if (is.null(model_name)) {
    return(names(model_data))
  }
  my_model <- model_data[[model_name]]
  if (is.null(my_model)) {
    stop(paste("Unknown model name:", model_name))
  }
  else {
    my_model
  }
}

# Define function to select Topt by groups (from literature, create citation)
xopt_func <- function(clade){
  if(clade == "I")
  {return(28)}
  if(clade == "II")
  {return(36.31)}
  else if(clade == "III")
  {return(39.27)}
  else if(clade == "IV")
  {return(38.735)}
  else if(clade == "V")
  {return(37.375)}
  else (clade == "VII")
  {return(42.35)}
}

# Define check_secondary_pars
check_secondary_pars <- function (starting_point, known_pars, sec_model_names, primary_pars = "mu_opt") 
{
  if (any(names(starting_point) %in% names(known_pars))) {
    stop("Parameters cannot be defined as both fixed and to be fitted.")
  }
  par_names <- c(names(starting_point), names(known_pars))
  missing_primary <- primary_pars[!primary_pars %in% par_names]
  if (length(missing_primary) > 0) {
    stop(paste("Parameter not defined:", missing_primary, 
               "\n"))
  }
  my_regex <- paste(primary_pars, collapse = "|")
  par_names <- par_names[!grepl(my_regex, par_names)]
  for (each_factor in names(sec_model_names)) {
    model_data <- secondary_model_data(sec_model_names[[each_factor]])
    req_pars <- paste0(each_factor, "_", model_data$pars)
    missing_pars <- req_pars[!req_pars %in% par_names]
    if (length(missing_pars) > 0) {
      stop(paste("Parameter not defined:", missing_pars, 
                 "\n"))
    }
    this_pars <- par_names[grepl(paste0(each_factor, "_"), 
                                 par_names)]
    unknown_pars <- this_pars[!this_pars %in% req_pars]
    if (length(unknown_pars) > 0) {
      stop(paste("Unknown parameter: ", unknown_pars, 
                 "\n"))
    }
  }
  my_regex <- paste(paste0(names(sec_model_names), "_"), collapse = "|")
  wtpars <- par_names[!grepl(my_regex, par_names)]
  if (length(wtpars) > 0) {
    stop(paste("Unknown parameter: ", wtpars, "\n"))
  }
}

# Define approx_env
approx_env <- function (env_conditions) 
{
  out <- lapply(names(env_conditions[-1]), function(this_col) {
    x <- env_conditions$time
    y <- env_conditions[[this_col]]
    approxfun(x, y, rule = 2)
  })
  names(out) <- names(env_conditions[-1])
  out
}

# Define predict_dynamic_growth
predict_dynamic_growth <- function (times, env_conditions, primary_pars, secondary_models, 
          ..., check = TRUE, logbase_logN = 10, logbase_mu = logbase_logN, 
          formula = . ~ time) 
{
  x_col <- formula.tools:::rhs(formula)
  env_conditions <- rename(env_conditions, time = x_col)
  if (isTRUE(check)) {
    sec_model_names <- secondary_models %>% map(~.$model) %>% 
      unlist()
    check_pars_names <- lapply(names(secondary_models), 
                               function(each_factor) {
                                 this_model <- secondary_models[[each_factor]]
                                 this_model$model <- NULL
                                 paste0(each_factor, "_", names(this_model))
                               }) %>% unlist()
    check_pars <- rep(1, length(check_pars_names))
    names(check_pars) <- check_pars_names
    check_secondary_pars(check_pars, unlist(primary_pars), 
                         sec_model_names, primary_pars = c("mu_opt", "N0", 
                                                           "Nmax", "Q0"))
  }
  primary_pars_calc <- primary_pars
  primary_pars_calc$mu_opt <- primary_pars_calc$mu_opt/log(10, 
                                                           base = logbase_mu)
  my_env <- approx_env(env_conditions)
  yini <- c(Q = primary_pars_calc$Q0, N = primary_pars_calc$N0)
  my_sim <- ode(yini, times, dBaranyi, primary_pars_calc, 
                env_func = my_env, sec_models = secondary_models, ...) %>% 
    as.data.frame() %>% as_tibble() %>% mutate(logN = log(.data$N, 
                                                          base = logbase_logN))
  gammas <- lapply(times, function(x) {
    calculate_gammas(x, my_env, secondary_models)
  })
  gammas <- do.call(rbind.data.frame, gammas)
  names(gammas) <- names(secondary_models)
  gammas <- bind_cols(time = times, gammas)
  out <- list(simulation = my_sim, gammas = as_tibble(gammas), 
              env_conditions = my_env, primary_pars = primary_pars, 
              sec_models = secondary_models, logbase_mu = logbase_mu, 
              logbase_logN = logbase_logN)
  class(out) <- c("DynamicGrowth", class(out))
  out
}
