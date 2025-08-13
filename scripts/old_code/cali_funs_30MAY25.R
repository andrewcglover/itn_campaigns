# cali_funs.R
# local cali functions required for cluster due to malariasimulation version
# conflict

cali_pf_2_10 <- function (eir, ft = 0.1) 
{
  eq <- malariaEquilibrium::human_equilibrium(eir, ft = ft, 
                                              p = malariaEquilibrium::load_parameter_set("Jamie_parameters.rds"), 
                                              age = 0:100)
  sum(eq$states[3:11, "pos_M"])/sum(eq$states[3:11, "prop"])
}

cali_pfpr_6_59_mo <- function (eir, ft = 0.1) 
{
  eq <- malariaEquilibrium::human_equilibrium(eir, ft = ft, 
                                              p = malariaEquilibrium::load_parameter_set("Jamie_parameters.rds"), 
                                              age = 0:100)
  sum(eq$states[2:6, "pos_M"], eq$states[1, "pos_M"]/2)/sum(eq$states[2:6, 
                                                                      "prop"], eq$states[1, "prop"]/2)
}

cali_eq_objective <- function (eir, target_pfpr, ft) 
{
  cali_pf_2_10(eir, ft) - target_pfpr
}

cali_eq_objective_6_59_mo <- function (eir, target_pfpr, ft) 
{
  cali_pfpr_6_59_mo(eir, ft) - target_pfpr
}

cali_get_eq_eir <- function (target_pfpr, ft = 0.1, pfpr_6_59_on = FALSE) 
{
  if (target_pfpr <= 0 | target_pfpr > 0.9) {
    stop("eq_pfpr must be between 0 and 0.9")
  }
  if (ft < 0 | ft >= 1) {
    stop("ft must be between 0 and 1")
  }
  if (pfpr_6_59_on) {
    opt <- stats::uniroot(f = cali_eq_objective_6_59_mo, lower = .Machine$double.eps, 
                          upper = 200, target_pfpr = target_pfpr, ft = ft, 
                          extendInt = "upX")
  }
  else {
    opt <- stats::uniroot(f = cali_eq_objective, lower = .Machine$double.eps, 
                          upper = 200, target_pfpr = target_pfpr, ft = ft, 
                          extendInt = "upX")
  }
  opt$root
}

cali_check_elimination <- function (output, target) 
{
  any(output == 0 & target != 0)
}

cali_check_fit <- function (proposed_eir, parameters, target, summary_function) 
{
  parameters <- malariasimulation::set_equilibrium(parameters, 
                                                   init_EIR = proposed_eir)
  simulation <- malariasimulation::run_simulation(timesteps = parameters$timesteps, 
                                                  parameters = parameters)
  output <- summary_function(simulation)
  elimination <- cali_check_elimination(output, target)
  objective <- ifelse(elimination, NA, sum(output - target))
  return(objective)
}

cali_to_real <- function (x, limits) 
{
  log(x - limits[1]) - log(limits[2] - x)
}

cali_from_real <- function (x, limits) 
{
  (limits[2] * exp(x) + limits[1])/(1 + exp(x))
}

cali_proposal <- function (current_eir, limits, direction, step = 0.5) 
{
  if (current_eir < limits[1] | current_eir > limits[2]) {
    stop("Current eir outside of limits")
  }
  r <- cali_to_real(current_eir, limits)
  r_prop <- ifelse(direction == "decrease", r - step, r + step)
  proposed_eir <- cali_from_real(r_prop, limits)
  return(proposed_eir)
}

cali_good_move <- function (objective, direction) 
{
  if (direction == "decrease") {
    cali_good_move <- objective[2] <= 0
  }
  else {
    cali_good_move <- objective[2] >= 0
  }
  return(cali_good_move)
}

cali_linear_interpolate <- function (eir, objective) 
{
  eir <- sort(eir)
  objective <- sort(objective)
  if (objective[1] == objective[2] | diff(sign(objective)) == 
      0) {
    return(mean(eir, na.rm = TRUE))
  }
  b <- (eir[2] - eir[1])/(objective[2] - objective[1])
  eir[1] - b * objective[1]
}

cali_calibrate <- function (parameters, target, summary_function, eq_prevalence, 
                            eq_ft = 0, human_population = c(1000, 10000, 1e+05), eir_limits = c(0, 
                                                                                                1500), max_attempts = 10, use_pfpr_6_59_mo = FALSE) 
{
  eir <- rep(0, 2)
  objective <- rep(NA, 2)
  parameters$human_population <- human_population[1]
  eir[1] <- cali_get_eq_eir(target_pfpr = eq_prevalence, ft = eq_ft, 
                       pfpr_6_59_on = use_pfpr_6_59_mo)
  if (eir[1] < eir_limits[1]) {
    eir[1] <- eir_limits[1] + 1
  }
  if (eir[1] > eir_limits[2]) {
    eir[1] <- eir_limits[2] - 1
  }
  min_eir <- 0
  attempts <- 0
  while (is.na(objective[1])) {
    objective[1] <- cali_check_fit(proposed_eir = eir[1], parameters = parameters, 
                              target = target, summary_function = summary_function)
    attempts <- attempts + 1
    if (is.na(objective[1])) {
      if (parameters$human_population < max(human_population)) {
        parameters$human_population <- human_population[which(human_population == 
                                                                parameters$human_population) + 1]
      }
      else {
        parameters$human_population <- human_population[1]
        min_eir <- eir[1]
        eir[1] <- cali_proposal(current_eir = eir[1], limits = eir_limits, 
                           direction = "increase", step = 2)
      }
    }
    if (attempts > max_attempts) {
      stop("Failure due to max attempts reached before successful run")
    }
  }
  direction <- ifelse(objective[1] > 0, "decrease", "increase")
  if (direction == "decrease") {
    eir_limits[1] <- min_eir
  }
  estimated_eir <- NA
  eir[2] <- eir[1]
  repeat {
    eir[2] <- cali_proposal(current_eir = eir[2], limits = eir_limits, 
                       direction = direction)
    objective[2] <- cali_check_fit(proposed_eir = eir[2], parameters = parameters, 
                              target = target, summary_function = summary_function)
    if (!is.na(objective[2])) {
      if (cali_good_move(objective, direction)) {
        estimated_eir <- cali_linear_interpolate(eir, objective)
        break
      }
      else {
        eir[1] <- eir[2]
        objective[1] <- objective[2]
      }
    }
    else {
      if (direction == "decrease") {
        eir_limits[1] <- eir[2]
        direction <- "increase"
      }
    }
    attempts <- attempts + 1
    if (attempts > max_attempts) {
      estimated_eir <- mean(eir, na.rm = TRUE)
      break
    }
  }
  return(estimated_eir)
}