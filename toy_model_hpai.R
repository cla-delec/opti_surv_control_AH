library(SimInf)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(purrr)
library(ggpubr)
library(nloptr)

##########################################
# MODEL DEFINITION

# define simple hpai toy model w SimInf

model_simulation <- function(cull_prop, vacc_prop) {
  
  compartments <- c("S", "I", "Icum", "C", "V")

  transitions <- c("@ -> Lambda*S -> S", "S -> mu*S -> @", 
                   "S -> beta*S*I/(S+I+C+V) -> I + Icum", 
                   "S -> p*S -> V", "V -> mu*V -> @", 
                   "V -> (1 - phi)*beta*V*I/(S+I+C+V) -> I + Icum",
                   "I -> mu*I -> @", "I -> Cu*I -> C")
  
  n <- 100
  u0 <- data.frame(S = rep(95, n), I = rep(5, n), Icum = rep(5, n), 
                   C = rep(0, n), V = rep(0, n))
  
  model <- mparse(transitions = transitions, compartments = compartments,
                  gdata = c(mu = 0.05, p = vacc_prop, 
                            Lambda = 0.05, beta = 0.8, 
                            phi = 0.5, Cu = cull_prop), 
                  u0 = u0, tspan = 1:50)
  
  
  result <- run(model = model)
  #plot(result)
  
  return(trajectory(result))
}



##########################################
# OBJECTIVE FUNCTION DEFINITION

#metric to minimize

cost_function <- function(x) {
  
  # decision variables to be optimized
  cull_prop <- x[1]
  vacc_coverage <- x[2]
  
  # run model
  result <- model_simulation(cull_prop, vacc_coverage)
  
  # Compute total infected, metric to minimize
  n_inf <- result %>% 
    group_by(node) %>% 
    summarise(n_inf = max(Icum))
  
  return(mean(n_inf$n_inf))  #what we minimize
  
}

##########################################
# OPTIMIZATION

# using package nloptr

result_nloptr <- nloptr(
  x0 = c(0.2, 0.2),
  eval_f = cost_function,
  
  # constraints on decision variables
  lb = c(0, 0), # lower bounds
  ub = c(0.3, 0.9), # upper bounds
  
  opts = list(
    algorithm = "NLOPT_GN_MLSL_LDS",
    maxeval = 500,
    xtol_rel = 1.0e-8,
    local_opts = list(
      algorithm = "NLOPT_LN_COBYLA",
      xtol_rel = 1.0e-6))
)

best_params <- result_nloptr$solution
result_nloptr$status

##########################################
# OBJECTIVE FUNCTION 2 
# based on costs

cost_function2 <- function(x, vacc_cost, cull_cost) {
  
  # decision variables to be optimized
  cull_prop <- x[1]
  vacc_coverage <- x[2]
  
  # run model
  result <- model_simulation(cull_prop, vacc_coverage)
  
  # Compute total infected
  n_inf <- result %>% 
    group_by(node) %>% 
    summarise(n_inf = max(Icum))
  
  n_individuals <- result %>%
    filter(time==1) %>% group_by(node) %>%
    summarise(n_ind=S+I+C+V)
    
  
  cost_inf <- 12

  #cost function
  tot_costs <- vacc_cost*vacc_coverage*mean(n_individuals$n_ind) + 
    cull_cost*cull_prop*mean(n_individuals$n_ind) + 
    cost_inf*mean(n_inf$n_inf)
  
  return(tot_costs)  
}

#tested costs
vacc_costs <- seq(1, 10, by = 2)
cull_costs <- seq(1, 10, by = 2)
param_grid <- expand.grid(vacc_cost = vacc_costs, cull_cost = cull_costs)

#wrapper function - to test optimal parameter combination based on costs
optimize_wrapper <- function(vacc_cost, cull_cost) {
  opt <- nloptr(
    x0 = c(0.2, 0.2),
    eval_f = function(x) cost_function2(x, vacc_cost = vacc_cost, cull_cost = cull_cost),
    lb = c(0, 0),
    ub = c(0.3, 0.9),
    opts = list(
      algorithm = "NLOPT_GN_MLSL_LDS",
      maxeval = 300,
      xtol_rel = 1.0e-8,
      local_opts = list(
        algorithm = "NLOPT_LN_COBYLA",
        xtol_rel = 1.0e-6
      )
    )
  )
  
  tibble(
    vacc_cost = vacc_cost,
    cull_cost = cull_cost,
    vacc_coverage = opt$solution[2],
    cull_prop = opt$solution[1]
  )
}


results <- param_grid %>%
  pmap_dfr(optimize_wrapper) 

g1 <- ggplot(results, aes(x = vacc_cost, y = cull_cost, fill = vacc_coverage)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Vaccination coverage", x = "Vaccination Cost", y = "Culling Cost") +
  theme_minimal()

g2 <- ggplot(results, aes(x = vacc_cost, y = cull_cost, fill = cull_prop)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Preventive culling", x = "Vaccination Cost", y = "Culling Cost") +
  theme_minimal()

ggarrange(g1, g2)

