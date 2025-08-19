library(SimInf)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(lubridate)
library(INLA)
library(purrr)
library(sf)
library(xlsx)
library(openxlsx)
library(stringr)
library(vegan)
library(readxl)
library(reshape2)
library(patchwork)
library(ggpubr)
library(DEoptim)
library(GenSA)
library(nloptr)

#hpai model
n <- 100

compartments <- c("S", "I", "C", "V")

transitions <- c("@ -> Lambda*S -> S", "S -> mu*S -> @", "S -> beta*S*I/(S+I+C+V) -> I", 
                 "S -> p*S -> V", "V -> mu*V -> @", "V -> (1 - phi)*beta*V*I/(S+I+C+V) -> I",
                 "I -> mu*I -> @", "I -> Cu*I -> C")

u0 <- data.frame(S = rep(95, n), I = rep(5, n), C = rep(0, n), V = rep(0, n))
model <- mparse(transitions = transitions, compartments = compartments,
                gdata = c(mu = 0.05, p = 0.2, Lambda = 0.05, beta = 0.8, phi = 0.5, Cu = 0.3), 
                u0 = u0, tspan = 1:50)


result <- run(model = model)
plot(result)

#model 2, with PE culling
transitions_CuPe <- c("@ -> Lambda*S -> S", "S -> mu*S -> @", "S -> beta*S*I/(S+I+C+V) -> I", 
                 "S -> p*S -> V", "V -> mu*V -> @", "V -> (1 - phi)*beta*V*I/(S+I+C+V) -> I",
                 "I -> mu*I -> @", "I -> Cu*I -> C", "S -> CuPe*I*S -> C")

u0 <- data.frame(S = rep(95, n), I = rep(5, n), C = rep(0, n), V = rep(0, n))
model_CuPe <- mparse(transitions = transitions_CuPe, compartments = compartments,
                gdata = c(mu = 0.05, p = 0.2, Lambda = 0.05, beta = 0.8, phi = 0.5, 
                          Cu = 0.5, CuPe=0.3), 
                u0 = u0, tspan = 1:50)


result <- run(model = model)
plot(result)

#model 3, dens dep transmission
compartments2 <- c("S", "I", "Icum", "C", "Cpe", "V", "Vcum")

transitions_densdep <- c("@ -> Lambda*S -> S", "S -> mu*S -> @", "S -> beta*S*I -> I + Icum", 
                      "S -> p*S -> V + Vcum", "V -> mu*V -> @", "V -> (1 - phi)*beta*V*I -> I + Icum",
                      "I -> mu*I -> @", "I -> Cu*I -> C", "S -> CuPe*I*S -> Cpe")

u0 <- data.frame(S = rep(99, n), I = rep(1, n), Icum=rep(0,n), C = rep(0, n), Cpe = rep(0, n), V = rep(0, n), Vcum = rep(0, n))
model_densdep <- mparse(transitions = transitions_densdep, compartments = compartments2,
                     gdata = c(mu = 0.05, p = 0.1, Lambda = 0.05, beta = 0.02, phi = 0.5, 
                               Cu = 0.2, CuPe=0.1), 
                     u0 = u0, tspan = 1:50)


result <- run(model = model_densdep)
plot(result)

traj <- trajectory(model = result, compartments = c("Icum", "Vcum", "C", "Cpe")) %>% 
  group_by(time) %>% 
  summarise(Icum_avg=mean(Icum), Vcum_avg=mean(Vcum), C_avg=mean(C), Cpe_avg=mean(Cpe))


#function to minimize
cost_function <- function(x) {
  cull_prop <- x[1]
  vacc_coverage <- x[2]
  
  model <- mparse(transitions = transitions, compartments = compartments,
                  gdata = c(mu = 0.05, p = vacc_coverage, Lambda = 5, beta = 0.8, phi = 0.5, Cu = cull_prop), 
                  u0 = u0, tspan = 1:180)
  
  result <- run(model = model)
  
  # Compute total infected over time
  id_I<-seq(2,dim(result@U)[1],4)
  tot_infected <- rowSums(result@U[id_I,])
  avg_tot_infected <- mean(tot_infected)
  
  return(avg_tot_infected)  #what we minimize
}

### TEST SEVERAL OPTI PACKAGES

# optim: gets trapped in local minima

opt <- optim(
  par = c(0.25, 0.6),         # initial values: culling = 20%, vaccination = 20%
  fn = cost_function,        # objective function
  method = "L-BFGS-B",
  lower = c(0, 0),           # min values
  upper = c(0.3, 0.9)            # max values
)

opt$par       # Optimal [culling, vaccination]
opt$value 

# DEoptim: long

result_deoptim <- DEoptim(
  fn = cost_function,
  lower = c(0, 0),
  upper = c(0.3, 0.9),
  control = DEoptim.control(
    NP = 40,         
    itermax = 200,   
    parallelType = "none" 
  )
)

best_params <- result_deoptim$optim$bestmem #best combination

#plot the trace
param_trace <- as.data.frame(result_deoptim$member$bestmemit)
colnames(param_trace) <- c("culling", "vaccination")

param_trace$iteration <- 1:nrow(param_trace)
param_trace$objective <- result_deoptim$member$bestvalit

ggplot(param_trace, aes(x = culling, y = vaccination)) +
  geom_path(aes(color = objective), size = 1) +
  geom_point(aes(color = objective), size = 2) +
  scale_color_viridis_c(option = "plasma", direction = -1) +
  geom_hline(yintercept=0.9, linetype="dashed", color = "red") +
  geom_vline(xintercept=0.3, linetype="dashed", color = "red") +
  theme_minimal() +
  labs(
    title = "Optimization Trace of DEoptim",
    x = "Culling Proportion",
    y = "Vaccination Coverage",
    color = "Objective"
  )


# GenSA

result_gensa <- GenSA(
  par = c(0.2, 0.2),
  fn = cost_function,
  lower = c(0, 0),
  upper = c(0.3, 0.9)
)

best_params <- result_gensa$par

# nloptr: gradient ?

result_nloptr <- nloptr(
  x0 = c(0.2, 0.2),
  eval_f = cost_function,
  lb = c(0, 0),
  ub = c(0.3, 0.9),
  opts = list(
    "algorithm" = "NLOPT_LN_COBYLA",  
    "xtol_rel" = 1.0e-8)
)

best_params <- result_nloptr$solution

result_nloptr <- nloptr(
  x0 = c(0.2, 0.2),
  eval_f = cost_function,  
  lb = c(0, 0),
  ub = c(0.3, 0.9),
  opts = list(
    algorithm = "NLOPT_GN_MLSL_LDS",
    maxeval = 500,
    xtol_rel = 1.0e-8,
    local_opts = list(
      algorithm = "NLOPT_LN_COBYLA",
      xtol_rel = 1.0e-6  
    )
  )
)

best_params <- result_nloptr$solution
result_nloptr$status


### New cost function

cost_function2 <- function(x, vacc_cost, cull_cost) {
  cull_prop <- x[1]
  vacc_coverage <- x[2]
  
  u0 <- data.frame(S = rep(95, n), I = rep(5, n), C = rep(0, n), V = rep(0, n))
  
  model <- mparse(transitions = transitions, compartments = compartments,
                  gdata = c(mu = 0.05, p = vacc_coverage, Lambda = 5, beta = 0.8, phi = 0.5, Cu = cull_prop), 
                  u0 = u0, tspan = 1:180)
  
  result <- run(model = model)
  
  
  # Compute total infected over time
  id_I<-seq(2,dim(result@U)[1],4)
  tot_infected <- rowSums(result@U[id_I,])
  avg_tot_infected <- mean(tot_infected)
  
  cost_inf <- 12
  
  obj<-vacc_cost*vacc_coverage*u0$S[1] + cull_cost*cull_prop*u0$S[1] + cost_inf*avg_tot_infected
  
  return(obj)  
}

#tested costs
vacc_costs <- seq(1, 10, by = 1)
cull_costs <- seq(1, 10, by = 1)
param_grid <- expand.grid(vacc_cost = vacc_costs, cull_cost = cull_costs)

#wrapper function
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
  pmap_dfr(optimize_wrapper) %>%
  mutate(ratio=cull_prop/vacc_coverage)

ggplot(results, aes(x = vacc_cost, y = cull_cost, fill = vacc_coverage)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Vaccination to Culling Ratio", x = "Vaccination Cost", y = "Culling Cost") +
  theme_minimal()

ggplot(results, aes(x = vacc_cost, y = cull_cost, fill = cull_prop)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Vaccination to Culling Ratio", x = "Vaccination Cost", y = "Culling Cost") +
  theme_minimal()

ggplot(results, aes(x = vacc_cost, y = cull_cost, fill = ratio)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Vaccination to Culling Ratio", x = "Vaccination Cost", y = "Culling cost") +
  theme_minimal()


### Cost function w PE culling

cost_function3 <- function(x, vacc_cost, cull_cost) {
  cull_prop <- x[1]
  pe_cull_prop <- x[2]
  vacc_coverage <- x[3]
  
  u0 <- data.frame(S = rep(99, n), I = rep(1, n), Icum=rep(0,n), C = rep(0, n), Cpe = rep(0, n), V = rep(0, n), Vcum = rep(0, n))
  model_densdep <- mparse(transitions = transitions_densdep, compartments = compartments2,
                          gdata = c(mu = 0.05, p = vacc_coverage, Lambda = 0.05, beta = 0.02, phi = 0.5, 
                                    Cu = cull_prop, CuPe=pe_cull_prop), 
                          u0 = u0, tspan = 1:50)
  
  
  result <- run(model = model_densdep)
  
  # Compute total infected over time
  traj <- trajectory(model = result, compartments = c("Icum", "Vcum", "C", "Cpe")) %>% 
    group_by(time) %>% 
    summarise(Icum_avg=mean(Icum), Vcum_avg=mean(Vcum), C_avg=mean(C), Cpe_avg=mean(Cpe))

  cost_inf <- 12
  
  obj<-vacc_cost*tail(traj$Vcum_avg, n=1) + cull_cost*pe_cull_prop*tail(traj$Cpe_avg, n=1) + cost_inf*tail(traj$Icum_avg, n=1)
  
  return(obj)  
}

#tested costs
vacc_costs <- seq(1, 10, by = 1)
cull_costs <- seq(1, 10, by = 1)
param_grid <- expand.grid(vacc_cost = vacc_costs, cull_cost = cull_costs)

#wrapper function
optimize_wrapper <- function(vacc_cost, cull_cost) {
  opt <- nloptr(
    x0 = c(0.2, 0.2, 0.2),
    eval_f = function(x) cost_function3(x, vacc_cost = vacc_cost, cull_cost = cull_cost),
    lb = c(0, 0, 0),
    ub = c(0.9, 0.3, 0.9),
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
    vacc_coverage = opt$solution[3],
    cull_prop = opt$solution[1],
    pe_cull_prop = opt$solution[2]
  )
}

start.time <- Sys.time()
results <- param_grid %>%
  pmap_dfr(optimize_wrapper) %>%
  mutate(ratio=cull_prop/vacc_coverage)
end.time <- Sys.time()

time.taken <- end.time - start.time
time.taken

ggplot(results, aes(x = vacc_cost, y = cull_cost, fill = vacc_coverage)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Vaccination coverage", x = "Vaccination Cost", y = "Culling Cost") +
  theme_minimal()
ggsave("./Figures/diff_costs_resultvacc_v2wPEcull.png", bg="white")

ggplot(results, aes(x = vacc_cost, y = cull_cost, fill = pe_cull_prop)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Preventive culling rate", x = "Vaccination Cost", y = "Culling Cost") +
  theme_minimal()
ggsave("./Figures/diff_costs_resultpecull_v2wPEcull.png", bg="white")

ggplot(results, aes(x = vacc_cost, y = cull_cost, fill = cull_prop)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Reactive culling rate", x = "Vaccination Cost", y = "Culling Cost") +
  theme_minimal()
ggsave("./Figures/diff_costs_resultcull_v2wPEcull.png", bg="white")

# ggplot(results, aes(x = vacc_cost, y = cull_cost, fill = ratio)) +
#   geom_tile() +
#   scale_fill_viridis_c() +
#   labs(title = "Vaccination to Culling Ratio", x = "Vaccination Cost", y = "Culling cost") +
#   theme_minimal()

#plot increase of costs
v1<-seq(0,1,by=0.1)
costs<-0.5*v1^2*200
plot(v1, costs)
