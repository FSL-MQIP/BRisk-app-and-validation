setwd("C:/Users/sujun/Documents/GitHub/Bacillus-cereus-exposure-assessment-model")

# load packages 
library(tibble)
library(EnvStats)         # to load rtri function 
library(truncnorm)        # to load rtruncnorm function
library(jmuOutlier)       # to load rlaplace function
library(dplyr)
library(tidyr)
library(formula.tools)    # to load 'rhs'
library(purrr)            # to load 'map'
library(deSolve)          # to load 'ode'
library(rlang)

# load utility functions
source("Exposure assessment model/UtilityFunctions_dynamic_growth.R")

## Set up dataframe for modeling 100 units of HTST milk 
n_sim = 1
data = data.frame(unit_id = rep(seq(1,n_sim)))

# Stage 1: facility storage 
## (a)  temperature
data$T_F <- 4
## (b) storage time (in days)
data$t_F <- 1.5
# Stage 2: transport from facility to retail store
## (a)  temperature
data$T_T <- 10
## (b) transportation time (in days)
data$t_T <- 5
# Stage 3: storage/display at retail store
## (a)  temperature
data$T_S <- 4
## (b) Sample the storage time (in days) distribution
data$t_S <- 2
## Stage 4: transportation from retail store to home
## (a)  temperature
data$T_T2 <- 10
## (b) transportation time (in days) 
data$t_T2 <- 0.04
## Stage 5: home storage 
## (a)  Sample the temperature distribution
data$T_H <- 10
## (b) transportation time (in days) 
data$t_H <- 14
# data$t_H <- 21
# data$t_H <- 35

## Model temperature profiles of n_sim units HTST milk 
env_cond_time <- matrix(c(rep(0,n_sim),
                          data$t_F, 
                          data$t_F+0.001,
                          data$t_F + data$t_T,
                          data$t_F + data$t_T+0.001,
                          data$t_F + data$t_T + data$t_S,
                          data$t_F + data$t_T + data$t_S+0.001,
                          data$t_F + data$t_T + data$t_S + data$t_T2,
                          data$t_F + data$t_T + data$t_S + data$t_T2+0.001,
                          data$t_F + data$t_T + data$t_S + data$t_T2 + data$t_H), ncol = 10)

env_cond_temp <- matrix(c(data$T_F, 
                          data$T_F,
                          data$T_T,
                          data$T_T,
                          data$T_S,
                          data$T_S,
                          data$T_T2,
                          data$T_T2,
                          data$T_H,
                          data$T_H), ncol = 10)

# Import data set
data_Q0 = read.csv("Exposure assessment model/InputFiles/Q0_h0_summary.csv")
data_Nmax = read.csv("Exposure assessment model/InputFiles/Nmax_new.csv")
data_sec_model = read.csv("Exposure assessment model/InputFiles/sec_model_new.csv")
Clade = c("I","II","VII","IV","IV","IV","IV","II","III","IV","II","VII","II","V","V","IV")
# Generate N0
N0 = 1000

# Generate simulation input
# Input Q0
simulation_input <- data.frame(isolate = data_Q0$isolate, Q0 = data_Q0$Q0, 
                               Nmax = data_Nmax$average_Nmax, b = data_sec_model$Slope, 
                               Tmin = data_sec_model$Tmin,Clade)

simulation_input$Topt = sapply(simulation_input$Clade, xopt_func)
simulation_input$mu_opt = (simulation_input$b*(simulation_input$Topt-simulation_input$Tmin))^2 

# Run simulation
final_conc <- simulation_input %>%
  rowwise() %>%
  mutate(final_conc_isolate = list(lapply(1:n_sim, function(j) {
    my_primary <- list(mu_opt = mu_opt, Nmax = Nmax, N0 = N0[j], Q0 = Q0)
    sec_temperature <- list(model = "reducedRatkowsky", xmin = Tmin, b = b, clade = Clade)
    my_secondary <- list(temperature = sec_temperature)
    growth <- predict_dynamic_growth(times = env_cond_time[j,],
                                     env_conditions = tibble(time = env_cond_time[j,],
                                                             temperature = env_cond_temp[j,]),
                                     my_primary,
                                     my_secondary)
    sim <- growth$simulation
    return(tail(sim$logN, 1))
  }))) %>%
  pull(final_conc_isolate)

# Convert the list of 16 elements into a matrix
matrix <- t(sapply(final_conc, function(x) sapply(x, tail, n=1)))
final_count_d14 <- as.list(matrix)
# final_count_d21 <- as.list(matrix)
# final_count_d35 <- as.list(matrix)
simulation_input$final_count_d14 = final_count_d14
# simulation_input$final_count_d21 = final_count_d21
# simulation_input$final_count_d35 = final_count_d35
simulation_input <- simulation_input[order(simulation_input$Clade), ]
Validation_1 <- simulation_input[,c(1,6,9:11)]
Validation_1 <- Validation_1 %>%
  unnest(final_count_d14) %>%
  unnest(final_count_d21) %>%
  unnest(final_count_d35)
write.csv(Validation_1,"Validation/Validation_1.csv")

