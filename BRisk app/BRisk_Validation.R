# Load packages
library(tidyverse)
library(tibble)
library(EnvStats)         # to load rtri function 
library(truncnorm)        # to load rtruncnorm function
library(jmuOutlier)       # to load rlaplace function
library(formula.tools)    # to load 'rhs'
library(purrr)            # to load 'map'
library(deSolve)          # to load 'ode'
library(rlang)
library(ggplot2)

# Load utility functions
source("UtilityFunctions_dynamic_growth.R")

# Generate database 
# BTyper data
BTyper3_input = read.csv("Btyper3_Results.csv")
colnames(BTyper3_input)[1] <- "Isolate.Name"
gp_input = read.csv("simulation_input.csv")
database = cbind(BTyper3_input,gp_input[,3:7])
database <- database %>% 
  separate(Closest_Type_Strain.ANI., into = c("species","ANI"), sep = "\\(") %>%
  mutate(ANI = gsub("\\)", "", ANI))

# Input BTyper3 result for a B cereus isolate 
df <- read.csv("Genomic data for detected isolates/PS00087_contigs_final_results.csv")
df <- df %>% 
  separate(Closest_Type_Strain.ANI., into = c("species","ANI"), sep = "\\(") %>%
  separate(Adjusted_panC_Group.predicted_species., into = c("panC_Group","predicted_species"), sep = "\\(") %>%
  mutate(ANI = gsub("\\)", "", ANI),
         panC_Group = gsub("\\)", "", panC_Group),
         predicted_species = gsub("\\)", "", predicted_species))

# Input BTyper3 result for ANI
ANI_file <- read.csv("Genomic data for detected isolates/PS00087_contigs_cytotoxicity_fastani.csv", header = FALSE)
colnames(ANI_file) <- c("query", "reference", "ANI", "matched_genes", "total_genes")
ANI_file <- ANI_file[, c("reference", "ANI")]
ANI_file$reference <- sub("^(PS\\d+).*", "\\1", ANI_file$reference)
ANI_file$num_id <- as.numeric(sub("PS", "", ANI_file$reference))
ANI_file <- ANI_file[order(ANI_file$num_id), ]
ANI_file <- ANI_file[, c("reference", "ANI")]

# Filter the database input for rows with the same species as the BTyper3 input
df$species <- trimws(df$species)
database$ANI_new = ANI_file$ANI #Adding new ANI
matching_species_df <- subset(database, species == df$species)
# Assigning ANI weight 
matching_species_df$ANI_wght <- matching_species_df$ANI_new / sum(matching_species_df$ANI_new)

# Simulate HTST milk products along the supply chain
## Set seed
set.seed(1)

## Assign isolate names to 10,000 units of HTST milk products
## Isolates from the same species are represented by weight determined by ANI
n_sim = 10000
matching_species_df$n_units <- round(n_sim * matching_species_df$ANI_wght)
sampled_isolates <- rep(matching_species_df$Isolate.Name, matching_species_df$n_units)
if(length(sampled_isolates) != n_sim){
  sampled_isolates <- sample(sampled_isolates, n_sim)
}

ModelData <- data.frame(
  unit_id = seq_len(n_sim),
  isolate = sampled_isolates
)

# Validation temp profiles 
# Stage 1: facility storage 
## (a)  Sample the temperature distribution
ModelData$T_F = 4
## (b) Sample the storage time (in days) distribution
ModelData$t_F = 1.5

# Stage 2: transport from facility to retail store
## (a)  Sample the temperature distribution
ModelData$T_T = 10
## (b) Sample the transportation time (in days) distribution
ModelData$t_T = 5

# Stage 3: storage/display at retail store
## (a)  Sample the temperature distribution
ModelData$T_S = 4
## (b) Sample the storage time (in days) distribution
ModelData$t_S = 2

## Stage 4: transportation from retail store to home
## (a)  Sample the temperature distribution
ModelData$T_T2 = 10
## (b) Sample the transportation time (in days) distribution 
ModelData$t_T2 = 1/24

## Stage 5: home storage 
## (a)  Sample the temperature distribution
ModelData$T_H <- 10
## (b) Define t_H as 35 days for all units
ModelData$t_H <- rep(35, each = n_sim)

## Model temperature profiles of 10,000 units HTST milk 
env_cond_time <- matrix(c(rep(0,n_sim),
                          ModelData$t_F, 
                          ModelData$t_F+0.001,
                          ModelData$t_F + ModelData$t_T,
                          ModelData$t_F + ModelData$t_T+0.001,
                          ModelData$t_F + ModelData$t_T + ModelData$t_S,
                          ModelData$t_F + ModelData$t_T + ModelData$t_S+0.001,
                          ModelData$t_F + ModelData$t_T + ModelData$t_S + ModelData$t_T2,
                          ModelData$t_F + ModelData$t_T + ModelData$t_S + ModelData$t_T2+0.001,
                          ModelData$t_F + ModelData$t_T + ModelData$t_S + ModelData$t_T2 + ModelData$t_H), ncol = 10)

env_cond_temp <- matrix(c(ModelData$T_F, 
                          ModelData$T_F,
                          ModelData$T_T,
                          ModelData$T_T,
                          ModelData$T_S,
                          ModelData$T_S,
                          ModelData$T_T2,
                          ModelData$T_T2,
                          ModelData$T_H,
                          ModelData$T_H), ncol = 10)

## Generate simulation input 
## Assign growth parameters to 10,000 units of HTST milk 
ModelData$index = match(ModelData$isolate, matching_species_df$Isolate.Name)
ModelData$Q0 = matching_species_df$Q0[ModelData$index]
ModelData$Nmax = matching_species_df$Nmax[ModelData$index]
ModelData$b = matching_species_df$b[ModelData$index]
ModelData$Tmin = matching_species_df$Tmin[ModelData$index]
ModelData$Clade = matching_species_df$Clade[ModelData$index]

# Validation N0 
N0 = 100
ModelData$N0 = N0 

ModelData$Topt = sapply(ModelData$Clade, xopt_func)
ModelData$mu_opt = (ModelData$b*(ModelData$Topt-ModelData$Tmin))^2

# Run simulation
for (i in 1:nrow(ModelData)){
  my_primary <- list(mu_opt = ModelData$mu_opt[i], Nmax = ModelData$Nmax[i], N0 = ModelData$N0[i], Q0 = ModelData$Q0[i])
  sec_temperature <- list(model = "reducedRatkowsky", xmin = ModelData$Tmin[i], b = ModelData$b[i], clade = ModelData$Clade[i])
  my_secondary <- list(temperature = sec_temperature)
  growth <- predict_dynamic_growth(times = env_cond_time[i,],
                                   env_conditions = tibble(time = env_cond_time[i,],
                                                           temperature = env_cond_temp[i,]),
                                   my_primary,
                                   my_secondary)
  sim <- growth$simulation
  ModelData$conc[i] = tail(sim$logN, 1)
}

quantile(ModelData$conc, 0.25)
median(ModelData$conc)
quantile(ModelData$conc, 0.75)

