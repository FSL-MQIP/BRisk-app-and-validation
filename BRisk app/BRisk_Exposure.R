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
df <- read.csv("Genomic data for detected isolates/PS00597_contigs_final_results.csv")
df <- df %>% 
  separate(Closest_Type_Strain.ANI., into = c("species","ANI"), sep = "\\(") %>%
  separate(Adjusted_panC_Group.predicted_species., into = c("panC_Group","predicted_species"), sep = "\\(") %>%
  mutate(ANI = gsub("\\)", "", ANI),
         panC_Group = gsub("\\)", "", panC_Group),
         predicted_species = gsub("\\)", "", predicted_species))

# Input BTyper3 result for ANI
ANI_file <- read.csv("Genomic data for detected isolates/PS00597_contigs_cytotoxicity_fastani.csv", header = FALSE)
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

# Stage 1: facility storage 
## (a)  Sample the temperature distribution
ModelData$T_F <- rep(runif(n_sim,min=3.5,max=4.5)) #uniform distribution
## (b) Sample the storage time (in days) distribution
ModelData$t_F <- rep(runif(n_sim,min=1,max=2)) #uniform distribution

# Stage 2: transport from facility to retail store
## (a)  Sample the temperature distribution
ModelData$T_T <- rep(rtri(n_sim,min=1.7,max=10.0,mode=4.4)) #triangular distribution
## (b) Sample the transportation time (in days) distribution
ModelData$t_T <- rep(rtri(n_sim,min=1,max=10,mode=5))

# Stage 3: storage/display at retail store
## (a)  Sample the temperature distribution
ModelData$T_S <- rep(rtruncnorm(n_sim,a=-1.4,b=5.4,mean=2.3,sd=1.8)) #truncated normal distribution
## (b) Sample the storage time (in days) distribution
ModelData$t_S <- rep(rtruncnorm(n_sim,a=0.042,b=10.0, mean=1.821,sd=3.3)) #truncated normal distribution

## Stage 4: transportation from retail store to home
## (a)  Sample the temperature distribution
ModelData$T_T2 <- rep(rtruncnorm(n_sim,a=0,b=10,mean=8.5,sd=1.0)) #truncated normal distribution
## (b) Sample the transportation time (in days) distribution 
ModelData$t_T2 <- rep(rtruncnorm(n_sim,a=0.01,b=0.24, mean=0.04,sd=0.02)) #truncated normal distribution

## Stage 5: home storage 
## (a)  Sample the temperature distribution
temps <- rep(NA, n_sim)
for (i in 1:n_sim){
number <- rlaplace(1,m=4.06,s=2.31)
while (number > 15 | number < -1) {
number <- rlaplace(1,m=4.06,s=2.31) #truncated laplace distribution 
}
temps[i] <- number
}
ModelData$T_H <- temps

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

## Generate N0 from a Poisson distribution 
set.seed(42)
N0 = rpois(n = n_sim, lambda = 100)
ModelData$N0 = N0 

ModelData$Topt = sapply(ModelData$Clade, xopt_func)
ModelData$mu_opt = (ModelData$b*(ModelData$Topt-ModelData$Tmin))^2

# Run simulation
my_times <- seq(0,35)
num_iterations <- nrow(ModelData)
all_simulations <- list()
for (i in 1:num_iterations) {
  my_primary <- list(mu_opt = ModelData$mu_opt[i], Nmax = ModelData$Nmax[i], N0 = ModelData$N0[i], Q0 = ModelData$Q0[i])
  sec_temperature <- list(model = "reducedRatkowsky", xmin = ModelData$Tmin[i], b = ModelData$b[i], clade = ModelData$Clade[i])
  my_secondary <- list(temperature = sec_temperature)
  growth <- predict_dynamic_growth(times = my_times,
                                   env_conditions = tibble(time = env_cond_time[i,],
                                                           temperature = env_cond_temp[i,]),
                                   my_primary,
                                   my_secondary)
  sim <- growth$simulation
  all_simulations[[i]] <- sim
}

final_conc <- do.call(rbind, all_simulations)
df <- final_conc

PS00399_df_day14 = subset(df, time == "14")
PS00399_df_day21 = subset(df, time == "21")
PS00399_df_day35 = subset(df, time == "35")

sum(PS00399_df_day35$logN>5)/10000
sum(PS00399_df_day35$logN<3)/10000
sum(PS00399_df_day35$logN>=3 & PS00399_df_day35$logN<=5)/10000

PS00399_df_day35$color<-ifelse(test = PS00399_df_day35$logN>=5,yes = "Above 5 log",no = 
                    ifelse(PS00399_df_day35$logN>=3,yes = "Between 3 and 5 log",no = "Below 3 log"))

PS00399_df_day35$isolate = "PS00399"
cytotoxicus_df_day35 = rbind(PS00399_df_day35,PS00399_df_day35)

cytotoxicus_df_day35$color <- factor(
  cytotoxicus_df_day35$color,
  levels = c("Below 3 log", "Between 3 and 5 log", "Above 5 log")
)

breaks <- seq(0, 10, by = 0.5)
finalhist <- ggplot(cytotoxicus_df_day35, aes(x = logN, fill = color)) +
  geom_histogram(
    breaks = breaks,
    color = NA
    ) +
  facet_wrap(~ isolate, ncol = 3) +
  scale_fill_manual(
    name = expression(italic(B~cereus) ~ "count per ml"),
    values = c(
      "Below 3 log" = "springgreen3",
      "Between 3 and 5 log" = "darkorange1",
      "Above 5 log" = "red3"
      ),
    breaks = c("Below 3 log", "Between 3 and 5 log", "Above 5 log")
    ) +
  xlab("log CFU per ml") +
  ylab("Number of Units (log scale)") +
  ggtitle(expression(italic(cytotoxicus))) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
  theme_classic() +   # removes background grid lines
  theme(
    plot.title = element_text(size = 28, face = "bold"),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    strip.text = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 22, face = "bold")
    )
cytotoxicus_d35_plot = finalhist

# Supplemental Figure 1
library(gridExtra)
grid.arrange(pseudomycoides_d21_plot,albus_d21_plot,mobilis_d21_plot)
grid.arrange(tropicus_d21_plot,pacificus_d21_plot,cereus_d21_plot)
grid.arrange(thuringiensis_d21_plot,toyonensis_d21_plot,cytotoxicus_d21_plot)

grid.arrange(pseudomycoides_d35_plot,albus_d35_plot,mobilis_d35_plot)
grid.arrange(tropicus_d35_plot,pacificus_d35_plot,cereus_d35_plot)
grid.arrange(thuringiensis_d35_plot,toyonensis_d35_plot,cytotoxicus_d35_plot)

# Figure 4
dat = read.csv("Figure 4 input.csv")
dat$Species <- factor(dat$Species, levels = unique(dat$Species))
dat$Day <- factor(dat$Day)
p <- ggplot(dat, aes(x = Species, y = log_5plus, fill = Day)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_continuous(limits = c(0, 5), expand = c(0, 0)) +
  labs(x = "Species",
       y = "Percent milk containers over 5 logs",
       fill = "Shelf-life day") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

p
