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

# Parallel setup
library(furrr)
library(future)
plan(multisession, workers = 2)

# Load utility functions
source("UtilityFunctions_dynamic_growth.R")

# Generate database 
# BTyper data
BTyper3_input = read.csv("Btyper3_Results.csv")
colnames(BTyper3_input)[1] <- "Isolate.Name"
gp_input = read.csv("simulation_input.csv")
database = cbind(BTyper3_input,gp_input[,3:10])
database <- database %>% 
  separate(Closest_Type_Strain.ANI., into = c("species","ANI"), sep = "\\(") %>%
  mutate(ANI = gsub("\\)", "", ANI))

# Input BRiskTyper intermediate outputs for all 22 B cereus group isolates
N0_df <- read.csv("Raw_data_validation.csv")
N0_df_sub <- subset(N0_df, Consumer.storage.day == 0) # just to get the order of isolates

cases <- data.frame(
  id = N0_df_sub$Isolate
) %>%
  mutate(
    final_file = sprintf(
      "Genomic data for detected isolates/%s_contigs_final_results.csv",
      id
    ),
    ani_file = sprintf(
      "Genomic data for detected isolates/%s_contigs_cytotoxicity_fastani.csv",
      id
    )
  )

# Parallel simulation function 
run_simulation <- function(final_file, ani_file, isolate_id) {

# Input BRiskTyper information on the type strain closest to the input strains  
df <- read.csv(final_file)
df <- df %>% 
  separate(Closest_Type_Strain.ANI., into = c("species","ANI"), sep = "\\(") %>%
  separate(Adjusted_panC_Group.predicted_species., into = c("panC_Group","predicted_species"), sep = "\\(") %>%
  mutate(ANI = gsub("\\)", "", ANI),
         panC_Group = gsub("\\)", "", panC_Group),
         predicted_species = gsub("\\)", "", predicted_species))

# Input BRiskTyper ANI similarity matrix
ANI_file <- read.csv(ani_file, header = FALSE)
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

set.seed(1)

# Simulation setup
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

# Sensitivity analysis to assess uncertainty of temperature 
# + 1.5C
# ModelData$T_F = ModelData$T_F + 1.5
# ModelData$T_T = ModelData$T_T + 1.5
# ModelData$T_S = ModelData$T_S + 1.5
# ModelData$T_T2 = ModelData$T_T2 + 1.5
# ModelData$T_H = ModelData$T_H + 1.5 

# - 1.5C
# ModelData$T_F = ModelData$T_F - 1.5
# ModelData$T_T = ModelData$T_T - 1.5
# ModelData$T_S = ModelData$T_S - 1.5
# ModelData$T_T2 = ModelData$T_T2 - 1.5
# ModelData$T_H = ModelData$T_H - 1.5 

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

## Assign growth parameters to 10,000 units of HTST milk 
ModelData$index = match(ModelData$isolate, matching_species_df$Isolate.Name)
ModelData$mean_LOG10Q0 = matching_species_df$mean_LOG10Q0[ModelData$index]
ModelData$sd_LOG10Q0 = matching_species_df$sd_LOG10Q0[ModelData$index]
ModelData <- ModelData %>%
  mutate(LOGQ0 = rnorm(n(),mean = mean_LOG10Q0,sd = sd_LOG10Q0))
ModelData$Q0 = 10^ModelData$LOGQ0
ModelData$mean_Nmax = matching_species_df$mean_Nmax[ModelData$index]
ModelData$sd_Nmax = matching_species_df$sd_Nmax[ModelData$index]
ModelData <- ModelData %>%
  mutate(LOGNmax = rnorm(n(),mean = mean_Nmax,sd = sd_Nmax))
ModelData$Nmax = 10^(ModelData$LOGNmax)
ModelData$b = matching_species_df$b[ModelData$index]
ModelData$Tmin = matching_species_df$Tmin[ModelData$index]
ModelData$Clade = matching_species_df$Clade[ModelData$index]

## Generate N0 from a Poisson distribution 
set.seed(42)
N0 = rpois(n = n_sim, lambda = 1*1900)
# Sensitivity analysis to assess uncertainty of lambda 
# N0 = rpois(n = n_sim, lambda = 0.5*1900) # lower level
# N0 = rpois(n = n_sim, lambda = 10*1900) # higher level
ModelData$N0 = N0/1900

ModelData$Topt = sapply(ModelData$Clade, xopt_func)
ModelData$mu_opt = (ModelData$b*(ModelData$Topt-ModelData$Tmin))^2

# Parallel unit simulation
my_times <- seq(0,35)
num_iterations <- nrow(ModelData)
all_simulations <- future_map_dfr(seq_len(num_iterations), function(i) {
  
  # Primary parameters
  primary <- list(
    mu_opt = ModelData$mu_opt[i],
    Nmax   = ModelData$Nmax[i],
    N0     = ModelData$N0[i],
    Q0     = ModelData$Q0[i]
  )
  
  # Secondary parameters
  secondary <- list(
    temperature = list(
      model = "reducedRatkowsky",
      xmin  = ModelData$Tmin[i],
      b     = ModelData$b[i],
      clade = ModelData$Clade[i]
    )
  )
  
  growth <- predict_dynamic_growth(
    times = my_times,
    env_conditions = tibble(
      time = env_cond_time[i, ],
      temperature = env_cond_temp[i, ]
    ),
    primary,
    secondary
  )
  
  sim <- growth$simulation
  sim$isolate_id <- isolate_id
  sim
},
.options = furrr_options(seed = TRUE)
)
return(all_simulations)
}

results_list <- future_pmap(
  cases,
  function(id, final_file, ani_file) {
    run_simulation(final_file = final_file, ani_file = ani_file,isolate_id = id)
    },
  .options = furrr_options(seed = TRUE)
)

final_results <- bind_rows(results_list) 

final_results_day14 = subset(final_results, time == "14")
final_results_day21 = subset(final_results, time == "21")
final_results_day35 = subset(final_results, time == "35")

# Figure 6
summary_df_14 <- final_results_day14 %>%
  group_by(isolate_id) %>%
  summarise(
    n_total = n(),
    n_gt_5 = sum(logN > 5),
    n_lt_3 = sum(logN < 3),
    n_3_5  = sum(logN >= 3 & logN <= 5),
    pct_gt_5 = n_gt_5 / n_total * 100,
    pct_lt_3 = n_lt_3 / n_total * 100,
    pct_3_5  = n_3_5  / n_total * 100
  )
summary_df_14 <- N0_df_sub %>%
  select(Isolate, Closest.Type.Strain) %>%
  left_join(summary_df_14, by = c("Isolate" = "isolate_id"))

summary_df_21 <- final_results_day21 %>%
  group_by(isolate_id) %>%
  summarise(
    n_total = n(),
    n_gt_5 = sum(logN > 5),
    n_lt_3 = sum(logN < 3),
    n_3_5  = sum(logN >= 3 & logN <= 5),
    pct_gt_5 = n_gt_5 / n_total * 100,
    pct_lt_3 = n_lt_3 / n_total * 100,
    pct_3_5  = n_3_5  / n_total * 100
  )
summary_df_21 <- N0_df_sub %>%
  select(Isolate, Closest.Type.Strain) %>%
  left_join(summary_df_21, by = c("Isolate" = "isolate_id"))

final_results_day35 <- final_results_day35[!is.na(final_results_day35$logN), ]
summary_df_35 <- final_results_day35 %>%
  group_by(isolate_id) %>%
  summarise(
    n_total = n(),
    n_gt_5 = sum(logN > 5),
    n_lt_3 = sum(logN < 3),
    n_3_5  = sum(logN >= 3 & logN <= 5),
    pct_gt_5 = n_gt_5 / n_total * 100,
    pct_lt_3 = n_lt_3 / n_total * 100,
    pct_3_5  = n_3_5  / n_total * 100
  )
summary_df_35 <- N0_df_sub %>%
  select(Isolate, Closest.Type.Strain) %>%
  left_join(summary_df_35, by = c("Isolate" = "isolate_id"))

dat <- list(
  "14" = summary_df_14,
  "21" = summary_df_21,
  "35" = summary_df_35
) %>%
  bind_rows(.id = "Day") %>%
  mutate(
    Day = as.numeric(Day)
  ) %>%
  select(Isolate, Closest.Type.Strain, pct_gt_5, Day)

colnames(dat) <- c("Isolate", "species", "log_5plus", "Day")

dat <- dat %>%
  group_by(Day) %>%
  distinct(species, .keep_all = TRUE) %>%
  ungroup()

dat$species <- factor(dat$species, levels = unique(dat$species))
dat$Day <- factor(dat$Day)
p <- ggplot(dat, aes(x = species, y = log_5plus, fill = Day)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0)) +
  scale_x_discrete(
    labels = function(x) parse(text = paste0("italic('B.')~italic('", x, "')"))
  ) +
  labs(x = "Species",
       y = "Percent milk containers over 5 logs",
       fill = "Shelf-life day") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
p
ggsave("Figure 6.png", p, width = 12, height = 6, dpi = 600)

# Supplemental Figure 3
final_results_day14 <- N0_df_sub %>%
  select(Isolate, Closest.Type.Strain) %>%
  left_join(final_results_day14, by = c("Isolate" = "isolate_id"))

final_results_day21 <- N0_df_sub %>%
  select(Isolate, Closest.Type.Strain) %>%
  left_join(final_results_day21, by = c("Isolate" = "isolate_id"))

final_results_day35 <- N0_df_sub %>%
  select(Isolate, Closest.Type.Strain) %>%
  left_join(final_results_day35, by = c("Isolate" = "isolate_id"))

final_results_day14<- final_results_day14[final_results_day14$logN != -Inf, ]
final_results_day14$color<-ifelse(test = final_results_day14$logN>=5,yes = "Above 5 log",no = 
                    ifelse(final_results_day14$logN>=3,yes = "Between 3 and 5 log",no = "Below 3 log"))

final_results_day21<- final_results_day21[final_results_day21$logN != -Inf, ]
final_results_day21$color<-ifelse(test = final_results_day21$logN>=5,yes = "Above 5 log",no = 
                    ifelse(final_results_day21$logN>=3,yes = "Between 3 and 5 log",no = "Below 3 log"))

final_results_day35<- final_results_day35[final_results_day35$logN != -Inf, ]
final_results_day35$color<-ifelse(test = final_results_day35$logN>=5,yes = "Above 5 log",no = 
                    ifelse(final_results_day35$logN>=3,yes = "Between 3 and 5 log",no = "Below 3 log"))

final_results_day14$color <- factor(
  final_results_day14$color,
  levels = c("Below 3 log", "Between 3 and 5 log", "Above 5 log")
)
breaks <- seq(0, 10, by = 1)
plot_list_day_14 <- split(final_results_day14, final_results_day14$Closest.Type.Strain)
plot_list_day_14 <- lapply(plot_list_day_14, function(data_subset) {
  ggplot(data_subset, aes(x = logN, fill = color)) +
    geom_histogram(
      breaks = breaks,
      color = NA
    ) +
    facet_wrap(~ Isolate, ncol = 3) +
    scale_fill_manual(
      name = expression(italic(italic(B~cereus) ~ "count per ml")),
      values = c(
        "Below 3 log" = "#779ECC",
        "Between 3 and 5 log" = "#F2C894",
        "Above 5 log" = "#FF985A"
      ),
      breaks = c("Below 3 log", "Between 3 and 5 log", "Above 5 log")
    ) +
    xlab("log CFU per ml") +
    ylab("Number of Units (log scale)") +
    ggtitle(bquote(italic("B.") ~ italic(.(unique(data_subset$Closest.Type.Strain)))))+
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 28, face = "bold"),
      axis.title.x = element_text(size = 24),
      axis.title.y = element_text(size = 24),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 22, face = "bold")
    )
})

final_results_day21$color <- factor(
  final_results_day21$color,
  levels = c("Below 3 log", "Between 3 and 5 log", "Above 5 log")
)
breaks <- seq(0, 10, by = 1)
plot_list_day_21 <- split(final_results_day21, final_results_day21$Closest.Type.Strain)
plot_list_day_21 <- lapply(plot_list_day_21, function(data_subset) {
  ggplot(data_subset, aes(x = logN, fill = color)) +
    geom_histogram(
      breaks = breaks,
      color = NA
    ) +
    facet_wrap(~ Isolate, ncol = 3) +
    scale_fill_manual(
      name = expression(italic(italic(B~cereus) ~ "count per ml")),
      values = c(
        "Below 3 log" = "#779ECC",
        "Between 3 and 5 log" = "#F2C894",
        "Above 5 log" = "#FF985A"
      ),
      breaks = c("Below 3 log", "Between 3 and 5 log", "Above 5 log")
    ) +
    xlab("log CFU per ml") +
    ylab("Number of Units (log scale)") +
    ggtitle(bquote(italic("B.") ~ italic(.(unique(data_subset$Closest.Type.Strain)))))+
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 28, face = "bold"),
      axis.title.x = element_text(size = 24),
      axis.title.y = element_text(size = 24),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 22, face = "bold")
    )
})

final_results_day35$color <- factor(
  final_results_day35$color,
  levels = c("Below 3 log", "Between 3 and 5 log", "Above 5 log")
)
breaks <- seq(0, 10, by = 1)
plot_list_day_35 <- split(final_results_day35, final_results_day35$Closest.Type.Strain)
plot_list_day_35 <- lapply(plot_list_day_35, function(data_subset) {
  ggplot(data_subset, aes(x = logN, fill = color)) +
    geom_histogram(
      breaks = breaks,
      color = NA
    ) +
    facet_wrap(~ Isolate, ncol = 3) +
    scale_fill_manual(
      name = expression(italic(italic(B~cereus) ~ "count per ml")),
      values = c(
        "Below 3 log" = "#779ECC",
        "Between 3 and 5 log" = "#F2C894",
        "Above 5 log" = "#FF985A"
      ),
      breaks = c("Below 3 log", "Between 3 and 5 log", "Above 5 log")
    ) +
    xlab("log CFU per ml") +
    ylab("Number of Units (log scale)") +
    ggtitle(bquote(italic("B.") ~ italic(.(unique(data_subset$Closest.Type.Strain))))) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 28, face = "bold"),
      axis.title.x = element_text(size = 24),
      axis.title.y = element_text(size = 24),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 22, face = "bold")
    )
})


# Supplemental Figure 3
library(gridExtra)

S3A1 = grid.arrange(plot_list_day_14$pseudomycoides,plot_list_day_14$albus, plot_list_day_14$mobilis)
ggsave("Supplemental 3A1.png", S3A1, width = 12, height = 14, dpi = 600)

S3A2 = grid.arrange(plot_list_day_14$tropicus,plot_list_day_14$pacificus, plot_list_day_14$cereus)
ggsave("Supplemental 3A2.png", S3A2, width = 12, height = 14, dpi = 600)

S3A3 = grid.arrange(plot_list_day_14$thuringiensis,plot_list_day_14$toyonensis, plot_list_day_14$cytotoxicus)
ggsave("Supplemental 3A3.png", S3A3, width = 12, height = 14, dpi = 600)

S3B1 = grid.arrange(plot_list_day_21$pseudomycoides,plot_list_day_21$albus, plot_list_day_21$mobilis)
ggsave("Supplemental 3B1.png", S3B1, width = 12, height = 14, dpi = 600)

S3B2 = grid.arrange(plot_list_day_21$tropicus,plot_list_day_21$pacificus, plot_list_day_21$cereus)
ggsave("Supplemental 3B2.png", S3B2, width = 12, height = 14, dpi = 600)

S3B3 = grid.arrange(plot_list_day_21$thuringiensis,plot_list_day_21$toyonensis, plot_list_day_21$cytotoxicus)
ggsave("Supplemental 3B3.png", S3B3, width = 12, height = 14, dpi = 600)

S3C1 = grid.arrange(plot_list_day_35$pseudomycoides,plot_list_day_35$albus, plot_list_day_35$mobilis)
ggsave("Supplemental 3C1.png", S3C1, width = 12, height = 14, dpi = 600)

S3C2 = grid.arrange(plot_list_day_35$tropicus,plot_list_day_35$pacificus, plot_list_day_35$cereus)
ggsave("Supplemental 3C2.png", S3C2, width = 12, height = 14, dpi = 600)

S3C3 = grid.arrange(plot_list_day_35$thuringiensis,plot_list_day_35$toyonensis, plot_list_day_35$cytotoxicus)
ggsave("Supplemental 3C3.png", S3C3, width = 12, height = 14, dpi = 600)

# Figure 7
# Uncertainty of N0
result_df_N0 <- summary_df_35_1900 %>%
  select(Isolate, Closest.Type.Strain, pct_gt_5) %>%
  rename(pct_gt_5_1900 = pct_gt_5) %>%
  left_join(
    summary_df_35_950 %>%
      select(Isolate, Closest.Type.Strain, pct_gt_5) %>%
      rename(pct_gt_5_950 = pct_gt_5),
    by = c("Isolate", "Closest.Type.Strain")
  ) %>%
  left_join(
    summary_df_35_19000 %>%
      select(Isolate, Closest.Type.Strain, pct_gt_5) %>%
      rename(pct_gt_5_19000 = pct_gt_5),
    by = c("Isolate", "Closest.Type.Strain")
  )

tornado_dfN0 <- result_df_N0 %>%
  group_by(Closest.Type.Strain) %>%
  summarise(
    pct_gt_5_950   = mean(pct_gt_5_950, na.rm = TRUE),
    pct_gt_5_1900  = mean(pct_gt_5_1900, na.rm = TRUE),
    pct_gt_5_19000 = mean(pct_gt_5_19000, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    decrease_950   = pct_gt_5_950 - pct_gt_5_1900,
    increase_19000 = pct_gt_5_19000 - pct_gt_5_1900
  )

order_levels <- tornado_dfN0 %>%
  arrange(increase_19000) %>%
  pull(Closest.Type.Strain) %>%
  unique()

tornado_dfN0 <- tornado_dfN0 %>%
  mutate(Closest.Type.Strain = factor(Closest.Type.Strain, levels = order_levels))

plot_dfN0 <- tornado_dfN0 %>%
  select(Closest.Type.Strain, decrease_950, increase_19000) %>%
  pivot_longer(
    cols = c(decrease_950, increase_19000),
    names_to = "scenario",
    values_to = "change"
  ) %>%
  mutate(
    scenario = recode(
      scenario,
      decrease_950 = "lambda = 950 CFU/container",
      increase_19000 = "lambda = 19000 CFU/container"
    )
  )

Figure_7A = ggplot(plot_dfN0,
       aes(
         y = Closest.Type.Strain,
         x = change,
         fill = scenario
       )) +
  geom_col(position = "identity") +
  scale_fill_manual(
    values = c(
      "lambda = 950 CFU/container" = "skyblue",
      "lambda = 19000 CFU/container" = "#C2185B"
    )
  ) +
  scale_y_discrete(
    labels = function(x) parse(text = paste0("italic('B.')~italic('", x, "')"))
  ) +
  labs(
    x = "Percentage point increases and decreases from the baseline",
    y = "Species",
    fill = NULL
  ) +
  theme_bw() +
  coord_cartesian(xlim = c(-1.5, 4))
ggsave("Figure 7A.png", Figure_7A, width = 12, height = 8, dpi = 600)

# Uncertainty of temperature profile
result_df_temp <- summary_df_35_1900 %>%
  select(Isolate, Closest.Type.Strain, pct_gt_5) %>%
  rename(pct_gt_5_base = pct_gt_5) %>%
  left_join(
    summary_df_35_minus1.5C %>%
      select(Isolate, Closest.Type.Strain, pct_gt_5) %>%
      rename(pct_gt_5_minus1.5C = pct_gt_5),
    by = c("Isolate", "Closest.Type.Strain")
  ) %>%
  left_join(
    summary_df_35_1.5C %>%
      select(Isolate, Closest.Type.Strain, pct_gt_5) %>%
      rename(pct_gt_5_1.5C = pct_gt_5),
    by = c("Isolate", "Closest.Type.Strain")
  )

tornado_dfTemp <- result_df_temp %>%
  group_by(Closest.Type.Strain) %>%
  summarise(
    pct_gt_5_base        = mean(pct_gt_5_base, na.rm = TRUE),
    pct_gt_5_minus1.5C   = mean(pct_gt_5_minus1.5C, na.rm = TRUE),
    pct_gt_5_1.5C        = mean(pct_gt_5_1.5C, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    decrease_minus1.5C = pct_gt_5_minus1.5C - pct_gt_5_base,
    increase_1.5C      = pct_gt_5_1.5C - pct_gt_5_base
  )

order_levels <- tornado_dfTemp %>%
  arrange(increase_1.5C) %>%
  pull(Closest.Type.Strain) %>%
  unique()

tornado_dfTemp <- tornado_dfTemp %>%
  mutate(Closest.Type.Strain = factor(Closest.Type.Strain, levels = order_levels))

plot_dfTemp <- tornado_dfTemp %>%
  select(Closest.Type.Strain, decrease_minus1.5C, increase_1.5C) %>%
  pivot_longer(
    cols = c(decrease_minus1.5C, increase_1.5C),
    names_to = "scenario",
    values_to = "change"
  ) %>%
  mutate(
    scenario = recode(
      scenario,
      decrease_minus1.5C = "-1.5C",
      increase_1.5C = "+1.5C"
    )
  )

Figure_7B = ggplot(plot_dfTemp,
       aes(
         y = Closest.Type.Strain,
         x = change,
         fill = scenario
       )) +
  geom_col(position = "identity") +
  scale_fill_manual(
    values = c(
      "-1.5C" = "skyblue",
      "+1.5C" = "#C2185B"
    )
  ) +
  scale_y_discrete(
    labels = function(x) parse(text = paste0("italic('B.')~italic('", x, "')"))
  ) +
  labs(
    x = "Percentage point change from baseline",
    y = "Species",
    fill = NULL
  ) +
  theme_bw() +
  coord_cartesian(xlim = c(-1.5, 4))
ggsave("Figure 7B.png", Figure_7B, width = 12, height = 8, dpi = 600)

# Figure 3
library(flextable)
tbl_core <- as.data.frame(
  matrix(
    c(
      "Negligible", "High", "High", "Very high", "Very high",
      "Negligible", "Medium", "High", "High", "Very high",
      "Negligible", "Medium", "Medium", "High", "Very high",
      "Negligible", "Low", "Medium", "Medium", "High",
      "Negligible", "Negligible", "Negligible", "Negligible", "Negligible"
    ),
    nrow = 5,
    byrow = TRUE
  )
)

col_headers <- c("", "Negligible", "Low", "Medium", "High", "Very high")

tbl_with_header <- rbind(
  col_headers,
  cbind(
    Severity = c("Very high", "High", "Medium", "Low", "Negligible"),
    tbl_core
  )
)

tbl_with_header <- as.data.frame(tbl_with_header, stringsAsFactors = FALSE)

colors <- c(
  "Very high"  = "#FF985A",
  "High"       = "#FFB347",
  "Medium"     = "#F2C894",
  "Low"        = "#9FC0CE",
  "Negligible" = "#779ECC"
)

ft <- flextable(tbl_with_header)

ft <- bg(ft, i = 1, bg = "white")   # top row
ft <- bg(ft, j = 1, bg = "white")   # left column

for (i in 2:nrow(tbl_with_header)) {
  for (j in 2:ncol(tbl_with_header)) {
    
    val <- tbl_with_header[i, j]
    
    ft <- bg(
      ft,
      i = i,
      j = j,
      bg = colors[val]
    )
  }
}

ft <- align(ft, align = "center", part = "all")
ft
save_as_image(
  ft,
  path = "Figure_3A.png",
  expand = 10
)

tbl_core <- as.data.frame(
  matrix(
    c(
      "Medium", "High", "High", "Very high", "Very high",
      "Medium", "Medium", "High", "High", "Very high",
      "Low", "Medium", "Medium", "High", "Very high",
      "Low", "Low", "Medium", "Medium", "High",
      "Negligible", "Negligible", "Negligible", "Negligible", "Negligible"
    ),
    nrow = 5,
    byrow = TRUE
  )
)

col_headers <- c("", "Very low", "Low", "Medium", "High", "Very high")

tbl_with_header <- rbind(
  col_headers,
  cbind(
    Severity = c("Very high", "High", "Medium", "Low", "Negligible"),
    tbl_core
  )
)

tbl_with_header <- as.data.frame(tbl_with_header, stringsAsFactors = FALSE)

colors <- c(
  "Very high"  = "#FF985A",
  "High"       = "#FFB347",
  "Medium"     = "#F2C894",
  "Low"        = "#9FC0CE",
  "Negligible" = "#779ECC"
)

ft1 <- flextable(tbl_with_header)

ft1 <- bg(ft1, i = 1, bg = "white")  # top row
ft1 <- bg(ft1, j = 1, bg = "white")  # left column

for (i in 2:nrow(tbl_with_header)) {
  for (j in 2:ncol(tbl_with_header)) {
    
    val <- tbl_with_header[i, j]
    
    ft1 <- bg(
      ft1,
      i = i,
      j = j,
      bg = colors[val]
    )
  }
}

ft1 <- align(ft1, align = "center", part = "all")
ft1
save_as_image(
  ft1,
  path = "Figure_3B.png",
  expand = 10
)
