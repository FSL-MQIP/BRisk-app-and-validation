# Load packages
library(tidyverse)
library(tibble)
library(EnvStats)
library(truncnorm)
library(jmuOutlier)
library(formula.tools)
library(purrr)
library(deSolve)
library(rlang)
library(ggplot2)

# Parallel setup
library(furrr)
library(future)
plan(multisession, workers = 2)

# Utility functions
source("UtilityFunctions_dynamic_growth.R")

# Generate database 
# Temperature profile data set
temp_data = read.csv("Validation Study Temp.csv", header = FALSE)

# BTyper data
BTyper3_input <- read.csv("Btyper3_Results.csv")
colnames(BTyper3_input)[1] <- "Isolate.Name"
gp_input <- read.csv("simulation_input.csv")
database <- cbind(BTyper3_input, gp_input[,3:10])
database <- database %>% 
  separate(Closest_Type_Strain.ANI., into = c("species","ANI"), sep = "\\(") %>%
  mutate(ANI = gsub("\\)", "", ANI))

# Input BRiskTyper intermediate outputs, N0 for all 22 B cereus group isolates, consumer storage days 14, 21, 35
# Organize N0 file 
N0_df <- read.csv("Raw_data_validation.csv")
N0_df_sub <- subset(N0_df, Consumer.storage.day == 0)
N0_df_sub$Count <- 10^N0_df_sub$Log.Count

# Define consumer storage time 
tH_values <- c(14, 21, 35)

# Combine BRiskTyper intermediate outputs, N0, tH in the same file
cases <- expand_grid(
  id = N0_df_sub$Isolate,
  t_H = tH_values) 

cases <- cases %>%
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

cases <- cases %>%
  left_join(
    N0_df_sub %>% select(Isolate, Count),
    by = c("id" = "Isolate")
  ) %>%
  rename(N0 = Count)

# Parallel simulation function 
run_simulation <- function(final_file, ani_file, isolate_id, N0_value, t_H) {
  
  # Input BRiskTyper information on the type strain closest to the input strains 
  df <- read.csv(final_file)
  df <- df %>% 
    separate(Closest_Type_Strain.ANI., into = c("species","ANI"), sep = "\\(") %>%
    separate(Adjusted_panC_Group.predicted_species., into = c("panC_Group","predicted_species"), sep = "\\(") %>%
    mutate(
      ANI = gsub("\\)", "", ANI),
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
  
  # Filter database input for rows with the same species as the BRiskTyper input
  df$species <- trimws(df$species)
  database$ANI_new <- ANI_file$ANI
  matching_species_df <- subset(database, species == df$species)
  
  # Assigning ANI weight
  matching_species_df$ANI_wght <- matching_species_df$ANI_new / sum(matching_species_df$ANI_new)
  
  set.seed(1)
  
  # Simulation setup
  n_sim <- 10000
  matching_species_df$n_units <- round(n_sim * matching_species_df$ANI_wght)
  sampled_isolates <- rep(matching_species_df$Isolate.Name, matching_species_df$n_units)
  if(length(sampled_isolates) != n_sim){
    sampled_isolates <- sample(sampled_isolates, n_sim)
  }
  
  ModelData <- data.frame(
    unit_id = seq_len(n_sim),
    isolate = sampled_isolates
  )
  
  ModelData$N0 <- rpois(n_sim, lambda = N0_value)
  
  # Validation temp profiles 
  # Stage 1: facility storage 
  ## (a)  Sample the temperature distribution
  ModelData$T_F <- sample(temp_data[1:1820,], size = 10000, replace = TRUE)
  
  ## (b) Sample the storage time (in days) distribution
  ModelData$t_F = 1.5
  
  # Stage 2: transport from facility to retail store
  ## (a)  Sample the temperature distribution
  ModelData$T_T <- sample(temp_data[1821:8554,], size = 10000, replace = TRUE)
  
  ## (b) Sample the transportation time (in days) distribution
  ModelData$t_T = 5
  
  # Stage 3: storage/display at retail store
  ## (a)  Sample the temperature distribution
  ModelData$T_S <- sample(temp_data[8555:11801,], size = 10000, replace = TRUE)
  
  ## (b) Sample the storage time (in days) distribution
  ModelData$t_S = 2
  
  ## Stage 4: transportation from retail store to home
  ## (a)  Sample the temperature distribution
  ModelData$T_T2 <- sample(temp_data[11802:11862,], size = 10000, replace = TRUE)
  
  ## (b) Sample the transportation time (in days) distribution 
  ModelData$t_T2 = 1/24
  
  ## Stage 5: home storage 
  ## (a)  Sample the temperature distribution
  ModelData$T_H <- sample(temp_data[11863:61844,], size = 10000, replace = TRUE)
  
  ## (b) Define t_H (in days) as consumer home storage days
  ModelData$t_H <- rep(t_H, each = n_sim)
  
  ## Model temperature profiles of 10,000 units HTST milk 
   env_cond_time <- matrix(c(
    rep(0, n_sim),
    ModelData$t_F,
    ModelData$t_F + 0.001,
    ModelData$t_F + ModelData$t_T,
    ModelData$t_F + ModelData$t_T + 0.001,
    ModelData$t_F + ModelData$t_T + ModelData$t_S,
    ModelData$t_F + ModelData$t_T + ModelData$t_S + 0.001,
    ModelData$t_F + ModelData$t_T + ModelData$t_S + ModelData$t_T2,
    ModelData$t_F + ModelData$t_T + ModelData$t_S + ModelData$t_T2 + 0.001,
    ModelData$t_F + ModelData$t_T + ModelData$t_S + ModelData$t_T2 + ModelData$t_H), ncol = 10)
  
  env_cond_temp <- matrix(c(
    ModelData$T_F,
    ModelData$T_F,
    ModelData$T_T,
    ModelData$T_T,
    ModelData$T_S,
    ModelData$T_S,
    ModelData$T_T2,
    ModelData$T_T2,
    ModelData$T_H,
    ModelData$T_H), ncol = 10)
  
  # Assign growth parameters to 10,000 units of HTST milk 
  ModelData$index <- match(ModelData$isolate, matching_species_df$Isolate.Name)
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
  ModelData$b     <- matching_species_df$b[ModelData$index]
  ModelData$Tmin  <- matching_species_df$Tmin[ModelData$index]
  ModelData$Clade <- matching_species_df$Clade[ModelData$index]
  
  ModelData$Topt <- sapply(ModelData$Clade, xopt_func)
  ModelData$mu_opt <- (ModelData$b * (ModelData$Topt - ModelData$Tmin))^2
  
  # Parallel unit simulation
  conc <- future_map_dbl(seq_len(n_sim), function(i) {
    
    my_primary <- list(
      mu_opt = ModelData$mu_opt[i],
      Nmax   = ModelData$Nmax[i],
      N0     = ModelData$N0[i],
      Q0     = ModelData$Q0[i]
    )
    
    my_secondary <- list(
      temperature = list(
        model = "reducedRatkowsky",
        xmin  = ModelData$Tmin[i],
        b     = ModelData$b[i],
        clade = ModelData$Clade[i]
      )
    )
    
    growth <- predict_dynamic_growth(
      times = env_cond_time[i, ],
      env_conditions = data.frame(
        time = env_cond_time[i, ],
        temperature = env_cond_temp[i, ]
      ),
      my_primary,
      my_secondary
    )
    
    tail(growth$simulation$logN, 1)
    
  }, .options = furrr_options(seed = TRUE))
  
  ModelData$conc <- conc
  ModelData$isolate_id <- isolate_id
  
  return(ModelData)
}

# Run all isolates (parallel)
results_list <- future_pmap(
  cases,
  function(id, t_H, final_file, ani_file, N0) {
    run_simulation(final_file, ani_file, id, N0, t_H)
  },
  .options = furrr_options(seed = TRUE)
)

# Combine results
final_results <- bind_rows(results_list)

iso_order <- c(
  "PS00134","PS00136","PS00123",
  "PS00161","PS00423","PS00048",
  "PS00906","PS00515",
  "PS00047","PS00452",
  "PS00436","PS00485","PS00513",
  "PS00087",
  "PS00107","PS00512","PS00597",
  "PS00598","PS00565","PS00809",
  "PS00539","PS00399"
)

# Summary stats
summary_stats = final_results %>%
  filter(t_H %in% c(14, 21, 35)) %>%
  group_by(isolate_id, t_H) %>%
  summarise(
    median_conc = median(conc, na.rm = TRUE),
    Q1_conc = quantile(conc, 0.25, na.rm = TRUE),
    Q3_conc = quantile(conc, 0.75, na.rm = TRUE),
    .groups = "drop"
  )%>%
  mutate(isolate_id = factor(isolate_id, levels = iso_order)) %>%
  arrange(isolate_id, t_H)

N0_df <- N0_df %>%
  filter(Consumer.storage.day != "0")

joined_summary_stats <- N0_df %>%
  mutate(Consumer.storage.day = as.numeric(Consumer.storage.day)) %>%
  left_join(
    summary_stats %>% mutate(t_H = as.numeric(t_H)),
    by = c(
      "Isolate" = "isolate_id",
      "Consumer.storage.day" = "t_H"
    )
  ) %>%
  mutate(Isolate = factor(Isolate, levels = iso_order)) %>%
  arrange(Isolate, Consumer.storage.day)

colnames(joined_summary_stats) = c("isolate","species","consumer storage day","observed count","predicted median","predicted Q1", "predicted Q3")
write.csv(joined_summary_stats,"validation_stat_summary.csv")

# Sensitivity analysis 
library(sensitivity)

input_vars <- c("N0","T_F","T_T","T_S","T_T2","T_H","Q0","Nmax")
output_var <- "conc"

prcc_results <- final_results %>%
  group_by(isolate_id, t_H) %>%
  group_modify(~{
    X <- .x[, input_vars]
    Y <- .x[[output_var]]
    # Run PRCC
    prcc <- pcc(
      X = X,
      y = Y,
      rank = TRUE
    )
    tibble(
      parameter = rownames(prcc$PRCC),
      PRCC = prcc$PRCC[,1],
      p_value = prcc$p.value[,1]
    )
  }) %>%
  ungroup()

# Supplemental Figure 2
# Filter for Q0 only
df_q0 <- prcc_results %>%
  filter(parameter == "Q0")
# Ensure ordering (optional but recommended)
df_q0$t_H <- factor(df_q0$t_H, levels = c(14, 21, 35))

mean_prcc_by_time <- df_q0 %>%
  group_by(t_H) %>%
  summarise(mean_PRCC = mean(PRCC, na.rm = TRUE))

df_Nmax <- prcc_results %>%
  filter(parameter == "Nmax")
# Ensure ordering (optional but recommended)
df_Nmax$t_H <- factor(df_Nmax$t_H, levels = c(14, 21, 35))

mean_prcc_by_time_Nmax <- df_Nmax %>%
  group_by(t_H) %>%
  summarise(mean_PRCC = mean(PRCC, na.rm = TRUE))


mean_prcc_by_time <- df_q0 %>%
  group_by(t_H) %>%
  summarise(mean_PRCC = mean(PRCC, na.rm = TRUE))

# Heatmap
prcc_results$parameter <- factor(
  prcc_results$parameter,
  levels = c("Q0", "Nmax", "N0", "T_H", "T_F", "T_S", "T_T", "T_T2")
)
heat_map <- ggplot(prcc_results,
                   aes(x = parameter,
                       y = isolate_id,
                       fill = PRCC)) +
  geom_tile(color = "white") +
  facet_wrap(~ t_H, nrow = 1,
             labeller = labeller(
               t_H = function(x)
                 paste0("Consumer storage day ", x)
             )) +
  scale_fill_gradient2(
    low = "white",
    high = "#b2182b",
    midpoint = 0,
    name = "PRCC"
  ) +
  labs(
    title = "PRCC Heatmap Across Parameters at Different Consumer Storage Days",
    x = "Parameter",
    y = "Isolate"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    panel.spacing = unit(1, "lines")
  )

# Figure 4 
library(patchwork)
dat  = read.csv("Raw_data_validation.csv")
Group <- c(rep("II", 12),rep("IV", 4),rep("V", 8),rep("VII", 8),rep("II", 8),rep("III", 12),rep("I", 12),rep("IV", 12),rep("V", 4),rep("II", 8))
dat$Group = Group

shape_palette <- c(16, 17, 15, 3, 7, 8, 0, 1, 2, 5, 6, 9)
make_group_plot <- function(df, group_name) {
  df_g <- df %>%
    filter(Group == group_name)
  df_g$Isolate <- factor(df_g$Isolate)
  iso_levels <- levels(df_g$Isolate)
  shape_map <- setNames(
    shape_palette[seq_along(iso_levels)],
    iso_levels
  )
  ggplot(df_g, aes(
    x = Consumer.storage.day,
    y = Log.Count,
    color = Closest.Type.Strain,
    shape = Isolate,
    group = Isolate
  )) +
    geom_point(size = 3) +
    geom_line() +
    scale_shape_manual(values = shape_map) +
    coord_cartesian(ylim = c(0, 6.5)) +
    labs(
      title = paste("Group", group_name),
      x = "Consumer Storage Day",
      y = "B. cereus concentration (log CFU/mL)",
      color = "Species",
      shape = "Isolate"
    ) +
    guides(
      color = guide_legend(order = 1),
      shape = guide_legend(order = 2)
    ) +
    theme_bw() +
    theme(
      text = element_text(size = 12),
      legend.position = "right"
    )
}
group_levels <- c("I", "II", "III", "IV", "V", "VII")
plots <- lapply(group_levels, function(g) {
  make_group_plot(dat, g)
})
names(plots) <- group_levels
wrap_plots(plots)

net_grow <- dat %>%
  filter(Consumer.storage.day %in% c(0, 35)) %>%
  pivot_wider(names_from = Consumer.storage.day, values_from = Log.Count) %>%
  mutate(LogDiff = `35` - `0`)

net_grow %>%
  filter(Closest.Type.Strain != "pseudomycoides") %>%
  summarise(
    min_LogDiff = min(LogDiff, na.rm = TRUE),
    max_LogDiff = max(LogDiff, na.rm = TRUE)
  )

 net_grow %>%
  group_by(Closest.Type.Strain) %>%
  summarise(
    avg_LogDiff = mean(LogDiff, na.rm = TRUE),
    n = n()
  )

  net_grow %>%
   group_by(Closest.Type.Strain) %>%
   summarise(
     min_LogDiff = min(LogDiff, na.rm = TRUE),
     max_LogDiff = max(LogDiff, na.rm = TRUE),
     diff_LogDiff = max_LogDiff - min_LogDiff,
     .groups = "drop"
   )
 
 
# Figure 5
species_map <- tibble::tibble(
  isolate_id = c(
    "PS00134","PS00136","PS00123",
    "PS00161","PS00423","PS00048",
    "PS00906","PS00515",
    "PS00047","PS00452",
    "PS00436","PS00485","PS00513",
    "PS00087",
    "PS00107","PS00512","PS00597",
    "PS00598","PS00565","PS00809",
    "PS00539","PS00399"
  ),
  species = c(
    rep("B. pseudomycoides", 3),
    rep("B. albus", 3),
    rep("B. mobilis", 2),
    rep("B. tropicus", 2),
    rep("B. pacificus", 3),
    rep("B. cereus", 1),
    rep("B. thuringiensis", 3),
    rep("B. toyonensis", 3),
    rep("B. cytotoxicus", 2)
  )
)

# Day 14
exp_count_d14 = subset(N0_df, Consumer.storage.day == "14")
Boxplot_d14 = final_results %>%
  filter(t_H == 14) %>%
  left_join(species_map, by = "isolate_id") %>%
  mutate(isolate_id = factor(isolate_id, levels = iso_order)) %>%
  ggplot(aes(x = isolate_id, y = conc, fill = species)) +
  geom_boxplot() +
  geom_point(
    data = exp_count_d14,
    aes(x = Isolate, y = Log.Count, color = "Observed"),
    inherit.aes = FALSE,
    size = 2
  ) +
  scale_color_manual(
    name = "Data",
    values = c("Observed" = "red")
  ) +
  labs(
    title = "Predicted vs Observed B. cereus group concentration (Day 14)",
    x = "Isolate",
    y = "Concentration (log CFU/mL)",
    fill = "Species"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
 
# Day 21
exp_count_d21 = subset(N0_df, Consumer.storage.day == "21")
Boxplot_d21 = final_results %>%
  filter(t_H == 21) %>%
  left_join(species_map, by = "isolate_id") %>%
  mutate(isolate_id = factor(isolate_id, levels = iso_order)) %>%
  ggplot(aes(x = isolate_id, y = conc, fill = species)) +
  geom_boxplot() +
  geom_point(
    data = exp_count_d21,
    aes(x = Isolate, y = Log.Count, color = "Observed"),
    inherit.aes = FALSE,
    size = 2
  ) +
  scale_color_manual(
    name = "Data",
    values = c("Observed" = "red")
  ) +
  labs(
    title = "Predicted vs Observed B. cereus group concentration (Day 21)",
    x = "Isolate",
    y = "Concentration (log CFU/mL)",
    fill = "Species"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Day 35
exp_count_d35 = subset(N0_df, Consumer.storage.day == "35")
Boxplot_d35 = final_results %>%
  filter(t_H == 35) %>%
  left_join(species_map, by = "isolate_id") %>%
  mutate(isolate_id = factor(isolate_id, levels = iso_order)) %>%
  ggplot(aes(x = isolate_id, y = conc, fill = species)) +
  geom_boxplot() +
  geom_point(
    data = exp_count_d35,
    aes(x = Isolate, y = Log.Count, color = "Observed"),
    inherit.aes = FALSE,
    size = 2
  ) +
  scale_color_manual(
    name = "Data",
    values = c("Observed" = "red")
  ) +
  labs(
    title = "Predicted vs Observed B. cereus group concentration (Day 35)",
    x = "Isolate",
    y = "Concentration (log CFU/mL)",
    fill = "Species"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Calculate IQR
iqr_by_isolate_time <- final_results %>%
  group_by(isolate_id, t_H) %>%
  summarise(
    Q1 = quantile(conc, 0.25, na.rm = TRUE),
    Q3 = quantile(conc, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    .groups = "drop"
  ) %>%
  arrange(isolate_id, t_H)

N0_lookup <- N0_df %>%
  distinct(Isolate, .keep_all = TRUE) %>%
  select(Isolate, Closest.Type.Strain)

iqr_with_strain <- iqr_by_isolate_time %>%
  left_join(N0_lookup, by = c("isolate_id" = "Isolate"))

mean_iqr <- iqr_with_strain %>%
  group_by(t_H) %>%
  summarise(mean_IQR = mean(IQR, na.rm = TRUE))

mean_iqr_df <- iqr_with_strain %>%
  group_by(Closest.Type.Strain, t_H) %>%
  summarise(
    mean_IQR = mean(IQR, na.rm = TRUE),
    .groups = "drop"
  )

mean_iqr_df %>%
  filter(t_H == 14) %>%
  summarise(
    min_mean_IQR = min(mean_IQR, na.rm = TRUE),
    min_strain = Closest.Type.Strain[which.min(mean_IQR)],
    max_mean_IQR = max(mean_IQR, na.rm = TRUE),
    max_strain = Closest.Type.Strain[which.max(mean_IQR)]
  )

mean_iqr_df %>%
  filter(t_H == 21) %>%
  summarise(
    min_mean_IQR = min(mean_IQR, na.rm = TRUE),
    min_strain = Closest.Type.Strain[which.min(mean_IQR)],
    max_mean_IQR = max(mean_IQR, na.rm = TRUE),
    max_strain = Closest.Type.Strain[which.max(mean_IQR)]
  )

iqr_diff_by_strain_time <- iqr_with_strain %>%
  group_by(Closest.Type.Strain, t_H) %>%
  summarise(
    IQR_difference = max(IQR, na.rm = TRUE) - min(IQR, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Closest.Type.Strain, t_H)

# Get percentile for observed data
percentile_results <- N0_df %>%
  rowwise() %>%
  mutate(
    percentile = mean(
      final_results$conc[
        final_results$isolate_id == Isolate &
          final_results$t_H == Consumer.storage.day
      ] <= Log.Count,
      na.rm = TRUE
    ) * 100
  ) %>%
  ungroup()

quartile_results <- N0_df %>%
  rowwise() %>%
  mutate(
    q1 = quantile(
      final_results$conc[
        final_results$isolate_id == Isolate &
          final_results$t_H == Consumer.storage.day
      ],
      0.25,
      na.rm = TRUE
    ),
    
    q3 = quantile(
      final_results$conc[
        final_results$isolate_id == Isolate &
          final_results$t_H == Consumer.storage.day
      ],
      0.75,
      na.rm = TRUE
    ),
    
    category = case_when(
      Log.Count < (q1 - 1) ~ "Below Q1 by >1 log",
      Log.Count < q1 ~ "Below Q1",
      Log.Count > q3 ~ "Above Q3",
      TRUE ~ "Between Q1 and Q3"
    )
  ) %>%
  ungroup()

