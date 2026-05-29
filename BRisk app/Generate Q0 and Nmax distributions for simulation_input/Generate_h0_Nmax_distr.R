library(dplyr)
source("UtilityFunctions_Q0.R")

simulation_input = read.csv("simulation_input_1.csv")

# Calculate Q0 
h0_Q0_table = read.csv("h0_table.csv")
h0_Q0_table$Q0 = sapply(h0_Q0_table$h0, Calculate_Q0)
h0_Q0_table$LOG10Q0 = log10(h0_Q0_table$Q0)
h0_Q0_table <- h0_Q0_table %>%
  group_by(isolate) %>%
  summarise(
    mean_LOG10Q0 = mean(LOG10Q0, na.rm = TRUE),
    sd_LOG10Q0 = sd(LOG10Q0, na.rm = TRUE)
  )


# Calculate Nmax 
Nmax_table = read.csv("Nmax_table.csv")
Nmax_table <- Nmax_table %>%
  group_by(isolate) %>%
  summarise(
    mean_Nmax = mean(LOG10Nmax, na.rm = TRUE),
    sd_Nmax   = sd(LOG10Nmax, na.rm = TRUE)
  )

# Simulation input 
simulation_input <- simulation_input %>%
  left_join(h0_Q0_table, by = "isolate") %>%
  left_join(Nmax_table, by = "isolate")

write.csv(simulation_input,"simulation_input.csv")
