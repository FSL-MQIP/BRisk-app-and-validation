library(ggplot2)
library(tidyverse)

dat = read.csv("Exposure data input.csv")

# Convert to long format
df_long <- dat %>%
  pivot_longer(cols = c(log_3_5, log_5plus),
               names_to = "Category",
               values_to = "Value")

df_long$Species <- factor(df_long$Species, levels = unique(dat$Species))

# --- Plot: 3–5 log ---
p1 <- df_long %>%
  filter(Category == "log_3_5") %>%
  ggplot(aes(x = Species, y = Value, fill = factor(Day))) +
  geom_col(position = position_dodge(width = 0.8)) +
  scale_y_continuous(limits = c(0, 5)) +  
  labs(title = "B",
       x = "Species",
       y = "Percent milk units 3 to 5 logs",
       fill = "Consumer storage day") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())

# --- Plot: >5 log ---
p2 <- df_long %>%
  filter(Category == "log_5plus") %>%
  ggplot(aes(x = Species, y = Value, fill = factor(Day))) +
  geom_col(position = position_dodge(width = 0.8)) +
  scale_y_continuous(limits = c(0, 5)) + 
  labs(title = "A",
       x = "Species",
       y = "Percent milk units over 5 logs",
       fill = "Consumer storage day") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())

# Show plots
p1
p2
