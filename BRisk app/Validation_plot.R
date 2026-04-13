library(ggplot2)
library(dplyr)
library(gridExtra)

df = read.csv("Validation dataset_1.csv")
df2 = read.csv("Validation dataset_2.csv")
df3 = read.csv("Validation dataset_3.csv")

# Plots (species with 1 isolate in the database)
p_pseudomycoides <- df %>%
  filter(Closest.Type.Strain == "pseudomycoides") %>%
  ggplot(aes(x = Consumer.storage.day,
             colour = Isolate)) +   # color by isolate
  geom_point(aes(y = Observed.Count, shape = "Observed"), size = 3) +
  geom_point(aes(y = Predicted.Count, shape = "Predicted"), 
             size = 3, colour = "black") +
  scale_shape_manual(values = c("Observed" = 16, "Predicted" = 17)) +
  scale_y_continuous(limits = c(0, 7)) +   
  labs(
    title = expression(italic("pseudomycoides")),
    x = "Consumer storage day",
    y = "Bacterial count (log CFU/mL)",
    colour = "Isolate",
    shape = ""
  ) +
  guides(
    shape = guide_legend(order = 1),
    colour = guide_legend(order = 3),   # isolate colors moved to bottom
    `Predicted points` = guide_legend(order = 2) # predicted colors
  ) +
  theme_bw(base_size = 14)
p_pseudomycoides

p_albus <- df %>%
  filter(Closest.Type.Strain == "albus") %>%
  ggplot(aes(x = Consumer.storage.day,
             colour = Isolate)) +   # color by isolate
  geom_point(aes(y = Observed.Count, shape = "Observed"), size = 3) +
  geom_point(aes(y = Predicted.Count, shape = "Predicted"), 
             size = 3, colour = "black") +
  scale_shape_manual(values = c("Observed" = 16, "Predicted" = 17)) +
  scale_y_continuous(limits = c(0, 7)) +   
  labs(
    title = expression(italic("albus")),
    x = "Consumer storage day",
    y = "Bacterial count (log CFU/mL)",
    colour = "Isolate",
    shape = ""
  ) +
  guides(
    shape = guide_legend(order = 1),
    colour = guide_legend(order = 3),   # isolate colors moved to bottom
    `Predicted points` = guide_legend(order = 2) # predicted colors
    ) +
  theme_bw(base_size = 14)
p_albus

p_mobilis <- df %>%
  filter(Closest.Type.Strain == "mobilis") %>%
  ggplot(aes(x = Consumer.storage.day,
             colour = Isolate)) +   # color by isolate
  geom_point(aes(y = Observed.Count, shape = "Observed"), size = 3) +
  geom_point(aes(y = Predicted.Count, shape = "Predicted"), 
             size = 3, colour = "black") +
  scale_shape_manual(values = c("Observed" = 16, "Predicted" = 17)) +
  scale_y_continuous(limits = c(0, 7)) +   
  labs(
    title = expression(italic("mobilis")),
    x = "Consumer storage day",
    y = "Bacterial count (log CFU/mL)",
    colour = "Isolate",
    shape = ""
  ) +
  guides(
    shape = guide_legend(order = 1),
    colour = guide_legend(order = 3),   # isolate colors moved to bottom
    `Predicted points` = guide_legend(order = 2) # predicted colors
  ) +
  theme_bw(base_size = 14)
p_mobilis

p_pacificus <- df %>%
  filter(Closest.Type.Strain == "pacificus") %>%
  ggplot(aes(x = Consumer.storage.day)) +
  # observed = colored by isolate
  geom_point(aes(y = Observed.Count, colour = Isolate, shape = "Observed"), size = 3) +
  # predicted = black triangles
  geom_point(aes(y = Predicted.Count, shape = "Predicted"),
             size = 3, colour = "black") +
  scale_shape_manual(values = c("Observed" = 16, "Predicted" = 17)) +
  scale_y_continuous(limits = c(0, 7)) +   
  labs(
    title = expression(italic("pacificus")),
    x = "Consumer storage day",
    y = "Bacterial Count (log CFU/mL)",
    colour = "Isolate",
    shape = ""
  ) +
  guides(
    shape = guide_legend(order = 1),
    colour = guide_legend(order = 2)
  ) +
  theme_bw(base_size = 14)
p_pacificus

# Plots (species with 2 isolates in the database)
p_tropicus <- df2 %>%
  filter(Closest.Type.Strain == "tropicus") %>%
  ggplot(aes(x = Consumer.storage.day)) +
  
  # Observed: colored by isolate
  geom_point(aes(y = Observed.Count,
                 colour = Isolate,
                 shape = "Observed"),
             size = 3) +
  
  # Predicted Q1: black triangle
  geom_point(aes(y = Pred_Q1,
                 shape = "Predicted (1st quartile)"),
             size = 3,
             colour = "lightgrey") +
  
  # Shape mapping (all predicted values share same triangle)
  scale_shape_manual(values = c(
    "Observed" = 16,         # solid circle
    "Predicted (1st quartile)" = 17     # triangle
  )) +
  
  scale_y_continuous(limits = c(0, 7.5)) +
  
  labs(
    title = expression(italic("tropicus")),
    x = "Consumer storage day",
    y = "Bacterial count (log CFU/mL)",
    colour = "Isolate",
    shape = ""
  ) +
  guides(
    shape = guide_legend(order = 1),
    colour = guide_legend(order = 3),   # isolate colors moved to bottom
    `Predicted points` = guide_legend(order = 2) # predicted colors
  ) +
  theme_bw(base_size = 14)

p_tropicus

p_thuringiensis <- df2 %>%
  filter(Closest.Type.Strain == "thuringiensis") %>%
  ggplot(aes(x = Consumer.storage.day)) +
  
  # Observed: colored by isolate
  geom_point(aes(y = Observed.Count,
                 colour = Isolate,
                 shape = "Observed"),
             size = 3) +
  
  # Predicted Q1: black triangle
  geom_point(aes(y = Pred_Q1,
                 shape = "Predicted (1st quartile)"),
             size = 3,
             colour = "lightgrey") +
  
  # Shape mapping (all predicted values share same triangle)
  scale_shape_manual(values = c(
    "Observed" = 16,         # solid circle
    "Predicted (1st quartile)" = 17     # triangle
  )) +
  
  scale_y_continuous(limits = c(0, 7.5)) +
  
  labs(
    title = expression(italic("thuringiensis")),
    x = "Consumer storage day",
    y = "Bacterial count (log CFU/mL)",
    colour = "Isolate",
    shape = ""
  ) +
  
  guides(
    shape = guide_legend(order = 1),
    colour = guide_legend(order = 3),   # isolate colors moved to bottom
    `Predicted points` = guide_legend(order = 2) # predicted colors
  ) +
  
  theme_bw(base_size = 14)

p_thuringiensis

p_toyonensis <- df2 %>%
  filter(Closest.Type.Strain == "toyonensis") %>%
  ggplot(aes(x = Consumer.storage.day)) +
  
  # Observed: colored by isolate
  geom_point(aes(y = Observed.Count,
                 colour = Isolate,
                 shape = "Observed"),
             size = 3) +
  
  # Predicted Q1: black triangle
  geom_point(aes(y = Pred_Q1,
                 shape = "Predicted (1st quartile)"),
             size = 3,
             colour = "lightgrey") +
  
  # Shape mapping (all predicted values share same triangle)
  scale_shape_manual(values = c(
    "Observed" = 16,         # solid circle
    "Predicted (1st quartile)" = 17      # triangle
  )) +
  
  scale_y_continuous(limits = c(0, 7.5)) +
  
  labs(
    title = expression(italic("toyonensis")),
    x = "Consumer storage day",
    y = "Bacterial count (log CFU/mL)",
    colour = "Isolate",
    shape = ""
  ) +
  
  guides(
    shape = guide_legend(order = 1),
    colour = guide_legend(order = 3),   # isolate colors moved to bottom
    `Predicted points` = guide_legend(order = 2) # predicted colors
  ) +
  
  theme_bw(base_size = 14)

p_toyonensis

p_cytotoxicus <- df2 %>%
  filter(Closest.Type.Strain == "cytotoxicus") %>%
  ggplot(aes(x = Consumer.storage.day)) +
  
  # Observed: colored by isolate
  geom_point(aes(y = Observed.Count,
                 colour = Isolate,
                 shape = "Observed"),
             size = 3) +
  
  # Predicted Q1: black triangle
  geom_point(aes(y = Pred_Q1,
                 shape = "Predicted (1st quartile)"),
             size = 3,
             colour = "lightgrey") +
  
  # Shape mapping (all predicted values share same triangle)
  scale_shape_manual(values = c(
    "Observed" = 16,         # solid circle
    "Predicted (1st quartile)" = 17      # triangle
  )) +
  
  scale_y_continuous(limits = c(0, 7.5)) +
  
  labs(
    title = expression(italic("cytotoxicus")),
    x = "Consumer storage day",
    y = "Bacterial count (log CFU/mL)",
    colour = "Isolate",
    shape = ""
  ) +
  
  guides(
    shape = guide_legend(order = 1),
    colour = guide_legend(order = 3),   # isolate colors moved to bottom
    `Predicted points` = guide_legend(order = 2) # predicted colors
  ) +
  
  theme_bw(base_size = 14)

p_cytotoxicus

p_cereus <- df3 %>%
  filter(Closest.Type.Strain == "cereus") %>%
  ggplot(aes(x = Consumer.storage.day)) +
  
  # Observed: colored by isolate
  geom_point(aes(y = Observed.Count,
                 colour = Isolate,
                 shape = "Observed"),
             size = 3) +
  
  # Predicted Q1: black triangle
  geom_point(aes(y = Pred_Q1,
                 shape = "Predicted (1st quartile)"),
             size = 3,
             colour = "lightgrey") +
  
  # Shape mapping (all predicted values share same triangle)
  scale_shape_manual(values = c(
    "Observed" = 16,         # solid circle
    "Predicted (1st quartile)" = 17    # triangle
  )) +
  
  scale_y_continuous(limits = c(0, 7)) +
  
  labs(
    title = expression(italic("cereus")),
    x = "Consumer storage day",
    y = "Bacterial count (log CFU/mL)",
    colour = "Isolate",
    shape = ""
  ) +
  
  guides(
    shape = guide_legend(order = 1),
    colour = guide_legend(order = 3),   # isolate colors moved to bottom
    `Predicted points` = guide_legend(order = 2) # predicted colors
  ) +
  
  theme_bw(base_size = 14)

p_cereus

library(gridExtra)
grid.arrange(p_pseudomycoides,p_albus,p_mobilis,p_pacificus, ncol = 2)
grid.arrange(p_tropicus,p_thuringiensis,p_toyonensis,p_cytotoxicus, ncol = 2)
