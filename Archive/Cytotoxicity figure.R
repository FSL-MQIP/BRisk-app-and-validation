library(ggplot2)
library(dplyr)

cytotoxicity_input = read.csv("Cytotoxicity_data.csv")
colnames(cytotoxicity_input)[1] <- "Isolate.Name"
colnames(cytotoxicity_input)[colnames(cytotoxicity_input) == "Average_Cell_Viability_F"] <- "Normalized_Cytotoxicity"
cytotoxicity_input$panC_Group <- trimws(cytotoxicity_input$panC_Group)

# Calculate quartiles and min/max
quartiles <- quantile(cytotoxicity_input$Normalized_Cytotoxicity, probs = c(0, 0.25, 0.75, 1))
quartile_labels <- c("Min", "Q1", "Q3", "Max")
label_text <- paste0(quartile_labels, ": ", round(quartiles, 2))

# Compute y_max dynamically for this group
hist_data <- ggplot_build(
  ggplot(subset(cytotoxicity_input, panC_Group == "Group_IV"), 
         aes(x = Normalized_Cytotoxicity)) + geom_histogram(binwidth = 0.05)
)$data[[1]]
y_max <- max(hist_data$count)

# Position labels slightly above the tallest bar (10% above)
label_y <- y_max * 1.10  

ggplot() +
  # Histogram of all isolates
  geom_histogram(data = cytotoxicity_input, 
                 aes(x = Normalized_Cytotoxicity, fill = "All Isolates"), 
                 binwidth = 0.05, alpha = 0.8) +
  # Histogram of Group X
  geom_histogram(data = subset(cytotoxicity_input, panC_Group == "Group_IV"), 
                 aes(x = Normalized_Cytotoxicity, fill = "Phylogenetic Group"), 
                 binwidth = 0.05) +
  # Vertical lines for min, Q1, Q3, max
  geom_vline(xintercept = quartiles, linetype = "dashed", color = "red", size = 1) +
  # Labels for the quartiles with real numbers, moved up and smaller
  geom_text(aes(x = quartiles, y = label_y, label = label_text),
            color = "red", angle = 90, vjust = -0.5, hjust = 0, size = 3) +
  xlab("Cytotoxicity Value") +
  ylab("Number of Isolates") +
  ggtitle("Phylogenetic Group IV") + 
  theme_minimal() +
  theme(plot.title = element_text(size = 24, face = "bold"),       
        axis.title.x = element_text(size = 22),                    
        axis.title.y = element_text(size = 22),  
        axis.text.x = element_text(size = 22),                    
        axis.text.y = element_text(size = 22),
        legend.text = element_text(size = 22)) +
  scale_fill_manual(values = c("All Isolates" = "lightblue", 
                               "Phylogenetic Group" = "yellow"),
                    labels = c("All Isolates", "Phylogenetic Group IV (n=82)")) +
  labs(fill = "")
