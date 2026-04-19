library(car)

cytotoxicity_input = read.csv("Cytotoxicity_data.csv")
colnames(cytotoxicity_input)[1] <- "Isolate.Name"
colnames(cytotoxicity_input)[colnames(cytotoxicity_input) == "Average_Cell_Viability_F"] <- "Normalized_Cytotoxicity"
cytotoxicity_input$panC_Group <- trimws(cytotoxicity_input$panC_Group)
cytotoxicity_input <- subset(cytotoxicity_input, !(panC_Group %in% c("Group_clarus","Group_VI")))

anova_model <- aov(Normalized_Cytotoxicity ~ panC_Group, data = cytotoxicity_input)

# test normality assumption 
shapiro.test(residuals (anova_model)) # p > 0.05
qqnorm(residuals(anova_model))
qqline(residuals(anova_model))

# test equal variance assumption
leveneTest(Normalized_Cytotoxicity ~ panC_Group, data = cytotoxicity_input) 

summary(anova_model)
stats = TukeyHSD(anova_model)