library(lme4)
library(lmerTest)

df = read.csv("Validation_combined.csv")

# factors
df$Isolate <- factor(df$Isolate)
df$Closest.Type.Strain <- factor(df$Closest.Type.Strain)
df$Consumer.storage.day <- factor(df$Consumer.storage.day)

# mixed model
model <- lmer(Diff ~ Closest.Type.Strain + Consumer.storage.day + (1 | Isolate), data = df)
anova(model)

# assumption check
# normality of residuals
shapiro.test(residuals(model))

# homoscedasticity
lm_tmp <- lm(residuals(model) ~ fitted(model))
lmtest::bptest(lm_tmp)
