library(lme4)
library(lmerTest)

df = read.csv("Validation_result_combined.csv")

# factors
df$isolate_id <- factor(df$isolate_id)
df$species <- factor(df$species)
df$consumer.storage.day <- factor(df$consumer.storage.day)

# mixed model
model <- lmer(Diff ~ species + consumer.storage.day + (1 | isolate_id), data = df)
anova(model)

# assumption check
library(DHARMa)
# residual distribution
sim_res <- simulateResiduals(model)
plot(sim_res)
