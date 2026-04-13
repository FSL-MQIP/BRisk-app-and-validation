library(car)

dat = read.csv("Validation_combined.csv")
dat$Number.of.strains = as.factor(dat$Number.of.strains)

anova_model_1 <- aov(Diff ~ Closest.Type.Strain, data = dat)
qqnorm(residuals(anova_model_1))
qqline(residuals(anova_model_1))
leveneTest(Diff ~ Closest.Type.Strain, data = dat)
kruskal.test(Diff ~  Closest.Type.Strain, data = dat)

anova_model_2 <- aov(Diff ~ Number.of.strains, data = dat)
qqnorm(residuals(anova_model_2))
qqline(residuals(anova_model_2))
leveneTest(Diff ~ Number.of.strains, data = dat)
kruskal.test(Diff ~  Number.of.strains, data = dat)
