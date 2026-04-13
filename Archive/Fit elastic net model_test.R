set.seed(2025)

n <- 20       # number of samples
num_snps <- 5 # number of SNPs
snp_levels <- 4

# ----- SNPs -----
# Each SNP is coded by 4 dummy variables (one-hot encoding)
SNP_list <- lapply(1:num_snps, function(s) {
  mat <- t(sapply(1:n, function(i) {
    x <- rep(0, snp_levels)
    x[sample(1:snp_levels, 1)] <- 1
    x
  }))
  colnames(mat) <- paste0("SNP", s, "_code", 1:snp_levels)
  mat
})

# Combine all SNPs
SNP_mat <- do.call(cbind, SNP_list)

# ----- Phylogenetic group -----
phylo_group <- factor(sample(paste0("G", 1:8), n, replace = TRUE))

# ----- Cytotoxicity outcome -----
# Simulate as a function of a couple SNP codes + phylo_group effects
cytotoxicity <- 1.5*SNP_mat[, "SNP2_code2"] - 
  2*SNP_mat[, "SNP4_code4"] +
  ifelse(phylo_group %in% c("G3","G6"), 2, 0) +
  rnorm(n, mean=0, sd=1)

# ----- Final dataframe -----
df <- data.frame(cytotoxicity, phylo_group, SNP_mat)

head(df)

library(glmnet)

# Create design matrix (dummy-encodes phylo_group automatically)
X <- model.matrix(cytotoxicity ~ . , data = df[, -1])[,-1]
y <- df$cytotoxicity

# Cross-validated elastic net (alpha=0.5 = halfway between ridge & lasso)
set.seed(123)
cvfit <- cv.glmnet(X, y, alpha = 0.5, family = "gaussian")

# Best lambda
best_lambda <- cvfit$lambda.min

# Final model at best lambda
final_model <- glmnet(X, y, alpha = 0.5, lambda = best_lambda, family = "gaussian")

# Coefficients
coef(final_model)

# New strain data
new_strain <- data.frame(
  phylo_groupG2 = 1, phylo_groupG3 = 0, phylo_groupG4 = 0, 
  phylo_groupG5 = 0, phylo_groupG6 = 0, phylo_groupG7 = 0, 
  phylo_groupG8 = 0, 
  SNP1_code1 = 0, SNP1_code2 = 0, SNP1_code3 = 0, SNP1_code4 = 0,
  SNP2_code1 = 1, SNP2_code2 = 0, SNP2_code3 = 0, SNP2_code4 = 0,
  SNP3_code1 = 0, SNP3_code2 = 0, SNP3_code3 = 1, SNP3_code4 = 0,
  SNP4_code1 = 0, SNP4_code2 = 1, SNP4_code3 = 0, SNP4_code4 = 0,
  SNP5_code1 = 0, SNP5_code2 = 0, SNP5_code3 = 0, SNP5_code4 = 1
)

# Make design matrix
X_new <- model.matrix(~ ., data = new_strain)

# Predict cytotoxicity
pred_cyto <- predict(final_model, newx = X_new[,-1])
pred_cyto
