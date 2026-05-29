## Utility Functions 
# Last edited: 04/02/22

# Source: 
# Baranyi, J., & Roberts, T. A. (1994). 
# A dynamic approach to predicting bacterial growth in food. International journal of food microbiology, 23(3-4), 277-294.

# Purpose: Calculate Q0 

# Note: 
# mumax is in natural (i.e. exp(1)) scale to calculate h0 in this formula

# Function: 
Calculate_Q0 = function(h0) {
  Q0 <- 1/(exp(h0)-1)
  return(Q0)
}