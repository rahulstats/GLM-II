edge_pairs <- t(combn(K, 2))
nu_true = nu_fit; nu_transitive = hat_nu_linear; sigma = as.numeric(sqrt(sigma2_hat))

# Calculate the true win probabilities under the cyclic model
tau <- pnorm(nu_true / sigma)

# Calculate the expected win probabilities under the baseline transitive model
omega <- pnorm(nu_transitive / sigma)

# Initialize the expected gain vector
Win = 0

# Loop through all K2 edges to calculate the expected gain for each match-up
for (ai in 1:K2) {
  # Case 1: The cyclic model predicts a higher win probability than the baseline.
  # Expected gain from betting 1 unit FOR the outcome.
  if(tau[ai] > omega[ai]) {
    Win[ai] = (tau[ai] - omega[ai]) / omega[ai]
  } 
  # Case 2: The cyclic model predicts a lower win probability than the baseline.
  # Expected gain from betting 1 unit AGAINST the outcome (for the complement).
  else if(tau[ai] < omega[ai]) {
    Win[ai] = (omega[ai] - tau[ai]) / (1 - omega[ai])
  } 
  # Case 3: Both models predict the exact same probability. No edge exists.
  else {
    Win[ai] = 0
  }
}

# Combine the edge indices, probabilities, and calculated expected gains into a single matrix
win_hat = cbind(edge_pairs, tau, omega, Win)

# Calculate the total expected gain over the entire round-robin tournament
sum(Win)

# Output the expected gain calculation specifically for items 2,5,19
win_hat[win_hat[,5] !=0,]
