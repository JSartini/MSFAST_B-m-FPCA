data {
  int N;      // Number of time series
  int L;      // Total number of observed data points
  int K;      // Number of Eigenfunctions
  
  vector[L] Y;                    // Observed data
  array[L] int<upper=N> ts_idx;   // Time series indices
  matrix[L, K] FPC;               // FPC basis evaluated at observed points
}

parameters {
  real<lower=0> sigma2; // Error in observation
  
  // Components/weights
  positive_ordered[K] lambda;        // Eigenvalues
  matrix[N, K] Scores;               // EF scores (xi)
}

model {

  // Error
  sigma2 ~ inv_gamma(0.001, 0.001);
  
  // Score priors 
  for(i in 1:K){
    to_vector(Scores[,i]) ~ normal(0, sqrt(lambda[i]));
  }
  
  // Likelihood
  Y ~ normal(rows_dot_product(Scores[ts_idx, ], FPC), sqrt(sigma2)); 
}
