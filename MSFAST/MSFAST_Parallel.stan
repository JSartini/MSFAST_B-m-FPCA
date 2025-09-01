functions {
  real partial_sum_lik(array[] int p_slice, int start, int end,
                      array[] int Tp_card, array[] int Subj, array[] int S,
                      matrix Psi, matrix Scores, matrix B, vector w_mu,
                      vector sigma2, array[] int start_indices, vector Y,
                      int Q, int M, int K) {
    real acc = 0;    // log density accumulator for this chunk

    for (i in 1:size(p_slice)) {
      int p = p_slice[i];
      int sdx = (p-1)*Q+1;
      int edx = p*Q;
      int Tp = Tp_card[p];

      array[Tp] int Subj_p = segment(Subj, start_indices[p], Tp);
      array[Tp] int S_p = segment(S, start_indices[p], Tp);
      matrix[M, K] Phi_mat = B * Psi[sdx:edx, ];
      vector[Tp] Theta = rows_dot_product(Scores[Subj_p, ], Phi_mat[S_p, ]);
      vector[M] mu = B * w_mu[sdx:edx];

      acc += normal_lpdf(segment(Y, start_indices[p], Tp) | mu[S_p] + Theta, sqrt(sigma2[p]));
    }
    return acc;
  }
}

data {
  int N;   // Number of time series/individuals
  int M;   // Cardinality of observation time set
  int L;   // Total number of observed data points
  int Q;   // Number of spline bases
  int K;   // Number of eigenfunctions
  int P;   // Number of covariates
  
  vector[L] Y;                         // Observed data for all variables stacked
  array[L] int<upper=N> Subj;          // Time series indices all variables stacked
  array[L] int<upper=M> S;             // Time point indices all variables stacked
  array[P] int Tp_card;                // Number of observed data points for each variable
  
  matrix[M, Q] B;               // Orthogonalized basis over T
  matrix[Q, Q] P_alpha;         // Penalty matrix for splines
}

transformed data{
  array[P] int start_indices;
  array[P] int p_vals;
  {
    int pos = 1;
    for(p in 1:P){
      p_vals[p] = p;
      start_indices[p] = pos;
      pos = pos + Tp_card[p];
    } 
  }
}

parameters {
  vector<lower=0>[P] sigma2; // Error in observation
  
  // Fixed-effect components
  vector[P*Q] w_mu;         // Population mean parameters
  vector<lower=0>[P] h_mu;  // Population mean smoothing parameter
  
  // Components/weights
  positive_ordered[K] lambda;        // Eigenvalues
  matrix<lower=0>[P,K] H;            // EF Smoothing parameters
  matrix[P*Q, K] X;                  // Unconstrained EF matrix
  matrix[N, K] Xi_Raw;               // EF scores unscaled
}

transformed parameters{
  // Orthogonal basis weights
  matrix[P*Q, K] Psi;
  matrix[N, K] Scores;
  
  // Polar decomposition
  {
    matrix[K,K] evec_XtX = eigenvectors_sym(crossprod(X)); 
    vector[K] eval_XtX = eigenvalues_sym(crossprod(X));
    Psi = X*evec_XtX*diag_matrix(1/sqrt(eval_XtX))*evec_XtX'; 
  }
  
  // Scaled scores
  Scores = Xi_Raw * diag_matrix(sqrt(lambda));
}

model {
  // Variance component priors
  lambda ~ inv_gamma(0.01, 0.01); 
  sigma2 ~ inv_gamma(0.01, 0.01);
  
  // Smoothing priors
  h_mu ~ gamma(0.01, 0.01); 
  to_vector(H) ~ gamma(0.01, 0.01); 
  
  int sdx;
  int edx;
  for(p in 1:P){
    sdx = (p-1)*Q+1;
    edx = p*Q;
    
    target += Q / 2.0 * log(h_mu[p]) - h_mu[p] / 2.0 * quad_form(P_alpha, w_mu[sdx:edx]);
    
    for(k in 1:K){
      target += Q / 2.0 * log(H[p,k]) -  H[p,k] / 2.0 * quad_form(P_alpha, Psi[sdx:edx,k]);
    }
  }
  
  // Uniform priors through matrix normals
  to_vector(X) ~ std_normal();
  
  // Score priors 
  to_vector(Xi_Raw) ~ std_normal();
  
  // Model likelihood
  target += reduce_sum(
    partial_sum_lik,      // likelihood over subset
    p_vals,               // split by covariates
    1,                    // grainsize of 1 (1 covariate each)
    Tp_card, Subj, S, Psi, Scores, B, w_mu, sigma2, start_indices, Y, Q, M, K
  );
}

