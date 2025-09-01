# Format input data for GFSR
GFSR_datalist <- function(Y, Subj, Arg, N, K, Q, Domain, alpha = 0.1){
  
  # FE matrix - just a population intercept
  X_mat = matrix(0, nrow = N, ncol = 1)
  X_mat[,1] = 1
  
  # Splines
  oB = bSpline(Domain, df = Q, intercept = T)
  BS = bSpline(Arg, df = Q, intercept = T)
  
  # Penalty matrix - taken from Goldsmith et al. in 2015
  D = length(Domain)
  diff0 = diag(1, D, D)
  diff2 = matrix(rep(c(1,-2,1, rep(0, D-2)), D-2)[1:((D-2)*D)], D-2, D, byrow = TRUE)
  P0 = t(oB) %*% t(diff0) %*% diff0 %*% oB
  P2 = t(oB) %*% t(diff2) %*% diff2 %*% oB
  P.mat = alpha * P0 + (1-alpha) * P2
  
  # Format for GFSR
  stan_dat = list(N = length(Y), I = N, D = length(Domain), p = 1,
                  Kt = Q, Kp = K, Y = Y, subjId = Subj,
                  X = X_mat, BS = BS, PenMat = P.mat, B = oB)
  return(stan_dat)
}

# Format GFSR outputs for evaluation
GFSR_extract <- function(mod_fit, BS, M, K){
  samples = extract(mod_fit)
  n_samp = length(samples$sigma2)
  
  smooth_list = map(1:n_samp, function(x){
    smooth_sample = BS %*% t(samples$beta_psi[x,,])
    return(svd(t(smooth_sample)))
  })
  
  EF_list = map(smooth_list, function(svd_x){
    return(svd_x$v[,1:K]*sqrt(M))
  })
  
  Score_list = map(1:n_samp, function(x){
    U = smooth_list[[x]]$u[,1:K]
    D = diag(smooth_list[[x]]$d[1:K])
    Score_sample = samples$c[x,,] %*% U %*% D / sqrt(M)
    return(Score_sample)
  })
  
  Mu_list = map(1:n_samp, function(x){
    return(BS %*% samples$beta[x,,])
  }) 
  
  return(list(EF = EF_list, Score = Score_list, Mu = Mu_list))
}
