# Format input data for PACE
PACE_smooths <- function(mu, scores, score_var, FPC, N_S, Domain){
  mu_df = data.frame(Mu = mu, Arg = Domain)
  
  
  
  return(smooth_samples)
}

# Format PACE outputs for evaluation
PACE_datalist <- function(Y, Subj, S, N, sub_domain){
  lY = map(1:N, function(x){
    return(Y[Subj == x])
  })
  lT = map(1:N, function(x){
    return(sub_domain[S[Subj == x]])
  })
  return(list(list_y = lY, list_t = lT))
}

PACE_extract <- function(mod, K, anchor, sim_domain, N_S = 500){
  # Align with Procrustes
  pace_phi = procrust_FPC(list(mod$phi), anchor)[[1]]
  
  Mu_df = data.frame(Arg = sim_domain, Est = mod$mu, LB = NA, UB = NA)
  EF_df = FPC_df(pace_phi, sim_domain) %>%
    mutate(LB = NA, UB = NA)
  Smooth_df = map(1:nrow(mod$xiEst), function(x){
    score_samples = mvrnorm(N_S, mu = mod$xiEst[x,], Sigma = mod$xiVar[[x]])
    sample_matrix = mod$phi %*% t(score_samples)
    sample_df = data.frame(sample_matrix)
    colnames(sample_df) = 1:N_S
    sample_df$Arg = sim_domain
    sample_df = sample_df %>%
      pivot_longer(-c(Arg), names_to = "Sample", values_to = "Smooth") %>%
      left_join(Mu_df, by = "Arg") %>%
      mutate(Smooth = Smooth + Est) %>%
      select(-c(Est, LB, UB)) %>%
      mutate(Curve = paste0("Curve ", x))
    return(sample_df)
  }) %>% list_rbind() %>%
    group_by(Curve, Arg) %>%
    summarize(Est = mean(Smooth), 
              LB = quantile(Smooth, probs = c(0.025)), 
              UB = quantile(Smooth, probs = c(0.975)))
  
  return(list(Mu = Mu_df, EF = EF_df, Smooth = Smooth_df))
}
