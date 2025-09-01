# Format input data for VMP
vmp_datalist <- function(df, N, P){
  time_obs = map(1:N, function(i){
    subj_time = map(1:P, function(p){
      vals = df %>%
        filter(Subj == i & Var == p) %>%
        arrange(Arg) %>%
        pull(Arg)
      return(vals)
    })
    return(subj_time)
  })
  
  Y_vals = map(1:N, function(i){
    subj_vals = map(1:P, function(p){
      vals = df %>%
        filter(Subj == i & Var == p) %>%
        arrange(Arg) %>%
        pull(Y)
      return(vals)
    })
    return(subj_vals)
  })
  
  return(list(t = time_obs, y = Y_vals))
}

# Format VMP outputs for evaluation
extract_vmp <- function(mod, domain, N, K, P, anchor){
  cond_phi_list = map(1:K, function(k){
    return(c(mod$list_Psi_hat[[k]]))
  })
  cond_phi = do.call(cbind, cond_phi_list)
  
  # Align with Procrustes
  align = procrust_FPC(list(cond_phi), anchor)
  cond_phi = align[[1]]
  
  # Extract other elements
  Mu_vec = c(mod$mu_hat)
  Mu_df = data.frame(Arg = rep(domain, P),
                     Var = rep(1:P, each = length(domain)),
                     Est = Mu_vec, LB = NA, UB = NA)
  EF_df = FPC_df(cond_phi, sim_domain, P) %>%
    mutate(LB = NA, UB = NA)
  Smooth_df = map(1:N, function(x){
    person_df = map(1:P, function(p){
      output = data.frame(Est = mod$Y_hat[[x]][[p]], 
                          LB = mod$Y_low[[x]][[p]], 
                          UB = mod$Y_upp[[x]][[p]], 
                          Arg = sim_domain, 
                          Curve = paste0("Curve ", x), 
                          Var = p)
      return(output)
    }) %>% list_rbind()
    return(person_df)
  }) %>% list_rbind()
  
  return(list(FE = Mu_df, EF = EF_df, Smooths = Smooth_df))
}
