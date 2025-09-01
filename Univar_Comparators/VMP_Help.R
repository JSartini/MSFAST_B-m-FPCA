# Format input data for VMP
VMP_datalist <- function(Y, Subj, S, N, sub_domain){
  lY = map(1:N, function(x){
    return(Y[Subj == x])
  })
  lT = map(1:N, function(x){
    return(sub_domain[S[Subj == x]])
  })
  return(list(list_y = lY, list_t = lT))
}

# Format VMP outputs for evaluation
VMP_extract <- function(mod, K, anchor, sim_domain, N){
  # Align with Procrustes
  vmp_phi = procrust_FPC(list(mod$list_Psi_hat), anchor)[[1]]
  
  Mu_df = data.frame(Arg = sim_domain, Est = mod$mu_hat, LB = NA, UB = NA)
  EF_df = FPC_df(vmp_phi, sim_domain) %>%
    mutate(LB = NA, UB = NA)
  Smooth_df = map(1:N, function(x){
    output = data.frame(Est = mod$Y_hat[[x]], 
                        LB = mod$Y_low[[x]], 
                        UB = mod$Y_upp[[x]], 
                        Arg = sim_domain, 
                        Curve = paste0("Curve ", x))
    return(output)
  }) %>% list_rbind()
  
  return(list(Mu = Mu_df, EF = EF_df, Smooth = Smooth_df))
}
