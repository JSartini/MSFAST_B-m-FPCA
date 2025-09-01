# Format input data for FACE
FACE_datalist <- function(sim_df){
  output_df = data.frame(subj = sim_df$Subj, 
                         argvals = sim_df$Arg, 
                         y = sim_df$Y)
  return(output_df)
}

# Format FACE outputs for evaluation
FACE_extract <- function(mod, sim_domain, anchor, N, M, K, fpca_df){
  # Align with Procrustes
  face_phi = procrust_FPC(list(mod$eigenfunctions), anchor)[[1]]
  
  Mu_df = data.frame(Arg = sim_domain, Est = mod$mu.new, LB = NA, UB = NA)
  EF_df = FPC_df(face_phi, sim_domain) %>%
    mutate(LB = NA, UB = NA)

  Smooth_df = map(1:N, function(chosen_subj){
    pred_df = fpca_df %>%
      rbind(data.frame(subj = rep(chosen_subj, length(sim_domain)), 
                       argvals = sim_domain, 
                       y = NA))
    face_pred = predict(mod, newdata = pred_df)
    sub_smooth = face_pred$newdata[-(1:nrow(fpca_df)),] %>%
      select(-c(y)) %>%
      dplyr::rename(Arg = argvals, Curve = subj) %>%
      mutate(Curve = paste0("Curve ", Curve))
    sub_smooth$Est = face_pred$y.pred[-(1:nrow(fpca_df))]
    sub_smooth$SE = face_pred$se.pred[-(1:nrow(fpca_df))]
    return(sub_smooth)
  }) %>% list_rbind() %>%
    mutate(UB = Est + 1.96*SE, LB = Est - 1.96*SE) %>%
    select(-c(SE))
  
  return(list(EFm = face_phi, Muv = mod$mu.new, Mu = Mu_df, EF = EF_df,
              Smooth = Smooth_df))
}
