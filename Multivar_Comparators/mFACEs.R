# Format input data for mFACEs
mfaces_datalist <- function(df, N, P, domain){
  
  for_fit = map(1:P, function(p){
    sub_df = df %>%
      filter(Var == p) %>%
      select(Arg, Subj, Y) %>%
      rename(argvals = Arg, 
             subj = Subj, y = Y) %>%
      mutate(subj = as.numeric(subj)) %>%
      as.data.frame()
    return(sub_df)
  })
  names(for_fit) = paste0("y", 1:P)
  
  new_data = map(1:P, function(p){
    append = expand_grid(1:N, domain)
    colnames(append) = c("subj", "argvals")
    append$y = NA
    append = append %>%
      arrange(subj, argvals) # %>%
      # anti_join(for_fit[[p]], by = c("subj", "argvals"))
    return(append)
  })
  names(new_data) = paste0("y", 1:P)

  return(list(fit = for_fit, new = new_data))
}

# Format mFACEs outputs for evaluation
extract_mfaces <- function(mod, domain, N, P, anchor, fit_data, pred_data){
  # Mean functions
  df_lens = rep(NA, P)
  sub_pred = map(1:P, function(p){
    data_df = fit_data[[p]] %>%
      filter(subj == 1)
    df_lens[p] <<- nrow(data_df)
    
    pred_df = pred_data[[p]] %>%
      filter(subj == 1)
    
    return(rbind(data_df,pred_df))
  })
  
  pred_mod = predict(mod, sub_pred)
  
  mu_df = map(1:P, function(p){
    est = pred_mod$mu.pred[[p]][-(1:df_lens[p])]
    return(data.frame(Arg = domain, Var = p, 
                      Est = est, LB = NA, UB = NA))
  }) %>% list_rbind()
  
  # Eigenfunctions
  align_FPC = procrust_FPC(list(mod$eigenfunctions[,1:ncol(anchor)]), anchor)[[1]]
  EF_df = FPC_df(align_FPC, domain, P) %>%
    mutate(LB = NA, UB = NA)
  
  # Predicted smooths
  size_group = 1
  sub_groups = split(1:N, ceiling(seq_along(1:N)/size_group))
  Smooth_df = map(sub_groups, function(idxs){
    
    df_lens = rep(NA, P)
    sub_pred = map(1:P, function(p){
      data_df = fit_data[[p]] %>%
        filter(subj %in% idxs)
      df_lens[p] <<- nrow(data_df)
      
      pred_df = pred_data[[p]] %>%
        filter(subj %in% idxs)
      
      return(rbind(data_df,pred_df))
    })
    
    pred_mod = predict(mod, sub_pred)
    
    smooth_df = map(1:P, function(p){
      est = pred_mod$y.pred[[p]][-(1:df_lens[p])]
      se = pred_mod$se.pred[[p]][-(1:df_lens[p])]
      out_df = data.frame(Arg = rep(domain, size_group), 
                          Var = p, 
                          Curve = rep(paste0("Curve ", idxs), each = length(domain)), 
                          Est = est, 
                          LB = est - 1.96*se, 
                          UB = est + 1.96*se)
      return(out_df)
    }) %>% list_rbind()
    
    return(smooth_df)
  }) %>% list_rbind()
  
  return(list(EF = EF_df, Smooths = Smooth_df, FE = mu_df))
}
