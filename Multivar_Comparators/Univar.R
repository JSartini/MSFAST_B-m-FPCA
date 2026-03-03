# Format input data for univariate FPCAs
ufpca_datalist <- function(df, N, P, domain){
  
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

# Extract univariate FPCA results
extract_ufpca <- function(mods, domain, N, P, anchor, datas){
  # Mean functions
  mu_df = map(1:P, function(p){
    return(data.frame(Arg = domain, Var = p, Est = mods[[p]]$mu.new,
                      LB = NA, UB = NA))
  }) %>% list_rbind()
  
  # Eigenfunctions
  FPC_stacked = map(1:P, function(p){
    return(mods[[p]]$eigenfunctions[,1:ncol(anchor)])
  })
  FPC_stacked = do.call(rbind, FPC_stacked)
  align_FPC = procrust_FPC(list(FPC_stacked), anchor)[[1]]
  EF_df = FPC_df(align_FPC, domain, P) %>%
    mutate(LB = NA, UB = NA)

  # Predicted smooths
  size_group = 1
  sub_groups = split(1:N, ceiling(seq_along(1:N)/size_group))
  Smooth_df = map(sub_groups, function(idxs){
    sub_df = map(1:P, function(p){
      for_predict = rbind(datas$fit[[p]], 
                          datas$new[[p]] %>% 
                            filter(subj %in% idxs))
      
      pred_obj = predict(mods[[p]], newdata = for_predict)
      est = pred_obj$y.pred[-(1:nrow(datas$fit[[p]]))]
      se = pred_obj$se.pred[-(1:nrow(datas$fit[[p]]))]
      return(data.frame(Arg = domain,
                        Var = p, 
                        Curve = rep(paste0("Curve ", idxs), each = length(domain)), 
                        Est = est, 
                        LB = est - 1.96*se, 
                        UB = est + 1.96*se))
    }) %>% list_rbind()
    
    return(sub_df)
  }) %>% list_rbind()
  
  return(list(EF = EF_df, Smooths = Smooth_df, FE = mu_df))
}
