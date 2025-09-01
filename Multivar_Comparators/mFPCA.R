# Function for aligning bootstrap-based FPC Confidence Intervals
align_CI_flips <- function(phi, anchor, LB, UB){
  K = ncol(anchor)
  out_LB = LB
  out_UB = UB
  for(k in 1:K){
    if(sum(phi[,k]*anchor[,k]) < 0){
      out_LB[,k] = -UB[,k]
      out_UB[,k] = -LB[,k]
    }
  }
  return(list(LB = out_LB, UB = out_UB))
}

# Format input data for mFPCA
mfpca_datalist <- function(df, N, P, domain){
  
  fun_dats = map(1:P, function(p){ 
    
    ppts = sort(unique(df$Subj))
    list_args = map(ppts, function(x){
      sargs = df %>%
               filter(Subj == x & Var == p) %>%
               pull(Arg)
      return(sargs)
    })
    list_obs = map(ppts, function(x){
      sy = df %>%
        filter(Subj == x & Var == p) %>%
        pull(Y)
      return(sy)
    })
      
    return(as.funData(irregFunData(argvals = list_args,
                                   X = list_obs)))
  })
  mf_fit = multiFunData(fun_dats)
  
  expan = map(1:P, function(p){
    return(list(type = "uFPCA"))
  })
  
  return(list(Data = mf_fit, Expansions = expan))
}

# Format mFPCA outputs for evaluation
extract_mfpca <- function(mod, N, P, domain, anchor, CI = F){
  
  # Where are the functions observed - handle indexing when observations missing
  M = length(domain)
  dom_idx = round(domain*(M-1))
  assignment = map(1:P, function(p){
    domain = mod$functions[[p]]@argvals[[1]]
    return(rep(p, length(domain)))
  }) %>% unlist()
  
  select_rows = map(1:P, function(p){
    sd_idx = round(mod$functions[[p]]@argvals[[1]]*(M-1))
    return(dom_idx %in% sd_idx)
  })
  
  # Fixed effects - all domains need updated indexing
  mu_df = map(1:P, function(p){
    est = mod$meanFunction[[p]]@X %>% as.vector()
    out_df = data.frame(Arg = domain[select_rows[[p]]], Var = p, 
                        Est = est, LB = NA, UB = NA)
    return(out_df)
  }) %>% list_rbind()
  
  # Eigenfunctions
  ef_mat_list = map(1:P, function(p){
    return(t(mod$functions[[p]]@X))
  })
  phi_est = do.call(rbind, ef_mat_list)
  
  sub_anchor = anchor[select_rows %>% unlist(),]
  aligned_est = procrust_FPC(list(phi_est), sub_anchor)[[1]]
  
  EF_df = map(1:P, function(p){
    out_df = FPC_df(aligned_est[assignment == p, ], domain[select_rows[[p]]])
    out_df$Var = p
    return(out_df)
  }) %>% list_rbind()
  
  if(CI){
    up_mat_list = map(1:P, function(p){
      return(t(mod$CI$alpha_0.05$upper[[p]]@X))
    })
    upp_est = do.call(rbind, up_mat_list)
    
    lo_mat_list = map(1:P, function(p){
      return(t(mod$CI$alpha_0.05$lower[[p]]@X))
    })
    low_est = do.call(rbind, lo_mat_list)
    
    aligned_CI = align_CI_flips(phi_est, sub_anchor, low_est, upp_est) 
    
    LB_df = map(1:P, function(p){
      out_df = FPC_df(aligned_CI$LB[assignment == p, ], domain[select_rows[[p]]]) %>%
        rename(LB = FPC_Val)
      out_df$Var = p
      return(out_df)
    }) %>% list_rbind()
    
    UB_df = map(1:P, function(p){
      out_df = FPC_df(aligned_CI$UB[assignment == p, ], domain[select_rows[[p]]]) %>%
        rename(UB = FPC_Val)
      out_df$Var = p
      return(out_df)
    }) %>% list_rbind()
    
    EF_df = EF_df %>%
      inner_join(LB_df) %>%
      inner_join(UB_df)
  }
  else{
    EF_df = EF_df %>%
      mutate(LB = NA, UB = NA)
  }
  
  # Predicted smooths
  Smooth_df = map(1:P, function(p){
    smooth_matrix = mod$fit[[p]]@X
    out_df = data.frame(t(smooth_matrix))
    colnames(out_df) = paste0("Curve ", 1:N)
    out_df$Arg = domain[select_rows[[p]]]
    out_df$Var = p
    out_df = out_df %>%
      pivot_longer(-c(Arg, Var), names_to = "Curve", values_to = "Est")
    return(out_df)
  }) %>% list_rbind() %>%
    mutate(LB = NA, UB = NA)
  
  return(list(FE = mu_df, EF = EF_df, Smooths = Smooth_df))
}
