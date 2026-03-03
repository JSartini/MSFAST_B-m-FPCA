# Produce P_alpha using functional basis - numerical integration
FAST_P <- function(derivs, Q, upp){
  P2 = matrix(0, ncol = Q, nrow = Q)
  for(i in 1:Q){
    for(j in 1:i){
      # Rescale for numerical precision purposes
      f2 <- function(x){ return(derivs[[i]](x) * derivs[[j]](x)) }
      P2[i,j] = stats::integrate(f2, lower = 0, upper = upp,
                                 subdivisions = 10000)$value
      if(i != j){
        P2[j,i] = P2[i,j]
      }
    }
  } 
  return(P2)
}

# Evaluate the orthogonal spline basis at specified time points
FAST_B <- function(basis_type = "B", Q, Domain){
  if(basis_type == "B"){
    B_f = OBasis(c(rep(0, 3), seq(0, 1, length.out = Q-2), rep(1, 3)))
    B = evaluate(B_f, Domain)
  }
  else{
    if(basis_type == "Fourier"){
      basis = Fourier_bases(Q)
    }
    else if(basis_type == "Legendre"){
      basis = Legendre_bases(Q)
    }
    else if(basis_type == "Splinet"){
      basis = Splinet_bases(Q)$B
    }
    else{
      stop("Basis not supported")
    }
    B = map(1:Q, function(q){
      return(basis[[q]](Domain))
    }) %>% abind(along = 2)
  }
  return(B)
}

# Produce inputs for FAST
MSFAST_datalist <- function(df, N, K, Q, Domain, ord = "Var", basis_type = "B",
                            alpha = 0.1, scale = T, threads = 1){
  
  # Generate spline basis with derivatives
  P0 = diag(Q)
  if(basis_type == "B"){
    B_f = OBasis(c(rep(0, 3), seq(0, 1, length.out = Q-2), rep(1, 3)))
    B = evaluate(B_f, Domain)
    P2 = OuterProdSecondDerivative(B_f)
  }
  else{
    B = FAST_B(basis_type, Q, Domain)
    upp = 1
    if(basis_type == "Fourier"){
      basis = Fourier_bases(Q)
      derivs = Fourier_d2(Q)
    }
    else if(basis_type == "Splinet"){
      bobj = Splinet_bases(Q)
      basis = bobj$B
      derivs = Splinet_d2(Q, bobj$cInt, bobj$cSlo)
      upp = 10 # Larger range for stability of numerical integration
    }
    else if(basis_type == "Legendre"){
      basis = Legendre_bases(Q)
      derivs = Legendre_d2(Q)
    }
    else{
      stop("Basis not supported")
    }
    P2 = FAST_P(derivs, Q, upp)
    P_alpha = alpha * diag(Q) + (1-alpha) * P2
  }
  P_alpha = alpha * P0 + (1-alpha) * P2
  
  # Scale the Y-values, storing the location and scale
  if(scale){
    const_df = df %>%
      group_by(Var) %>%
      summarize(mu_Y = mean(Y), 
                sd_Y = sd(Y))
  }
  else{
    const_df = df %>%
      group_by(Var) %>%
      summarize(mu_Y = 0, 
                sd_Y = 1)
  }
  M_val = length(Domain)
  
  if(ord == "Var"){
    s_df = data.frame(S_new = 1:length(Domain), 
                      Arg = sort(Domain))
    
    arranged_df = df %>%
      arrange(Var, Arg) %>%
      left_join(const_df, by = "Var") %>%
      left_join(s_df) %>%
      mutate(Y = (Y - mu_Y)/sd_Y) %>%
      select(-c(mu_Y, sd_Y))  
    
    TPC = arranged_df %>%
      group_by(Var) %>%
      summarize(m_c = n()) %>%
      pull(m_c)
    dim(TPC) = c(P)
    
    # Format for FAST
    fast_list = list(N = N, M = M_val, L = nrow(arranged_df), Q = Q, K = K, P = length(TPC), 
                     Y = arranged_df$Y, Subj = arranged_df$Subj, S = arranged_df$S_new,
                     Tp_card = TPC, B = B, P_alpha = P_alpha, consts = const_df)
  }
  else if(ord == "ID"){
    s_df = data.frame(S_new = 1:length(Domain), 
                      Arg = sort(Domain))
    
    arranged_df = df %>%
      arrange(Subj, Var, Arg) %>%
      left_join(const_df, by = "Var") %>%
      left_join(s_df) %>%
      mutate(Y = (Y - mu_Y)/sd_Y) %>% #, 
             # S_new = S_new + (Var - 1)*M_val) %>%
      select(-c(mu_Y, sd_Y))  
    
    IPC = matrix(NA, nrow = N, ncol = P)
    for(p in 1:P){
      IPC[,p] = arranged_df %>%
        filter(Var == p) %>%
        group_by(Subj) %>%
        summarize(m_c = n()) %>%
        pull(m_c)
    }
    
    # Format for FAST
    fast_list = list(N = N, M = M_val, L = nrow(arranged_df), Q = Q, K = K, 
                     P = n_distinct(arranged_df$Var), Y = arranged_df$Y, 
                     S = arranged_df$S_new, ID_card = IPC, B = B, 
                     P_alpha = P_alpha, consts = const_df)
    fast_list$Th = threads
  }
  else{
    stop("Ordering for parallelization not supported")
  }
  return(fast_list)
}

# Place EF and scores in standard form for post-processing
FAST_extract <- function(mod_fit, B, DL){
  samples = extract(mod_fit)
  n_samp = dim(samples$sigma2)[1]
  kB = kronecker(diag(DL$P), B)
  
  Weight_list = map(1:n_samp, function(x){
    return(samples$Psi[x,,])
  })
  
  EF_list = map(Weight_list, function(weights){
    EF_sample = kB %*% weights
    return(EF_sample)
  })
  
  Score_list = map(1:n_samp, function(x){
    Score_sample = samples$Scores[x,,]
    return(Score_sample)
  })
  
  Mu_list = map(1:n_samp, function(x){
    return((kB %*% samples$w_mu[x,]))
  })
  
  return(list(Weights = Weight_list, EF = EF_list, 
              Score = Score_list, Mu = Mu_list))
}

# Extract EF and Scores by-chain
FAST_byChain <- function(mod_fit, N, M, Q, K, P){
  chain_samples = extract(mod_fit, permuted = F)
  dimen_names = dimnames(chain_samples)$parameters
  
  n_samp = dim(chain_samples)[1]
  n_chain = dim(chain_samples)[2]
  
  Psi_idx = grepl("Psi", dimen_names)
  Score_idx = grepl("Scores", dimen_names)
  
  Psi_byChain = list()
  Score_byChain = list()
  
  for(j in 1:n_chain){
    Psi_byChain[[j]] = map(1:n_samp, function(i){
      Psi = chain_samples[i,j,Psi_idx]
      dim(Psi) = c(P*Q, K)
      return(Psi)
    })
    Score_byChain[[j]] = map(1:n_samp, function(i){
      Xi = chain_samples[i,j,Score_idx]
      dim(Xi) = c(N, K)
      return(Xi)
    })
  }
  return(list(Psi = Psi_byChain, Score = Score_byChain))
}
