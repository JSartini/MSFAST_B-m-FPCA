
cond_datalist <- function(N, K, Y, Subj, S, Mu, FPCm){
  Y_res = Y - Mu[S]
  FPC_eval = FPCm[S,]
  return(list(N = N, L = length(Y), K = K, Y = Y_res, 
              ts_idx = Subj, FPC = FPC_eval))
}

cond_extract <- function(mod, N, K, FPCm, sim_domain){
  samples = extract(mod)
  n_sample = dim(samples$Scores)[1]
  
  smooth_df = map(1:n_sample, function(x){
    smooths = data.frame(FPCm %*% t(samples$Scores[x,,]))
    colnames(smooths) = paste0("Curve ", 1:N)
    smooths$Arg = sim_domain
    smooths$iter = x
    smooths = smooths %>%
      pivot_longer(-c(Arg, iter), names_to = "Curve", values_to = "Func")
    return(smooths)
  }) %>% list_rbind() %>%
    group_by(Curve, Arg) %>%
    summarize(Est = mean(Func), 
              LB = quantile(Func, probs = c(0.025)), 
              UB = quantile(Func, probs = c(0.975)))
  
  return(list(Smooth = smooth_df))
}
