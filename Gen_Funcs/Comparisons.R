# Function to calculate equally-spaced quadrature weights
booleQuad <- function(dx){
  n = length(dx)
  h = diff(dx)[1]
  weights = rep(NA, n)
  weights[1] = 7
  weights[n] = 7
  weights[seq.int(1, n-1, 2)] = 32
  weights[seq.int(2, n-2, 4)] = 12
  weights[seq.int(4, n-4, 4)] = 14
  return(2*h/45*weights)
}

# Summary of FE performance
FE_comp <- function(trueFE_df, sampleFE_df, method){
  M_t = n_distinct(trueFE_df$Arg)
  M_o = n_distinct(sampleFE_df$Arg)
  
  jdf = inner_join(trueFE_df, sampleFE_df, by = c("Var", "Arg")) %>%
    mutate(Cov = LB <= Func & UB >= Func, 
           SE = (Func - Est)^2) %>%
    group_by(Var) %>%
    summarize(Cov = mean(Cov), ISE = sum(wei * SE) * (M_t/M_o))
  jdf$Method = method
  
  return(jdf)
}

# Summary of EF performance
EF_comp <- function(trueEF_df, sampleEF_df, method){
  M_t = n_distinct(trueEF_df$Arg)
  M_o = n_distinct(sampleEF_df$Arg)
  
  jdf = inner_join(trueEF_df, sampleEF_df, by = c("Var", "Arg", "FPC_Num")) %>%
    mutate(Cov = LB <= Func & UB >= Func, 
           SE = (Func - FPC_Val)^2) %>%
    group_by(Var, FPC_Num) %>%
    summarize(Cov = mean(Cov), ISE = sum(wei * SE) * (M_t/M_o))
  jdf$Method = method
  return(jdf)
}

# Summary of Smooth performance
Smooth_comp <- function(true_Smooth, sample_Smooth, method){
  M_t = n_distinct(true_Smooth$Arg)
  M_o = n_distinct(sample_Smooth$Arg)
  
  jdf = inner_join(true_Smooth, sample_Smooth, by = c("Var", "Arg", "Curve")) %>%
    mutate(Cov = LB <= Func & UB >= Func, 
           SE = (Func - Est)^2) %>%
    group_by(Var, Curve) %>%
    summarize(Cov = mean(Cov), ISE = sum(wei * SE) * (M_t/M_o))
  jdf$Method = method
  return(jdf)
}
