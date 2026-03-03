# General Helper functions
library(cmdstanr)
source("Gen_Funcs/Libs.R")
source("Gen_Funcs/Bases.R")
source("Gen_Funcs/Generate.R")
source("Gen_Funcs/PostProcess.R")
source("Gen_Funcs/Comparisons.R")
source("Gen_Funcs/Convergence.R")

# Specific to each method
source("MSFAST/MSFAST_Help.R")

# For cmdstanr
extract <- function(fit_obj) {
  vars = fit_obj$metadata()$stan_variables
  draws = posterior::as_draws_rvars(fit_obj$draws())
  
  lapply(vars, \(var_name){  
    posterior::draws_of(draws[[var_name]], with_chains = FALSE)
  }) |> setNames(vars)
}

# Constants
{
  obs_ran = 3:7
  P = 3
  diff_obs = T
  M = 2000
  Q = 20
  K = 3
  sim_domain = seq(0, 1, length.out = M)
  quad_weights = booleQuad(sim_domain) 
  
  func_objs = Multivar_F(P, K, 4)
  MuFun = func_objs$Mu
  FPCs = func_objs$FPC
  EVals = func_objs$EV
  sigma2 = func_objs$Sig2
}

time_out = data.frame()
n_sim = 1#5
N_vals = c(100)#, 200, 300, 400, 500)

sinthr_mod = cmdstan_model("MSFAST/MSFAST.stan")
parvar_mod = cmdstan_model("MSFAST/MSFAST_Parallel.stan",
                           cpp_options = list(stan_threads = TRUE))
                         

for(N in N_vals){
  
  for(i in 1:n_sim){
    
    # Generate the data
    sim_dataset = gen_mFPCA(FPCs, N, P, sim_domain, MuFun, EVals, sigma2, obs_ran, diff_obs)
    sub_domain = sim_dataset$new_domain
    {
      mObjs = true_funcs(MuFun, FPCs, P, quad_weights, sim_domain)
      
      trueSmooth = data.frame(t(sim_dataset$Y_true))
      colnames(trueSmooth) = paste0("Curve ", 1:ncol(trueSmooth))
      trueSmooth$Arg = rep(sim_domain, P)
      trueSmooth$Var = rep(1:P, each = M)
      trueSmooth$wei = rep(quad_weights, P)
      trueSmooth = trueSmooth %>%
        pivot_longer(-c(Arg, Var, wei), names_to = "Curve", values_to = "Func") %>%
        left_join(mObjs$Mu %>% rename(Mu = Func)) %>%
        mutate(Func = Func + Mu) %>%
        select(-c(Mu))
    }
    
    # MSFAST - 2 cores
    {
      start.time = Sys.time()
      
      data_list = MSFAST_datalist(sim_dataset$ddf, N, K, Q, sub_domain, 
                                  basis_type = "B", scale = T)
      
      # Fit model
      mod = sinthr_mod$sample(data_list, chains = 2, parallel_chains = 2,
                              iter_warmup = 250, refresh = 0, 
                              iter_sampling = 1000, show_exceptions=F,
                              max_treedepth = 12, show_messages=F)
    
      # Complete timing
      end.time = Sys.time()
    }
    core2 = difftime(end.time, start.time, units = "secs") %>% as.numeric()
    
    # MSFAST - 4 cores
    {
      start.time = Sys.time()
      
      data_list = MSFAST_datalist(sim_dataset$ddf, N, K, Q, sub_domain, 
                                  basis_type = "B", scale = T)
      
      # Fit model
      mod = sinthr_mod$sample(data_list, chains = 4, parallel_chains = 4,
                              iter_warmup = 250, refresh = 0, 
                              iter_sampling = 500, show_exceptions=F,
                              max_treedepth = 12, show_messages=F)
      
      # Complete timing
      end.time = Sys.time()
    }
    core4 = difftime(end.time, start.time, units = "secs") %>% as.numeric()
    
    # MSFAST - 2 cores with 2 threads each
    {
      start.time = Sys.time()
      
      data_list = MSFAST_datalist(sim_dataset$ddf, N, K, Q, sub_domain, 
                                  basis_type = "B", scale = T)
      
      # Fit model
      mod = parvar_mod$sample(data_list, chains = 2, parallel_chains = 2,
                              refresh = 0, threads_per_chain = 2, 
                              iter_warmup = 250, iter_sampling = 1000, 
                              max_treedepth = 12, show_messages=F, 
                              show_exceptions=F)
      
      # Complete timing
      end.time = Sys.time()
    }
    c2t2 = difftime(end.time, start.time, units = "secs") %>% as.numeric()
    
    new_df = data.frame(Method = c("4 Cores", "2 Cores", "Threading"), 
                        Times = c(core4, core2, c2t2), 
                        Iter = i, N = N)
    time_out = rbind(time_out, new_df)
    write.csv(time_out, "Timing_Results/Parallel.csv")
  }
}

# RHat calculations

# Extract
# new_B = FAST_B(basis_type = "B", Q, sim_domain)
# objects = FAST_extract(mod, new_B, data_list)

# psi_hat = pmap(list(objects$Weights, objects$Score), function(w, s){
# return(s %*% t(w))
# }) %>% abind(along = 3) %>%
#   apply(c(1,2), mean)
# kB = kronecker(diag(1, P), new_B)
# phi_est = kB %*% svd(psi_hat, nv = K)$v

# RHat_FAST(mod, data_list, new_B, phi_est)



