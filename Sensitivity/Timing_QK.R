# General Helper functions
library(cmdstanr)
source("../Gen_Funcs/Libs.R")
source("../Gen_Funcs/Bases.R")
source("../Gen_Funcs/Generate.R")
source("../Gen_Funcs/PostProcess.R")
source("../Gen_Funcs/Comparisons.R")

# Specific to each method
source("../MSFAST/MSFAST_Help.R")

results_directory = "Results"

create_dir <- function(direct){
  if (!file.exists(direct)){
    dir.create(direct, recursive = T)
  }
}

extract <- function(fit_obj) {
  vars = fit_obj$metadata()$stan_variables
  draws = posterior::as_draws_rvars(fit_obj$draws())
  
  lapply(vars, \(var_name){  
    posterior::draws_of(draws[[var_name]], with_chains = FALSE)
  }) |> setNames(vars)
}

args = commandArgs(trailingOnly=TRUE)

n_sim = as.numeric(args[3])
SNR = as.numeric(args[6])
obs_ran = 3:7
P = 3
diff_obs = T
N = 100
M = 100
K = 3
sim_domain = seq(0, 1, length.out = M)
quad_weights = booleQuad(sim_domain)

func_objs = Multivar_F(P, K, SNR)
MuFun = func_objs$Mu
FPCs = func_objs$FPC
EVals = func_objs$EV
sigma2 = func_objs$Sig2

out_dir = paste0(args[1], "_Q", args[4], "_K", args[5])

# Create output directories
{
  Time_dir = paste0(results_directory, "/Time/", out_dir)
  
  create_dir(Time_dir)
}

# Compile the cmdstan files
sinthr_mod = cmdstan_model("../MSFAST/MSFAST.stan")

message("Output directory: ", args[1])
message("Output Directory: ", out_dir)
message("Output file: ", args[2])
message("Number of simulations: ", args[3])
message("Number of splines Q: ", args[4])
message("Number of FPCs fit K: ", args[5])

fit_Q = as.numeric(args[4])
fit_K = as.numeric(args[5])

Time_out = data.frame(Seconds = c(), Method = c(), Sample = c(), K = c(), 
                      Q = c())

for(x in 1:n_sim){
  print(paste0("Iteration: ", x))
  
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
  
  # MSFAST
  {
    start.time = Sys.time()
    
    # Collate data
    data_list = MSFAST_datalist(sim_dataset$ddf, N, fit_K, fit_Q, sub_domain, 
                                basis_type = "B", scale = T)
    
    # Fit model
    cmd_mod = sinthr_mod$sample(data_list, chains = 1, parallel_chains = 1,
                                iter_warmup = 1000, refresh = 0, 
                                iter_sampling = 1000, show_exceptions=F,
                                max_treedepth = 12, show_messages=F)
    
    # Extract
    new_B = FAST_B(basis_type = "B", fit_Q, sim_domain)
    objects = FAST_extract(cmd_mod, new_B, data_list)
    align = procrust_WEI(objects$Weights, new_B, P, mObjs$mFPC)
    
    Mu_df = FE_Summary(objects$Mu, sim_domain, P) %>%
      left_join(data_list$consts, by = "Var") %>%
      mutate(Est = Est * sd_Y + mu_Y, 
             LB = LB * sd_Y + mu_Y, 
             UB = UB * sd_Y + mu_Y)
    
    EF_ests = FPC_Est_WEI(align$Weights, new_B, P, sim_domain, mObjs$mFPC)
    EF_df = inner_join(FPC_CI(align$EF, sim_domain, P), 
                       EF_ests)
    
    Smooth_df = Smooth_Summary(N, P, objects$Mu, objects$EF, 
                               objects$Score, sim_domain) %>%
      left_join(data_list$consts, by = "Var") %>%
      mutate(Est = Est * sd_Y + mu_Y, 
             LB = LB * sd_Y + mu_Y, 
             UB = UB * sd_Y + mu_Y) %>% 
      select(-c(mu_Y, sd_Y))
    
    # Complete timing
    end.time = Sys.time()
    fastS_time = difftime(end.time, start.time, units = "secs") %>% as.numeric()
  }
  
  # Collate
  {
    Time_df = data.frame(Seconds = fastS_time,
                         Method = "FAST",
                         Sample = x, K = fit_K, Q = fit_Q)
  }
  
  # Add most recent result
  Time_out = rbind(Time_out, Time_df)
  
  # Write to storage
  write.csv(Time_out, paste0(Time_dir, "/", args[2], ".csv"))
}