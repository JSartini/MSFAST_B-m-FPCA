# General Helper functions
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

args = commandArgs(trailingOnly=TRUE)

n_sim = as.numeric(args[3])
SNR = as.numeric(args[5])
obs_ran = 3:7
P = 3
diff_obs = T
N = 100
M = 100
Q = 20
K = 3
sim_domain = seq(0, 1, length.out = M)
quad_weights = booleQuad(sim_domain)

func_objs = Multivar_F(P, K, SNR)
MuFun = func_objs$Mu
FPCs = func_objs$FPC
EVals = func_objs$EV
sigma2 = func_objs$Sig2

out_dir = paste0(args[1], "_K", args[4])

# Create output directories
{
  FE_dir = paste0(results_directory, "/FE/", out_dir)
  EF_dir = paste0(results_directory, "/EF/", out_dir)
  Smooth_dir = paste0(results_directory, "/Smooth/", out_dir)
  
  create_dir(FE_dir)
  create_dir(EF_dir)
  create_dir(Smooth_dir)
}

message("Output directory: ", args[1])
message("Output Directory: ", out_dir)
message("Output file: ", args[2])
message("Number of simulations: ", args[3])
message("Number of FPCs fit: ", args[4])

fit_K = as.numeric(args[4])
sub_K = min(K, fit_K)

FE_out = data.frame(Method = c(), Var = c(), Cov = c(), ISE = c(), Sample = c(),
                    P = c(), NObs = c())
EF_out = data.frame(Method = c(), Var = c(), FPC_Num = c(), Cov = c(), 
                    ISE = c(), Sample = c(), P = c(), NObs = c())
Smooth_out = data.frame(Method = c(), Var = c(), RISE = c(), Cov = c(), 
                        Sample = c(), P = c(), NObs = c())

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
    # Collate data
    data_list = MSFAST_datalist(sim_dataset$ddf, N, K = fit_K, Q,
                                sub_domain, basis_type = "B", scale = T)
    
    # Fit model
    mod = stan(
      file = "../MSFAST/MSFAST.stan",
      data = data_list, 
      chains = 4, 
      cores = 4, 
      warmup = 500, 
      iter = 1000, 
      # init = init_list,
      control = list(max_treedepth = 12,
                     adapt_delta = 0.85),
      verbose = F,
      refresh = 0
    )
    
    # Extract
    new_B = FAST_B(basis_type = "B", Q, sim_domain)
    objects = FAST_extract(mod, new_B, data_list)
    align = procrust_WEI(objects$Weights, new_B, P, mObjs$mFPC[,1:sub_K])
    
    Mu_df = FE_Summary(objects$Mu, sim_domain, P) %>%
      left_join(data_list$consts, by = "Var") %>%
      mutate(Est = Est * sd_Y + mu_Y, 
             LB = LB * sd_Y + mu_Y, 
             UB = UB * sd_Y + mu_Y)
    
    # EF_ests = FPC_Est_WEI(objects$Weights, new_B, P, sim_domain, mObjs$mFPC)
    EF_ests = Latent_RSV(objects$Weights, objects$Score, new_B, P, sim_domain, mObjs$mFPC[,1:sub_K])
    EF_df = inner_join(FPC_CI(align$EF, sim_domain, P), 
                       EF_ests)
    
    Smooth_df = Smooth_Summary(N, P, objects$Mu, objects$EF, 
                               objects$Score, sim_domain) %>%
      left_join(data_list$consts, by = "Var") %>%
      mutate(Est = Est * sd_Y + mu_Y, 
             LB = LB * sd_Y + mu_Y, 
             UB = UB * sd_Y + mu_Y) %>% 
      select(-c(mu_Y, sd_Y))
    
    # Calculate performance
    FE_FAST =  FE_comp(mObjs$Mu, Mu_df, "MSFAST") %>%
      mutate(Sample = x)
    EF_FAST = EF_comp(mObjs$FPC, EF_df, "MSFAST") %>%
      mutate(Sample = x)
    Smooth_FAST = Smooth_comp(trueSmooth, Smooth_df, "MSFAST") %>%
      mutate(Sample = x)
  }
  
  # Simple mean functions as smooth
  {
    Smooth_df = map(1:N, function(x){
      person_data = map(1:P, function(p){
        sub_data = sim_dataset$ddf %>%
          filter(Subj == x & Var == p) %>%
          pull(Y)
        
        return(data.frame(Curve = paste0("Curve ", x),
                          Var = p,
                          Arg = sim_domain, 
                          Est = mean(sub_data),
                          UB = NA, LB = NA))
      }) %>% list_rbind()
      return(person_data)
    }) %>% list_rbind()
    
    Smooth_mean = Smooth_comp(trueSmooth, Smooth_df, "Mean") %>%
      mutate(Sample = x)
  }
  
  # Collate model component results
  {
    FE_df = FE_FAST %>%
      mutate(P = P, NObs = "3-7", K = fit_K)
    EF_df = EF_FAST %>%
      mutate(P = P, NObs = "3-7", K = fit_K)
  }
  
  # Collate prediction results (RISE)
  {
    Smooth_df = left_join(Smooth_FAST, Smooth_mean %>% 
                            rename(BaseISE = ISE) %>%
                            select(Var, Curve, BaseISE)) %>%
      mutate(RISE = ISE/BaseISE) %>%
      select(-c(BaseISE, ISE)) %>%
      group_by(Method, Var) %>%
      summarize(RISE = mean(RISE), Cov = mean(Cov), 
                Sample = first(Sample)) %>%
      mutate(P = P, NObs = "3-7", K = fit_K)
  }
  
  # Add most recent result
  FE_out = rbind(FE_out, FE_df)
  EF_out = rbind(EF_out, EF_df)
  Smooth_out = rbind(Smooth_out, Smooth_df)
  
  # Write to storage
  write.csv(FE_out, paste0(FE_dir, "/", args[2], ".csv"))
  write.csv(EF_out, paste0(EF_dir, "/", args[2], ".csv"))
  write.csv(Smooth_out, paste0(Smooth_dir, "/", args[2], ".csv"))
}


