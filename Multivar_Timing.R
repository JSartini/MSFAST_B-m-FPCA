# General Helper functions
source("Gen_Funcs/Libs.R")
source("Gen_Funcs/Bases.R")
source("Gen_Funcs/Generate.R")
source("Gen_Funcs/PostProcess.R")
source("Gen_Funcs/Comparisons.R")

# Specific to each method
source("MSFAST/MSFAST_Help.R")
source("Multivar_Comparators/mFACEs.R")
source("Multivar_Comparators/mFPCA.R")
source("Multivar_Comparators/vmp.R")

results_directory = "Results"

create_dir <- function(direct){
  if (!file.exists(direct)){
    dir.create(direct, recursive = T)
  }
}

args = commandArgs(trailingOnly=TRUE)

n_sim = as.numeric(args[3])
obs_ran = 3:7
N = as.numeric(args[4])
diff_obs = T
M = 2000
Q = 20
K = 3
P = 3
sim_domain = seq(0, 1, length.out = M)
quad_weights = booleQuad(sim_domain)
Sys.setenv(STAN_NUM_THREADS=3) # Enable STAN threading

func_objs = Multivar_F(P, K, 4)
MuFun = func_objs$Mu
FPCs = func_objs$FPC
EVals = func_objs$EV
sigma2 = func_objs$Sig2

# Create output directories
{
  out_dir = paste0(args[1], "_N", N , "_M", M)
  time_dir = paste0(results_directory, "/Time/", out_dir)
  create_dir(time_dir)
}

message("Output directory: ", args[1])
message("Output file: ", args[2])
message("Number of simulations: ", args[3])
message("Number of time series: ", args[4])

time_out = data.frame(Method = c(), Timing = c(), Sample = c())

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
    data_list = MSFAST_datalist(sim_dataset$ddf, N, K, Q, sub_domain, 
                                basis_type = "B", scale = T)
    
    # Fit model
    mod = stan(
      file = "MSFAST/MSFAST_Parallel.stan",
      data = data_list, 
      chains = 1, 
      cores = 1, 
      warmup = 1000, 
      iter = 2000, 
      control = list(max_treedepth = 12),
      verbose = F,
      refresh = 0
    )
    
    # Extract
    new_B = FAST_B(basis_type = "B", Q, sim_domain)
    objects = FAST_extract(mod, new_B, data_list)
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
    fast_time = difftime(end.time, start.time, units = "secs") %>% as.numeric()
  }
  
  # mFACEs (Li and Xiao 2020)
  {
    start.time = Sys.time()
    
    DL = mfaces_datalist(sim_dataset$ddf, N, P, sim_domain)
    
    face_fit = mface.sparse(DL$fit, argvals.new = sim_domain,
                            knots.option = "quantile", center = T, pve = 0.99) 
    
    FACE_objs = extract_mfaces(face_fit, sim_domain, N, P, mObjs$mFPC, DL$fit, DL$new)
    
    # Complete timing
    end.time = Sys.time()
    face_time = difftime(end.time, start.time, units = "secs") %>% as.numeric()
  }
  
  # mFPCA (Happ and Greven 2018)
  {
    start.time = Sys.time()
    
    DL = mfpca_datalist(sim_dataset$ddf, N, P, sim_domain)
    
    fpca_fit = MFPCA(DL$Data, M = K, DL$Expansions, fit = T)
    FPCA_objs = extract_mfpca(fpca_fit, N, P, sim_domain, mObjs$mFPC, CI = F)
    
    # Complete timing
    end.time = Sys.time()
    mfpca_time = difftime(end.time, start.time, units = "secs") %>% as.numeric()
  }
  
  # VMP (T.H. Nolan 2025)
  {
    start.time = Sys.time()
    
    DL = vmp_datalist(sim_dataset$ddf, N, P)
    
    vmp_fit = run_vmp_fpca(DL$t, DL$y, K, Q, verbose = F, time_g = sim_domain)
    
    VMP_objs = extract_vmp(vmp_fit, sim_domain, N, K, P, mObjs$mFPC)
    
    # Complete timing
    end.time = Sys.time()
    vmp_time = difftime(end.time, start.time, units = "secs") %>% as.numeric()
  }
  
  # Collate timing results
  {
    time_df = data.frame(Method = c("MSFAST", "mFACEs", "VMP", "mFPCA"),
                         Timing = c(fast_time, face_time, vmp_time, mfpca_time), 
                         Sample = x)
  }
  
  # Add most recent result
  time_out = rbind(time_out, time_df)
  
  # Write to storage
  write.csv(time_out, paste0(time_dir, "/", args[2], ".csv"))
}
