# General Helper functions
library(cmdstanr)
source("Gen_Funcs/Libs.R")
source("Gen_Funcs/Bases.R")
source("Gen_Funcs/Generate.R")
source("Gen_Funcs/PostProcess.R")
source("Gen_Funcs/Comparisons.R")

# Specific to each method
source("MSFAST/MSFAST_Help.R")
source("Multivar_Comparators/Univar.R")
source("Multivar_Comparators/mFACEs.R")
source("Multivar_Comparators/mFPCA.R")
source("Multivar_Comparators/vmp.R")

results_directory = "Timing_Results"

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
obs_ran = 3:7
N = as.numeric(args[4])
P = as.numeric(args[5])
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

# Create output directories
{
  out_dir = paste0(args[1], "_N", N , "_P", P)
  time_dir = paste0(results_directory, "/Time/", out_dir)
  create_dir(time_dir)
}

# Compile the cmdstan files
sinthr_mod = cmdstan_model("MSFAST/MSFAST.stan")

# Messages
message("Output directory: ", args[1])
message("Output file: ", args[2])
message("Number of simulations: ", args[3])
message("Number of time series: ", args[4])
message("Number of variables: ", args[5])

time_out = data.frame(Method = c(), Timing = c(), Sample = c(), N = c(), P = c())

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
  
  # Univariate FACE separately (Xiao and Li 2017)
  {
    start.time = Sys.time()

    DL = ufpca_datalist(sim_dataset$ddf, N, P, sim_domain)

    fits = map(1:P, function(p){
      fit = face.sparse(data = DL$fit[[p]], argvals.new = sim_domain,
                        pve = 0.999, center = T, calculate.scores = T)
      return(fit)
    })

    uFACE_objs = extract_ufpca(fits, sim_domain, N, P, mObjs$mFPC, DL)

    # Complete timing
    end.time = Sys.time()
    univar_time = difftime(end.time, start.time, units = "secs") %>% as.numeric()
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
  
  # Mean Field Approximation (T.H. Nolan 2025)
  {
    start.time = Sys.time()
    
    DL = vmp_datalist(sim_dataset$ddf, N, P)
    
    mf_fit = run_mfvb_fpca(DL$t, DL$y, K, Q, verbose = F, time_g = sim_domain)
    
    MF_objs = extract_vmp(mf_fit, sim_domain, N, K, P, mObjs$mFPC)
    
    # Complete timing
    end.time = Sys.time()
    mf_time = difftime(end.time, start.time, units = "secs") %>% as.numeric()
  }
  
  # MSFAST - single core
  {
    start.time = Sys.time()
    
    # Collate data
    data_list = MSFAST_datalist(sim_dataset$ddf, N, K, Q, sub_domain, 
                                basis_type = "B", scale = T)
    
    # Fit model
    cmd_mod = sinthr_mod$sample(data_list, chains = 1, parallel_chains = 1,
                                iter_warmup = 1000, refresh = 0, 
                                iter_sampling = 1000, show_exceptions=F,
                                max_treedepth = 12, show_messages=F)
    
    # Extract
    new_B = FAST_B(basis_type = "B", Q, sim_domain)
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
                               objects$Score, sim_domain, parallel = F) %>%
      left_join(data_list$consts, by = "Var") %>%
      mutate(Est = Est * sd_Y + mu_Y, 
             LB = LB * sd_Y + mu_Y, 
             UB = UB * sd_Y + mu_Y) %>% 
      select(-c(mu_Y, sd_Y))
    
    # Complete timing
    end.time = Sys.time()
    fast_time = difftime(end.time, start.time, units = "secs") %>% as.numeric()
  }
  
  # Collate timing results
  {
    time_df = data.frame(Method = c("MSFAST", "MSFAST Parallel - Variates", "mFACEs", "uFPCA",
                                     "mFPCA", "VMP", "MF"), # "MSFAST Parallel - Subjects",
                         Timing = c(fast_time, fast_parvar_time, face_time, univar_time,
                                    mfpca_time, vmp_time, mf_time), # fast_parsub_time, 
                         Sample = x, N = N, P = P)
  }
  
  # Add most recent result
  time_out = rbind(time_out, time_df)
  
  # Write to storage
  write.csv(time_out, paste0(time_dir, "/", args[2], ".csv"))
}
