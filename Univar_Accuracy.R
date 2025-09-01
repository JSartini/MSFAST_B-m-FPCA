# General Helper functions
source("Gen_Funcs/Libs.R")
source("Gen_Funcs/Bases.R")
source("Gen_Funcs/Generate.R")
source("Gen_Funcs/PostProcess.R")
source("Gen_Funcs/Comparisons.R")

# Specific to each method
source("MSFAST/MSFAST_Help.R")
source("Univar_Comparators/GFSR/GFSR_Help.R")
source("Univar_Comparators/Conditional/Cond_Help.R")
source("Univar_Comparators/FACE_Help.R")
source("Univar_Comparators/PACE_Help.R")
source("Univar_Comparators/VMP_Help.R")

results_directory = "Results"

create_dir <- function(direct){
  if (!file.exists(direct)){
    dir.create(direct, recursive = T)
  }
}

args = commandArgs(trailingOnly=TRUE)

n_sim = as.numeric(args[3])
obs_ran = as.numeric(args[4]):as.numeric(args[5])
snr = as.numeric(args[6])
N = 100
M = 100
Q = 20
K = 3
P = 1
sim_domain = seq(0, 1, length.out = M)
quad_weights = booleQuad(sim_domain)

func_objs = Univar_F(snr)
MuFun = func_objs$Mu
FPCs = func_objs$FPC
EVals = func_objs$EV
sigma2 = func_objs$Sig2

# Create output directories
{
  out_dir = args[1]
  
  FE_dir = paste0(results_directory, "/FE/", out_dir)
  EF_dir = paste0(results_directory, "/EF/", out_dir)
  Smooth_dir = paste0(results_directory, "/Smooth/", out_dir)
  
  create_dir(FE_dir)
  create_dir(EF_dir)
  create_dir(Smooth_dir)
}

message("Output directory: ", args[1])
message("Output file: ", args[2])
message("Number of simulations: ", args[3])
message("Range: " ,args[4], " - ", args[5])
message("SNR: ", args[6])

FE_out = data.frame(Method = c(), Cov = c(), MSE = c(), Sample = c())
EF_out = data.frame(Method = c(), FPC_Num = c(), Cov = c(), MSE = c(), 
                    Sample = c())
Smooth_out = data.frame(Method = c(), RISE = c(), Cov = c(), Sample = c())

for(x in 1:n_sim){
  print(paste0("Iteration: ", x))
  
  # Generate the data
  sim_dataset = gen_mFPCA(FPCs, N, P, sim_domain, MuFun, EVals, sigma2, obs_ran, F)
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
    data_list = MSFAST_datalist(sim_dataset$ddf, N, K, Q, sub_domain, 
                                basis_type = "B", scale = T)
    
    # Fit model
    mod = stan(
      file = "MSFAST/MSFAST.stan",
      data = data_list, 
      chains = 4, 
      cores = 4, 
      warmup = 1000, 
      iter = 1500, 
      control = list(max_treedepth = 12, 
                     adapt_delta = 0.85),
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
    
    # Calculate performance
    FE_FAST =  FE_comp(mObjs$Mu, Mu_df, "MSFAST") %>%
      mutate(Sample = x)
    EF_FAST = EF_comp(mObjs$FPC, EF_df, "MSFAST") %>%
      mutate(Sample = x)
    Smooth_FAST = Smooth_comp(trueSmooth, Smooth_df, "MSFAST") %>%
      mutate(Sample = x)
  }
  
  # GFSR (Gertheis and Goldsmith 2018)
  {
    # Collate data
    data_list = GFSR_datalist(sim_dataset$ddf$Y, sim_dataset$ddf$Subj, 
                              sim_dataset$ddf$Arg, N, K, Q, sim_domain)
    
    # Unable to install package in 2024, had to take STAN code directly
    mod = stan(
      file = "Univar_Comparators/GFSR/GFSR.stan",
      data = data_list,
      chains = 4, 
      cores = 4, 
      warmup = 1000, 
      iter = 1500, 
      control = list(max_treedepth = 12,
                     adapt_delta = 0.85),
      verbose = F,
      refresh = 0
    )
    
    # Extract
    new_B = bSpline(sim_domain, df = Q, intercept = T)
    objects = GFSR_extract(mod, new_B, M, K)
    align = procrust_FPC(objects$EF, mObjs$mFPC)
    
    Mu_df = FE_Summary(objects$Mu, sim_domain, P)
    EF_ests = FPC_Est_RSV(align, M, P, sim_domain, mObjs$mFPC)
    EF_df = inner_join(FPC_CI(align, sim_domain, P), 
                       EF_ests)
    Smooth_df = Smooth_Summary(N, P, objects$Mu, objects$EF, 
                                objects$Score, sim_domain)
    
    # Calculate performance
    FE_GFSR = FE_comp(mObjs$Mu, Mu_df, "GFSR") %>%
      mutate(Sample = x)
    EF_GFSR = EF_comp(mObjs$FPC, EF_df, "GFSR") %>%
      mutate(Sample = x)
    Smooth_GFSR = Smooth_comp(trueSmooth, Smooth_df, "GFSR") %>%
      mutate(Sample = x)
  }
  
  # FACE (Fast Covariance Smoothing - Xiao 2018)
  {
    # Collate data
    fpca_df = FACE_datalist(sim_dataset$ddf)
    
    # Find pve for 3 components with binary search
    {
      flag = FALSE
      lend = 0
      rend = 1
      var_exp = 0.5
      while(flag == FALSE){
        test_mod = face.sparse(data = fpca_df, argvals.new = sim_domain,
                               knots = Q-3, pve = var_exp)
        if(test_mod$npc == K){
          flag = TRUE
        }
        else if(test_mod$npc > K){
          rend = var_exp
          var_exp = (var_exp + lend)/2
        }
        else{
          lend = var_exp
          var_exp = (var_exp + rend)/2
        }
      } 
    }
    
    # Fit the model
    mod = face.sparse(data = fpca_df, argvals.new = sim_domain, 
                      knots = Q-3, calculate.scores = T, pve = var_exp)
    
    # Extract
    FACE_objs = FACE_extract(mod, sim_domain, mObjs$mFPC, N, M, K, fpca_df)
    
    # Calculate performance
    FE_FACE = FE_comp(mObjs$Mu, FACE_objs$Mu %>% mutate(Var = 1), "FACE") %>%
      mutate(Sample = x)
    EF_FACE = EF_comp(mObjs$FPC, FACE_objs$EF %>% mutate(Var = 1), "FACE") %>%
      mutate(Sample = x)
    Smooth_FACE = Smooth_comp(trueSmooth, FACE_objs$Smooth %>% mutate(Var = 1),
                             "FACE") %>%
      mutate(Sample = x)
  }
  
  # Conditional on eigenfunction estimates from FACE
  {
    # Collate data
    data_list = cond_datalist(N, K, sim_dataset$ddf$Y, sim_dataset$ddf$Subj, 
                              sim_dataset$ddf$S, FACE_objs$Muv, FACE_objs$EFm)
    
    # Fit conditional model
    mod = stan(
      file = "Univar_Comparators/Conditional/Cond.stan",
      data = data_list, 
      chains = 4, 
      cores = 4, 
      warmup = 1000, 
      iter = 1500,
      verbose = F,
      refresh = 0
    )
    
    # Extract
    objects = cond_extract(mod, N, K, FACE_objs$EFm, sim_domain)
    
    Mu_df = data.frame(Arg = sim_domain, Est = FACE_objs$Muv, LB = NA, UB = NA)
    EF_df = FPC_df(FACE_objs$EFm, sim_domain) %>%
      mutate(LB = NA, UB = NA)
    Smooth_df = objects$Smooth
    
    # Calculate performance
    FE_COND = FE_comp(mObjs$Mu, Mu_df %>% mutate(Var = 1), "COND") %>%
      mutate(Sample = x)
    EF_COND = EF_comp(mObjs$FPC, EF_df %>% mutate(Var = 1), "COND") %>%
      mutate(Sample = x)
    Smooth_COND = Smooth_comp(trueSmooth, Smooth_df %>% mutate(Var = 1),
                             "COND") %>%
      mutate(Sample = x)
  }
  
  # fdapace (Wang and Mueller last updated 2024)
  {
    # Collate data
    pace_data = PACE_datalist(sim_dataset$ddf$Y, sim_dataset$ddf$Subj, 
                              sim_dataset$ddf$S, N, sub_domain)
    
    # Fit the model
    mod = FPCA(Ly = pace_data$list_y, Lt = pace_data$list_t, 
               optns = list(dataType = 'Sparse', maxK = K, 
                            nRegGrid = length(sim_domain)))
    
    # Extract
    objects = PACE_extract(mod, K, mObjs$mFPC, sim_domain)
    
    # Calculate performance
    FE_PACE = FE_comp(mObjs$Mu, objects$Mu %>% mutate(Var = 1), "PACE") %>%
      mutate(Sample = x)
    EF_PACE = EF_comp(mObjs$FPC, objects$EF %>% mutate(Var = 1), "PACE") %>%
      mutate(Sample = x)
    Smooth_PACE = Smooth_comp(trueSmooth, objects$Smooth %>% mutate(Var = 1),
                             "PACE") %>%
      mutate(Sample = x)
  }
  
  # Variational message passing (Nolan 2025)
  {
    # Collate data
    vmp_data = VMP_datalist(sim_dataset$ddf$Y, sim_dataset$ddf$Subj, 
                             sim_dataset$ddf$S, N, sub_domain)
    
    # Fit the model
    mod = run_vmp_fpca(time_obs = vmp_data$list_t, Y = vmp_data$list_y,
                       time_g = sim_domain, L = K, verbose = F)
    
    # Extract
    objects = VMP_extract(mod, K, mObjs$mFPC, sim_domain, N)
    
    # Calculate performance
    FE_VMP = FE_comp(mObjs$Mu, objects$Mu %>% mutate(Var = 1), "VMP") %>%
      mutate(Sample = x)
    EF_VMP = EF_comp(mObjs$FPC, objects$EF %>% mutate(Var = 1), "VMP") %>%
      mutate(Sample = x)
    Smooth_VMP = Smooth_comp(trueSmooth, objects$Smooth %>% mutate(Var = 1),
                            "VMP") %>%
      mutate(Sample = x)
  }
  
  # Simple mean functions as smooth
  {
    Smooth_df = map(sort(unique(fpca_df$subj)), function(x){
      sub_data = fpca_df %>%
        filter(subj == x) %>%
        pull(y)
      
      return(data.frame(Curve = paste0("Curve ", x),
                        Arg = sim_domain, 
                        Est = mean(sub_data),
                        UB = NA, LB = NA))
    }) %>% list_rbind()
    
    Smooth_mean = Smooth_comp(trueSmooth, Smooth_df %>% mutate(Var = 1),
                             "Mean") %>%
      mutate(Sample = x)
  }
  
  # Collate model component results
  {
    FE_df = rbind(FE_FAST, FE_GFSR, FE_FACE, FE_COND, FE_PACE, FE_VMP)
    EF_df = rbind(EF_FAST, EF_GFSR, EF_FACE, EF_COND, EF_PACE, EF_VMP)
  }
  
  # Collate prediction results (RISE)
  {
    Smooth_df = rbind(Smooth_FAST, Smooth_GFSR, Smooth_FACE, 
                      Smooth_COND, Smooth_PACE, Smooth_VMP)
    Smooth_df = left_join(Smooth_df, Smooth_mean %>% 
                            rename(BaseISE = ISE) %>%
                            select(Var, Curve, BaseISE)) %>%
      mutate(RISE = ISE/BaseISE) %>%
      select(-c(BaseISE, ISE)) %>%
      group_by(Method) %>%
      summarize(RISE = mean(RISE), Cov = mean(Cov), 
                Sample = first(Sample))
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
