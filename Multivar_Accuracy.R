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
obs_ran = as.numeric(args[4]):as.numeric(args[5])
diff_obs = T
N = 100
M = 100
Q = 20
K = 3
P = 3
sim_domain = seq(0, 1, length.out = M)
quad_weights = booleQuad(sim_domain)

func_objs = Multivar_F(P, K, 4)
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

FE_out = data.frame(Method = c(), Var = c(), Cov = c(), ISE = c(), Sample = c())
EF_out = data.frame(Method = c(), Var = c(), FPC_Num = c(), Cov = c(), 
                    ISE = c(), Sample = c())
Smooth_out = data.frame(Method = c(), Var = c(), RISE = c(), Cov = c(), 
                        Sample = c())

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
    data_list = MSFAST_datalist(sim_dataset$ddf, N, K, Q, sub_domain, 
                                basis_type = "B", scale = T)
    
    # Fit model
    mod = stan(
      file = "MSFAST/MSFAST.stan",
      data = data_list, 
      chains = 4, 
      cores = 4, 
      warmup = 1000, 
      iter = 2000, 
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
  
  # mFACEs (Li and Xiao 2020)
  {
    DL = mfaces_datalist(sim_dataset$ddf, N, P, sim_domain)
    
    face_fit = mface.sparse(DL$fit, argvals.new = sim_domain,
                            knots.option = "quantile", center = T, pve = 0.99) 
    
    FACE_objs = extract_mfaces(face_fit, sim_domain, N, P, mObjs$mFPC, DL$fit, DL$new)
    
    # Calculate performance
    FE_FACE = FE_comp(mObjs$Mu, FACE_objs$FE, "mFACEs") %>%
      mutate(Sample = x)
    EF_FACE = EF_comp(mObjs$FPC, FACE_objs$EF, "mFACEs") %>%
      mutate(Sample = x)
    Smooth_FACE = Smooth_comp(trueSmooth, FACE_objs$Smooth, "mFACEs") %>%
      mutate(Sample = x)
  }
  
  # mFPCA (Happ and Greven 2018)
  {
    DL = mfpca_datalist(sim_dataset$ddf, N, P, sim_domain)
    
    fpca_fit = MFPCA(DL$Data, M = K, DL$Expansions, fit = T)
    FPCA_objs = extract_mfpca(fpca_fit, N, P, sim_domain, mObjs$mFPC, CI = F)
    
    # Calculate performance
    FE_FPCA = FE_comp(mObjs$Mu, FPCA_objs$FE, "mFPCA") %>%
      mutate(Sample = x)
    EF_FPCA = EF_comp(mObjs$FPC, FPCA_objs$EF, "mFPCA") %>%
      mutate(Sample = x)
    Smooth_FPCA = Smooth_comp(trueSmooth, FPCA_objs$Smooth, "mFPCA") %>%
      mutate(Sample = x)
  }
  
  # VMP (T.H. Nolan 2025)
  {
    DL = vmp_datalist(sim_dataset$ddf, N, P)
    
    vmp_fit = run_vmp_fpca(DL$t, DL$y, K, Q, verbose = F, time_g = sim_domain)
    
    VMP_objs = extract_vmp(vmp_fit, sim_domain, N, K, P, mObjs$mFPC)
    
    # Calculate performance
    FE_VMP = FE_comp(mObjs$Mu, VMP_objs$FE, "VMP") %>%
      mutate(Sample = x)
    EF_VMP = EF_comp(mObjs$FPC, VMP_objs$EF, "VMP") %>%
      mutate(Sample = x)
    Smooth_VMP = Smooth_comp(trueSmooth, VMP_objs$Smooth, "VMP") %>%
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
    FE_df = rbind(FE_FAST, FE_FACE, FE_FPCA, FE_VMP)
    EF_df = rbind(EF_FAST, EF_FACE, EF_FPCA, EF_VMP)
  }
  
  # Collate prediction results (RISE)
  {
    Smooth_df = rbind(Smooth_FAST, Smooth_FACE, 
                      Smooth_FPCA, Smooth_VMP)
    Smooth_df = left_join(Smooth_df, Smooth_mean %>% 
                            rename(BaseISE = ISE) %>%
                            select(Var, Curve, BaseISE)) %>%
      mutate(RISE = ISE/BaseISE) %>%
      select(-c(BaseISE, ISE)) %>%
      group_by(Method, Var) %>%
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
