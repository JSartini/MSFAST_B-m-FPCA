# Function to evaluate FPC bases
eval_FPCs <- function(func_list, domain){
  M = length(domain)
  K = length(func_list)
  out = matrix(0, nrow = M, ncol = K)
  for(k in 1:K){
    out[,k] = func_list[[k]](domain)
  }
  return(out)
}

# Function to generate random GP deviations
gen_dev <- function(FPC_func, domain, Lambdas){
  P = length(FPC_func)
  K = length(FPC_func[[1]])
  
  Xi = rep(0, K)
  for(k in 1:K){
    Xi[k] = rnorm(1, sd = sqrt(Lambdas[k]))
  }
  output = map(1:P, function(p){
    Phi = eval_FPCs(FPC_func[[p]], domain[[p]])
    return(list(Delta = Phi %*% Xi, Scores = Xi, Points = domain[[p]]))
  })
  return(output)
}

# Generative based upon KKL decomposition
gen_mFPCA <- function(FPCs, N, P, domain, mu_func, Lambdas, sigma_eps, obs_ran, 
                      diff = F){
  M = length(domain)
  
  # Latent smooths
  GP_devs = map(1:N, function(n){
    if(diff){
      n_obs = sample(obs_ran, size = P)
      obs_points = map(1:P, function(p){
        pts = sample(domain, n_obs[p])
        return(pts)
      })
    }
    else{
      n_obs = sample(obs_ran, 1)
      pts = sample(domain, n_obs)
      obs_points = map(1:P, function(p){
        return(pts)
      })
    }
    
    latent = gen_dev(FPCs, obs_points, Lambdas)
    return(latent)
  })
  
  data_df = map(1:N, function(n){
    df_vals = map(1:P, function(p){
      latent = GP_devs[[n]][[p]]
      denoised = latent$Delta + mu_func[[p]](latent$Points)
      y_vals = denoised + rnorm(length(denoised), sd = sigma_eps[p])
      df_vals = data.frame(Y = y_vals, Arg = latent$Points, Subj = n, 
                           Var = p)
    }) %>% list_rbind()
    return(df_vals)
  }) %>% 
    list_rbind() %>%
    arrange(Arg) %>%
    group_by(Arg) %>%
    mutate(S = cur_group_id()) %>%
    ungroup()
  
  scores = do.call(rbind, map(GP_devs, function(latent){
    return(latent[[1]]$Scores) # Scores are identical over outcomes
  }))
  
  overall_phi = do.call(rbind, map(1:P, function(p){
    return(eval_FPCs(FPCs[[p]], domain))
  }))
  smooth = scores %*% t(overall_phi)
  
  return(list(ddf = data_df, Y_true = smooth, Score_true = scores, 
              new_domain = sort(unique(data_df$Arg))))
}

# Function for extracting Truth at various granularities
true_funcs <- function(MuFun, FPCs, P, weights, domain){
  mfpc = do.call(rbind, map(1:P, function(p){
    return(eval_FPCs(FPCs[[p]], domain))
  })) 
  tEFM = data.frame(mfpc)
  colnames(tEFM) = paste0("FPC ", 1:ncol(tEFM))
  tEFM$Arg = rep(domain, P)
  tEFM$Var = rep(1:P, each = length(domain))
  tEFM$wei = rep(weights, P)
  tEFM = tEFM %>%
    pivot_longer(-c(Arg, Var, wei), names_to = "FPC_Num", values_to = "Func")
  
  tMUM = data.frame(Arg = rep(domain, P),
                    Func = map(1:P, function(p){
                      return(MuFun[[p]](domain))
                      }) %>% unlist(),
                    Var = rep(1:P, each = length(domain)),
                    wei = rep(weights, P))
  
  return(list(FPC = tEFM, Mu = tMUM, mFPC = mfpc))
}

# Simulation adapted from Xiao 2018
Univar_F = function(SNR){
  MuFun = list(function(x){
    return(rep(0, length.out = length(x)))
  })
  
  FPCs = list(list(
    FPC1 = function(x) {
      return(sqrt(2) * (sin(2 * pi * x)))
    },
    FPC2 = function(x) {
      return(sqrt(2) * (cos(4 * pi * x)))
    },
    FPC3 = function(x) {
      return(sqrt(2) * (sin(4 * pi * x)))
    }
  ))
  
  EVals = 2^(0:-2)
  Sigma = sum(EVals)/SNR
  return(list(Mu = MuFun, FPC = FPCs, EV = EVals, Sig2 = Sigma))
}

# Simulation adapted from Nolan 2025
Multivar_F = function(P, K, SNR){
  MuFun = map(1:P, function(p){
    func = function(t){
      return((-1)^p*2*sin((2*pi + p)*t))
    }
    return(func)
  })
  
  FPCs = map(1:P, function(p){
    var_FPC = map(1:K, function(k){
      func = function(t){
        khat = (k-1) %/% 2 + 1
        if(k %% 2 == 1){
          return((-1)^p*sqrt(2/P)*cos(2*khat*pi*t))
        }
        else{
          return((-1)^p*sqrt(2/P)*sin(2*khat*pi*t))
        }
      }
      return(func)
    })
    return(var_FPC)
  })
  
  EVals = 2^(-(0:K-1))
  Sigma = rep(sum(EVals)/SNR, P)
  
  return(list(Mu = MuFun, FPC = FPCs, EV = EVals, Sig2 = Sigma))
}


