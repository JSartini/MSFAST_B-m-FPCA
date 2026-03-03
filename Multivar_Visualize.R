source("Gen_Funcs/Libs.R")
library(ggridges)
library(ggforce)
library(scales)
library(patchwork)

# Helper function for reading in simulation results
read_directory <- function(directory){
  files = list.files(directory, full.names = T)
  max_sample = 0
  all_results = map(files, function(x){
    ind_data = read.csv(x)
    ind_data$Sample = ind_data$Sample + max_sample
    max_sample <<- max(ind_data$Sample)
    return(ind_data)
  }) %>% list_rbind()
  return(all_results)
}

# Set colors for the methods
levord = c("MSFAST", "uFPCA", "VMP", "MF", "mFACEs", "mFPCA")
colors = pal_hue()(6)
names(colors) = levord

#--------------------------------------------#
# Timing Simulations                         #
#--------------------------------------------#

# Visualization by subjects
N_vals = c(100, 200, 300, 400, 500, 1000)
P_vals = c(3)

Timing_N = data.frame()
for(N in N_vals){
  for(P in P_vals){
    sub_name = paste0("Multi_Time_N", N, "_P", P)
    file_name = paste0("output_", P, "_", N, ".csv")
    new_df = read.csv(paste0("Timing_Results/Time/", sub_name, "/", file_name)) %>%
      select(-c(X)) %>%
      mutate(N = N, P = P)
    Timing_N = rbind(Timing_N, new_df)
  }
}

N_vals = c(100)
P_vals = c(3, 4, 5, 6, 7, 8)
Timing_P = data.frame()
for(N in N_vals){
  for(P in P_vals){
    sub_name = paste0("Multi_Time_N", N, "_P", P)
    file_name = paste0("output_", P, "_", N, ".csv")
    new_df = read.csv(paste0("Timing_Results/Time/", sub_name, "/", file_name)) %>%
      select(-c(X)) %>%
      mutate(N = N, P = P)
    Timing_P = rbind(Timing_P, new_df)
  }
}

p1 = Timing_N %>%
  filter(Method != "MSFAST Parallel - Variates") %>%
  mutate(Time = Timing/60) %>%
  mutate(Method = factor(Method, levels = levord)) %>%
  group_by(N, Method) %>%
  summarize(med_val = median(Time, na.rm = T), 
            min_val = min(Time, na.rm = T), 
            max_val = max(Time, na.rm = T)) %>%
  ggplot(aes(x = N, group = Method, color = Method)) + 
  geom_line(aes(y = med_val)) +
  geom_errorbar(aes(ymin = min_val, ymax = max_val), width = 50) + 
  theme_bw() + 
  scale_color_manual(values = colors) +
  labs(x = "Number of Participants I", y = "Computation Time (min)")
p1

p2 = Timing_P %>%
  filter(Method != "MSFAST Parallel - Variates") %>%
  mutate(Time = Timing/60) %>%
  mutate(Method = factor(Method, levels = levord)) %>%
  group_by(P, Method) %>%
  summarize(med_val = median(Time, na.rm = T), 
            min_val = min(Time, na.rm = T), 
            max_val = max(Time, na.rm = T)) %>%
  ggplot(aes(x = P, group = Method, color = Method)) + 
  geom_line(aes(y = med_val)) +
  geom_errorbar(aes(ymin = min_val, ymax = max_val), width = 0.2) + 
  theme_bw() + 
  scale_color_manual(values = colors) +
  labs(x = "Number of Variates P", y = "Computation Time (min)")
p2

(p1 | p2) + plot_layout(guides = "collect", axes = "collect")
ggsave("Figs/Timing_MV.png", height = 5, width = 9)

#--------------------------------------------#
# Parallel timing                            #
#--------------------------------------------#

read.csv("Timing_Results/Parallel.csv") %>%
  group_by(N, Method) %>%
  summarize(atime = median(Times/60),
            ltime = min(Times/60),
            utime = max(Times/60)) %>%
  mutate(Method = case_when(Method == "2 Cores" ~ "2 Chains", 
                            Method == "4 Cores" ~ "4 Chains", 
                            Method == "Threading" ~ "2 Chains\n4 Threads"), 
         Method = factor(Method, levels = c("2 Chains", "4 Chains", 
                                            "2 Chains\n4 Threads")))%>%
  ggplot(aes(x = N, group = Method, color = Method)) +
  geom_pointrange(aes(y = atime, ymin = ltime, ymax = utime),
                  position = position_dodge(width = 10)) +
  geom_line(aes(y = atime, group = Method),
            position = position_dodge(width = 10)) +
  theme_bw() +
  labs(x = "Number of Subjects I", y = "Compute Time (min)")
ggsave("Figs/Parallel_Compute.png", width = 9, height = 5)

read.csv("Timing_Results/Parallel.csv") %>%
  group_by(N, Method) %>%
  summarize(atime = median(Times/60)) %>%
  pivot_wider(names_from = Method, values_from = atime) %>%
  mutate(rel_thread = Threading/`2 Cores`,
         rel_4 = `4 Cores`/`2 Cores`)

#--------------------------------------------#
# Properly specified simulations             #
#--------------------------------------------#

# Collate cluster performance data
direct = "Results"

Smooths = read_directory(paste0(direct, "/Smooth/Multi")) %>%
  select(-c(X)) %>%
  mutate(Data = case_when(NObs == "3-7" ~ "~5 Obs.",
                          TRUE ~ "~10 Obs."))

EF = read_directory(paste0(direct, "/EF/Multi")) %>%
  select(-c(X)) %>%
  mutate(Data = case_when(NObs == "3-7" ~ "~5 Obs.",
                          TRUE ~ "~10 Obs."))

FE = read_directory(paste0(direct, "/FE/Multi")) %>%
  select(-c(X)) %>%
  mutate(Data = case_when(NObs == "3-7" ~ "~5 Obs.",
                          TRUE ~ "~10 Obs."))

num_components = 3
snr = 4

# Accuracy of Smooth components
p1 = Smooths %>%
  filter(Method != "MSFAST_vb") %>%
  filter(P == num_components & SNR == snr) %>%
  mutate(Data = factor(Data, levels = c("~5 Obs.", "~10 Obs.")),
         Var = factor(paste0("Var ", Var), levels = c("Var 1", "Var 2", "Var 3", "Var 4")),
  Method = factor(Method, levels = levord)) %>%
  filter(RISE < 1.5) %>%
  ggplot() + 
  geom_boxplot(aes(x = Method, y = RISE, fill = Method), outliers = F) + 
  theme_bw() + 
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  facet_grid(Data~Var, scales = "free_y") + 
  scale_fill_manual(values = colors) +
  scale_y_continuous(trans = "log10") + 
  labs(y = "Relative ISE")
p1

# Point-wise coverage of Smooth components by 95% credible intervals
p2 = Smooths %>%
  filter(Method != "MSFAST_vb") %>%
  filter(P == num_components & SNR == snr) %>%
  mutate(Data = factor(Data, levels = c("~5 Obs.", "~10 Obs.")),
         Var = factor(paste0("Var ", Var), levels = c("Var 1", "Var 2", "Var 3", "Var 4")),
         Method = factor(Method, levels = levord)) %>% 
  rename(COV = Cov) %>%
  ggplot() + 
  geom_density_ridges(aes(y = Method, fill = Method, x = COV, color = Method), alpha = 0.2, 
                      quantile_lines=TRUE, quantile_fun=function(Mean_Cov,...)mean(Mean_Cov)) + 
  geom_vline(xintercept = 0.95, linetype = "dashed") + 
  scale_x_continuous(limits = c(0.5, 1)) + 
  theme_bw() + 
  coord_flip() + # xlim = c(0.5, 1) 
  theme(legend.position = "none",
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(x = "Average Coverage") + 
  scale_fill_manual(values = colors) +
  facet_grid(Data~Var, scales = "free_y")
p2

# Combine RISE and Coverage into a single, two-panel figure
(p1 | p2) + plot_annotation(tag_levels = 'A', tag_suffix = ")")
ggsave("Figs/Smooths_all_MV.png", width = 12, height = 6)

# Functional model component summaries
func_df = FE %>%
  filter(P == num_components & SNR == snr) %>%
  mutate(EFNum = "mu(t)") %>%
  rbind(EF %>% 
          filter(P == num_components) %>%
          mutate(EFNum = paste0("phi[", substring(FPC_Num, 4), "](t)")) %>%
          select(-c(FPC_Num))) %>%
  mutate(Data = factor(Data, levels = c("~5 Obs.", "~10 Obs.")),
         Method = factor(Method, levels = levord),
         Var = factor(paste0("Var ", Var), levels = c("Var 1", "Var 2", "Var 3", "Var 4"))) 

# Accuracy - ISE
func_df %>%
  filter(Method != "uFPCA") %>%
  ggplot() + 
  geom_boxplot(aes(x = Method, y = ISE, fill = Method), outliers = F) + 
  theme_bw() + 
  facet_grid(EFNum ~ Data + Var, labeller = labeller(EFNum = label_parsed),
             scales = "free_y") +
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values = colors) +
  labs(y = "Eigenfunction ISE")
ggsave("Figs/Functional_ISE.png", width = 9, height = 8)

# Inference - Point-wise coverage of 95% credible intervals
func_df %>%
  filter(Method != "uFPCA") %>%
  rename(COV = Cov) %>%
  drop_na() %>%
  ggplot() + 
  geom_density_ridges(aes(y = Method, fill = Method, x = COV, color = Method), alpha = 0.2, 
                      quantile_lines=TRUE, quantile_fun=function(Mean_Cov,...)mean(Mean_Cov)) + 
  geom_vline(xintercept = 0.95, linetype = "dashed") + 
  theme_bw() + 
  coord_flip(xlim = c(0.5, 1), expand = F) +
  theme(legend.position = "none",
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(x = "Average Coverage") + 
  scale_fill_manual(values = colors) + 
  facet_grid(EFNum ~ Data + Var, labeller = labeller(EFNum = label_parsed),
             scales = "free_y")
ggsave("Figs/Functional_COV.png", width = 9, height = 8)




