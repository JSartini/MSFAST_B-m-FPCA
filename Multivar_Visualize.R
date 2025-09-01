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
levord = c("MSFAST", "b1", "VMP", "b2", "mFACEs", "mFPCA")
colors = pal_hue()(6)
names(colors) = levord

# Visualization of computation time scaling 
Timing_N = rbind(read.csv("Results/Time/Timing_N100_M2000/output.csv") %>%
                   select(-c(X)) %>%
                   mutate(N = 100),
                 read.csv("Results/Time/Timing_N200_M2000/output.csv") %>%
                   select(-c(X)) %>%
                   mutate(N = 200),
                 read.csv("Results/Time/Timing_N300_M2000/output.csv") %>%
                   select(-c(X)) %>%
                   mutate(N = 300),
                 read.csv("Results/Time/Timing_N400_M2000/output.csv") %>%
                   select(-c(X)) %>%
                   mutate(N = 400),
                 read.csv("Results/Time/Timing_N500_M2000/output.csv") %>%
                   select(-c(X)) %>%
                   mutate(N = 500),
                 read.csv("Results/Time/Timing_N1000_M2000/output.csv") %>%
                   select(-c(X)) %>%
                   mutate(N = 1000))

Timing_N %>%
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
  labs(x = "Number of Participants N", y = "Computation Time (min)")
ggsave("Figs/Timing_MV.png", height = 5, width = 9)

# Collate cluster performance data
Smooths = rbind(read_directory("Results/Smooth/Multi_37") %>%
                  select(-c(X)) %>%
                  mutate(Data = "~5 Obs."),
                read_directory("Results/Smooth/Multi_515") %>% 
                  select(-c(X)) %>%
                  mutate(Data = "~10 Obs."))

EF = rbind(read_directory("Results/EF/Multi_37") %>%
             select(-c(X)) %>%
             mutate(Data = "~5 Obs."),
           read_directory("Results/EF/Multi_515") %>%
             select(-c(X)) %>%
             mutate(Data = "~10 Obs."))

FE = rbind(read_directory("Results/FE/Multi_37") %>%
             select(-c(X)) %>%
             mutate(Data = "~5 Obs."),
           read_directory("Results/FE/Multi_515") %>%
             select(-c(X)) %>%
             mutate(Data = "~10 Obs."))

# Accuracy of Smooth components
p1 = Smooths %>%
  mutate(Data = factor(Data, levels = c("~5 Obs.", "~10 Obs.")),
         Var = factor(paste0("Var ", Var), levels = c("Var 1", "Var 2", "Var 3")),
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
ggsave("Figs/Smooth_RISE.png", width = 9, height = 8)

# Point-wise coverage of Smooth components by 95% credible intervals
p2 = Smooths %>%
  mutate(Data = factor(Data, levels = c("~5 Obs.", "~10 Obs.")),
         Var = factor(paste0("Var ", Var), levels = c("Var 1", "Var 2", "Var 3")),
         Method = factor(Method, levels = levord)) %>%
  rename(COV = Cov) %>%
  ggplot() + 
  geom_density_ridges(aes(y = Method, fill = Method, x = COV, color = Method), alpha = 0.2, 
                      quantile_lines=TRUE, quantile_fun=function(Mean_Cov,...)mean(Mean_Cov)) + 
  geom_vline(xintercept = 0.95, linetype = "dashed") + 
  theme_bw() + 
  coord_flip() + # xlim = c(0.5, 1) 
  theme(legend.position = "none",
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(x = "Average Coverage") + 
  scale_fill_manual(values = colors) + 
  facet_grid(Data~Var, scales = "free_y")
p2
ggsave("Figs/Smooth_COV.png", width = 9, height = 8)

# Combine RISE and Coverage into a single, two-panel figure
(p1 | p2) + plot_annotation(tag_levels = 'A', tag_suffix = ")")
ggsave("Figs/Smooths_all_MV.png", width = 12, height = 6)

# Functional model component summaries
func_df = FE %>%
  mutate(EFNum = "mu(t)") %>%
  rbind(EF %>% 
          mutate(EFNum = paste0("phi[", substring(FPC_Num, 4), "](t)")) %>%
          select(-c(FPC_Num))) %>%
  mutate(Data = factor(Data, levels = c("~5 Obs.", "~10 Obs.")),
         Var = factor(paste0("Var ", Var), levels = c("Var 1", "Var 2", "Var 3")),
         Method = factor(Method, levels = levord)) 

# Accuracy - ISE
func_df %>%
  ggplot() + 
  geom_boxplot(aes(x = Method, y = ISE, fill = Method), outliers = F) + 
  theme_bw() + 
  facet_grid(EFNum ~ Data + Var, labeller = labeller(EFNum = label_parsed),
             scales = "free_y") +
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(y = "Eigenfunction ISE")
ggsave("Figs/Functional_ISE.png", width = 9, height = 8)

# Inference - Point-wise coverage of 95% credible intervals
func_df %>%
  rename(COV = Cov) %>%
  drop_na() %>%
  ggplot() + 
  geom_density_ridges(aes(y = Method, fill = Method, x = COV, color = Method), alpha = 0.2, 
                      quantile_lines=TRUE, quantile_fun=function(Mean_Cov,...)mean(Mean_Cov)) + 
  geom_vline(xintercept = 0.95, linetype = "dashed") + 
  theme_bw() + 
  coord_flip(xlim = c(0.5, 1)) + 
  theme(legend.position = "none",
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(x = "Average Coverage") + 
  scale_fill_manual(values = colors) + 
  facet_grid(EFNum ~ Data + Var, labeller = labeller(EFNum = label_parsed),
             scales = "free_y")
ggsave("Figs/Functional_COV.png", width = 9, height = 8)



