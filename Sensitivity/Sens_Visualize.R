setwd("Sensitivity")

library(ggplot2)
library(ggridges)
library(ggforce)
library(scales)
library(tidyverse)
library(ggpubr)

#---------------------------------#
# Visualize Q simulations         #
#---------------------------------#

Q_vals = c(5, 10, 20, 30, 40)

EF_Q = data.frame()
for(Q in Q_vals){
  sub_name = paste0("Results/EF/MS_Q", Q, "/")
  file_name = list.files(sub_name)
  new_df = read.csv(paste0(sub_name, file_name[1])) %>%
    select(-c(X)) 
  EF_Q = rbind(EF_Q, new_df)
}

EF_Q %>%
  mutate(EFNum = paste0("phi[", substring(FPC_Num, 5), "](t)"),
         Var = paste0("Var ", Var)) %>%
  select(-c(FPC_Num)) %>%
  ggplot() + 
  geom_boxplot(aes(x = Q, y = ISE, group = Q), outliers = F) + 
  theme_bw() + 
  facet_grid(EFNum ~ Var, labeller = labeller(EFNum = label_parsed)) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(y = "Eigenfunction ISE", x = "Spline Basis Dimension Q")
ggsave("Figures/Q_EF_ISE.png", height = 8, width = 9)

EF_Q %>%
  mutate(EFNum = paste0("phi[", substring(FPC_Num, 5), "](t)"),
         Var = paste0("Var ", Var)) %>%
  select(-c(FPC_Num)) %>% 
  drop_na() %>%
  rename(COV = Cov) %>%
  ggplot() + 
  geom_density_ridges(aes(y = Q, group = Q, x = COV), alpha = 0.2, 
                      quantile_lines=TRUE, quantile_fun=function(Mean_Cov,...)mean(Mean_Cov)) + 
  geom_vline(xintercept = 0.95, linetype = "dashed") + 
  theme_bw() + 
  coord_flip(xlim = c(0.7, 1)) + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(x = "Eigenfunction Coverage", y = "Spline Basis Dimension Q") + 
  facet_grid(EFNum ~ Var, labeller = labeller(EFNum = label_parsed),
             scales = "free_y")
ggsave("Figures/Q_EF_Cov.png", height = 8, width = 9)

Smooth_Q = data.frame()
for(Q in Q_vals){
  sub_name = paste0("Results/Smooth/MS_Q", Q, "/")
  file_name = list.files(sub_name)
  new_df = read.csv(paste0(sub_name, file_name[1])) %>%
    select(-c(X)) 
  Smooth_Q = rbind(Smooth_Q, new_df)
}

Smooth_Q %>%
  mutate(Var = paste0("Var ", Var)) %>%
  ggplot() + 
  geom_boxplot(aes(x = Q, y = RISE, group = Q), outliers = F) + 
  theme_bw() + 
  facet_grid(. ~ Var, labeller = labeller(EFNum = label_parsed)) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(y = "Smooth Function ISE", x = "Spline Basis Dimension Q")
ggsave("Figures/Q_Smooth_ISE.png", height = 5, width = 9)

Smooth_Q %>%
  drop_na() %>%
  rename(COV = Cov) %>%
  mutate(Var = paste0("Var ", Var)) %>%
  ggplot() + 
  geom_density_ridges(aes(y = Q, group = Q, x = COV), alpha = 0.2, 
                      quantile_lines=TRUE, quantile_fun=function(Mean_Cov,...)mean(Mean_Cov)) + 
  geom_vline(xintercept = 0.95, linetype = "dashed") + 
  theme_bw() + 
  coord_flip(xlim = c(0.7, 1)) + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(x = "Smooth Function Coverage", y = "Spline Basis Dimension Q") + 
  facet_grid(. ~ Var, labeller = labeller(EFNum = label_parsed))
ggsave("Figures/Q_Smooth_Cov.png", height = 5, width = 9)

#---------------------------------#
# Visualize K simulations         #
#---------------------------------#

K_vals = c(2, 3, 4, 5, 6)

EF_K = data.frame()
for(K in K_vals){
  sub_name = paste0("Results/EF/MS_K", K, "/")
  file_name = list.files(sub_name)
  new_df = read.csv(paste0(sub_name, file_name[1])) %>%
    select(-c(X)) 
  EF_K = rbind(EF_K, new_df)
}

EF_K %>%
  mutate(EFNum = paste0("phi[", substring(FPC_Num, 5), "](t)"),
         Var = paste0("Var ", Var)) %>%
  select(-c(FPC_Num)) %>%
  ggplot() + 
  geom_boxplot(aes(x = K, y = ISE, group = K), outliers = F) + 
  theme_bw() + 
  facet_grid(EFNum ~ Var, labeller = labeller(EFNum = label_parsed)) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(y = "Eigenfunction ISE", x = "FPC Basis Dimension K")
ggsave("Figures/K_EF_ISE.png", height = 8, width = 9)

EF_K %>%
  mutate(EFNum = paste0("phi[", substring(FPC_Num, 5), "](t)"),
         Var = paste0("Var ", Var)) %>%
  select(-c(FPC_Num)) %>% 
  drop_na() %>%
  rename(COV = Cov) %>%
  ggplot() + 
  geom_density_ridges(aes(y = K, group = K, x = COV), alpha = 0.2, 
                      quantile_lines=TRUE, quantile_fun=function(Mean_Cov,...)mean(Mean_Cov)) + 
  geom_vline(xintercept = 0.95, linetype = "dashed") + 
  theme_bw() + 
  coord_flip(xlim = c(0.7, 1)) + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(x = "Eigenfunction Coverage", y = "FPC Basis Dimension K") +
  facet_grid(EFNum ~ Var, labeller = labeller(EFNum = label_parsed),
             scales = "free_y")
ggsave("Figures/K_EF_Cov.png", height = 8, width = 9)

Smooth_K = data.frame()
for(K in K_vals){
  sub_name = paste0("Results/Smooth/MS_K", K, "/")
  file_name = list.files(sub_name)
  new_df = read.csv(paste0(sub_name, file_name[1])) %>%
    select(-c(X)) 
  Smooth_K = rbind(Smooth_K, new_df)
}

Smooth_K %>%
  mutate(Var = paste0("Var ", Var)) %>%
  ggplot() + 
  geom_boxplot(aes(x = K, y = RISE, group = K), outliers = F) + 
  theme_bw() + 
  facet_grid(. ~ Var, labeller = labeller(EFNum = label_parsed)) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(y = "Smooth Function ISE",  x = "FPC Basis Dimension K")
ggsave("Figures/K_Smooth_ISE.png", height = 5, width = 9)

Smooth_K %>%
  drop_na() %>%
  rename(COV = Cov) %>%
  mutate(Var = paste0("Var ", Var)) %>%
  ggplot() + 
  geom_density_ridges(aes(y = K, group = K, x = COV), alpha = 0.2, 
                      quantile_lines=TRUE, quantile_fun=function(Mean_Cov,...)mean(Mean_Cov)) + 
  geom_vline(xintercept = 0.95, linetype = "dashed") + 
  theme_bw() + 
  coord_flip(xlim = c(0.8, 1)) + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(x = "Smooth Function Coverage", y = "FPC Basis Dimension K") + 
  facet_grid(. ~ Var, labeller = labeller(EFNum = label_parsed))
ggsave("Figures/K_Smooth_Cov.png", height = 5, width = 9)
  
#---------------------------------#
# Visualize alpha simulations     #
#---------------------------------#

A_vals = c(0.01, 0.05, 0.1, 0.2, 0.3)

EF_A = data.frame()
for(A in A_vals){
  sub_name = paste0("Results/EF/MS_Alpha", A, "/")
  file_name = list.files(sub_name)
  new_df = read.csv(paste0(sub_name, file_name[1])) %>%
    select(-c(X)) 
  EF_A = rbind(EF_A, new_df)
}

EF_A %>%
  mutate(EFNum = paste0("phi[", substring(FPC_Num, 5), "](t)"),
         Var = paste0("Var ", Var)) %>%
  select(-c(FPC_Num)) %>%
  ggplot() + 
  geom_boxplot(aes(x = Alpha, y = ISE, group = Alpha), outliers = F, 
               width = 0.03) + 
  theme_bw() + 
  facet_grid(EFNum ~ Var, labeller = labeller(EFNum = label_parsed)) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(y = "Eigenfunction ISE", x = "Penalty Hyper-parameter Alpha") + 
  scale_fill_manual(values = colors)
ggsave("Figures/Alpha_EF_ISE.png", height = 8, width = 9)

EF_A %>%
  mutate(EFNum = paste0("phi[", substring(FPC_Num, 5), "](t)"),
         Var = paste0("Var ", Var)) %>%
  select(-c(FPC_Num)) %>% 
  drop_na() %>%
  rename(COV = Cov) %>%
  ggplot() + 
  geom_density_ridges(aes(y = Alpha, group = Alpha, x = COV), alpha = 0.2, 
                      quantile_lines=TRUE, quantile_fun=function(Mean_Cov,...)mean(Mean_Cov)) + 
  geom_vline(xintercept = 0.95, linetype = "dashed") + 
  theme_bw() + 
  coord_flip(xlim = c(0.7, 1)) + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(x = "Eigenfunction Coverage", y = "Penalty Hyper-parameter Alpha") + 
  scale_fill_manual(values = colors) + 
  facet_grid(EFNum ~ Var, labeller = labeller(EFNum = label_parsed),
             scales = "free_y")
ggsave("Figures/Alpha_EF_Cov.png", height = 8, width = 9)

Smooth_A = data.frame()
for(A in A_vals){
  sub_name = paste0("Results/Smooth/MS_Alpha", A, "/")
  file_name = list.files(sub_name)
  new_df = read.csv(paste0(sub_name, file_name[1])) %>%
    select(-c(X)) 
  Smooth_A = rbind(Smooth_A, new_df)
}

Smooth_A %>%
  mutate(Var = paste0("Var ", Var)) %>%
  ggplot() + 
  geom_boxplot(aes(x = Alpha, y = RISE, group = Alpha), 
               outliers = F, width = 0.03) + 
  theme_bw() + 
  facet_grid(. ~ Var, labeller = labeller(EFNum = label_parsed)) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(y = "Smooth Function ISE",  x = "Penalty Hyper-parameter Alpha") + 
  scale_fill_manual(values = colors)
ggsave("Figures/Alpha_Smooth_ISE.png", height = 5, width = 9)

Smooth_A %>%
  drop_na() %>%
  rename(COV = Cov) %>%
  mutate(Var = paste0("Var ", Var)) %>%
  ggplot() + 
  geom_density_ridges(aes(y = Alpha, group = Alpha, x = COV), alpha = 0.2, 
                      quantile_lines=TRUE, quantile_fun=function(Mean_Cov,...)mean(Mean_Cov)) + 
  geom_vline(xintercept = 0.95, linetype = "dashed") + 
  theme_bw() + 
  coord_flip(xlim = c(0.8, 1)) + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(x = "Smooth Function Coverage", y = "Penalty Hyper-parameter Alpha") + 
  scale_fill_manual(values = colors) + 
  facet_grid(. ~ Var, labeller = labeller(EFNum = label_parsed))
ggsave("Figures/Alpha_Smooth_Cov.png", height = 5, width = 9)

#---------------------------------#
# Timing table                    #
#---------------------------------#

Q_vals = c(5, 10, 20, 30, 40)
K_vals = c(3, 4, 5, 6)

timing_sens = data.frame()
for(Q in Q_vals){
  for(K in K_vals){
    sub_name = paste0("Results/TIme/MS_Q", Q, "_K", K,  "/")
    file_name = list.files(sub_name)
    new_df = read.csv(paste0(sub_name, file_name[1])) %>%
      select(-c(X)) 
    timing_sens = rbind(timing_sens, new_df)
  }
}

timing_sens %>%
  group_by(K, Q) %>%
  summarize(med_time = median(Seconds)) %>%
  pivot_wider(names_from = Q, values_from = med_time)
