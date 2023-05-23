library(readxl)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(xtable)

setwd("H:/.shortcut-targets-by-id/1hRMzSqAcE5AcmrsOu5k3nw1rmUY1NHyO/QA&COVID19/VisitingLUH2021/Code/funHDGM_sim/Application")


####################################
########## Models summary ##########
####################################
Models_summary_Paolo <- read_excel("Models_summary_auto.xlsx", n_max = 17)

Models_summary_Paolo <- Models_summary_Paolo %>%
  rename(SpatPart = `Spatial partitioning`,
         Basis = `Basis number (p=x10+4)`,
         'DSTEM time' = `Computation time - DSTEM (Minutes)`,
         'VarCov time' = `Computation time - VarCov (Minutes)`,
         'PenLik time' = `Computation time - PenLik (Minutes)`,
         'Total time' = `Computation time - Total (Minutes)`) %>%
  mutate(SpatPart = case_when(SpatPart == 1 ~ "k = 1",
                              SpatPart == 2 ~ "k = 2",
                              SpatPart == 3 ~ "k = 3",
                              SpatPart == 4 ~ "k = 4",
                              SpatPart == 5 ~ "k = 5"))

Models_summary_Paolo_red <- Models_summary_Paolo %>%
  select(VarCov,SpatPart,Basis,`DSTEM time`,`VarCov time`,`PenLik time`,`Total time`,
         `Storage dimension (MB)`,`min RMSE`,`1-SE min RMSE`,`min MAE`,`1-SE min MAE`)
xtable(Models_summary_Paolo_red,
       caption = "Summary of the empirical models.",
       label = "Tab_appl_models_sum")

p1 <- Models_summary_Paolo %>%
  filter(SpatPart %in% c("k = 1","k = 2")) %>%
  select(VarCov,SpatPart,Basis,
         'DSTEM time', 'VarCov time','PenLik time','Total time') %>%
  pivot_longer(cols = 4:7,names_to = "Timing",values_to = "Time") %>%
  ggplot() + 
  geom_line(mapping = aes(x = Basis,
                          y = Time,
                          col = VarCov),
            size = 2) + 
  facet_grid(rows = vars(Timing), cols = vars(SpatPart), scales = "free") + 
  labs(x = "Number of basis (b)",
       y = "Computation time (minutes)",
       title = "Computational cost by partition number and model complexity") + 
  scale_x_discrete(limits=c(5,7,9))
ggsave("ComputCost1.png", plot = p1, width = 7.25, height = 4.25, dpi=1200)


p2 <- Models_summary_Paolo %>%
  filter(Basis == 9,
         VarCov == "Approx") %>%
  select(SpatPart,Basis,
         'DSTEM time', 'PenLik time' , 'VarCov time', 'Total time') %>%
  pivot_longer(cols = 3:6,names_to = "Timing",values_to = "Time") %>%
  mutate(SpatPart = ordered(SpatPart, levels=c("k = 1","k = 2","k = 3","k = 4","k = 5"))) %>%
  ggplot(mapping = aes(x = SpatPart,
                       y = Time,
                       col = Timing)) + 
  geom_point(size = 4) + 
  labs(x = "Number of groups (k)",
       y = "Computation time (minutes)",
       title = "Computational cost by partition number",
       subtitle = "b=9 bases")
ggsave("ComputCost2.png", plot = p2, width = 7.25, height = 4.25, dpi=1200)

p3 <- Models_summary_Paolo %>%
  filter(VarCov %in% c("Approx")) %>%
  select(VarCov,SpatPart,Basis,
         'DSTEM time', 'VarCov time','PenLik time','Total time') %>%
  pivot_longer(cols = 4:7,names_to = "Timing",values_to = "Time") %>%
  ggplot() + 
  geom_point(mapping = aes(x = Basis,
                           y = Time,
                           col = SpatPart),
             size = 4) + 
  facet_wrap( ~ Timing, scales = "free") + 
  labs(x = "Number of basis (b)",
       y = "Computation time (minutes)",
       title = "Computational cost by partition number and model complexity") + 
  scale_x_discrete(limits=c(5,7,9))
ggsave("ComputCost3.png", plot = p3, width = 7.25, height = 4.25, dpi=1200)



p4 <- Models_summary_Paolo %>%
  filter(VarCov %in% c("Approx")) %>%
  select(SpatPart,Basis,
         'min RMSE', '1-SE min RMSE', '1-SE min MAE', 'min MAE') %>%
  pivot_longer(cols = 3:6,names_to = "Error",values_to = "Error_val") %>%
  mutate(Error = factor(Error,levels=c('min RMSE','1-SE min RMSE',
                                       'min MAE','1-SE min MAE'))) %>%
  ggplot() + 
  geom_point(mapping = aes(x = Basis,
                           y = Error_val,
                           col = SpatPart),
             size = 4) + 
  facet_wrap(. ~ Error, scales = "free") + 
  labs(x = "Number of basis (b)",
       y = "RMSE or MAE",
       title = "RMSE and MAE by partition number and model complexity",
       subtitle = "minimum and 1-SE rule values") + 
  scale_x_discrete(limits=c(5,7,9))
ggsave("ComputCost4.png", plot = p4, width = 7.25, height = 4.25, dpi=1200)


pcomb <- ggarrange(plotlist = list(p1,p2,p3,p4), ncol = 2, nrow = 2)
ggsave("ComputCosts.png", plot = pcomb, width = 12.25, height = 9.25, dpi=800)



setwd("H:/.shortcut-targets-by-id/1hRMzSqAcE5AcmrsOu5k3nw1rmUY1NHyO/QA&COVID19/VisitingLUH2021/Code/funHDGM_sim/Application")
Coefs <- read_excel("nullcoefs.xlsx")
Coefs_summary <- Coefs %>%
  pivot_longer(cols = 3:8, names_to = "Coef", values_to = "Beta_val") %>%
  group_by(Coef) %>%
  mutate(nCoefs = n()) %>%
  ungroup() %>%
  filter(Beta_val == 0) %>%
  group_by(Coef) %>%
  summarise(NullCoefs = n(), nCoefs = mean(nCoefs), PercNullCoefs = NullCoefs / nCoefs * 100) %>%
  mutate(Basis = case_when(grepl("nb5",Coef) ~ 5,
                           grepl("nb7",Coef) ~ 7,
                           TRUE ~ 9),
         Error = case_when(grepl("1se",Coef) ~ "1-SE min RMSE",
                           TRUE ~ "min RMSE"))

p5 <- Coefs_summary %>%
  mutate(Error = factor(Error,levels=c('min RMSE','1-SE min RMSE'))) %>%
  ggplot(mapping = aes(x = Basis,
                       y = PercNullCoefs,
                       col = Error),) + 
  geom_point(size = 4) + 
  labs(x = "Number of basis (b)",
       y = "% of zero values",
       title = "Percentage of zero coefficients by model complexity",
       subtitle = "Spatial partitioning with k=2 groups") + 
  scale_x_discrete(limits=c(5,7,9)) + 
  theme(legend.position = "bottom")


Coefs_summary2 <- Coefs %>%
  pivot_longer(cols = 3:8, names_to = "Coef", values_to = "Beta_val") %>%
  mutate(Error = case_when(grepl("1se",Coef) ~ "1-SE min RMSE",
                           TRUE ~ "min RMSE")) %>%
  filter(!is.na(Beta_val)) %>%
  group_by(Beta,Error) %>%
  mutate(nCoefs = n()) %>%
  ungroup() %>%
  filter(Beta_val == 0) %>%
  group_by(Beta,Error) %>%
  summarise(NullCoefs = n(), nCoefs = mean(nCoefs), PercNullCoefs = NullCoefs / nCoefs * 100)

beta_names <- c("\u03b2[1]","\u03b2[2]","\u03b2[3]","\u03b2[4]",
                "\u03b2[5]","\u03b2[6]","\u03b2[7]","\u03b2[8]","\u03b2[9]")
p6 <- Coefs_summary2 %>%
  mutate(Error = factor(Error,levels=c('min RMSE','1-SE min RMSE')),
         Beta = factor(Beta)) %>%
  ggplot(mapping = aes(x = Beta,
                       y = PercNullCoefs,
                       col = Error)) + 
  geom_point(size = 4) + 
  labs(x = "Coefficient (basis)",
       y = "% of zero values",
       title = "Percentage of zero coefficients by basis function",
       subtitle = "Spatial partitioning with k=2 groups") + 
  scale_x_discrete(labels = beta_names) + 
  theme(legend.position = "bottom")

pcomb <- ggarrange(plotlist = list(p5,p6), ncol = 2, nrow = 1)
ggsave("ComputCosts2.png", plot = pcomb, width = 12.25, height = 9.25, dpi=800)


