##### Working directory
setwd("H:/.shortcut-targets-by-id/1hRMzSqAcE5AcmrsOu5k3nw1rmUY1NHyO/QA&COVID19/VisitingLUH2021/Code/funHDGM_sim/VarDecomposition")

##### Packages
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(readxl)

##### Load simulations
Sims_G <- read_excel("Sims_G.xlsx")
Sims_Vz <- read_excel("Sims_Vz.xlsx")
Sims_Theta <- read_excel("Sims_Theta.xlsx")

##### Study the behavior of G
p1 <- Sims_G %>%
  group_by(g) %>%
  summarise(m = mean(VAR_Z), std = sd(VAR_Z), n = n()) %>%
  ggplot(mapping = aes(x = g, y = m)) +
  geom_ribbon(aes(ymin = m - qnorm(0.975)*std/sqrt(n),
                  ymax = m + qnorm(0.975)*std/sqrt(n)), fill = "grey70") +
  geom_line(aes(y = m)) + 
  labs(x = "g", y = "Simulated random effects variance",
       title = "Random effects variance",
       subtitle = latex2exp::TeX("n = 10 & $V_{z}$ = 1 & $\\theta_{z}$ = 0.50"))
p2 <- Sims_G %>%
  group_by(g) %>%
  summarise(m = mean(sigma2_y_sim), std = sd(sigma2_y_sim), n = n()) %>%
  ggplot(mapping = aes(x = g, y = m)) +
  geom_ribbon(aes(ymin = m - qnorm(0.975)*std/sqrt(n),
                  ymax = m + qnorm(0.975)*std/sqrt(n)), fill = "grey70") +
  geom_line(aes(y = m)) + 
  labs(x = "g", y = "Simulated Y total variance",
       title = "Total variance",
       subtitle = latex2exp::TeX("n = 10 & $V_{z}$ = 1 & $\\theta_{z}$ = 0.50"))
p3 <- Sims_G %>%
  group_by(g) %>%
  summarise(m = mean(sigma2_eps_sim), std = sd(sigma2_eps_sim), n = n()) %>%
  ggplot(mapping = aes(x = g, y = m)) +
  geom_ribbon(aes(ymin = m - qnorm(0.975)*std/sqrt(n),
                  ymax = m + qnorm(0.975)*std/sqrt(n)), fill = "grey70") +
  geom_line(aes(y = m)) + 
  labs(x = "g", y = "Simulated Y total variance",
       title = "Total variance",
       subtitle = latex2exp::TeX("n = 10 & $V_{z}$ = 1 & $\\theta_{z}$ = 0.50"))
p12 <- ggarrange(p1,p2, ncol = 2)
p12 <- annotate_figure(p12, 
                top = text_grob("Increasing temporal autocorrelation (g)",
                                color = "red", face = "bold", size = 14))
ggsave(filename = "SimG.png", plot = p12, device = "png", dpi = 300, width = 10)



##### Study the behavior of Vz
p3 <- Sims_Vz %>%
  group_by(v_z) %>%
  summarise(m = mean(VAR_Z), std = sd(VAR_Z), n = n()) %>%
  ggplot(mapping = aes(x = v_z, y = m)) +
  geom_ribbon(aes(ymin = m - qnorm(0.975)*std/sqrt(n),
                  ymax = m + qnorm(0.975)*std/sqrt(n)), fill = "grey70") +
  geom_line(aes(y = m)) + 
  labs(x = latex2exp::TeX("$V_{z}$"), y = "Simulated random effects variance",
       title = "Random effects variance",
       subtitle = latex2exp::TeX("n = 10 & $g$ = 0.85 & $\\theta_{z}$ = 0.50"))
p4 <- Sims_Vz %>%
  group_by(v_z) %>%
  summarise(m = mean(sigma2_y_sim), std = sd(sigma2_y_sim), n = n()) %>%
  ggplot(mapping = aes(x = v_z, y = m)) +
  geom_ribbon(aes(ymin = m - qnorm(0.975)*std/sqrt(n),
                  ymax = m + qnorm(0.975)*std/sqrt(n)), fill = "grey70") +
  geom_line(aes(y = m)) + 
  labs(x = latex2exp::TeX("$V_{z}$"), y = "Simulated Y total variance",
       title = "Total variance",
       subtitle = latex2exp::TeX("n = 10 & $g$ = 0.85 & $\\theta_{z}$ = 0.50"))
p34 <- ggarrange(p3,p4, ncol = 2)
p34 <- annotate_figure(p34, 
                       top = text_grob(latex2exp::TeX("Increasing marginal variance ($V_{Z}$)"),
                                       color = "red", face = "bold", size = 14))
ggsave(filename = "SimVz.png", plot = p34, device = "png", dpi = 300, width = 10)



##### Study the behavior of Theta_z
p5 <- Sims_Theta %>%
  group_by(theta_z) %>%
  summarise(m = mean(VAR_Z), std = sd(VAR_Z), n = n()) %>%
  ggplot(mapping = aes(x = theta_z, y = m)) +
  geom_ribbon(aes(ymin = m - qnorm(0.975)*std/sqrt(n),
                  ymax = m + qnorm(0.975)*std/sqrt(n)), fill = "grey70") +
  geom_line(aes(y = m)) + 
  labs(x = latex2exp::TeX("$\\theta_{z}$ (km)"), y = "Simulated random effects variance",
       title = "Random effects variance",
       subtitle = latex2exp::TeX("n = 10 & $g$ = 0.85 & $V_{z}$ = 1"))
p6 <- Sims_Theta %>%
  group_by(theta_z) %>%
  summarise(m = mean(sigma2_y_sim), std = sd(sigma2_y_sim), n = n()) %>%
  ggplot(mapping = aes(x = theta_z, y = m)) +
  geom_ribbon(aes(ymin = m - qnorm(0.975)*std/sqrt(n),
                  ymax = m + qnorm(0.975)*std/sqrt(n)), fill = "grey70") +
  geom_line(aes(y = m)) + 
  labs(x = latex2exp::TeX("$\\theta_{z}$ (km)"), y = "Simulated Y total variance",
       title = "Total variance",
       subtitle = latex2exp::TeX("n = 10 & $g$ = 0.85 & $V_{z}$ = 1"))
p56 <- ggarrange(p5,p6, ncol = 2)
p56 <- annotate_figure(p56, 
                       top = text_grob(latex2exp::TeX("Increasing range ($\\theta_{Z}$)"),
                                       color = "red", face = "bold", size = 14))
ggsave(filename = "SimTheta.png", plot = p56, device = "png", dpi = 300, width = 10)


##### Unified plot
p16 <- ggarrange(p12,p34,p56, ncol = 1)
ggsave(filename = "SimVarRE.png", plot = p16, device = "png", dpi = 300, width = 10, height = 10)

