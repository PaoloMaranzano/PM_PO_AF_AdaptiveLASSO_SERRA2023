#########################################
########## Functional boxplots ##########
#########################################

library(readr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

setwd("H:/.shortcut-targets-by-id/1hRMzSqAcE5AcmrsOu5k3nw1rmUY1NHyO/QA&COVID19/VisitingLUH2021/Code/funHDGM_sim/Application")

LongData <- read_csv("LongData.csv")

SummData <- LongData %>%
  group_by(h) %>%
  summarise(across(c(NO2,Pressure,RelHumid,Temperature,Rainfall,WindU,WindV), 
                   list(mean = ~ mean(.x, na.rm=T),
                        median = ~ median(.x, na.rm=T),
                        q1 = ~ quantile(.x, probs = 0.01, na.rm=T),
                        q99 = ~ quantile(.x, probs = 0.99, na.rm=T),
                        q5 = ~ quantile(.x, probs = 0.05, na.rm=T),
                        q95 = ~ quantile(.x, probs = 0.95, na.rm=T),
                        q25 = ~ quantile(.x, probs = 0.25, na.rm=T),
                        q75 = ~ quantile(.x, probs = 0.75, na.rm=T),
                        sd = ~ sd(.x, na.rm=T)),
                   .names = "{.col}_{.fn}")) %>%
  pivot_longer(cols = 2:last_col(),
               names_to = c("Var","Stat"), names_sep = "_",
               values_to = "Value") %>%
  pivot_wider(names_from = Stat, values_from = Value) %>%
  mutate(Var2 = factor(Var,
                       levels = c("NO2","Pressure","Rainfall","RelHumid","Temperature","WindU","WindV"),
                       labels = c(attributes(latex2exp::TeX("NO2 ($\\mu$g/$m^3$)"))$plotmath,
                                  attributes(latex2exp::TeX("Pressure (hPa)"))$plotmath,
                                  attributes(latex2exp::TeX("Rainfall (mm)"))$plotmath,
                                  attributes(latex2exp::TeX("RelHumid (%)"))$plotmath,
                                  attributes(latex2exp::TeX("Temperature (°)"))$plotmath,
                                  attributes(latex2exp::TeX("WindU (m/s)"))$plotmath,
                                  attributes(latex2exp::TeX("WindV (m/s)"))$plotmath)))


cols <- c("q01-99" = "grey80", "q05-95" = "grey60", "q25-75" = "grey40",
          "median" = "black", "mean" = "blue", "std.dev." = "red")


p1 <- SummData %>%
  filter(Var != "NO2") %>%
  ggplot(mapping = aes(x = h)) + 
  geom_ribbon(mapping = aes(ymin = q1, ymax = q99, fill = "q01-99")) + 
  geom_ribbon(mapping = aes(ymin = q5, ymax = q95, fill = "q05-95")) + 
  geom_ribbon(mapping = aes(ymin = q25, ymax = q75, fill = "q25-75")) + 
  geom_line(mapping = aes(y = median, col = "median"), size = 1.2) + 
  geom_line(mapping = aes(y = mean, col = "mean"), size = 1.2) + 
  facet_wrap(~ Var2, scales = "free", labeller = label_parsed) + 
  theme(legend.position = "bottom") + 
  scale_fill_manual(values = cols, name="") + 
  scale_color_manual(values = cols, name="") + 
  theme_minimal() + 
  theme(plot.title = element_text(size = 14,face = "bold")) + 
  scale_x_continuous(breaks = c(0,6,12,18,23)) + 
  labs(y = "", x = "Hour of the day", title = "Functional boxplots of the covariates",
       subtitle = "Hourly data from 1st March to 31st May 2020 for Lombardy (84 stations)")
ggpubr::ggexport(p1,width = 1500, height = 1000, res = 220, filename = "FunBoxPlot_covs.png")

p2 <- SummData %>%
  filter(Var == "NO2") %>%
  ggplot(mapping = aes(x = h)) + 
  geom_ribbon(mapping = aes(ymin = q1, ymax = q99, fill = "q01-99")) + 
  geom_ribbon(mapping = aes(ymin = q5, ymax = q95, fill = "q05-95")) + 
  geom_ribbon(mapping = aes(ymin = q25, ymax = q75, fill = "q25-75")) + 
  geom_line(mapping = aes(y = median, col = "median"), size = 1.2) + 
  geom_line(mapping = aes(y = mean, col = "mean"), size = 1.2) + 
  geom_line(mapping = aes(y = sd, col = "std.dev."), size = 1.2) + 
  facet_wrap(~ Var2, scales = "free") + 
  theme(legend.position = "bottom") + 
  scale_fill_manual(values = cols, name="") + 
  scale_color_manual(values = cols, name="") + 
  scale_x_continuous(breaks = c(0,6,12,18,23)) + 
  scale_y_continuous(breaks = seq(from = 0, to = 80, by = 5)) + 
  labs(y = latex2exp::TeX("$\\mu$g/$m^3$"), x = "Hour of the day",
       title = latex2exp::TeX("Functional boxplots of NO2 concentrations ($\\mu$g/$m^3$)"),
       subtitle = "Hourly data from 1st March to 31st May 2020 for Lombardy (84 stations)") + 
  theme_minimal() + 
  theme(strip.text.x = element_blank(),
        plot.title = element_text(size = 14,face = "bold"))
p2
ggpubr::ggexport(p2,width = 1500, height = 1000, res = 220, filename = "FunBoxPlot_no2.png")


p3 <- SummData %>%
  filter(Var == "NO2") %>%
  ggplot(mapping = aes(x = h)) + 
  geom_line(mapping = aes(y = sd, col = "std.dev."), size = 1.2) + 
  geom_line(mapping = aes(y = mean, col = "mean"), size = 1.2) + 
  geom_line(mapping = aes(y = median, col = "median"), size = 1.2) + 
  facet_wrap(~ Var2, scales = "free") + 
  theme(legend.position = "bottom") + 
  scale_fill_manual(values = cols, name="") + 
  scale_color_manual(values = cols, name="") + 
  scale_x_continuous(breaks = c(0,6,12,18,23)) + 
  labs(y = latex2exp::TeX("$\\mu$g/$m^3$"), x = "Hour of the day",
       title = latex2exp::TeX("Average, median and variability of NO2 concentrations ($\\mu$g/$m^3$)"),
       subtitle = "Hourly data from 1st March to 31st May 2020 for Lombardy (84 stations)") + 
  theme_minimal() + 
  theme(strip.text.x = element_blank(),
        plot.title = element_text(size = 14,face = "bold"))
p3
ggexport(p3,width = 1500, height = 1000, res = 220, filename = "NO2_dynamics.png")


p23 <- ggarrange(p2,p3,ncol = 2)
p23
