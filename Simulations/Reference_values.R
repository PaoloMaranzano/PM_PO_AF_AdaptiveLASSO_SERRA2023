library(xtable)
library(readxl)
library(tidyverse)

setwd("H:/.shortcut-targets-by-id/1hRMzSqAcE5AcmrsOu5k3nw1rmUY1NHyO/QA&COVID19/VisitingLUH2021/Code/funHDGM_sim")

SimRefVals <- read_excel("SimRefVals.xlsx")

xtable(x = SimRefVals,
       label = "Ref_values",
       caption = "Simulated (n=100 simulations) reference values for each setting",
       digits = 4)
