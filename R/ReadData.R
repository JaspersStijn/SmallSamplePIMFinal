library(readxl)
library(tidyverse)
library(MASS)

source("./R/Main_Function.R")

## Read in data
dat = read_xlsx("./Data/Pyridine.xlsx")
dat = dat %>% mutate(Dose = as.factor(Dose),
                     Sex = factor(Sex,labels=c("M","F")),
                     weight_use = as.numeric(weight_use))
