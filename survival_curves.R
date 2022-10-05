## set WD
# setwd() # set the working directory

## load packages
library(survival)
library(survminer)
library(car) # logit transforming percentage data
library(GGally) # data explorartion
library(dplyr) # select function

## load data
therm_surv <- read.csv("gammroes_temp_survival_LF.csv", header = T, sep = ",", row.names = c(1) ,stringsAsFactors = FALSE)
therm_surv$pop <- as.factor(therm_surv$pop)
therm_surv$motu <- as.factor(therm_surv$motu)
therm_surv$rep <- as.factor(therm_surv$rep)

## data exploration
# change plot size (optional)
options(repr.plot.width = 20, repr.plot.height = 10)

therm_surv %>% 
  select("surv", "surv_perc", "temp", "pop", "motu", "rep") %>%
  ggpairs(mapping = aes(color = therm_surv$motu, alpha = 0.5))

## GLM with binomial distribution
plot(therm_surv$surv_perc ~ therm_surv$temp,
     pch = 19,
     col = factor(therm_surv$pop))
mod <- glm(surv_perc ~ temp, data = therm_surv, family = "binomial")
summary(mod)

fun.gen <- function(awd) exp(mod$coef[1] + mod$coef[2] * awd)
fun.acd <- function(awd) exp(mod$coef[1] + mod$coef[2] * awd + mod$coef[3])
fun.voc <- function(awd) exp(mod$coef[1] + mod$coef[2] * awd + mod$coef[4])

## generate prediction frame
ggplot(therm_surv, aes(temp, surv_perc, col = motu)) +
  geom_point() +
  geom_smooth(method = "glm", se = TRUE,
              method.args = list(family = "binomial"), linetype = "dashed")  ## use prediction data here
