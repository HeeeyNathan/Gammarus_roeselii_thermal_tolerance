# Read in the survival data
survdat <- read.csv("survival_own_data.csv")
survdat$popf <- as.factor(survdat$pop)

# Inspect distribution of LOS variable
hist(survdat$temp)
hist(subset(survdat, surv == 0)$temp)
hist(subset(survdat, surv == 1)$temp)

# Kaplan-Meier Analysis - life table
library(survival)

# specify KM analysis model
km_fit1 <- survfit(Surv(temp, surv) ~ 1,    # first argument = time to event variable, second 
                   data = survdat,          # argument = survival function. Adding 1 after tilda 
                   type = "kaplan-meier")   # makes null model, if you want a categorical variable, add it there.
print(km_fit1) # print model results
str(km_fit1) # Look at model structure

# for entire span of study survival rate
## 2020 observations
## 1388 is the number of gammarids that experienced the events (i.e., died at a specific temp)
## 28 is the median temperature before dying.

# summarise km analysis results using default time intervals and create life table
summary(km_fit1)
# summarise km analysis results be prespecified time intervals and create life table
summary(km_fit1, times = c(22, 25, 28, 31, 34)) # specifies time intervals of interest (3 degree increments)
# plot cumulative survival rates (probabilities)
plot(km_fit1) # x-axis is temp, y axis is survival.

# Aestically pleasing plots
library(survminer)
library(ggsurvfit)
library(lubridate)
library(gtsummary)
library(tidycmprsk)
library(ggplot2)
library(condsurv)

# Plot km output
ggsurvplot(km_fit1, data = survdat, risk.table = T, conf.int = T,
           ggtheme = theme_minimal()) # strata = non-independent groups of people

survfit2(Surv(temp, surv) ~ 1, data = survdat) %>% 
# build Kaplan-Meier plot ----------------------------------------------------
  ggsurvfit(size = 1) +
  add_confidence_interval() +
  add_risktable()+
  # add_quantile(y_value = 0.6, color = "gray50", size = 0.75)
  # use ggplot2 functions to style the plot and update the labels --------------
# limit plot to only show relevant temperatures
coord_cartesian(xlim = c(22, 34)) +
  # update figure labels/titles
  labs(
    y = "Percentage Survival",
    title = "Temperature (in degrees Celcius)",
  ) +
  # reduce padding on edges of figure and format axes
  scale_y_continuous(label = scales::percent, 
                     breaks = seq(0, 1, by = 0.1),
                     expand = c(0.02, 0)) +
  scale_x_continuous(breaks = seq(22, 34, by = 3), 
                     expand = c(0.04, 0))

# specify KM analysis model with a categorical covariate
km_fit2 <- survfit(Surv(temp, surv) ~ motu, # first argument = time to event variable, second 
                   data = survdat,          # argument = survival function. Adding 1 after tilda 
                   type = "kaplan-meier")   # makes null model, if you want a categorical variable, add it there.
print(km_fit2) # print model results
str(km_fit2) # Look at model structure

# summarise km analysis results using default time intervals and create life table
summary(km_fit2)

# summarise km analysis results be prespecified time intervals and creat life table
summary(km_fit2, times = c(22, 25, 28, 31, 34)) # specifies time intervals of interest (3 degree increments)

# plot cumulative survival rates (probabilities)
plot(km_fit2) # x-axis is temperature, y axis is survival.

# Aestically pleasing plots
ggsurvplot(km_fit2, data = survdat, risk.table = F, conf.int = T, # are the cumulative survival rate curves 
           pval = T, pval.method = T,                             # different between the races, pval.method = T = default is log rank methd 
           ggtheme = theme_minimal())

survfit2(Surv(temp, surv) ~ motu, data = survdat) %>% 
# build Kaplan-Meier plot ----------------------------------------------------
ggsurvfit(size = 1) +
  add_confidence_interval() +
  # add_risktable()+
add_quantile(x_value = 28, color = "gray50", size = 0.75) + # median survival time
# use ggplot2 functions to style the plot and update the labels --------------
  # limit plot to only show relevant temperatures
coord_cartesian(xlim = c(22, 34)) +
  # update figure labels/titles
  labs(
    y = "Percentage Survival",
    title = "Temperature (in degrees Celcius)",
  ) +
  # reduce padding on edges of figure and format axes
  scale_y_continuous(label = scales::percent, 
                     breaks = seq(0, 1, by = 0.1),
                     expand = c(0.02, 0)) +
  scale_x_continuous(breaks = seq(22, 34, by = 3), 
                     expand = c(0.04, 0))

# estimating x-temperature survival
summary(km_fit2, times = 32) # surival probability at 32 degrees

# comparing survival between groups
survdiff(Surv(time, surv) ~ motu, data = survdat)

# estimate cox proportional hazards (ph) model (cox regression model)
cox_reg1 <- coxph(Surv(temp, surv) ~ motu,
                  data = survdat) # dont interpret the second lines of information if the p-values of the first lines are not significant
print(cox_reg1)
summary(cox_reg1)
# Concordance score = the performance of the model = how well did the model predict the survival rate
# = proportion of pairs of gammarids that died and how well the model predicted it
# a concordance of 0.5 means that model does no better than chance (50% chance of model getting this correct)
# concordance of 0.511 = that the race categorical variable slighly inproves the accuracy of our 
# predictions but not by too much (better than nothing)

## LOG RANK TEST IS RECOMMENED FOR SMALLER SAMPLE SIZES!!!!!!!!!!!!!!!
# the model with the motu covariate does not perform significantly better than a model with no covariates.

# estimate cox proportional hazards (ph) model (cox regression model)
cox_reg2 <- coxph(Surv(temp, surv) ~ pop + run + rep + motu + popf,
                  data = survdat) # don't interpret the second lines of information if the p-values of the first lines are not significant
print(cox_reg2)
summary(cox_reg2) # pay satisfactions data had 38 missing observations

# what is the hazard ratio?
# Harzard ratio = exp(-0.35287) = 0.70267
# 1 - HR = 1 - 0.70267 =  0.29733 = the risk of an individual dying will decrease by 29.7% 
# in motu G compared with motu A (which is 0).

# reciprocal of HR = 1 / 0.70267 = 1.423143 = the chance of an individual staying alive will increase 1.423143 times 
# for motu G compared with motu A (which is based on 0). If you do present the reciprocal, you should defintiely 
# report the HR because the confidence intervals are related to the HR and not the reciprocal.

## CHOOSE WHAT IS EASIEST TO EXPLAIN TO PEOPLE

# Nested model comparisons (using anova)
cox_reg1 # 701 individuals in the sample size
cox_reg2 # sample is 663 (38 individuals removed)

anova(cox_reg1, cox_reg2) # cox_reg2 = model with motu is better than the fuller model with many covariates



