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
cox_reg2 <- coxph(Surv(temp, surv) ~ motu + popf,
                  data = survdat) # don't interpret the second lines of information if the p-values of the first lines are not significant
print(cox_reg2)
summary(cox_reg2) # pay satisfactions data had 38 missing observations

# what is the hazard ratio?
# Harzard ratio = exp(-0.10971) = 0.8961
# 1 - HR = 1 - 0.8961 =  0.1039 = the risk of an individual quitting their job will decrease by 10.4% 
# for every additional dollar in hourly pay that the individual earned.

# reciprocal of HR = 1 / 0.8961 = 1.1160 = risk of staying in the organisation will increase 1.116 times 
# for every additional dollar of hourly payed earned. If you do present the reciprocal, you should defintiely 
# report the HR because the confidence intervals are related to the HR and not the reciprocal.

## CHOOSE WHAT IS EASIEST TO EXPLAIN TO PEOPLE

# write out the above in equation format
# hypothetical values for a hispranic latino person who makes 16$ per hour and is satisfied with their pay
RaceHL <- 1
RaceW <- 0
PayHour <- 16.00
PaySat <- 4.00
log_overallrisk <- .121*RaceHL -.044*RaceW -.110*PayHour + .141*PaySat
print(log_overallrisk)
exp(log_overallrisk) # the overall risk of quitting their job is 65.1% (1 - 0.3412978 = 0.6587022)
                     # less likely to quit their job when compared to an individual with scores of zero
                     # in each of the covariates in the model.

# making covariates more meaningful
# grand-mean centre continuous covariates
survdat$c_Pay_hourly <- scale(survdat$Pay_hourly, center = T, scale = F) # now model will compare to the average not 0 which is meaningless
survdat$c_Pay_sat <- scale(survdat$Pay_sat, center = T, scale = F)# now model will compare to the average not 0 which is meaningless

# Change reference group to HispanicLatino for categorical covariate by reordering levels
survdat$HL_Race <- factor(survdat$Race, levels = c("HispanicLatino", "Black", "White"))

# Estimate Cox proportional hazard model with categorical & continuous covariates
cox_reg3 <- coxph(Surv(LOS, censored) ~ HL_Race + c_Pay_hourly + c_Pay_sat,
                  data = survdat)
summary(cox_reg3)

# hypothetical values for a hispranic latino person who makes 16$ per hour and is satisfied with their pay
RaceB <- 0 # HL is now the reference group - they identify as HL, not W or B
RaceW <- 0
PayHour <- 16.00 - mean(survdat$Pay_hourly, na.rm = T) # drop people with missing data (na.rm)
PayHour # this person has above average pay - the number is positive
PaySat <- 4.00 - mean(survdat$Pay_sat, na.rm = T) # drop people with missing data (na.rm)
PaySat # this person has above average pay satisfaction - the number is positive

log_overallrisk <- -.121*RaceB -.165*RaceW -.110*PayHour + .141*PaySat
print(log_overallrisk) # -0.1623714
exp(log_overallrisk) # the overall risk of quitting their job is 15% (1 - 0.8501254 = 0.1498746)
# less likely to quit their job when compared to an individual with scores of zero
# in each of the covariates in the model (i.e., compared to hispanicLatinos with average hourly pay 
# and have average levels of pay satisfaction).

# Nested model comparisons (using anova)
cox_reg1 # 701 individuals in the sample size
cox_reg2 # sample is 663 (38 individuals removed)

library(tidyr)

# estimate cox proportional hazards (ph) model (cox regression model)
cox_reg1 <- coxph(Surv(LOS, censored) ~ Race,
                  data = drop_na(survdat, LOS, censored, Race, Pay_hourly, Pay_sat)) # dont interpret the second lines of information if the p-values of the first lines are not significant
summary(cox_reg1)
# estimate cox proportional hazards (ph) model (cox regression model)
cox_reg2 <- coxph(Surv(LOS, censored) ~ Race + Pay_hourly + Pay_sat,
                  data = drop_na(survdat, LOS, censored, Race, Pay_hourly, Pay_sat)) # don't interpret the second lines of information if the p-values of the first lines are not significant
summary(cox_reg2)

anova(cox_reg1, cox_reg2) # cox_reg2 = full model = is better than the smaller nested model
