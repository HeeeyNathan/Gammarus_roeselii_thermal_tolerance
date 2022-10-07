# Read in the survival data
survdat <- read.csv("Survival_example.csv")

# create sensoring variable (right censoring)
survdat$censored[survdat$Turnover == 1] <- 1
survdat$censored[survdat$Turnover == 0 | survdat$Turnover == 2] <- 0

# Inspect distribution of LOS variable
hist(survdat$LOS)
hist(subset(survdat, Turnover == 0)$LOS)
hist(subset(survdat, Turnover == 1)$LOS)

# Kaplan-Meier Analysis - life table
library(survival)

# specify KM analysis model
km_fit1 <- survfit(Surv(LOS, censored) ~ 1, # first argument = time to event variable, second 
                   data = survdat,          # argument = survival function. Adding 1 after tilda 
                   type = "kaplan-meier")   # makes null model, if you want a categorical variable, add it there.
print(km_fit1)
# for entire span of study survival rate
## 701 observations
## 463 is the number of people who quit the organisation
## 1059 is the median number of days before quitting the organisation

# summarise km analysis results using default time intervals and create life table
summary(km_fit1)

# summarise km analysis results be prespecified time intervals and creat life table
summary(km_fit1, times = c(30, 60, 90*(1:30))) # specifies time intervals of interest (90 day increments)

# plot cumulative survival rates (probabilities)
plot(km_fit1) # x-axis is LOS, y axis is survival.

# Aestically pleasing plots
library(survminer)

ggsurvplot(km_fit1, data = survdat, risk.table = T, conf.int = T,
           ggtheme = theme_minimal()) # strata = non-independent groups of people

# Kaplan-Meier Analysis - life table
library(survival)

# specify KM analysis model with a categorical covariate
km_fit2 <- survfit(Surv(LOS, censored) ~ Race, # first argument = time to event variable, second 
                   data = survdat,          # argument = survival function. Adding 1 after tilda 
                   type = "kaplan-meier")   # makes null model, if you want a categorical variable, add it there.
print(km_fit2)

# summarise km analysis results using default time intervals and create life table
summary(km_fit2)

# summarise km analysis results be prespecified time intervals and creat life table
summary(km_fit2, times = c(30, 60, 90*(1:30))) # specifies time intervals of interest (90 day increments)

# plot cumulative survival rates (probabilities)
plot(km_fit2) # x-axis is LOS, y axis is survival.

# Aestically pleasing plots
ggsurvplot(km_fit2, data = survdat, risk.table = T, conf.int = T, # are the cumulative survival rate curves 
           pval = T, pval.method = T,                             # different between the races, pval.method = T = default is log rank methd 
           ggtheme = theme_minimal())

# estimate cox proportional hazards (ph) model (cox regression model)
cox_reg1 <- coxph(Surv(LOS, censored) ~ Race,
                  data = survdat) # dont interpret the second lines of information if the p-values of the first lines are not significant
print(cox_reg1)
summary(cox_reg1)
# Concordance score = the performance of the model = how well did the model predict the survival rate
# = proportion of pairs of individuals that quit job and how well the model predicted it
# a concordance of 0.5 means that model does no better than chance (50% chance of model getting this correct)
# concordance of 0.522 = that the race categorical variable slighly inproves the accuracy of our 
# predictions but not by too much (better than nothing)

## LOG RANK TEST IS RECOMMENED FOR SMALLER SAMPLE SIZES!!!!!!!!!!!!!!!
# the model with the race covariate does not perform significantly better than a model with no covariates.

# estimate cox proportional hazards (ph) model (cox regression model)
cox_reg1 <- coxph(Surv(LOS, censored) ~ Race + Pay_hourly + Pay_sat,
                  data = survdat) # dont interpret the second lines of information if the p-values of the first lines are not significant
print(cox_reg1)
summary(cox_reg1) # pay satisfactions data had 38 missing observations
