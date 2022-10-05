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
mod <- glm(surv_perc ~ temp, data = therm_surv, family = binomial)
summary(mod)

mod2 <- glm(surv_perc ~ temp, data = therm_surv, family = quasibinomial)
summary(mod2) # dispersion parametre = 0.1644116

## generate prediction frame
ggplot(therm_surv, aes(temp, surv_perc, col = motu)) +
  geom_point() +
  xlim(22, 34) +
  geom_smooth(method = "glm", se = T,
              method.args = list(family = "quasibinomial"), linetype = "dashed")  ## use prediction data here

## visualising the glm
MyData <- data.frame(temp = 
                       seq(from = min(therm_surv$temp),
                           to = max(therm_surv$temp), by = 4))
P1 <- predict(mod2, newdata = MyData, type = "link", se = TRUE)
plot(MyData$temp, exp(P1$fit) / (1+exp(P1$fit)),
     type = "l", ylim = c(0, 1),
     xlab = "Temperature (degrees Celcius)",
     ylab = "Probability of survival")
lines(MyData$temp, exp(P1$fit + 1.96 * P1$se.fit) / 
      (1 + exp(P1$fit + 1.96 * P1$se.fit)), lty = 2)
lines(MyData$temp, exp(P1$fit - 1.96 * P1$se.fit) / 
        (1 + exp(P1$fit - 1.96 * P1$se.fit)), lty = 2)
points(therm_surv$time, therm_surv$surv_perc)

## Model validation
EP <- resid(mod2, type = "pearson")
ED <- resid(mod2, type = "deviance")
mu = predict(mod2, type = "response")
