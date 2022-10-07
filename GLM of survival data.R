## set WD
# setwd() # set the working directory

## load packages
library(survival)
library(survminer)
library(car) # logit transforming percentage data
library(GGally) # data exploration
library(dplyr) # select function

## load data
therm_surv_LF <- read.csv("gammroes_temp_survival_LF_HD.csv", header = T, sep = ",", row.names = c(1) ,stringsAsFactors = FALSE)
therm_surv_LF$popf <- as.factor(therm_surv_LF$pop)
therm_surv_LF$motuf <- as.factor(therm_surv_LF$motu)
therm_surv_LF$repf <- as.factor(therm_surv_LF$rep)
therm_surv_LF$runf <- as.factor(therm_surv_LF$run)

## Quick & dirty GLM
ggplot(therm_surv_LF, aes(temp, surv_perc, col = motu)) +
  geom_point() +
  geom_smooth(method = "glm", se = T,
              method.args = list(family = "quasibinomial"), linetype = "dashed")

ggplot(therm_surv_LF, aes(temp, death_perc, col = motu)) +
  geom_point() +
  geom_smooth(method = "glm", se = T,
              method.args = list(family = "quasibinomial"), linetype = "dashed")

### All data - long format
## data exploration
# change plot size (optional)
options(repr.plot.width = 20, repr.plot.height = 10)
therm_surv_LF %>% 
  select("death_perc", "surv_perc", "temp", "popf", "motuf", "repf", "runf") %>%
  ggpairs(mapping = aes(color = therm_surv_LF$motu, alpha = 0.5))
### All data combined
## GLM with quasibinomial distribution to account for overdispersion
plot(therm_surv_LF$surv_perc ~ therm_surv_LF$temp, pch = 19, col = factor(therm_surv_LF$motuf))
mod <- glm(surv_perc ~ temp, data = therm_surv_LF, family = binomial)
summary(mod)
mod2 <- glm(surv_perc ~ temp, data = therm_surv_LF, family = quasibinomial)
summary(mod2) # dispersion parametre = 0.1644116
## visualising the glm
MyData <- data.frame(temp = 
                       seq(from = min(therm_surv_LF$temp),
                           to = max(therm_surv_LF$temp), by = 4))
P1 <- predict(mod2, newdata = MyData, type = "link", se = TRUE)
plot(MyData$temp, exp(P1$fit) / (1+exp(P1$fit)),
     type = "l", ylim = c(0, 1),
     xlab = "Temperature (degrees Celcius)",
     ylab = "Probability of survival")
lines(MyData$temp, exp(P1$fit + 1.96 * P1$se.fit) / 
      (1 + exp(P1$fit + 1.96 * P1$se.fit)), lty = 2)
lines(MyData$temp, exp(P1$fit - 1.96 * P1$se.fit) / 
        (1 + exp(P1$fit - 1.96 * P1$se.fit)), lty = 2)
points(therm_surv_LF$temp, therm_surv_LF$surv_perc)
ggplot(therm_surv_LF, aes(temp, surv_perc)) +
  geom_point() +
  geom_smooth(method = "glm", se = T,
              method.args = list(family = "quasibinomial"), linetype = "dashed")
## Model validation
# standardized model outputs
op <- par(mfrow = c(2,2))
plot(mod2, which = 1:4)
# residuals
EP <- resid(mod2, type = "pearson")
ED <- resid(mod2, type = "deviance")
mu <- predict(mod2, type = "response")
E <- therm_surv_LF$surv_perc - mu
EP2 <- E / sqrt(0.3328222 * mu)
plot(x = mu, y = E, main = "Response residuals")
abline(0, 0, col = "red")
plot(x = mu, y = EP, main = "Pearson residuals")
abline(0, 0, col = "red")
plot(x = mu, y = EP2, main = "Pearson residuals scaled")
abline(0, 0, col = "red")
plot(x = mu, y = ED, main = "Deviance residuals")
abline(0, 0, col = "red")
par(op)

### At motu level
## Motu G
## GLM with quasibinomial distribution to account for overdispersion
plot(therm_surv_LF[c(1:135),]$surv_perc ~ therm_surv_LF[c(1:135),]$temp, pch = 19, col = factor(therm_surv_LF[c(1:135),]$motuf))
mod3 <- glm(surv_perc ~ temp, data = therm_surv_LF[c(1:135),], family = binomial)
summary(mod3)
mod4 <- glm(surv_perc ~ temp, data = therm_surv_LF[c(1:135),], family = quasibinomial)
summary(mod4) # dispersion parametre = 0.4378245
## generate prediction frame
ggplot(therm_surv_LF[c(1:135),], aes(temp, surv_perc, col = motuf)) +
  geom_point() +
  xlim(22, 34) +
  geom_smooth(method = "glm", se = T,
              method.args = list(family = "quasibinomial"), linetype = "dashed")  ## use prediction data here
## visualising the glm
MyData <- data.frame(temp = 
                       seq(from = min(therm_surv_LF[c(1:135),]$temp),
                           to = max(therm_surv_LF[c(1:135),]$temp), by = 4))
P1 <- predict(mod4, newdata = MyData, type = "link", se = TRUE)
plot(MyData$temp, exp(P1$fit) / (1+exp(P1$fit)),
     type = "l", ylim = c(0, 1),
     xlab = "Temperature (degrees Celcius)",
     ylab = "Probability of survival")
lines(MyData$temp, exp(P1$fit + 1.96 * P1$se.fit) / 
        (1 + exp(P1$fit + 1.96 * P1$se.fit)), lty = 2)
lines(MyData$temp, exp(P1$fit - 1.96 * P1$se.fit) / 
        (1 + exp(P1$fit - 1.96 * P1$se.fit)), lty = 2)
points(therm_surv_LF$temp, therm_surv_LF$surv_perc)
ggplot(therm_surv_LF[c(1:135),], aes(temp, surv_perc)) +
  geom_point() +
  geom_smooth(method = "glm", se = T,
              method.args = list(family = "quasibinomial"), linetype = "dashed")
## Model validation
# standardized model outputs
op <- par(mfrow = c(2,2))
plot(mod4, which = 1:4)
# residuals
EP <- resid(mod4, type = "pearson")
ED <- resid(mod4, type = "deviance")
mu <- predict(mod4, type = "response")
E <- therm_surv_LF[c(1:135),]$surv_perc - mu
EP2 <- E / sqrt(0.4378245 * mu)
plot(x = mu, y = E, main = "Response residuals")
abline(0, 0, col = "red")
plot(x = mu, y = EP, main = "Pearson residuals")
abline(0, 0, col = "red")
plot(x = mu, y = EP2, main = "Pearson residuals scaled")
abline(0, 0, col = "red")
plot(x = mu, y = ED, main = "Deviance residuals")
abline(0, 0, col = "red")
par(op)

## Motu A
## GLM with quasibinomial distribution to account for overdispersion
plot(therm_surv_LF[c(136:270),]$surv_perc ~ therm_surv_LF[c(136:270),]$temp, pch = 19, col = factor(therm_surv_LF[c(136:270),]$motuf)) # population level
mod5 <- glm(surv_perc ~ temp, data = therm_surv_LF[c(136:270),], family = binomial)
summary(mod5)
mod6 <- glm(surv_perc ~ temp, data = therm_surv_LF[c(136:270),], family = quasibinomial)
summary(mod6) # dispersion parametre = 0.3562261
## generate prediction frame
ggplot(therm_surv_LF[c(136:270),], aes(temp, surv_perc, col = motu)) +
  geom_point() +
  xlim(22, 34) +
  geom_smooth(method = "glm", se = T,
              method.args = list(family = "quasibinomial"), linetype = "dashed")  ## use prediction data here
## visualising the glm
MyData <- data.frame(temp = 
                       seq(from = min(therm_surv_LF[c(136:270),]$temp),
                           to = max(therm_surv_LF[c(136:270),]$temp), by = 4))
P1 <- predict(mod6, newdata = MyData, type = "link", se = TRUE)
plot(MyData$temp, exp(P1$fit) / (1+exp(P1$fit)),
     type = "l", ylim = c(0, 1),
     xlab = "Temperature (degrees Celcius)",
     ylab = "Probability of survival")
lines(MyData$temp, exp(P1$fit + 1.96 * P1$se.fit) / 
        (1 + exp(P1$fit + 1.96 * P1$se.fit)), lty = 2)
lines(MyData$temp, exp(P1$fit - 1.96 * P1$se.fit) / 
        (1 + exp(P1$fit - 1.96 * P1$se.fit)), lty = 2)
points(therm_surv_LF$temp, therm_surv_LF$surv_perc)
ggplot(therm_surv_LF[c(136:270),], aes(temp, surv_perc)) +
  geom_point() +
  geom_smooth(method = "glm", se = T,
              method.args = list(family = "quasibinomial"), linetype = "dashed")
## Model validation
# standardized model outputs
op <- par(mfrow = c(2,2))
plot(mod6, which = 1:4)
# residuals
EP <- resid(mod6, type = "pearson")
ED <- resid(mod6, type = "deviance")
mu <- predict(mod6, type = "response")
E <- therm_surv_LF[c(136:270),]$surv_perc - mu
EP2 <- E / sqrt(0.3562261 * mu)
plot(x = mu, y = E, main = "Response residuals")
abline(0, 0, col = "red")
plot(x = mu, y = EP, main = "Pearson residuals")
abline(0, 0, col = "red")
plot(x = mu, y = EP2, main = "Pearson residuals scaled")
abline(0, 0, col = "red")
plot(x = mu, y = ED, main = "Deviance residuals")
abline(0, 0, col = "red")
par(op)

## Motu L
## GLM with quasibinomial distribution to account for overdispersion
plot(therm_surv_LF[c(271:405),]$surv_perc ~ therm_surv_LF[c(271:405),]$temp, pch = 19, col = factor(therm_surv_LF[c(271:405),]$motuf))
mod7 <- glm(surv_perc ~ temp, data = therm_surv_LF[c(271:405),], family = binomial)
summary(mod5)
mod8 <- glm(surv_perc ~ temp, data = therm_surv_LF[c(271:405),], family = quasibinomial)
summary(mod6) # dispersion parametre = 0.3562261
## generate prediction frame
ggplot(therm_surv_LF[c(271:405),], aes(temp, surv_perc, col = motu)) +
  geom_point() +
  xlim(22, 34) +
  geom_smooth(method = "glm", se = T,
              method.args = list(family = "quasibinomial"), linetype = "dashed")  ## use prediction data here
## visualising the glm
MyData <- data.frame(temp = 
                       seq(from = min(therm_surv_LF[c(271:405),]$temp),
                           to = max(therm_surv_LF[c(271:405),]$temp), by = 4))
P1 <- predict(mod8, newdata = MyData, type = "link", se = TRUE)
plot(MyData$temp, exp(P1$fit) / (1+exp(P1$fit)),
     type = "l", ylim = c(0, 1),
     xlab = "Temperature (degrees Celcius)",
     ylab = "Probability of survival")
lines(MyData$temp, exp(P1$fit + 1.96 * P1$se.fit) / 
        (1 + exp(P1$fit + 1.96 * P1$se.fit)), lty = 2)
lines(MyData$temp, exp(P1$fit - 1.96 * P1$se.fit) / 
        (1 + exp(P1$fit - 1.96 * P1$se.fit)), lty = 2)
points(therm_surv_LF$temp, therm_surv_LF$surv_perc)
ggplot(therm_surv_LF[c(271:405),], aes(temp, surv_perc)) +
  geom_point() +
  geom_smooth(method = "glm", se = T,
              method.args = list(family = "quasibinomial"), linetype = "dashed")
## Model validation
# standardized model outputs
op <- par(mfrow = c(2,2))
plot(mod8, which = 1:4)
# residuals
EP <- resid(mod8, type = "pearson")
ED <- resid(mod8, type = "deviance")
mu <- predict(mod8, type = "response")
E <- therm_surv_LF[c(271:405),]$surv_perc - mu
EP2 <- E / sqrt(0.3562261 * mu)
plot(x = mu, y = E, main = "Response residuals")
abline(0, 0, col = "red")
plot(x = mu, y = EP, main = "Pearson residuals")
abline(0, 0, col = "red")
plot(x = mu, y = EP2, main = "Pearson residuals scaled")
abline(0, 0, col = "red")
plot(x = mu, y = ED, main = "Deviance residuals")
abline(0, 0, col = "red")
par(op)
