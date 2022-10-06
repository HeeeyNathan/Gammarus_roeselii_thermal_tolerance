## set WD
# setwd() # set the working directory

## load packages
library(survival)
library(survminer)
library(car) # logit transforming percentage data
library(GGally) # data exploration
library(dplyr) # select function

## load data
therm_surv_LF <- read.csv("gammroes_temp_survival_LF.csv", header = T, sep = ",", row.names = c(1) ,stringsAsFactors = FALSE)
therm_surv_LF$popf <- as.factor(therm_surv_LF$pop)
therm_surv_LF$motuf <- as.factor(therm_surv_LF$motu)
therm_surv_LF$repf <- as.factor(therm_surv_LF$rep)

### All data - long format
## data exploration
# change plot size (optional)
options(repr.plot.width = 20, repr.plot.height = 10)
therm_surv_LF %>% 
  select("surv", "surv_perc", "temp", "pop", "motu", "rep") %>%
  ggpairs(mapping = aes(color = therm_surv_LF$motu, alpha = 0.5))
### All data combined
## GLM with quasibinomial distribution to account for overdispersion
plot(therm_surv_LF$surv_perc ~ therm_surv_LF$temp, pch = 19, col = "black")
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
points(therm_surv_LF$time, therm_surv_LF$surv_perc)
## Model validation
# standardized model outputs
op <- par(mfrow = c(2,2))
plot(mod2, which = 1:4)
# residuals
EP <- resid(mod2, type = "pearson")
ED <- resid(mod2, type = "deviance")
mu <- predict(mod2, type = "response")
E <- therm_surv_LF$surv_perc - mu
EP2 <- E / sqrt(0.1644116 * mu)
plot(x = mu, y = E, main = "Response residuals")
abline(0, 0, col = "red")
plot(x = mu, y = EP, main = "Pearson residuals")
abline(0, 0, col = "red")
plot(x = mu, y = EP2, main = "Pearson residuals scaled")
abline(0, 0, col = "red")
plot(x = mu, y = ED, main = "Deviance residuals")
abline(0, 0, col = "red")
par(op)

### At population level
## GLM with quasibinomial distribution to account for overdispersion
plot(therm_surv_LF$surv_perc ~ therm_surv_LF$temp, pch = 19, col = factor(therm_surv_LF$popf))
mod3 <- glm(surv_perc ~ temp + popf + motuf + repf, data = therm_surv_LF, family = quasibinomial)
summary(mod3) # dispersion parametre = 0.1644116
## generate prediction frame
ggplot(therm_surv_LF, aes(temp, surv_perc, col = popf)) +
  geom_point() +
  geom_smooth(method = "glm", se = F,
              method.args = list(family = "quasibinomial"), linetype = "dashed")
ggplot(therm_surv_LF, aes(temp, surv_perc, col = popf)) +
  geom_point() +
  geom_smooth(method = "glm", se = T,
              method.args = list(family = "quasibinomial"), linetype = "dashed")
## visualising the glm
MyData <- data.frame(temp = 
                       seq(from = min(therm_surv_LF$temp),
                           to = max(therm_surv_LF$temp), by = 4))
P1 <- predict(mod3, newdata = MyData, type = "link", se = TRUE)
plot(MyData$temp, exp(P1$fit) / (1+exp(P1$fit)),
     type = "l", ylim = c(0, 1),
     xlab = "Temperature (degrees Celcius)",
     ylab = "Probability of survival")
lines(MyData$temp, exp(P1$fit + 1.96 * P1$se.fit) / 
        (1 + exp(P1$fit + 1.96 * P1$se.fit)), lty = 2)
lines(MyData$temp, exp(P1$fit - 1.96 * P1$se.fit) / 
        (1 + exp(P1$fit - 1.96 * P1$se.fit)), lty = 2)
points(therm_surv_LF$time, therm_surv_LF$surv_perc)

MyData2 <- data.frame(pop = 
                       seq(from = min(therm_surv_LF$pop),
                           to = max(therm_surv_LF$pop), by = 1))
P1 <- predict(mod3, newdata = MyData2, type = "link", se = TRUE)
plot(MyData$temp, exp(P1$fit) / (1+exp(P1$fit)),
     type = "l", ylim = c(0, 1),
     xlab = "Temperature (degrees Celcius)",
     ylab = "Probability of survival")
lines(MyData$temp, exp(P1$fit + 1.96 * P1$se.fit) / 
        (1 + exp(P1$fit + 1.96 * P1$se.fit)), lty = 2)
lines(MyData$temp, exp(P1$fit - 1.96 * P1$se.fit) / 
        (1 + exp(P1$fit - 1.96 * P1$se.fit)), lty = 2)
points(therm_surv_LF$time, therm_surv_LF$surv_perc)
## Model validation
# standardized model outputs
op <- par(mfrow = c(2,2))
plot(mod2, which = 1:4)
# residuals
EP <- resid(mod2, type = "pearson")
ED <- resid(mod2, type = "deviance")
mu <- predict(mod2, type = "response")
E <- therm_surv_LF$surv_perc - mu
EP2 <- E / sqrt(0.1644116 * mu)
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
## GLM with quasibinomial distribution to account for overdispersion
plot(therm_surv_LF$surv_perc ~ therm_surv_LF$temp, pch = 19, col = factor(therm_surv_LF$motu)) # population level
plot(therm_surv_LF$surv_perc ~ therm_surv_LF$temp,
     pch = 19,
     col = factor(therm_surv_LF$motu)) # motu level
mod <- glm(surv_perc ~ temp, data = therm_surv_LF, family = binomial)
summary(mod)
mod2 <- glm(surv_perc ~ temp, data = therm_surv_LF, family = quasibinomial)
summary(mod2) # dispersion parametre = 0.1644116
## generate prediction frame
ggplot(therm_surv_LF, aes(temp, surv_perc, col = motu)) +
  geom_point() +
  xlim(22, 34) +
  geom_smooth(method = "glm", se = T,
              method.args = list(family = "quasibinomial"), linetype = "dashed")  ## use prediction data here
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
points(therm_surv_LF$time, therm_surv_LF$surv_perc)
## Model validation
# standardized model outputs
op <- par(mfrow = c(2,2))
plot(mod2, which = 1:4)
# residuals
EP <- resid(mod2, type = "pearson")
ED <- resid(mod2, type = "deviance")
mu <- predict(mod2, type = "response")
E <- therm_surv_LF$surv_perc - mu
EP2 <- E / sqrt(0.1644116 * mu)
plot(x = mu, y = E, main = "Response residuals")
abline(0, 0, col = "red")
plot(x = mu, y = EP, main = "Pearson residuals")
abline(0, 0, col = "red")
plot(x = mu, y = EP2, main = "Pearson residuals scaled")
abline(0, 0, col = "red")
plot(x = mu, y = ED, main = "Deviance residuals")
abline(0, 0, col = "red")
par(op)