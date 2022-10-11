## set WD
# setwd() # set the working directory

## load data
therm_surv_LF <- read.csv("Survival_own_data.csv", header = T, sep = ",", row.names = c(1) ,stringsAsFactors = FALSE)
therm_surv_LF$popf <- as.factor(therm_surv_LF$pop)
therm_surv_LF$motuf <- as.factor(therm_surv_LF$motu)
therm_surv_LF$repf <- as.factor(therm_surv_LF$rep)
therm_surv_LF$runf <- as.factor(therm_surv_LF$run)

## data exploration
boxplot(surv ~ temp, data = therm_surv_LF)
plot(surv ~ temp,  data = therm_surv_LF, main = "Scatterplot", xlab = "Temperature in degrees Celcius", ylab = "Survival", pch = 19)
abline(lm(surv ~ temp, data = therm_surv_LF), col = "red") # regression line (y~x)
lines(lowess(therm_surv_LF$temp, therm_surv_LF$surv), col = "blue") # lowess line (x,y)

nrow(therm_surv_LF[therm_surv_LF$surv == 1,]) # 1388 alive individuals (1)
nrow(therm_surv_LF[therm_surv_LF$surv == 0,]) # 632 dead individuals (0)

## GLM with binomial distribution - all data
GLM_bi1 <- glm(surv ~ temp, data = therm_surv_LF, family = "binomial"(link = "logit")) # assumes ~ equal number of 1s and 0s
summary(GLM_bi1, corr = T) # in a Bernoulli GLM (binomial GLM) because the data are code as 0 and 1, overdispersion cannot occur.
GLM_bi2 <- glm(surv ~ temp, data = therm_surv_LF, family = "binomial"(link = "probit")) # assumes ~ equal number of 1s and 0s
summary(GLM_bi2, corr = T) # in a Bernoulli GLM (binomial GLM) because the data are code as 0 and 1, overdispersion cannot occur.
GLM_bi3 <- glm(surv ~ temp, data = therm_surv_LF, family = "binomial"(link = "cloglog")) # used for considerable more 1s than 0s or vice versa
summary(GLM_bi3, corr = T) # in a Bernoulli GLM (binomial GLM) because the data are code as 0 and 1, overdispersion cannot occur.

# deviance table of best model
anova(GLM_bi1, test = "Chisq")

# model selection
drop1(GLM_bi1, test = "Chi")
# model validation
op <- par(mfrow = c(2,2))
plot(GLM_bi1, which = 1:4)
par(op)

# using the predict command
MyData <- data.frame(temp = 
                       seq(from = min(therm_surv_LF$temp),
                           to = max(therm_surv_LF$temp), by = 0.1))
pred <- predict(GLM_bi1, newdata = MyData, type = "response")
plot(therm_surv_LF$temp, therm_surv_LF$surv,
     xlab = "Temperature (degrees Celcius)",
     ylab = "Probability of survival (%)")
lines(MyData$temp, pred, lty = 1)

## plotting GLM
plot(therm_surv_LF$temp, fitted(glm(surv ~ temp, family = "binomial"(link = "logit"), data = therm_surv_LF)))
GLM_bi1
temp1 <- seq(22,34,0.1)
newtemp <- data.frame(temp=temp1)
predicted.probability <- predict(GLM_bi1, + newtemp,type="resp")
plot(predicted.probability ~ temp1, type="l")

# all data
ggplot(therm_surv_LF, aes(temp, surv)) +
  #geom_point(aes(col = motu)) +
  geom_smooth(method = "glm", se = T,
              method.args = list(family = "binomial"(link = "logit")), linetype = "dashed") +
  coord_cartesian(xlim = c(22, 34)) +
  # update figure labels/titles
  labs(
    y = "Percentage Survival",
    x = "Temperature in degrees Celcius",
    title = "Survival curve - GLM"
  ) +
  # reduce padding on edges of figure and format axes
  scale_y_continuous(label = scales::percent, 
                     breaks = seq(0, 1, by = 0.1),
                     expand = c(0.02, 0)) +
  scale_x_continuous(breaks = seq(22, 34, by = 3), 
                     expand = c(0.04, 0))

## GLM with binomial distribution - temp + motu
GLM_bi4 <- glm(surv ~ temp + motuf, data = therm_surv_LF, family = "binomial"(link = "logit")) # assumes ~ equal number of 1s and 0s
summary(GLM_bi4, corr = T) # in a Bernoulli GLM (binomial GLM) because the data are code as 0 and 1, overdispersion cannot occur.
GLM_bi4.1 <- glm(surv ~ motuf + temp, data = therm_surv_LF, family = "binomial"(link = "logit")) # assumes ~ equal number of 1s and 0s
summary(GLM_bi4.1, corr = T) # in a Bernoulli GLM (binomial GLM) because the data are code as 0 and 1, overdispersion cannot occur.
GLM_bi5 <- glm(surv ~ temp + motuf, data = therm_surv_LF, family = "binomial"(link = "probit")) # assumes ~ equal number of 1s and 0s
summary(GLM_bi5, corr = T) # in a Bernoulli GLM (binomial GLM) because the data are code as 0 and 1, overdispersion cannot occur.
GLM_bi6 <- glm(surv ~ temp + motuf, data = therm_surv_LF, family = "binomial"(link = "cloglog")) # used for considerable more 1s than 0s or vice versa
summary(GLM_bi6, corr = T) # in a Bernoulli GLM (binomial GLM) because the data are code as 0 and 1, overdispersion cannot occur.

## model selection
# deviance table of best model
anova(GLM_bi4, test = "Chisq")
anova(GLM_bi4.1, test = "Chisq") #swopped variables around

# using drop1 function is easier
drop1(GLM_bi4, test = "Chi")

## model validation
op <- par(mfrow = c(2,2))
plot(GLM_bi4, which = 1:4)
par(op)

## Presentation as odds-ratio estimates
exp(cbind(OR=coef(GLM_bi4), confint(GLM_bi4)))

# here, the (Intercept) is really the odds of dying at high temperatures (for not motu A) and not an odds ratio.

# Hardin & Hilde 2018 (4th ed.) - Clog-log and log-log links are asymmetrically sigmoidal. 
# For clog-log models, the upper part of the sigmoid is more elongated or stretched out than 
# the logit or probit. Log-log models are based on the converse. The bottom of the sigmoid is
# elongated or skewed to the left.

# The clog-log inverse link function, which defines the fitted value for the model,
# is also used to predict the probability of death in data that are otherwise modeled
# using Cox proportional hazards models. In a Cox model, given that death is the
# failure event, the probability of death given predictors is a function of time.
# However, if we interpret the data only as those observations having reached the
# end point (death in this example) by a specific time, then we can estimate the
# likelihood of death by that specific time given the various predictors in the
# model as Pr(death) = 1 - exp{-exp(xß)}. One can thus use the log-log inverse link 
# to predict the probability of survival over the same time period.

# There is good reason to suspect that the log-log model is preferable to other binomial links.
# Modeling the grouped data with the log-log link assumes that the data have
# significantly more zeros than ones in the binary form. Accordingly, the resulting
# log-log model should have a lower deviance statistic than that of other links.

# motuf with standard error
ggplot(therm_surv_LF, aes(temp, surv, col = motuf)) +
  #geom_point(aes(col = motu)) +
  geom_smooth(method = "glm", se = T,
              method.args = list(family = "binomial"(link = "logit")), linetype = "dashed") +
  coord_cartesian(xlim = c(22, 34)) +
  # update figure labels/titles
  labs(
    y = "Percentage Survival",
    x = "Temperature in degrees Celcius",
    title = "Survival curve - GLM"
  ) +
  # reduce padding on edges of figure and format axes
  scale_y_continuous(label = scales::percent, 
                     breaks = seq(0, 1, by = 0.1),
                     expand = c(0.02, 0)) +
  scale_x_continuous(breaks = seq(22, 34, by = 3), 
                     expand = c(0.04, 0))

# compare models
anova(GLM_bi1, GLM_bi4, test = "Chi")

## GLM with binomial distribution - temp + motuf + temp*motuf
# Interaction between a factor and a numeric variable. In this case, the
# model with interaction contains linear effects of the continuous variable
# but with different slopes within each group defined by the
# factor.

GLM_bi7 <- glm(surv ~ temp + motuf + temp*motuf, data = therm_surv_LF, family = "binomial"(link = "logit")) # assumes ~ equal number of 1s and 0s
summary(GLM_bi7) # in a Bernoulli GLM (binomial GLM) because the data are code as 0 and 1, overdispersion cannot occur.
GLM_bi8 <- glm(surv ~ temp + motuf + temp*motuf, data = therm_surv_LF, family = "binomial"(link = "probit")) # assumes ~ equal number of 1s and 0s
summary(GLM_bi5) # in a Bernoulli GLM (binomial GLM) because the data are code as 0 and 1, overdispersion cannot occur.
GLM_bi9 <- glm(surv ~ temp + motuf + temp*motuf, data = therm_surv_LF, family = "binomial"(link = "cloglog")) # used for considerable more 1s than 0s or vice versa
summary(GLM_bi6) # in a Bernoulli GLM (binomial GLM) because the data are code as 0 and 1, overdispersion cannot occur.

# cell means
tapply(therm_surv_LF$surv, list(therm_surv_LF$temp, therm_surv_LF$motuf), mean)

# The difference between motus decreases with increasing temperature, making an additive model 
# inadequate. When this is the case, the individual tests for the two factors make no sense. 
# If the interaction had not been significant, then we would have been able to perform separate F
# tests for the two factors.

# model selection
drop1(GLM_bi7, test = "Chi")
# model validation
op <- par(mfrow = c(2,2))
plot(GLM_bi4, which = 1:4)
par(op)

# Hardin & Hilde 2018 (4th ed.) - Clog-log and log-log links are asymmetrically sigmoidal. 
# For clog-log models, the upper part of the sigmoid is more elongated or stretched out than 
# the logit or probit. Log-log models are based on the converse. The bottom of the sigmoid is
# elongated or skewed to the left.

# The clog-log inverse link function, which defines the fitted value for the model,
# is also used to predict the probability of death in data that are otherwise modeled
# using Cox proportional hazards models. In a Cox model, given that death is the
# failure event, the probability of death given predictors is a function of time.
# However, if we interpret the data only as those observations having reached the
# end point (death in this example) by a specific time, then we can estimate the
# likelihood of death by that specific time given the various predictors in the
# model as Pr(death) = 1 - exp{-exp(xß)}. One can thus use the log-log inverse link 
# to predict the probability of survival over the same time period.

# There is good reason to suspect that the log-log model is preferable to other binomial links.
# Modeling the grouped data with the log-log link assumes that the data have
# significantly more zeros than ones in the binary form. Accordingly, the resulting
# log-log model should have a lower deviance statistic than that of other links.

# motuf with standard error
ggplot(therm_surv_LF, aes(temp, surv, col = motuf)) +
  #geom_point(aes(col = motu)) +
  geom_smooth(method = "glm", se = T,
              method.args = list(family = "binomial"(link = "logit")), linetype = "dashed") +
  coord_cartesian(xlim = c(22, 34)) +
  # update figure labels/titles
  labs(
    y = "Percentage Survival",
    x = "Temperature in degrees Celcius",
    title = "Survival curve - GLM"
  ) +
  # reduce padding on edges of figure and format axes
  scale_y_continuous(label = scales::percent, 
                     breaks = seq(0, 1, by = 0.1),
                     expand = c(0.02, 0)) +
  scale_x_continuous(breaks = seq(22, 34, by = 3), 
                     expand = c(0.04, 0))

# compare models
anova(GLM_bi1, GLM_bi4, test = "Chi")

# all pops
ggplot(therm_surv_LF, aes(temp, surv, col = popf)) +
  #geom_point(aes(col = motu)) +
  geom_smooth(method = "glm", se = F,
              method.args = list(family = "quasibinomial"), linetype = "dashed") +
  coord_cartesian(xlim = c(22, 34)) +
  # update figure labels/titles
  labs(
    y = "Percentage Survival",
    x = "Temperature in degrees Celcius",
    title = "Survival curve - GLM"
    ) +
  # reduce padding on edges of figure and format axes
  scale_y_continuous(label = scales::percent, 
                     breaks = seq(0, 1, by = 0.1),
                     expand = c(0.02, 0)) +
  scale_x_continuous(breaks = seq(22, 34, by = 3), 
                     expand = c(0.04, 0))

## all pops with standard error
ggplot(therm_surv_LF, aes(temp, surv, col = popf)) +
  #geom_point(aes(col = motu)) +
  geom_smooth(method = "glm", se = T,
              method.args = list(family = "quasibinomial"), linetype = "dashed") +
  coord_cartesian(xlim = c(22, 34)) +
  # update figure labels/titles
  labs(
    y = "Percentage Survival",
    x = "Temperature in degrees Celcius",
    title = "Survival curve - GLM"
  ) +
  # reduce padding on edges of figure and format axes
  scale_y_continuous(label = scales::percent, 
                     breaks = seq(0, 1, by = 0.1),
                     expand = c(0.02, 0)) +
  scale_x_continuous(breaks = seq(22, 34, by = 3), 
                     expand = c(0.04, 0))

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
