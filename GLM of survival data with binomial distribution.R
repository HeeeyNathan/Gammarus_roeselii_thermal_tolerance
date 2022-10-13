## set WD
# setwd() # set the working directory

## Load require packages
library(ggplot2)
library(gridExtra)
library(grid)
library(gridtext)
library(tidyverse)

## load data
therm_surv_LF <- read.csv("Survival_own_data.csv", header = T, sep = ",", row.names = c(1) ,stringsAsFactors = FALSE)
therm_surv_LF$popf <- as.factor(therm_surv_LF$pop)
therm_surv_LF$motuf <- as.factor(therm_surv_LF$motu)
therm_surv_LF$motu_numf <- as.factor(therm_surv_LF$motu_num)
therm_surv_LF$motusf <- as.factor(therm_surv_LF$motus)
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
# compare models
anova(GLM_bi1, GLM_bi2, GLM_bi3)
AIC(GLM_bi1, GLM_bi2, GLM_bi3)
BIC(GLM_bi1, GLM_bi2, GLM_bi3)

## model selection
# deviance table of best model
anova(GLM_bi1, test = "Chisq")
# using drop1 function is easier
drop1(GLM_bi1, test = "Chi")

## model validation
op <- par(mfrow = c(2,2))
plot(GLM_bi1, which = 1:4)
par(op)

## Plotting the model
# using the predict command
MyData <- data.frame(temp = 
                       seq(from = min(therm_surv_LF$temp),
                           to = max(therm_surv_LF$temp), by = 0.1))
pred <- predict(GLM_bi1, newdata = MyData, type = "response")
plot(therm_surv_LF$temp, therm_surv_LF$surv,
     xlab = "Temperature (degrees Celcius)",
     ylab = "Probability of survival (%)")
lines(MyData$temp, pred, lty = 1)
# base plot
plot(therm_surv_LF$temp, fitted(glm(surv ~ temp, family = "binomial"(link = "logit"), data = therm_surv_LF)))
GLM_bi1
temp1 <- seq(22,34,0.1)
newtemp <- data.frame(temp=temp1)
predicted.probability <- predict(GLM_bi1, + newtemp, type = "resp")
plot(predicted.probability ~ temp1, type = "l")

# ggplot
p1 <- ggplot(therm_surv_LF, aes(temp, surv, col = "All MOTUs")) +
  #geom_point(aes(col = motu)) +
  geom_smooth(method = "glm", se = T,
              method.args = list(family = "binomial"(link = "logit")), linetype = "dashed") +
  coord_cartesian(xlim = c(22, 34)) +
  scale_color_manual(values = "black") +
  theme(legend.position = "none") + 
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
p1

## GLM with binomial distribution - temp + motu
GLM_bi4 <- glm(surv ~ temp + motuf, data = therm_surv_LF, family = "binomial"(link = "logit")) # assumes ~ equal number of 1s and 0s
summary(GLM_bi4, corr = T) # in a Bernoulli GLM (binomial GLM) because the data are code as 0 and 1, overdispersion cannot occur.
GLM_bi5 <- glm(surv ~ temp + motuf, data = therm_surv_LF, family = "binomial"(link = "probit")) # assumes ~ equal number of 1s and 0s
summary(GLM_bi5, corr = T) # in a Bernoulli GLM (binomial GLM) because the data are code as 0 and 1, overdispersion cannot occur.
GLM_bi6 <- glm(surv ~ temp + motuf, data = therm_surv_LF, family = "binomial"(link = "cloglog")) # used for considerable more 1s than 0s or vice versa
summary(GLM_bi6, corr = T) # in a Bernoulli GLM (binomial GLM) because the data are code as 0 and 1, overdispersion cannot occur.
# compare models
anova(GLM_bi4, GLM_bi5, GLM_bi6)
AIC(GLM_bi4, GLM_bi5, GLM_bi6)
BIC(GLM_bi4, GLM_bi5, GLM_bi6)

## posthoc testing to determine factor comparisons
library(multcomp)
summary(glht(GLM_bi4, mcp(motuf = "Tukey"), alternative = "two.sided"))

## model selection
# deviance table of best model
anova(GLM_bi4, test = "Chisq")
# using drop1 function is easier
drop1(GLM_bi4, test = "Chi")

## model validation
op <- par(mfrow = c(2,2))
plot(GLM_bi4, which = 1:4)
par(op)

## Presentation as odds-ratio estimates
exp(cbind(OR = coef(GLM_bi4), confint(GLM_bi4)))
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
p2 <- ggplot(therm_surv_LF, aes(temp, surv, col = motuf)) +
  geom_smooth(method = "glm", se = T,
              method.args = list(family = "binomial"(link = "logit")), linetype = "dashed") +
  coord_cartesian(xlim = c(22, 34)) +
  scale_color_manual(values = c("red", "green", "blue")) +
  theme(legend.position = "none") + 
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
p2

# compare models
anova(GLM_bi1, GLM_bi4, test = "Chi") # fuller model is better than model with just temp

## Create aestically pleasing plots
p3 <- ggplot(therm_surv_LF, aes(temp, surv)) + 
  theme_bw() +
  # Add reference lines
  # geom_vline(xintercept = c(25, 28, 31), linetype = "dashed", color = "gray50", size = 0.5) +
  # geom_hline(yintercept = c(0.5), linetype = "dashed", color = "gray50", size = 0.5) +
  geom_segment(x = 21, y = 0.5, xend = 30.24, yend = 0.5, color = "#464646", linetype = "dashed") +
  geom_segment(x = 30.24, y = 0.5, xend = 30.24, yend = -0.5, color = "#464646", linetype = "dashed") + # motu L
  geom_segment(x = 30.24, y = 0.5, xend = 31, yend = 0.5, color = "#AB82FF", linetype = "dashed") +
  geom_segment(x = 31, y = 0.5, xend = 31, yend = -0.5, color = "#AB82FF", linetype = "dashed") + # motu A
  geom_segment(x = 31, y = 0.5, xend = 31.25, yend = 0.5, color = "#9ACD32", linetype = "dashed") +
  geom_segment(x = 31.25, y = 0.5, xend = 31.25, yend = -0.5, color = "#9ACD32", linetype = "dashed") + # motu G
  # add smoothers
  geom_smooth(data = therm_surv_LF, aes(temp, surv, col = motuf, fill = motuf), method = "glm", se = T,
              method.args = list(family = "binomial"(link = "logit")), linetype = "dashed") + 
  # geom_smooth(data = therm_surv_LF, aes(temp, surv, col = "All MOTUs"), method = "glm", se = T,
  #             method.args = list(family = "binomial"(link = "logit")), linetype = "dashed", col = "black") +
  coord_cartesian(xlim = c(22, 34)) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_fill_manual(values = c("#AB82FF", "#9ACD32", "#464646")) +
  scale_color_manual(values = c("#AB82FF", "#9ACD32", "#464646")) +
  theme(legend.position = "bottom") + 
  # add text at LC50
  geom_label(aes(label = "LT50 at 30.24\u00B0C"), x = 30.24, y = 0.48, colour = "#464646", fontface = "bold", size = 3) +
  geom_label(aes(label = "LT50 at 31\u00B0C"), x = 31, y = 0.52, colour = "#AB82FF", fontface = "bold", size = 3) +
  geom_label(aes(label = "LT50 at 31.25\u00B0C"), x = 31.25, y = 0.44, colour = "#9ACD32", fontface = "bold", size = 3) +
  # update figure labels/titles
  labs(
    y = "Probability of Survival (%)",
    x = "Temperature (\u00B0C)",
    title = "Survival curve (GLM) - Split MOTUs"
  ) +
  # reduce padding on edges of figure and format axes
  scale_y_continuous(label = scales::percent, 
                     breaks = seq(0, 1, by = 0.05),
                     expand = c(0.02, 0)) +
  scale_x_continuous(breaks = seq(22, 34, by = 1), 
                     expand = c(0.04, 0))
p3

p4 <- ggplot(therm_surv_LF, aes(temp, surv)) + 
  theme_bw() +
  # Add reference lines
  geom_segment(x = 21, y = 0.5, xend = 30.83, yend = 0.5, color = "#b22222", linetype = "dashed") +
  geom_segment(x = 30.83, y = 0.5, xend = 30.83, yend = -0.5, color = "#b22222", linetype = "dashed") + # motu L
  # add smoothers
  geom_smooth(data = therm_surv_LF, aes(temp, surv, col = motusf, fill = motusf), method = "glm", se = T,
              method.args = list(family = "binomial"(link = "logit")), linetype = "dashed") +
  coord_cartesian(xlim = c(22, 34)) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_fill_manual(values = c("#b22222")) +
  scale_color_manual(values =  c("#b22222")) +
  theme(legend.position = "bottom") + 
  # add text at LC50
  geom_label(aes(label = "LT50 at 30.83\u00B0C"), x = 30.83, y = 0.52, colour = "#b22222", fontface = "bold", size = 3) +
  # update figure labels/titles
  labs(
    y = "Probability of Survival (%)",
    x = "Temperature (\u00B0C)",
    title = "Survival curve (GLM) - All MOTUs"
  ) +
  # reduce padding on edges of figure and format axes
  scale_y_continuous(label = scales::percent, 
                     breaks = seq(0, 1, by = 0.05),
                     expand = c(0.02, 0)) +
  scale_x_continuous(breaks = seq(22, 34, by = 1), 
                     expand = c(0.04, 0))
p4

## join the two plots
# Remove axis titles from all plots
p5 = list(p4, p3) %>% map(~.x + labs(x = NULL, y = NULL))
# create new axis names
yleft = richtext_grob("Probability of Survival (%)", rot = 90)
bottom = richtext_grob("Temperature (\u00B0C)")
# Lay out plots
grid.arrange(grobs = p5, ncol = 2, nrow = 1, left = yleft, bottom = bottom)

## save plots
# combined motus glm
svg(filename = "Combined_motus_glm.svg",
    width = 7, height = 7, pointsize = 12,
    onefile = FALSE, family = "sans", bg = "white")
p4
dev.off()
# seperate motus glm
svg(filename = "Seperate_motus_glm.svg",
    width = 7, height = 7, pointsize = 12,
    onefile = FALSE, family = "sans", bg = "white")
p3
dev.off()
# panel plot
svg(filename = "Panel_motus_glm.svg",
    width = 14, height = 7, pointsize = 12,
    onefile = FALSE, family = "sans", bg = "white")
grid.arrange(grobs = p5, ncol = 2, nrow = 1, left = yleft, bottom = bottom)
dev.off()
# updated panel plot
pdf("Panel_motus_glm.pdf", onefile = T, width = 14, height = 7)
grid.arrange(grobs = p5, ncol = 2, nrow = 1, left = yleft, bottom = bottom)
dev.off()
