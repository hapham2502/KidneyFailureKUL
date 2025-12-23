library(dplyr)
library(MASS)
library(multcomp)

days <- c(
  0, 2, 1, 3, 0, 2, 0, 5, 6, 8, 2, 4, 7, 12, 15, 4, 3, 1, 5, 20, 15, 10, 8, 5, 25, 16, 7, 30, 3, 27, 0, 1, 1, 0, 4, 2, 7, 4, 0, 3, 5, 3, 2, 0, 1, 1, 3, 6, 7, 9, 10, 8, 12, 3, 7, 15, 4, 9, 6, 1
)

duration <- c(rep(1, 30), rep(2, 30)) 
weight   <- rep(c(rep(1, 10), rep(2, 10), rep(3, 10)), 2)

# Combine into a data frame
KidneyFailure <- data.frame(
  Days = days,
  Duration = factor(duration, labels = c("Short", "Long")),
  WeightGain = factor(weight, labels = c("Slight", "Moderate", "Substantial"))
)

# Checking cell balance
table(KidneyFailure$Duration, KidneyFailure$WeightGain) 
# 10 in every cell --> balanced

summary.stats <- KidneyFailure %>%
  group_by(Duration, WeightGain) %>%
  summarise(
    Mean = mean(Days),
    Var = var(Days),
    SD = sd(Days)
  )

print(summary.stats)
#SD for short/slight = 2.7, for short/substantial = 9.7 --> red flag for heteroscedasticity

# boxplot (HA)
KidneyFailure$Group <- interaction(KidneyFailure$Duration, 
                                   KidneyFailure$WeightGain,
                                   sep=":")
KidneyFailure$Group <- factor(KidneyFailure$Group,
                              levels=c("Short:Slight", "Short:Moderate", 
                                       "Short:Substantial", "Long:Slight",
                                       "Long:Moderate", "Long:Substantial"))

par(mar = c(8, 4, 4, 2) + 0.1)  # Bottom margin = 8 lines

boxplot(Days ~ Group, data = KidneyFailure,
        main = "Hospital Days by Duration & Weight Gain",
        ylab = "Days Hospitalized",
        xlab = "",
        las = 2,  # Makes x-axis labels vertical
        col = rep(c("lightblue", "lightgreen"), each = 3),
        border = "black",
        ylim = c(0, max(KidneyFailure$Days) * 1.1))

# means per factor (HA)
plot.design(Days ~ Duration + WeightGain, data = KidneyFailure,
            main = "Mean Hospitalization Days by Factors",
            xlab = "Factors",
            ylab = "Mean Days",
            col = "black",
            pch = 16,
            cex = 1.2)
# interaction plot (HA)
interaction.plot(x.factor = KidneyFailure$WeightGain,  
                 trace.factor = KidneyFailure$Duration, 
                 response = KidneyFailure$Days,         
                 fun = mean,                            
                 xlab = "Weight Gain",                  
                 ylab = "Mean Hospitalization Days",    
                 main = "Interaction Plot: Duration × Weight Gain",
                 col = c("blue", "red"),               
                 lwd = 2,                               
                 trace.label = "Duration",              
                 type = "b")                           

# ANALYSIS (HA)
# 1. Test interaction 
fit_full <- aov(Days ~ Duration * WeightGain, data = KidneyFailure)
summary(fit_full)

# 2. Test main effects (additive model)
fit_additive <- aov(Days ~ Duration + WeightGain, data = KidneyFailure)
summary(fit_additive)

# MODEL AND ASSUMPTIONS

model.full <- lm(Days ~ Duration * WeightGain, data = KidneyFailure)
model.full

par(mfrow = c(2, 2))
plot(model.full) 
par(mfrow = c(1, 1))
#Residuals plot shows funnel shape (bigger spread in residuals for higher values)

shapiro.test(residuals(model.full))
#significant --> abnormal distribution of residuals

library(car)
# Homogeneity test (HA)
levene_result <- leveneTest(Days ~ Duration * WeightGain, 
                            data = KidneyFailure, center = mean)

print(levene_result)
#p-value < 0.001 indicates unequal variances across treatment duration and
# weight gain group.

# TRANSFORMATIONS

#box cox to find best transform
#use days + 1 since log(0) won't work
model.plus.1 <- lm(Days + 1 ~ Duration * WeightGain, data = KidneyFailure)
bc <- boxcox(model.plus.1, lambda = seq(-2, 2, 0.1))

# Find best Lambda
lambda.max <- max(bc$y)
best.lambda <- bc$x[bc$y == lambda.max]
print(best.lambda)

# Since lambda approaches 0 we use log-transformation
KidneyFailure$Days.Log <- log(KidneyFailure$Days + 1)

# TRANSFORMED MODEL

fit.full <- aov(Days.Log ~ Duration * WeightGain, data = KidneyFailure)

par(mfrow = c(2, 2))
plot(fit.full) 
par(mfrow = c(1, 1))

shapiro.test(residuals(fit.full))
#transformation seems to have fixed homoscedastisticy

summary(fit.full)
# interaction is not significant --> leave it out of model

fit.main <- aov(Days.Log ~ Duration + WeightGain, data = KidneyFailure)
summary(fit.main)

# POST HOC ANALYSIS (Weight gain)

#Tukey for WeightGain
tukey.weight <- glht(fit.main, linfct = mcp(WeightGain = "Tukey"))
tukey.ci <- confint(tukey.weight)
summary(tukey.weight)
plot(tukey.weight)

#back transofrmation for Weight Gain
weightgain.bt <- data.frame(
  Comparison = rownames(tukey.ci$confint),
  Log.Diff   = tukey.ci$confint[, "Estimate"],
  Log.Lwr    = tukey.ci$confint[, "lwr"],
  Log.Upr    = tukey.ci$confint[, "upr"]
)

#add p-values (HA)
tukey.summary <- summary(tukey.weight)
weightgain.bt$p.value <- round(tukey.summary$test$pvalues, 4)
weightgain.bt$p.value
## conclusion: p < 0.05 for all three comparisons => all weight groups differ
## significantly from each other. Every increase in weight gain leads to 
## significantly longer hospitalized days. 

# the ratio is e^(log_diff)
weightgain.bt$Ratio <- exp(weightgain.bt$Log.Diff)
weightgain.bt$Ratio.Lwr <- exp(weightgain.bt$Log.Lwr)
weightgain.bt$Ratio.Upr <- exp(weightgain.bt$Log.Upr)

# Percentage change (Formula: (Ratio - 1) * 100)
weightgain.bt$Pct.Change <- (weightgain.bt$Ratio - 1) * 100
weightgain.bt$Pct.Lwr    <- (weightgain.bt$Ratio.Lwr - 1) * 100
weightgain.bt$Pct.Upr    <- (weightgain.bt$Ratio.Upr - 1) * 100

weightgain.bt <- weightgain.bt %>%
  dplyr::select(Comparison, Ratio, Pct.Change, Pct.Lwr, Pct.Upr) %>%
  mutate(across(where(is.numeric), round, 1))

print(weightgain.bt)


#patients with a substantial weight gain have a significanly higher stay.
#3.7x times higher compared to slight weight (267%)
#moderate --> sligt 1.9 (88%)
#substantial --> moderate 2.0 (95%)

# ANALYSIS OF DURATION

#Duration has only two groups, so we already know that the effect is longer stay --> shorter hospitalization

#confidence intervals
ci.log <- confint(fit.main, parm="DurationLong", level = .95)

#coefficient
coef.log <- coef(fit.main)["DurationLong"]
#Back transformation for Duration
duration.bt <- data.frame(
  Comparison = "Long - Short" ,
  Ratio      = exp(coef(fit.main)["DurationLong"]),
  Pct.Change = (exp(coef.log) - 1) * 100,
  Pct.Lwr    = (exp(ci.log[1])-1) * 100 ,
  Pct.Upr    = (exp(ci.log[2])-1) * 100
)

duration.bt <- duration.bt %>%
  dplyr::select(Comparison, Ratio, Pct.Change, Pct.Lwr, Pct.Upr) %>%
  mutate(across(where(is.numeric), round, 1))

print(duration.bt)
#patients with a longer treatment had a 32,6% shorter stay in the hospital

# 7. Visualization & Conclusion

interaction.plot(x.factor = KidneyFailure$WeightGain, 
                 trace.factor = KidneyFailure$Duration, 
                 response = KidneyFailure$Days.Log, 
                 fun = mean, 
                 type = "b", 
                 legend = TRUE, 
                 ylab = "Mean Log(Days + 1)", 
                 xlab = "Weight Gain",
                 pch = c(1, 19),
                 col = c("black", "red"))

#Weight gain between treatment has a very strong impact on days in the hospital. 
#The extremer the weight gain, leads to longer hospital stays
#The time of dialysis impacts hospitalization, regardless of weight gain. 
#No interaction effect --> effects are additive
#Have to keep in mind log-transformation, the difference in *actual* amount of days between
#weight gain classes, but also withing for example substantial weight gain between
#short and long treatment within substantial weight gain is huge