library(nortest)
library(lmtest)
library(car)
library(fpp)
library(plyr)
library(coin)
library(PMCMRplus)
library(mblm)

RiskFire = read.csv("Focos_2022-01-01_2022-05-31.csv")
View(RiskFire)

#Ajust Dataset
RiskFire$pais <- NULL
RiskFire$estado <- NULL
RiskFire$municipio <- NULL
RiskFire$bioma <- NULL
RiskFire$latitude <- NULL
RiskFire$longitude <- NULL
RiskFire$datahora <- NULL
RiskFire$satelite <- NULL
RiskFire <-RiskFire[!(RiskFire$riscofogo==-999 | is.na(RiskFire$riscofogo) | is.na(RiskFire$frp)),]

RiskFire <- rename(RiskFire, c("riscofogo"="riskfire", "diasemchuva"="daywithoutrain", "precipitacao"="precipitation")) 

RiskFire$RiskFireLevel[RiskFire$riskfire < 0.15] <- "Minimum" 
RiskFire$RiskFireLevel[RiskFire$riskfire >= 0.15 & RiskFire$riskfire < 0.4] <- "Low"
RiskFire$RiskFireLevel[RiskFire$riskfire >= 0.4 & RiskFire$riskfire < 0.7] <- "Medium"
RiskFire$RiskFireLevel[RiskFire$riskfire >= 0.7 & RiskFire$riskfire < 0.95] <- "High" 
RiskFire$RiskFireLevel[RiskFire$riskfire >= 0.95] <- "Critical"
RiskFire$RiskFireLevel=factor(RiskFire$RiskFireLevel, levels=c("Minimum", "Low", "Medium", "High", "Critical"))

summary(RiskFire)

#boxplot(RiskFire$daywithoutrain~RiskFire$riskfire, main="Comparison: Risk of Fire x Days without rain", col=c("red", "blue"))
#boxplot(RiskFire$precipitation~RiskFire$riskfire, main="Comparison: Fire Risk x Precipitation", col=c("red", "blue"))
#boxplot(RiskFire$precipitation~RiskFire$daywithoutrain, main="Comparison: Days without rain x Precipitation", col=c("red", "blue"))
#boxplot(RiskFire$frp~RiskFire$riskfire, main="Comparison: Fire Risk x FRP", col=c("red", "blue"))
#boxplot(RiskFire$riskfire~RiskFire$RiskFireLevel, main="Comparison: Fire Risk x Fire Risk Scale");

#Mean by level
(mean_daywithoutrain <- tapply(RiskFire$daywithoutrain, RiskFire$RiskFireLevel, mean))
(mean_precipitation <- tapply(RiskFire$precipitation, RiskFire$RiskFireLevel, mean))
(mean_frp <- tapply(RiskFire$frp, RiskFire$RiskFireLevel, mean))

#Standard Deviation by level
(sd_daywithoutrain <- tapply(RiskFire$daywithoutrain, RiskFire$RiskFireLevel, sd))
(sd_precipitation <- tapply(RiskFire$precipitation, RiskFire$RiskFireLevel, sd))
(sd_frp <- tapply(RiskFire$frp, RiskFire$RiskFireLevel, sd))

boxplot(RiskFire$daywithoutrain~RiskFire$RiskFireLevel, main="Comparative : Risk Levels Fire x Days without rain", col=c("red", "blue"))
boxplot(RiskFire$precipitation~RiskFire$RiskFireLevel, main="Comparison: Risk Levels Fire x Precipitation", col=c("red", "blue"))
boxplot(RiskFire$frp~RiskFire$RiskFireLevel, main="Comparison: Fire Risk Levels x FRP", col=c("red", "blue"))

boxplot(RiskFire$precipitation~RiskFire$daywithoutrain, main="Comparison: Days without rain x Precipitation", col=c("red", "blue"))
boxplot(RiskFire$riskfire~RiskFire$RiskFireLevel, main="Comparison: Risk Levels x Fire Risk");


#Test the null hypothesis that the data are normal
hist(RiskFire$riskfire)
hist(RiskFire$daywithoutrain)
hist(RiskFire$precipitation)
hist(RiskFire$frp)

#Anderson-Darlin normality test
ad.test(RiskFire$riskfire)
ad.test(RiskFire$daywithoutrain)
ad.test(RiskFire$precipitation)
ad.test(RiskFire$frp)

qqPlot(RiskFire$riskfire)
qqPlot(RiskFire$daywithoutrain)
qqPlot(RiskFire$precipitation)
qqPlot(RiskFire$frp)


#Do not apply, as this is not a normal distribution.
#m=aov(riskfire~daywithoutrain*precipitation*frp, data=RiskFire)
#anova(m)

#Anderson-Darlin normality test residuals
#ad.test(m$residuals)
#qqPlot(m$residuals)
# test for the residuals homoscedasticity
#bptest(m)
#summary.lm(m)

#-------------------------------------------------------------------------------
#Transform in log1p - disregard the zero value
logRF <- log1p(RiskFire$riskfire)
ad.test(logRF)
qqPlot(logRF)

#Transform in sqrt
sqrtRF = sqrt(RiskFire$riskfire)
ad.test(sqrtRF)
qqPlot(sqrtRF)

#Transform cubic
cubica <- function(x){
  return(x^(1/3))
}
raizCRF <- cubica(RiskFire$riskfire)
ad.test(raizCRF) 
qqPlot(raizCRF)

#Use Box-Cox
##Go through each row and determine if a value is zero
RF <- RiskFire$riskfire +1
lambda <- BoxCox.lambda (RF, method=c("loglik"), lower=-5, upper= 5)
lambda
BoxCoxRF <- BoxCox (RF, lambda)
ad.test(BoxCoxRF)
qqPlot(BoxCoxRF)

# Nonparametric equivalent of independent-samples t-test
kruskal.test(daywithoutrain ~ RiskFireLevel, data=RiskFire)
kruskal.test(precipitation ~RiskFireLevel, data=RiskFire)
kruskal.test(frp ~ RiskFireLevel, data=RiskFire)

#posthoc
kwAllPairsConoverTest(daywithoutrain ~ RiskFireLevel, data=RiskFire, p.adjust.method="holm")
kwAllPairsConoverTest(precipitation ~ RiskFireLevel, data=RiskFire, p.adjust.method="holm")
kwAllPairsConoverTest(frp ~ RiskFireLevel, data=RiskFire, p.adjust.method="holm")


#Correlation 
#Close to zero are considered independent
knitr::kable(cor(RiskFire[, c("riskfire","daywithoutrain","precipitation","frp")], method = "spearman"))

#---------------------------------------------------------------
#Prediction in riskfire
pred = lm(riskfire ~ daywithoutrain*precipitation*frp, data=RiskFire)
summary(pred)
plot(pred)

#multicolinearity
# large values are candidates to drop from the model
vif(pred)

