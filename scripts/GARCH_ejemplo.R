setwd("C:/Users/zaira/Documents/Doctorado PCF/3. Cursos/2. Primavera 2021/3. Series de tiempo/3. Tareas/Tarea 13")
rm(list=ls())
options(max.print = 10000)

rendHD<- read.csv("rend_HD.csv",header = TRUE)
View(rendHD)
resHD<-rendHD$r_HOME_DEPOT

par(mfrow=c(1,2),mar=c(5,4,3,1),cex=0.6)
plot.ts(resHD)
hist(resHD, main="", breaks=20, freq=FALSE, col="grey")
dev.off()

par(mfrow=c(2,2),mar=c(5,4,3,1),cex=0.6)
acf(resHD)
pacf(resHD)
acf(resHD^2)
pacf(resHD^2)


###ARCH

library("dynlm")
resHD_mean<-dynlm(resHD~1)
summary(resHD_mean)
ehat2<-ts(resid(resHD_mean)^2)
HD_arch<-dynlm(ehat2~L(ehat2))
summary(HD_arch)

##Estimación de parametros
library(broom)
T <- 500
q <- length(coef(HD_arch))-1
Rsq <- glance(HD_arch)[[1]]
LM <- (T-q)*Rsq
alpha <- 0.05
Chicr <- qchisq(1-alpha, q)
#Otra forma
library(FinTS)
(HD_archtest<-ArchTest(resHD,lags = 1,demean = T))


#ARCH (0,1)
library(tseries)
HD.arch_m2 <- garch(resHD,c(0,1))
summary(HD.arch_m2)

hhat <- ts(2*HD.arch_m2$fitted.values[-1,1]^2)
plot.ts(hhat)

##Modelo GARCH
library(rugarch)
garch_Spec <- ugarchspec(variance.model=list(model="sGARCH",
                      garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0)),
                      distribution.model="std")
garchFit <- ugarchfit(spec=garch_Spec, data=resHD)
(coef(garchFit))
rhat <- garchFit@fit$fitted.values

par(mfrow=c(1,2),mar=c(5,4,3,1),cex=0.6)
plot.ts(rhat)
hhat2 <- ts(garchFit@fit$sigma^2)
plot.ts(hhat2)
