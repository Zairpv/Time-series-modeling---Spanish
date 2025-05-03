#TAREA 3
#Presenta: Zaira Rosario Pérez Vázquez

rm(list=ls())
options(max.print = 10000)

##1. Librerias
library(TSA)
library(tseries)
library(ggfortify)
library(ggpubr)

###2. Simulación de un proceso estacionario

set.seed(98765)
?arima.sim

##Simulación AR(1)
Sim1_AR<-arima.sim(list(order=c(1,0,0),ar=0.8),n=1000)
Sim2_AR<-arima.sim(list(order=c(1,0,0),ar=0.9),n=100)
#Gráfica de las simulaciones AR(1)
Fig1_AR_S1<-autoplot(Sim1_AR, ts.colour = "dark blue",xlab = "Tiempo",ylab = "x")
Fig1_AR_S2<-autoplot(Sim2_AR, ts.colour = "brown",xlab = "Tiempo",ylab = "x")
ggarrange(Fig1_AR_S1,Fig1_AR_S2,ncol = 1, nrow = 2,labels = c("A","B"))

##Simulación MA(1)
Sim3_MA<-arima.sim(list(order=c(0,0,1),ma=0.7),n=1000)
Sim4_MA<-arima.sim(list(order=c(0,0,1),ma=-0.5),n=100) 
#Gráfica de las simulaciones MA(1)
Fig2_MA_S1<-autoplot(Sim3_MA, ts.colour = "dark blue",xlab = "Tiempo",ylab = "x")
Fig2_MA_S2<-autoplot(Sim4_MA, ts.colour = "brown",xlab = "Tiempo",ylab = "x")
ggarrange(Fig2_MA_S1,Fig2_MA_S2,ncol = 1, nrow = 2,labels = c("A","B"))

##Simulación ARMA(1,0,1)
Sim5_ARMA<-arima.sim(list(order=c(1,0,1),ar=0.8,ma=0.7),n=1000)
Sim6_ARMA<-arima.sim(list(order=c(1,0,1),ar=0.9,ma=-0.5),n=100)
#Gráfica de las simulaciones MA(1,0,1)
Fig3_ARMA_S1<-autoplot(Sim5_ARMA, ts.colour = "dark blue",xlab = "Tiempo",ylab = "x")
Fig3_ARMA_S2<-autoplot(Sim6_ARMA, ts.colour = "brown",xlab = "Tiempo",ylab = "x")
ggarrange(Fig3_ARMA_S1,Fig3_ARMA_S2,ncol = 1, nrow = 2,labels = c("A","B"))


###3. Identificación del modelo

##Autorregresivo AR(1)
#ACF y PACF
par(mfrow=c(2,2),cex=0.6)
#Sim1_AR
acf(Sim1_AR,main=(expression(ACF: AR(1)~~~phi==0.8)),ylim=c(-0.1,1.1)) 
pacf(Sim1_AR,main=(expression(PACF: AR(1)~~~phi==0.8)),ylim=c(-0.1,1.1))
#Sim2_AR
acf(Sim2_AR,main=(expression(ACF: AR(1)~~~phi==0.9)),ylim=c(-0.4,1.1)) 
pacf(Sim2_AR,main=(expression(PACF: AR(1)~~~phi==0.9)),ylim=c(-0.1,1.1))
dev.off()

##Promedios móviles AM(1)
par(mfrow=c(2,2),cex=0.6)
#Sim3_MA
acf(Sim3_MA,main=(expression(ACF: MA(1)~~~theta==0.7))) 
pacf(Sim3_MA,main=(expression(PACF: MA(1)~~~theta==0.7)))
#Sim4_MA
acf(Sim4_MA,main=(expression(ACF: MA(1)~~~theta==-0.5))) 
pacf(Sim4_MA,main=(expression(PACF: MA(1)~~~theta==-0.5)))
dev.off()

##ARMA (1,0,1)
par(mfrow=c(2,2),cex=0.6)
#Sim5_ARMA
acf(Sim5_ARMA,main=(expression(ACF: ARMA(1,0,1)~~~phi==0.8~~~theta==0.7))) 
pacf(Sim5_ARMA,main=(expression(PACF: ARMA(1,0,1)~~~phi==0.8~~~theta==0.7)))
#Sim6_ARMA
acf(Sim6_ARMA,main=(expression(ACF: ARMA(1,0,1)~~~phi==0.9~~~theta==-0.5))) 
pacf(Sim6_ARMA,main=(expression(PACF: ARMA(1,0,1)~~~phi==0.9~~~theta==-0.5)))
dev.off()

###4. Ajuste del modelo
?arima
##Autorregresivos AR(1)
#Sim1_AR
S1_AR<-arima(Sim1_AR,order=c(1,0,0)) #Con intercepto
S1_AR
S2_AR<-arima(Sim1_AR,order=c(1,0,0),include.mean = FALSE) #Sin intercepto
S2_AR
#Sim2_AR
S3_AR<-arima(Sim2_AR,order=c(1,0,0))
S3_AR
S4_AR<-arima(Sim2_AR,order=c(1,0,0),include.mean = FALSE)
S4_AR

##Promedios móviles MA(1)
#Sim3_MA
S1_MA<-arima(Sim3_MA,order=c(0,0,1))
S1_MA
S2_MA<-arima(Sim3_MA,order=c(0,0,1),include.mean = FALSE)
S2_MA
#Sim4_MA
S3_MA<-arima(Sim4_MA,order=c(0,0,1))
S3_MA
S4_MA<-arima(Sim4_MA,order=c(0,0,1),include.mean = FALSE)
S4_MA

##ARMA (1,1)
#Sim5_ARMA
S1_ARMA<-arima(Sim5_ARMA,order=c(1,0,1))
S1_ARMA
S2_ARMA<-arima(Sim5_ARMA,order=c(1,0,1),include.mean = FALSE)
S2_ARMA
#Sim6_ARMA
S3_ARMA<-arima(Sim6_ARMA,order=c(1,0,1))
S3_ARMA
S4_ARMA<-arima(Sim6_ARMA,order=c(1,0,1),include.mean = FALSE)
S4_ARMA

###5. Diagnostico de los modelos
##Autorregresivo
#Sim1_AR
AIC(S1_AR)
BIC(S1_AR)
AIC(S2_AR)
BIC(S2_AR)

#Sim2_AR
AIC(S3_AR)
BIC(S3_AR)
AIC(S4_AR)
BIC(S4_AR)

##Promedios moviles
#Sim3_MA
AIC(S1_MA)
BIC(S1_MA)
AIC(S2_MA)
BIC(S2_MA)
#Sim4_MA
AIC(S3_MA)
BIC(S3_MA)
AIC(S4_MA)
BIC(S4_MA)

##ARMA 
#Sim5_ARMA
AIC(S1_ARMA)
BIC(S1_ARMA)
AIC(S2_ARMA)
BIC(S2_ARMA)
#Sim6_ARMA
AIC(S3_ARMA)
BIC(S3_ARMA)
AIC(S4_ARMA)
BIC(S4_ARMA)

#Por lo anterior, se grafican los residuales de los 
#modelos ajustados sin intercepto para cada una de las simulaciones

##Autorregresivos
ggtsdiag(S2_AR)
ggtsdiag(S4_AR)
##Promedios móviles
ggtsdiag(S2_MA)
ggtsdiag(S4_MA)
##ARMA
ggtsdiag(S2_ARMA)
ggtsdiag(S4_ARMA)
