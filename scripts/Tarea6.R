#TAREA 6
#Presenta: Zaira Rosario Pérez Vázquez

rm(list=ls())
options(max.print = 10000)
setwd("C:/Users/zaira/Documents/Doctorado PCF/3. Cursos/2. Primavera 2021/3. Series de tiempo/3. Tareas/Tarea 6")


#Librerias
library(forecast)
library(tseries)
library(ggfortify)
library(car)
library(nortest)
library(fUnitRoots)
library(astsa)
library(ggpubr)
library(lmtest)
library(TSPred)


### 1. Conjunto de datos
#Datos de temperatura
#01/1985 a 12/2019
datos_tma<-read.csv("Temperatura_Serie.csv",header=FALSE)
datos_tma
tma_serie<-ts(datos_tma,frequency = 12,start = c(1985,1),names = "Temperatura")
tma_serie

### 2. Reconocimiento de la serie
summary(tma_serie)
sd(tma_serie)
autoplot(tma_serie,xlab = "Tiempo",ylab = "Temperatura media (°C)", main = "a)")
par(mfrow=c(1,3),cex=0.9)
hist(tma_serie,main = "b)",xlab = "Temperatura media (°C)",ylab = "Frecuencia")
plot(density(tma_serie), main = "c)",ylab = "Densidad")
qqPlot(tma_serie,dist="norm",ylab = "Temperatura media (°C)",main = "d)")
dev.off() 


### 3. Prueba de normalidad
jarque.bera.test(tma_serie) #La serie no sigue una distribución normal


### 4. Prueba de estacionariedad y raiz unitaria. 
acf2(tma_serie) #Serie estacionaria
adf.test(tma_serie) #Se rechaza Ho, con un pvalue=0.01
pp.test(tma_serie)


### 4. Componentes de la serie de tiempo
Componentes=decompose(tma_serie)
plot(Componentes, cex=0.9)
seasonplot(tma_serie, col=rainbow(12), year.labels=TRUE,cex=0.5)
#Grafica polar
ggseasonplot(tma_serie,polar=TRUE,main="Temperatura" )

### 5. Identificación de un modelo ARIMA

#Eliminar estacionalidad
Serie_tma_ajustada<-tma_serie-Componentes$seasonal
Fig1<-autoplot(Serie_tma_ajustada,xlab = "",ylab = "(°C)", main = "Tma-season")
#Aplicando la primer diferencia
Serie_tma_estacionaria<-diff(Serie_tma_ajustada,differences=1)
Fig2<-autoplot(Serie_tma_estacionaria,xlab = "Tiempo",ylab = "(°C)", main = "d=1")
ggarrange(Fig1,Fig2,ncol = 1, nrow = 2)


#Identificación de parametros
acf(Serie_tma_estacionaria)
pacf(Serie_tma_estacionaria)
adf.test(Serie_tma_estacionaria) #función del paquete tseries
acf2(Serie_tma_estacionaria)


### 6. Ajuste del modelo
#Partición de datos
tma_train = window(tma_serie, start = c(1985,1), end = c(2012,12)) 
plot(tma_train)
dim(as.matrix(tma_train))
tma_test = window(tma_serie, start = c(2013,1))
plot(tma_test)
dim(as.matrix(tma_test))

frequency(tma_train)
Arima1<- arima(tma_train, order=c(1,0,0),seasonal = list(order = c(1,1,0), period = 12),method="ML")
Arima2<- arima(tma_train, order=c(1,0,0),seasonal = list(order = c(0,1,1), period = 12),method="ML")
Arima3<- arima(tma_train, order=c(1,0,0),seasonal = list(order = c(1,1,1), period = 12),method="ML")
Arima4<- arima(tma_train, order=c(0,0,1),seasonal = list(order = c(1,1,0), period = 12),method="ML")
Arima5<- arima(tma_train, order=c(0,0,1),seasonal = list(order = c(0,1,1), period = 12),method="ML")
Arima6<- arima(tma_train, order=c(0,0,1),seasonal = list(order = c(1,1,1), period = 12),method="ML")
Arima7<- arima(tma_train, order=c(1,0,1),seasonal = list(order = c(1,1,0), period = 12),method="ML")
Arima8<- arima(tma_train, order=c(1,0,1),seasonal = list(order = c(0,1,1), period = 12),method="ML")
Arima9<- arima(tma_train, order=c(1,0,1),seasonal = list(order = c(1,1,1), period = 12),method="ML")
Arima10<- arima(tma_train, order=c(2,0,2),seasonal = list(order = c(1,1,1), period = 12),method="ML")
Arima11<- arima(tma_train, order=c(1,0,2),seasonal = list(order = c(2,1,2), period = 12),method="ML")
Arima12<- arima(tma_train, order=c(0,0,3),seasonal = list(order = c(2,1,2), period = 12),method="ML")
Arima13<- arima(tma_train, order=c(2,0,0),seasonal = list(order = c(2,1,1), period = 12),method="ML")
Arima14<- arima(tma_serie, order=c(3,0,0),seasonal = list(order = c(2,1,2), period = 12),method="ML")

AIC(Arima1,Arima2,Arima3,Arima4,Arima5,Arima6,Arima7,Arima8,Arima9,Arima10,Arima11,Arima12,Arima13,Arima14)
#<AIC=Arima11
BIC(Arima1,Arima2,Arima3,Arima4,Arima5,Arima6,Arima7,Arima8,Arima9,Arima10,Arima11,Arima12,Arima13,Arima14)
#<BIC=Arima8

coeftest(Arima8)
coeftest(Arima11)

confint(Arima8)
coeftest(Arima11)

### 7. Diagnóstico del modelo
#Arima8
ggtsdiag(Arima8)
Box.test(Arima8$residuals) # Test de Box-Pierce #Independencia
Box.test(Arima8$residuals, type="Ljung-Box") #Independencia
jarque.bera.test(Arima8$residuals) #Normalidad
par(mfrow=c(1,2),cex=0.7)
hist(Arima8$residuals, main="a)",xlab = "Residuales",ylab = "Frecuencia")
qqPlot(Arima8$residuals,distribution = "norm",main = "b)",ylab = "Residuales")
dev.off()
#Arima11
ggtsdiag(Arima11)
Box.test(Arima11$residuals) # Test de Box-Pierce #Independencia
Box.test(Arima11$residuals, type="Ljung-Box") #Independencia
jarque.bera.test(Arima11$residuals) #Normalidad
par(mfrow=c(1,2),cex=0.7)
hist(Arima11$residuals, main="a)",xlab = "Residuales",ylab = "Frecuencia")
qqPlot(Arima11$residuals,distribution = "norm",main = "b)",ylab = "Residuales")
dev.off()

### 8. Prediccion del modelo
##Predicción (con el mejor modelo: ARIMA11)
#Tres periodos
#forecast::forecast()

pred1_arima11<-forecast::forecast(Arima11,h=12, level=c(99.5))
Fig3<-autoplot(pred1_arima11, main = "1 año")

pred2_arima11<-forecast::forecast(Arima11,h=24, level=c(99.5))
Fig4<-autoplot(pred2_arima11, main = "2 años")

pred3_arima11<-forecast::forecast(Arima11,h=36, level=c(99.5))
Fig5<-autoplot(pred3_arima11, main = "3 años")

ggarrange(Fig3,Fig4,Fig5,ncol = 1, nrow = 3)

### 9. Validación de la predicción
forecast::accuracy(pred1_arima11, tma_test)
forecast::accuracy(pred2_arima11, tma_test)
forecast::accuracy(pred3_arima11, tma_test)
plotarimapred(tma_test, Arima11, xlim=c(2013, 2019), range.percent = 0.05)

