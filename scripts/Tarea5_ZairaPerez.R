#TAREA 5
#Presenta: Zaira Rosario Pérez Vázquez

rm(list=ls())
options(max.print = 10000)
setwd("C:/Users/zaira/Documents/Doctorado PCF/3. Cursos/2. Primavera 2021/3. Series de tiempo/3. Tareas/Tarea 5")

#Lectura de datos
Datos=scan('https://bit.ly/3afYHIK')
View(Datos)

#Librerias
library(tseries)
library(forecast)
library(ggplot2)
library(fUnitRoots)
library(ggpubr)
library(astsa)
library(lmtest)
library(ggfortify)
library(car)

### 1. Reconocimiento de la serie
Serie_PM = ts(Datos, start = c(2006,1), frequency = 12,names = "PM 2.5")
Serie_PM
autoplot(Serie_PM,xlab = "Tiempo",ylab = "Particulas suspendidas", main = "a)")
summary(Serie_PM)
var(Serie_PM)
sd(Serie_PM)
par(mfrow=c(1,3),cex=0.9)
hist(Serie_PM,main = "b)",xlab = "Partículas Suspendidas",ylab = "Frecuencia")
plot(density(Serie_PM), main = "c)",ylab = "Densidad")
qqPlot(Serie_PM,dist="norm",ylab = "Partículas Suspendidas",main = "d)")
dev.off() #800x320px



### 2. Prueba de normalidad
jarque.bera.test(Serie_PM) #La serie sigue una distribución normal



### 3. Componentes de la serie de tiempo
Componentes=decompose(Serie_PM)
plot(Componentes)



### 4. Prueba de estacionariedad
adfTest(Serie_PM) #Serie no estacionaria, con raiz unitaria
adfTest(log(Serie_PM))
acf2(Serie_PM) #Decaimiento lento
acf2(log(Serie_PM))

### 5. Prueba de raiz unitaria
#Ajuste de un modelo AR a la serie diferenciada
?adfTest()
plot(diff(log(Serie_PM)))

#Nota: la serie diferenciada no presenta tendencia
acf2(diff(log(Serie_PM)))
mod1=ar(diff(log(Serie_PM)),method="mle")
mod1$order #Fue 12

#Prueba ADF ajustada a 12 lags
adfTest(log(Serie_PM),lags=12,type=c("c")) 

#¿Hay pendiente? ¿Cuánto vale la constante?
t=seq(2:length(log(Serie_PM)))
mod2=lm(diff(log(Serie_PM))~t)
summary(mod2)
#Nota: No es significativo el intercepto ni la pendiente
#entonces no hay deriva

#Prueba ADF ajustada sin constante
adfTest(log(Serie_PM),lags=12,type=c("nc")) 
#Nota: la serie tiene  raiz unitaria
#Requiere de una diferenciación
#Comprobación
ndiffs(Serie_PM)
ndiffs(log(Serie_PM))

#Prueba con unitrootTest
unitrootTest(log(Serie_PM),lags=1,type="nc")


### 5. Identificación de un modelo ARIMA

#Eliminar no estacionariedad
Serie_PM_diferenciada<-diff(Serie_PM,differences=1)
Fig1<-autoplot(Serie_PM_diferenciada,xlab = "Tiempo",ylab = "Particulas suspendidas", main = "a)")

#Eliminar estacionalidad
Serie_PM_ajustada<-Serie_PM-Componentes$seasonal
Fig2<-autoplot(Serie_PM_ajustada,xlab = "Tiempo",ylab = "Particulas suspendidas", main = "b)")

Serie_PM_estacionaria<-diff(Serie_PM_ajustada,differences=1)
Fig3<-autoplot(Serie_PM_estacionaria,xlab = "Tiempo",ylab = "Particulas suspendidas", main = "c)")

ggarrange(Fig1,Fig2,Fig3,ncol = 1, nrow = 3)

#Identificación de los valores p, d, q. 
acf(Serie_PM_estacionaria)
pacf(Serie_PM_estacionaria)
acf2(Serie_PM_estacionaria)
frequency(Serie_PM)



### 6. Ajuste del modelo ARIMA estacional
#Arima1: ARIMA (1, 1, 1) (1, 0, 0) (12) 
#Arima2: ARIMA (2, 0, 0) (0, 1, 2) (12)
#Arima3: ARIMA (2, 0, 1) (0, 1, 2) (12)
#Arima4: ARIMA (2, 0, 2) (1, 1, 2) (12)
#Arima5: ARIMA (2, 0, 2) (0, 1, 2) (12)
#Arima6: ARIMA (0, 0, 1) (0, 1, 1) (12)
#Arima7: ARIMA (0, 1, 1) (0, 1, 1) (12)

Arima1<- arima(Serie_PM, order=c(1,1,1),seasonal = list(order = c(1,0,0), period = 12),method="ML")
coeftest(Arima1)
Arima2<- arima(Serie_PM, order=c(2,0,0),seasonal = list(order = c(0,1,2), period = 12),method="ML")
coeftest(Arima2) #Coeficientes significantes 
Arima3<- arima(Serie_PM, order=c(2,0,1),seasonal = list(order = c(0,1,2), period = 12),method="ML")
coeftest(Arima3)
Arima4<- arima(Serie_PM, order=c(2,0,2),seasonal = list(order = c(1,1,2), period = 12),method="ML")
coeftest(Arima4)
Arima5<- arima(Serie_PM, order=c(2,0,2),seasonal = list(order = c(0,1,2), period = 12),method="ML")
coeftest(Arima5)
Arima6<- arima(Serie_PM, order=c(0,0,1),seasonal = list(order = c(0,1,1), period = 12),method="ML")
coeftest(Arima6)
Arima7<- arima(Serie_PM, order=c(0,1,1),seasonal = list(order = c(0,1,1), period = 12),method="ML")
coeftest(Arima7)
Arima8<- arima(Serie_PM, order=c(1,1,0),seasonal = list(order = c(1,1,0), period = 12),method="ML")
coeftest(Arima8)


AIC(Arima1,Arima2,Arima3,Arima4,Arima5,Arima6,Arima7,Arima8) #Best: Arima4
BIC(Arima1,Arima2,Arima3,Arima4,Arima5,Arima6,Arima7,Arima8) #Best: Arima4

auto.arima(Serie_PM, stepwise = FALSE, approximation = FALSE)
#Nota: Sugiere modelo Arima2



### 7. Diagnóstico de los modelos
ggtsdiag(Arima8)
ggtsdiag(Arima2)

#Pruebas adicionales a los residuales
Box.test(Arima2$residuals) # Test de Box-Pierce #Independencia
Box.test(Arima2$residuals, type="Ljung-Box") #Independencia
jarque.bera.test(Arima2$residuals) #Normalidad
par(mfrow=c(1,2),cex=0.7)
hist(Arima2$residuals, main="a)",xlab = "Residuales",ylab = "Frecuencia")
qqPlot(Arima2$residuals,distribution = "norm",main = "b)",ylab = "Residuales")
dev.off()

Box.test(Arima8$residuals) # Test de Box-Pierce #Independencia
Box.test(Arima8$residuals, type="Ljung-Box") #Independencia
jarque.bera.test(Arima8$residuals) #Normalidad
par(mfrow=c(1,2),cex=0.7)
hist(Arima8$residuals, main="a)",xlab = "Residuales",ylab = "Frecuencia")
qqPlot(Arima8$residuals,distribution = "norm",main = "b)",ylab = "Residuales")
dev.off()



### 8. Prediccion del modelo
##Predicción (con el mejor modelo: ARIMA2)
forecast::forecast()
#Arima2
pred1_arima2<-forecast::forecast(Arima2,h=12, level=c(99.5))
Fig4<-autoplot(pred1_arima2, main = "2018")

pred2_arima2<-forecast::forecast(Arima2,h=24, level=c(99.5))
Fig5<-autoplot(pred2_arima2, main = "2019")

pred3_arima2<-forecast::forecast(Arima2,h=36, level=c(99.5))
Fig6<-autoplot(pred3_arima2, main = "2020")

ggarrange(Fig4,Fig5,Fig6,ncol = 1, nrow = 3)

#Arima 8
pred4_arima8<-forecast::forecast(Arima8,h=12, level=c(99.5))
Fig4<-autoplot(pred4_arima8, main = "2018")

pred5_arima8<-forecast::forecast(Arima8,h=24, level=c(99.5))
Fig5<-autoplot(pred5_arima8, main = "2019")

pred6_arima8<-forecast::forecast(Arima8,h=36, level=c(99.5))
Fig6<-autoplot(pred6_arima8, main = "2020")

ggarrange(Fig4,Fig5,Fig6,ncol = 1, nrow = 3)
