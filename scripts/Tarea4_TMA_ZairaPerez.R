#TAREA 4
#Presenta: Zaira Rosario Pérez Vázquez

#PRIMER SERIE DE TIEMPO: TEMPERATURA

rm(list=ls())
options(max.print = 10000)
setwd("C:/Users/zaira/Documents/Doctorado PCF/3. Cursos/2. Primavera 2021/3. Series de tiempo/3. Tareas/Tarea 4")

#Librerias
library(ggfortify)
library(tseries)
library(car)
library(nortest)
library(forecast)
library(ggpubr)
library(strucchange)
library(fUnitRoots)
library(astsa)
library(lmtest)

### 1. Primer conjunto de datos - Temperatura media anual
datos_TMA<-read.csv("Temperatura_SerieAnual.csv",header = FALSE)
View(datos_TMA)
TMA_serie<-ts(datos_TMA,start = 1945,frequency = 1,names = "Temperatura")
TMA_serie 

##1.1 Reconocimiento de la serie

#1.1.1. Estadistica descriptiva de la serie
autoplot(TMA_serie,xlab = "Tiempo",ylab = "Temperatura (°C)")
summary(TMA_serie)
var(TMA_serie)
sd(TMA_serie)
par(mfrow=c(1,3),cex=0.9)
hist(TMA_serie,main = "a)",xlab = "Temperatura media anual (°C)",ylab = "Frecuencia")
plot(density(TMA_serie), main = "b)",ylab = "Densidad")
qqPlot(TMA_serie,dist="norm",ylab = "Temperatura",main = "c)")
dev.off() #800x320px

#1.1.2 Prueba de normalidad
jarque.bera.test(TMA_serie) #La serie sigue una distribución normal
shapiro.test(TMA_serie)
#Nota: Se puede aplicar la metodología de Box-Jenkins


##1.2. Identificación de la señal de la serie
tt_tma <- 1:length(TMA_serie)
fit_tma <- ts(loess(TMA_serie ~ tt_tma, span = .2)$fitted, start = 1945, frequency = 1)
plot.ts(fit_tma, type='l')
grid()
lines(TMA_serie, col = "red") 


##1.3. Estimación del nivel de la serie (y=bo+e)
fit_level_tma<-lm(TMA_serie~1)
summary(fit_level_tma)
plot.ts(TMA_serie, main="Nivel de la serie")
lines(ts(fitted(fit_level_tma), start=1945, frequency = 1), col = "red")
#Nota: el nivel de la serie se encuentra en un 
#valor de temperatura igual a 13.508°C (Valor muy cercano a la media y mediana de la serie)


##1.4 Identificación de cambios estructurales de nivel
#1.4.1. Breakpoints
temp_brk <- breakpoints(TMA_serie ~ 1, h = 0.1)
summary(temp_brk)
plot(temp_brk)
#Nota: De acuerdo con AIC y BIC, aproximadamente existen 6 cambios de nivel

breakdates(temp_brk, breaks = 6)
plot(TMA_serie,col="dark gray")
lines(fitted(temp_brk, breaks = 6), col = 4)
lines(confint(temp_brk, breaks = 6))
coef(temp_brk, breaks = 6)

#1.5. Tendencias
#1.5.1. Ajuste del modelo de regresión lineal
tt_tma
trend_fit_tma <- lm(TMA_serie ~ tt_tma)
summary(trend_fit_tma)
#Nota: hay presencia de tendencia significativa
plot(TMA_serie)
lines(ts(fitted(trend_fit_tma), start=1945, frequency = 1), col = "red")
#La tendencia es ascendente - la serie no es estacionaria en media.

#1.5.2. Obtención de breakpoints
temp_brk_trend <- breakpoints(TMA_serie ~ tt_tma, h = 0.1)
summary(temp_brk_trend)
plot(temp_brk_trend)
#Nota: De acuerdo con AIC y BIC, existen seis cambios estructurales de tendencia
breakdates(temp_brk_trend, breaks = 6)
coef(temp_brk_trend, breaks = 6)
plot(TMA_serie,col="dark gray")
lines(fitted(temp_brk_trend, breaks = 6), col = 4)
lines(confint(temp_brk_trend, breaks = 6))


##1.6. Identificación de un modelo ARIMA

#1.6.1. Estacionariedad
adf.test(TMA_serie) #La serie no es estacionaria
#Otra forma: 
adfTest(TMA_serie,lags=0,type=c("ct"))
#Nota: se selecciona el tipo ct debido a la tendencia existente
#Nota: La serie es no estacionaria.

#1.6.2. Función acf y pacf
acf(TMA_serie)
pacf(TMA_serie)
acf2(TMA_serie)

#1.6.3. Diferenciación de la serie
ndiffs(TMA_serie) #d = 1
nsdiffs(TMA_serie) #No existen datos estacionales
TMA_serie_diff1 <- diff(TMA_serie, lag = 1)

#Comparación entre serie original y diferenciada en lag=1
both<-cbind(TMA_serie, TMA_serie_diff1)
head(both)
plot(both) 
#Otra forma de ver la gráfica
autoplot(both,xlab = "Tiempo",ylab = "Temperatura (°C)")

#Verificar estacionariedad
plot(TMA_serie_diff1)
acf(TMA_serie_diff1)
pacf(TMA_serie_diff1)
acf2(TMA_serie_diff1)

adf.test(TMA_serie_diff1) #La serie es estacionaria
#Gráficas
par(mfrow=c(2,2))
acf(TMA_serie,main="ACF - Original")
pacf(TMA_serie,main="PACF - Original")
acf(TMA_serie,main="ACF - Diff (1)")
pacf(TMA_serie_diff1,main="PACF - Diff(1)")
dev.off()


##1.7. Ajuste de un modelo ARIMA
#Propuestas: 
#Arima1: ARIMA (0,1,0)
#Arima2: ARIMA (1,1,0)
#Arima3: ARIMA (0,1,14)
#Arima4: ARIMA (1,1,14)

Arima1<-arima(TMA_serie, c(0, 1, 0))
Arima2<-arima(TMA_serie, c(1, 1, 0))
Arima3<-arima(TMA_serie, c(0, 1, 14)) #No cumplen el principio de parsimonia (ilustrativo)
Arima4<-arima(TMA_serie, c(1, 1, 14))

AIC(Arima1,Arima2,Arima3,Arima4)
BIC(Arima1,Arima2,Arima3,Arima4)

?auto.arima
auto.arima(TMA_serie, stepwise = FALSE, approximation = FALSE)


##1.8. Diagnostico de un modelo ARIMA
ggtsdiag(Arima1) #Modelo seleccionado
ggtsdiag(Arima2)
ggtsdiag(Arima3)
ggtsdiag(Arima4)

#Pruebas adicionales a los residuales
Box.test(Arima1$residuals) # Test de Box-Pierce #Independencia
Box.test(Arima1$residuals, type="Ljung-Box") #Independencia
jarque.bera.test(Arima1$residuals) #Normalidad
shapiro.test(Arima1$residuals) #Normalidad
par(mfrow=c(1,2))
hist(Arima1$residuals, main="a)",xlab = "Residuales",ylab = "Frecuencia")
qqPlot(Arima1$residuals,distribution = "norm",main = "b)",ylab = "Residuales")
dev.off()

#Verificar normalidad en Arima2,3 y 4
#Arima2
jarque.bera.test(Arima2$residuals) #Sin Normalidad
shapiro.test(Arima2$residuals) #Sin Normalidad
#Arima3
jarque.bera.test(Arima3$residuals) #Sin Normalidad
shapiro.test(Arima3$residuals) #Normalidad
#Arima4
jarque.bera.test(Arima4$residuals) #Sin Normalidad
shapiro.test(Arima4$residuals) #Normalidad

#Nota: al parecer el modelo Arima3 es el más cercano a una distribución normal
Box.test(Arima3$residuals) # Test de Box-Pierce #Independencia
Box.test(Arima3$residuals, type="Ljung-Box") #Independencia
par(mfrow=c(1,2))
hist(Arima3$residuals, main="a)",xlab = "Residuales",ylab = "Frecuencia")
qqPlot(Arima3$residuals,distribution = "norm",main = "b)",ylab = "Residuales")
dev.off()
#Nota: Este modelo cumple con todos los supuestos de ruido blanco.


##1.9. Predicción (con el mejor modelo: ARIMA6)
#2025
tma_pred_arima1<-forecast(Arima1,level = c(95), h = 6)
autoplot(tma_pred_arima1)
tma_pred_arima3<-forecast(Arima3,level = c(95), h = 6)
autoplot(tma_pred_arima3)

#2030
tma_pred2_arima1<-forecast(Arima1,level = c(95), h = 11)
autoplot(tma_pred2_arima1)
tma_pred2_arima3<-forecast(Arima3,level = c(95), h = 11)
autoplot(tma_pred2_arima3)

