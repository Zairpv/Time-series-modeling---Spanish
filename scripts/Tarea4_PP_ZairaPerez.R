#TAREA 4
#Presenta: Zaira Rosario Pérez Vázquez

#SEGUNDA SERIE DE TIEMPO: PRECIPITACIÓN
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

### 1. Segundo conjunto de datos - Precipitación total anual
datos_PP<-read.csv("Precipitacion_SerieAnual.csv",header = FALSE)
View(datos_PP)
PP_serie<-ts(datos_PP,start = 1945,frequency = 1,names = "Precipitacion")
PP_serie

##1.1 Reconocimiento de la serie

#1.1.1. Estadistica descriptiva de la serie
autoplot(PP_serie,xlab = "Tiempo",ylab = "Precipitación (mm)")
summary(PP_serie)
var(PP_serie)
sd(PP_serie)
par(mfrow=c(1,3),cex=0.9)
hist(PP_serie,main = "a)",xlab = "Precipitación total anual (mm)",ylab = "Frecuencia")
plot(density(PP_serie), main = "b)",ylab = "Densidad")
qqPlot(PP_serie,dist="norm",ylab = "Precipitación",main = "c)")
dev.off() #800x320px

#1.1.2 Prueba de normalidad
jarque.bera.test(PP_serie) #La serie no sigue una distribución normal
shapiro.test(PP_serie)
#Nota: No se puede aplicar la metodología de Box-Jenkins
#Se requiere de la transformación de los datos. 

#1.1.3 Transformación de la serie

#Transformación logaritmica
PP_serie_log<-log(PP_serie)
par(mfrow=c(1,3),cex=0.9)
hist(PP_serie_log,main = "a)",xlab = "Log Precipitación",ylab = "Frecuencia")
plot(density(PP_serie_log), main = "b)",ylab = "Densidad")
qqPlot(PP_serie_log,main = "c)",ylab = "Log Precipitación")
dev.off()
jarque.bera.test(PP_serie_log)
lillie.test(PP_serie_log) #solo esta prueba sugiere normalidad
ad.test(PP_serie_log)
serie_orig_Fig<-autoplot(PP_serie,xlab = "Tiempo",ylab = "Precipitación")
serie_log_Fig<-autoplot(PP_serie_log,xlab = "Tiempo",ylab = "Log Precipitación")
ggarrange(serie_orig_Fig,serie_log_Fig,ncol = 1, nrow = 2,labels = c("A","B"))

#Transformación BoxCox
lambda<-BoxCox.lambda(PP_serie)
lambda
PP_serie_boxcox<-BoxCox(PP_serie, lambda = lambda)
PP_serie_boxcox
par(mfrow=c(1,3),cex=0.9)
hist(PP_serie_boxcox,main = "a)",xlab = "BoxCox - Precipitación",ylab = "Frecuencia")
plot(density(PP_serie_boxcox), main = "b)",ylab = "Densidad")
qqPlot(PP_serie_boxcox,main = "c)",ylab = "BoxCox - Precipitación")
dev.off()
jarque.bera.test(PP_serie_boxcox) #Normalidad
lillie.test(PP_serie_boxcox)
ad.test(PP_serie_boxcox) #Sin normalidad
serie_bc_Fig<-autoplot(PP_serie_log,xlab = "Tiempo",ylab = "Log Precipitación")
ggarrange(serie_orig_Fig,serie_bc_Fig,ncol = 1, nrow = 2,labels = c("A","B"))
#Nota: la transformación Boxcox es la más adecuada para obtener normalidad
#en la serie de precipitación


##1.2. Identificación de la señal de la serie
#Serie original
tt_pp <- 1:length(PP_serie)
fit_pp <- ts(loess(PP_serie ~ tt_pp, span = .2)$fitted, start = 1945, frequency = 1)
plot.ts(fit_pp, type='l')
grid()
lines(PP_serie, col = "red") 

#Serie transformada
tt_pp_bc <- 1:length(PP_serie_boxcox)
fit_pp_bc <- ts(loess(PP_serie_boxcox ~ tt_pp_bc, span = .2)$fitted, start = 1945, frequency = 1)
plot.ts(fit_pp_bc, type='l')
grid()
lines(PP_serie_boxcox, col = "red") 


##1.3. Estimación del nivel de la serie (y=bo+e)
#serie original
fit_level_pp<-lm(PP_serie~1)
summary(fit_level_pp)
plot.ts(PP_serie, main="")
lines(ts(fitted(fit_level_pp), start=1945, frequency = 1), col = "red")
#Nota: el nivel de la serie se encuentra en un 
#valor de precipitación igual a 1486.41 mm  

#serie original transformada
fit_level_pp_bc<-lm(PP_serie_boxcox~1)
summary(fit_level_pp_bc)
plot.ts(PP_serie_boxcox, main="")
lines(ts(fitted(fit_level_pp_bc), start=1945, frequency = 1), col = "red")


##1.4 Identificación de cambios estructurales de nivel
#1.4.1. Serie original
#Breakpoints
pp_brk <- breakpoints(PP_serie ~ 1, h = 0.1)
summary(pp_brk)
plot(pp_brk)
#Nota: De acuerdo con  BIC, aproximadamente existen 2 cambios de nivel
breakdates(pp_brk, breaks = 2)
plot(PP_serie,col="dark gray")
lines(fitted(pp_brk, breaks = 2), col = 4)
lines(confint(pp_brk, breaks = 2))
coef(pp_brk, breaks = 2)

#1.4.2. Serie transformada
#Breakpoints
pp_brk_bc <- breakpoints(PP_serie_boxcox ~ 1, h = 0.1)
summary(pp_brk_bc)
plot(pp_brk_bc)
#Nota: De acuerdo con  BIC, aproximadamente existen 4 cambios de nivel
breakdates(pp_brk_bc, breaks = 4)
plot(PP_serie_boxcox,col="dark gray")
lines(fitted(pp_brk_bc, breaks = 4), col = 4)
lines(confint(pp_brk_bc, breaks = 4))
coef(pp_brk_bc, breaks = 4)


#1.5. Tendencias

#1.5.1. Serie original

#Ajuste del modelo de regresión lineal
tt_pp
trend_fit_pp <- lm(PP_serie ~ tt_pp)
summary(trend_fit_pp)
#Nota: hay presencia de tendencia significativa
plot(PP_serie)
lines(ts(fitted(trend_fit_pp), start=1945, frequency = 1), col = "red")
#La tendencia es descendente - la serie no es estacionaria en media.

#Obtención de breakpoints
pp_brk_trend <- breakpoints(PP_serie ~ tt_pp, h = 0.1)
summary(pp_brk_trend)
plot(pp_brk_trend)
#Nota: De acuerdo con BIC, existen 2 cambios estructurales de tendencia
breakdates(pp_brk_trend, breaks = 2)
coef(pp_brk_trend, breaks = 2)
plot(PP_serie,col="dark gray")
lines(fitted(pp_brk_trend, breaks = 2), col = 4)
lines(confint(pp_brk_trend, breaks = 2))

#1.5.1. Serie transformada

#Ajuste del modelo de regresión lineal
tt_pp_bc
trend_fit_pp_bc <- lm(PP_serie_boxcox ~ tt_pp_bc)
summary(trend_fit_pp_bc)
#Nota: hay presencia de tendencia significativa
plot(PP_serie_boxcox)
lines(ts(fitted(trend_fit_pp_bc), start=1945, frequency = 1), col = "red")
#La tendencia es descendente - la serie no es estacionaria en media.

#Obtención de breakpoints
pp_brk_trend_bc <- breakpoints(PP_serie_boxcox ~ tt_pp_bc, h = 0.1)
summary(pp_brk_trend_bc)
plot(pp_brk_trend_bc)
#Nota: De acuerdo con BIC, existen 2 cambios estructurales de tendencia
breakdates(pp_brk_trend_bc, breaks = 2)
coef(pp_brk_trend_bc, breaks = 2)
plot(PP_serie_boxcox,col="dark gray")
lines(fitted(pp_brk_trend_bc, breaks = 2), col = 4)
lines(confint(pp_brk_trend_bc, breaks = 2))

##1.6. Identificación de un modelo ARIMA

#1.6.1. Estacionariedad
#Serie original
adf.test(PP_serie) #La serie es estacionaria
#Otra forma: 
adfTest(PP_serie,lags=0,type=c("ct"))
#Nota: se selecciona el tipo ct debido a la tendencia existente
#Nota: La serie es  estacionaria.

#Serie transformada
adfTest(PP_serie_boxcox,lags=0,type=c("ct")) #La serie es estacionaria
#Nota: se selecciona el tipo ct debido a la tendencia existente


#1.6.2. Función acf y pacf
#Original
acf(PP_serie)
pacf(PP_serie)
acf2(PP_serie)
#Transformada
acf(PP_serie_boxcox)
pacf(PP_serie_boxcox)
acf2(PP_serie_boxcox)

#1.6.3. ¿Diferenciación de la serie?
ndiffs(PP_serie_boxcox) #d = 1
nsdiffs(PP_serie_boxcox) #No existen datos estacionales
#¿Y sin transformación?
ndiffs(PP_serie)
PP_serie_diff1 <- diff(PP_serie, lag = 1)
ad.test(PP_serie_diff1) 
adf.test(PP_serie_diff1)
#Nota. La serie diferenciada si es estacionaria, pero no
#es normal. 

##1.7. Ajuste de un modelo ARIMA
#Propuestas: 
#Arima1: ARIMA (2,0,3)
#Arima2: ARIMA (2,0,0)
#Arima3: ARIMA (0,0,3)

Arima1<-arima(PP_serie_boxcox, c(2, 0, 3)) #ARMA
Arima2<-arima(PP_serie_boxcox, c(2, 0, 0)) #AR
Arima3<-arima(PP_serie_boxcox, c(0, 0, 3)) #MA


AIC(Arima1,Arima2,Arima3)
BIC(Arima1,Arima2,Arima3)

?auto.arima
auto.arima(PP_serie_boxcox, stepwise = FALSE, approximation = FALSE)


##1.8. Diagnostico de un modelo ARIMA
ggtsdiag(Arima3) #Modelo seleccionado
ggtsdiag(Arima1) 
ggtsdiag(Arima2)

#Pruebas adicionales a los residuales
Box.test(Arima3$residuals) # Test de Box-Pierce #Independencia
Box.test(Arima3$residuals, type="Ljung-Box") #Independencia
jarque.bera.test(Arima3$residuals) #Normalidad
shapiro.test(Arima3$residuals) #Normalidad
par(mfrow=c(1,2))
hist(Arima3$residuals, main="a)",xlab = "Residuales",ylab = "Frecuencia")
qqPlot(Arima3$residuals,distribution = "norm",main = "b)",ylab = "Residuales")
dev.off()


##1.9. Predicción (con el mejor modelo: ARIMA3)
#2025
pp_pred_arima3<-forecast(Arima3,level = c(95), h = 6)
autoplot(pp_pred_arima3)

#2030
pp_pred2_arima3<-forecast(Arima3,level = c(95), h = 11)
autoplot(pp_pred2_arima3) 

#Con los otros modelos tampoco se alcanza a predecir hasta 2030
#Arima1
#2025
pp_pred_arima1<-forecast(Arima1,level = c(95), h = 6)
autoplot(pp_pred_arima1)
#2030
pp_pred2_arima1<-forecast(Arima1,level = c(95), h = 11)
autoplot(pp_pred2_arima1) 
#Arima2
#2025
pp_pred_arima2<-forecast(Arima2,level = c(95), h = 6)
autoplot(pp_pred_arima2)
#2030
pp_pred2_arima2<-forecast(Arima2,level = c(95), h = 11)
autoplot(pp_pred2_arima2) 


