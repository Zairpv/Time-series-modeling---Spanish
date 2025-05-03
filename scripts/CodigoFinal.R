rm(list=ls())
options(max.print = 10000)
setwd("C:/Users/zaira/Documents/Doctorado PCF/3. Cursos/2. Primavera 2021/3. Series de tiempo/3. Tareas/Tarea 8")
getwd()

### 1. Leer datos. 
library(readr)
datos_petroleo <- read_csv("datos_petroleo.csv")
View(datos_petroleo)

#Librerias
library(tseries)
library(ggfortify)
library(fBasics)
library(car)
library(tidyverse)
library(ggplot2)
library(Rmisc)
library(forecast)
library(nortest)
library(strucchange)
library(tsoutliers)
library(fUnitRoots)
library(astsa)
library(dplyr)
library(lmtest)
library(aTSA)
library(TSA)


#Definir variables
gasolinaY_ts <- ts(datos_petroleo$gasolina,start = c(2013,01),frequency = 12)
brentX_ts <- ts(datos_petroleo$brent,start = c(2013,01),frequency = 12)


### 1. Análisis exploratorio de Yt y Xt
basicStats(gasolinaY_ts)
basicStats(brentX_ts)

## a)  Análisis exploratorio de Y
autoplot(gasolinaY_ts,xlab = "Tiempo",ylab = "Gasolina (US$)",main = "a)")
par(mfrow=c(1,3),mar=c(5,4,3,1),cex=0.6)
hist(gasolinaY_ts,main = "b)",xlab = "Gasolina (US$)",ylab = "Frecuencia")
plot(density(gasolinaY_ts), main = "c)",ylab = "Densidad")
qqPlot(gasolinaY_ts,dist="norm",ylab = "Gasolina (US$)",main = "d)")
dev.off() 
jarque.bera.test(gasolinaY_ts)

## b) Transformación de Yt 
#Transformación logaritmica
jarque.bera.test(log(gasolinaY_ts))
#Transformación BoxCox
(lambdaY<-BoxCox.lambda(gasolinaY_ts))
gasolinaY_BC<-BoxCox(gasolinaY_ts,lambda = lambdaY)
jarque.bera.test(gasolinaY_BC)
#Comparación entre transformaciones
par(mfrow=c(2,3),mar=c(5,4,3,1),cex=0.6)
hist(gasolinaY_ts,main = "Original",xlab = "Gasolina (US$)",ylab = "Frecuencia")
hist(log(gasolinaY_ts),main = "Logarítmica",xlab = "Gasolina (US$)",ylab = "Frecuencia")
hist(gasolinaY_BC,main = expression(paste("BoxCox (", lambda, " = -0.4164)")),xlab = "Gasolina (BoxCox)",ylab = "Frecuencia")
qqPlot(gasolinaY_ts,dist="norm",ylab = "Gasolina (US$)",main = "")
qqPlot(log(gasolinaY_ts),dist="norm",ylab = "Gasolina (US$)",main = "")
qqPlot(gasolinaY_BC,dist="norm",ylab = "Gasolina (US$)",main = "")
dev.off()

## c) Análisis exploratorio de X1
autoplot(brentX_ts,xlab = "Tiempo",ylab = "Brent (US$)",main = "a)")
par(mfrow=c(1,3),mar=c(5,4,3,1),cex=0.6)
hist(brentX_ts,main = "b)",xlab = "Brent (US$)",ylab = "Frecuencia")
plot(density(brentX_ts), main = "c)",ylab = "Densidad")
qqPlot(brentX_ts,dist="norm",ylab = "Gasolina (US$)",main = "d)")
dev.off() 
jarque.bera.test(brentX_ts)

## d) Transformación BoxCox de Xt
(lambdaX1<-BoxCox.lambda(brentX_ts))
brentX1_BC<-BoxCox(brentX_ts,lambda = lambdaX1)
hist(brentX1_BC,main = "",xlab = "Brent (BC)",ylab = "Frecuencia")
qqPlot(brentX1_BC,dist="norm",ylab = "Brent (BC)",main = "")
jarque.bera.test(brentX1_BC)



### 2. Identificación de la señal y nivel de la serie Y
## a) Señal de la serie
par(mfrow=c(2,1),cex=0.9,mar=c(5,4,3,1),cex=0.6)
tt_gas<-1:length(gasolinaY_BC)
fit_gas<-ts(loess(gasolinaY_BC~tt_gas,span = 0.2)$fitted,start = 2013,frequency=12)
plot.ts(fit_gas,type="l",main="a)",ylab=expression(paste("Gasolina (", lambda, " = -0.4164)")))
grid()
lines(gasolinaY_BC,col="Red")

## b) Nivel de la serie
fitlevel_gas<-lm(gasolinaY_BC~1)
summary(fitlevel_gas)
plot.ts(gasolinaY_BC,main="b)",ylab=expression(paste("Gasolina (", lambda, " = -0.4164)")))
lines(ts(fitted(fitlevel_gas),start = 2013,frequency = 12),col="red")
dev.off()


### 3. Cambios estructurales Y
## a) Cambios de nivel
gasBC_brk<-breakpoints(gasolinaY_BC~1,h=0.1)
summary(gasBC_brk)
par(mfrow=c(1,1),cex=0.9,mar=c(5,4,3,1))
plot(gasBC_brk) #BIC=6 brk
breakdates(gasBC_brk,breaks = 6)
plot(gasolinaY_BC,col="dark gray",ylab=expression(paste("Gasolina (", lambda, " = -0.4164)")),xlab="Tiempo")
lines(fitted(gasBC_brk,breaks = 6),col=4)
lines(confint(gasBC_brk,breaks = 6))
coef(gasBC_brk,breaks = 6)

## b) Cambios de tendencia
trend_fit_gas<-lm(gasolinaY_BC~tt_gas)
summary(trend_fit_gas)
plot(gasolinaY_BC,ylab=expression(paste("Gasolina (", lambda, " = -0.4164)")),xlab="Tiempo")
lines(ts(fitted(trend_fit_gas),start = 2013,frequency = 12),col="red")
gasBC_brk_trend<-breakpoints(gasolinaY_BC~tt_gas,h=0.1)
summary(gasBC_brk_trend)
plot(gasBC_brk_trend) #6 Brk
breakdates(gasBC_brk_trend,breaks = 6)
coef(gasBC_brk_trend,breaks = 6)
plot(gasolinaY_BC,col="dark gray",ylab=expression(paste("Gasolina (", lambda, " = -0.4164)")),xlab="Tiempo")
lines(fitted(gasBC_brk_trend,breaks = 6),col=4)
lines(confint(gasBC_brk_trend,breaks = 6))
coef(gasBC_brk_trend,breaks = 6)


### 4. Raiz unitaria de las series Yt, Xt
##Serie Y
par(mfrow=c(1,1),cex=0.9,mar=c(5,4,3,1))
acf(gasolinaY_BC,main="a)")
F1<-autoplot(gasolinaY_BC,ylab = "Gasolina",main="Y: BoxCox")
F2<-autoplot(diff(gasolinaY_BC),ylab = "Gasolina",main="Y: BoxCox + d=1")
acf(diff(gasolinaY_BC),main="b)")
mod1<-ar(diff(gasolinaY_BC),method = "mle")
mod1$order #Fue = 3
#Hay constante? es significativa?
t_Y=seq(2:length(gasolinaY_BC))
mod2<-lm(diff(gasolinaY_BC)~t_Y)
summary(mod2) #No es significativa
#Ajustar a 3 lags
adfTest(gasolinaY_BC,lags = 3,type = "nc")
adfTest(diff(gasolinaY_BC),lags = 3,type = "nc")
#La serie transformada tiene raiz unitaria y por lo tanto necesita diferenciarse
dev.off()

##Serie X1
par(mfrow=c(1,1),cex=0.9,mar=c(5,4,3,1))
acf(brentX1_BC,main="a)")
F3<-autoplot(brentX1_BC,ylab = "Brent",main="X1: BoxCox")
F4<-autoplot(diff(brentX1_BC),ylab = "Brent",main="X1: BoxCox + d=1")
acf(diff(brentX1_BC),main="b)")
mod3<-ar(diff(brentX1_BC),method = "mle")
mod3$order #Fue = 1
#Hay constante? es significativa?
t_X1=seq(2:length(brentX1_BC))
mod4<-lm(diff(brentX1_BC)~t_X1)
summary(mod4) #No es significativa
#Ajustar a 1 lag
adfTest(brentX1_BC,lags = 1,type = "nc")
adfTest(diff(brentX1_BC),lags = 1,type = "nc")
dev.off()
#La serie transformada tiene raiz unitaria y por lo tanto necesita diferenciarse

multiplot(F1,F2,F3,F4)


### 5. Correlación entre pares de variables 
#La correlacion debe ser entre variables estacionarias
par(mfrow=c(1,1),cex=0.9,mar=c(5,4,3,1))
ccf(diff(brentX1_BC),diff(gasolinaY_BC),ylab="Cross correlation",main="Precio Brent & Precio Gasolina",type = "correlation")
dev.off()


### 6. Test de causalidad Granger -- Ho: X no causa a Y
grangertest(diff(brentX1_BC),diff(gasolinaY_BC),order = 3)
#¿Existe una relación inversa?
grangertest(diff(gasolinaY_BC),diff(brentX1_BC))


### 7. Preblanqueamiento Xt
acf(diff(brentX1_BC))
pacf(diff(brentX1_BC))
acf2(diff(brentX1_BC),main = "Serie: Brent (d=1)")

## a) Modelos arima propuestos
Arima1_x<-Arima(brentX1_BC,order = c(0,1,0),seasonal = list(order = c(1,0,0),period = 12),method = "ML",lambda = lambdaX1)
coeftest(Arima1_x)
Arima2_x<-Arima(brentX1_BC,order = c(0,1,0),seasonal = list(order = c(0,0,1),period = 12),method = "ML",lambda = lambdaX1)
coeftest(Arima2_x)
#Modelo sin estacionalidad
Arima3_x<-Arima(gasolinaY_BC,order = c(0,1,0),method = "ML",lambda = lambdaX1)
Arima4_x<-Arima(gasolinaY_BC,order = c(1,1,0),method = "ML",lambda = lambdaX1)
coeftest(Arima4_x)
Arima5_x<-Arima(gasolinaY_BC,order = c(0,1,1),method = "ML",lambda = lambdaX1)
coeftest(Arima5_x)

BIC(Arima1_x,Arima2_x,Arima3_x,Arima4_x,Arima5_x)
AIC(Arima1_x,Arima2_x,Arima3_x,Arima4_x,Arima5_x)

## b) Diagnostico del modelo
ggtsdiag(Arima5_x)
Box.test(Arima5_x$residuals) # Test de Box-Pierce #Independencia
Box.test(Arima5_x$residuals, type="Ljung-Box") #Independencia
jarque.bera.test(Arima5_x$residuals) #Normalidad
par(mfrow=c(1,2),cex=0.9,mar=c(5,4,3,1))
hist(Arima5_x$residuals, main="a)",xlab = "Residuales",ylab = "Frecuencia")
qqPlot(Arima5_x$residuals,distribution = "norm",main = "b)",ylab = "Residuales")

## c) Una forma de preblanqueamiento general
prewhiten(x=diff(brentX1_BC),y=diff(gasolinaY_BC))
#Nota: la gráfica es la misma que la obtenida en el paso 5. 


## d) Identificación de outliers
tso(brentX1_BC)
plot(tso(brentX1_BC))
#Se identifico un outlier LS (cambio de nivel), 
#¿y en residuales?
res_Arima5_x<-residuals(Arima5_x)
polin_X<-coefs2poly(Arima5_x)
(outliersX1<-locate.outliers(res_Arima5_x,polin_X)) 
#Nota. esta función solo se utiliza con residuales
plot(tso(res_Arima5_x))


## e) Incorporar intervención al modelo de X1 
#Crear vector binario
length(brentX1_BC)
LS_X1<-outliers("LS",71)
(LSef_X1<-outliers.effects(LS_X1,length(brentX1_BC)))

library(TSA)
#Omega
Arima_xLS<-arimax(brentX1_BC,order = c(0,1,1),method = "ML",
                    xreg = LSef_X1)
coeftest(Arima_xLS)

#Omega y delta
Arima2_xLS<-arimax(brentX1_BC,order = c(0,1,1),method = "ML",
                     xtrans=LSef_X1,transfer = list(c(1,0)))
coeftest(Arima2_xLS)
#vALOR Delta no significativo. 

## f) Diagnostico del modelo
ggtsdiag(Arima2_xLS)
Box.test(Arima2_xLS$residuals) # Test de Box-Pierce #Independencia
Box.test(Arima2_xLS$residuals, type="Ljung-Box") #Independencia
hist(Arima2_xLS$residuals, main="a)",xlab = "Residuales",ylab = "Frecuencia")
qqPlot(Arima2_xLS$residuals,distribution = "norm",main = "b)",ylab = "Residuales")
dev.off()
plot(tso(Arima_x3_LS2$residuals))


### 8. Preblanqueamiento Yt
## a) Modelo arima univariado
acf(diff(gasolinaY_BC))
pacf(diff(gasolinaY_BC))
acf2(diff(gasolinaY_BC))

#Propuestas
Arima1_y<-Arima(gasolinaY_BC,order = c(0,1,1),method = "ML",lambda = lambdaY)
coeftest(Arima1_y) #Mismo modelo que el seleccionado por X
BIC(Arima1_y)
AIC(Arima1_y)

## b) Diagnostico
ggtsdiag(Arima1_y)
Box.test(Arima1_y$residuals) # Test de Box-Pierce #Independencia
Box.test(Arima1_y$residuals, type="Ljung-Box") #Independencia
jarque.bera.test(Arima1_y$residuals) #Normalidad
par(mfrow=c(1,2),cex=0.9,mar=c(5,4,3,1))
hist(Arima1_y$residuals, main="a)",xlab = "Residuales",ylab = "Frecuencia")
qqPlot(Arima1_y$residuals,distribution = "norm",main = "b)",ylab = "Residuales")
dev.off()

## c) Identificación de outliers en Yt
tso(gasolinaY_BC)
plot(tso(gasolinaY_BC))
#No se identificaron outliers, ¿y en residuales?
resid_Arima1_y<-residuals(Arima1_y)
polin_y<-coefs2poly(Arima1_y)
(outliersY1<-locate.outliers(resid_Arima1_y,polin_y))
plot(tso(resid_Arima1_y))


### 9. Correlación cruzada entre series preblanqueadas
## Se utilizaron residuales de modelos X y Y 
#Correlación cruzada entre series preblanqueadas
par(mfrow=c(1,1),cex=0.9,mar=c(5,4,3,1))
ccf(Arima2_xLS$residuals,Arima1_y$residuals,ylab="Cross correlation",main="Precio Brent & Precio Gasolina",type = "correlation",na.action = na.omit)
dev.off()

### 10. Función de transferencia

### a). Modelo bivariado sin intervención
Arima1_xy_si<-arimax(gasolinaY_BC,order = c(0,1,1),xreg = brentX1_BC)
BIC(Arima1_xy_si)

coeftest(Arima1_xy_si)
ggtsdiag(Arima1_xy_si)
checkresiduals(Arima1_xy_si)
plot(rstandard(Arima1_xy_si))
resid_Arima1_xy<-residuals(Arima1_xy_si)
polin1_XY<-coefs2poly(Arima1_xy_si)
(outliersX1<-locate.outliers(resid_Arima1_xy,polin1_XY))
tso(residuals(Arima1_xy_si))

### b) Modelo bivariado con intervención
Arima2_xy_ci<-arimax(gasolinaY_BC,order = c(0,1,1),xreg = data.frame(brentX1_BC,LSef_X1))
BIC(Arima2_xy_ci)
coeftest(Arima2_xy_ci)
ggtsdiag(Arima2_xy_ci)
checkresiduals(Arima2_xy_ci)
plot(rstandard(Arima2_xy_ci))
resid_Arima2_xy<-residuals(Arima2_xy_ci)
polin2_XY<-coefs2poly(Arima2_xy_ci)
(outliersX2<-locate.outliers(resid_Arima2_xy,polin2_XY))
tso(residuals(Arima2_xy_ci))


### c). Modelo bivariado con intervención y efecto
Arima3_xy_ci<-arimax(gasolinaY_BC,order = c(0,1,1),xreg = data.frame(brentX1_BC),xtrans=LSef_X1,transfer = list(c(1,0)))
AIC(Arima3_xy_ci)
coeftest(Arima3_xy_ci)
ggtsdiag(Arima3_xy_ci)
checkresiduals(Arima3_xy_ci)
plot(rstandard(Arima3_xy_ci))
resid_Arima3_xy<-residuals(Arima3_xy_ci)
polin3_XY<-coefs2poly(Arima3_xy_ci)
(outliersX3<-locate.outliers(resid_Arima3_xy,polin3_XY))
tso(residuals(Arima3_xy_ci))
AIC(Arima1_xy_si,Arima2_xy_ci,Arima3_xy_ci)

