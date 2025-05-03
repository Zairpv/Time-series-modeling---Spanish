rm(list=ls())
options(max.print = 10000)
setwd("C:/Users/zaira/Documents/Doctorado PCF/3. Cursos/2. Primavera 2021/3. Series de tiempo/3. Tareas/Tarea 10")

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

### 1. Leer datos. 
datos<- read.csv("DATOS_HGO.csv",header = TRUE)
View(datos)


#Definir variables
PFM_Y_ts <- ts(datos$PFM,start = c(2014,01),frequency = 12)
IF_X1_ts <- ts(datos$IF,start = c(2014,01),frequency = 12)

### 1. Análisis exploratorio de Y1 y X1 
basicStats(PFM_Y_ts)
basicStats(IF_X1_ts)


## a)  Análisis exploratorio de Y
plot(PFM_Y_ts,xlab = "Tiempo",ylab=expression("PFM "~"(m"^"3"~")"),main = "a)")
grid()
par(mfrow=c(1,3),mar=c(5,4,3,1),cex=0.6)
hist(PFM_Y_ts,main = "b)",xlab =expression("PFM "~"(m"^"3"~")"),ylab = "Frecuencia")
plot(density(PFM_Y_ts), main = "c)",ylab = "Densidad")
qqPlot(PFM_Y_ts,dist="norm",ylab = expression("PFM "~"(m"^"3"~")"),main = "d)")
dev.off() 
jarque.bera.test(PFM_Y_ts)

## b) Transformación de Yt 
#Transformación logaritmica
jarque.bera.test(log(PFM_Y_ts))
#Transformación BoxCox
(lambdaY<-BoxCox.lambda(PFM_Y_ts))
PFMY_BC<-BoxCox(PFM_Y_ts,lambda = lambdaY)
jarque.bera.test(PFMY_BC)
#Comparación entre transformaciones
par(mfrow=c(2,3),mar=c(5,4,3,1),cex=0.6)
hist(PFM_Y_ts,main = "Original",xlab = expression("PFM "~"(m"^"3"~")"),ylab = "Frecuencia")
hist(log(PFM_Y_ts),main = "Logarítmica",xlab = expression("PFM "~"(m"^"3"~")"),ylab = "Frecuencia")
hist(PFMY_BC,main = expression(paste("BoxCox (", lambda, " = -0.4164)")),xlab = "PFM (BoxCox)",ylab = "Frecuencia")
qqPlot(PFM_Y_ts,dist="norm",ylab = expression("PFM "~"(m"^"3"~")"),main = "")
qqPlot(log(PFM_Y_ts),dist="norm",ylab = expression("PFM "~"(m"^"3"~")"),main = "")
qqPlot(PFMY_BC,dist="norm",ylab = expression("PFM "~"(m"^"3"~")"),main = "")
dev.off()

## c) Análisis exploratorio de X1: num. árboles infectados
plot(IF_X1_ts,xlab = "Tiempo",ylab = expression("Incendios forestales "~"(m"^"2"~")"),main = "a)")
grid()
par(mfrow=c(1,3),mar=c(5,4,3,1),cex=0.6)
hist(IF_X1_ts,main = "b)",xlab = expression("Incendios forestales "~"(m"^"2"~")"),ylab = "Frecuencia")
plot(density(IF_X1_ts), main = "c)",ylab = "Densidad")
qqPlot(IF_X1_ts,dist="norm",ylab = expression("Incendios forestales "~"(m"^"2"~")"),main = "d)")
dev.off() 
jarque.bera.test(IF_X1_ts)
## d) Transformación BoxCox de Xt
(lambdaX1<-BoxCox.lambda(IF_X1_ts))
IFX1_BC<-BoxCox(IF_X1_ts,lambda = lambdaX1)
hist(IFX1_BC,main = "",xlab = expression("Incendios forestales "~"(m"^"2"~")"),ylab = "Frecuencia")
qqPlot(IFX1_BC,dist="norm",ylab = expression("Incendios forestales "~"(m"^"2"~")"),main = "")
jarque.bera.test(IFX1_BC)


### 2. Identificación de la señal y nivel de la serie Y
## a) Señal de la serie
par(mfrow=c(2,1),mar=c(4.5,4.6,1.7,1),cex=0.6)
tt_pfm<-1:length(PFMY_BC)
fit_pfm<-ts(loess(PFMY_BC~tt_pfm,span = 0.2)$fitted,start = 2014,frequency=12)
plot.ts(fit_pfm,type="l",main="a)",ylab=expression(paste("PFM (", lambda, " = -0.4164)")))
grid()
lines(PFMY_BC,col="Red")

## b) Nivel de la serie
fitlevel_pfm<-lm(PFMY_BC~1)
summary(fitlevel_pfm)
plot.ts(PFMY_BC,main="b)",ylab=expression(paste("PFM", lambda, " = -0.4164)")))
lines(ts(fitted(fitlevel_pfm),start = 2014,frequency = 12),col="red")
dev.off()


### 3. Cambios estructurales Y
## a) Cambios de nivel
PFMY_BC_brk<-breakpoints(PFMY_BC~1,h=0.1)
summary(PFMY_BC_brk)
par(mfrow=c(1,1),cex=0.9,mar=c(5,4,3,1))
plot(PFMY_BC_brk) #BIC=6 brk
breakdates(PFMY_BC_brk,breaks = 6)
plot(PFMY_BC,col="dark gray",ylab=expression(paste("PMF (", lambda, " = -0.4164)")),xlab="Tiempo")
lines(fitted(PFMY_BC_brk,breaks = 6),col=4)
lines(confint(PFMY_BC_brk,breaks = 6))
coef(PFMY_BC_brk,breaks = 6)

## b) Cambios de tendencia
trend_fit_pfm<-lm(PFMY_BC~tt_pfm)
summary(trend_fit_pfm)
plot(PFMY_BC,ylab=expression(paste("PFM (", lambda, " = -0.4164)")),xlab="Tiempo")
lines(ts(fitted(trend_fit_pfm),start = 2014,frequency = 12),col="red")
pfmBC_brk_trend<-breakpoints(PFMY_BC~tt_pfm,h=0.1)
summary(pfmBC_brk_trend)
plot(pfmBC_brk_trend) #6 Brk
breakdates(pfmBC_brk_trend,breaks = 6)
coef(pfmBC_brk_trend,breaks = 6)
plot(PFMY_BC,col="dark gray",ylab=expression(paste("PFM ", lambda, " = -0.4164)")),xlab="Tiempo")
lines(fitted(pfmBC_brk_trend,breaks = 6),col=4)
lines(confint(pfmBC_brk_trend,breaks = 6))
coef(pfmBC_brk_trend,breaks = 6)


### 4. Raiz unitaria de las series Yt, X1, X2
##Serie Y
par(mfrow=c(2,1),mar=c(5.1,4.9,1.7,1),cex=0.6)
plot(PFMY_BC,ylab = expression(paste("PFM (", lambda, " = -0.4164)")))
plot(diff(PFMY_BC),ylab = expression(paste("PFM (", lambda, " = -0.4164)")),sub="d=1")
acf(PFMY_BC,main="a)")
acf(diff(PFMY_BC),main="b)",sub="d=1")
mod1<-ar(diff(PFMY_BC),method = "mle")
mod1$order #Fue = 3
#Hay constante? es significativa?
t_Y=seq(2:length(PFMY_BC))
mod2<-lm(diff(PFMY_BC)~t_Y)
summary(mod2) #No es significativa
#Ajustar a 3 lags
adfTest(PFMY_BC,lags = 3,type = "nc")
adfTest(diff(PFMY_BC),lags = 3,type = "nc")
#La serie transformada tiene raiz unitaria y por lo tanto necesita diferenciarse
ndiffs(PFMY_BC)
dev.off()

##Serie X1: arboles infectados
par(mfrow=c(2,1),mar=c(5.1,4.9,1.7,1),cex=0.6)
plot(IFX1_BC,ylab = expression(paste("árboles infectados (", lambda, " = -0.047)")))
plot(diff(IFX1_BC),ylab = expression(paste("Sup. afect (", lambda, " = -0.047)")),sub="d=1")
acf(IFX1_BC,main="a)")
acf(diff(IFX1_BC),main="b)")
mod3<-ar(diff(IFX1_BC),method = "mle")
mod3$order #Fue = 1
#Hay constante? es significativa?
t_X1=seq(2:length(IFX1_BC))
mod4<-lm(diff(IFX1_BC)~t_X1)
summary(mod4) #No es significativa
#Ajustar a 1 lag
adfTest(IFX1_BC,lags = 1,type = "nc")
adfTest(diff(IFX1_BC),lags = 1,type = "nc")
dev.off()
#La serie transformada tiene raiz unitaria y por lo tanto necesita diferenciarse


### 5. Correlación entre pares de variables 
#La correlacion debe ser entre variables estacionarias
par(mfrow=c(1,1),mar=c(5.1,4.9,3,1),cex=0.6)
ccf(diff(IFX1_BC),diff(PFMY_BC),ylab="Cross correlation",main="",type = "correlation")
dev.off()


### 6. Test de causalidad Granger -- Ho: X no causa a Y
grangertest(diff(IFX1_BC),diff(PFMY_BC),order = 3)

#¿Existe una relación inversa?
grangertest(diff(PFMY_BC),diff(IFX1_BC))
#No existe


### 7. Preblanqueamiento X1
acf(diff(IFX1_BC))
pacf(diff(IFX1_BC))
acf2(diff(IFX1_BC),main = "Serie: Incendios forestales (d=1)")

## a) Modelos arima 
Arima1_x1_model<- arima(IFX1_BC, order=c(0,1,0),seasonal = list(order = c(0,1,1), period = 12),method="ML")
coeftest(Arima1_x1_model)

## b) Diagnostico del modelo
ggtsdiag(Arima1_x1_model)
Box.test(Arima1_x1_model$residuals) # Test de Box-Pierce #Independencia
Box.test(Arima1_x1_model$residuals, type="Ljung-Box") #Independencia
jarque.bera.test(Arima1_x1_model$residuals) #Normalidad
par(mfrow=c(1,2),cex=0.9,mar=c(5,4,3,1))
hist(Arima1_x1_model$residuals, main="a)",xlab = "Residuales",ylab = "Frecuencia")
qqPlot(Arima1_x1_model$residuals,distribution = "norm",main = "b)",ylab = "Residuales")
dev.off()

## c) Una forma de preblanqueamiento general
prewhiten(x=diff(IFX1_BC),y=diff(PFMY_BC))
#Nota: la gráfica es la misma que la obtenida en el paso 5. 


## d) Identificación de outliers
tso(IFX1_BC,types = c("AO","LS","TC"))
coefhat_LSX1=-0.166
plot(tso(IFX1_BC,types = c("AO","LS","TC")))
#Se identifico un outlier LS (cambio de nivel), 
#¿y en residuales?
res_Arima1_x<-residuals(Arima1_x1_model)
polin_X1<-coefs2poly(Arima1_x1_model)
(outliersX1<-locate.outliers(res_Arima1_x,polin_X1)) 
#Nota. esta función solo se utiliza con residuales
plot(tso(res_Arima1_x))

## e) Incorporar intervención al modelo de X1 
#Crear vector binario
LS_X1<-outliers("LS",71)
xreg_lsx1<-outliers.effects(LS_X1,length(IFX1_BC))

ls_efec<-xreg_lsx1*coefhat_LSX1
ls_effect_ts <- ts(ls_efec, frequency = frequency(IFX1_BC), start = start(IFX1_BC))
if_wots<-IFX1_BC-ls_effect_ts
plot(cbind(IFX1_BC, if_wots, ls_effect_ts))
plot(IFX1_BC, type ='l', ylab = "IF",col="light gray")
lines(if_wots, col = 'red', lty = 3, type ='l')

#Una revisión rápida de los residuos de la serie de tiempo
#sin el efecto de cambio transitorio confirma la validez del mod ARIMA (0,1,0)
sarima(if_wots, p=0, d=1, q=0,P=0,D=1,Q=1,S=12)

arima2_x1_model <- arimax(IFX1_BC,
                       order = c(0,1,0),
                       seasonal = list(order = c(0,1,1), period = 12),
                       xtransf = data.frame(LS = (1*(seq(IFX1_BC) == 72))),                                                                                      
                       transfer = list(c(1,0)),
                       method='ML')
summary(arima2_x1_model)
coeftest(arima2_x1_model)
plot(IFX1_BC)
lines(fitted(arima2_x1_model), col = 'blue')

AIC(Arima1_x1_model,arima2_x1_model)
BIC(Arima1_x1_model,arima2_x1_model)

## f) Diagnostico del modelo
ggtsdiag(arima2_x1_model)
Box.test(arima2_x1_model$residuals) # Test de Box-Pierce #Independencia
Box.test(arima2_x1_model$residuals, type="Ljung-Box") #Independencia
hist(arima2_x1_model$residuals, main="a)",xlab = "Residuales",ylab = "Frecuencia")
qqPlot(arima2_x1_model$residuals,distribution = "norm",main = "b)",ylab = "Residuales")
dev.off()
plot(tso(arima2_x1_model$residuals))


### 8. Preblanqueamiento Yt
## a) Modelo arima univariado
acf(diff(PFMY_BC))
pacf(diff(PFMY_BC))
acf2(diff(PFMY_BC))

#Propuestas
Arima1_y_model<- arima(PFMY_BC, order=c(1,1,0),seasonal = list(order = c(1,1,0), period = 12),method="ML")
summary(Arima1_y_model)
coeftest(Arima1_y_model)

## b) Diagnostico
ggtsdiag(Arima1_y_model)
Box.test(Arima1_y_model$residuals) # Test de Box-Pierce #Independencia
Box.test(Arima1_y_model$residuals, type="Ljung-Box") #Independencia
jarque.bera.test(Arima1_y_model$residuals) #Normalidad
par(mfrow=c(1,2),cex=0.9,mar=c(5,4,3,1))
hist(Arima1_y_model$residuals, main="a)",xlab = "Residuales",ylab = "Frecuencia")
qqPlot(Arima1_y_model$residuals,distribution = "norm",main = "b)",ylab = "Residuales")
dev.off()

## c) Identificación de outliers en Yt
tso(PFMY_BC)
plot(tso(PFMY_BC))
#No se identificaron outliers, ¿y en residuales?
resid_Arima1_y<-residuals(Arima1_y_model)
polin_y<-coefs2poly(Arima1_y_model)
(outliersY1<-locate.outliers(resid_Arima1_y,polin_y))
plot(tso(resid_Arima1_y))

## e) Incorporar intervención al modelo de X1 
#Crear vector binario
AO_Y<-outliers("AO",13)
xreg_AOY1<-outliers.effects(AO_Y,length(PFMY_BC))
coefhat_AOY=-0.007739
AO_efec<-xreg_AOY1*coefhat_AOY
AO_effect_ts <- ts(AO_efec, frequency = frequency(PFMY_BC), start = start(PFMY_BC))
pfm_wots<-PFMY_BC-AO_effect_ts
plot(cbind(PFMY_BC, pfm_wots, AO_effect_ts))
plot(PFMY_BC, type ='l', ylab = "PFM",col="light gray")
lines(pfm_wots, col = 'red', lty = 3, type ='l')

#Una revisión rápida de los residuos de la serie de tiempo
#sin el efecto de cambio transitorio confirma la validez del mod ARIMA (0,1,0)
sarima(pfm_wots, p=1, d=1, q=0,P=1,D=1,Q=0,S=12)

arima2_y_model <- arimax(Arima1_y_model$residuals,
                          order = c(1,1,0),
                          seasonal = list(order = c(1,1,0), period = 12),
                          xreg = AO_efec,                                                                                      
                          method='ML')
summary(arima2_y_model)
coeftest(arima2_y_model)
plot(Arima1_y_model$residuals)
lines(fitted(arima2_y_model), col = 'blue')

AIC(Arima1_y_model,arima2_y_model)
BIC(Arima1_y_model,arima2_y_model)


### 9. Correlación cruzada entre series preblanqueadas
## Se utilizaron residuales de modelos X y Y 
#Correlación cruzada entre series preblanqueadas
par(mfrow=c(1,1),cex=0.9,mar=c(5,4,3,1))
ccf(arima2_x1_model$residuals,Arima1_y_model$residuals,ylab="Cross correlation",main="Precio Brent & Precio Gasolina",type = "correlation",na.action = na.omit)
dev.off()

### 10. Función de transferencia

### a). Modelo bivariado sin intervención
Arima1_xy_model_si<- arima(PFMY_BC, order=c(1,1,0),seasonal = list(order = c(1,1,0), period = 12),method="ML",xreg = IFX1_BC)
summary(Arima1_xy_model_si)
coeftest(Arima1_xy_model_si)
BIC(Arima1_xy_model_si)

### b) Modelo bivariado con intervención
Arima2_xy_model_ci<-arima(PFMY_BC, order=c(1,1,0),seasonal = list(order = c(1,1,0), period = 12),method="ML",xreg = data.frame(IFX1_BC,ls_efec))
coeftest(Arima2_xy_model_ci)
ggtsdiag(Arima2_xy_model_ci)


### c). Modelo bivariado con intervención y efecto
Arima3_xy_model_ci<-arima(PFMY_BC, order=c(1,1,0),seasonal = list(order = c(1,1,0), period = 12),method="ML",xreg = IFX1_BC,xtrans=ls_efec,transfer = list(c(1,0)))
summary(Arima3_xy_model_ci)
coeftest(Arima3_xy_model_ci)
ggtsdiag(Arima3_xy_model_ci)
jarque.bera.test(Arima3_xy_model_ci$residuals)

AIC(Arima1_xy_model_si,Arima2_xy_model_ci,Arima3_xy_model_ci)

############# Figuras del articulo
## e) Figura 1
par(mfrow=c(2,1),mar=c(4.5,4.6,1.7,1),cex=0.8)
plot(PFM_Y_ts,xlab = "",ylab=expression("Vol. "~"(m"^"3"~")"),main = "Producción forestal maderable")
grid()
plot(IF_X1_ts,xlab = "",ylab = expression("Infestación "~"(árboles ha"^"-1"~")"),main = "Árboles infectados por muérdago")
grid()

## e) Figura 2
par(mfrow=c(2,1),mar=c(4.5,4.6,1.7,1),cex=0.6)
plot(diff(PFMY_BC),xlab = "",ylab=expression("Vol. "~"(m"^"3"~") ("~lambda~" = -0.4164)"),main = "Producción forestal maderable")
grid()
plot(diff(IFX1_BC),xlab = "",ylab = expression("Infestación "~"(árb ha"^"-1"~") ("~lambda~" = -0.0470)"),main = "Árboles infectados por muérdago")
grid()

## Figura 4
par(mfrow=c(1,1),mar=c(4.5,4.6,1.7,1),cex=0.6)
plot(tso(IFX1_BC,types = c("AO","LS","TC")),cex=0.6)
plot(IFX1_BC, type ='l', ylab =expression("Infestación "~"(árb ha"^"-1"~") ("~lambda~" = -0.0470)"),col="light gray")
grid()
lines(if_wots, col = 'red', lty = 3, type ='l')


########################################
#### Autocorrelación en residuales
#Residuales

#Modelo ajustado en la linea 338. 
ggtsdiag(Arima3_xy_model_ci)
checkresiduals(Arima3_xy_model_ci)

residxy<-residuals(Arima3_xy_model_ci)
acf2(residxy) #Hay presencia de algunos picos en acf y pacf a 2 lags
library(car)
library(lmtest)
residxy<-as.vector(residxy)
durbinWatsonTest(residxy)
#El valor cercano a 2 es indicativo de ausencia de autocorrelación de primer orden
shapiro.test(residxy) #Son normales
Box.test(residxy, type="Ljung-Box") #son independientes

#Definir tiempo
trend=time(residxy)
reg_model<-lm(residxy~trend)
summary(reg_model)
acf2(residuals(reg_model))
#Ho: no hay autocorrelacion 
dwtest(residxy~trend)


################# Tarea 12
#Heterocedasticidad ############ TAREA 12
bptest(reg_model)
# Ho: Los residuales son homocedasticos
# No se rechaza Ho
hist(reg_model$residuals)
plot(density(reg_model$residuals))
shapiro.test(reg_model$residuals)



##Afinar el modelo
adjreg = sarima (residuals(reg_model), 0,0,1, xreg=trend) # This is the adjustment regression with MA(1) residuals
adjreg2 = sarima (residuals(reg_model), 2,0,3, xreg=trend) # This is the adjustment regression with MA(1) residuals
##Este modelo es mas cercano a ruido blanco

acf2(resid(adjreg$fit))
acf2(resid(adjreg2$fit))

resid_adjreg2<-resid(adjreg2$fit)
resid_adjreg2<-as.vector(resid_adjreg2)
durbinWatsonTest(resid_adjreg2)
shapiro.test(resid_adjreg2) #son normales
Box.test(resid_adjreg2, type="Ljung-Box") #Son independientes


#Poner el modelo sarima en objeto arima
adjreg22 = arima(residuals(reg_model), order=c(2,0,3), xreg=trend) # This is the adjustment regression with MA(1) residuals
checkresiduals(adjreg22) 
ggtsdiag(adjreg22)
arch.test(adjreg22) #probar homocedasticidad , h0 es homocedast..  
res2=residuals(adjreg22)^2
plot(res2,type="l")
acf2(res2,ci.type="ma")
