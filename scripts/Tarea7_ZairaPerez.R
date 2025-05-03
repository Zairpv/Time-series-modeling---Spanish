  #TAREA 7
  #Presenta: Zaira Rosario Pérez Vázquez
  
  #Librerias
  library(ggplot2)
  library(forecast)
  library(astsa)
  library(fUnitRoots)
  library(lmtest)
  library(strucchange)
  library(reshape)
  library(fBasics)
  library(FitARMA)
  library(Rmisc)
  library(tidyverse)
  library(magrittr)
  library(ggfortify)
  library(TSA)
  library(dygraphs)
  library(dplyr)
  library(car)
  library(tseries)
  
  rm(list=ls())
  options(max.print = 10000)
  setwd("C:/Users/zaira/Documents/Doctorado PCF/3. Cursos/2. Primavera 2021/3. Series de tiempo/3. Tareas/Tarea 7")
  
  nacimientos <- read.csv("londres.csv", header=TRUE)
  View(nacimientos)
  
  ####1. Exploración de la serie de tiempo
  nacimientos_tot_ts<-ts(nacimientos$boys + nacimientos$girls, frequency = 1 , 
                       start = 1629)
  exceso_boys<-ts(nacimientos$excess,frequency = 1 , 
                  start = 1629)
  
  #a) Nacimientos por género
  nacimientos_rs <- melt(nacimientos[1:3], id = c("year"))
  F1<-ggplot(data = nacimientos_rs, aes(x = year)) + geom_line(aes(y = value, colour = variable)) +
    scale_colour_manual(values = c("blue", "red"))
  #b) Nacimientos totales (boys + girls)
  F2<-autoplot(nacimientos_tot_ts,xlab = "",ylab = "births")
  #c) Exceso de nacimientos masculinos
  F3<-autoplot(exceso_boys,xlab = "year",ylab = "excess boys")
  multiplot(F1,F2,F3)
  
  basicStats(exceso_boys)
  par(mfrow=c(1,3),cex=0.9)
  hist(exceso_boys,main = "a)",xlab = "Nacimientos",ylab = "Frecuencia")
  plot(density(exceso_boys), main = "b)",ylab = "Densidad")
  qqPlot(exceso_boys,dist="norm",ylab = "Nacimientos",main = "c)")
  dev.off() 
  jarque.bera.test(exceso_boys)
  
  ### 2. Analisis de datos atipicos
  bp<-boxplot(exceso_boys)
  bp$out
  
  library(tsoutliers)
  ?tso
  (outliers_excesoboys<-tso(exceso_boys))
  plot(outliers_excesoboys)
  #¿Es igual con la serie transformada a logaritmo?
  exceso_boys_log<-log(exceso_boys)
  autoplot(exceso_boys_log)
  (outliers_excesoboys_log<-tso(exceso_boys_log))
  plot(outliers_excesoboys_log)
  
  
  ### 3. Cambios estructurales
  # a) Identificación de la señal de la serie
  tt_exceso <- 1:length(exceso_boys)
  fit_exceso <- ts(loess(exceso_boys ~ tt_exceso, span = .2)$fitted, start = 1629, frequency = 1)
  
  plot.ts(fit_exceso, type='l',ylab="values")
  grid()
  lines(exceso_boys, col = "red") 
  
  # b) Estimación del nivel de la serie (y=bo+e)
  fit_level_exceso<-lm(exceso_boys~1)
  summary(fit_level_exceso)
  plot.ts(exceso_boys, main="",ylab="values")
  grid()
  lines(ts(fitted(fit_level_exceso), start=1629, frequency = 1), col = "red")
  
  # c) Cambios estructurales de tendencia
  excesos_brk_trend <- breakpoints(exceso_boys ~ tt_exceso, h = 0.1)
  summary(excesos_brk_trend)
  plot(excesos_brk_trend)
  #Nota: De acuerdo con BIC, no existen cambios estructurales de tendencia
  
  
  # d) Identificación de cambios estructurales de nivel
  #Breakpoints
  excesos_brk <- breakpoints(exceso_boys ~ 1, h = 0.1)
  summary(excesos_brk)
  plot(excesos_brk)
  #Nota: De acuerdo con  BIC, aproximadamente existe 1 cambios de nivel
  breakdates(excesos_brk, breaks = 1)
  plot(exceso_boys,col="dark gray",ylab="values")
  lines(fitted(excesos_brk, breaks = 1), col = 4)
  lines(confint(excesos_brk, breaks = 1))
  coef(excesos_brk, breaks = 1)
  
  # e) Tasa del cambio de nivel
  fitted(excesos_brk)[1]
  fitted(excesos_brk)[length(exceso_boys)]
  #La tasa de cambio fue de 0.08 a 0.05
  
  # f) Significancia en las medias 
  (break_date<-breakdates(excesos_brk))
  pre_interv<-window(exceso_boys,end=break_date)
  post_interv<-window(exceso_boys,start=break_date+1)
  t.test(pre_interv,post_interv)
  
  ### 4. Prueba de raiz unitaria
  autoplot(exceso_boys)
  autoplot(exceso_boys_log)
  acf(exceso_boys) 
  acf(exceso_boys_log) #Son estacionarias, no se necesita diferenciar
  
  mod1=ar(exceso_boys_log,method="mle")
  mod1$order
  # ¿Hay pendiente? ¿Cuánto vale la constante? 
  t_exceso=seq(1:length(exceso_boys_log))
  mod2=lm((exceso_boys_log)~t_exceso) 
  summary(mod2) #La constante es significativa
  adfTest(exceso_boys_log,lags=1,type=c("c")) #serie estacionaria
  ndiffs(exceso_boys_log)
  
  
  ### 5. Funciones acf y pacf
  acf(exceso_boys_log)
  pacf(exceso_boys_log)
  
  acf2(exceso_boys_log)
  #S=10
  
  
  ### 6. Modelo ARIMA sin intervencion
  #Modelo ARIMA (p,d,q)(P,D,Q)(S) #S=10
  
  #a) ARIMA1(1, 0, 0) (1, 0, 0) (10). 
  ARIMA1 <- Arima(exceso_boys_log, order = c(1,0,0), 
                   seasonal = list(order = c(1,0,0), 
                                   period = 10),include.mean = TRUE)
  summary(ARIMA1)
  coeftest(model_1)
  
  # b) ARIMA2(0, 0, 1) (0, 0, 1) (10). 
  ARIMA2 <- Arima(exceso_boys_log, order = c(0,0,1), 
                  seasonal = list(order = c(0,0,1), 
                                  period = 10),include.mean = TRUE)
  summary(ARIMA2)
  coeftest(model_2)
  
  # c) ARIMA3(1, 0, 1) (1, 0, 1) (10). 
  ARIMA3 <- Arima(exceso_boys_log, order = c(1,0,1), 
                  seasonal = list(order = c(1,0,1), 
                                  period = 10),include.mean = TRUE)
  summary(ARIMA3)
  coeftest(ARIMA3)
  
  AIC(ARIMA1,ARIMA2,ARIMA3)
  BIC(ARIMA1,ARIMA2,ARIMA3)
  #Mejor: ARIMA1
  
  
  ### 7. Creación de variable regresora de intervencion
  nivel <- c(rep(0, excesos_brk$breakpoints), 
             rep(1, length(exceso_boys_log) - excesos_brk$breakpoints))
  plot(data.frame(nacimientos$year,nivel),pch=20,ylim=c(0,1),ylab="Value",xlab="Year")
  
  
  ### 8. Modelo ARIMA con intervención
  #ARIMA4 (1, 0, 0) (1, 0, 0) (10).
  ARIMA4 <- Arima(exceso_boys_log, order = c(1,0,0), 
                   seasonal = list(order = c(1,0,0), period = 10), 
                   xreg = nivel, include.mean = TRUE)
  summary(ARIMA4)
  coeftest(ARIMA4)
  
  #ARIMA5 (0, 0, 1) (0, 0, 1) (10).
  ARIMA5 <- Arima(exceso_boys_log, order = c(0,0,1), 
                  seasonal = list(order = c(0,0,1), period = 10), 
                  xreg = nivel, include.mean = TRUE)
  summary(ARIMA5)
  coeftest(ARIMA5)
  
  #ARIMA6 (1, 0, 1) (1, 0, 1) (10).
  ARIMA6 <- Arima(exceso_boys_log, order = c(1,0,1), 
                  seasonal = list(order = c(1,0,1), period = 10), 
                  xreg = nivel, include.mean = TRUE)
  summary(ARIMA6)
  coeftest(ARIMA6)
  
  AIC(ARIMA4,ARIMA5,ARIMA6)
  BIC(ARIMA4,ARIMA5,ARIMA6)
  #Mejor: ARIMA4
  
  
  ### 9. Diagnóstico del modelo
  #¿Comparando todos?
  AIC(ARIMA1,ARIMA2,ARIMA3,ARIMA4,ARIMA5,ARIMA6)
  BIC(ARIMA1,ARIMA2,ARIMA3,ARIMA4,ARIMA5,ARIMA6)
  #Mejor: ARIMA4
  #Nota. Los valores de AIC y BIC disminuyeron al incorporar 
  #... la variable regresora de intervención
  
  ggtsdiag(ARIMA4)
  Box.test(ARIMA4$residuals) # Test de Box-Pierce #Independencia
  Box.test(ARIMA4$residuals, type="Ljung-Box") #Independencia
  jarque.bera.test(ARIMA4$residuals) #Normalidad
  par(mfrow=c(1,2),cex=0.7)
  hist(ARIMA4$residuals, main="a)",xlab = "Residuales",ylab = "Frecuencia")
  qqPlot(ARIMA4$residuals,distribution = "norm",main = "b)",ylab = "Residuales")
  dev.off()
  
  
  ### 10. Predicción del modelo
  pred1<-forecast(ARIMA4, h = 5, xreg = rep(1, 5))
  F1<-autoplot(pred1,main="a)")
  
  pred2<-forecast(ARIMA4, h = 10, xreg = rep(1, 10))
  F2<-autoplot(pred2,main="b)")
  
  pred3<-forecast(ARIMA4, h = 20, xreg = rep(1, 20))
  F3<-autoplot(pred3,main="c)",xlab="Year")
  
  multiplot(F1,F2,F3)
  

