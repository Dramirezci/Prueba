#Taller 2 - Dasometria 2023-2s

#Cargar datos y revisar datos
bd <- read.csv2("taller-2_Dasometria.csv")
str(bd)
summary(bd)

#Grafico general de los datos: como hay que seleccionar el 80% esto se hace despues de la seleccion
#windows()
#plot(bd$DAP_cm,bd$H_m, xlab="DAP (cm)", ylab="Altura (m)")
#Hay 3 ptos que no se si seran outliers pero llaman la atencion que tienen DAP de los mas grandes(mayores a 40), pero altura de 26

#############################PRIMER PUNTO################################
# Selección del 80% de los datos
set.seed(99)
muestra <- sample(1:nrow(bd),round(nrow(bd)*0.8))
muestra <- muestra[order(muestra)]

#Base de datos con el 80% de los datos
bd.m <- bd[muestra,]
summary(bd.m)
str(bd.m)

#Base de datos para la validacion
bd.v <- bd[-muestra,]
range(bd.v$DAP_cm)
range(bd.m$DAP_cm)
#Con range miro que el rango del diametro este entre los datos de validacion
#Si eso no me da, entonces cambio la semilla

###################SEGUNDO PUNTO#########################
#Grafico general de los datos, para ver el comportamiento de los mismos.
windows()
plot(bd.m$DAP_cm,bd.m$H_m, xlab="DAP (cm)", ylab="Altura (m)")

#Evaluacion de modelos:

#diferencias entre modelos lineales y no lineales:
#Los modelos no lineales no tienen ecuacion definida, toca escribirla
#es decir, se define el modelo no lineal específico que se va a usar. En lugar de
#indicar y ~ x, se debe escribir el modelo completo
#Se usa la función nls ("non-linear squares")

#Modelo lineal: H=b0+b1DAP
ml <- lm(H_m~DAP_cm,data=bd.m)
summary(ml)
anova(ml)                       #DAP*** significativo
#Ecuacion modelo lineal: H=2.45+0.71*DAP
#Valor F: 1762***
#R cuadrado ajustado: 0.8374

#Graficos en los que reviso homocedasticidad y normalidad
windows()
par(mfrow=c(2,2))
plot(ml)

#Supuestos:
#Normalidad de los residuales
shapiro.test(residuals(ml))
#W=0.843 Como el p-value es menor a 0.05 los residuales no siguen una distribucion normal

#Prueba Kolmogorov:solo se hace si shapiro no dio normalidad, porque shapiro es la mas exigente.
ks.test(residuals(ml),"pnorm",mean(residuals(ml)),sd(residuals(ml)))
#D:0.14  Como el p-value es menor a 0.05 los residuales no siguen una distribucion normal

#Modelo cuadratico: H=b0+b1DAP+b2(DAP)^2

mc <- lm(H_m~DAP_cm+I(DAP_cm^2),data=bd.m)
summary(mc)
anova(mc)                                 #DAP y I(DAP^2)*** significativos
#Ecuacion modelo cuadratico: H=1.16+1.13*DAP-0.013DAP^2
#Valor F:1454***
#R cuadrado ajustado:0.8947

#Grafico para revisar homocedasticidad y normalidad:
windows()
par(mfrow=c(2,2))
plot(mc)

#Supuestos:
#Normalidad de los residuales:
shapiro.test(residuals(mc))
#W:0.88*** Como el p-value es menor a 0.05 los residuales no siguen una distribucion normal

#Prueba Kolmogorov:
ks.test(residuals(mc),"pnorm",mean(residuals(mc)),sd(residuals(mc)))
#D:0.15 Como el p-value es menor a 0.05 los residuales no siguen una distribucion normal

#Modelo exponencial: ln(H)=b0+b1*ln(DAP) esta es la ecuacion linealizada (modelo logaritmico)
#no se puede usar AIC en este modelo para compararlo con los otros modelos
#Porque la variable x no esta en las mismas unidades que las de los otros
#Muchas veces este es mejor para hacer Diametro vs altura

mexp <- lm(log(H_m)~log(DAP_cm),data=bd.m)
summary(mexp)
anova(mexp)

#Factor de correccion:
fc <- sigma(mexp)^2/2       #error estandar de los residuales que me da el anova...La nueva pendiente es el factor de correccion
b1 <- fc+(coef(mexp)[1]) #Intercepto

#Ecuacion del modelo: para obtener el valor de H se saca exponencial a ambos lados
#exp(ln(H))= exp(0.71+0.74*log(DAP))

#Graficos para revisar normalidad y homocedasticidad:
windows()
par(mfrow=c(2,2))
plot(mexp)

#Supuestos:
#Normalidad:
shapiro.test(residuals(mexp))
#W:0.98 Como el p-value es menor a 0.05 los residuales no siguen una distribucion normal

#Prueba Kolmogorov:
ks.test(residuals(mexp),"pnorm",mean(residuals(mexp)),sd(residuals(mexp)))
#D=0.05 Como el p-value es mayor a 0.05 los residuales siguen una distribucion normal

#RES: #El error estandar de los residuales que da el summary no lo podemos utilizar para este modelo
rse_mexp <- sqrt((sum((bd.m$H_m-exp(predict(mexp)+fc))^2))/mexp$df.residual)
#Error estandar de los residuales: 1.96

#Michaeles: H = a*DAP/(b + DAP)

#Funcion de inicializacion
m_mm <- nls(H_m ~ SSmicmen(DAP_cm, a, b), data = bd.m)
summary(m_mm)

# Incluyo los parámetros obtenidos para el modelo, para validar que el resultado es igual
m_mm1<- nls(H_m ~ (a * DAP_cm)/(b + DAP_cm),
            star = list(a = 40.2, b = 24.2), data = bd.m)
summary(m_mm1)
#Ecuacion del modelo: H=40.2*DAP/(24.2 + DAP)

#Supuestos:
#Normalidad:
shapiro.test(residuals(m_mm1))
#W:0.86 Como el p-value es menor a 0.05 los residuales no siguen una distribucion normal

#Prueba Kolmogorov:
ks.test(residuals(m_mm1),"pnorm",mean(residuals(m_mm1)),sd(residuals(m_mm1)))
#D=0.15 Como el p-value es mayor a 0.05 los residuales siguen una distribucion normal

#Weibull: H	=	a x [1	-	exp(-b x DAP^c)]

weibull_f <- function(a, b, c, DAP) {a*(1-exp(-b*DAP^c))}

m_weibull <- nls(H_m ~ weibull_f(a,b,c,DAP_cm), start=list(a=38,b=0.05,c=0.9), data=bd.m)
summary(m_weibull) # H = 35.9*(1-exp(-0.05*DAP^0.86))

#Supuestos:
#Normalidad:
shapiro.test(residuals(m_weibull))
#W:0.87 Como el p-value es menor a 0.05 los residuales no siguen una distribucion normal

#Prueba Kolmogorov:
ks.test(residuals(m_weibull),"pnorm",mean(residuals(m_weibull)),sd(residuals(m_weibull)))
#D=0.14 Como el p-value es mayor a 0.05 los residuales siguen una distribucion normal

#Se grafica la linea de tendencia de los 5 modelos:
#Modelo lineal
windows()
plot(bd.m$DAP_cm, bd.m$H_m, xlab="DAP (cm)", ylab="Altura (m)")
abline(ml,lwd=2)
#Modelo cuadratico
hmc <- predict(mc,newdata = list(DAP_cm=seq(0,50,.5)))
lines(seq(0,50,.5),hmc,col=2,lwd=2)
#Modelo exponencial
hmexp <- exp(predict(mexp,newdata=list(DAP_cm=seq(0,50,.5))))*exp(fc)
lines(seq(0,50,.5),hmexp,col=3,lwd=2)
#Michaeles
hmm <- predict(m_mm1,newdata = list(DAP_cm= seq(0,50,.5)))
lines(seq(0,50,.5),hmm,col=4, lwd= 2)
#Weibull
hmw <- predict(m_weibull, newdata = list(DAP_cm= seq(0,50,.5)))
lines(seq(0,50,.5),hmw,col=5, lwd= 2)
legend("bottomright",lty=1,lwd=3,col=1:5,legend = c("ML","MC","MEXP","Michaeles-M","Weibull"))

#PARAMETROS DE DESICION
Modelos <- c("Lineal", "Cuadrático", "Exponencial", "Michaeles-M", "Weibull")
RSE <- c(summary(ml)$sigma, summary(mc)$sigma, rse_mexp,
         summary(m_mm1)$sigma, summary(m_weibull)$sigma)

AIC <- c(AIC(ml), AIC(mc), AIC(mexp), AIC(m_mm1), AIC(m_weibull))

Tabla <- data.frame(Modelos, RSE, AIC)
Tabla$Sel <- c("", "", "X", "", "X")
show(Tabla)

ks.test(residuals(mexp),"pnorm",mean(residuals(mexp)),sd(residuals(mexp)))
#D:0.054 Como el p-value es mayor a 0.05 los residuales siguen una distribucion normal

ks.test(residuals(m_weibull),"pnorm", mean(residuals(m_weibull)), sd(residuals(m_weibull)))
# D=0.14 segun Kolmogorov los residuales del modelo weibull no siguen una distribucion normal

windows(10,5)
par(mfrow = c(1,2))
plot(predict(mexp), residuals(mexp))  #tiende a ser homocedastico
abline(h=0, col=2, lwd=2)
plot(predict(m_weibull), residuals(m_weibull)) #presenta problemas de heterocedasticidad
abline(h=0, col=2, lwd=2)

#A pesar que el error estandar de lo residuales es mayor en el modelo exponencial (1.962)
#comparado con el del modelo Weibull (1.799), el modelo que mejor predice los datos es el exponencial
#ya que cumple con los supuestos de homocedasticidad y normalidad.

##########################TERCER PUNTO######################
#Graficos:
windows(10,5)
par(mfrow=c(1,2))
plot(bd.m$DAP_cm, bd.m$H_m, xlab="DAP (cm)", ylab="Altura (m)", main = "Modelo DAP vs Altura", cex.main=0.8)
#Linea de tendencia
hmexp <- exp(predict(mexp,newdata=list(DAP_cm=seq(0,50,.5))))*exp(fc)
lines(seq(0,50,.5),hmexp,col="purple",lwd=2)
legend("bottomright",lty=1,lwd=3,col="purple",legend = c("H= exp(0.71+0.74*log(10))"))

#Estimados vs Predichos
plot(bd.m$H_m,exp(predict(mexp))*exp(fc),xlab = "Observados", ylab = "Predichos", main = "Observados vs Predichos", cex.main=0.8)
abline(0,1,col="blue",lwd=2)

##########################CUARTO PUNTO######################
#Revisar datos
summary(bd.v)
str(bd.v)
#Observados
obs <- bd.v$H_m
#Alturas predichas con el 20%
hmexp<- exp(predict(mexp,newdata = bd.v)+fc)
#Error:
error <- mean((obs-hmexp)/obs)
#Sesgo %
error*100
#Incertidumbre:
I <- sd(hmexp/obs)
#Incertidumbre %
I*100

#Como el sesgo es negativo el modelo exponencial esta sobrestimando en un 3,49% los datos de altura
#El porcentaje de incertidumbre fue del 22.40%, esto quiere decir, que los valores predichos por el modelo
#podrían variar un 22.40% respecto a los valores reales dentro del rango de los datos. 
