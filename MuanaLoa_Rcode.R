rm(list=ls())
mauna <- read.csv("mauna_loa.csv",sep=",",dec=".")
plot(mauna$X,mauna$CO2,xlab = "Time", 
     ylab=expression(CO[2] ~ "(ppm)"),type="l")

par(mfrow=c(4,4))
lag.plot(mauna$CO2,lags = 12)

par(mfrow=c(1,1))
f <- 1/12
A <- 2
phi <- pi/2
t <- seq(0,24,length=48*2)
mu.t <- A*sin(2*pi*f*t+phi)
plot(t,mu.t,type="l",xlab="t",ylab=expression(E(Y[t])))
abline(h=0,lty="dashed")

segments(0,0,0,2,col=2)
text(1,1,"A",col=2)

segments(0,0,12,0,col=3)
text(6,-0.25,"1/f",col=3)

f <- 1/12

mauna$time <- mauna$X
mauna$sin.t <- sin(2*pi*f*mauna$time)
mauna$cos.t <- cos(2*pi*f*mauna$time)
lm.fit <- lm(CO2 ~ sin.t+cos.t+time,data=mauna)
plot(mauna$time,mauna$CO2,xlab = "Time", ylab=expression(CO[2] ~ "(ppm)"),type="l")
lines(predict(lm.fit),col=2)

plot(lm.fit$fitted.values,residuals(lm.fit),xlab="Fitted values",ylab="Residuals")
abline(h=0,lty="dashed")

acf(residuals(lm.fit))
plot(mauna$time,residuals(lm.fit),xlab="Time",ylab="Residuals",type="l")
abline(h=0,lty="dashed")

f2 <- 1/6
mauna$sin.t2 <- sin(2*pi*f2*mauna$time)
mauna$cos.t2 <- cos(2*pi*f2*mauna$time)
lm.fit2 <- lm(CO2 ~ sin.t+cos.t+
                sin.t2+cos.t2+time
              ,data=mauna)
summary(lm.fit2)

plot(mauna$time,residuals(lm.fit2),xlab="Time",ylab="Residuals",type="l")
abline(h=0,lty="dashed")

beta.season <- coef(lm.fit2)[2:5]
curve(beta.season[1]*sin(2*pi*f*x)+beta.season[2]*cos(2*pi*f*x)+
        beta.season[3]*sin(2*pi*f2*x)+beta.season[4]*cos(2*pi*f2*x),
      xlim=c(0,48),xlab = "t", ylab="",main="Seasonal trend")


acf(residuals(lm.fit2),
    main="Correlogram of the residuals")
ar.fit <- arima(mauna$CO2,
                xreg=mauna[,c("sin.t","cos.t",
                              "sin.t2","cos.t2",
                              "time")],order=c(1,0,0),
                method="ML")
acf(ar.fit$residuals,lag.max = 40)

plot(mauna$X,mauna$CO2,xlab = "Time", ylab=expression(CO[2] ~ "(ppm)"),type="l")
predicted.CO2 <- mauna$CO2-ar.fit$residuals
lines(predicted.CO2,col=2)
plot(ar.fit$residuals,type="l", xlab="Time",
     ylab="",main="AR(1) Residuals")
abline(h=0,lty="dashed")

acf(ar.fit$residuals,main="Series AR(1) residuals")

