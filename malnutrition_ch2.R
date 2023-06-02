rm(list=ls())

library(PrevMap)
Ghana.bndrs <- read.csv("Ghana_bndrs.csv")

maln <- read.csv("malnutrition.csv")

ID.coords <- create.ID.coords(maln,~utm_x+utm_y)
HAZ.avg <- tapply(maln$HAZ,ID.coords,mean)
coords <- unique(maln[,c("utm_x","utm_y")])

par(mfrow=c(2,2),mar=c(4,4,2,2))
plot(Ghana.bndrs, type = "l", asp = 1, xlab = "", ylab = "",
     main=c("(a)"))
points(coords,pch=20,cex=0.5)
plot(HAZ ~ age, data = maln,main="(b)",xlab="Age",pch=20,cex=0.5)
plot(HAZ ~ factor(edu), data = maln,main=c("(c)"),xlab="Maternal education")
plot(HAZ ~ factor(wealth), data = maln,main="(d)",xlab="Wealth index")
par(mfrow=c(1,1))


library(splines)
max.vec <- function(val,x) sapply(x,function(i) max(0,i-val))
lm.fit <- lm(HAZ ~ age+I(max.vec(1,age))+I(max.vec(2,age))+edu+wealth,data=maln)
summary(lm.fit)

beta.hat <- coef(lm.fit)
age.set <- seq(0,5,length=1000)
age.vars <- cbind(age.set,max.vec(1,age.set),max.vec(2,age.set))
broken.sticks <- as.numeric(age.vars%*%beta.hat[2:4])
std.errors <- sqrt(diag(age.vars%*%vcov(lm.fit)[2:4,2:4]%*%t(age.vars)))
ci.95 <- cbind(broken.sticks-qnorm(0.975)*std.errors,
               broken.sticks+qnorm(0.975)*std.errors)

matplot(age.set,cbind(ci.95,broken.sticks),type="l",xlab="Age (years)",
        ylab="",
        lty=c("dashed","dashed","solid"),col = 1)

