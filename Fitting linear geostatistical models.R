rm(list=ls())
library(PrevMap)

maln <- read.csv("malnutrition.csv")
ghana <- read.csv("Ghana_bndrs.csv")

max.vec <- function(x,c) sapply(x,function(x) max(x,c))

lm.fit <- lm(HAZ ~ age + max.vec(age,2),data=maln)

summary(lm.fit)

library(lme4)

maln$ID.utm <- create.ID.coords(maln,~utm_x+utm_y)
lmer.fit <- lmer(HAZ ~ age + max.vec(age,2)+(1|ID.utm),data=maln)

summary(lmer.fit)

residuals <- data.frame(utm_x=unique(maln$utm_x),
                        utm_y=unique(maln$utm_y),
                        Z=ranef(lmer.fit)$ID.utm[,1])
vari <- variogram(residuals,~Z,~utm_x+utm_y)
plot(vari,type="l")
out.eyefit <- eyefit(vari)

theta.cov.start <- c(phi=19.34,var.ratio=1.4137/0.04)

geo.lm <- linear.model.MLE(HAZ ~ age + max.vec(age,2),
                 coords=~utm_x+utm_y,
                 ID.coords = maln$ID.utm,
                 data=maln,kappa=0.5,
                 fixed.rel.nugget = 0,
                 start.cov.pars = theta.cov.start,
                 method="nlminb")

summary(geo.lm)


# Create a prediction grid

library(splancs)
plot(ghana,type="l",asp=1)
ghana.grid.pred <- gridpts(as.matrix(ghana),xs=10,ys=10)
points(ghana.grid.pred,pch=20,cex=0.1)


predictors.maln <- data.frame(age=rep(5,nrow(ghana.grid.pred)))

maln.pred <- spatial.pred.linear.MLE(geo.lm, predictors = predictors.maln,
                        grid.pred = ghana.grid.pred,
                        scale.predictions = "logit",
                        thresholds = -2,
                        scale.thresholds = "logit",
                        standard.errors = TRUE)

plot(maln.pred,"logit","predictions")
lines(ghana)

plot(maln.pred,"logit","standard.errors")
lines(ghana)

plot(maln.pred,summary="exceedance.prob")
lines(ghana)

