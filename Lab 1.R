rm(list=ls())

library(PrevMap)
library(lme4)

chikwawa <- read.csv("lab_data.csv")

str(chikwawa)

# Number of sampled people per household
table(chikwawa$ID)

attach(chikwawa) # attach variables to use

# Boxplot Hb against sex
boxplot(Hb ~ sex)

# Linear regression of Hb against maternal education (write down then model)
lm.fit <- lm(log(Hb) ~ sex,data=chikwawa)

summary(lm.fit) # interpet the results

# Accounting for clustering effects (write down the model)
lmer.fit <- lmer(log(Hb) ~ sex + (1|ID), data=chikwawa)

summary(lmer.fit) # interpret the results

# Extraction of the household random effects
Z.hat <- ranef(lmer.fit)$ID[,1]

# Coordinates
coords <- unique(chikwawa[,c("web_x","web_y")]) # Never use longitude/latitude!

# Empirical variogram
args(variogram)
residuals <- data.frame(Z.hat=Z.hat, web_x=coords[,1],
                        web_y=coords[,2])
EV <- variogram(residuals,var.name=~Z.hat,coords=~web_x+web_y)
plot(EV,type="b") # What do you see?

# Testing correlation 
check.corr <- spat.corr.diagnostic(log(Hb) ~ sex, coords = ~web_x+web_y,data=chikwawa,
                     ID=chikwawa$ID,likelihood = "Gaussian",plot.results = FALSE)
plot(check.corr,xlab="Spatial distance",ylab="Variogram") # What do you conclude?


