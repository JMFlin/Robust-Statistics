library(MASS)

#This data set contains the effect of different forms of therapy on the body
#weight of girls suffering from anorexia.

data(anorexia)
str(anorexia)
summary(anorexia)


anorexia$TREAT <- relevel(anorexia$Treat, ref="Cont")
par(mfrow = c(1,2))
boxplot(Prewt~TREAT, data=anorexia, ylab="preweight")
boxplot(Postwt~TREAT, data=anorexia, ylab="postweight")

Prewt.mean <- with(anorexia, tapply(Prewt, TREAT, mean))
Prewt.sd <- with(anorexia, tapply(Prewt, TREAT, sd))
Postwt.mean <- with(anorexia, tapply(Postwt, TREAT, mean))
Postwt.sd <- with(anorexia, tapply(Postwt, TREAT, sd))
n.group <- with(anorexia, tapply(Prewt, TREAT, length))
ANOREX.SUMMARY <- cbind(n.group, Prewt.mean, Prewt.sd, Postwt.mean, Postwt.sd)
ANOREX.SUMMARY

anfit1 <- lm(Postwt~TREAT, data=anorexia)
sixth <- seq(from = 6, to = 72, by = 6)
model.matrix(anfit1)[sixth,]#If one is interested how the design matrix looks one can use the function model.matrix.
#This function returns for an lm object the design matrix where one for
#example can see which contrast was used for a factor and so on. Especially
#when there are factors in your model it might be a good idea to check this
#matrix so that you know how to interpret the result.

anfit1b <- lm(Postwt~TREAT-1, data=anorexia)
model.matrix(anfit1b)[sixth,]

summary(anfit1b)

anfit1c <- lm(Postwt~TREAT, data=anorexia, contrast=list(TREAT="contr.sum"))
model.matrix(anfit1c)[sixth,]

summary(anfit1c)

anfit1d <- lm(Postwt~TREAT, data=anorexia, contrast=list(TREAT="contr.helmert"))
model.matrix(anfit1d)[sixth,]

summary(anfit1d)

par(mfrow=c(2,2))
plot(Postwt~Prewt, data=anorexia, col=as.numeric(TREAT), pch=as.numeric(TREAT),
        main="All groups")
plot(Postwt~Prewt, data=anorexia, col=as.numeric(TREAT), pch=as.numeric(TREAT),
        main="Cont", subset=TREAT=="Cont")
plot(Postwt~Prewt, data=anorexia, col=as.numeric(TREAT), pch=as.numeric(TREAT),
        main="CBT", subset=TREAT=="CBT")
plot(Postwt~Prewt, data=anorexia, col=as.numeric(TREAT), pch=as.numeric(TREAT),
        main="FT", subset=TREAT=="FT")

anfit2 <-lm(Postwt~TREAT + Prewt, data=anorexia)
coef(anfit2)

par(mfrow = c(1,1))
plot(Postwt~Prewt, data=anorexia, col=as.numeric(TREAT), pch=as.numeric(TREAT))
abline(coef(anfit2)[1],coef(anfit2)[4], col=1)
abline(coef(anfit2)[1]+coef(anfit2)[2],coef(anfit2)[4], col=2)
abline(coef(anfit2)[1]+coef(anfit2)[3],coef(anfit2)[4], col=3)
with(anorexia,points(Prewt,fitted(anfit2),pch=15, col=as.numeric(TREAT)))


anfit3 <- lm(Postwt~TREAT*Prewt, data=anorexia)
summary(anfit3)
model.matrix(anfit3)[sixth,]
coef(anfit3)

plot(Postwt~Prewt, data=anorexia, col=as.numeric(TREAT), pch=as.numeric(TREAT))
abline(coef(anfit3)[1],coef(anfit3)[4], col=1)
abline(coef(anfit3)[1]+coef(anfit3)[2],coef(anfit3)[4]+coef(anfit3)[5], col=2)
abline(coef(anfit3)[1]+coef(anfit3)[3],coef(anfit3)[4]+coef(anfit3)[6], col=3)
with(anorexia,points(Prewt,fitted(anfit3),pch=15, col=as.numeric(TREAT)))

min(anorexia$Prewt)

anorexia$Prewt2 <- anorexia$Prewt - min(anorexia$Prewt)
anfit4 <- lm(Postwt~TREAT*Prewt2, data=anorexia)
summary(anfit4)
coef(anfit4)

par(mfrow=c(1,1))
plot(Postwt~Prewt2, data=anorexia, col=as.numeric(TREAT), pch=as.numeric(TREAT))
abline(coef(anfit4)[1],coef(anfit4)[4], col=1)
abline(coef(anfit4)[1]+coef(anfit4)[2],coef(anfit4)[4]+coef(anfit4)[5], col=2)
abline(coef(anfit4)[1]+coef(anfit4)[3],coef(anfit4)[4]+coef(anfit4)[6], col=3)
with(anorexia,points(Prewt2,fitted(anfit4),pch=15, col=as.numeric(TREAT)))

par(mfrow=c(2,2))
plot(anfit4)

IM.fit4 <- influence.measures(anfit4)
INFobs <- apply(IM.fit4$is.inf, 1, any)
round(IM.fit4$infmat[INFobs,],3)

