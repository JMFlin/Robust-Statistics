library(MASS)

#As a first example consider the trees data set from the MASS package.
#The data set contains the girth, height and volume of 31 felled black
#cherry trees. The aim is to obtain a model which can be used to predict
#the volume of a tree based on its height and girth.

data(trees)
head(trees)
summary(trees)
plot(trees)

#Let us first fit a marginal model for the two explaining variables.

fit.girth <- lm(Volume ~ Girth, data = trees)
summary(fit.girth)

fit.height <- lm(Volume ~ Height, data = trees)
summary(fit.height)

#Model containing both explaining variables.

fit.both <- lm(Volume ~ Girth + Height, data = trees)
summary(fit.both)
coef(fit.both)
confint(fit.both)
anova(fit.both)

#A full model might here rather a model with a second degree polynomial
#for both variables.

fit.full <- lm(Volume ~ Girth + I(Girth^2) + Height + I(Height^2), data = trees)
summary(fit.full)

coef(fit.full)
confint(fit.full)
anova(fit.full)

#The polynomials make the parameters difficult to interpret.

with(trees, cor(Girth, Girth^2))
with(trees, cor(Height, Height^2))

m.Girth <- with(trees, mean(Girth))
m.Height <- with(trees, mean(Height))
with(trees, cor(Girth-m.Girth, (Girth-m.Girth)^2))

#So lets us the centered variables.

fit.full.c <- lm(Volume ~ I(Girth-m.Girth)+ I((Girth-m.Girth)^2)+ I(Height-m.Height) + I((Height-m.Height)^2), data = trees)
summary(fit.full.c)#Last one is not stats signif


fit.2 <- lm(Volume ~ I(Girth-m.Girth)+ I((Girth-m.Girth)^2) + I(Height-m.Height), data = trees)
summary(fit.2)#better without last

#Some model comparison if we would need the dropped term

anova(fit.full.c, fit.2)#model fit.full is not significantly different from fit.2
AIC(fit.full.c)
AIC(fit.2)#fit 2 has a smaller AIC

par(mfrow=c(2,2))
plot(fit.2)

influence.measures(fit.2)#leverage values for each observation are called "hat"

#The model in this form might not be fully satisfactory.