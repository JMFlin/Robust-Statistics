library(quantreg)
library(robustbase)

#1. Consider the data set Cars available in R
#(a) Compute the mean stopping distance for each initial speed.
#(b) Make a scatter plot for the data and add the mean stopping distances computed above as a line.
#(c) Add to the plot the lines when you fit only a linear model and a quadratic model using least squares to the figure.
#(d) Decide what is the best LS fit.
#(e) Similar choose the best LAD fit for the data.
#(f ) Make a scatter plot and add the regression lines for = 0.05, 0.25, 0.5, 0.75, 0.95 to the data where the model is based on the best LAD model. Add a legend to the figure to explain the curves.

#2. Make yourself familiar with the R data set stackloss.
#(a) Fit an appropriate least squares model using stack.loss as the response variable.
#(b) Are there influential observations? If yes, in which way are they influential?
#(c) Fit the corresponding LAD model for the same data. How does it differ from the least squares fit?
#(d) Estimate robustly the scale based on the LAD fit.
#(e) Using the function Lmrob in the package robustbase to fit an M-estimate as regression for the same model. 
#Use as initial coefficients those from the LAD fit and as initial scale estimate the one estimated above. 
#You need to specify therefore the argument method and init in the lmrob call. Compare it to the other two fits.



#-------------1

my_data <- cars
summary(my_data)#Klaus looked at this and said that if you want to standardize or center you can only subtract
#the minimum so 4. If it was 1000 then you should center or standardize. Look at the other data analysis examples
str(my_data)
head(my_data)

#Plot the mean, lm and poly for the data
my_mean <- aggregate(dist ~ speed, my_data, function(x) mean(as.numeric(x)))

par(mfrow = c(1,1))
plot(my_data)
lines(my_mean)

my_model <- lm(dist~speed, my_data)
summary(my_model)
abline(my_model, col = "blue")

my_model_q <- lm(dist~speed+I(speed^2), my_data)
summary(my_model_q)

plot2 <- function(x) my_model_q$coefficient[3]*x^2 + my_model_q$coefficient[2]*x + my_model_q$coefficient[1] 
curve(plot2, lwd=1, add=T, col="red")

anova(my_model, my_model_q)#best fit is model my_model_q
AIC(my_model)
AIC(my_model_q)#better

#You have to perform quantile regression for the 0.5 quantile (median) 
#to have least absolute deviations instead of least squares 

rq.1 <- rq(my_data$dist~my_data$speed, 0.5)
rq.2 <- rq(my_data$dist~my_data$speed+I(my_data$speed^2), 0.5)

anova(rq.1, rq.2)#best fit is rq.2, but it is not significant
AIC(rq.1)
AIC(rq.2)#better

#Do quantile regresion
rq.a <- rq(my_data$dist~my_data$speed+I(my_data$speed^2), 0.05)
rq.b <- rq(my_data$dist~my_data$speed+I(my_data$speed^2), 0.25)
rq.c <- rq(my_data$dist~my_data$speed+I(my_data$speed^2), 0.5)
rq.d <- rq(my_data$dist~my_data$speed+I(my_data$speed^2), 0.75)
rq.e <- rq(my_data$dist~my_data$speed+I(my_data$speed^2), 0.95)

plot(my_data)

plota <- function(x) rq.a$coefficient[3]*x^2 + rq.a$coefficient[2]*x + rq.a$coefficient[1] 
curve(plota, lwd=1, add=T, col=1)

plotb <- function(x) rq.b$coefficient[3]*x^2 + rq.b$coefficient[2]*x + rq.b$coefficient[1] 
curve(plotb, lwd=1, add=T, col=2)

plotc <- function(x) rq.c$coefficient[3]*x^2 + rq.c$coefficient[2]*x + rq.c$coefficient[1] 
curve(plotc, lwd=1, add=T, col=3)

plotd <- function(x) rq.d$coefficient[3]*x^2 + rq.d$coefficient[2]*x + rq.d$coefficient[1] 
curve(plotd, lwd=1, add=T, col=4)

plote <- function(x) rq.e$coefficient[3]*x^2 + rq.e$coefficient[2]*x + rq.e$coefficient[1] 
curve(plote, lwd=1, add=T, col=5)

names <- c("alpha = 0.05","alpha = 0.25","alpha = 0.5","alpha = 0.75","alpha = 0.95")
legend("topleft", names, col = c(1,2,3,4,5), lty = c(1,1))

#-------------2

my_data <- stackloss
summary(my_data)
str(my_data)
head(my_data)
plot(my_data)

boxplot(my_data)#some outliers?

model <- lm(stack.loss~Air.Flow + Water.Temp + Acid.Conc., my_data)
summary(model)
model.matrix(model)

#Are there influential observations? If yes, which ones?
par(mfrow = c(2,2))
plot(model)#Some outliers in cook's distance and QQ-plot. Observation 21 and maybe observation 4 since they are both
#marked on all 4 plots.

IM.fit4 <- influence.measures(model)
round(IM.fit4$infmat,3)#no influence points??

#Fit the corresponding LAD estimate
rq.loss <- rq(stack.loss~Air.Flow + Water.Temp + Acid.Conc., data = my_data, 0.5)
summary(rq.loss)#Warning is not a problem. The answer can be non-unique.

AIC(model)
AIC(rq.loss)#better

#Estimate the robust scale based on LAD

#As the residuals of the LAD estimate have median zero and many zero
#residuals it is custom to estimate the scale as
(1/0.675)*median(abs(coef(rq.loss)))#All coef are non-zero

#fit an M-estimate as regression for the same modeluse the inital value and scale
#from above.
LAD_scale <- (1/0.675)*median(abs(coef(rq.loss)))

model_robust <- lmrob(stack.loss~Air.Flow+Water.Temp+Acid.Conc., data = my_data, method = "M", 
      init = list(coefficients = rq.loss$coef, scale = LAD_scale))# this is a monotone regression
#it defaults to tukeys estimate. Next practicals we will find better initial values so look there
summary(model_robust)#4 observations c(1,3,4,21) are outliers
#Outlier residuals are downweighted!
#This is iteratively reweighted least squares (IRWLS)

#ANOVA saman mallin eri muotoja. Ei eri metoedien välilla!

#compare residuals and intercept. Can't use model fit. If the parameter estimates differ a lot
#Use robust fit if they differ.

#Quantile regression is used for heteroscedastic data! Growth of child. If you start
#small then you usually follow your quantile.

#These are models to test for prediction.

#in RQ you can only compare the different curves with each other! Not to other models
#LAD can be compared with other models though




#------------------TEACHERS

data(cars)
attach(cars)
#a)
#Means of distances by group
mdist <- aggregate(dist, by = list(speed),  FUN = mean)

#b)
plot(speed, dist)
lines(mdist)

#c)
#Linear
fit1 <- lm(dist ~ speed)
summary(fit1)
lines(speed,fit1$fitted, col = "red")

#You could also center variables in such way that predictors have mean value of 0
# Then intercept means means the expected value of the response,
#   when all predictors are at their means
#Scaling could be good, as some variables with very
# high variablility can produce annoyingly small coefficients
# Change the scale: like area unit of a country could be 1000 square km's 
#   instead of km's


#Quadratic
#Also linear term is included in the quadratic model!
fit2 <- lm(dist ~ speed + I(speed^2)) 
summary(fit2)
lines(speed, fit2$fitted, col = "blue")


#d)
#Quadratic seems cooler

#You could test it:
anova(fit1, fit2) # No significant difference, but:
#AIC smaller in quadratic one (Smaller -> better)
AIC(fit1, fit2) 

# Also linear form may lead to negative values of response dist with very
#   low speed -> bad!
#   -> Quadratic form included -> no such problem

#Decide yourself on how you do the comparison and stick with it here.

#e)
library(quantreg)
fit3 <- rq(dist ~ speed)
summary(fit3)  #Solution may be nonunique
lines(speed, fit3$fitted, col = "green")
#"Solutions cannot be given in closed forms and are usually not unique."

fit4 <- rq(dist ~ speed + I(speed^2))
summary(fit4)
lines(speed, fit4$fitted, col = "gray")
#speed^2 seems cooler

AIC(fit3)
AIC(fit4) #quadratic gives lower AIC -> better

#Also with anova function, you get the same result:
anova(fit3, fit4) #With speed^2 significantly better!
#speed^2 confidence limit does not include zero!


#f)
plot(speed, dist)
alpha <- c(0.05, 0.25, 0.5, 0.75, 0.95)
fit3a <- rq(dist ~ speed + I(speed^2), tau = alpha[1])
summary(fit3a)
lines(speed, fit3a$fitted, col = "blue")
fit3b <- rq(dist ~ speed + I(speed^2), tau = alpha[2])
summary(fit3b)
lines(speed, fit3b$fitted, col = "green")
fit3c <- rq(dist ~ speed + I(speed^2), tau = alpha[3])
summary(fit3c) #Solution may be nonunique
lines(speed, fit3c$fitted, col = "yellow")
fit3d <- rq(dist ~ speed + I(speed^2), tau = alpha[4])
summary(fit3d)
lines(speed, fit3d$fitted, col = "orange")
fit3e <- rq(dist ~ speed + I(speed^2), tau = alpha[5])
summary(fit3e)
lines(speed, fit3e$fitted, col = "red")
legend(5, 120, alpha, lty=c(1,1), lwd=c(2.5,2.5),
       col = c("blue", "green", "yellow", "orange", "red"), title = "alpha")
title("Quantile regression lines")

detach(cars)


#2
#a)
data("stackloss")
summary(stackloss)
head(stackloss)
plot(stackloss) #Scatterplot of the variables

#Some initial thoughts about data:
# - Airflow & water temp: some connection?
# - stackloss: vs. acid & water: ^2; vs. air.flow: linear?

res0a <- lm(stack.loss ~ (Air.Flow + Water.Temp + Acid.Conc.)^2 +
                I(Air.Flow^2) + I(Water.Temp^2) + I(Acid.Conc.^2),
            data = stackloss)
summary(res0a) #Starting here from big model

#Choose here some model choice criterion

#I choose here automatic selection based on step function
res1a <- step(res0a) #Model selection based on AIC (backward selection)
summary(res1a)

#b)
influence.measures(res1a)

#NOTE: Your results may vary according to how you have modelled

plot(res1a)
#Observations 1, 12, 17, 19 and 21 are influential
#12 and 21: very high leverage:
# High-leverage points are those observations, if any, made at extreme or
#   outlying values of the independent variables such that the lack of
#   neighboring observations means that the fitted regression model will
#   pass close to that particular observation. (Wikipedia: Leverage (statistics))
# If a leverage h_i is large, it must be due to extreme values in x_i.

#12 and 21: very high covariance ratio:
#1, 17 and 19: quite high covariance ratio:
# The covariance ratio measures the change in the determinant of the
#   covariance matrix of the estimates by deleting the ith observation.
# - From SAS documentation

#DFBETA: Values larger than 1 should be investigated
#   Assume we are interested in parameter beta_k and omit observation i.
#   Again, 21 has three of these values

#DFFIT: impact on fitted individual values: values larger than 2/sqrt(p/n) i.e.
#   here 3.24:  No obs above this

#Cook's distance which measures the impact of the ith observation on _ALL_
# fitted values.
# - No observations larger than 1

#More detailed version of influence measures:
# - shows more closely which measures are worth to look at
infl1 <- influence.measures(res1a)
summary(infl1) 

#This following produces the same table than before, but
#   it does not show the column where there are the *'s that mark
#   influential observations, so be careful!
infl1$infmat


#c)
#LAD model
#Automatic selection:
#res1b <- rq(stack.loss ~ (Air.Flow + Water.Temp + Acid.Conc.)^2 +
#              I(Air.Flow^2) + I(Water.Temp^2) + I(Acid.Conc.^2),
#            data = stackloss)
#res2b <- step(res1b) #Model selection based on AIC (backward selection)
#summary(res2b)

#..or I use just the same parameters than in 2a)
#   (little difference compared to automatic)
library(quantreg)
res1b <- rq(stack.loss ~ Air.Flow + Water.Temp + Acid.Conc. + Air.Flow:Water.Temp +
                Acid.Conc.:Air.Flow + I(Air.Flow^2) + I(Acid.Conc.^2),
            data = stackloss)

summary(res1b)



#d)
#As r gives zero residuals sometimes in the form like -1.4E-15,
# we need to do something in order to make them be actually zero
resid1 <- round(res1b$residuals, 7) #Residuals with 7 decimals

resid2 <- resid1[resid1 != 0] #Which of those are != 0
#Scale estimate:
sigma <- (1/0.675)*median(abs(resid2))

#e)
library(robustbase)

# - Initial coefficients from LAD model (so keep same parameters than in LAD)
# - Initial scale estimaste = sigma
res1c <- lmrob(stack.loss ~ Air.Flow + Water.Temp + Acid.Conc. +
                   Air.Flow:Water.Temp + Acid.Conc.:Air.Flow + I(Air.Flow^2) +
                   I(Acid.Conc.^2), method = "M",
               init = list(coefficients = res1b$coef, scale = sigma),
               data = stackloss, setting = "KS2014")
#Setting KS2014 is not necessary here, but in categorical data it is necessary
summary(res1c)


#Comparing coefficients of parameter estimates:
summary(res1a)
summary(res1b)
summary(res1c)

#In THIS chosen model:
# - Coefficients of least squares and M-estimates are quite close together
# - In LAD some parameter estimates are much different


#As LAD does not require a scale estimate it is the natural "first" estimate
#to obtain the residuals.
# - It may result to non-unique results.
# - LAD is more like a provider of initial estimates for M-regression.


#Checking some residuals..
plot(NA, xlim = c(0, 22), ylim = c(-4, 8), xaxs = "i")
abline(1,0)
points(res1a$residuals, col = "gray")
points(res1b$residuals, col = "red")
points(res1c$residuals, col = 'blue')

#Note that results depend on the model you ended up with!
