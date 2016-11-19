library(quantreg)
library(robustbase)

#1. Reconsider again the R data set stackloss.
#(a) Fit an MM regression to the stackloss data. Does it differ from the fits from the previous practicals?
#(b) Are there outliers? Use the weights to decide how many and which?
#(c) Look at the diagnostic plots of the MM- t. Are you satisfied?

#2. Make yourself familiar with the coleman data in the robustbase package.
#(a) Plot the scatter plot matrix of the data. Are there any suspicious observations?
#(b) Fit a LS model to the data where all variables enter the model linearly. Perform a diagnostic check of the model. Are you satisfied?
#(c) Fit the same model as a MM regression to the data using the default bisquare but with 85% efficiency 
#(the tuning constant should be then 3.44, or more precisely 3.443689). 
#Check lmrob.control for help. Check the fit.
#(d) Do you think now that there are outliers in the data? 
#If yes, identify them. Plot the scatter matrix again and mark them
#- if there are outliers - using color. 
#If there are outliers - do they stick anywhere esp ecially out?



#-------------------------1
my_data <- stackloss

#Fit an MM regression to the stackloss data. Does it differ from
#the fts from the previous practicals?

r_model <- lmrob(stack.loss~Air.Flow+Water.Temp+Acid.Conc., data = my_data, method = "MM")

#Are there outliers? Use the weights to decide how many and which?

summary(r_model)#tunin.psi is 95% efficiency here.
#The higher the efficiency is chosen the more suffers the robustness and bias might be introduced.

#Adjusted R-squared:  0.9521

#observation 21 is an outlier

#Look at the diagnostic plots of the MM-fit. Are you satisfied?

par(mfrow = c(2,3))
plot(r_model)
#4 and 21 seem to be influential here. Maybe observation 3

#----------------2

#Plot the scatter plot matrix of the data. Are there any suspicious observations?

my_data1 <- coleman
pairs(my_data1)

#Fit a LS model to the data where all variables enter the model
#linearly. Perform a diagnostic check of the model. Are you satisfied?

fit1 <- lm(Y ~., data = my_data1)
summary(fit1)
par(mfrow = c(2,2))
plot(fit1)#18 is influential, maybe 3
influence.measures(fit1)

#Fit the same model as a MM regression to the data using the
#default bisquare but with 85% effciency 
#(the tuning constant should be then 3.44, or more precisely 3.443689). 
#Check lmrob.control for help. Check the fit.
#The higher the efficiency is chosen the more suffers the robustness
#and bias might be introduced. 85% efficiency are usually considered reasonable.

r_model1 <- lmrob(Y ~., data = my_data1, method = "MM", tuning.psi = 3.44)
summary(r_model1)
par(mfrow = c(2,3))
plot(r_model1)#18 and 3 are influential

#Do you think now that there are outliers in the data? If yes,
#identify them. Plot the scatter matrix again and mark them - if
#there are outliers - using color. If there are outliers - do they stick
#anywhere especially out?

summary(r_model1)#2 observations c(3,18) are outliers. These are the row numbers in your data set!
pairs(my_data1, col=ifelse(rownames(my_data1)==3, "red", ifelse(rownames(my_data1)==18, "red", "black")))
#They don't look like outliers, which is exactly the point!

#At work compare method = "MM" to setting = "KS2014" for model comparison.
test <- lmrob(Y ~., data = my_data1, setting = "KS2014")
summary(test)

#influence measure = outlier? Not always. 

#What diagnostics to believe? Summary outlier detection.



#----------------TEACHER


#Question 1

data("stackloss")

# a)
#Using the same model than last time!
library(robustbase)
res1d <- lmrob(stack.loss ~ Air.Flow + Water.Temp + Acid.Conc. + Air.Flow:Water.Temp +
                   Acid.Conc.:Air.Flow + I(Air.Flow^2) + I(Acid.Conc.^2), method = "MM",
               data = stackloss, setting = "KS2014")

res1a$coefficients #Previous exercise
res1b$coefficients #Previous exercise
res1c$coefficients #Previous exercise
res1d$coefficients

# b)
round(weights(res1d, type = "r"), 3) #Any < 0.2?

# c)
par(mfrow = c(3,2))
plot(res1d) #Hit <Return> to see next plot
#Obs 4 might have some issues
# Low weight (< 0.3)
# -> Quite high residual, sqrt(|residual|), not too good in QQplot,
# -> high robust std. residual



#Question 2
data(coleman)
# a)
# Checking the data
summary(coleman)
head(coleman) #See top of the dataset
plot(coleman)
#Some observations a bit further from others


# b)
#Least squares fit
# . includes all variables in to the model as main effects (except response; Y here)
lmfit1 <- lm(Y ~ ., data = coleman) 

plot(lmfit1) #Hit <Return> to see next plot
#Obs 18 seems to be influential
#Obs 3 likely as well
#Maybe 11 too; 12?

#Method MM is the default; no need to add it there
lmrobfit1 <- lmrob(Y ~ ., data = coleman, tuning.psi = 3.443689)
plot(lmrobfit1) #Hit <Return> to see next plot
summary(lmrobfit1) #tuning.psi looks as it should be according to our settings
#2 observations c(3,18) are outliers with |weight| = 0 ( < 0.005)
w1 <- round(weights(lmrobfit1, type = 'r'), 3)
w1

# d)
#Observations 3 and 18 seem to be outliers: weight zero!
#Others: at least 0.8

#Plot with different colours for observations with low weight
plot(coleman, col = ifelse(w1 >= 0.2, 1, 2), pch = 15)

#They do not seem stand out that clearly from the plot