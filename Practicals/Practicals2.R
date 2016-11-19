#Robust statistics  (Fall 2015)
#Practical Session 2

#Question 1
#Compare normal dist to cauchy
#Assume Cauchy distribution with location = 0 and scale = 1

m <- 2000
set.seed(1)
x1 <- replicate(m, mean(rcauchy(100)))
set.seed(1)
x2 <- replicate(m, mean(rcauchy(1000)))
set.seed(1)
x3 <- replicate(m, mean(rcauchy(10000)))

#Just for fun:
mean(x1); mean(x2); mean(x3)
min(x1); min(x2); min(x3)
max(x1); max(x2); max(x3)

#Mean squared errors again with respect to location mu = 0
#MSE = mean((mu_i - mu)^2)
mean((x1 - 0)^2)
mean((x2 - 0)^2)
mean((x3 - 0)^2)
#Basically any random value was able to be there..
# - Depends on in which cases there are really high or low observations

#Just for fun, let's check if MSE were calculated
#   with trimmed mean or with median (not in relaity; not needed):
mean((x1 - 0)^2, trim = 0.1)
mean((x2 - 0)^2, trim = 0.1)
mean((x3 - 0)^2, trim = 0.1)
median((x1 - 0)^2, 0.1)
median((x2 - 0)^2, 0.1)
median((x3 - 0)^2, 0.1)
#Much lower, as the highest values to power 2 are no longer present.


#If X, Y ~ N(0,1), then X/Y ~ Cauchy(0,1).
# - Imagine from X we get every once in a while 3,
#     and from Y in such case sometimes 0.001. Then X/Y would be 3000...
#In few cases one value of mean is so high or low
#   that it might take the whole mean of the means further away from 0.
#  - Though if there are such high values in both sides,
#     it still might look good enough.
# - But: In calculating MSE's we calculate all mean^2's,
#     so all high absolute values ^2 means bad stuff

#QQplot compared to normal distribution
quantnorm <- qnorm(1:m/(m + 1), 0, 1)
qqplot(quantnorm, x3)
qqline(x3)
#Or
qqnorm(x3)
qqline(x3)
#Tails are clearly not in line

#QQplot compared to cauchy distribution
quantcauchy <- qcauchy(1:m/(m + 1), 0, 1)
qqplot(quantcauchy, x3)
qqline(x3)
#Even own quantiles are not exactly in line,
#   as some values can be really big compared to quantiles

#Question 2
#Plot the influence function of the winzorized mean
alpha <- c(0.15, 0.3, 0.45)
mu <- 0; sd <- 1
qal <- qnorm(alpha, mu, sd)
q1al <- qnorm(1 - alpha, mu, sd)
dqal <- dnorm(qal, mu, sd)
dq1al <- dnorm(q1al, mu, sd)
twmf <- 0 + alpha*qal + alpha*q1al #T_WM(F) = 0 (symmetric around zero!)
cf <- twmf - alpha^2/dqal + alpha^2/dq1al # C(F) = 0 (symmetric around zero!)
x1 <- qal - alpha/dqal - cf #For x < qnorm(alpha)
x2 <- q1al + alpha/dq1al - cf #For x > qnorm(1 - alpha)

for (i in 1:3){
    #Winsorized mean
    #Create an empty plot:
    plot(NA, xlim = c(-10, 10), ylim = c(-10, 10), xaxs = "i",
         main = paste("alpha = ", alpha[i]))
    #For qnorm(alpha) < x < qnorm(1 - alpha):
    curve((x - cf[i]), from = qal[i], to = q1al[i], add = T, lwd = 2) 
    #For x < qnorm(alpha):
    curve(0*x + x1[i], from = -10, to = qal[i], add = T, lwd = 2)
    #For x > qnorm(1 - alpha):
    curve(0*x + x2[i], from = q1al[i], to = 10, add = T, lwd = 2)
    #One plot at a time; press enter while in Console window to get next one:
    readline() 
    #Trimmed mean (just for comparison)
    curve((1/(1 - 2*alpha[i]))*(x - twmf[i]),
          from = qal[i], to = q1al[i], add = T, col = "blue", lwd = 2)
    curve(0*x + (1/(1 - 2*alpha[i]))*(qal[i] - twmf[i]),
          from = -10, to = qal[i], add = T, col = "blue", lwd = 2)
    curve(0*x + (1/(1 - 2*alpha[i]))*(q1al[i] - twmf[i]),
          from = q1al[i], to = 10, add = T, col = "blue", lwd = 2)
    readline()
    #Mean (just for comparison)
    curve(1*x, from = -10, to = 10, add = T, col = "red", lwd = 2)
    readline()
}

#When observation x is winsorized, it same influence than all the others winsorizad!
# Within nonwinsorized area, if observation grows, so does the winsorized mean
#   - Same are true for the trimmed mean

#This is still valid even if we change means or variances
# - After some point influence anyway stays the same.
# - Also normal distribution is symmetric; it will just change the
#     location of the function, or rescale it a bit, but no essential change.
#   - Think about location-scale model!


#Sensitivity curves
set.seed(1)
scnorm <- rnorm(50, 0, 1)
x <- seq(-10, 10, 0.1)
#Making new datasets with one value added
scnorm1 <- data.frame(matrix(NA, nrow = 51, ncol = length(x)))
for(i in 1:length(x)) {
    scnorm1[i] <- c(x[i], scnorm)
}
#Curves:
#For winsorized mean
for (i in 1:3){
    plot1 <- apply(scnorm1, 2, winmean, alpha[i]) - winmean(scnorm, alpha[i])
    readline()
    plot(x, plot1, type = 'l')
}
#For trimmed mean (just for comparison..)
for (i in 1:3){
    plot1 <- apply(scnorm1, 2, trimmean, alpha[i]) - trimmean(scnorm, alpha[i])
    readline()
    plot(x, plot1, type = 'l')
}


#Question 3: Huber's M-estimate of location
#Write your own HUber's function
# w=u(r) is psi(r)/(r)
huberM <- function(x, c, eps = 1e-06, init = 0)
{
    mu <- init
    diff <- Inf
    while (diff > eps)
    {
        r <- x - mu
        psi.r <- r #In general value is just r
        psi.r <- ifelse(r <= -c, -c, psi.r) #If < -c then -c, otherwise psi.r (= r)
        psi.r <- ifelse(r >=  c,  c, psi.r) #If > c then c, otherwise psi.r (= r)
        #r can be zero in case of median -> eps added (Not the most elegant way):
        w <- psi.r/(r + eps)
        mu.new  <- sum(w * x)/sum(w) #New value for mu
        #Check if it differs much from previous estimate:
        diff <- abs(mu - mu.new) #If diff < eps, then convergence is achieved
        mu <- mu.new
    }
    mu.new #Return mu value
} 

x <- rnorm(101,4)
huberM(x, 1.5)
huberM(x, 1.5, init = median(x)) #Starting with robust location estimate

huberM(x, c = 1.345, eps = 1E-06)

#huber-function in MASS package in R gives different results,
#   as it takes also scale into account.



#Question 4
#Asymptotic relative efficiency (ARE)
m <- 1000
set.seed(1)
rep1 <- replicate(m, rt(1000, 3))
tdm1 <- colMeans(rep1)
tdm2 <- apply(rep1, 2, median)
var(tdm1)/var(tdm2) 
#median vs. mean -> Median more efficient asymptotically

tdm3 <- apply(rep1, 2, huberM, c = 1.345, eps = 1E-06)
var(tdm1)/var(tdm3) 
#Huber vs. mean -> Huber more efficient asymptotically

#Extra check:
var(tdm2)/var(tdm3)
#Huber vs. median -> Huber more efficient asymptotically

#Efficiency of Huber M estimate
#ARE(huber, mean) = AVAR(mean)/AVAR(huber)

#Median is more efficient than mean
#Huber M estimate is more efficient than mean (and median)