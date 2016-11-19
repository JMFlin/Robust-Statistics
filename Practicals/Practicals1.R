#Robust statistics  (Fall 2015)
#Practical Session 1

#Question 1: Function for contaminated normal distribution
#
#G is the "good" one, H causes the contamination

#1. a random variable is selected by chance from the collection
#   according to given probabilities of selection
#2. The value of the selected random variable is drawn
rcnorm <- function(n, mu1, sd1, mu2, sd2, eps) {
    mu <- c(mu1, mu2)
    sd <- c(sd1, sd2)
    #Which normal model to take:
    choice <- sample(1:2, prob = c(1 - eps, eps), size = n, replace = T) 
    rnorm(n, mu[choice], sd[choice]) 
}
#To get a result x from function, either DO NOT assign it to anything, i.e.
# just use        x  (or return(x))
# and __NOT__     something <- x

set.seed(1); ans1 <- rcnorm(9997, 0, 1, 5, 1, 0.1) #Testing
plot(density(ans1))

#Wrong way to do it:
plot(density(0.9*rnorm(100,0,1) + 0.1*rnorm(100, 5, 1)))
#This draws values from both distributions and just gives them weights!


#Question 2: Function for trimmed mean
trimmean <- function(x, alpha) {
    n <- length(x)
    orx <- x[order(x)]
    m <- floor(alpha * n)
    mean(orx[(m + 1):(n - m)])
}
trimmean(ans1, 0.1) #Testing
mean(ans1, trim = 0.1) #Compare to this (should be exatly this)


#Question 3: Function for winsorized mean
winmean <- function(x, alpha) {
    n <- length(x)
    orx <- x[order(x)]
    m <- floor(alpha * n)
    (1/n) * (sum(orx[(m + 1):(n - m)])
             + sum(rep(orx[(m + 1)], m)) + sum(rep(orx[(n - m)], m)))
}
winmean(ans1, 0.1) #Testing

library(psych)
winsor.mean(ans1, trim = 0.1) #Compare to this (should be very close to this)
#NOTE: In r winsorized mean is defined differently,
#   i.e. it will give a bit different results


#Question 4

# Write       a       function sim.rc
#which       simulates       for       given       values       of
#n,mu1,sd1,mu2,sd2,eps and alpha from      the      contaminated      normal      mo del and     
#computes      for      the       sample      the       mean,       median,       trimmed      mean       and
#winsorized       mean       with       b oth alpha=0.1.           The       function       should       rep eat
#this m times     and     return     a     dataframe     with     the     values     of     the     parameters
#submitted        and        the        four        lo cation        estimates.       

#With for loop (not that efficient; can be slow if samples are large):
sim.rc <- function(m, n, mu1, sd1, mu2, sd2, eps, alpha = 0.1){
    dfr <- data.frame(matrix(NA, nrow = m, ncol = 4))
    colnames(dfr) <- c("Mean", "Median", "Trim. mean", "Wins. mean")
    for (i in 1:m) {
        cn <- rcnorm(n, mu1, sd1, mu2, sd2, eps)
        dfr[i, ] <- data.frame(mean(cn), median(cn),
                               trimmean(cn, alpha), winmean(cn, alpha))
    }
    dfr
}
#Or with the replicate function (more efficient in R):
sim.rc2 <- function(m, n, mu1, sd1, mu2, sd2, eps, alpha = 0.1){
    rep1 <- replicate(m, rcnorm(n, mu1, sd1, mu2, sd2, eps))
    dfr <- data.frame(cbind(colMeans(rep1), apply(rep1, 2, median),
                            apply(rep1, 2, trimmean, alpha), apply(rep1, 2, winmean, alpha)))
    colnames(dfr) <- c("Mean", "Median", "Trim. mean", "Wins. mean")
    dfr
}

set.seed(24); sim.rc(10, 997, 0, 1, 5, 1, 0.1)
set.seed(24); sim.rc2(10, 997, 0, 1, 5, 1, 0.1)
#Both give the exact same result!

#Question 5
set.seed(1); s11 <- sim.rc2(200, 100, 0, 1, 1.5, 1, 0)
set.seed(1); s21 <- sim.rc2(200, 200, 0, 1, 1.5, 1, 0)
set.seed(1); s31 <- sim.rc2(200, 500, 0, 1, 1.5, 1, 0)
set.seed(1); s12 <- sim.rc2(200, 100, 0, 1, 1.5, 1, 0.05)
set.seed(1); s22 <- sim.rc2(200, 200, 0, 1, 1.5, 1, 0.05)
set.seed(1); s32 <- sim.rc2(200, 500, 0, 1, 1.5, 1, 0.05)
set.seed(1); s13 <- sim.rc2(200, 100, 0, 1, 1.5, 1, 0.15)
set.seed(1); s23 <- sim.rc2(200, 200, 0, 1, 1.5, 1, 0.15)
set.seed(1); s33 <- sim.rc2(200, 500, 0, 1, 1.5, 1, 0.15)

par(mfrow=c(1,3))
#
boxplot(s11); title("n = 100, eps = 0")
boxplot(s21); title("n = 200, eps = 0")
boxplot(s31); title("n = 500, eps = 0")
#
boxplot(s12); title("n = 100, eps = 0.05")
boxplot(s22); title("n = 200, eps = 0.05")
boxplot(s32); title("n = 500, eps = 0.05")
#
boxplot(s13); title("n = 100, eps = 0.15")
boxplot(s23); title("n = 200, eps = 0.15")
boxplot(s33); title("n = 500, eps = 0.15")
par(mfrow=c(1,1))
#Use zoom to see all the boxplot clearer
#This could be done more automatically as well, as some of you did


#Mean squared errors with respect to location of good data
#   - In good data location mu = 0
#   - MSE = mean((mu_i - mu)^2)

#For mean in case of n = 100, eps = 0: Diferent ways to calculate
1/length(s11[ , 1]) * sum((s11[ , 1] - 0)^2)
mean((s11[ , 1] - 0)^2)
library(Metrics)
mse(rep(0, length(s11[ , 1])), s11[ , 1]) #Function in R


#For all methods in all settings:
apply((s11 - 0)^2, 2, mean) #In apply function: 1 = row, 2 = column
apply((s21 - 0)^2, 2, mean)
apply((s31 - 0)^2, 2, mean)
#When eps = 0, then no big differences between methods
#   - All symmetric around zero
#   - Median has the largest MSE,as its values change more,
#       compared to other methods (see boxplots)

apply((s12 - 0)^2, 2, mean)
apply((s22 - 0)^2, 2, mean)
apply((s32 - 0)^2, 2, mean)
#When eps = 0.05, then mean starts to deviate a bit more from 0
#   and has larger MSE value compared to median when n is large
#Median values still the most spread (boxplot), but contamination starts
#   to affect mean more to make it worse than median

apply((s13 - 0)^2, 2, mean)
apply((s23 - 0)^2, 2, mean)
apply((s33 - 0)^2, 2, mean)
#When eps = 0.15, median is the best of all,
#   as trimmed mean and winsorized mean starts to get affected as well
#More contamination also means higher MSE's, as the location estimates
#   are further away from zero an MSE basically uses here their squared values


#When n grows, MSE's get smaller as samples get larger (biggest change in eps = 0):
#   location estimates get closer together
#   contamination reduces this effect: MSE's do not get smaller that fast
# For example if eps = 0 we just follow normal distribution
#   With large n, location estimates in each repetition are much closer
#     to true value 0