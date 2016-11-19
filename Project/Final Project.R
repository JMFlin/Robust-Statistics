#Janne Flinck
#79428

ID <- 79428
set.seed(ID)

Hampel <- function(x, a=1.2, b=3.5, c=8.0, eps = 1e-06){
    mu <- median(x)
    s <- mad(x)
    diff <- Inf
    while(diff > eps*s){
        r <- (x - mu)/s
        psi.r <- r
        psi.r <- ifelse(0 <= abs(x) | abs(x) < a, abs(x)*sign(x), psi.r)
        psi.r <- ifelse(a <= abs(x) | abs(x) <= b, a*sign(x), psi.r) 
        psi.r <- ifelse(b <= abs(x) | abs(x) <= c, ((c-abs(x))/(c-b)) *a*sign(x), psi.r)
        psi.r <- ifelse(abs(x) > c, 0, psi.r) 

        w <- psi.r/r
        w[r == 0] <- 1
        mu.new  <- sum(w * x)/sum(w)
        diff <- abs(mu - mu.new)
        mu <- mu.new
    }
    mu.new #Return mu value
} 
set.seed(ID);x2 <- rnorm(3,4)
Hampel(x2)

#it is a redescender

H.38 <- chgDefaults(hampelPsi, k = c(1.2, 3.5, 8))
H.38
par(mfrow = c(1,1))
plot(H.38)

#--------------------------------------3

alpha <- c(1.2, 3.5, 8)

set.seed(ID)
scnorm <- rnorm(50, 0, 1)
x <- seq(-10, 10, 0.1)
#Making new datasets with one value added
scnorm1 <- data.frame(matrix(NA, nrow = 51, ncol = length(x)))
for(i in 1:length(x)) {
    scnorm1[i] <- c(x[i], scnorm)
}
#Curves:
for (i in 1:1){
    plot1 <- apply(scnorm1, 2, Hampel, alpha[1]) - Hampel(scnorm, alpha[1])
    readline()
    plot(x, plot1, type = 'l')
}

#--------------------------4
DATA <- read.csv("Data.txt", sep = "")

ID <- 79428
set.seed(ID)
DATA.ANA <- DATA[sort(sample(1:nrow(DATA),800)), ]

#Descriptive statistics
library(psych)
str(DATA.ANA)
head(DATA.ANA)
describe(DATA.ANA)[3:12]
#I first briefly examine the data. I notice that the variables have vastly different means.
#Not surprisingly, the variables also have vastly different variances
table(DATA.ANA$DIET)
#Visuals

#PAIRS
pairs(DATA.ANA[,2:6], gap=0, pch=19, cex=0.4, col="darkblue")

#pairs for men and women separately
par(mfrow = c(1,1))
a <- subset(DATA.ANA, SEX == "female")
plot(jitter(a$AGE), a$HOR_LEVEL, col="darkblue", 
     xlab = "Age", ylab = "Hormone level", main = "Female")
b <- subset(DATA.ANA, SEX == "male")
plot(jitter(b$AGE), b$HOR_LEVEL, col="darkblue",
     xlab = "Age", ylab = "Hormone level", main = "Male")

#HISTOGRAMS
library(ggplot2)
library(gridExtra)
my_data <- DATA.ANA[,2:6]
my_data$DIET <- as.numeric(my_data$DIET)
my_data$SEX <- as.numeric(my_data$SEX)

for(i in 1:5){
    breaks <- pretty(range(my_data[,i]), n = nclass.Sturges(my_data[,i]), min.n = 1)
    bwidth <- breaks[2]-breaks[1]
    df <- data.frame(my_data[,i])
    a <- ggplot(df,aes(my_data[,i]))+geom_histogram(binwidth=bwidth,fill="white",colour="black")+
        labs(title = paste(names(my_data)[i],"- Sturges' Method"))+
        theme(plot.title = element_text(size = rel(1.5)))
    
    
    breaks <- pretty(range(my_data[,i]), n = nclass.scott(my_data[,i]), min.n = 1)
    bwidth <- breaks[2]-breaks[1]
    df <- data.frame(my_data[,i])
    b <- ggplot(df,aes(my_data[,i]))+geom_histogram(binwidth=bwidth,fill="white",colour="black")+
        labs(title = paste(names(my_data)[i],"- Scott's Method"))+
        theme(plot.title = element_text(size = rel(1.5)))
    
    
    breaks <- pretty(range(my_data[,i]), n = nclass.FD(my_data[,i]), min.n = 1)
    bwidth <- breaks[2]-breaks[1]
    df <- data.frame(my_data[,i])
    c <- ggplot(df,aes(my_data[,i]))+geom_histogram(binwidth=bwidth,fill="white",colour="black")+
        labs(title = paste(names(my_data)[i],"- Freedman-Diaconis' Method"))+
        theme(plot.title = element_text(size = rel(1.5)))
    
    grid.arrange(a, b, c, ncol=1)
}#Data is not normally distributed

#CORRELATIONS
vector <- c(1:6)
remove <- c(1,3,4)
vector <- vector[! vector %in% remove]#remove non-continuous variables
round(cor(DATA.ANA[,vector], use="complete.obs", method="kendall"), 3)
round(cor(DATA.ANA[,vector], use="complete.obs", method="pearson"), 3)

#PCA
pca1 <- prcomp(my_data[,1:5], scale = TRUE, center = TRUE)

pr.var1 <- pca1$sdev^2
pve1 <- pr.var1/sum(pr.var1)

variances <- data.frame(variances=pca1$sdev**2, pcomp=1:length(pca1$sdev))
varPlot <- ggplot(variances, aes(pcomp, variances))+ geom_bar(stat="identity", fill="gray") + geom_line() + 
    labs(title = "Scree plot with normalization")+
    theme(plot.title = element_text(size = rel(1.5)))
varPlot
par(mfrow = c(1,1))
plot(cumsum(pve1), main = "Normalized", xlab="Principal  Component", 
     ylab="Cumulative  Proportion  of Variance  Explained", ylim=c(0,1),type="b")

library(grid)
PCbiplot <- function(PC, x="PC1", y="PC2", mm) {
    # PC being a prcomp object
    # data <- data.frame(obsnames=row.names(PC$x), PC$x)
    data <- data.frame(obsnames="", PC$x)#data.frame(obsnames=my_data[,12], PC$x)
    plot <- ggplot(data, aes_string(x=x, y=y))# + geom_text(alpha=.4, size=3, aes(label=obsnames))
    #geom_text is redundant here because we have obsnames as blank
    plot <- plot + geom_hline(aes(0), size=.2) + geom_vline(aes(0), size=.2) + geom_point(alpha = 1,size = 2)
    datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
    mult <- min(
        (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
        (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
    )
    datapc <- transform(datapc,
                        v1 = .7 * mult * (get(x)),
                        v2 = .7 * mult * (get(y))
    )
    plot <- plot + coord_equal() + geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 6.5, vjust=1, color="red")
    plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="red")+
        labs(title = mm)+
        theme(plot.title = element_text(size = rel(1.5)))
    plot
}
fit.1 <- prcomp(my_data[,1:5], scale=TRUE, center = TRUE)#Zero mean and unit variance

PCbiplot(fit.1, mm = "PCA with normalization")

#2D MDS
library(ggfortify)
a <- scale(my_data, scale = TRUE, center = TRUE)

d <- dist(a, method = "euclidean") # euclidean distances between the rows
fit <- cmdscale(d, eig=TRUE, k=2) # k is the number of dimensions
plot <- autoplot(fit, xlab="Coordinate 1", ylab="Coordinate 2", main = "Metric MDS")
plot <- plot +theme(plot.title = element_text(size = rel(1.5))) +geom_hline(aes(0), size=.2) + geom_vline(aes(0), size=.2)
plot

#--------------------MODELING

#Basic model
fit <- lm(HOR_LEVEL~.-ID, DATA.ANA)
summary(fit)
confint(fit)

# Look for assumption violation for continuous variables
# y,c, i.i.d?
ylag <- DATA.ANA$HOR_LEVEL[-1]     #shift y_1 by one
ynew <- DATA.ANA$HOR_LEVEL[-(nrow(DATA.ANA))]   #need equally sized vectors...
clag <- DATA.ANA$AGE[-1]     #do the same for c..
cnew <- DATA.ANA$AGE[-nrow(DATA.ANA)]

#Is there a relation between the ith and the i+1th observation of y? (Note that the lag might be larger, i.e. the ith and i+kth observation might be correlated)
cor(ylag,ynew)
cor(cnew,clag)
#does not seem to be the case.

#do residuals and yhats correlate? 
#Is the relationship linear? Plotting y as a function of c...
par(mfrow = c(2,2))
plot(fit)
#first plot: visual inspection -->no relation between residuals and fitted values--> Cov(X'u)=0 OK
#qq plot: some departures from the line -- doubts about a normal distribution.
par(mfrow = c(1,1))
hist(resid(lm(HOR_LEVEL~.-ID, DATA.ANA)), "Scott")#Residuals are not normally distributed
round(cbind(coef = coef(fit), confint(fit)),3)

###we could perform a shapiro-wilk test for normality. H0: normality H1:non-normality
shapiro.test(resid(lm(HOR_LEVEL~.-ID, DATA.ANA)))#H0 is rejected, becase the distribution of our data is not normal


#Basic model with interaction terms
fit.1 <- lm(HOR_LEVEL~(AGE + SEX + DIET + SPORT)^2, DATA.ANA)
summary(fit.1)
par(mfrow = c(2,2))
plot(fit.1)
round(cbind(coef = coef(fit.1), confint(fit.1)),3)

anova(fit, fit.1)

#Checking polynomials for AGE
library(boot)
set.seed(ID)
cv.error.10 <- rep(0,10)#We begin by initializing the vector.
for(i in 1:10){
    glm.fit <- glm(HOR_LEVEL~poly(AGE,i), data = DATA.ANA)
    cv.error.10[i] <- cv.glm(DATA.ANA, glm.fit, K = 10)$delta[1]
}#10-fold cross validation
cv.error.10
par(mfrow = c(1,1))
plot(1:10, cv.error.10, xlab="Degree", ylab="CV error", type="o", lwd=1)
min.point <- min(cv.error.10)
sd.points <- sd(cv.error.10)
abline(h=min.point + 0.2 * sd.points, col="red", lty="dashed")
abline(h=min.point - 0.2 * sd.points, col="red", lty="dashed")
legend("topleft", "0.2-standard deviation lines", lty="dashed", col="red")
which.min(cv.error.10)


#Fitting the poly model
fit.2 <- lm(HOR_LEVEL~.-ID + I(AGE^2) + (AGE + SEX + DIET + SPORT)^2, DATA.ANA)
summary(fit.2)
par(mfrow = c(2,2))
plot(fit.2)
round(cbind(coef = coef(fit.2), confint(fit.2)),3)

#Fitting a model where I have removed the insignificant interaction terms
fit.3 <- lm(HOR_LEVEL~.-ID + I(AGE^2) + AGE:SPORT + DIET:SPORT, DATA.ANA)
summary(fit.3)
par(mfrow = c(2,2))
plot(fit.3)
round(cbind(coef = coef(fit.3), confint(fit.3)),3)

anova(fit, fit.1, fit.2, fit.3)
AIC(fit)
AIC(fit.1)
AIC(fit.2)
AIC(fit.3)#best one
summary(influence.measures(fit.2))

#Model Selection
library(leaps)
#FORWARD AND BACKWARD STEPWISE SELECTION
regfit.fwd <- regsubsets(HOR_LEVEL~.-ID + I(AGE^2) + (AGE + SEX + DIET + SPORT)^2,data = DATA.ANA, nvmax=(ncol(DATA.ANA)-1), method ="forward")
summary (regfit.fwd)
regfit.bwd <- regsubsets(HOR_LEVEL~.-ID + I(AGE^2) + (AGE + SEX + DIET + SPORT)^2, data = DATA.ANA, nvmax = (ncol(DATA.ANA)-1), method = "backward")
summary(regfit.bwd)

par(mfrow = c(2,2))
#plots for forward step wise
plot(regfit.fwd, scale="r2")
plot(regfit.fwd, scale="adjr2")
plot(regfit.fwd, scale="Cp")
plot(regfit.fwd, scale="bic")

#plots for backwards step wise
plot(regfit.bwd, scale="r2")
plot(regfit.bwd, scale="adjr2")
plot(regfit.bwd, scale="Cp")
plot(regfit.bwd, scale="bic")

library(robustbase)
rob_fit <- lmrob(HOR_LEVEL~.-ID + I(AGE^2) + (AGE + SEX + DIET + SPORT)^2, method = "MM",
               data = DATA.ANA, setting = "KS2014")
summary(rob_fit)#c(44,242,364,385,405,479,500,645) are outliers
par(mfrow = c(3,2)) 
plot(rob_fit, sub = NULL)
w1 <- round(weights(rob_fit, type = 'r'), 3)

plot(DATA.ANA, col = ifelse(w1 >= 0.01, 1, 2))
