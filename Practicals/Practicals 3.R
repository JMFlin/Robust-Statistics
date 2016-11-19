library(MASS)
#-------------c
k <- 1.34

plot(NA, xlim=c(-8,8), ylim=c(-4,4), xaxs = "i")

curve(sin(x/k), from = -k*pi, to= k*pi, add = T)
curve(0*x + 0, from = -8, to = -k*pi, add = T)
curve(0*x + 0, from = k*pi, to= 8, add = T)

#-------------d
andrewM <- function(x, c = 1, eps = 1e-06, init = 0){
  mu <- init
  diff <- Inf
  while (diff > eps){
    r <- x - mu
    psi.r <- r
    psi.r <- ifelse(abs(psi.r) <= c*pi, psi.r, ifelse(abs(psi.r) > c*pi,0,psi.r))
    #psi.r <- ifelse(,,)
    w <- psi.r/(r + eps)
    mu.new  <- sum(w * x)/sum(w)
    diff <- abs(mu - mu.new)
    mu <- mu.new
  }
  mu.new
} 

#-------------e
HuberM <- function(x, c = 1, init = 0, eps = 1e-06){
  n <- length(x)
  mu <- median(x)
  s <- mad(x)
  repeat{
    xx <- pmin(pmax(mu - c * s, x), mu + c * s)
    mu1 <- sum(xx)/n
    if (abs(mu - mu1) < eps * s){
      break
    }
    mu <- mu1
  }
  mu
}

#-------------f
my_data <- read.csv("HTEST.txt", header = T)

andrewM(my_data$x)
HuberM(my_data$x)


#-------------2a

my_function <- function(x, mu, sigma, k){
  
}

#-------------2b

plot(my_data$x)

x_mean <- mean(my_data$x)
s_sd <- sd(my_data$x)

error <- qt(0.975,df=length(my_data$x)-1)*sd(my_data$x)/sqrt(length(my_data$x))
left <- mean(my_data$x)-error
right <- mean(my_data$x)+error
cbind(left, right)

#-------------2c

huber1 <- huber(my_data$x, k=1)

error <- qt(0.975,df=length(my_data$x)-1)*huber1$s/sqrt(length(my_data$x))
left <- huber1$mu-error
right <- huber1$mu+error
cbind(left, right)

huber1.5 <- huber(my_data$x, k=1.5)

error <- qt(0.975,df=length(my_data$x)-1)*huber1.5$s/sqrt(length(my_data$x))
left <- huber1.5$mu-error
right <- huber1.5$mu+error
cbind(left, right)

huber2 <- huber(my_data$x, k=2)

error <- qt(0.975,df=length(my_data$x)-1)*huber2$s/sqrt(length(my_data$x))
left <- huber2$mu-error
right <- huber2$mu+error
cbind(left, right)

#-------------2d

t.test(my_data$x, conf.level = 0.95)

#-------------2e

t.test(my_data$x, mu = huber1$mu, s = huber1$s, conf.level = 0.95)

t.test(my_data$x, mu = huber1.5$mu, s = huber1.5$s, conf.level = 0.95)

t.test(my_data$x, mu = huber2$mu, s = huber2$s, conf.level = 0.95)
