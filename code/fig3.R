library(survival)
library(survminer)
library(dplyr)
library(pracma)
library(lognorm)
library(tikzDevice)
library(gbutils)

clr1<-rgb(27/255,158/255,119/255)
clr2<-rgb(217/255,95/255,2/255)
clr3<-rgb(117/255,112/255,179/255)
clr4<-rgb(231/255,41/255,138/255)
clr5<-rgb(102/255,166/255,30/255)
clr6<-rgb(230/255,171/255,2/255)
clr7<-rgb(166/255,118/255,29/255)

clr1b<-rgb(27/255,158/255,119/255,0.2)
clr2b<-rgb(217/255,95/255,2/255,0.2)
clr3b<-rgb(117/255,112/255,179/255,0.2)
clr4b<-rgb(231/255,41/255,138/255,0.2)
clr5b<-rgb(102/255,166/255,30/255,0.2)
clr6b<-rgb(230/255,171/255,2/255,0.2)
clr7b<-rgb(166/255,118/255,29/255,0.2)

tikz("fig3.tex",width=6,height=5)
#pdf("fig3.pdf")

par(mfrow = c(1,1),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)

lambda_ = 1/60 # Weibull scale
k = 3 # Weibull shape
delta = 1 # Gamma variance
num = 10^7

## Empirical theta ##
U = runif(num, min = 0, max = 1)
U0 = rgamma(num,shape=1/delta,scale=1/delta)
T0 = (-log(U)/(lambda_*U0))^{1/k} 
F0 = ecdf(T0)
S0 = function(t){1-F0(t)}

eval_ts0 = quantile(F0,c(1e-4, seq(0.001,0.999,0.001), 1 - 1e-4))
eval_S0s = S0(eval_ts0)
inverse_S0 = function(t){interp1(eval_S0s, eval_ts0, t)}

theta <-function(t, S1_t){
  return(inverse_S0(S1_t)/t)
}

###################################
####### E[U_1]=(1/3)^(1/3) ########
###################################

mu = (1/3)**(1/3) ###=E[U_1]

### T0 ~ Weibull-Gamma, U_1 ~ Gamma, unit variance ###
rho_1 = 1 ###=var[U_1]

k1 = mu**2/rho_1 #shape
theta1= rho_1/mu #scale

U1 = rgamma(num, shape=k1, scale=theta1)

# print(paste0("var U: ", var(U1)))
# print(paste0("E[U]: ", mean(U1)))

T1 = T0*(1/U1)

print(round(c(mu, mean(T0)/mean(T1), exp(mean(T0)-mean(T1))),3))

F1 = ecdf(T1)
S1 = function(t){1-F1(t)}

eval_ts1 = sapply(seq(0.001,0.999,0.001), function(q){cdf2quantile(q, F1, interval = c(0.001, 1000))})
eval_S1s = S1(eval_ts1)
theta_ = mapply(theta, eval_ts1, eval_S1s)

plot((1-eval_S1s), theta_, type="l", xlab="$q$", ylab='$\\theta$', xlim=c(0,1), ylim=c(0,2.6), col=clr2, axes=F)
abline(h=(1/3)**(1/3),lwd=2,col=clr2b)
abline(h=1,lwd=2,col=rgb(220/255,220/255,220/255,0.2))

### T0 ~ Weibull-Gamma, U_1 ~ Gamma, small variance ###
rho_1 = 0.5  ###=var[U_1]

k1 = mu**2/rho_1 #shape
theta1= rho_1/mu #scale

U1 = rgamma(num, shape=k1, scale=theta1)

# print(paste0("var U: ", var(U1)))
# print(paste0("E[U]: ", mean(U1)))

T1 = T0*(1/U1)

print(round(c(mu, mean(T0)/mean(T1), exp(mean(T0)-mean(T1))),3))


#S1 = function(t){S0((mu_1)*t)*p_1 + S0((mu_2)*t)*p_2 + S0(t)*(1-p_1-p_2)}
F1 = ecdf(T1)
S1 = function(t){1-F1(t)}
eval_ts1 = sapply(seq(0.001,0.999,0.001), function(q){cdf2quantile(q, F1, interval = c(0.001, 1000))})
eval_S1s = S1(eval_ts1)
theta_ = mapply(theta, eval_ts1, eval_S1s)

lines((1-eval_S1s), theta_, col=clr2,lty=2)

### T0 ~ Weibull-Gamma, U_1 ~ Gamma, large variance ###
rho_1 = 2  ###=var[U_1]

k1 = mu**2/rho_1 #shape
theta1= rho_1/mu #scale

U1 = rgamma(num, shape=k1, scale=theta1)

# print(paste0("var U: ", var(U1)))
# print(paste0("E[U]: ", mean(U1)))

T1 = T0*(1/U1)

print(round(c(mu, mean(T0)/mean(T1), exp(mean(T0)-mean(T1))),3))

#S1 = function(t){S0((mu_1)*t)*p_1 + S0((mu_2)*t)*p_2 + S0(t)*(1-p_1-p_2)}
F1 = ecdf(T1)
S1 = function(t){1-F1(t)}
eval_ts1 = sapply(seq(0.001,0.999,0.001), function(q){cdf2quantile(q, F1, interval = c(0.001, 1000))})
eval_S1s = S1(eval_ts1)

theta_ = mapply(theta, eval_ts1, eval_S1s)

lines((1-eval_S1s), theta_, col=clr2,lty=3)

##############################
####### E[U_1]=3^(1/3) #######
##############################

mu = 3**(1/3) ###=E[U_1]
p_1 = 0.05
mu_1 = 0.5

### T0 ~ Weibull-Gamma, U_1 ~ Gamma, unit variance ###
rho_1 = 1 ###=var[U_1]

k1 = mu**2/rho_1 #shape
theta1= rho_1/mu #scale

U1 = rgamma(num, shape=k1, scale=theta1)

# print(paste0("var U: ", var(U1)))
# print(paste0("E[U]: ", mean(U1)))

T1 = T0*(1/U1)

print(round(c(mu, mean(T0)/mean(T1), exp(mean(T0)-mean(T1))),3))


F1 = ecdf(T1)
S1 = function(t){1-F1(t)}

eval_ts1 = sapply(seq(0.001,0.999,0.001), function(q){cdf2quantile(q, F1, interval = c(0.001, 1000))})
eval_S1s = S1(eval_ts1)
theta_ = mapply(theta, eval_ts1, eval_S1s)

lines((1-eval_S1s), theta_, col=clr1)
abline(h=3**(1/3),lwd=2,col=clr1b)

### T0 ~ Weibull-Gamma, U_1 ~ Gamma, small variance ###
rho_1 = 0.5 ###=var[U_1]

k1 = mu**2/rho_1 #shape
theta1= rho_1/mu #scale

U1 = rgamma(num, shape=k1, scale=theta1)

# print(paste0("var U: ", var(U1)))
# print(paste0("E[U]: ", mean(U1)))

T1 = T0*(1/U1)

print(round(c(mu, mean(T0)/mean(T1), exp(mean(T0)-mean(T1))),3))


F1 = ecdf(T1)
S1 = function(t){1-F1(t)}

#S1 = function(t){S0((1/mu_1)*t)*p_1 + S0((1/mu_2)*t)*p_2 + S0(t)*(1-p_1-p_2)}
eval_ts1 = sapply(seq(0.001,0.999,0.001), function(q){cdf2quantile(q, F1, interval = c(0.001, 1000))})
eval_S1s = S1(eval_ts1)
theta_ = mapply(theta, eval_ts1, eval_S1s)

lines((1-eval_S1s), theta_, col=clr1,lty=2)

### T0 ~ Weibull-Gamma, U_1 ~ Gamma, large variance ###
rho_1 = 2  ###=var[1/U_1]

k1 = mu**2/rho_1 #shape
theta1= rho_1/mu #scale

U1 = rgamma(num, shape=k1, scale=theta1)

# print(paste0("var U: ", var(U1)))
# print(paste0("E[U]: ", mean(U1)))

T1 = T0*(1/U1)

print(round(c(mu, mean(T0)/mean(T1), exp(mean(T0)-mean(T1))),3))

F1 = ecdf(T1)
S1 = function(t){1-F1(t)}

eval_ts1 = sapply(seq(0.001,0.999,0.001), function(q){cdf2quantile(q, F1, interval = c(0.001, 1000))})
eval_S1s = S1(eval_ts1)
theta_ = mapply(theta, eval_ts1, eval_S1s)

lines((1-eval_S1s), theta_, col=clr1,lty=3)

axis(side = 1,
     labels = TRUE)

axis(side = 2,
     labels = TRUE)
box(which = "plot", bty = "l")

title(xlab = "$F_{T^a}(t)$",
      ylab = "$\\theta(t)$",
      outer = TRUE, line = 3)

legend(0.6,0.75,legend=c('$0.5$','$1$','$2$'),lty=c(2,1,3),seg.len=4,title="$\\mathrm{Var}[U_1]$", cex=0.9)
dev.off()








