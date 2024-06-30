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

tikz("fig2.tex",width=6,height=3)
#pdf("fig2.pdf")
par(mfrow = c(1,2),
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
p_1 = 0.7
mu_1 = 0.3

### T0 ~ Weibull-Gamma, U_1 ~ BHN[unit variance rho_1] ###
rho_1 = 1 ###=var[U_1]

p_2 = (-1 + mu + p_1 - mu_1*p_1)^2/(1+rho_1-2*mu+mu^2-p_1+2*mu_1*p_1-mu_1^2*p_1)
mu_2 = (rho_1-mu+mu^2+mu_1*p_1-mu_1^2*p_1)/(-1+mu+p_1-mu_1*p_1)

U1 = sapply(runif(num,0,1),function(x){if(x<p_1){return(mu_1)}else{
  if(x<(1-p_2)){return(1)}else{return(mu_2)}}})

print(paste0("var U: ", var(U1)))
print(paste0("E[U]: ", mean(U1)))

print(paste0("p1 + p2: ", p_1 + p_2))
print(paste0("mu1: ", mu_1))
print(paste0("mu2: ", mu_2))
print(paste0("p1: ", p_1))
print(paste0("p2: ", p_2))


T1 = T0*(1/U1)

print(round(c(mu, mean(T0)/mean(T1), exp(mean(T0)-mean(T1))),3))

F1 = ecdf(T1)
S1 = function(t){1-F1(t)}

eval_ts1 = sapply(seq(0.001,0.999,0.001), function(q){cdf2quantile(q, F1, interval = c(0.001, 1000))})
eval_S1s = S1(eval_ts1)
theta_ = mapply(theta, eval_ts1, eval_S1s)

plot((1-eval_S1s), theta_, type="l", xlab="$q$", ylab='$\\theta$', xlim=c(0,1), ylim=c(0,2.1), col=clr2, axes=F)
#plot(eval_ts1, theta_, type="l", xlab="$q$", ylab='$\\theta$', xlim=c(0,max(eval_ts1)), ylim=c(0,2.1), col=clr2, axes=F)
abline(h=(1/3)**(1/3),lwd=2,col=clr2b)
abline(h=1,lwd=2,col=rgb(220/255,220/255,220/255,0.2))

# ### T0 ~ Weibull-Gamma, U_1 ~ BHN[small variance rho_1] ###
# rho_1 = 0.5  ###=var[U_1]
# 
# p_2 = (-1 + mu + p_1 - mu_1*p_1)^2/(1+rho_1-2*mu+mu^2-p_1+2*mu_1*p_1-mu_1^2*p_1)
# mu_2 = (rho_1-mu+mu^2+mu_1*p_1-mu_1^2*p_1)/(-1+mu+p_1-mu_1*p_1)
# 
# U1 = sapply(runif(num,0,1),function(x){if(x<p_1){return(mu_1)}else{
#   if(x<(1-p_2)){return(1)}else{return(mu_2)}}})
# 
# print(paste0("var U: ", var(U1)))
# print(paste0("E[U]: ", mean(U1)))
# 
# print(paste0("p1 + p2: ", p_1 + p_2))
# print(paste0("mu1: ", mu_1))
# print(paste0("mu2: ", mu_2))
# 
# T1 = T0*(1/U1)
# 
# print(c(mu, mean(T0/T1)))
# 
# F1 = ecdf(T1)
# S1 = function(t){1-F1(t)}
# eval_ts1 = sapply(seq(0.001,0.999,0.001), function(q){cdf2quantile(q, F1, interval = c(0.001, 1000))})
# eval_S1s = S1(eval_ts1)
# theta_ = mapply(theta, eval_ts1, eval_S1s)
# 
# lines((1-eval_S1s), theta_, col=clr2,lty=2)
# 
# ### T0 ~ Weibull-Gamma, U_1 ~ BHN[large variance rho_1] ###
# rho_1 = 2  ###=var[U_1]
# 
# p_2 = (-1 + mu + p_1 - mu_1*p_1)^2/(1+rho_1-2*mu+mu^2-p_1+2*mu_1*p_1-mu_1^2*p_1)
# mu_2 = (rho_1-mu+mu^2+mu_1*p_1-mu_1^2*p_1)/(-1+mu+p_1-mu_1*p_1)
# 
# U1 = sapply(runif(num,0,1),function(x){if(x<p_1){return(mu_1)}else{
#   if(x<(1-p_2)){return(1)}else{return(mu_2)}}})
# 
# print(paste0("var U: ", var(U)))
# print(paste0("E[U]: ", mean(U)))
# 
# print(paste0("p1 + p2: ", p_1 + p_2))
# print(paste0("mu1: ", mu_1))
# print(paste0("mu2: ", mu_2))
# 
# T1 = T0*(1/U1)
# 
# print(c(mu, mean(T0/T1)))
# 
# F1 = ecdf(T1)
# S1 = function(t){1-F1(t)}
# eval_ts1 = sapply(seq(0.001,0.999,0.001), function(q){cdf2quantile(q, F1, interval = c(0.001, 1000))})
# eval_S1s = S1(eval_ts1)
# 
# theta_ = mapply(theta, eval_ts1, eval_S1s)
# 
# lines((1-eval_S1s), theta_, col=clr2,lty=3)

##############################
####### E[U_1]=3^(1/3) #######
##############################

mu = 3**(1/3) ###=E[U_1]
p_1 = 0.05
mu_1 = 0.5

### T0 ~ Weibull-Gamma, U_1 ~ BHN[unit variance rho_1] ###
rho_1 = 1 ###=var[U_1]

p_2 = (-1 + mu + p_1 - mu_1*p_1)^2/(1+rho_1-2*mu+mu^2-p_1+2*mu_1*p_1-mu_1^2*p_1)
mu_2 = (rho_1-mu+mu^2+mu_1*p_1-mu_1^2*p_1)/(-1+mu+p_1-mu_1*p_1)

U1 = sapply(runif(num,0,1),function(x){if(x<p_1){return(mu_1)}else{
  if(x<(1-p_2)){return(1)}else{return(mu_2)}}})

# print(paste0("var U: ", var(U1)))
# print(paste0("E[U]: ", mean(U1)))
# 
# print(paste0("p1 + p2: ", p_1 + p_2))
# print(paste0("mu1: ", mu_1))
# print(paste0("mu2: ", mu_2))
# print(paste0("p1: ", p_1))
# print(paste0("p2: ", p_2))


T1 = T0*(1/U1)

print(round(c(mu, mean(T0)/mean(T1), exp(mean(T0)-mean(T1))),3))

F1 = ecdf(T1)
S1 = function(t){1-F1(t)}

eval_ts1 = sapply(seq(0.001,0.999,0.001), function(q){cdf2quantile(q, F1, interval = c(0.001, 1000))})
eval_S1s = S1(eval_ts1)
theta_ = mapply(theta, eval_ts1, eval_S1s)

lines((1-eval_S1s), theta_, col=clr1)
#lines(eval_ts1, theta_, col=clr1)
abline(h=3**(1/3),lwd=2,col=clr1b)

# ### T0 ~ Weibull-Gamma, U_1 ~ BHN[small variance rho_1] ###
# rho_1 = 0.5 ###=var[U_1]
# 
# p_2 = (-1 + mu + p_1 - mu_1*p_1)^2/(1+rho_1-2*mu+mu^2-p_1+2*mu_1*p_1-mu_1^2*p_1)
# mu_2 = (rho_1-mu+mu^2+mu_1*p_1-mu_1^2*p_1)/(-1+mu+p_1-mu_1*p_1)
# 
# U1 = sapply(runif(num,0,1),function(x){if(x<p_1){return(mu_1)}else{
#   if(x<(1-p_2)){return(1)}else{return(mu_2)}}})
# 
# print(paste0("var U: ", var(U)))
# print(paste0("E[U]: ", mean(U)))
# 
# print(paste0("p1 + p2: ", p_1 + p_2))
# print(paste0("mu1: ", mu_1))
# print(paste0("mu2: ", mu_2))
# 
# T1 = T0*(1/U1)
# 
# print(c(mu, mean(T0/T1)))
# 
# F1 = ecdf(T1)
# S1 = function(t){1-F1(t)}
# 
# eval_ts1 = sapply(seq(0.001,0.999,0.001), function(q){cdf2quantile(q, F1, interval = c(0.001, 1000))})
# eval_S1s = S1(eval_ts1)
# theta_ = mapply(theta, eval_ts1, eval_S1s)
# 
# lines((1-eval_S1s), theta_, col=clr1,lty=2)

# ### T0 ~ Weibull-Gamma, U_1 ~ BHN[large variance rho_1] ###
# rho_1 = 2  ###=var[1/U_1]
# 
# p_2 = (-1 + mu + p_1 - mu_1*p_1)^2/(1+rho_1-2*mu+mu^2-p_1+2*mu_1*p_1-mu_1^2*p_1)
# mu_2 = (rho_1-mu+mu^2+mu_1*p_1-mu_1^2*p_1)/(-1+mu+p_1-mu_1*p_1)
# 
# U1 = sapply(runif(num,0,1),function(x){if(x<p_1){return(mu_1)}else{
#   if(x<(1-p_2)){return(1)}else{return(mu_2)}}})
# 
# print(paste0("var U: ", var(U1)))
# print(paste0("E[U]: ", mean(U1)))
# 
# print(paste0("p1 + p2: ", p_1 + p_2))
# print(paste0("mu1: ", mu_1))
# print(paste0("mu2: ", mu_2))
# 
# 
# T1 = T0*(1/U1)
# 
# print(c(mu, mean(T0/T1)))
# 
# F1 = ecdf(T1)
# S1 = function(t){1-F1(t)}
# 
# eval_ts1 = sapply(seq(0.001,0.999,0.001), function(q){cdf2quantile(q, F1, interval = c(0.001, 1000))})
# eval_S1s = S1(eval_ts1)
# theta_ = mapply(theta, eval_ts1, eval_S1s)
# 
# lines((1-eval_S1s), theta_, col=clr1,lty=3)

axis(side = 1,
     labels = TRUE)

axis(side = 2,
     labels = TRUE)
box(which = "plot", bty = "l")

legend(0.6,2,legend=c("$(1/3)^{1/3}$","$3^{1/3}$"), col=c(clr2,clr1), lty=c(1,1) ,title="$E[U_1]$")


###################################
############ Mixture ##############
###################################
###################################
####### E[U_1]=(1/3)^(1/3) ########
###################################
mu = (1/3)**(1/3) ###=E[U_1]
p_1 = 0.7
mu_1 = 0.3

### C. T~ Weibull mixture, U_1 ~ BHN[with unit variance rho_1] ###
n = 10^7

labels = sample(c(1,2),n,TRUE)
means = c(1,10)

T0 = rweibull(n, shape =2, scale = means[labels]/gamma(1+1/2))

rho_1 = 1 

p_2 = (-1 + mu + p_1 - mu_1*p_1)^2/(1+rho_1-2*mu+mu^2-p_1+2*mu_1*p_1-mu_1^2*p_1)
mu_2 = (rho_1-mu+mu^2+mu_1*p_1-mu_1^2*p_1)/(-1+mu+p_1-mu_1*p_1)

# print(paste0("p1 + p2: ", p_1 + p_2))
# print(paste0("mu1: ", mu_1))
# print(paste0("mu2: ", mu_2))

U1 = sapply(runif(n,0,1),function(x){if(x<p_1){return(mu_1)}else{
  if(x<(1-p_2)){return(1)}else{return(mu_2)}}})
T1 = T0*(1/U1)

print(round(c(mu, mean(T0)/mean(T1), exp(mean(T0)-mean(T1))),3))

F0 = ecdf(T0)
F1 = ecdf(T1)
S0 = function(t){1-F0(t)}
S1 = function(t){1-F1(t)}

eval_ts0 = quantile(T0,seq(0.001,0.999,0.001))
eval_ts1 = quantile(T1,seq(0.001,0.999,0.001))
eval_S0s = S0(eval_ts0)
eval_S1s = S1(eval_ts1)
inverse_S0 = function(t){interp1(eval_S0s, eval_ts0, t)}

theta_ = mapply(theta, eval_ts1, eval_S1s)

#plot(eval_ts1, theta_, type="l", xlab="$q$", ylab='$\\theta$', xlim=c(0,max(eval_ts1)), ylim=c(0,2.1), col=clr2, axes=F)
plot((1-eval_S1s), theta_, type="l", xlab="$q$", ylab='$\\theta$', xlim=c(0,1), ylim=c(0,2.1), col=clr2, axes=F)
abline(h=(1/3)**(1/3),lwd=2,col=clr2b)
abline(h=1,lwd=2,col=rgb(220/255,220/255,220/255,0.2))

# ### C. T~ Weibull mixture, U_1 ~ BHN[with small variance rho_1] ###
# rho_1 = 0.5 
# 
# p_2 = (-1 + mu + p_1 - mu_1*p_1)^2/(1+rho_1-2*mu+mu^2-p_1+2*mu_1*p_1-mu_1^2*p_1)
# mu_2 = (rho_1-mu+mu^2+mu_1*p_1-mu_1^2*p_1)/(-1+mu+p_1-mu_1*p_1)
# 
# # print(paste0("p1 + p2: ", p_1 + p_2))
# # print(paste0("mu1: ", mu_1))
# # print(paste0("mu2: ", mu_2))
# 
# U1 = sapply(runif(n,0,1),function(x){if(x<p_1){return(mu_1)}else{
#   if(x<(1-p_2)){return(1)}else{return(mu_2)}}})
# T1 = T0*(1/U1)
# 
# print(c(mu, mean(T0/T1)))
# 
# F0 = ecdf(T0)
# F1 = ecdf(T1)
# S0 = function(t){1-F0(t)}
# S1 = function(t){1-F1(t)}
# 
# eval_ts0 = quantile(T0,seq(0.01,0.99,0.001))
# eval_ts1 = quantile(T1,seq(0.01,0.99,0.001))
# eval_S0s = S0(eval_ts0)
# eval_S1s = S1(eval_ts1)
# inverse_S0 = function(t){interp1(eval_S0s, eval_ts0, t)}
# 
# theta_ = mapply(theta, eval_ts1, eval_S1s)
# lines((1-eval_S1s), theta_, col=clr2, lty=2)
# 
# ### C. T~ Weibull mixture, U_1 ~ BHN[with large variance rho_1] ###
# rho_1 = 5 
# 
# p_2 = (-1 + mu + p_1 - mu_1*p_1)^2/(1+rho_1-2*mu+mu^2-p_1+2*mu_1*p_1-mu_1^2*p_1)
# mu_2 = (rho_1-mu+mu^2+mu_1*p_1-mu_1^2*p_1)/(-1+mu+p_1-mu_1*p_1)
# 
# # print(paste0("p1 + p2: ", p_1 + p_2))
# # print(paste0("mu1: ", mu_1))
# # print(paste0("mu2: ", mu_2))
# 
# U1 = sapply(runif(n,0,1),function(x){if(x<p_1){return(mu_1)}else{
#   if(x<(1-p_2)){return(1)}else{return(mu_2)}}})
# T1 = T0*(1/U1)
# 
# print(c(mu, mean(T0/T1)))
# 
# F0 = ecdf(T0)
# F1 = ecdf(T1)
# S0 = function(t){1-F0(t)}
# S1 = function(t){1-F1(t)}
# 
# eval_ts0 = quantile(T0,seq(0.01,0.99,0.001))
# eval_ts1 = quantile(T1,seq(0.01,0.99,0.001))
# eval_S0s = S0(eval_ts0)
# eval_S1s = S1(eval_ts1)
# inverse_S0 = function(t){interp1(eval_S0s, eval_ts0, t)}
# 
# theta_ = mapply(theta, eval_ts1, eval_S1s)
# lines((1-eval_S1s), theta_, col=clr2, lty=3)

##############################
####### E[U_1]=3^(1/3) #######
##############################

mu = 3**(1/3) ###=E[U_1]
p_1 = 0.05
mu_1 = 0.5

### C. T~ Weibull mixture, U_1 ~ BHN[with unit variance rho_1] ###
rho_1 = 1 

p_2 = (-1 + mu + p_1 - mu_1*p_1)^2/(1+rho_1-2*mu+mu^2-p_1+2*mu_1*p_1-mu_1^2*p_1)
mu_2 = (rho_1-mu+mu^2+mu_1*p_1-mu_1^2*p_1)/(-1+mu+p_1-mu_1*p_1)

# print(paste0("p1 + p2: ", p_1 + p_2))
# print(paste0("mu1: ", mu_1))
# print(paste0("mu2: ", mu_2))

U1 = sapply(runif(n,0,1),function(x){if(x<p_1){return(mu_1)}else{
  if(x<(1-p_2)){return(1)}else{return(mu_2)}}})
T1 = T0*(1/U1)

print(round(c(mu, mean(T0)/mean(T1), exp(mean(T0)-mean(T1))),3))

F0 = ecdf(T0)
F1 = ecdf(T1)
S0<-function(t){1-F0(t)}
S1<-function(t){1-F1(t)}

eval_ts0 = quantile(T0,seq(0.01,0.99,0.001))
eval_ts1 = quantile(T1,seq(0.01,0.99,0.001))
eval_S0s = S0(eval_ts0)
eval_S1s = S1(eval_ts1)

inverse_S0 = function(t){interp1(eval_S0s, eval_ts0, t)}

theta_ = mapply(theta, eval_ts1, eval_S1s)
lines((1-eval_S1s), theta_, col=clr1)
#lines(eval_ts1, theta_, col=clr1)
abline(h=3**(1/3),lwd=2,col=clr1b)

# ### C. T~ Weibull mixture, U_1 ~ BHN[with small variance rho_1] ###
# 
# rho_1 = 5 ###=var[U_1]
# 
# p_2 = (-1 + mu + p_1 - mu_1*p_1)^2/(1+rho_1-2*mu+mu^2-p_1+2*mu_1*p_1-mu_1^2*p_1)
# mu_2 = (rho_1-mu+mu^2+mu_1*p_1-mu_1^2*p_1)/(-1+mu+p_1-mu_1*p_1)
# 
# # print(paste0("p1 + p2: ", p_1 + p_2))
# # print(paste0("mu1: ", mu_1))
# # print(paste0("mu2: ", mu_2))
# 
# U1 = sapply(runif(n,0,1),function(x){if(x<p_1){return(mu_1)}else{
#   if(x<(1-p_2)){return(1)}else{return(mu_2)}}})
# T1 = T0*(1/U1)
# 
# print(c(mu, mean(T0/T1)))
# 
# F0 = ecdf(T0)
# F1 = ecdf(T1)
# S0 = function(t){1-F0(t)}
# S1 = function(t){1-F1(t)}
# 
# eval_ts0 = quantile(T0,seq(0.01,0.99,0.001))
# eval_ts0 = quantile(T1,seq(0.01,0.99,0.001))
# eval_S0s = S0(eval_ts0)
# eval_S1s = S1(eval_ts1)
# inverse_S0 = function(t){interp1(eval_S0s, eval_ts0, t)}
# 
# theta_ = mapply(theta, eval_ts1, eval_S1s)
# 
# lines((1-eval_S1s), theta_, col=clr1, lty=2)
# 
# ### C. T~ Weibull mixture, U_1 ~ BHN[with large variance rho_1] ###
# 
# rho_1 = 1.1 ###=var[U_1]
# 
# p_2 = (-1 + mu + p_1 - mu_1*p_1)^2/(1+rho_1-2*mu+mu^2-p_1+2*mu_1*p_1-mu_1^2*p_1)
# mu_2 = (rho_1-mu+mu^2+mu_1*p_1-mu_1^2*p_1)/(-1+mu+p_1-mu_1*p_1)
# 
# # print(paste0("p1 + p2: ", p_1 + p_2))
# # print(paste0("mu1: ", mu_1))
# # print(paste0("mu2: ", mu_2))
# 
# U1 = sapply(runif(n,0,1),function(x){if(x<p_1){return(mu_1)}else{
#   if(x<(1-p_2)){return(1)}else{return(mu_2)}}})
# T1 = T0*(1/U1)
# 
# print(c(mu, mean(T0/T1)))
# 
# F0 = ecdf(T0)
# F1 = ecdf(T1)
# S0 = function(t){1-F0(t)}
# S1 = function(t){1-F1(t)}
# 
# eval_ts0 = quantile(T0,seq(0.01,0.99,0.001))
# eval_ts1 = quantile(T1,seq(0.01,0.99,0.001))
# 
# eval_S0s = S0(eval_ts0)
# eval_S1s = S1(eval_ts1)
# inverse_S0 = function(t){interp1(eval_S0s, eval_ts0, t)}
# 
# theta_ = mapply(theta, eval_ts1, eval_S1s)
# lines((1-eval_S1s), theta_, col=clr1, lty=3)

axis(side = 1,
     labels = TRUE)

axis(side = 2,
     labels = FALSE)
box(which = "plot", bty = "l")

title(xlab = "$F_{T^a}(t)$",
      ylab = "$\\theta(t)$",
      outer = TRUE, line = 3)

dev.off()






