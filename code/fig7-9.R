library("extRemes")
library("survival")
library("statmod")
library("gamma")
library("pracma")
library("tikzDevice")


library(copula)
library(gumbel)
library(GoFKernel)
library(stabledist)

clr1<-rgb(27/255,158/255,119/255)
clr2<-rgb(217/255,95/255,2/255)
clr3<-rgb(117/255,112/255,179/255)
clr4<-rgb(231/255,41/255,138/255)

simcop<-function(tau0,tau1,beta_LA){
  
  #Gaussian Copula
  rho<- c(0, sin(pi*tau0/2), sin(pi*tau1/2))
  
  nc<- normalCopula(rho,dim=3,dispstr = "un")
  NU  <- rCopula(num, copula = nc)
  
  U0<-sapply(NU[,1],function(x){qgamma(x, shape=1/theta0,scale=theta0)})
  
  p_1 = 0.05
  mu_1 = 0.5
  
  p_2 <- (-1 + HR0 + p_1 - mu_1*p_1)^2/(1+theta1-2*HR0+HR0^2-p_1+2*mu_1*p_1-mu_1^2*p_1)
  mu_2 <- (theta1-HR0+HR0^2+mu_1*p_1-mu_1^2*p_1)/(-1+HR0+p_1-mu_1*p_1)
  
  U1<-sapply(NU[,2],function(x){if(x<p_1){return(mu_1)}else{
    
    if(x<(1-p_2)){return(1)}
    
    else{
      
      return(mu_2)}
    
  }
    
  })
  
  pL<-0.5
  L<-sapply(NU[,3],function(x){ifelse(x<pL,1,0)})
  #print('hei')
  p_0 = mean(L)
  
  NT<-runif(num,0,1)
  
  T0<-(-(60/U0)*log(NT))^(1/3)
  T1<-T0/U1
  
  pA<- 0.5 + beta_LA*(2*L-1)
  A<-sapply(pA, function(x){rbinom(1,1,x)}) 
  
  T<-A*T1+(1-A)*T0
  
  df = data.frame(time=T, status=1, log_time=log(T), a=A, l=L)
  return(df)
}

AF_causal_empirical <- function() { 
  #########################
  ### unconfounded data ###
  #########################
  
  df = sim_data()
  
  survobj = Surv(time = df$time, event = df$status)
  fit = survfit(survobj~a, data = df)
  
  fit_0 = fit[1]
  fit_a = fit[2]
  
  S0<-function(t){summary(fit_0, times=t, extend=TRUE)$surv}
  S1<-function(t){summary(fit_a, times=t, extend=TRUE)$surv}
  
  max_ts = max(df$time)
  
  eval_ts = seq(0.5, max_ts, 1e-1)
  eval_S0s = S0(eval_ts)
  eval_S1s = S1(eval_ts)
  
  inverse_S0 = function(t){interp1(eval_S0s, eval_ts, t)}
  
  theta <-function(t, S1_t){
    return(inverse_S0(S1_t)/t)
  }
  
  max_idx = tail(which(eval_S1s > tail(eval_S0s, 1)), 1)
  theta_ = mapply(theta, eval_ts[1:max_idx], eval_S1s[1:max_idx])
  
  return(
    list(
      "ts"=eval_ts[1:max_idx],
      "theta"=theta_
    )
  )
}

AF_causal <- function(ts) {
  lambda_ = 1/60 #scale
  k = 3 #shape
  
  p_1 = 0.05
  mu_1 = 0.5
  p_2 <- (-1 + HR0 + p_1 - mu_1*p_1)^2/(1+theta1-2*HR0+HR0^2-p_1+2*mu_1*p_1-mu_1^2*p_1)
  mu_2 = (theta1-HR0+HR0^2+mu_1*p_1-mu_1^2*p_1)/(-1+HR0+p_1-mu_1*p_1)
  
  S0_gamma = function(t){(1 + theta0*((t**k)*lambda_))**(-1/theta0)}
  inverse_S0_gamma = function(p){((p**(-theta0) - 1)/(lambda_*theta0))**(1/k)}
  
  theta_BHN = function(t){inverse_S0_gamma(S0_gamma((mu_1)*t)*p_1 + S0_gamma((mu_2)*t)*p_2 + S0_gamma(t)*(1-p_1-p_2))/t}
  
  theta_BHN_ = sapply(ts, function(t){theta_BHN(t)})
  
  return(theta_BHN_)
}

num = 1e5
theta0 = 1 
theta1 = 1
HR0 = 3**(1/3)

beta_LA_values = c(0, 0.25, 0.45)
tau0_values = c(0, 0.25, 0.45)
tau1_values = c(0, 0.25, 0.45)

num_beta_LA = length(beta_LA_values)
num_tau0 = length(beta_LA_values)
num_tau1 = length(beta_LA_values)
#uncorr_res = AF_causal_empirical()
#eval_ts_uncorr = uncorr_res$ts
#theta_uncorr = uncorr_res$theta

par(mfrow=c(3, 3))

# q0 <- function(ts) {
#   lambda_ = 1/60 #scale
#   k = 3 #shape
#   
#   S0_gamma = function(t){(1 + theta0*((t**k)*lambda_))**(-1/theta0)}
#   
#   q_ = sapply(ts, function(t){1-S0_gamma(t)})
#   
#   return(q_)
# }

q0 <- function(ts) {
  lambda_ = 1/60 #scale
  k = 3 #shape
  
  S0_gamma = function(t){(1 + theta0*((t**k)*lambda_))**(-1/theta0)}
  
  q_ = sapply(ts, function(t){1-S0_gamma(HR0*t)})
  
  return(q_)
}

#pdf("confounders_XXX.pdf")
tikz("confounders_XXX.tex",width=6,height=3)
par(mfrow = c(3,3),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)

for (beta_idx in 1:num_beta_LA) {
  print(beta_idx)
  for (tau0_idx in 1:num_tau0) {
    print(tau0_idx)
    #print(c(tau0_idx, beta_idx))  
    
    tau0 = tau0_values[tau0_idx]
    #tau1 = 0
    beta_LA = beta_LA_values[beta_idx]
    
    print(c(tau0,0,beta_LA))
    
    df = simcop(tau0, tau0, beta_LA)
    
    
    p_0 = mean(df$l==0)
    p_1 = mean(df$l==1)
    
    survobj = Surv(time = df$time, event = df$status)
    fit_marg = survfit(survobj~a, data = df)
    fit_cond = survfit(survobj~a+l, data = df)
    
    fit_0 = fit_marg[1]
    fit_a = fit_marg[2]
    
    fit_00 = fit_cond[1]
    fit_10 = fit_cond[3]
    fit_01 = fit_cond[2]
    fit_11 = fit_cond[4]
    
    S0<-function(t){summary(fit_0, times=t, extend=TRUE)$surv}
    S1<-function(t){summary(fit_a, times=t, extend=TRUE)$surv}
    
    S00<-function(t){summary(fit_00, times=t, extend=TRUE)$surv}
    S10<-function(t){summary(fit_10, times=t, extend=TRUE)$surv}
    S01<-function(t){summary(fit_01, times=t, extend=TRUE)$surv}
    S11<-function(t){summary(fit_11, times=t, extend=TRUE)$surv}
    
    S0_adj<-function(t){S00(t)*p_0 + S01(t)*p_1}
    S1_adj<-function(t){S10(t)*p_0 + S11(t)*p_1}
    
    max_ts = max(df$time)
    
    eval_ts = seq(0.5, max_ts, 1e-1)
    eval_S0s = S0(eval_ts)
    eval_S1s = S1(eval_ts)
    
    eval_S0s_adj = S0_adj(eval_ts)
    eval_S1s_adj = S1_adj(eval_ts)
    
    eval_S00s = S00(eval_ts)
    eval_S10s = S10(eval_ts)
    eval_S01s = S01(eval_ts)
    eval_S11s = S11(eval_ts)
    
    inverse_S0 = function(t){interp1(eval_S0s, eval_ts, t)}
    inverse_S00 = function(t){interp1(eval_S00s, eval_ts, t)}
    inverse_S01 = function(t){interp1(eval_S01s, eval_ts, t)}
    inverse_S0_adj = function(t){interp1(eval_S0s_adj, eval_ts, t)}
    
    theta_marg <-function(t, S1_t){
      return(inverse_S0(S1_t)/t)
    }
    
    theta_adj <-function(t, S1_t){
      return(inverse_S0_adj(S1_t)/t)
    }
    
    max_idx_marg = tail(which(eval_S1s > tail(eval_S0s, 1)), 1)
    theta_marg_ = mapply(theta_marg, eval_ts[1:max_idx_marg], eval_S1s[1:max_idx_marg])
    
    max_idx_adj = tail(which(eval_S1s_adj > tail(eval_S0s_adj, 1)), 1)
    theta_adj_ = mapply(theta_adj, eval_ts[1:max_idx_adj], eval_S1s_adj[1:max_idx_adj])
    
    AF_causal_ = sapply(eval_ts,function(ts){AF_causal(ts)})
    
    #plot(eval_ts_uncorr, theta_uncorr, type="l", xlab="t", ylab='y', ylim=c(0.5,2.5), col=clr1, lty=2)
    plot(q0(eval_ts[1:max_idx_marg]), theta_marg_, type="l", xlab="", ylab="Estimand", xlim=c(0.15,0.95), ylim=c(0.7,1.6), col=clr1, lty=2, axes=F)
    lines(q0(eval_ts[1:max_idx_adj]), theta_adj_, col=clr2, lty=2)
    lines(q0(eval_ts), AF_causal_, col=clr3, lty=2)
    
    if (tau0_idx==1 & beta_idx==3) {
      axis(side = 1,
           labels = TRUE)
      axis(side = 2,
           labels = TRUE)
      box(which = "plot", bty = "l")
    } else if (tau0_idx==1) {
      axis(side = 1,
           labels = FALSE)
      axis(side = 2,
           labels = TRUE)
      box(which = "plot", bty = "l")
    } else if (beta_idx==3) {
      axis(side = 1,
           labels = TRUE)
      axis(side = 2,
           labels = FALSE)
    } else {
      axis(side = 1,
           labels = FALSE)
      axis(side = 2,
           labels = FALSE)
      box(which = "plot", bty = "l")
    }
  }
}
title(xlab = "$tau_0$",
      ylab = "$beta_LA$",
      outer = TRUE, line = 3)

legend("topright",c("$theta_m$","$theta_adj$", "$theta$"), lty=c(2,2,2), col=c(clr1,clr2,clr3),seg.len=4,title="Estimands", cex=0.5)
dev.off()
