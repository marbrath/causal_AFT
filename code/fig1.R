clr1<-rgb(27/255,158/255,119/255)
clr2<-rgb(217/255,95/255,2/255)
clr3<-rgb(117/255,112/255,179/255)
clr4<-rgb(231/255,41/255,138/255)
clr5<-rgb(102/255,166/255,30/255)
clr6<-rgb(230/255,171/255,2/255)
clr7<-rgb(166/255,118/255,29/255)
clr8<-rgb(158/255,27/255,66/255) #cox purple
clr9<-rgb(30/255,159/255,80/255) #OH green


clr1b<-rgb(27/255,158/255,119/255,0.2)
clr2b<-rgb(217/255,95/255,2/255,0.2)
clr3b<-rgb(117/255,112/255,179/255,0.2)
clr4b<-rgb(231/255,41/255,138/255,0.2)
clr5b<-rgb(102/255,166/255,30/255,0.2)
clr6b<-rgb(230/255,171/255,2/255,0.2)
clr7b<-rgb(166/255,118/255,29/255,0.2)
clr8b<-rgb(233/255,125/255,157/255) #cox purple pale
clr9b<-rgb(144/255,210/255,169/255) #OH green pale

clr1c<-rgb(27/255,158/255,119/255,0.5)
clr2c<-rgb(217/255,95/255,2/255,0.5)
clr3c<-rgb(117/255,112/255,179/255,0.5)
clr4c<-rgb(231/255,41/255,138/255,0.5)
clr5c<-rgb(102/255,166/255,30/255,0.5)
clr6c<-rgb(230/255,171/255,2/255,0.5)
clr7c<-rgb(166/255,118/255,29/255,0.5)

##########
###COX####
##########

set.seed(4112024)
n<-1e6

library(survival)
library(pracma)
library(tikzDevice)

par(mfrow = c(1,1),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)


###Discrete + Gamma ####

OH_1f<-function(t,HR0,theta0){
  HR0*(1/theta0+t^3/60)/(1/theta0+t^3/(60/HR0))
}

q0 <- function(ts) {
  lambda_ = 1/60 #scale
  k = 3 #shape
  theta0 = 1
  
  S0_gamma = function(t){(1 + theta0*((t**k)*lambda_))**(-1/theta0)}
  
  q_ = sapply(ts, function(t){1-S0_gamma(HRR*t)})
  
  return(q_)
}

#pdf("fig1.pdf")
tikz("fig1",width=6,height=4)

var0=1
HR_ = 1/3
HRR = (1/3)**(1/3)
#plot(seq(0.1,10,0.1),sapply(seq(0.1,10,0.1),function(x){OH_1f(x,HR=HR_,theta0=1)}),type="l",lty=1,col=clr5c,ylim=c(0,3),ylab="Hazard ratio",xlab="t",axes=F)
plot(q0(seq(0.1,10,0.1)),rep(HRR, length(seq(0.1,10,0.1))),type="l",lty=1,col=clr9b,ylim=c(0.2,1),xlab="t",ylab="Estimands",axes=F)
abline(h=HR_,lwd=1,col=clr8b)
#abline(h=HR_**(1/3),lwd=1,col=clr5b)
abline(h=1,lwd=2,col=rgb(220/255,220/255,220/255,0.2))

count<-0
for(k in c(0, 0.5, 1, 2)){
  count<-count+1
  HR<-HR_

  var0<-1
  U0<-rgamma(n,shape=1/var0,scale=var0)

  NT<-runif(n,0,1)
  T0<-(-(60/U0)*log(NT))^(1/3)
  #T1<-(-(60/(U0*U1))*log(NT))^(1/3)
  T1<-(-(60/(U0*HR))*log(NT))^(1/3)

  pA<-0.5
  A<-rbinom(n,1,pA)

  T<-A*T1+(1-A)*T0

  if(k==0){
    TC<-rep(Inf,n)
    C<-rep(0,n)
  }
  else{
    TC<-rexp(n,1/(k*min(mean(T1),mean(T0))))

    C<-(TC<T)
  }
  COXn<-function(FT){
    a<-(mean(sapply(T[which(T<=FT & !C)],function(x){log(OH_1f(x,HR=HR_,theta0=var0))})))
    #b<-as.numeric(exp(coef(coxph(Surv(T, E) ~ A, data = data.frame(T=apply(cbind(T,TC,rep(FT,n)),1,function(x){min(x)}), E=as.numeric((T<=FT & !C)), A=A),init=a))))
    #return(c(exp(a),b))
    return(exp(a))
    #return(b)
  }

  theta_n<-function(FT){
    df = data.frame(T=apply(cbind(T,TC,rep(FT,n)),1,function(x){min(x)}), E=as.numeric((T<=FT & !C)), A=A)

    survobj = Surv(time = df$T, event = df$E)
    fit = survfit(survobj~A, data = df)

    fit_0 = fit[1]
    fit_a = fit[2]

    S0<-function(t){summary(fit_0, times=t, extend=TRUE)$surv}
    S1<-function(t){summary(fit_a, times=t, extend=TRUE)$surv}

    #max_ts = quantile(df$T,0.9)
    #max_ts = max(df$T)
    
    #eval_ts = seq(0.1, max_ts + 1e0, 1e0)
    eval_ts = seq(0.1, FT, 0.1)
    
    eval_S0s = S0(eval_ts)
    eval_S1s = S1(eval_ts)
    inverse_S0 = function(t){interp1(eval_S0s, eval_ts, t)}

    theta <-function(t, S1_t){
      return(inverse_S0(S1_t)/t)
    }

    #max_idx = tail(which(eval_S1s > tail(eval_S0s, 1)), 1)

    #theta_ = mapply(theta, eval_ts[1:max_idx], eval_S1s[1:max_idx])
    #target.index <- which(abs(eval_ts[1:max_idx] - FT) == min(abs(eval_ts[1:max_idx] - FT)))

    #We need a very large N if we want to start at 0.1
    res_ = sapply(seq(0.5,FT,0.1),function(s){inverse_S0(S1(s))/s})
    return(res_)
  }


  #test<-cbind(sapply(seq(2,10,0.5),function(x){COXn(x)}))
  #lines(seq(2,10,0.5),(test[2,]),type="l",col=clr3,lty=count)
  #lines(seq(2,10,0.5),(test[1,]),type="l",col=clr4,lty=count)

  #lines(seq(1,10,0.5),sapply(seq(1,10,0.5),function(x){COXn(x)}),col=clr4,lty=count)
  
  points(q0(seq(1,10,0.5)),sapply(seq(1,10,0.5),function(x){COXn(x)}),col=clr8,pch=14+count,cex=0.5)
  for(FT in seq(2,10,0.5)){
    Y<-theta_n(FT)
  lines(q0(seq(0.5,FT,0.1)),Y,type="l",col=clr9,lty=count)
  points(q0(FT),tail(Y,1),col=clr9,cex=0.5,pch=14+count)
  #sapply(seq(1,10,0.5),function(x){theta_n(x)})
  }
}

 axis(side = 1,
      labels = TRUE)

 axis(side = 2,
      labels = TRUE)
 box(which = "plot", bty = "l")

 title(xlab = "$F_{T^a}(t)$",
       ylab = "Estimand",
       outer = TRUE, line = 3)

 legend(0.6,0.98,c("-", "0.5 E[T1]","E[T1]","2  E[T1]"),lty=c(1,2,3,4),pch=c(15,16,17,18),seg.len=4,title="E[Tc]", cex=0.5)
 dev.off()
