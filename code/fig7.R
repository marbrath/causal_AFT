library(survival)
library(pracma)
library(tikzDevice)
library(rstpm2)

clr1<-rgb(27/255,158/255,119/255)
clr2<-rgb(217/255,95/255,2/255)
clr3<-rgb(117/255,112/255,179/255)
clr4<-rgb(231/255,41/255,138/255)
clr5<-rgb(147/255,112/255,219/255)
clr6<-rgb(210/255,105/255,30/255) #chocolate
clr7<-rgb(85/255,107/255,47/255) #dark olive green
clr8<-rgb(58/255, 40/255, 74/255) #purple

tikz("case_study_fig7.tex",width=6,height=5)
#pdf("case_study_fig7.pdf")

################################################################
### Functions to extend predict.aft() to present theta(t)+CI ###
################################################################
{grad1<-
  function (func, x, ..., log.transform = FALSE) 
  {
    ny <- length(func(x, ...))
    if (ny == 0L) 
      stop("Length of function equals 0")
    if (log.transform) {
      h <- .Machine$double.eps^(1/3)
      value <- (func(x * exp(h), ...) - func(x * exp(-h), ...))/2/h/x
    }
    else {
      h <- .Machine$double.eps^(1/3) * ifelse(abs(x) > 1, abs(x), 
                                              1)
      value <- (func(x + h, ...) - func(x - h, ...))/2/h
    }
    return(value)
  }

lpmatrix.lm<-function (object, newdata, na.action = na.pass) 
{
  tt <- terms(object)
  if (!inherits(object, "lm")) 
    warning("calling predict.lm(<fake-lm-object>) ...")
  if (missing(newdata) || is.null(newdata)) {
    X <- model.matrix(object)
  }
  else {
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, na.action = na.action, 
                     xlev = object$xlevels)
    if (!is.null(cl <- attr(Terms, "dataClasses"))) 
      .checkMFClasses(cl, m)
    X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
  }
  X
}

setMethod("predict", "aft", function (object, ...) 
{
  .local <- function (object, newdata = NULL, type = c("surv", 
                                                       "cumhaz", "hazard", "density", "hr", 
                                                       "sdiff", "hdiff", "loghazard", "link", 
                                                       "meansurv", "meansurvdiff", "odds", 
                                                       "or", "meanhaz", "af", "fail", 
                                                       "accfac", "theta", "gradh"), grid = FALSE, seqLength = 300, 
                      level = 0.95, se.fit = FALSE, link = NULL, exposed = incrVar(var), 
                      var = NULL, keep.attributes = TRUE, ...) 
  {
    type = match.arg(type)
    if (object@args$tvc.integrated) 
      predict.aft_integrated2(object, newdata, type, grid, 
                              seqLength, level, se.fit, link, exposed, var, 
                              keep.attributes, ...)
    else predict.aft_mixture2(object, newdata, type, grid, 
                              seqLength, level, se.fit, link, exposed, var, keep.attributes, 
                              ...)
  }
  .local(object, ...)
})

predict.aft_mixture2<-function (object, newdata = NULL, type = c("surv", "cumhaz", 
                                                                 "hazard", "density", "hr", "sdiff", 
                                                                 "hdiff", "loghazard", "link", "meansurv", 
                                                                 "meansurvdiff", "odds", "or", "meanhaz", 
                                                                 "af", "fail", "accfac", "theta", "gradh"), 
                                grid = FALSE, seqLength = 300, level = 0.95, se.fit = FALSE, 
                                link = NULL, exposed = incrVar(var), var = NULL, keep.attributes = TRUE, 
                                ...) 
{
  type <- match.arg(type)
  args <- object@args
  if (type %in% c("fail")) {
    out <- 1 - predict(object, newdata = newdata, type = "surv", 
                       grid, seqLength, se.fit, link, exposed, var, keep.attributes, 
                       ...)
    if (se.fit) {
      temp <- out$lower
      out$lower <- out$upper
      out$upper <- temp
    }
    return(out)
  }
  if (is.null(exposed) && is.null(var) & type %in% c("hr", 
                                                     "sdiff", "hdiff", "meansurvdiff", "or", 
                                                     "af", "accfac", "theta")) 
    stop("Either exposed or var required for type in (\"hr\",\"sdiff\",\"hdiff\",\"meansurvdiff\",\"or\",\"af\",\"accfac\")")
  if (is.null(newdata) && type %in% c("hr", "sdiff", 
                                      "hdiff", "meansurvdiff", "or", "af", 
                                      "accfac", "theta")) 
    stop("Prediction using type in ('hr','sdiff','hdiff','meansurvdiff','or','af','accfac') requires newdata to be specified.")
  calcX <- !is.null(newdata)
  time <- NULL
  if (is.null(newdata)) {
    X <- args$X
    XD <- args$XD
    y <- args$y
    time <- as.vector(y[, ncol(y) - 1])
    newdata <- as.data.frame(args$data)
  }
  loglpfunc <- function(x, newdata, ...) {
    newdata[[object@args$timeVar]] <- exp(x)
    lpmatrix.lm(object@args$lm.obj, newdata = newdata)
  }
  if (grid) {
    Y <- args$y
    event <- Y[, ncol(Y)] == 1
    time <- args$data[[args$timeVar]]
    eventTimes <- time[event]
    tt <- seq(min(eventTimes), max(eventTimes), length = seqLength)[-1]
    data.x <- data.frame(tt)
    names(data.x) <- args$timeVar
    newdata[[args$timeVar]] <- NULL
    newdata <- merge(newdata, data.x)
    calcX <- TRUE
  }
  if (calcX) {
    X <- lpmatrix.lm(args$lm.obj, newdata)[, args$X.index, 
                                           drop = FALSE]
    XD <- grad1(loglpfunc, log(newdata[[object@args$timeVar]]), 
                log.transform = FALSE, newdata = newdata)
    XD <- matrix(XD, nrow = nrow(X))[, args$X.index, drop = FALSE]
    Xc <- lpmatrix.lm(args$glm.cure.obj, newdata)
    time <- eval(args$timeExpr, newdata)
  }
  if (type %in% c("hr", "sdiff", "hdiff", 
                  "meansurvdiff", "or", "af", "accfac", "theta")) {
    newdata2 <- exposed(newdata)
    X2 <- lpmatrix.lm(args$lm.obj, newdata2)[, args$X.index, 
                                             drop = FALSE]
    XD2 <- grad1(loglpfunc, log(newdata2[[object@args$timeVar]]), 
                 log.transform = FALSE, newdata = newdata2)
    XD2 <- matrix(XD2, nrow = nrow(X))[, args$X.index, drop = FALSE]
    time2 <- eval(args$timeExpr, newdata2)
    Xc2 = model.matrix(args$glm.cure.obj, newdata2)
  }
  if (type %in% c("grad")) {
    return(predict.aft.ext(object, type = type, time = time, 
                           X = X, XD = XD))
  }
  local <- function(object, newdata = NULL, type = "surv", 
                    exposed, ...) {
    args <- object@args
    betafull <- coef(object)
    beta <- betafull[1:ncol(args$X)]
    betas <- betafull[-(1:ncol(args$X))]
    tt <- args$terms
    eta <- as.vector(X %*% beta)
    #eta = -log(theta(t)) theta(t) = exp(-eta)
    logtstar <- log(time) - eta
    etas <- as.vector(predict(args$design, logtstar) %*% 
                        betas)
    H <- exp(etas)
    S <- exp(-H)
    if (type == "cumhaz") 
      return(H)
    if (type == "surv") 
      return(S)
    if (type == "fail") 
      return(1 - S)
    if (type == "odds") 
      return((1 - S)/S)
    if (type == "meansurv") 
      return(tapply(S, newdata[[object@args$timeVar]], 
                    mean))
    etaDs <- as.vector(predict(args$designD, logtstar) %*% 
                         betas)
    etaD <- as.vector(XD %*% beta)
    h <- H * etaDs * (1/time - etaD)
    Sigma = vcov(object)
    if (type == "link") 
      return(eta)
    if (type == "density") 
      return(S * h)
    if (type == "hazard") 
      return(h)
    if (type == "loghazard") 
      return(log(h))
    if (type == "meanhaz") 
      return(tapply(S * h, newdata[[object@args$timeVar]], 
                    sum)/tapply(S, newdata[[object@args$timeVar]], 
                                sum))
    eta2 <- as.vector(X2 %*% beta)
    logtstar2 <- log(time2) - eta2
    etas2 <- as.vector(predict(args$design, logtstar2) %*% 
                         betas)
    H2 <- exp(etas2)
    S2 <- exp(-H2)
    if (type == "sdiff") 
      return(S2 - S)
    if (type == "or") 
      return((1 - S2)/S2/((1 - S)/S))
    if (type == "meansurvdiff") 
      return(tapply(S2, newdata[[object@args$timeVar]], 
                    mean) - tapply(S, newdata[[object@args$timeVar]], 
                                   mean))
    etaDs2 <- as.vector(predict(args$designD, logtstar2) %*% 
                          betas)
    etaD2 <- as.vector(XD2 %*% beta)
    h2 <- H2 * etaDs2 * (1/time2 - etaD2)
    if (type == "hdiff") 
      return(h2 - h)
    if (type == "hr") 
      return(h2/h)
    if (type == "af") {
      meanS <- tapply(S, newdata[[object@args$timeVar]], 
                      mean)
      meanS2 <- tapply(S2, newdata[[object@args$timeVar]], 
                       mean)
      return((meanS2 - meanS)/(1 - meanS))
    }
    if (type == "accfac") {
      accfac <- -eta + log(1 - etaD)
      accfac2 <- -eta2 + log(1 - etaD2)
      return(exp(accfac2 - accfac))
    }
    if (type == "theta") {
      theta <- -eta 
      theta2 <- -eta2 
      return(exp(theta2 - theta))
      #return(exp(-theta))
      
    }
  }
  local2 <- function(object, newdata = NULL, type = "surv", 
                     exposed, ...) {
    args <- object@args
    betafull <- coef(object)
    beta <- betafull[1:ncol(args$X)]
    betac <- betafull[(ncol(args$X) + 1):(ncol(args$X) + 
                                            ncol(args$Xc))]
    betas <- betafull[-(1:(ncol(args$X) + ncol(args$Xc)))]
    tt <- args$terms
    eta <- as.vector(X %*% beta)
    etac <- as.vector(Xc %*% betac)
    cure_frac <- exp(etac)/(1 + exp(etac))
    logtstar <- log(time) - eta
    etas <- as.vector(predict(args$design, logtstar) %*% 
                        betas)
    Hu <- exp(etas)
    Su <- exp(-Hu)
    S <- cure_frac + (1 - cure_frac) * Su
    if (type == "cumhaz") 
      return(-log(cure_frac + (1 - cure_frac) * exp(-Hu)))
    if (type == "surv") 
      return(S)
    if (type == "fail") 
      return(1 - S)
    if (type == "odds") 
      return((1 - S)/S)
    if (type == "meansurv") 
      return(tapply(S, newdata[[object@args$timeVar]], 
                    mean))
    etaDs <- as.vector(predict(args$designD, logtstar) %*% 
                         betas)
    etaD <- as.vector(XD %*% beta)
    hu <- Hu * etaDs * (1/time - etaD)
    h <- (1 - cure_frac) * exp(-Hu) * hu/(cure_frac + (1 - 
                                                         cure_frac) * exp(-Hu))
    Sigma = vcov(object)
    if (type == "link") 
      return(eta)
    if (type == "density") 
      return(S * h)
    if (type == "hazard") 
      return(h)
    if (type == "loghazard") 
      return(log(h))
    if (type == "meanhaz") 
      return(tapply(S * h, newdata[[object@args$timeVar]], 
                    sum)/tapply(S, newdata[[object@args$timeVar]], 
                                sum))
    eta2 <- as.vector(X2 %*% beta)
    logtstar2 <- log(time2) - eta2
    etas2 <- as.vector(predict(args$design, logtstar2) %*% 
                         betas)
    etac2 <- as.vector(Xc2 %*% betac)
    cure_frac2 <- exp(etac2)/(1 + exp(etac2))
    Hu2 <- exp(etas2)
    Su2 <- exp(-Hu2)
    S2 <- cure_frac2 + (1 - cure_frac2) * Su2
    if (type == "sdiff") 
      return(S2 - S)
    if (type == "or") 
      return((1 - S2)/S2/((1 - S)/S))
    if (type == "meansurvdiff") 
      return(tapply(S2, newdata[[object@args$timeVar]], 
                    mean) - tapply(S, newdata[[object@args$timeVar]], 
                                   mean))
    etaDs2 <- as.vector(predict(args$designD, logtstar2) %*% 
                          betas)
    etaD2 <- as.vector(XD2 %*% beta)
    hu2 <- Hu2 * etaDs2 * (1/time2 - etaD2)
    h2 <- (1 - cure_frac2) * exp(-Hu2) * hu2/(cure_frac2 + 
                                                (1 - cure_frac2) * exp(-Hu2))
    if (type == "hdiff") 
      return(h2 - h)
    if (type == "hr") 
      return(h2/h)
    if (type == "af") {
      meanS <- tapply(S, newdata[[object@args$timeVar]], 
                      mean)
      meanS2 <- tapply(S2, newdata[[object@args$timeVar]], 
                       mean)
      return((meanS2 - meanS)/(1 - meanS))
    }
    if (type == "accfac") {
      accfac <- -eta + log(1 - etaD)
      accfac2 <- -eta2 + log(1 - etaD2)
      return(exp(accfac2 - accfac))
    }
    if (type == "theta") {
      theta <- -eta 
      theta2 <- -eta2 
      return(exp(theta2 - theta))
      #return(exp(-theta))
      
    }
    
  }
  local <- if (args$mixture) 
    local2
  else local
  out <- if (!se.fit) {
    local(object, newdata, type = type, exposed = exposed, 
          ...)
  }
  else {
    if (is.null(link)) 
      link <- switch(type, surv = "cloglog", cumhaz = "log", 
                     hazard = "log", hr = "log", sdiff = "I", 
                     hdiff = "I", loghazard = "I", link = "I", 
                     odds = "log", or = "log", meansurv = "I", 
                     meanhaz = "I", af = "I", accfac = "log", theta = "log")
    invlinkf <- switch(link, I = I, log = exp, cloglog = cexpexp, 
                       logit = expit)
    pred <- predictnl(object, local, link = link, newdata = newdata, 
                      type = type, gd = NULL, exposed = exposed, ...)
    ci <- confint.predictnl(pred, level = level)
    out <- data.frame(Estimate = pred$fit, lower = ci[, 1], 
                      upper = ci[, 2])
    if (link == "cloglog") 
      out <- data.frame(Estimate = out$Estimate, lower = out$upper, 
                        upper = out$lower)
    invlinkf(out)
  }
  if (keep.attributes) 
    attr(out, "newdata") <- newdata
  return(out)
}
}

whole_df = read.csv("./GASTRIC advanced data.csv")
study_ids<-unique(whole_df$id_study)
i<-14
df = subset(whole_df, id_study==study_ids[i])
df$arm = df$arm-1

knots<-c(quantile(df$pfs[df$arm==0],0.1),quantile(df$pfs[df$arm==0],0.2))

# is the time-invariant model? Why ~1 and not ~a?
aft14<-aft(Surv(pfs,PFS_status_0cens_1ProgDeath==1)~1,data=df,df=3,smooth.formula=~arm:nsx(log(pfs),knots=log(knots)))

#Realize S_0(t) = S_1(t*theta(t)) scale
theta<-predict(aft14,newdata=data.frame(arm=1, pfs=seq(1,600,1)),type="theta",exposed=function(data) transform(data,arm=0))

S0<-function(t){predict(aft14,newdata=data.frame(arm=1, pfs=t),type="surv")[1]} ####???????????

S1spline<-function(t){predict(aft14,newdata=data.frame(arm=0, pfs=t),type="surv")[1]} #########?????
S1<-function(t,mu1,mu2){0.5*predict(aft14,newdata=data.frame(arm=1, pfs=mu1*t),type="surv")[1]+0.5*predict(aft14,newdata=data.frame(arm=1, pfs=t*mu2),type="surv")[1]}

mu1<-0.9
mu2<-0.45

plot(survfit(Surv(pfs,PFS_status_0cens_1ProgDeath==1)~arm,data=df),col=c(clr2,clr8), ylab="S(t)", xlab='t', lty=2) #a=0=orange, a=1=green
#plot(aft14,newdata=data.frame(arm=1), add=TRUE, line.col="green",ci=F)
#plot(aft14,newdata=data.frame(arm=0), add=TRUE, line.col="blue",ci=F)
plot(aft14,newdata=data.frame(arm=0), add=TRUE, line.col=clr2,ci=F)
plot(aft14,newdata=data.frame(arm=1), add=TRUE, line.col=clr8,ci=F)
legend(420,0.9,legend=c(1,0),lty=c(1,1), col=c(clr2,clr8),title="a")
dev.off()
