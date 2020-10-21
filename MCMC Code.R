##### Load Packages
library(invgamma)
library(truncnorm)
library(mvtnorm)
library(deSolve)
library(Rmpfr)
library(splines2)
library(MASS)
library(HDInterval)

############################################################################
#####  Specify the path to the folder to which output will be written  #####
############################################################################
out.path = "out/path/"

##############################
#####  Read in the Data  #####
##############################

###See Example Data for the Appropriate Formatting
df = read.csv(file.choose())
pop.df = read.csv(file.choose())

##############################
#####  Input Parameters  #####
##############################

### The names of the counties to be included in the Model
county.list = c('County1', 'County2')

###First day in the model (usually the first day for which data is available)
start.day = as.Date("3/15/2020",tryFormat = "%m/%d/%Y")

###Last day in the model (usually the last day on which any data is observed)
end.day = as.Date("6/8/2020",tryFormat = "%m/%d/%Y")

###Date through which predictions are desired
end.pred.day = as.Date("6/22/2020",tryFormat = "%m/%d/%Y")

###First Date that reported COVID-19 Incidence data is available
start.county.incidence.day = as.Date("3/15/2020",tryFormat = "%m/%d/%Y")

###Last Date that reported COVID-19 Incidence data is available
end.county.incidence.day = as.Date("6/8/2020",tryFormat = "%m/%d/%Y")

###First Date that COVID-19 inpatient census data is available
start.census.day = as.Date("3/15/2020",tryFormat = "%m/%d/%Y")

###Last Date that COVID-19 inpatient census data is available
end.census.day = as.Date("6/8/2020",tryFormat = "%m/%d/%Y")

###First Date that COVID-19 daily hospital admissions data is available
start.admits.day = as.Date("3/15/2020",tryFormat = "%m/%d/%Y")

###Last Date that COVID-19 daily hospital admissions data is available
end.admits.day = as.Date("6/8/2020",tryFormat = "%m/%d/%Y")

###First day for which social distancing data is available
start.soc.dist.day = as.Date("3/1/2020",tryFormat = "%m/%d/%Y")

###Last day for which social distancing data is available
end.soc.dist.day = as.Date("6/8/2020",tryFormat = "%m/%d/%Y")

###Assumed delay between changes in the amount of social distancing and changes in the reported COVID-19 incidence
soc.dist.lag = 14

###Degrees of Freedom for Spline Basis for the Transmission Rate
beta.df = 4
hr.df = 3

######################################################
#####  Variables Defined using input parameters  #####
######################################################

###Time Variables
###Number of days of observed data
tm = length(seq.Date(from = start.day,to = end.day,by = 1))
###Number of days of observed data plus forecasted days
tm.pred = length(seq.Date(from = start.day, to = end.pred.day,by = 1))
t <- 1:tm
t.pred <- 1:tm.pred

###Data Index Variables
###Should match the formatting of the date feild in the data file
date.format = "%Y-%m-%d"
census.start.ind = which(as.Date(df$Date,format = date.format) == start.census.day)
census.end.ind = which(as.Date(df$Date,format = date.format) == end.census.day)
admits.start.ind = which(as.Date(df$Date,format = date.format) == start.admits.day)
admits.end.ind = which(as.Date(df$Date,format = date.format) == end.admits.day)
county.incidence.start.ind = which(as.Date(df$Date,format = date.format) == start.county.incidence.day)
county.incidence.end.ind = which(as.Date(df$Date,format = date.format) == end.county.incidence.day)
#dhec.pred.end.ind = which(as.Date(df$Date,format = date.format) == max(as.Date(dhec.df$Date)))

###The number of counties
n.mod.cty = length(county.list)
###The positions of those counties in df
county.ind = c()
for(i in 1:n.mod.cty){
  county.ind = c(county.ind, which(names(df) == paste0("New_Reported_Infections_",county.list[i])))
}

#############################
#####  Format the Data  #####
#############################

###Format the Population Data
N.mat = matrix(NA,tm.pred,n.mod.cty)
for(i in 1:n.mod.cty){
  ind = which(pop.df$County.Name == county.list[i])
  N.mat[,i] = pop.df$Population[ind]
}

###Format the DHEC and Prisma Data
Y <- df[county.incidence.start.ind:county.incidence.end.ind,county.ind]
if(start.county.incidence.day > start.day){
  for(i in 1:length(seq.Date(from = as.Date(start.day,format = date.format),to = as.Date((start.county.incidence.day -1),format = date.format),1))){
    Y = rbind(NA,Y)
  }
}
if(end.county.incidence.day < end.day){
  for(i in 1:length(seq.Date(from = as.Date(end.county.incidence.day,format = date.format)+1,to = as.Date(end.day,format = date.format),1))){
    Y = rbind(Y,NA)
  }
}
Y.ind = !is.na(Y)

###Format the non-ventilated census data
Y.c <- df$Non_Ventitlated_Census[census.start.ind:census.end.ind]
if(start.day <start.census.day){
  Y.c = c(rep(NA,length(seq.Date(from = as.Date(start.day,format = date.format),to = (as.Date(start.census.day,format = date.format)-1),by = 1))),Y.c)
}
if(end.census.day < end.day){
  Y.c = c(Y.c,rep(NA,length(seq.Date(from = as.Date(end.census.day,format = date.format)+1,to = as.Date(end.day,format = date.format),by = 1))))
}
Y.c.ind = !is.na(Y.c)

###Format the ventilated census data
Y.v <- df$Ventilated_Census[census.start.ind:census.end.ind]
if(start.day <start.census.day){
  Y.v = c(rep(NA,length(seq.Date(from = as.Date(start.day,format = date.format),to = as.Date(start.census.day,format = date.format)-1,by = 1))),Y.v)
}
if(end.census.day < end.day){
  Y.v = c(Y.v,rep(NA,length(seq.Date(from = as.Date(end.census.day,format = date.format)+1,to = as.Date(end.day,format = date.format),by = 1))))
}
Y.v.ind = !is.na(Y.v)

###Format the admissions data
Y.h <- df$New_Admissions[admits.start.ind:admits.end.ind]
if(start.day <start.admits.day){
  Y.h = c(rep(NA,length(seq.Date(from = as.Date(start.day,format = date.format),to = as.Date(start.admits.day,format = date.format)-1,by = 1))),Y.h)
}
if(end.admits.day < end.day){
  Y.h = c(Y.h,rep(NA,length(seq.Date(from = as.Date(end.admits.day,format = date.format)+1,to = as.Date(end.day,format = date.format),by = 1))))
}
Y.h.ind = !is.na(Y.h)

###Examine plots to ensure the data is formatted correctly
for(i in 1:n.mod.cty){
  plot(seq.Date(from = as.Date(start.day,format = date.format),to = as.Date(end.day,format = date.format),by = 1),Y[,i],type = 'o')
}
plot(seq.Date(from = as.Date(start.day,format = date.format),to = as.Date(end.day,format = date.format),by = 1),Y.c,type = 'o')
plot(seq.Date(from = as.Date(start.day,format = date.format),to = as.Date(end.day,format = date.format),by = 1),Y.h,type = 'o')
plot(seq.Date(from = as.Date(start.day,format = date.format),to = as.Date(end.day,format = date.format),by = 1),Y.v,type = 'o')


###Format the Social Distancing Data
soc.dist.date.seq = seq.Date(from = as.Date((start.day-soc.dist.lag),format = date.format),to = end.soc.dist.day,by = 1)
###If we are missing social distancing data from early days in the model, use the first observed value
soc.dist.date.seq[soc.dist.date.seq < start.soc.dist.day] = start.soc.dist.day
###Forecast future social distancing data using last value carried forward
while(length(soc.dist.date.seq)< tm.pred){
  soc.dist.date.seq = c(soc.dist.date.seq,soc.dist.date.seq[length(soc.dist.date.seq)])
}
soc.dist.mat = matrix(NA,tm.pred,n.mod.cty)
for(c in 1:n.mod.cty){
  col.ind = which(names(df) == paste0('Social_Distancing_',county.list[c]))
  for(j in 1:dim(soc.dist.mat)[1]){
    row.ind = which(as.Date(df$Date,format = date.format) == soc.dist.date.seq[j])
    soc.dist.mat[j,c] = df[row.ind,col.ind]
  }
}

###Examine Plots to ensure data is formatted correctly
for(c in 1:n.mod.cty){
  plot(soc.dist.date.seq[1:dim(soc.dist.mat)[1]],soc.dist.mat[,c],type = 'o')
}

##########################################
#####  Create Spline Basis Matrices  #####
##########################################

#Make Matrix for Transmission Parameter
B.sd = bSpline(t,df = beta.df, Boundary.knots = c(-1,(tm+1)))
B.sd = cbind(1,B.sd)
numb.b = dim(B.sd)[2] + 1
b.sdi = numb.b
###Predict Out in Time
for(i in (tm + 1):(tm.pred)){
  B.sd = rbind(B.sd,B.sd[tm,])
}

###Make Matrix for Hospitalization Rate Parameter
B.hr = bSpline(t,df = hr.df, Boundary.knots = c(-1,(tm+1)))
B.hr = cbind(1,B.hr)
numb.hr = dim(B.hr)[2]# + 1
###Predict Out in Time
for(i in (tm + 1):(tm.pred)){
  B.hr = rbind(B.hr,B.hr[tm,])
}

#Make Matrix for Vent Rate Parameter
B.vr = matrix(1, tm.pred, 1)
numb.vr = dim(B.vr)[2]

###Make a matrix to format the RE inside the SIR Function
re.mat.one = matrix(1,n.mod.cty,tm.pred)

###############################################
#####  Sampling and Likelihood Functions  #####
###############################################

###DE Function
SIR <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    S = state[1:n.mod.cty]
    I = state[(n.mod.cty + 1):(2*n.mod.cty)]
    R = state[(2*n.mod.cty +1):(3*n.mod.cty)]
    H = state[(3*n.mod.cty + 1):(4*n.mod.cty)]
    V = state[(4*n.mod.cty + 1):(5*n.mod.cty)]
    dS <- -exp(as.vector(B.sd[time,]%*%b[-c(b.sdi)]) + b[b.sdi]*soc.dist.mat[time,] + t(re*re.mat.one[,time]))*I * S/N.mat[time,]
    dI <- exp(as.vector(B.sd[time,]%*%b[-c(b.sdi)]) + b[b.sdi]*soc.dist.mat[time,] + t(re*re.mat.one[,time]))*I * S/N.mat[time,] - gam* I - as.vector(exp(B.hr[time,]%*%hr)/(1 + exp(B.hr[time,]%*%hr)))*I
    dR <- gam * I + gam.hr*H
    dH <- as.vector(exp(B.hr[time,]%*%hr)/(1 + exp(B.hr[time,]%*%hr)))*I - gam.hr*H - as.vector(exp(B.vr[time,]%*%vr)/(1 + exp(B.vr[time,]%*%vr)))*H + gam.vr*V
    dV <- as.vector(exp(B.vr[time,]%*%vr)/(1 + exp(B.vr[time,]%*%vr)))*H - gam.vr*V
    list(c(dS, dI, dR, dH, dV))
  })
}

###The log likelihood function
llh<-function(lambda, lambda.h, lambda.c,lambda.v, nb,nb.h, nb.c,nb.v){
  res = sum(dnbinom(Y[Y.ind],mu =lambda[Y.ind],size = nb,log = TRUE))
  res = res + sum(dnbinom(Y.h[Y.h.ind],mu = lambda.h[Y.h.ind],size = nb.h,log = TRUE))
  res = res + sum(dnbinom(Y.c[Y.c.ind],mu = lambda.c[Y.c.ind],size = nb.c,log = TRUE))
  res = res + sum(dnbinom(Y.v[Y.v.ind],mu = lambda.v[Y.v.ind],size = nb.v,log = TRUE))
  return(res)
}

###Function to calculate the number of new infections
lambda.calc<-function(fit,par,t){
  res = (exp(as.vector(B.sd[t,]%*%par$b[-c(b.sdi)]) + par$b[b.sdi]*soc.dist.mat[t,] + t(par$re*re.mat.one[,t]))* as.matrix(fit[t,I.ind])*as.matrix(fit[t,S.ind])/N.mat[t,])
  return(res)
}

###Function to calculate the number of new admissions
lambda.h.calc<-function(fit,par,t){
  res = as.vector(exp(B.hr[t,]%*%par$hr)/(1 + exp(B.hr[t,]%*%par$hr)))*apply(fit[t,I.ind],1,sum)
  return(res)
}

###Function to calculate the non-ventilated census
lambda.c.calc<-function(fit,par,t){
  res = apply(fit[t,H.ind],1,sum)
  return(res)
}

###Function to calculate the ventilated census
lambda.v.calc<-function(fit,par,t){
  res = apply(fit[t,V.ind],1,sum)
  return(res)
}

###Sample the global transmssion parameters
b.samp<-function(par.g,b.i){
  b.i.prop = rnorm(1,unlist(par.g)[b.i],sd.b[b.i])
  par.prop = par.g
  par.prop[[1]][b.i] = b.i.prop
  fit.prop <- data.frame(ode(y = init.g, times = t, func = SIR, parms = par.prop,method = 'euler'))
  lambda.prop = lambda.calc(fit.prop,par.prop,t)
  lambda.hr.prop = lambda.h.calc(fit.prop,par.prop,t)
  lambda.c.prop = lambda.c.calc(fit.prop,par.prop,t)
  lambda.v.prop = lambda.v.calc(fit.prop,par.prop,t)
  llh.p = llh(lambda.prop,lambda.hr.prop,lambda.c.prop,lambda.v.prop,nb.g,nb.h.g,nb.c.g,nb.v.g)
  prior.g = log(dtruncnorm(x = unlist(par.g)[b.i],mean = m.b, sd = ps.b, a = l.b,b = u.b))
  prior.p = log(dtruncnorm(x = b.i.prop,mean = m.b, sd = ps.b, a = l.b,b = u.b))
  if(prior.p != -Inf){
    e.a =exp( mpfr(llh.p + prior.p,100) - mpfr(llh.g + prior.g,100))
    a = min(1,as.numeric(e.a))
    r = rbinom(1,1,a)
  }else{
    r = 0
  }
  res = unlist(par.prop)[b.i]*r + unlist(par.g)[b.i]*(1-r)
  llh.res = llh.p*r + llh.g*(1-r)
  if(is.na(r)){
    r = 0
    res = unlist(par.g)[b.i]
    llh.res = llh.g
  }
  return(list(res,llh.res,r))
}

###Sample the County level random effects
re.samp<-function(par.g,r.i){
  re.g = par.g$re
  b.i.prop = rnorm(1,re.g[r.i],sd.re[r.i])
  par.prop = par.g
  par.prop$re[r.i] = b.i.prop
  fit.prop <- data.frame(ode(y = init.g, times = t, func = SIR, parms = par.prop,method = 'euler'))
  lambda.prop = lambda.calc(fit.prop,par.prop,t)
  lambda.hr.prop = lambda.h.calc(fit.prop,par.prop,t)
  lambda.c.prop = lambda.c.calc(fit.prop,par.prop,t)
  lambda.v.prop = lambda.v.calc(fit.prop,par.prop,t)
  llh.p = llh(lambda.prop,lambda.hr.prop,lambda.c.prop,lambda.v.prop,nb.g,nb.h.g,nb.c.g,nb.v.g)
  prior.g = log(dtruncnorm(x = re.g[r.i],mean = m.b, sd = ps.b, a = l.b,b = u.b))
  prior.p = log(dtruncnorm(x = b.i.prop,mean = m.b, sd = ps.b, a = l.b,b = u.b))
  if(prior.p != -Inf){
    e.a =exp( mpfr(llh.p + prior.p,100) - mpfr(llh.g + prior.g,100))
    a = min(1,as.numeric(e.a))
    r = rbinom(1,1,a)
  }else{
    r = 0
  }
  res = b.i.prop*r + re.g[r.i]*(1-r)
  llh.res = llh.p*r + llh.g*(1-r)
  if(is.na(r)){
    r = 0
    res = re.g[r.i]*(1-r)
    llh.res = llh.g
  }
  return(list(res,llh.res,r))
}

###Sample the recovery rate
gam.samp<-function(par.g){
  log.b.g = log(par.g$gam/(1 - par.g$gam))
  log.b.p = rnorm(1,log.b.g,sd.gam)
  b.i.prop = exp(log.b.p)/(1 + exp(log.b.p))
  par.prop = par.g
  par.prop$gam = b.i.prop
  fit.prop <- data.frame(ode(y = init.g, times = t, func = SIR, parms = par.prop,method = 'euler'))
  lambda.prop = lambda.calc(fit.prop,par.prop,t)
  lambda.hr.prop = lambda.h.calc(fit.prop,par.prop,t)
  lambda.c.prop = lambda.c.calc(fit.prop,par.prop,t)
  lambda.v.prop = lambda.v.calc(fit.prop,par.prop,t)
  llh.p = llh(lambda.prop,lambda.hr.prop,lambda.c.prop,lambda.v.prop,nb.g,nb.h.g,nb.c.g,nb.v.g)
  prior.g = 0#log(dtruncnorm(x = par.g$gam,mean = m.gam, sd = ps.gam, a = l.gam,b = u.gam))
  prior.p = 0#log(dtruncnorm(x = b.i.prop,mean = m.gam, sd = ps.gam, a = l.gam,b = u.gam))
  if(prior.p != -Inf){
    e.a =exp( mpfr(llh.p + prior.p + log(1/(par.g$gam*(1-par.g$gam))),100) - mpfr(llh.g + prior.g + log(1/(b.i.prop*(1-b.i.prop))),100))
    a = min(1,as.numeric(e.a))
    r = rbinom(1,1,a)
  }else{
    r = 0
  }
  res = par.prop$gam*r + par.g$gam*(1-r)
  llh.res = llh.p*r + llh.g*(1-r)
  if(is.na(r)){
    r = 0
    res =  par.g$gam
    llh.res = llh.g
  }
  return(list(res,llh.res,r))
}

###Sample the discharge rate
gam.hr.samp<-function(par.g){
  log.b.g = log(par.g$gam.hr/(1 - par.g$gam.hr))
  log.b.p = rnorm(1,log.b.g,sd.gam.hr)
  b.i.prop = exp(log.b.p)/(1 + exp(log.b.p))
  par.prop = par.g
  par.prop$gam.hr = b.i.prop
  fit.prop <- data.frame(ode(y = init.g, times = t, func = SIR, parms = par.prop,method = 'euler'))
  lambda.prop = lambda.calc(fit.prop,par.prop,t)
  lambda.hr.prop = lambda.h.calc(fit.prop,par.prop,t)
  lambda.c.prop = lambda.c.calc(fit.prop,par.prop,t)
  lambda.v.prop = lambda.v.calc(fit.prop,par.prop,t)
  llh.p = llh(lambda.prop,lambda.hr.prop,lambda.c.prop,lambda.v.prop,nb.g,nb.h.g,nb.c.g,nb.v.g)
  prior.g = 0#log(dtruncnorm(x = par.g$gam.hr,mean = m.gam, sd = ps.gam, a = l.gam,b = u.gam))
  prior.p = 0#log(dtruncnorm(x = b.i.prop,mean = m.gam, sd = ps.gam, a = l.gam,b = u.gam))
  if(prior.p != -Inf){
    e.a =exp( mpfr(llh.p + prior.p + log(1/(par.g$gam.hr*(1-par.g$gam.hr))),100) - mpfr(llh.g + prior.g + log(1/(b.i.prop*(1-b.i.prop))),100))
    a = min(1,as.numeric(e.a))
    r = rbinom(1,1,a)
  }else{
    r = 0
  }
  res = par.prop$gam.hr*r + par.g$gam.hr*(1-r)
  llh.res = llh.p*r + llh.g*(1-r)
  if(is.na(r)){
    r = 0
    res =  par.g$gam.hr
    llh.res = llh.g
  }

  return(list(res,llh.res,r))
}

###Sample the vent discharge rate
gam.vr.samp<-function(par.g){
  log.b.g = log(par.g$gam.vr/(1 - par.g$gam.vr))
  log.b.p = rnorm(1,log.b.g,sd.gam.vr)
  b.i.prop = exp(log.b.p)/(1 + exp(log.b.p))
  par.prop = par.g
  par.prop$gam.vr = b.i.prop
  fit.prop <- data.frame(ode(y = init.g, times = t, func = SIR, parms = par.prop,method = 'euler'))
  lambda.prop = lambda.calc(fit.prop,par.prop,t)
  lambda.hr.prop = lambda.h.calc(fit.prop,par.prop,t)
  lambda.c.prop = lambda.c.calc(fit.prop,par.prop,t)
  lambda.v.prop = lambda.v.calc(fit.prop,par.prop,t)
  llh.p = llh(lambda.prop,lambda.hr.prop,lambda.c.prop,lambda.v.prop,nb.g,nb.h.g,nb.c.g,nb.v.g)
  prior.g = 0#log(dtruncnorm(x = par.g$gam.vr,mean = m.gam, sd = ps.gam, a = l.gam,b = u.gam))
  prior.p = 0#log(dtruncnorm(x = b.i.prop,mean = m.gam, sd = ps.gam, a = l.gam,b = u.gam))
  if(prior.p != -Inf){
    e.a =exp( mpfr(llh.p + prior.p + log(1/(par.g$gam.vr*(1-par.g$gam.vr))),100) - mpfr(llh.g + prior.g + log(1/(b.i.prop*(1-b.i.prop))),100))
    a = min(1,as.numeric(e.a))
    r = rbinom(1,1,a)
  }else{
    r = 0
  }
  res = par.prop$gam.vr*r + par.g$gam.vr*(1-r)
  llh.res = llh.p*r + llh.g*(1-r)
  if(is.na(r)){
    r = 0
    res =  par.g$gam.vr
    llh.res = llh.g
  }

  return(list(res,llh.res,r))
}

###Sample the number of initially infected individuals
init.inf.samp<-function(init.g,ii){
  log.init.inf.g = log(init.g[(I.ind[ii]-1)])
  log.init.inf.p = rnorm(1,log.init.inf.g,sd.inf[[ii]])
  init.inf.prop = exp(log.init.inf.p)
  init.prop = init.g
  init.prop[(I.ind[ii]-1)] = init.inf.prop
  init.prop[(S.ind[ii]-1)] = init.prop[(S.ind[ii]-1)] + init.g[(I.ind[ii]-1)] - init.inf.prop
  fit.prop <- data.frame(ode(y = init.prop, times = t, func = SIR, parms = parameters.g,method = 'euler'))
  lambda.prop = lambda.calc(fit.prop,parameters.g,t)
  lambda.hr.prop = lambda.h.calc(fit.prop,parameters.g,t)
  lambda.c.prop = lambda.c.calc(fit.prop,parameters.g,t)
  lambda.v.prop = lambda.v.calc(fit.prop,parameters.g,t)
  llh.p = llh(lambda.prop,lambda.hr.prop,lambda.c.prop,lambda.v.prop,nb.g,nb.h.g,nb.c.g,nb.v.g)
  prior.g = log(dtruncnorm(x = init.g[(I.ind[ii]-1)],mean = m.init.inf, sd = ps.init.inf, a = l.init.inf,b = u.init.inf[ii]))
  prior.p = log(dtruncnorm(x = init.prop[(I.ind[ii]-1)],mean = m.init.inf, sd = ps.init.inf, a = l.init.inf,b = u.init.inf[ii]))
  if(prior.p != -Inf){
    e.a =exp( mpfr(llh.p + prior.p + log(1/init.g[(I.ind[ii]-1)]),100) - mpfr(llh.g + prior.g + log(1/init.inf.prop),100))
    a = min(1,as.numeric(e.a))
    r = rbinom(1,1,a)
  }else{
    r = 0
  }
  res = init.prop*r + init.g*(1-r)
  llh.res = llh.p*r + llh.g*(1-r)
  if(is.na(r)){
    r = 0
    res =  init.g
    llh.res = llh.g
  }
  return(list(res,llh.res,r))
}

### Sample the hosptialization rate
hr.samp<-function(par.g,hr.i){
  hr = par.g$hr[hr.i]
  hr.i.prop = rnorm(1,hr,sd.hr[hr.i])
  par.prop = par.g
  par.prop$hr[hr.i] = hr.i.prop
  fit.prop <- data.frame(ode(y = init.g, times = t, func = SIR, parms = par.prop,method = 'euler'))
  lambda.prop = lambda.calc(fit.prop,par.prop,t)
  lambda.hr.prop = lambda.h.calc(fit.prop,par.prop,t)
  lambda.c.prop = lambda.c.calc(fit.prop,par.prop,t)
  lambda.v.prop = lambda.v.calc(fit.prop,par.prop,t)
  llh.p = llh(lambda.prop,lambda.hr.prop,lambda.c.prop,lambda.v.prop,nb.g,nb.h.g,nb.c.g,nb.v.g)
  prior.hr.g = log(dtruncnorm(x = par.g$hr[hr.i],mean = m.b, sd = ps.b, a = l.b,b = u.b))
  prior.hr.p = log(dtruncnorm(x = hr.i.prop,mean = m.b, sd = ps.b, a = l.b,b = u.b))
  e.a =exp( mpfr(llh.p + prior.hr.p,100) - mpfr(llh.g + prior.hr.g ,100))
  a = min(1,as.numeric(e.a))
  r = rbinom(1,1,a)
  res = hr.i.prop*r + hr*(1-r)
  llh.res = llh.p*r + llh.g*(1-r)
  if(is.na(r)){
    r = 0
    res =  hr
    llh.res = llh.g
  }
  return(list(res,llh.res,r))
}

### Sample the vent rate
vr.samp<-function(par.g,vr.i){
  vr = par.g$vr[vr.i]
  vr.i.prop = rnorm(1,vr,sd.vr[vr.i])
  par.prop = par.g
  par.prop$vr[vr.i] = vr.i.prop
  fit.prop <- data.frame(ode(y = init.g, times = t, func = SIR, parms = par.prop,method = 'euler'))
  lambda.prop = lambda.calc(fit.prop,par.prop,t)
  lambda.hr.prop = lambda.h.calc(fit.prop,par.prop,t)
  lambda.c.prop = lambda.c.calc(fit.prop,par.prop,t)
  lambda.v.prop = lambda.v.calc(fit.prop,par.prop,t)
  llh.p = llh(lambda.prop,lambda.hr.prop,lambda.c.prop,lambda.v.prop,nb.g,nb.h.g,nb.c.g,nb.v.g)
  prior.vr.g = log(dtruncnorm(x = par.g$vr[vr.i],mean = m.b, sd = ps.b, a = l.b,b = u.b))
  prior.vr.p = log(dtruncnorm(x = vr.i.prop,mean = m.b, sd = ps.b, a = l.b,b = u.b))
  e.a =exp( mpfr(llh.p + prior.vr.p,100) - mpfr(llh.g + prior.vr.g ,100))
  a = min(1,as.numeric(e.a))
  r = rbinom(1,1,a)
  res = vr.i.prop*r + vr*(1-r)
  llh.res = llh.p*r + llh.g*(1-r)
  if(is.na(r)){
    r = 0
    res =  vr
    llh.res = llh.g
  }
  return(list(res,llh.res,r))
}


###Sample the variance of the new admissions
nb.samp<-function(nb){
  log.nb.g = log(nb)
  log.nb.p = rnorm(1,log.nb.g,sd.nb)
  nb.prop = exp(log.nb.p)
  fit.g <- data.frame(ode(y = init.g, times = t, func = SIR, parms = parameters.g,method = 'euler'))
  lambda.g = lambda.calc(fit.g,parameters.g,t)
  lambda.hr.g = lambda.h.calc(fit.g,parameters.g,t)
  lambda.c.g = lambda.c.calc(fit.g,parameters.g,t)
  lambda.v.g = lambda.v.calc(fit.g,parameters.g,t)
  llh.p = llh(lambda.g,lambda.hr.g,lambda.c.g,lambda.v.g,nb.prop,nb.h.g,nb.c.g,nb.v.g)
  prior.g = log(dtruncnorm(x = nb.prop,mean = m.nb, sd = ps.nb, a = l.nb,b = u.nb))
  prior.p = log(dtruncnorm(x = nb.prop,mean = m.nb, sd = ps.nb, a = l.nb,b = u.nb))
  if(prior.p != -Inf){
    e.a =exp( mpfr(llh.p + prior.p + log(1/nb),100) - mpfr(llh.g + prior.g + log(1/nb.prop),100))
    a = min(1,as.numeric(e.a))
    r = rbinom(1,1,a)
  }else{
    r = 0
  }
  res = nb.prop*r + nb*(1-r)
  llh.res = llh.p*r + llh.g*(1-r)
  if(is.na(r)){
    r = 0
    res = nb
    llh.res =llh.g
  }
  return(list(res,llh.res,r))
}

###Sample the variance for the new admissions
nb.h.samp<-function(nb.h){
  log.nb.h.g = log(nb.h)
  log.nb.h.p = rnorm(1,log.nb.h.g,sd.nb.h)
  nb.h.prop = exp(log.nb.h.p)
  fit.g <- data.frame(ode(y = init.g, times = t, func = SIR, parms = parameters.g,method = 'euler'))
  lambda.g = lambda.calc(fit.g,parameters.g,t)
  lambda.hr.g = lambda.h.calc(fit.g,parameters.g,t)
  lambda.c.g = lambda.c.calc(fit.g,parameters.g,t)
  lambda.v.g = lambda.v.calc(fit.g,parameters.g,t)
  llh.p = llh(lambda.g,lambda.hr.g,lambda.c.g,lambda.v.g,nb.g,nb.h.prop,nb.c.g,nb.v.g)
  prior.g = log(dtruncnorm(x = nb.h,mean = m.nb, sd = ps.nb, a = l.nb,b = u.nb))
  prior.p = log(dtruncnorm(x = nb.h.prop,mean = m.nb, sd = ps.nb, a = l.nb,b = u.nb))
  if(prior.p != -Inf){
    e.a =exp( mpfr(llh.p + prior.p + log(1/nb.h),100) - mpfr(llh.g + prior.g + log(1/nb.h.prop),100))
    a = min(1,as.numeric(e.a))
    r = rbinom(1,1,a)
  }else{
    r = 0
  }
  res = nb.h.prop*r + nb.h*(1-r)
  llh.res = llh.p*r + llh.g*(1-r)
  if(is.na(r)){
    r = 0
    res = nb.h
    llh.res = llh.g
  }
  return(list(res,llh.res,r))
}

###Sample the variance for the census
nb.c.samp<-function(nb.c){
  log.nb.c.g = log(nb.c)
  log.nb.c.p = rnorm(1,log.nb.c.g,sd.nb.c)
  nb.c.prop = exp(log.nb.c.p)
  fit.g <- data.frame(ode(y = init.g, times = t, func = SIR, parms = parameters.g,method = 'euler'))
  lambda.g = lambda.calc(fit.g,parameters.g,t)
  lambda.hr.g = lambda.h.calc(fit.g,parameters.g,t)
  lambda.c.g = lambda.c.calc(fit.g,parameters.g,t)
  lambda.v.g = lambda.v.calc(fit.g,parameters.g,t)
  llh.p = llh(lambda.g,lambda.hr.g,lambda.c.g,lambda.v.g,nb.g,nb.h.g,nb.c.prop,nb.v.g)
  prior.g = log(dtruncnorm(x = nb.c,mean = m.nb, sd = ps.nb, a = l.nb,b = u.nb))
  prior.p = log(dtruncnorm(x = nb.c.prop,mean = m.nb, sd = ps.nb, a = l.nb,b = u.nb))
  if(prior.p != -Inf){
    e.a =exp( mpfr(llh.p + prior.p + log(1/nb.c),100) - mpfr(llh.g + prior.g + log(1/nb.c.prop),100))
    a = min(1,as.numeric(e.a))
    r = rbinom(1,1,a)
  }else{
    r = 0
  }
  res = nb.c.prop*r + nb.c*(1-r)
  llh.res = llh.p*r + llh.g*(1-r)
  if(is.na(r)){
    r = 0
    res = nb.c
    llh.res = llh.g
  }
  return(list(res,llh.res,r))
}

###Sample the variance for the census
nb.v.samp<-function(nb.v){
  log.nb.v.g = log(nb.v)
  log.nb.v.p = rnorm(1,log.nb.v.g,sd.nb.v)
  nb.v.prop = exp(log.nb.v.p)
  fit.g <- data.frame(ode(y = init.g, times = t, func = SIR, parms = parameters.g,method = 'euler'))
  lambda.g = lambda.calc(fit.g,parameters.g,t)
  lambda.hr.g = lambda.h.calc(fit.g,parameters.g,t)
  lambda.c.g = lambda.c.calc(fit.g,parameters.g,t)
  lambda.v.g = lambda.v.calc(fit.g,parameters.g,t)
  llh.p = llh(lambda.g,lambda.hr.g,lambda.c.g,lambda.v.g,nb.g,nb.h.g,nb.c.g,nb.v.prop)
  prior.g = log(dtruncnorm(x = nb.v,mean = m.nb, sd = ps.nb, a = l.nb,b = u.nb))
  prior.p = log(dtruncnorm(x = nb.v.prop,mean = m.nb, sd = ps.nb, a = l.nb,b = u.nb))
  if(prior.p != -Inf){
    e.a =exp( mpfr(llh.p + prior.p + log(1/nb.v),100) - mpfr(llh.g + prior.g + log(1/nb.v.prop),100))
    a = min(1,as.numeric(e.a))
    r = rbinom(1,1,a)
  }else{
    r = 0
  }
  res = nb.v.prop*r + nb.v*(1-r)
  llh.res = llh.p*r + llh.g*(1-r)
  if(is.na(r)){
    r = 0
    res = nb.v
    llh.res = llh.g
  }
  return(list(res,llh.res,r))
}

################################
#####  Initilize the MCMC  #####
################################
re.g = c(0,0.1)
gam.g = 1/14
hr.g = rep(-3,numb.hr)
vr.g = rep(-3,numb.vr)
gam.hr.g = 0.1
gam.vr.g = 0.1
gam.e.g = 0.1
parameters.g = list(c(-1,-1,-1,-1.5,-1.5,-1),gam.g,re.g,hr.g, gam.hr.g,vr.g,gam.vr.g)
names(parameters.g) = c('b','gam','re','hr','gam.hr','vr','gam.vr')
p = length(unlist(parameters.g))
I.0 = rep(1,n.mod.cty)
R.0 = rep(0,n.mod.cty)
H.0 = rep(Y.c[!is.na(Y.c)][1]/n.mod.cty,n.mod.cty)
V.0 = rep(Y.v[which(!is.na(Y.v))[1]]/n.mod.cty,n.mod.cty)
S.0 = N.mat[1,] - I.0 -R.0 - H.0 - V.0
init.g = c(S = S.0,I = I.0, R = R.0, H = H.0, V = V.0)
fit.g <- data.frame(ode(y = init.g, times = t, func = SIR, parms = parameters.g,method = 'euler'))
nb.g = 1
nb.h.g = 1
nb.c.g = 1
nb.v.g = 1
S.ind = c(2:(n.mod.cty + 1))
I.ind = c((2 + n.mod.cty):(2*n.mod.cty + 1))
R.ind = c((2 + 2*n.mod.cty):(3*n.mod.cty + 1))
H.ind = c((2+ 3*n.mod.cty):(4*n.mod.cty+1))
V.ind = c((2+ 4*n.mod.cty):(5*n.mod.cty+1))
lambda.g = lambda.calc(fit.g,parameters.g,t)
lambda.hr.g = lambda.h.calc(fit.g,parameters.g,t)
lambda.c.g = lambda.c.calc(fit.g,parameters.g,t)
lambda.v.g = lambda.v.calc(fit.g,parameters.g,t)
llh.g = llh(lambda.g,lambda.hr.g,lambda.c.g,lambda.v.g,nb.g,nb.h.g,nb.c.g,nb.v.g)


####################################
#####  Specify Hyperparamters  #####
####################################

###Support, mean, and standard deviation for spline coefficients and social distancing coefficients
l.b = -Inf
u.b = Inf
m.b = 0
ps.b = 1000

###Support, mean, and standard deviation for recovery, discharge, and vent discharge rates
l.gam = 0
u.gam = 1
m.gam = 1/7
ps.gam = 1000

###Support, mean, and standard deviation for number of initially infected individuals
m.init.inf = 4
l.init.inf = 0
u.init.inf = N.mat[1,]
ps.init.inf = 10000

###Support, mean, and standard deviation for variance parameters
l.nb = 0
u.nb = Inf
m.nb = 1
ps.nb = 10

###############################################################################
#####  Set up MCMC Mechanics (iterations, autotuning, parameter storage)  #####
###############################################################################

### Number of iterations and burn in
G = 5000
burn = 2500

###Track Acceptance rates 
acc.b = rep(0,(p-1))
sd.b = rep(0.1,(p-1))
acc.gam = 0
sd.gam = 0.5
acc.gam.hr = 0
sd.gam.hr = 0.5
acc.gam.vr = 0
sd.gam.vr = 0.5
acc.inf = rep(0,n.mod.cty)
sd.inf = rep(0.5,n.mod.cty)
acc.hr = rep(0,numb.hr)
sd.hr = rep(0.5,numb.hr)
acc.vr = rep(0,numb.vr)
sd.vr = rep(0.5,numb.vr)
acc.nb = 0
sd.nb = 0.5
acc.nb.h = 0
sd.nb.h = 0.5
acc.nb.c = 0
sd.nb.c = 0.5
acc.nb.v = 0
sd.nb.v = 0.5
acc.re = rep(0,n.mod.cty)
sd.re = rep(0.5, n.mod.cty)
par.array = matrix(NA,p,G)
inf.array = matrix(NA,5*(n.mod.cty),G)
hr.array = matrix(NA,numb.hr,G)
nb.array = matrix(NA,1,G)
nb.c.array = matrix(NA,1,G)
nb.h.array = matrix(NA,1,G)
nb.v.array = matrix(NA,1,G)

for(g in 1:G){
  if(is.na(llh.g)){
    break
  }
  ###Sample the transmission parameters
  for(i in 1:numb.b){
    b.res = b.samp(parameters.g,i)
    parameters.g[[1]][i] = b.res[[1]]
    llh.g = b.res[[2]]
    acc.b[i] = acc.b[i] + b.res[[3]]
  }
  ###Sample the Recovery Rate
  gam.res = gam.samp(parameters.g)
  parameters.g[[2]] = gam.res[[1]]
  llh.g = gam.res[[2]]
  acc.gam = acc.gam + gam.res[[3]]
  ###Sample the Discharge Rate
  gam.hr.res = gam.hr.samp(parameters.g)
  parameters.g$gam.hr = gam.hr.res[[1]]
  llh.g = gam.hr.res[[2]]
  acc.gam.hr = acc.gam.hr + gam.hr.res[[3]]
  ###Sample the Vent Discharge Rate
  gam.vr.res = gam.vr.samp(parameters.g)
  parameters.g$gam.vr = gam.vr.res[[1]]
  llh.g = gam.vr.res[[2]]
  acc.gam.vr = acc.gam.vr + gam.vr.res[[3]]
  ###Sample the Random Effects
  for(i in 1:n.mod.cty){
    re.res = re.samp(parameters.g,i)
    parameters.g$re[i] = re.res[[1]]
    llh.g = re.res[[2]]
    acc.re[i] = acc.re[i] + re.res[[3]]
  }
  ###Sample the initial number of infected individuals
  for(i in 1:n.mod.cty){
    init.res = init.inf.samp(init.g,i)
    init.g = init.res[[1]]
    llh.g = init.res[[2]]
    acc.inf[i] = acc.inf[i] + init.res[[3]]
  }
  ###Sample the hospitalization rate
  for(i in 1:numb.hr){
    hr.res = hr.samp(parameters.g,i)
    parameters.g$hr[i] = hr.res[[1]]
    llh.g = hr.res[[2]]
    acc.hr[i] = acc.hr[i] + hr.res[[3]]
  }
  ###Sample the vent rate
  for(i in 1:numb.vr){
    vr.res = vr.samp(parameters.g,i)
    parameters.g$vr[i] = vr.res[[1]]
    llh.g = vr.res[[2]]
    acc.vr[i] = acc.vr[i] + vr.res[[3]]
  }
  ###Sample the infection variance
  nb.res = nb.samp(nb.g)
  nb.g = nb.res[[1]]
  llh.g = nb.res[[2]]
  acc.nb = acc.nb + nb.res[[3]]
  ###Sample the admit variance
  nb.h.res = nb.h.samp(nb.h.g)
  nb.h.g = nb.h.res[[1]]
  llh.g = nb.h.res[[2]]
  acc.nb.h = acc.nb.h + nb.h.res[[3]]
  ###Sample the census variance
  nb.c.res = nb.c.samp(nb.c.g)
  nb.c.g = nb.c.res[[1]]
  llh.g = nb.c.res[[2]]
  acc.nb.c = acc.nb.c + nb.c.res[[3]]
  ###Sample the vent variance
  nb.v.res = nb.v.samp(nb.v.g)
  nb.v.g = nb.v.res[[1]]
  llh.g = nb.v.res[[2]]
  acc.nb.v = acc.nb.v + nb.v.res[[3]]
  ###Store the parameters
  par.array[,g] = unlist(parameters.g)
  inf.array[,g] = init.g
  nb.array[,g] = nb.g
  nb.h.array[,g] = nb.h.g
  nb.c.array[,g] = nb.c.g
  nb.v.array[,g] = nb.v.g
  ###Tune the Proposal Distributions
  if(g%%100 == 0){
    acc.b = acc.b/100
    ind = acc.b <=0.15
    sd.b[ind] = sd.b[ind]*0.9
    ind = acc.b >0.7
    sd.b[ind] = sd.b[ind]*1.1
    acc.b = rep(0,p-1)
    acc.gam = acc.gam/100
    if(acc.gam <=0.15){
      sd.gam = sd.gam*0.9
    }
    if(acc.gam >0.7){
      sd.gam = sd.gam*1.1
    }
    acc.gam = 0
    acc.gam.hr = acc.gam.hr/100
    if(acc.gam.hr <=0.15){
      sd.gam.hr = sd.gam.hr*0.9
    }
    if(acc.gam.hr >0.7){
      sd.gam.hr = sd.gam.hr*1.1
    }
    acc.gam.hr = 0
    acc.gam.vr = acc.gam.vr/100
    if(acc.gam.vr <=0.15){
      sd.gam.vr = sd.gam.vr*0.9
    }
    if(acc.gam.vr >0.7){
      sd.gam.vr = sd.gam.vr*1.1
    }
    acc.gam.vr = 0
    acc.inf = acc.inf/100
    ind = acc.inf <=0.15
    sd.inf[ind] = sd.inf[ind]*0.9
    ind = acc.inf >0.7
    sd.inf[ind] = sd.inf[ind]*1.1
    acc.inf = rep(0,n.mod.cty)
    acc.re = acc.re/100
    ind = acc.re <=0.15
    sd.re[ind] = sd.re[ind]*0.9
    ind = acc.re >0.7
    sd.re[ind] = sd.re[ind]*1.1
    acc.re = rep(0,n.mod.cty)
    acc.hr = acc.hr/100
    ind = acc.hr <=0.15
    sd.hr[ind] = sd.hr[ind]*0.9
    ind = acc.hr >0.7
    sd.hr[ind] = sd.hr[ind]*1.1
    acc.hr = rep(0,numb.hr)
    acc.vr = acc.vr/100
    ind = acc.vr <=0.15
    sd.vr[ind] = sd.vr[ind]*0.9
    ind = acc.vr >0.7
    sd.vr[ind] = sd.vr[ind]*1.1
    acc.vr = rep(0,numb.vr)
    acc.nb = acc.nb/100
    if(acc.nb <=0.15){
      sd.nb = sd.nb*0.9
    }
    if(acc.nb >0.7){
      sd.nb = sd.nb*1.1
    }
    acc.nb = 0
    acc.nb.h = acc.nb.h/100
    if(acc.nb.h <=0.15){
      sd.nb.h = sd.nb.h*0.9
    }
    if(acc.nb.h >0.7){
      sd.nb.h = sd.nb.h*1.1
    }
    acc.nb.h = 0
    acc.nb.c = acc.nb.c/100
    if(acc.nb.c <=0.15){
      sd.nb.c = sd.nb.c*0.9
    }
    if(acc.nb.c >0.7){
      sd.nb.c = sd.nb.c*1.1
    }
    acc.nb.c = 0
    acc.nb.v = acc.nb.v/100
    if(acc.nb.v <=0.15){
      sd.nb.v = sd.nb.v*0.9
    }
    if(acc.nb.v >0.7){
      sd.nb.v = sd.nb.v*1.1
    }
    acc.nb.v = 0
  }
  print(g)
}

############################################################################
#####  Generate the Sample from the Posterior Predictive Distribution  #####
############################################################################
lambda.array =  array(NA,c(tm.pred,n.mod.cty,length(burn:G)))
lambda.pred.array = array(NA,c(length(t.pred),n.mod.cty,length(burn:G)))
lambda.h.array =  array(NA,c(tm.pred,length(burn:G)))
lambda.h.pred.array = array(NA,c(length(t.pred),length(burn:G)))
lambda.c.array =  array(NA,c(tm.pred,length(burn:G)))
lambda.c.pred.array = array(NA,c(length(t.pred),length(burn:G)))
lambda.v.array =  array(NA,c(tm.pred,length(burn:G)))
lambda.v.pred.array = array(NA,c(length(t.pred),length(burn:G)))

for(g in burn:G){
  par.pred = list(par.array[1:numb.b,g],par.array[(numb.b + 1),g], par.array[(numb.b + 2):(numb.b + n.mod.cty + 1),g], par.array[(numb.b + n.mod.cty + 2):(numb.b + n.mod.cty + numb.hr + 1),g], par.array[(numb.b + n.mod.cty + numb.hr + 2),g],par.array[(numb.b + n.mod.cty + numb.hr + 3):(numb.b + n.mod.cty + numb.hr + numb.vr + 2),g], par.array[(numb.b + n.mod.cty + numb.hr + numb.vr + 3),g])
  names(par.pred) = c('b','gam','re','hr','gam.hr','vr','gam.vr')
  init.pred =inf.array[,g]
  fit.pred <- data.frame(ode(y = init.pred, times = t.pred, func = SIR, parms = par.pred,method = 'euler'))
  for(c in 1:n.mod.cty){
    lambda.mean = rnbinom(tm.pred,mu = lambda.calc(fit.pred,par.pred,t.pred)[,c], size = nb.array[,g])
    lambda.pred.array[,c,(g - burn +1)] = lambda.mean
    lambda.array[,c,(g-burn+1)] = lambda.calc(fit.pred,par.pred,t.pred)[1:tm.pred,c]
  }
  lambda.h.mean = lambda.h.calc(fit.pred,par.pred,t.pred)
  lambda.h.array[,(g - burn+1)] = lambda.h.mean#[1:tm]
  lambda.h.pred.array[,(g -burn + 1)] = rnbinom(tm.pred,mu = lambda.h.mean, size = nb.h.array[,g])
  lambda.c.mean = lambda.c.calc(fit.pred,par.pred,t.pred)
  lambda.c.array[,(g - burn+1)] = lambda.c.mean#[1:tm]
  lambda.c.pred.array[,(g -burn + 1)] = rnbinom(tm.pred,mu = lambda.c.mean, size = nb.c.array[,g])
  lambda.v.mean = lambda.v.calc(fit.pred,par.pred,t.pred)
  lambda.v.array[,(g - burn+1)] = lambda.v.mean#[1:tm]
  lambda.v.pred.array[,(g -burn + 1)] = rnbinom(tm.pred,mu = lambda.v.mean, size = nb.v.array[,g])
  print(g)
}
save.image(paste0(out.path,"Image No Res.RData"))

#############################################################################
############  Compute Point Estimates and HPD Intervals  ####################
#############################################################################

lambda.hpd = apply(lambda.pred.array,c(1,2),hdi)
lambda.median = apply(lambda.array,c(1,2),median)
lambda.c.hpd = apply(lambda.c.pred.array,1,hdi)
lambda.c.median = apply(lambda.c.array,1,median)
lambda.v.hpd = apply(lambda.v.pred.array,1,hdi)
lambda.v.median = apply(lambda.v.array,1,median)
lambda.h.hpd = apply(lambda.h.pred.array,1,hdi)
lambda.h.median = apply(lambda.h.array,1,median)

#############################
#####  Write Out Plots  #####
#############################

date.list = seq.Date(from = start.day,to = end.pred.day,by = 1)
for(c in 1:n.mod.cty){
  pdf(paste0(out.path,county.list[c],".pdf"))
  plot(date.list[t.pred],lambda.median[t.pred,c],type = 'l',ylim = range(lambda.hpd[,t.pred,c]), ylab = "Number of Cases", xlab = "Date", main = paste0("Cases: ", county.list[c]), col = 'blue',cex.axis = 1.25,cex.lab = 1.5)
  polygon(c(date.list,rev(date.list)),c(lambda.hpd[1,t.pred,c],rev(lambda.hpd[2,t.pred,c])),col = 'pink',border = 'red')
  points(date.list[t],Y[,c],col= 'red',pch = 16, cex = 1.5)
  points(date.list[t.pred],lambda.median[t.pred,c],type = 'l',ylim = range(lambda.hpd[,t.pred,c]), ylab = "Number of Cases", xlab = "Date", main = paste0("Cases: ", county.list[c]), col = 'blue')
  dev.off()
}

pdf(paste0(out.path,"Non Ventilated Census HPD.pdf"))
plot(date.list[t.pred],lambda.c.median[t.pred],type = 'l',ylim = range(lambda.c.hpd[,t.pred]), ylab = "Number of Cases", xlab = "Date", main = paste0("Non-Ventilated Census"), col = 'blue',cex.axis = 1.25,cex.lab = 1.5)
polygon(c(date.list,rev(date.list)),c(lambda.c.hpd[1,t.pred],rev(lambda.c.hpd[2,t.pred])),col = 'pink',border = 'red')
points(date.list[t],Y.c,col= 'red',pch = 16, cex = 1.5)
points(date.list[t.pred],lambda.c.median[t.pred],type = 'l',col = 'blue')
dev.off()

pdf(paste0(out.path,"Ventilated Census HPD.pdf"))
plot(date.list[t.pred],lambda.v.median[t.pred],type = 'l',ylim = range(lambda.v.hpd[,t.pred]), ylab = "Number of Cases", xlab = "Date", main = paste0("Ventilated Census"), col = 'blue',cex.axis = 1.25,cex.lab = 1.5)
polygon(c(date.list,rev(date.list)),c(lambda.v.hpd[1,t.pred],rev(lambda.v.hpd[2,t.pred])),col = 'pink',border = 'red')
points(date.list[t],Y.v,col= 'red',pch = 16, cex = 1.5)
points(date.list[t.pred],lambda.v.median[t.pred],type = 'l',col = 'blue')
dev.off()

pdf(paste0(out.path,"Admissions.pdf"))
plot(date.list[t.pred],lambda.h.median[t.pred],type = 'l',ylim = range(lambda.h.hpd[,t.pred]), ylab = "Number of Cases", xlab = "Date", main = paste0("Admissions"), col = 'blue',cex.axis = 1.25,cex.lab = 1.5)
polygon(c(date.list,rev(date.list)),c(lambda.h.hpd[1,t.pred],rev(lambda.h.hpd[2,t.pred])),col = 'pink',border = 'red')
points(date.list[t], Y.h,col= 'red',pch = 16, cex = 1.5)
points(date.list[t.pred],lambda.h.median[t.pred],type = 'l',col = 'blue')
dev.off()

