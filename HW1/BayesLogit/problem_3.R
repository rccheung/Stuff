#Problem 3
#
# Logistic regression
# 
# Y_{i} | \beta \sim \textrm{Bin}\left(n_{i},e^{x_{i}^{T}\beta}/(1+e^{x_{i}^{T}\beta})\right)
# \beta \sim N\left(\beta_{0},\Sigma_{0}\right)
#
##
library(mvtnorm)
library(coda)
library(MCMCpack)

########################################################################################

log.accept.prob <- function(X,y,m,Sigma.0.inv,beta.old,beta.new){
  n = length(y)
  sum = 0
  
  for(i in 1:n){
    sum = sum + m[i]*log( (1+exp(X[i,]%*%beta.old)) / (1+exp(X[i,]%*%beta.new)) ) + y[i]*X[i,]%*%(beta.new - beta.old)
  }
  return(sum + 1/2*(t(beta.old) %*% Sigma.0.inv %*% beta.old - t(beta.new) %*% Sigma.0.inv %*% beta.new))
}

"bayes.logregMG" <- function(m,y,X,beta.0,Sigma.0.inv,type, niter=50000,burnin=10000,
                           print.every=5000,retune=1000,verbose=TRUE)
{
  n = length(y)
  if(type == "LS"){
    int.beta = matrix(lm(y~X-1)$coefficients)      #Initialize Beta (using LS)
    propose.sigma = vcov(lm(y~X-1))
  }else{
    int.beta = matrix(glm(y~X-1, family=binomial(logit))$coefficients)
    propose.sigma = vcov(glm(y~X-1, family = binomial(logit)))
  }
  small = c(1,5,6,9,10,11); med = c(3,4); high = c(2,7,8)
  smallvar = propose.sigma[small,small]; medvar = propose.sigma[med,med]; highvar = propose.sigma[high,high]
  accept = numeric(3) #acceptance rate is accept/niter  
  k = 1  #Index for updating acceptance rate
  MH.samples = matrix(0, nrow = (niter+burnin), ncol = 11)
  
  #burnin and tuning part
  #do each as a single marginal
  current.beta = int.beta
  propose.beta = current.beta
  for(i in 1:burnin){
    
    #Small
    propose.beta[small] = t(rmvnorm(1,mean = current.beta[small], sigma = smallvar))
    log.u = log(runif(1))
    log.alpha = log.accept.prob(X,y,m,Sigma.0.inv,current.beta,propose.beta)
    if(log.u < log.alpha){
      current.beta = propose.beta   #update
      accept[k] = accept[k] + 1
    }else{propose.beta = current.beta}   #don't update
    k = k+1
    
    #Mediam
    propose.beta[med] = t(rmvnorm(1,mean = current.beta[med], sigma = medvar))
    log.u = log(runif(1))
    log.alpha = log.accept.prob(X,y,m,Sigma.0.inv,current.beta,propose.beta)
    if(log.u < log.alpha){
      current.beta = propose.beta   #update
      accept[k] = accept[k] + 1
    }else{propose.beta = current.beta}   #don't update
    k = k+1
    
    #high
    propose.beta[high] = t(rmvnorm(1,mean = current.beta[high], sigma = highvar))
    log.u = log(runif(1))
    log.alpha = log.accept.prob(X,y,m,Sigma.0.inv,current.beta,propose.beta)
    if(log.u < log.alpha){
      current.beta = propose.beta   #update
      accept[k] = accept[k] + 1
    }else{propose.beta = current.beta}   #don't update
    k = 1
    
    MH.samples[i,] = current.beta
    
    if((i %% retune) == 0){       #every 100, tune parameter
      accept.rate = accept/i      #tune for intercept
      cat(paste('The current acceptance rate at ',i,'is ',100*round(accept/i,2),'%\n',sep=""))
      
      #Small
      if(accept.rate[1] < 0.2){
        smallvar = 2.825*accept.rate[1] * smallvar
      }else if(accept.rate[1] > 0.45){
        diag(smallvar) = (100*accept.rate[1]) * diag(smallvar)
      }
      
      #Perimater and radius
      if(accept.rate[2] < 0.2){
        medvar = 2.825*accept.rate[2] * medvar
      }else if(accept.rate[2] > 0.45){
        diag(medvar) = (50*accept.rate[2]) * diag(medvar)
      }
      
      #Perimater and radius
      if(accept.rate[3] < 0.2){
        highvar = 2.825*accept.rate[3] * highvar
      }else if(accept.rate[3] > 0.45){
        diag(highvar) = (10*accept.rate[3]) * diag(highvar)
      }
    }   
  }     #end of burnin period
  
  print(smallvar); print(medvar); print(highvar);
  #Full sampling
  for(i in (burnin+1):(niter+burnin)){
    
    #Small
    propose.beta[small] = t(rmvnorm(1,mean = current.beta[small], sigma = smallvar))
    log.u = log(runif(1))
    log.alpha = log.accept.prob(X,y,m,Sigma.0.inv,current.beta,propose.beta)
    if(log.u < log.alpha){
      current.beta = propose.beta   #update
      accept[k] = accept[k] + 1
    }else{propose.beta = current.beta}   #don't update
    k = k+1
    
    #Mediam
    propose.beta[med] = t(rmvnorm(1,mean = current.beta[med], sigma = medvar))
    log.u = log(runif(1))
    log.alpha = log.accept.prob(X,y,m,Sigma.0.inv,current.beta,propose.beta)
    if(log.u < log.alpha){
      current.beta = propose.beta   #update
      accept[k] = accept[k] + 1
    }else{propose.beta = current.beta}   #don't update
    k = k+1
    
    #high
    propose.beta[high] = t(rmvnorm(1,mean = current.beta[high], sigma = highvar))
    log.u = log(runif(1))
    log.alpha = log.accept.prob(X,y,m,Sigma.0.inv,current.beta,propose.beta)
    if(log.u < log.alpha){
      current.beta = propose.beta   #update
      accept[k] = accept[k] + 1
    }else{propose.beta = current.beta}   #don't update
    k = 1
    
    MH.samples[i,] = current.beta
    
    if((i %% print.every) == 0)
      cat(paste('The current acceptance rate at ',i,'is ',100*round(accept/i,2),'%\n',sep=""))
  }
  return(list(MH.samples = MH.samples, smallvar = smallvar, med = medvar, highvar = highvar))
}

"bayes.logregMH" <- function(m,y,X,beta.0,Sigma.0.inv,type,niter=20000,burnin=2000,
                           print.every=2000,retune=200,verbose=TRUE)
{
  n = length(y)
  if(type == "LS"){
    int.beta = matrix(lm(y~X-1)$coefficients)      #Initialize Beta (using LS)
    propose.sigma = vcov(lm(y~X-1))
  }else{
    int.beta = matrix(glm(y~X-1, family=binomial(logit))$coefficients)
    propose.sigma = vcov(glm(y~X-1, family = binomial(logit)))
  }
  accept = 0 #acceptance rate is accept/niter  
  logs = matrix(0,nrow=(burnin+niter), ncol = 2)
  MH.samples = matrix(0, nrow = (niter+burnin), ncol = 11)
  
  #burnin and tuning part
  #do each as a single marginal
  current.beta = int.beta
  propose.beta = current.beta
  for(i in 1:burnin){
    
    propose.beta = t(rmvnorm(1,mean = current.beta, sigma = propose.sigma))
    log.u = log(runif(1))
    log.alpha = log.accept.prob(X,y,m,Sigma.0.inv,current.beta,propose.beta)
    logs[i,] = c(log.u, log.alpha)
    if(log.u < log.alpha){
      current.beta = propose.beta   #update
      accept = accept + 1
    }else{propose.beta = current.beta}   #don't update
    
    MH.samples[i,] = current.beta
    
    if((i %% retune) == 0){       #every 100, tune parameter
      accept.rate = accept/i      #tune for intercept
      cat(paste('The current acceptance rate at ',i,'is ',100*round(accept/i,2),'%\n',sep=""))
      
      if(accept.rate < 0.2){
        propose.sigma = 2.825*accept.rate * propose.sigma
      }else if(accept.rate > 0.45){
        diag(propose.sigma) = (10*accept.rate) * diag(propose.sigma)
      }
    }
  }
  
  print(propose.sigma)
  #Full sampling
  for(i in (burnin+1):(niter+burnin)){
    
    propose.beta = t(rmvnorm(1,mean = current.beta, sigma = propose.sigma))
    log.u = log(runif(1))
    log.alpha = log.accept.prob(X,y,m,Sigma.0.inv,current.beta,propose.beta)
    logs[i,] = c(log.u, log.alpha)
    if(log.u < log.alpha){
      current.beta = propose.beta   #update
      accept = accept + 1
    }else{propose.beta = current.beta}   #don't update
    
    MH.samples[i,] = current.beta
    
    if((i %% print.every) == 0)
      cat(paste('The current acceptance rate at ',i,'is ',100*round(accept/i,2),'%\n',sep=""))
  }
  return(list(MH.samples, logs, propose.sigma))
}

#################################################
# Set up the specifications:
beta.0 <- matrix(rep(0,11))
Sigma.0.inv <- diag(rep(0.001,11))
#################################################

#Data process
data = read.table('C:/Users/user/School Work/Fall 2013/STA 250/HW1/breast_cancer.txt', header = T)
# Extract X and y:
n = nrow(data)
y = data$diagnosis; m = rep(1,n); X = cbind(rep(1,n),as.matrix(data[,1:10]))
y = as.numeric(y=='M')
newx = matrix(0, nrow = 569, ncol = 11)
for(i in 2:11){
  newx[,i] = (X[,i] - mean(X[,i]))/sd(X[,i])
}
newx[,1] = rep(1,569)
X = newx
# Fit the Bayesian model:
fitMG = bayes.logregMG(m,y,X,beta.0,Sigma.0.inv, ,type = "GLM",burnin = 5000, retune = 500, niter = 100000)
fitMH = bayes.logregMH(m,y,X,beta.0,Sigma.0.inv, type = "GLM", burnin = 5000, retune = 500, niter = 100000)

##################################################
#Analysis
effectiveSize(fitMG[[1]][5001:105000,])
effectiveSize(fitMH[[1]][5001:105000,])
MG.sample = fitMG[[1]][5001:105000,]
MH.sample = fitMH[[1]][5001:105000,]
MG.mcmc = mcmc(MG.sample)
MH.mcmc = mcmc(MH.sample)
plot(MG.mcmc); plot(MH.mcmc);

#Compute lag1
lag1.MG = unlist(sapply(1:11, function(i) acf(MG.sample[,i],plot=F,lag.max=1))[1,])
lag1.MH = unlist(sapply(1:11, function(i) acf(MH.sample[,i],plot=F,lag.max=1))[1,])

#Thinning
lag1.MG = unlist(sapply(1:11, function(i) acf(MG.sample[seq(10,100000,by=10),i],plot=F,lag.max=1))[1,])
lag1.MH = unlist(sapply(1:11, function(i) acf(MH.sample[seq(10,100000,by=10),i],plot=F,lag.max=1))[1,])

##################
#Check covariates
quants = c(0.025, 0.975)
cov.res = sapply(1:11, function(i) quantile(MG.sample[,i],probs = quants))
cov.res = sapply(1:11, function(i) quantile(MH.sample[,i],probs = quants))

#Predictive Check
#MG
tail.MG = MG.sample[90000:100000,]
s = split(tail.MG, tail.MG[,1])
cc = length(s)
samplebetas = matrix(0,nrow = 11, ncol = cc)
for(i in 1:cc){
	samplebetas[,i] = s[[i]][1:11]
}
probabilities = X %*% samplebetas
probabilities = exp(probabilities)/(1+exp(probabilities))
for(i in 1:cc){    #replace NaN with 1
	probabilities[ which( is.nan(probabilities[,i]) ),i] = 1
}
simulations = matrix(0,nrow = n, ncol = ncol(probabilities))
for(i in 1:cc){
	simulations[,i] = rbinom(n,1,prob = probabilities[,i])
}
means = apply(simulations,2,mean)
jpeg('predictivecheck_MG.jpg')
truemean = mean(y)
hist(means, main = "PPC for MG results"); abline(v = truemean)
dev.off()

#MH
tail.MH = MH.sample[90000:100000,]
s = split(tail.MH, tail.MH[,1])
cc = length(s)
samplebetas = matrix(0,nrow = 11, ncol = cc)
for(i in 1:cc){
  samplebetas[,i] = s[[i]][1:11]
}
probabilities = X %*% samplebetas
probabilities = exp(probabilities)/(1+exp(probabilities))
for(i in 1:cc){    #replace NaN with 1
  probabilities[ which( is.nan(probabilities[,i]) ),i] = 1
}
simulations = matrix(0,nrow = n, ncol = ncol(probabilities))
for(i in 1:cc){
  simulations[,i] = rbinom(n,1,prob = probabilities[,i])
}
means = apply(simulations,2,mean)
jpeg('predictivecheck_MH.jpg')
truemean = mean(y)
hist(means, main = "PPC for MH results"); abline(v = truemean)
dev.off()
