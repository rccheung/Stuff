
##
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
########################################################################################
## Handle batch job arguments:

# 1-indexed version is used now.
args <- commandArgs(TRUE)

cat(paste0("Command-line arguments:\n"))
print(args)

####
# sim_start ==> Lowest simulation number to be analyzed by this particular batch job
###

#######################
sim_start <- 1000
length.datasets <- 200
#######################

if (length(args)==0){
  sinkit <- FALSE
  sim_num <- sim_start + 1
  set.seed(1330931)
} else {
  # Sink output to file?
  sinkit <- TRUE
  # Decide on the job number, usually start at 1000:
  sim_num <- sim_start + as.numeric(args[1])
  # Set a different random seed for every job number!!!
  set.seed(762*sim_num + 1330931)
}

# Simulation datasets numbered 1001-1200
#sim_num is from 1001 to 1200
########################################################################################
########################################################################################

log.accept.prob <- function(X,y,m,beta.old,beta.new){
  n = length(y)
  sum = 0
  
  for(i in 1:n){
    sum = sum + m[i]*log( (1+exp(X[i,]%*%beta.old)) / (1+exp(X[i,]%*%beta.new)) ) + y[i]*X[i,]%*%(beta.new - beta.old)
  }
  return(sum + 1/2*(t(beta.old) %*% beta.old - t(beta.new) %*% beta.new))
}

"bayes.logreg" <- function(m,y,X,beta.0,Sigma.0.inv,niter=10000,burnin=1000,
                           print.every=1000,retune=100,verbose=TRUE)
{
  n = length(y)
  int.beta = solve(t(X) %*% X) %*% t(X) %*% y      #Initialize Beta (using LS)
  propose.sigma = solve(t(X) %*% X) * sum((y - X%*%int.beta)^2)/(n - 2)  #Use LS estimate for initial variance
  accept = 0 #acceptance rate is accept/niter  
  
  MH.samples = matrix(0, nrow = (niter+burnin), ncol = 2)
  
  #burnin and tuning part
  current.beta = int.beta
  for(i in 1:burnin){
    propose.beta = t(rmvnorm(1,mean = current.beta, sigma = propose.sigma))
    log.u = log(runif(1))
    log.alpha = log.accept.prob(X,y,m,current.beta,propose.beta)
    if(log.u < log.alpha){
      current.beta = propose.beta
      accept = accept + 1
    }
    MH.samples[i,] = current.beta
    
    if((i %% retune) == 0){       #every 100, tune parameter
      accept.rate = accept/i
      if(accept.rate < 0.3){
        propose.sigma = 2.825*accept.rate * propose.sigma
      }else if(accept.rate > 0.6){
        diag(propose.sigma) = (10*accept.rate) * diag(propose.sigma)
      }
    }
  }
  
  #Full sampling
  for(i in (burnin+1):(niter+burnin)){
    propose.beta = t(rmvnorm(1,mean = current.beta, sigma = propose.sigma))
    log.u = log(runif(1))
    log.alpha = log.accept.prob(X,y,m,current.beta,propose.beta)
    if(log.u < log.alpha){
      current.beta = propose.beta
      accept = accept + 1
    }
    MH.samples[i,] = current.beta
    
    if((i %% print.every) == 0)
      cat(paste('The current acceptance rate is ',accept/i,'%\n',sep=""))
  }
  return(MH.samples[(burnin+1):(burnin+niter),])
}


#################################################
# Set up the specifications:
beta.0 <- matrix(c(0,0))
Sigma.0.inv <- diag(rep(1.0,2))
niter <- 10000
#################################################
#effectiveSize(test)
#tt = mcmc(test)
#plot(tt)


# Read data corresponding to appropriate sim_num:
#data = read.csv('C:/Users/user/School Work/Fall 2013/STA 250/HW1/blr_data_1001.csv')
#true.para = read.csv('C:/Users/user/School Work/Fall 2013/STA 250/HW1/blr_pars_1001.csv')
indir = "data/"
in_file = paste0(indir,"blr_data_",sim_num,".csv")
data = read.csv(in_file)

# Extract X and y:
y = data$y; m = data$n; X = as.matrix(data[,3:4])

# Fit the Bayesian model:
test = bayes.logreg(m,y,X)

# Extract posterior quantiles...
quants = seq(0.01, 0.99, by = 0.01)
beta.0.quant = quantile(test[,1], probs = quants)
beta.1.quant = quantile(test[,2], probs = quants)

# Write results to a (99 x p) csv file...
beta.results = matrix(cbind(beta.0.quant, beta.1.quant), ncol = 2)
outdir = "results/"
out_file = paste0(outdir,"blr_res_",sim_num,".csv")
write.table(beta.results, file = out_file, row.names = FALSE, col.names = FALSE, sep = ",")







