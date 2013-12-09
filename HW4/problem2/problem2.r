#Problem 2
library(msm)
library(MASS)
library(compiler)
library(RCUDA)

###############################################################
"compute_grid" <- function(N,sqrt_threads_per_block=16L,grid_nd=1)
{
  # if...
  # N = 1,000,000
  # => 1954 blocks of 512 threads will suffice
  # => (62 x 32) grid, (512 x 1 x 1) blocks
  # Fix block dims:
  block_dims <- c(as.integer(sqrt_threads_per_block), as.integer(sqrt_threads_per_block), 1L)
  threads_per_block <- prod(block_dims)
  if (grid_nd==1){
    grid_d1 <- as.integer(max(1L,ceiling(N/threads_per_block)))
    grid_d2 <- 1L
  } else {
    grid_d1 <- as.integer(max(1L, floor(sqrt(N/threads_per_block))))
    grid_d2 <- as.integer(ceiling(N/(grid_d1*threads_per_block)))
  }
  grid_dims <- c(grid_d1, grid_d2, 1L)
  return(list("grid_dims"=grid_dims,"block_dims"=block_dims))
}
###############################################################

#Part a
probit_mcmc_cpu = function(y, X, beta_0, Sigma_0_inv, niter, burnin){
  #Initializations
  n = nrow(X); p = ncol(X); total = burnin + niter;
  y_1 = sum(y == 1); y_0 = sum(y == 0)
  var = solve(Sigma_0_inv + t(X) %*% X); priors = Sigma_0_inv %*% beta_0; 
  samples = matrix(0, nrow = p, ncol = total)
  samples[,1] = beta_0;
  z = matrix(0, nrow = n);
  for(i in 2:total){
    means = X %*% samples[,(i-1)]
    z[y == 0] = rtnorm(y_0, mean = means[y == 0], sd = rep(1, y_0), lower = -Inf, upper = 0)
    z[y == 1] = rtnorm(y_1, mean = means[y == 1], sd = rep(1, y_1), lower = 0, upper = Inf)
    samples[,i] = mvrnorm(1, mu = var %*% (priors + t(X) %*% z), Sigma = var)
  }
  return(samples)
}
probit_mcmc_cpu = cmpfun(probit_mcmc_cpu)


#Part b
m <- loadModule('truncnorm2.ptx')
mykernel <- m$truncnorm_kernel

probit_mcmc_gpu = function(y, X, beta_0, Sigma_0_inv, niter, burnin)
{
  #Initializations
  n = nrow(X); p = ncol(X); total = burnin + niter;
  tmp = solve(Sigma_0_inv + t(X) %*% X); sigma = rep(1, n)
  samples = matrix(0, nrow = p, ncol = total)
  samples[,1] = beta_0;
  z = matrix(0, nrow = n);
  y_1 = sum(y == 1); y_0 = sum(y == 0);
  a = rep(0, n); a[y == 0] = -Inf;
  b = rep(0, n); b[y == 1] = Inf; 
  var = solve(Sigma_0_inv + t(X) %*% X); priors = Sigma_0_inv %*% beta_0; 
  samples = matrix(0, nrow = p, ncol = total)
  samples[,1] = beta_0;
  z = numeric(n);
  bg <- compute_grid(N=n); grid_dims <- bg$grid_dims; block_dims <- bg$block_dims
  for(i in 2:total){
    means = X %*% samples[,(i-1)]
    z = .cuda(mykernel, "z" = z, n, y, means, sigma, a, b, gridDim=grid_dims, blockDim=block_dims, outputs = "z")
    z = matrix(z)
    samples[,i] = mvrnorm(1, mu = var %*% (priors + t(X) %*% z), Sigma = var)
  }
  return(samples)
}
probit_mcmc_gpu = cmpfun(probit_mcmc_gpu)

#Analyzing small data set.
beta_0 = matrix(0,nrow = 8); Sigma_inv = diag(8)

data = read.table('mini_data.txt',header = T)
y = data[,1]; X = as.matrix(data[,2:9])
testcpu = probit_mcmc_cpu(y,X,beta_0, Sigma_inv, 2000, 500)
testgpu = probit_mcmc_gpu(y,X,beta_0, Signa_inv, 2000, 500)
GLM_beta = glm(y ~ X_2 + X_3 + X_4 + X_5 + X_6 + X_7 + X_8, data = data, family = binomial(link = 'probit') )

#Analyzing big data sets.
data = read.table('data_01.txt', header = T)
y = data[,1]; X = as.matrix(data[,2:9])
gpu_time1 = system.time(gpu_beta.sample_01 <- probit_mcmc_gpu(y, X, beta_0, Sigma_inv, 2000,500))
cpu_time1 = system.time(cpu_beta.sample_01 <- probit_mcmc_cpu(y, X, beta_0, Sigma_inv, 2000, 500))

data = read.table('data_02.txt', header = T)
y = data[,1]; X = as.matrix(data[,2:9])
gpu_time2 = system.time(gpu_beta.sample_02 <- probit_mcmc_gpu(y, X, beta_0, Sigma_inv, 2000,500))
cpu_time2 = system.time(cpu_beta.sample_02 <- probit_mcmc_cpu(y, X, beta_0, Sigma_inv, 2000, 500))

data = read.table('data_03.txt',header = T)
y = data[,1]; X = as.matrix(data[,2:9])
gpu_time3 = system.time(gpu_beta.sample_03 <- probit_mcmc_gpu(y, X, beta_0, Sigma_inv, 2000,500))
cpu_time3 = system.time(cpu_beta.sample_03 <- probit_mcmc_cpu(y, X, beta_0, Sigma_inv, 2000, 500))

data = read.table('data_04.txt', header = T)
y = data[,1]; X = as.matrix(data[,2:9])
gpu_time4 = system.time(gpu_beta.sample_04 <- probit_mcmc_gpu(y, X, beta_0, Sigma_inv, 20,5))
cpu_time4 = system.time(cpu_beta.sample_04 <- probit_mcmc_cpu(y, X, beta_0, Sigma_inv, 2000, 500))

data = read.table('data_05.txt', header = T)
y = data[,1]; X = as.matrix(data[,2:9])
gpu_time5 = system.time(gpu_beta.sample_05 <- probit_mcmc_gpu(y, X, beta_0, Sigma_inv, 1,1))
cpu_time5 = system.time(cpu_beta.sample_05 <- probit_mcmc_cpu(y, X, beta_0, Sigma_inv, 2000, 500))