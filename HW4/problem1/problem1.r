#Problem 1
library(RCUDA)
library(msm)
################################################
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
###############################################


m <- loadModule('truncnorm_skeleton.ptx')
mykernel <- m$truncnorm_kernel

#Part c
N = 10000
threads_per_block <- 512L
block_dims <- c(16L, 16L, 1L)
grid_dims <- c(40L, 1L, 1L)

cat("Grid size:\n")
print(grid_dims)
cat("Block size:\n")
print(block_dims)

nthreads <- prod(grid_dims)*prod(block_dims) 
cat("Total number of threads to launch = ",nthreads,"\n")
if (nthreads < N){
  stop("Grid is not large enough...!")
}


mu = 2; sigma = 1; a = 0; b = 1.5;
x = rep(0.0, N);
x_mem <- copyToDevice(x)
.cuda(mykernel, x_mem, N, mu, sigma, a, b, gridDim=grid_dims, blockDim=block_dims)
gpu_1c = copyFromDevice(obj=x_mem,nels=x_mem$nels,type="float")

#Comparing sample means with theoretical mean
mean(gpu_1c)
z = pnorm((b-mu)/sigma) - pnorm((a-mu)/sigma)
theory_mean = mu + (dnorm((a-mu)/sigma) - dnorm((b-mu)/sigma))/z*sigma


#Part d
r_1 = rtnorm(N, mean = mu, sd = sigma, lower = a, upper = b)
mean(r_1)

#Part e
#run time for each n
n = 10^(1:8)
todevicetime = vector('list', 8)
cudatime = vector('list',8)
fromdevicetime = vector('list', 8)
cputime = vector('list',8)
for(i in 1:8){
  #GPU
  x = rep(0.0, n[i]);
  bg <- compute_grid(N=n[i])
  grid_dims <- bg$grid_dims
  block_dims <- bg$block_dims
  todevicetime[[i]] = system.time(x_mem <- copyToDevice(x))
  cudatime[[i]] = system.time(.cuda(mykernel, x_mem, n[i], mu, sigma, a, b, gridDim=grid_dims, blockDim=block_dims))
  fromdevicetime[[i]] = system.time(gpu_1c <- copyFromDevice(obj=x_mem,nels=x_mem$nels,type="float"))
  cputime[[i]] = system.time(rtnorm(n[i], mean = mu, sd = sigma, lower = a, upper = b))
}

#Part f
n = 10000
x = rep(0.0, n)
a = -Inf; b = Inf;
bg <- compute_grid(N=n)
grid_dims <- bg$grid_dims
block_dims <- bg$block_dims
inf_gpu = .cuda(mykernel, 'result' = x, n, mu, sigma, a, b, gridDim=grid_dims, blockDim=block_dims, outputs = 'result')
inf_cpu = rtnorm(n, mean = mu, sd = sigma, lower = a, upper = b)

#Part g
a = -Inf; b = -10; mu = 0; sigma = 1;
x = rep(0.0, n)
lowtrun_gpu = .cuda(mykernel, 'result' = x, n, mu, sigma, a, b, gridDim=grid_dims, blockDim=block_dims, outputs = 'result')
lowtrun_cpu = rtnorm(n, mean = mu, sd = sigma, lower = a, upper = b)
