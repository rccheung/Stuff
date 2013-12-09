#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>
#include <curand_kernel.h>
#include <math_constants.h>

extern "C"
{
__global__ void
 truncnorm_kernel(
      float *z,      // Vector to contain returned samples
      int n,         // Number of samples to return
      float *y,     // Vector of y's
      float *mu,     // Vector of mu's
      float *sigma,  // Vector of sigma's
      float *a,      // Vector of lower-truncation values
      float *b)      // Vector of upper-truncation values
{
    // Usual block/thread indexing...
    int myblock = blockIdx.x + blockIdx.y * gridDim.x;
    int blocksize = blockDim.x * blockDim.y * blockDim.z;
    int subthread = threadIdx.z*(blockDim.x * blockDim.y) + threadIdx.y*blockDim.x + threadIdx.x;
    int idx = myblock * blocksize + subthread;

    float tmp;      // Trying to sample from normal distribution
    float phi;      // Calculated probability
    int accept = 0;   // Initialize accept or not
    int maxitr = 2000;  // Maximum try for just sampling in basic normal distribution

  	// Setting RNG
  	curandState state;
  	curand_init(9131 + idx*17, idx, 0, &state);

    // Draw sample
    if (idx < n){
        int i = 0;

        if (isfinite(a[idx]) && !isfinite(b[idx])){       //Left truncate
            float mu_ = (a[idx] - mu[idx]) / sigma[idx];    //standardizing, gives left cut point
            float alpha = ( mu_ + sqrt(mu_ * mu_ + 4)) / 2;  //optimal alpha

            while(!accept && i < maxitr){
                tmp = mu[idx] + sigma[idx] * curand_normal(&state);  //standard sampling
                if(tmp > a[idx]) {
                    z[idx] = tmp;
                    accept = 1;
                }
                else i++;
            }
            while(!accept){
                tmp = mu_ - log(curand_uniform(&state))/alpha;

                if (mu_ < alpha){
                    phi = exp( -(alpha - tmp) * (alpha - tmp) /2);
                }
                else {
                    phi = exp( -(mu_ - alpha) * (mu_ - alpha) / 2 ) * exp(-(alpha - tmp)*(alpha - tmp)/2);

                }

                if( curand_uniform(&state) < phi ){
                    z[idx] = tmp * sigma[idx] + mu[idx];
                    accept = 1;
                }
            }
        }
        else if (isfinite(b[idx]) && !isfinite(a[idx])){   //Right truncate
            float mu_new = -mu[idx];
            float mu_ = - (b[idx] - mu_new) / sigma[idx];
            float alpha = ( mu_ + sqrt(mu_ * mu_ + 4)) / 2;  //optimal alpha

            while(!accept && i < maxitr){
                tmp = mu[idx] + sigma[idx] * curand_normal(&state);  //standard sampling
                if(tmp < b[idx]) {
                    z[idx] = tmp;
                    accept = 1;
                }
                else i++;
            }
            while(!accept){
                tmp = mu_ - log(curand_uniform(&state))/alpha;

                if (mu_ < alpha){
                    phi = exp( -(alpha - tmp) * (alpha - tmp) /2);
                }
                else {
                    phi = exp( -(mu_ - alpha) * (mu_ - alpha) / 2 ) * exp(-(alpha - tmp)*(alpha - tmp)/2);

                }

                if( curand_uniform(&state) < phi ){
                    z[idx] = -(tmp * sigma[idx] + mu_new);
                    accept = 1;
                }
            }

        }
    }
    return;
}
} // END extern "C"

