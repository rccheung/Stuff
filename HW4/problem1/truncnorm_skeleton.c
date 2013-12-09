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
      float *x,      // Vector to contain returned samples
      int n,         // Number of samples to return
      float mu,     // Vector of mu's
      float sigma,  // Vector of sigma's
      float a,      // Vector of lower-truncation values
      float b)      // Vector of upper-truncation values
{
    // Usual block/thread indexing...
    int myblock = blockIdx.x + blockIdx.y * gridDim.x;
    int blocksize = blockDim.x * blockDim.y * blockDim.z;
    int subthread = threadIdx.z*(blockDim.x * blockDim.y) + threadIdx.y*blockDim.x + threadIdx.x;
    int idx = myblock * blocksize + subthread;

    float tmp;      // Trying to sample from normal distribution
    float z;        // Exponential random variable
    float phi;      // Calculated probability
    int accept = 0;   // Initialize accept or not
    int maxitr = 1000;  // Maximum try for just sampling in basic normal distribution

  	// Setting RNG
  	curandState state;
  	curand_init(9131 + idx*17, idx, 0, &state);

    // Draw sample
    if (idx < n){
        int i = 0;

        if (isfinite(a) && !isfinite(b)){       //Left truncate
            float mu_ = (a - mu) / sigma;    //standardizing, gives left cut point
            float alpha = ( mu_ + sqrt(mu_ * mu_ + 4)) / 2;  //optimal alpha

            while(!accept && i < maxitr){
                tmp = mu + sigma * curand_normal(&state);  //standard sampling
                if(tmp > a) {
                    x[idx] = tmp;
                    accept = 1;
                }
                else i++;
            }
            while(!accept){
                z = mu_ - log(curand_uniform(&state))/alpha;

                if (mu_ < alpha){
                    phi = exp( -(alpha - z) * (alpha - z) /2);
                }
                else {
                    phi = exp( -(mu_ - alpha) * (mu_ - alpha) / 2 ) * exp(-(alpha - z)*(alpha - z)/2);

                }

                if( curand_uniform(&state) < phi ){
                    x[idx] = z * sigma + mu;
                    accept = 1;
                }
            }
        }
        else if (isfinite(b) && !isfinite(a)){   //Right truncate
            float mu_new = -mu;
            float mu_ = - (b - mu_new) / sigma;
            float alpha = ( mu_ + sqrt(mu_ * mu_ + 4)) / 2;  //optimal alpha

            while(!accept && i < maxitr){
                tmp = mu + sigma * curand_normal(&state);  //standard sampling
                if(tmp < b) {
                    x[idx] = tmp;
                    accept = 1;
                }
                else i++;
            }
            while(!accept){
                z = mu_ - log(curand_uniform(&state))/alpha;

                if (mu_ < alpha){
                    phi = exp( -(alpha - z) * (alpha - z) /2);
                }
                else {
                    phi = exp( -(mu_ - alpha) * (mu_ - alpha) / 2 ) * exp(-(alpha - z)*(alpha - z)/2);

                }

                if( curand_uniform(&state) < phi ){
                    x[idx] = -(z * sigma + mu_new);
                    accept = 1;
                }
            }

        }
        else if (!isfinite(a) && !isfinite(b)){ //No truncation
            x[idx] = mu + sigma * curand_normal(&state);
        }
        else if (isfinite(a) && isfinite(b)){  //Finite truncation
            float mu_ = (a - mu) / sigma;
            float mu_plus = (b - mu) / sigma;

            while(!accept && i < maxitr){
                    tmp = mu + sigma * curand_normal(&state);  //standard sampling
                    if(tmp <= b && tmp >= a) {
                        x[idx] = tmp;
                        accept = 1;
                    }
                    else i++;
                }
            while(!accept){
                float g;
                z = mu_ + (mu_plus-mu_)*curand_uniform(&state);

                if ( 0 >= mu_ && 0 <= mu_plus)
                g = exp ( -z * z / 2);
                else if ( mu_plus < 0)
                g = exp( ( (mu_plus * mu_plus) - (z * z) )/2 );
                else if ( 0 < mu_ )
                g = exp( ( (mu_ * mu_) - (z * z) )/2 );

                if (curand_uniform(&state) < g){
                    x[idx] = z * sigma + mu;
                    accept = 1;
                }

            }

        }
    }
    return;
}
} // END extern "C"

