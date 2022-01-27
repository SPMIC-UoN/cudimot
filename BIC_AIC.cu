/* BIC_AIC.cu

   Moises Hernandez-Fernandez - FMRIB Image Analysis Group
   
   Copyright (C) 2005 University of Oxford */

/*  Part of FSL - FMRIB's Software Library
    http://www.fmrib.ox.ac.uk/fsl
    fsl@fmrib.ox.ac.uk

    Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
    Imaging of the Brain), Department of Clinical Neurology, Oxford
    University, Oxford, UK


    LICENCE

    FMRIB Software Library, Release 6.0 (c) 2018, The University of
    Oxford (the "Software")

    The Software remains the property of the Oxford University Innovation
    ("the University").

    The Software is distributed "AS IS" under this Licence solely for
    non-commercial use in the hope that it will be useful, but in order
    that the University as a charitable foundation protects its assets for
    the benefit of its educational and research purposes, the University
    makes clear that no condition is made or to be implied, nor is any
    warranty given or to be implied, as to the accuracy of the Software,
    or that it will be suitable for any particular purpose or for use
    under any specific conditions. Furthermore, the University disclaims
    all responsibility for the use which is made of the Software. It
    further disclaims any liability for the outcomes arising from using
    the Software.

    The Licensee agrees to indemnify the University and hold the
    University harmless from and against any and all claims, damages and
    liabilities asserted by third parties (including claims for
    negligence) which arise directly or indirectly from the use of the
    Software or the sale of any products based on the Software.

    No part of the Software may be reproduced, modified, transmitted or
    transferred in any form or by any means, electronic or mechanical,
    without the express permission of the University. The permission of
    the University is not required if the said reproduction, modification,
    transmission or transference is done without financial return, the
    conditions of this Licence are imposed upon the receiver of the
    product, and all original and amended source code is included in any
    transmitted product. You may be held legally responsible for any
    copyright infringement that is caused or encouraged by your failure to
    abide by these terms and conditions.

    You are not permitted under this Licence to use this Software
    commercially. Use for which any financial return is received shall be
    defined as commercial use, and includes (1) integration of all or part
    of the source code or the Software into a product for sale or license
    by or on behalf of Licensee to third parties or (2) use of the
    Software or any derivative of it for research with the final aim of
    developing software products for sale or license to a third party or
    (3) use of the Software or any derivative of it for research with the
    final aim of developing non-software products for sale or license to a
    third party, or (4) use of the Software to provide any service to an
    external organisation for which payment is received. If you are
    interested in using the Software commercially, please contact Oxford
    University Innovation ("OUI"), the technology transfer company of the
    University, to negotiate a licence. Contact details are:
    fsl@innovation.ox.ac.uk quoting Reference Project 9564, FSL.*/

// BIC & AIC calculation: Bayesian & Akaike Information Criteria

#include "BIC_AIC.h"
#include "functions_gpu.h"
#include "modelparameters.h"
#include "macro_numerical.h"
#include "modelfunctions.h"

namespace Cudimot{

#define VOXELS_BLOCK 8
#define THREADS_VOXEL 32 // Multiple of 32: Threads collaborating to compute a voxel. Do not change this, otherwise Synchronization will be needed and shuffles cannot be used
  
  // Returns the natural log of the 0th order modified Bessel function of first kind for an argument x
  // Follows the exponential implementation of the Bessel function in Numerical Recipes, Ch. 6
  __device__ inline float logIo(const float x){
    float y,b;
    b=fabsf(x);
    if (b<3.75f){
      float a=x/3.75f;
      a*=a;
      //Bessel function evaluation
      y=1.0f+a*(3.5156229f+a*(3.0899424f+a*(1.2067492f+a*(0.2659732f+a*(0.0360768f+a*0.0045813f)))));
      y=logf(y);
    }else{
      float a=3.75f/b; 
      //Logarithm of Bessel function
      y=b+logf((0.39894228f+a*(0.01328592f+a*(0.00225319f+a*(-0.00157565f+a*(0.00916281f+a*(-0.02057706f+a*(0.02635537f+a*(-0.01647633f+a*0.00392377f))))))))/sqrt(b));
    }
    return y;
  }
  
  // Version for Double precision
  __device__ inline double logIo(const double x){
    double y,b;
    b=fabs(x);
    if (b<3.75){
      double a=x/3.75;
      a*=a;
      //Bessel function evaluation
      y=1.0+a*(3.5156229+a*(3.0899424+a*(1.2067492+a*(0.2659732+a*(0.0360768+a*0.0045813)))));
      y=log(y);
    }
    else{
      double a=3.75/b; 
      //Logarithm of Bessel function
      y=b+log((0.39894228+a*(0.01328592+a*(0.00225319+a*(-0.00157565+a*(0.00916281+a*(-0.02057706+a*(0.02635537+a*(-0.01647633+a*0.00392377))))))))/sqrt(b));
    }
    return y;
  }


  template <typename T, bool RICIAN_NOISE>
  __device__ inline  void Compute_BIC_AIC(int idSubVOX,
					  int nmeas,
					  int CFP_Tsize,
					  T* measurements,
					  T* parameters,
					  T* tau,
					  T* CFP,
					  T* FixP,
					  T* BIC,
					  T* AIC)
  {
    int idMeasurement=idSubVOX;
    int nmeas2compute = nmeas/THREADS_VOXEL;
    if (idSubVOX<(nmeas%THREADS_VOXEL)) nmeas2compute++;
    
    T accumulated_error=(T)0.0;
    for(int dir=0;dir<nmeas2compute;dir++){
      T* myCFP = &CFP[idMeasurement*CFP_Tsize];
      T pred_error=Predicted_Signal(NPARAMS,parameters,myCFP,FixP);

      if(RICIAN_NOISE){
	T meas = measurements[idMeasurement];
	pred_error=log_gpu(meas)+(-(T)0.5*(*tau)*(meas*meas+pred_error*pred_error)+logIo((*tau)*pred_error*meas));
	accumulated_error+=pred_error;
      }else{
	pred_error=pred_error-measurements[idMeasurement];
	accumulated_error+=pred_error*pred_error;
      }

      idMeasurement+=THREADS_VOXEL;
    }
    
    #pragma unroll
    for(int offset=THREADS_VOXEL/2; offset>0; offset>>=1){
      accumulated_error+= shfl_down(accumulated_error,offset);
    }
        
    int np=NPARAMS;
    if(idSubVOX==0){
      if(RICIAN_NOISE){
	np++;
	*BIC = nmeas*log_gpu(*tau)+accumulated_error;
      }else{
	*BIC = -(nmeas/(T)2.0)*log_gpu(accumulated_error/(T)2.0);
      }
      *AIC = (T)-2.0 * (*BIC) + np * (T)2.0;
      *BIC = (T)-2.0 * (*BIC) + np * log_gpu((T)nmeas);
    }
  }

  template <typename T, bool RICIAN_NOISE>
  __global__ void bic_aic_kernel(
			     int nmeas, // num measurements
			     int nsamples,
			     int CFP_Tsize, // common fixed params: size*M-measurements
			     int FixP_Tsize, // fixed params: size*N-voxels
			     T* meas, // measurements
			     T* samples, // samples of estimated parameters 
			     T* CFP_global, // common fixed model parameters
			     T* FixP, // fixed model parameters
			     T* BIC_values, // to record BIC values
			     T* AIC_values, // to record AIC values
			     T* tau_samples) // TAU values of each voxel for Rician noise
  {
    // 1 block of threads process several voxels
    // Each warp processes 1 voxel
    int idVOX= (blockIdx.x*VOXELS_BLOCK)+int(threadIdx.x/THREADS_VOXEL);
    int idVOX_inBlock =  threadIdx.x/THREADS_VOXEL;
    int idSubVOX= threadIdx.x%THREADS_VOXEL;
    bool leader = (idSubVOX==0);  // Some steps are performed by only one thread of the warp

    ////////// DYNAMIC SHARED MEMORY ///////////
    extern __shared__ double shared[];		     		// Size: 
    T* CFP = (T*) shared;			 	       	// nmeas*CFP_Tsize
    T* meanSamples = (T*) &CFP[nmeas*CFP_Tsize]; 		// NPARAMS*VOXELS_BLOCK 
    T* TAU = (T*) &meanSamples[NPARAMS*VOXELS_BLOCK];	        // VOXELS_BLOCK 
    T* bic = (T*) &TAU[VOXELS_BLOCK]; 				// VOXELS_BLOCK
    T* aic = (T*) &bic[VOXELS_BLOCK]; 				// VOXELS_BLOCK 
    ////////////////////////////////////////////
    
    /// Copy common fixed model parameters to Shared Memory ///
    if(threadIdx.x==0){ // only one thread of the whole block. Common to all voxels
      for(int i=0;i<nmeas*CFP_Tsize;i++){
	CFP[i]=CFP_global[i];
      }
    }
    ///////////////////////////////////////////////////////////
    
    ///////// each voxel/warp of the block points to its data///////////
    meas = &meas[idVOX*nmeas]; // Global memory
    samples = &samples[idVOX*NPARAMS*nsamples]; // Global memory
    FixP = &FixP[idVOX*FixP_Tsize]; // Global memory
    meanSamples = &meanSamples[idVOX_inBlock*NPARAMS];
    TAU = &TAU[idVOX_inBlock];
    bic = &bic[idVOX_inBlock];
    aic = &aic[idVOX_inBlock];
        
    /// Ititialise shared values of each voxel: only the leader///
    if(leader){ 
      *TAU=(T)0.0;
      if(nsamples>1){
	for(int par=0;par<NPARAMS;par++){
	  T value=(T)0.0;
	  for(int samp=0;samp<nsamples;samp++){
	    value+= samples[par*nsamples+samp];
	  }
	  meanSamples[par]=value/nsamples;
	  if(RICIAN_NOISE){
	    // TAU provided
	    value=(T)0.0;
	    for(int samp=0;samp<nsamples;samp++){
	      value+=tau_samples[idVOX*nsamples+samp];
	    }
	    *TAU=value/nsamples;
	  }
	}
      }else{
        #pragma unroll
	for(int par=0;par<NPARAMS;par++){
	  meanSamples[par]=samples[par];
	}
	if(RICIAN_NOISE){
	  // TAU provided
	  *TAU=tau_samples[idVOX];
	}
      }
    }
    __syncthreads();
    ///////////////////////////////////////////
    
    Compute_BIC_AIC<T,RICIAN_NOISE>(idSubVOX,nmeas,CFP_Tsize,meas,meanSamples,TAU,CFP,FixP,bic,aic);
   
    __syncthreads();

    // Write in global memory the result
    if(leader){
      BIC_values[idVOX]=*bic;
      AIC_values[idVOX]=*aic;
    }
  }
  
  template <typename T>
  BIC_AIC<T>::BIC_AIC(int nvoxFitpart)
  {
    cudimotOptions& opts = cudimotOptions::getInstance();
    nvoxFit_part=nvoxFitpart;
    RicianNoise=(opts.rician.value() && opts.runMCMC.value());
    // if not MCMC executed, ignote --rician, there are not tau samples
  }
  
  template <typename T>
  void BIC_AIC<T>::run(
		   int nvox, int nmeas, int nsamples,
		   int CFP_size, int FixP_size,
		   T* meas,
		   T* samples,
		   T* CFP, T* FixP,
		   T* BIC,
		   T* AIC,
		   T* tau_samples) 
  {
    long int amount_shared_mem = 0;
    amount_shared_mem += (nmeas*CFP_size)*sizeof(T); // CFP
    amount_shared_mem += (NPARAMS*VOXELS_BLOCK)*sizeof(T); // mean_samples
    amount_shared_mem += (1*VOXELS_BLOCK)*sizeof(T); // TAU
    amount_shared_mem += (1*VOXELS_BLOCK)*sizeof(T); // bic
    amount_shared_mem += (1*VOXELS_BLOCK)*sizeof(T); // aic
    
    cout << "Shared Memory used in BIC_AIC kernel: " << amount_shared_mem << endl;
    
    int threads_block = VOXELS_BLOCK * THREADS_VOXEL;
    int nblocks=(nvox/VOXELS_BLOCK);
    if(nvox%VOXELS_BLOCK) nblocks++;
      
    if(RicianNoise){	
      bic_aic_kernel<T,true><<<nblocks,threads_block,amount_shared_mem>>>(nmeas,nsamples,CFP_size,FixP_size,meas,samples,CFP,FixP,BIC,AIC,tau_samples);
    }else{
      bic_aic_kernel<T,false><<<nblocks,threads_block,amount_shared_mem>>>(nmeas,nsamples,CFP_size,FixP_size,meas,samples,CFP,FixP,BIC,AIC,tau_samples);
    }
    
    sync_check("BIC_AIC Kernel");   
  }
 
  // Explicit Instantiations of the template
  template class BIC_AIC<float>;
  template class BIC_AIC<double>;
  
}
