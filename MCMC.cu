/* MCMC.cu

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

// Markoc Chain Monte Carlo method

#include "MCMC.h"
#include "functions_gpu.h"
#include "modelparameters.h"
#include "macro_numerical.h"
#include "modelfunctions.h"

namespace Cudimot{

#define VOXELS_BLOCK 8
#define THREADS_VOXEL 32 // Multiple of 32: Threads collaborating to compute a voxel. Do not change this, otherwise Synchronization will be needed and shuffles cannot be used

#define maxfloat 1e10
  
  __constant__ int MCbound_types [NPARAMS];
  __constant__ float MCbounds_min [NPARAMS];
  __constant__ float MCbounds_max [NPARAMS];
  __constant__ int MCprior_types [NPARAMS];
  __constant__ float MCpriors_a [NPARAMS];
  __constant__ float MCpriors_b [NPARAMS];
  __constant__ int MCfixed [NPARAMS];
  
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


  template <typename T, bool DEBUG>
  __device__ inline void Propose(int par, T* params, T* old, T* propSD, curandState* localrandState,int debugVOX){
    *old=params[par];
    if(!MCfixed[par]){
      params[par] = params[par] + curand_normal(localrandState)*propSD[par];
    }
      
    if(DEBUG){
      int idVOX= (blockIdx.x*VOXELS_BLOCK)+int(threadIdx.x/THREADS_VOXEL);
      if(idVOX==debugVOX){
	      printf("----\n",par,params[par]);
	      printf("Proposing Value for Parameter_%i: %f\n",par,params[par]);	
      }
    }
  }

   template <typename T>
   __device__ inline void Compute_TotalPrior(T* priors, T* TotalPrior){
     *TotalPrior=(T)0.0;
     #pragma unroll
     for(int p=0;p<NPARAMS;p++){
       *TotalPrior += priors[p]; 
     }
   }


  template <typename T>
  __device__ inline void Initialize_priors_params(T* params, T* priors, int nmeas, T* CFP, T* FixP){
    #pragma unroll
    for(int p=0;p<NPARAMS;p++){
      if(MCbound_types[p]==BMIN){
	      // Bounded with only min
	      if (params[p] < MCbounds_min[p]) 
	        params[p]=MCbounds_min[p];
      }else if(MCbound_types[p]==BMAX){
	      // Bounded with only max
	      if (params[p] > MCbounds_max[p])
	        params[p]=MCbounds_max[p];
      }else if(MCbound_types[p]==BMINMAX){
	      // Bounded with min & max
	      if (params[p] < MCbounds_min[p])
	        params[p]=MCbounds_min[p];
	      else if (params[p] > MCbounds_max[p])
	      params[p]=MCbounds_max[p];
      }
      // Initialise priors
      if(MCprior_types[p]==1){
	      // Gaussian(mean,std)
	      T std2= MCpriors_b[p]*MCpriors_b[p]; //variance*variance
	      priors[p]=(params[p]-MCpriors_a[p])*(params[p]-MCpriors_a[p])/(2*std2);
	      //prior=(param-mean)*(param-mean)/(2*var*var)
      }else if(MCprior_types[p]==2){
	      // Gamma(alpha,beta)
	      priors[p]= ((T)1.0-MCpriors_a[p])* log_gpu(params[p]) + MCpriors_b[p]*params[p];
      }else if(MCprior_types[p]==3){
	      // ARD(fudge_factor)
	      priors[p]=MCpriors_a[p]*log_gpu(params[p]);
	      //fudgefactor*log(param)
      }else if(MCprior_types[p]==4){
	      // sin()
        priors[p]=-log_gpu(fabs_gpu(sin_gpu(params[p])/(T)2.0));
      }else if(MCprior_types[p]==5){
	      // custom() .. defined by the user in modelfunctions.h
	      priors[p]=custom_priors(p,params,nmeas,CFP,FixP);  
      }else{
	      priors[p]=0;
      }
    }
  }
  
  template <typename T>
  __device__ inline int Check_bounds_constraints(int idpar, T* params){
    // Check Designer Constraints
    if(!ConstraintsMCMC(NPARAMS,params)) return 0;
    if(MCbound_types[idpar]==BMIN){
      // Bounded with only min
      if (params[idpar] < MCbounds_min[idpar]) return 0;
      return 1;
    }else if(MCbound_types[idpar]==BMAX){
      // Bounded with only max
      if (params[idpar] > MCbounds_max[idpar]) return 0;
      return 1;
    }else if(MCbound_types[idpar]==BMINMAX){
      // Bounded with min & max
      if (params[idpar] < MCbounds_min[idpar]) return 0;
      if (params[idpar] > MCbounds_max[idpar]) return 0;
      return 1;
    } 
    return 1;
   }

  template <typename T>
  __device__ inline void Compute_prior(int idpar, T* params, T* priors, T* old_prior, int nmeas, T* CFP, T* FixP){
    *old_prior = priors[idpar];
    if(MCprior_types[idpar]==GAUSSPRIOR){
      // Gaussian(mean,std)
      T std2= MCpriors_b[idpar]*MCpriors_b[idpar]; //std*std
      priors[idpar]=(params[idpar]-MCpriors_a[idpar])*(params[idpar]-MCpriors_a[idpar])/(2*std2);
      //prior=(param-mean)*(param-mean)/(2*std*std)
     
    }else if(MCprior_types[idpar]==GAMMAPRIOR){
      // Gamma(alpha,beta)
      priors[idpar]= ((T)1.0-MCpriors_a[idpar])* log_gpu(params[idpar]) + MCpriors_b[idpar]*params[idpar];

    }else if(MCprior_types[idpar]==ARDPRIOR){
      // ARD(fudge_factor)
      priors[idpar]=MCpriors_a[idpar]*log_gpu(params[idpar]);
      //fudgefactor*log(param)
   
    }else if(MCprior_types[idpar]==SINPRIOR){
      // sin()
      priors[idpar]=-log_gpu(fabs_gpu(sin_gpu(params[idpar])/(T)2.0));

    }else if(MCprior_types[idpar]==CUSTOM){
      // custom() .. defined by the user in modelfunctions.h
      priors[idpar]=custom_priors(idpar,params,nmeas,CFP,FixP);  

    }else{
      priors[idpar]=0;
    }
  }
  
  template <typename T>
  __device__ inline void Initialise_TauRician(int idSubVOX,
					      int nmeas,
					      int CFP_Tsize,
					      T* measurements,
					      T* parameters,
					      T* CFP,
					      T* FixP,
					      T* TAU,
					      T* TAUpropSD)
  {
    int idMeasurement=idSubVOX;
    int nmeas2compute = nmeas/THREADS_VOXEL;
    if (idSubVOX<(nmeas%THREADS_VOXEL)) nmeas2compute++;
    
    T accumulated_error=(T)0.0;
    //Calculate Mean of Residuals
    for(int dir=0;dir<nmeas2compute;dir++){
      T* myCFP = &CFP[idMeasurement*CFP_Tsize];
      T pred_error=Predicted_Signal(NPARAMS,parameters,myCFP,FixP);
      pred_error=pred_error-measurements[idMeasurement];
      accumulated_error+=pred_error;
      idMeasurement+=THREADS_VOXEL;
    }
    #pragma unroll
    for(int offset=THREADS_VOXEL/2; offset>0; offset>>=1){
      accumulated_error+= shfl_down(accumulated_error,offset);
    }
    if(idSubVOX==0){
      *TAU=accumulated_error/nmeas; // this is the mean
    }

    //__threadfence_block(); // TAU changed in shared
    __syncthreads();

    idMeasurement=idSubVOX;
    accumulated_error=(T)0.0;
    //Calculate the square of each difference from the mean 
    for(int dir=0;dir<nmeas2compute;dir++){
      T* myCFP = &CFP[idMeasurement*CFP_Tsize];
      T pred_error=Predicted_Signal(NPARAMS,parameters,myCFP,FixP);
      pred_error=pred_error-measurements[idMeasurement];
      pred_error=(*TAU)-pred_error;
      accumulated_error+=pred_error*pred_error;
      idMeasurement+=THREADS_VOXEL;
    }
    #pragma unroll
    for(int offset=THREADS_VOXEL/2; offset>0; offset>>=1){
      accumulated_error+= shfl_down(accumulated_error,offset);
    }
    //__threadfence_block(); // All completed before change TAU
    __syncthreads();

    if(idSubVOX==0){
      *TAU= accumulated_error/nmeas; //this is the variance
      *TAU = (T)1.0/(*TAU);
      *TAUpropSD = *TAU/(T)2.0;
    }
  }


  template <typename T, bool RICIAN_NOISE, bool DEBUG>
  __device__ inline  void Compute_Likelihood(int idSubVOX,
					     int nmeas,
					     int CFP_Tsize,
					     T* measurements,
					     T* parameters,
					     T* tau,
					     T* CFP,
					     T* FixP,
					     T* likelihood,
					     int debugVOX)
  {
    int idMeasurement=idSubVOX;
    int nmeas2compute = nmeas/THREADS_VOXEL;
    if (idSubVOX<(nmeas%THREADS_VOXEL)) nmeas2compute++;
    
    T accumulated_error=(T)0.0;
    for(int dir=0;dir<nmeas2compute;dir++){
      T* myCFP = &CFP[idMeasurement*CFP_Tsize];
      T pred_error=Predicted_Signal(NPARAMS,parameters,myCFP,FixP);

      if(DEBUG){
	      int idVOX= (blockIdx.x*VOXELS_BLOCK)+int(threadIdx.x/THREADS_VOXEL);
	      if(idVOX==debugVOX && idSubVOX==0){
	        printf("PredictedSignal[%i]: %f\n",idMeasurement,pred_error);
	      }
      }

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
        
    if(idSubVOX==0){
      if(RICIAN_NOISE){
	      *likelihood = -nmeas*log_gpu(*tau)-accumulated_error;
      }else{
	      *likelihood = (nmeas/(T)2.0)*log_gpu(accumulated_error/(T)2.0);
      }
    }
  }
  
  template <typename T, bool DEBUG>
  __device__ inline int Compute_test_energy(T* new_energy, T* old_energy, T* prior, T* likelihood, curandState* localrandState, int debugVOX){
    (*old_energy) = (*new_energy);
    (*new_energy) = (*prior)+ (*likelihood);
    
    T tmp=exp_gpu((*old_energy)-(*new_energy));

    if(DEBUG){
      int idVOX= (blockIdx.x*VOXELS_BLOCK)+int(threadIdx.x/THREADS_VOXEL);
      if(idVOX==debugVOX){
	      printf("OldEnergy(%f), Likelihood(%f) Prior(%f) NewEnergy(%f)\n",*old_energy,*likelihood,*prior,*new_energy);
      }
    }

    return (tmp>curand_uniform(localrandState));
  }
  

  template <typename T, bool RECORDING, bool UPDATE_PROP, bool RICIAN_NOISE, bool DEBUG>
  __global__ void mcmc_kernel(
			      curandState* randstate, // to generate random numbers
			      int nmeas, // num measurements
			      int CFP_Tsize, // common fixed params: size*M-measurements
			      int FixP_Tsize, // fixed params: size*N-voxels
			      int niters,
			      int nsamples, // num samples per parameter
			      int sampleevery, // record a sample every x iterations
			      int updateproposalevery, // update SD proposals every x iters
			      T* meas, // measurements
			      T* parameters, // model parameters 
			      T* propSD_global, // std of proposals
			      T* CFP_global, // common fixed model parameters
			      T* FixP, // fixed model parameters
			      T* samples, // to record parameters samples
			      T* tau_samples, // TAU values of each voxel for Rician noise
			      T* tau_propSD_global, // std of tau proposals
			      int debugVOX)
  {
    // 1 block of threads process several voxels
    // Each warp processes 1 voxel
    int idVOX= (blockIdx.x*VOXELS_BLOCK)+int(threadIdx.x/THREADS_VOXEL);
    int idVOX_inBlock =  threadIdx.x/THREADS_VOXEL;
    int idSubVOX= threadIdx.x%THREADS_VOXEL;
    bool leader = (idSubVOX==0);  // Some steps are performed by only one thread of the warp

    ////////// DYNAMIC SHARED MEMORY ///////////
    extern __shared__ double shared[];		     			// Size: 
    curandState* localrandState = (curandState*)shared;		

    T* CFP = (T*) &localrandState[VOXELS_BLOCK]; 	       	        // nmeas*CFP_Tsize
    T* params = (T*) &CFP[nmeas*CFP_Tsize]; 				// NPARAMS*VOXELS_BLOCK 
    T* priors = (T*) &params[NPARAMS*VOXELS_BLOCK]; 			// NPARAMS*VOXELS_BLOCK
    T* propSD = (T*) &priors[NPARAMS*VOXELS_BLOCK]; 			// NPARAMS*VOXELS_BLOCK
    
    T* TAU = (T*) &propSD[NPARAMS*VOXELS_BLOCK];	              	// VOXELS_BLOCK 
    T* TAUpropSD = (T*) &TAU[VOXELS_BLOCK]; 		    		// VOXELS_BLOCK

    T* likelihood = (T*) &TAUpropSD[VOXELS_BLOCK]; 			// VOXELS_BLOCK 
    T* TotalPrior = (T*) &likelihood[VOXELS_BLOCK]; 			// VOXELS_BLOCK 
    T* energy = (T*) &TotalPrior[VOXELS_BLOCK];				// VOXELS_BLOCK 
    
    T* old_param = (T*) &energy[VOXELS_BLOCK];		 	      	// VOXELS_BLOCK
    T* old_prior  =  (T*) &old_param[VOXELS_BLOCK];		       	// VOXELS_BLOCK 
    T* old_energy =  (T*) &old_prior[VOXELS_BLOCK];		       	// VOXELS_BLOCK 

    int* naccepted = (int*) &old_energy[VOXELS_BLOCK];			// NPARAMS*VOXELS_BLOCK
    int* nrejected = (int*) &naccepted[NPARAMS*VOXELS_BLOCK];		// NPARAMS*VOXELS_BLOCK
    int* TAU_accepted = (int*) &nrejected[NPARAMS*VOXELS_BLOCK];       	// VOXELS_BLOCK
    int* TAU_rejected = (int*) &TAU_accepted[VOXELS_BLOCK];       	// VOXELS_BLOCK
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
    FixP = &FixP[idVOX*FixP_Tsize]; // Global memory
    if(RECORDING){
      samples = &samples[idVOX*NPARAMS*nsamples]; //Global memory
    }
    localrandState = (curandState*)&localrandState[idVOX_inBlock];
    params = &params[idVOX_inBlock*NPARAMS];
    priors = &priors[idVOX_inBlock*NPARAMS];
    propSD = &propSD[idVOX_inBlock*NPARAMS];
    TAU = &TAU[idVOX_inBlock];
    TAUpropSD = &TAUpropSD[idVOX_inBlock];
    likelihood = &likelihood[idVOX_inBlock];
    TotalPrior = &TotalPrior[idVOX_inBlock];
    energy = &energy[idVOX_inBlock];
    old_param = &old_param[idVOX_inBlock];
    old_prior = &old_prior[idVOX_inBlock];
    old_energy = &old_energy[idVOX_inBlock];
    naccepted = &naccepted[idVOX_inBlock*NPARAMS];
    nrejected = &nrejected[idVOX_inBlock*NPARAMS];
    TAU_accepted = &TAU_accepted[idVOX_inBlock];
    TAU_rejected = &TAU_rejected[idVOX_inBlock];
    
    /// Ititialise shared values of each voxel: only the leader///
    if(leader){ 
      *localrandState = randstate[idVOX];
      *TAU=(T)0.0;
      *TAUpropSD=(T)0.0;
      *TAU_accepted=0;
      *TAU_rejected=0;
      #pragma unroll
      for(int par=0;par<NPARAMS;par++){
        params[par]=parameters[idVOX*NPARAMS+par];
        naccepted[par]=0;
        nrejected[par]=0;
        priors[par]=(T)0.0;
        if(RECORDING){
          // already have std of the proposals from previous iterations
          propSD[par]=propSD_global[idVOX*NPARAMS+par];
        }else{
          propSD[par]=params[par]/(T)10.0;
        }
      }
      
      if(RICIAN_NOISE && RECORDING){
        // TAU has been already initializated
        *TAU=tau_samples[idVOX*nsamples];
        *TAUpropSD=tau_propSD_global[idVOX];
      }
      if(DEBUG){
	      if(idVOX==debugVOX){
          printf("\n ----- MCMC GPU algorithm: voxel %i -----\n",idVOX);
          for(int i=0;i<NPARAMS;i++){
            printf("Initial Parameter[%i]: %f\n",i,params[i]);
          }
          for(int i=0;i<CFP_Tsize;i++){
            printf("Commonn Fixed Params[%i]: ",i);
            for(int j=0;j<nmeas;j++){
              printf("%f ",CFP_global[j*CFP_Tsize+i]);
            }
            printf("\n");
          }
          printf("Fix Parameters: ");
          for(int i=0;i<FixP_Tsize;i++){
            printf("%f, ",FixP[i]);
          }
          printf("\n--------------------------------------------------------\n",idVOX);  
	      }
      }
    }
    
    // __threadfence_block();
    __syncthreads();
    ///////////////////////////////////////////
    
    if(leader){
      Initialize_priors_params(params,priors,nmeas,CFP,FixP);
      Compute_TotalPrior(priors,TotalPrior);
    }
    // __threadfence_block();
    __syncthreads(); 

    if(RICIAN_NOISE && !RECORDING){
      // TAU needs to be initializated 
      Initialise_TauRician<T>(idSubVOX,nmeas,CFP_Tsize,meas,params,CFP,FixP,TAU,TAUpropSD);
      // __threadfence_block();
      __syncthreads(); 
    }
  
    Compute_Likelihood<T,RICIAN_NOISE,DEBUG>(idSubVOX,nmeas,CFP_Tsize,meas,params,TAU,CFP,FixP,likelihood,debugVOX);
   
    // __threadfence_block();
    __syncthreads();

    if(leader){ 
      *energy=(*TotalPrior)+(*likelihood);

      if(DEBUG){
        if(idVOX==debugVOX){
          printf("Initial Point: Likelihood(%f) Prior(%f) Energy(%f)\n",*likelihood,*TotalPrior,*energy);
          printf("--------------------------------------------------------\n");  
        }
      }
    }
    
    // __threadfence_block();
    __syncthreads();

    int criteria=0;
    
    for(int iter=0; iter<niters; iter++){

      if(DEBUG){
        if(idVOX==debugVOX&&leader){
          printf("---------------------- Iteration %i ---------------------\n",iter);
        }
      }

      // Propose Tau if Rician Noise
      if(RICIAN_NOISE){
        criteria=0;
        if(leader){
          *old_param=*TAU;
          *TAU = (*TAU) + curand_normal(localrandState)*(*TAUpropSD);
          if(DEBUG){
            if(idVOX==debugVOX){
              printf("Proposing Value for TAU: %f",*TAU);	
            }
	        }
	        criteria= (*TAU>(T)0.0);
	      }
	      criteria = shfl(criteria,0);
	      // __threadfence_block(); // TAU modified by leader
	      __syncthreads();
	
        if(criteria){
          Compute_Likelihood<T,1,DEBUG>(idSubVOX,nmeas,CFP_Tsize,meas,params,TAU,CFP,FixP,likelihood,debugVOX);
        }
        //__threadfence_block(); // all threads must have finished before modify TAU
        __syncthreads();
      
        if(leader){
          if(criteria){
            Compute_TotalPrior(priors,TotalPrior);
            criteria=Compute_test_energy<T,DEBUG>(energy,old_energy,TotalPrior,likelihood,localrandState,debugVOX);
            if(criteria){
              (*TAU_accepted)++;
              if(DEBUG){
                if(idVOX==debugVOX){
                  printf("Accepted TAU\n");	
                }
              }
            }else{
              (*TAU_rejected)++;
              *TAU=(*old_param);
              *energy=*old_energy;
              if(DEBUG){
                if(idVOX==debugVOX){
                  printf("Rejected TAU\n");	
                }
              }
            }
          }else{
            (*TAU_rejected)++;
            *TAU=(*old_param);
            if(DEBUG){
              if(idVOX==debugVOX){
                printf("Rejected TAU\n");	
              }
            }
          }
        }	
      }

      // Propose the rest of Parameters
      for(int par=0; par<NPARAMS; par++){
	      criteria=0;
	
        if(leader){
          Propose<T,DEBUG>(par,params,old_param,propSD,localrandState,debugVOX);
          criteria=Check_bounds_constraints(par,params);
        }	
        criteria = shfl(criteria,0);
        // __threadfence_block(); //params is modified by leader
        __syncthreads();

        if(criteria){
          Compute_Likelihood<T,RICIAN_NOISE,DEBUG>(idSubVOX,nmeas,CFP_Tsize,meas,params,TAU,CFP,FixP,likelihood,debugVOX);
        }  
	      // __threadfence_block(); // params cannot be modify until all threads finish
	      __syncthreads();
 
        if(leader){
          if(criteria){
            Compute_prior(par,params,priors,old_prior,nmeas,CFP,FixP);
            Compute_TotalPrior(priors,TotalPrior);
            criteria=Compute_test_energy<T,DEBUG>(energy,old_energy,TotalPrior,likelihood,localrandState,debugVOX);
            if(criteria){
              naccepted[par]++;
              if(DEBUG){
                if(idVOX==debugVOX){
                  printf("Accepted Parameter_%i\n",par);	
                }
	            }
	          }else{
              nrejected[par]++;
              params[par]=(*old_param);
              priors[par]=(*old_prior);
              *energy=*old_energy;
	            if(DEBUG){
                if(idVOX==debugVOX){
                  printf("Rejected Parameter_%i\n",par);	
                }
	            }
	          }
	        }else{
            nrejected[par]++;
            params[par]=(*old_param);
	          if(DEBUG){
	            if(idVOX==debugVOX){
                printf("Rejected Parameter_%i\n",par);	
              }
	          }
	        }
	      }
      }
      
      // Record Samples
      if(RECORDING){
	      if((!(iter%sampleevery))&&(leader)){
	        int nsamp=iter/sampleevery;
	        #pragma unroll
	        for(int par=0; par<NPARAMS; par++){
	          samples[par*nsamples+nsamp]=params[par];
	        }
	        if(RICIAN_NOISE){
	          tau_samples[idVOX*nsamples+nsamp]=*TAU;
	        }
	      }
      }

      // Update propsals Std
      if(!RECORDING || UPDATE_PROP){  // deactivated when not recording if --no_updateproposal
	      if((iter>0)&&(!(iter%updateproposalevery))&&(leader)){
          #pragma unroll
	        for(int par=0; par<NPARAMS; par++){
            propSD[par]*=sqrt((naccepted[par]+(T)1.0)/(nrejected[par]+(T)1.0));
            propSD[par]=min_gpu(propSD[par],(T)maxfloat);
            naccepted[par]=0;
            nrejected[par]=0;
	        }
	        if(RICIAN_NOISE){
            (*TAUpropSD)*=sqrt(((*TAU_accepted)+(T)1.0)/(*(TAU_rejected)+(T)1.0));
            (*TAUpropSD)=min_gpu(*TAUpropSD,(T)maxfloat);
            *TAU_accepted=0;
            *TAU_rejected=0;
	        }
	      }
      }

      if(DEBUG){
	      if(idVOX==debugVOX&&leader){
          for(int i=0;i<NPARAMS;i++){
            printf("Parameter[%i]: %f\n",i,params[i]);
          }
	        printf("--------------------------------------------------------\n");  
	      }
      }

    }  // end Iterations

    // Write in global memory before finishing the data that is needded next call
    if(leader){
      randstate[idVOX]=*localrandState; 
      // save state, otherwise random numbers will be repeated (start at the same point)
      #pragma unroll
      for(int par=0;par<NPARAMS;par++){
      	parameters[idVOX*NPARAMS+par]=params[par];
        if(!RECORDING){
          propSD_global[idVOX*NPARAMS+par]=propSD[par];
        }
      }
      if(DEBUG){
	      if(idVOX==debugVOX){
          for(int i=0;i<NPARAMS;i++){
            printf("Final Parameter[%i]: %f\n",i,params[i]);
          }
	      }
	    }
      if(RICIAN_NOISE && !RECORDING){
        //save TAU and TAUpropSD
        tau_samples[idVOX*nsamples]=*TAU;
        tau_propSD_global[idVOX]=*TAUpropSD;
      }
    }
  }
  
  __global__ void setup_randoms_kernel(curandState* randstate, double seed){
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    curand_init(seed,id,0,&randstate[id]);
  }
  
  template <typename T>
  MCMC<T>::MCMC(int nvoxFitpart,
		vector<int> bou_types, 
		vector<T> bou_min, 
		vector<T> bou_max,
		vector<int> pri_types, 
		vector<T> pri_a, 
		vector<T> pri_b,
		vector<int> fix)
  {
    cudimotOptions& opts = cudimotOptions::getInstance();
    nvoxFit_part=nvoxFitpart;
    nburnin=opts.nburn.value();
    if(nburnin<1) nburnin=1; // to set the std of the proposals
    njumps=opts.njumps.value();
    nsamples=opts.njumps.value()/opts.sampleevery.value();
    sampleevery=opts.sampleevery.value();
    if(njumps<sampleevery){
      cerr << "CUDIMOT Error: If MCMC is run, the method must record at least one sample. The number of jumps (--nj) must be greater than or equal to space between samples (--se)" << endl; 
      exit(-1);
    }
    updateproposalevery=opts.sampleevery.value();
    updateproposal=true;
    if(opts.no_updateproposal.value()) updateproposal=false;
    RicianNoise=opts.rician.value();

    DEBUG=false;
    if(opts.debug.set()){
      DEBUG=true;
      debugVOX= atoi(opts.debug.value().data());
    }

    // Set bounds - priors
    bound_types_host = new int[NPARAMS];
    bounds_min_host = new float[NPARAMS];
    bounds_max_host = new float[NPARAMS];
    for(int p=0;p<NPARAMS;p++){
      bound_types_host[p]=bou_types[p];
      bounds_min_host[p]=bou_min[p];
      bounds_max_host[p]=bou_max[p];
    }
    prior_types_host = new int[NPARAMS];
    priors_a_host = new float[NPARAMS];
    priors_b_host = new float[NPARAMS];
    for(int p=0;p<NPARAMS;p++){
      prior_types_host[p]=pri_types[p];
      priors_a_host[p]=pri_a[p];
      priors_b_host[p]=pri_b[p];
    }
    fixed_host = new int[NPARAMS];
    for(int p=0;p<NPARAMS;p++){
      fixed_host[p]=fix[p];
    }

    cudaMemcpyToSymbol(MCbound_types,bound_types_host,NPARAMS*sizeof(int));
    cudaMemcpyToSymbol(MCbounds_min,bounds_min_host,NPARAMS*sizeof(float));
    cudaMemcpyToSymbol(MCbounds_max,bounds_max_host,NPARAMS*sizeof(float));
    
    cudaMemcpyToSymbol(MCprior_types,prior_types_host,NPARAMS*sizeof(int));
    cudaMemcpyToSymbol(MCpriors_a,priors_a_host,NPARAMS*sizeof(float));
    cudaMemcpyToSymbol(MCpriors_b,priors_b_host,NPARAMS*sizeof(float));

    cudaMemcpyToSymbol(MCfixed,fixed_host,NPARAMS*sizeof(int));
    sync_check("MCMC: Setting Bounds - Priors - Fixed");

    //Allocate mem for proposal SD on GPU
    cudaMalloc((void**)&propSD, nvoxFit_part*NPARAMS*sizeof(T));
    cudaMalloc((void**)&tau_propSD, nvoxFit_part*sizeof(T));

    // Initialise Randoms
    int blocks_Rand = nvoxFit_part/256;
    if(nvoxFit_part%256) blocks_Rand++;
    cudaMalloc((void**)&randStates, blocks_Rand*256*sizeof(curandState));
    dim3 Dim_Grid_Rand(blocks_Rand,1);
    dim3 Dim_Block_Rand(256,1);
    srand(opts.seed.value()+opts.idPart.value());  //randoms seed
    setup_randoms_kernel<<<Dim_Grid_Rand,Dim_Block_Rand>>>(randStates,rand());
    sync_check("Setup_Randoms_kernel");
  }
  
  template <typename T>
  void MCMC<T>::run(
		    int nvox, int nmeas,
		    int CFP_size, int FixP_size,
		    T* meas,
		    T* params,
		    T* CFP, T* FixP,
		    T* samples,
		    T* tau_samples) 
  {
    long int amount_shared_mem = 0;
    amount_shared_mem += VOXELS_BLOCK*sizeof(curandState); // curandState
    amount_shared_mem += (nmeas*CFP_size)*sizeof(T); // CFP
    amount_shared_mem += (NPARAMS*VOXELS_BLOCK)*sizeof(T); // Parameters
    amount_shared_mem += (NPARAMS*VOXELS_BLOCK)*sizeof(T); // Priors
    amount_shared_mem += (NPARAMS*VOXELS_BLOCK)*sizeof(T); // PropSD
    amount_shared_mem += (2*VOXELS_BLOCK)*sizeof(T); // TAU, TAU_PropSD
    amount_shared_mem += (3*VOXELS_BLOCK)*sizeof(T); // Likelihod,TPrior,Energy
    amount_shared_mem += (3*VOXELS_BLOCK)*sizeof(T); // old_param, old_prior, old_energy
    amount_shared_mem += (NPARAMS*VOXELS_BLOCK)*sizeof(int); // naccepted
    amount_shared_mem += (NPARAMS*VOXELS_BLOCK)*sizeof(int); // nrejected
    amount_shared_mem += (2*VOXELS_BLOCK)*sizeof(int); // TAU_accepted, TAU_rejected
    
    cout << "Shared Memory used in MCMC kernel: " << amount_shared_mem << endl;
    
    int threads_block = VOXELS_BLOCK * THREADS_VOXEL;
    int nblocks=(nvox/VOXELS_BLOCK);
    if(nvox%VOXELS_BLOCK) nblocks++;
    
    // Burn-In   ... always update_proposals
    if(RicianNoise){
      if(DEBUG){
	      mcmc_kernel<T,false,true,true,true><<<nblocks,threads_block,amount_shared_mem>>>(randStates,nmeas,CFP_size,FixP_size,nburnin,nsamples,sampleevery,updateproposalevery,meas,params,propSD,CFP,FixP,samples,tau_samples,tau_propSD,debugVOX);
      }else{
	      mcmc_kernel<T,false,true,true,false><<<nblocks,threads_block,amount_shared_mem>>>(randStates,nmeas,CFP_size,FixP_size,nburnin,nsamples,sampleevery,updateproposalevery,meas,params,propSD,CFP,FixP,samples,tau_samples,tau_propSD,debugVOX);
      }
    }else{
      if(DEBUG){
	      mcmc_kernel<T,false,true,false,true><<<nblocks,threads_block,amount_shared_mem>>>(randStates,nmeas,CFP_size,FixP_size,nburnin,nsamples,sampleevery,updateproposalevery,meas,params,propSD,CFP,FixP,samples,tau_samples,tau_propSD,debugVOX);
      }else{
	      mcmc_kernel<T,false,true,false,false><<<nblocks,threads_block,amount_shared_mem>>>(randStates,nmeas,CFP_size,FixP_size,nburnin,nsamples,sampleevery,updateproposalevery,meas,params,propSD,CFP,FixP,samples,tau_samples,tau_propSD,debugVOX);
      }
    }
    sync_check("MCMC Kernel: burnin step");
    
    
    // Recordig
    if(updateproposal){
      if(RicianNoise){
	      if(DEBUG){
	        mcmc_kernel<T,true,true,true,true><<<nblocks,threads_block,amount_shared_mem>>>(randStates,nmeas,CFP_size,FixP_size,njumps,nsamples,sampleevery,updateproposalevery,meas,params,propSD,CFP,FixP,samples,tau_samples,tau_propSD,debugVOX);
	      }else{
	        mcmc_kernel<T,true,true,true,false><<<nblocks,threads_block,amount_shared_mem>>>(randStates,nmeas,CFP_size,FixP_size,njumps,nsamples,sampleevery,updateproposalevery,meas,params,propSD,CFP,FixP,samples,tau_samples,tau_propSD,debugVOX);
	      }
      }else{
	      if(DEBUG){
	        mcmc_kernel<T,true,true,false,true><<<nblocks,threads_block,amount_shared_mem>>>(randStates,nmeas,CFP_size,FixP_size,njumps,nsamples,sampleevery,updateproposalevery,meas,params,propSD,CFP,FixP,samples,tau_samples,tau_propSD,debugVOX);
	      }else{
	        mcmc_kernel<T,true,true,false,false><<<nblocks,threads_block,amount_shared_mem>>>(randStates,nmeas,CFP_size,FixP_size,njumps,nsamples,sampleevery,updateproposalevery,meas,params,propSD,CFP,FixP,samples,tau_samples,tau_propSD,debugVOX);
	      }
      }
      
    }else{ // no_updateproposal
      if(RicianNoise){
	      if(DEBUG){
	        mcmc_kernel<T,true,false,true,true><<<nblocks,threads_block,amount_shared_mem>>>(randStates,nmeas,CFP_size,FixP_size,njumps,nsamples,sampleevery,updateproposalevery,meas,params,propSD,CFP,FixP,samples,tau_samples,tau_propSD,debugVOX);
	      }else{
	        mcmc_kernel<T,true,false,true,false><<<nblocks,threads_block,amount_shared_mem>>>(randStates,nmeas,CFP_size,FixP_size,njumps,nsamples,sampleevery,updateproposalevery,meas,params,propSD,CFP,FixP,samples,tau_samples,tau_propSD,debugVOX);
	      }
      }else{
	      if(DEBUG){
	        mcmc_kernel<T,true,false,false,true><<<nblocks,threads_block,amount_shared_mem>>>(randStates,nmeas,CFP_size,FixP_size,njumps,nsamples,sampleevery,updateproposalevery,meas,params,propSD,CFP,FixP,samples,tau_samples,tau_propSD,debugVOX);
	      }else{
	        mcmc_kernel<T,true,false,false,false><<<nblocks,threads_block,amount_shared_mem>>>(randStates,nmeas,CFP_size,FixP_size,njumps,nsamples,sampleevery,updateproposalevery,meas,params,propSD,CFP,FixP,samples,tau_samples,tau_propSD,debugVOX);
	      }
      }
    }
    sync_check("MCMC Kernel: recording step"); 
  }

  // Explicit Instantiations of the template
  template class MCMC<float>;
  template class MCMC<double>;
}
