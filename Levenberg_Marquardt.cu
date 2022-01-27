/* Levenberg_Marquardt.cu

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

#include "Levenberg_Marquardt.h"
#include "functions_gpu.h"
#include "modelparameters.h"
#include "macro_numerical.h"
#include "modelfunctions.h"

namespace Cudimot{

#define VOXELS_BLOCK 2
#define THREADS_VOXEL 32 // Multiple of 32: Threads collaborating to compute a voxel. Do not change this, otherwise Synchronization will be needed and shuffles cannot be used
//Dynamic

#define CFTOL 1.0e-8
#define LTOL 1.0e20
#define EPS_gpu 2.0e-16 //Losely based on NRinC 20.1
#define TwoDivPI 0.636619772367581 //0.6353  // 0.636619772367581
#define avoidErrors 1e-4 // Avoid maximum and minimums because tan(pi/2) is undefined to inf
#define SpeedFactor 10 // Speed of Transformations-function

  __constant__ int LMbound_types [NPARAMS];
  __constant__ float LMbounds_min [NPARAMS];
  __constant__ float LMbounds_max [NPARAMS];
  __constant__ int LMfixed [NPARAMS];
  
  __device__ inline bool zero_cf_diff_conv(double* cfold, double* cfnew){
    return(2.0*fabs(*cfold-*cfnew) <= CFTOL*(fabs(*cfold)+fabs(*cfnew)+EPS_gpu));
  }

  template <typename T>
  __device__ inline void MinMaxInvTransform(int idpar, T* params, T* transfParams){

    if(params[idpar]<=((T)LMbounds_min[idpar]+(T)avoidErrors)) params[idpar]+=(T)avoidErrors;
    if(params[idpar]>=((T)LMbounds_max[idpar]-(T)avoidErrors)) params[idpar]-=(T)avoidErrors;
    // Avoid maximum and minimum because tan(pi/2) is undefined to inf

    transfParams[idpar] = params[idpar]-(((T)LMbounds_min[idpar]+(T)LMbounds_max[idpar])/(T)2.0);
    transfParams[idpar] *= ((T)2.0/(((T)LMbounds_max[idpar]-(T)LMbounds_min[idpar])*(T)TwoDivPI));
    transfParams[idpar] = tan_gpu(transfParams[idpar])*SpeedFactor;
  }
  
  template <typename T>
  __device__ inline void MinMaxTransform(int idpar, T* params, T* transfParams){ 
    params[idpar] = (((T)LMbounds_max[idpar]-(T)LMbounds_min[idpar])/(T)2.0);
    params[idpar] *= atan_gpu(transfParams[idpar]/SpeedFactor);
    params[idpar] *=  (T)TwoDivPI;
    params[idpar] += (((T)LMbounds_min[idpar]+(T)LMbounds_max[idpar])/(T)2.0);

    if(params[idpar]<=((T)LMbounds_min[idpar]+(T)avoidErrors)){ 
      params[idpar]+=(T)avoidErrors;
      transfParams[idpar] = params[idpar]-(((T)LMbounds_min[idpar]+(T)LMbounds_max[idpar])/(T)2.0);
      transfParams[idpar] *= ((T)2.0/(((T)LMbounds_max[idpar]-(T)LMbounds_min[idpar])*(T)TwoDivPI));
      transfParams[idpar] = tan_gpu(transfParams[idpar])*SpeedFactor;
    }
    if(params[idpar]>=((T)LMbounds_max[idpar]-(T)avoidErrors)){ 
      params[idpar]-=(T)avoidErrors;
      transfParams[idpar] = params[idpar]-(((T)LMbounds_min[idpar]+(T)LMbounds_max[idpar])/(T)2.0);
      transfParams[idpar] *= ((T)2.0/(((T)LMbounds_max[idpar]-(T)LMbounds_min[idpar])*(T)TwoDivPI));
      transfParams[idpar] = tan_gpu(transfParams[idpar])*SpeedFactor;
    }
    // Avoid maximum and minimum because tan(pi/2) is undefined to inf
  }
  
  template <typename T>
  __device__ inline void MinInvTransform(int idpar, T* params, T* transfParams){
    transfParams[idpar] = log_gpu(params[idpar]-(T)LMbounds_min[idpar]);
  }
  
  template <typename T>
  __device__ inline void MinTransform(int idpar, T* params, T* transfParams){
    params[idpar] = exp_gpu(transfParams[idpar]) + (T)LMbounds_min[idpar];
  }
  template <typename T>
  __device__ inline void MaxInvTransform(int idpar, T* params, T* transfParams){
    transfParams[idpar] = log_gpu((T)LMbounds_max[idpar]-params[idpar]);
  }
  
  template <typename T>
  __device__ inline void MaxTransform(int idpar, T* params, T* transfParams){
    params[idpar] = (T)LMbounds_max[idpar] - exp_gpu(transfParams[idpar]);
  }

  template <typename T>
  __device__ inline void invtransformAll(T* params, T* transfParams){
    #pragma unroll
    for(int p=0; p<NPARAMS; p++){	
      if(LMbound_types[p]==BMIN)
	MinInvTransform(p,params,transfParams);
      else if(LMbound_types[p]==BMAX)
	MaxInvTransform(p,params,transfParams);
      else if(LMbound_types[p]==BMINMAX)
	MinMaxInvTransform(p,params,transfParams);
      else
	transfParams[p]=params[p];
    }
  }

  template <typename T>
  __device__ inline void transformAll(T* params, T* transfParams){
    #pragma unroll
    for(int p=0; p<NPARAMS; p++){	
      if(LMbound_types[p]==BMIN)
	MinTransform(p,params,transfParams);
      else if(LMbound_types[p]==BMAX)
	MaxTransform(p,params,transfParams);
      else if(LMbound_types[p]==BMINMAX)
	MinMaxTransform(p,params,transfParams);
      else
	params[p]=transfParams[p];
    }
  }
  
  template <typename T>
  __device__ inline void derivative_transf(T* transfParams, T* derivatives){
    #pragma unroll
    for(int p=0; p<NPARAMS; p++){
      //if(threadIdx.x==0) printf("Input derivatve %.20f\n",derivatives[p]);
      if(LMbound_types[p]==BMIN)
	derivatives[p]=derivatives[p]*exp_gpu(transfParams[p]);
      else if(LMbound_types[p]==BMAX)
	derivatives[p]=derivatives[p]*(-exp_gpu(transfParams[p]));
      else if(LMbound_types[p]==BMINMAX){
	derivatives[p]=derivatives[p]*((LMbounds_max[p]-LMbounds_min[p])/(T)2.0)*(T)TwoDivPI;
	derivatives[p]=derivatives[p]* ((T)SpeedFactor/ ((T)SpeedFactor*(T)SpeedFactor + transfParams[p]*transfParams[p]));
      }
      //else keep the same
    }
  }
  
  template <typename T, bool DEBUG>
  __device__ inline void Cost_Function(
				       int idSubVOX,
				       int nmeas,
				       int CFP_Tsize,
				       T* measurements,
				       T* parameters,
				       T* CFP,
				       T* FixP,
				       double* result,
				       int debugVOX)
  {
    int idMeasurement=idSubVOX;
    T accumulated_error=(T)0.0;
    
    int nmeas2compute = nmeas/THREADS_VOXEL;
    if (idSubVOX<(nmeas%THREADS_VOXEL)) nmeas2compute++;
    
    for(int iter=0;iter<nmeas2compute;iter++){
      T* myCFP = &CFP[idMeasurement*CFP_Tsize];
      T pred_error=Predicted_Signal(NPARAMS,parameters,myCFP,FixP);
      
      if(DEBUG){
	int idVOX= (blockIdx.x*VOXELS_BLOCK)+int(threadIdx.x/THREADS_VOXEL);
	if(idVOX==debugVOX && idSubVOX==0){
	  printf("PredictedSignal[%i]: %f\n",idMeasurement,pred_error);
	}
      }
      
      pred_error=pred_error-measurements[idMeasurement];
      accumulated_error+=pred_error*pred_error;
      idMeasurement+=THREADS_VOXEL;
    }

    #pragma unroll
    for(int offset=THREADS_VOXEL/2; offset>0; offset>>=1){
      accumulated_error+= shfl_down(accumulated_error,offset);
    }
    if(idSubVOX==0){
      *result=accumulated_error;
      if(DEBUG){
	int idVOX= (blockIdx.x*VOXELS_BLOCK)+int(threadIdx.x/THREADS_VOXEL);
	if(idVOX==debugVOX && idSubVOX==0){
	   printf("COST FUNTION: %f\n",*result);
	}
      }
    }
  }
  
  template <typename T, bool DEBUG>
  __device__ inline void Calculate_Gradient(
					    int idSubVOX,
					    int nmeas,
					    int CFP_Tsize,
					    T* measurements,
					    T* parameters,
					    T* params_transf,
					    T* CFP,
					    T* FixP,
					    T* Gradient,
					    int debugVOX)
  {
    int idMeasurement=idSubVOX;
    T myderivatives[NPARAMS];
    
    int max_iters = nmeas/THREADS_VOXEL;
    if(nmeas%THREADS_VOXEL) max_iters++;

    if(idSubVOX==0){
      #pragma unroll
      for(int p=0;p<NPARAMS;p++){
	Gradient[p]=(T)0.0;
      }
    }
    
    for(int iter=0;iter<max_iters;iter++){
      #pragma unroll
      for(int p=0;p<NPARAMS;p++){
	// Maybe idMeasurement > nmeas, so set to 0
     	myderivatives[p]=(T)0.0; 
      }
      T* myCFP = &CFP[idMeasurement*CFP_Tsize];
      T pred_error=(T)0.0;
      
      if(idMeasurement<nmeas){
	pred_error=Predicted_Signal(NPARAMS,parameters,myCFP,FixP);
	pred_error=pred_error-measurements[idMeasurement];
	
	Partial_Derivatives(NPARAMS,parameters,myCFP,FixP,myderivatives);
	derivative_transf(params_transf,myderivatives);

	if(DEBUG){
	  int idVOX= (blockIdx.x*VOXELS_BLOCK)+int(threadIdx.x/THREADS_VOXEL);
	  if(idVOX==debugVOX && idSubVOX==0){
	    for(int i=0;i<NPARAMS;i++){
	      printf("Derivatives Measurement[%i]_Parameter[%i]: %f\n",idMeasurement,i,myderivatives[i]);
	    }
	  }
	}
      }
      
      #pragma unroll 
      for(int p=0;p<NPARAMS;p++){
	myderivatives[p]=(T)2.0*pred_error*myderivatives[p];
	#pragma unroll 
	for (int offset=THREADS_VOXEL/2; offset>0; offset>>=1){
	  myderivatives[p]+= shfl_down(myderivatives[p],offset);
	}
      }
      if(idSubVOX==0){
        #pragma unroll
	for(int p=0;p<NPARAMS;p++){
	  Gradient[p]+=myderivatives[p];
	  if(LMfixed[p]){
	    Gradient[p]=0;
	  }
	}
      }
      idMeasurement+=THREADS_VOXEL;
    }  
  }

  template <typename T>
  __device__ inline void Calculate_Hessian(
					   int idSubVOX,
					   int nmeas,
					   int CFP_Tsize,
					   T* measurements,
					   T* parameters,
					   T* params_transf,
					   T* CFP,
					   T* FixP,
					   T* Hessian)
  {
    int idMeasurement=idSubVOX;
    T myderivatives[NPARAMS];
    
    int max_dir = nmeas/THREADS_VOXEL;
    if(nmeas%THREADS_VOXEL) max_dir++;
    
    if(idSubVOX==0){
      #pragma unroll
      for(int p=0;p<NPARAMS;p++){
	#pragma unroll
	for(int p2=0;p2<NPARAMS;p2++){
	  Hessian[p*NPARAMS+p2]=(T)0.0;
	}
      }
    }

    for(int iter=0;iter<max_dir;iter++){

      T* myCFP = &CFP[idMeasurement*CFP_Tsize];

      #pragma unroll
      for(int p=0;p<NPARAMS;p++){
	myderivatives[p]=(T)0.0;
      }
      T pred_error=(T)0.0;

      if(idMeasurement<nmeas){
	pred_error=Predicted_Signal(NPARAMS,parameters,myCFP,FixP);
	pred_error=pred_error-measurements[idMeasurement];
	Partial_Derivatives(NPARAMS,parameters,myCFP,FixP,myderivatives);
	derivative_transf(params_transf,myderivatives);
      }
 
      #pragma unroll
      for(int p=0;p<NPARAMS;p++){
	#pragma unroll
	for(int p2=0;p2<NPARAMS;p2++){
	  T element = (T)2.0 * myderivatives[p] * myderivatives[p2];
	  #pragma unroll
	  for (int offset=THREADS_VOXEL/2; offset>0; offset>>=1){
	    element+= shfl_down(element,offset);
	  }
	  if(idSubVOX==0){
	    Hessian[p*NPARAMS+p2]+=element;
	  }
	}
      }
      idMeasurement+=THREADS_VOXEL;
    }
    if(idSubVOX==0){
      #pragma unroll
      for(int p=0;p<NPARAMS;p++){
	#pragma unroll
	for(int p2=0;p2<NPARAMS;p2++){
	  if(LMfixed[p] || LMfixed[p2]){
	    if(p==p2) Hessian[p*NPARAMS+p2]=(T)1.0;
	    else Hessian[p*NPARAMS+p2]=(T)0.0;
	  }
	}
      }
    }
  }
  
  // Hessian * X  =  Gradient -> Calculate X
  template <typename T>
  __device__ inline void LUsolver(int idSubVOX,
				  T* Hessian,
				  T* Gradient,
				  T* Solution){
    
    // If NPARAMS > 32 the current version of the method will fail !!
    // Need to generalise
    
    T col_elems[NPARAMS];
    T pivot;
    
    // Initialise Matrix. Each thread contains a column of the Hessian and one thread the Gradient column:   Matrix = [Hessian | Gradient]
    if (idSubVOX<NPARAMS){
      #pragma unroll
      for(int p=0;p<NPARAMS;p++){
	col_elems[p] = Hessian[p*NPARAMS+idSubVOX];
      }
    }else if(idSubVOX==NPARAMS){
      #pragma unroll
      for(int p=0;p<NPARAMS;p++){
	col_elems[p] = Gradient[p];
      }
    }else{
      #pragma unroll
      for(int p=0;p<NPARAMS;p++){
	col_elems[p] = (T)0.0;
      }
    }
    
    // Solve in two steps: 
    // Forward step: Zero's under diagonal
    // Backward step: Zero's above diagonal

    // Forward step
    #pragma unroll
    for (int col=0; col<NPARAMS; col++){
      // Divide row by diagonal element (1 in the diagonal)
      // Cannot have negative numbers in the diagonal -r/-r = +1
      pivot = col_elems[col];
      pivot = shfl(pivot,col);
      col_elems[col] = col_elems[col]/pivot; 
      
      // Eliminate all terms under diagonal element (1)
      // Pivot is the element to make zero, Pivot-Pivot*1 = 0
      // This_row = This_row - Pivot * row_of_diagonal_element
      #pragma unroll
      for (int row=col+1; row<NPARAMS; row++){
	pivot  = col_elems[row];
	pivot  = shfl(pivot,col);
	col_elems[row] -= pivot*col_elems[col];
      }
    }

    // Backward step
    #pragma unroll
    for (int col=NPARAMS-1; col>0; col--) {
      // Eliminate all terms above diagonal element
      for (int row=0; row<col; row++) {
	pivot  = col_elems[row];
	pivot  = shfl(pivot,col);
	col_elems[row] -= pivot*col_elems[col];
      }
    }

    if(idSubVOX==NPARAMS){
      #pragma unroll
      for(int p=0;p<NPARAMS;p++){
	Solution[p] = col_elems[p];
      }
    }
    
  }
  
  template <typename T, bool MARQUARDT, bool DEBUG>
  __global__ void levenberg_kernel(
				   int nmeas, // nmeasurements
				   int CFP_Tsize, // common fixed params: size*M-measurements
				   int FixP_Tsize, // fixed params: size*Nvoxels 
				   T* meas, // measurements
				   T* parameters, // model parameters 
				   T* CFP_global, // common fixed model parameters
				   T* FixP, // fixed model parameters
				   int nmax_iters,
				   int debugVOX)
{
  // 1 block of threads process several voxels
  // Each warp processes 1 voxel
  int idVOX= (blockIdx.x*VOXELS_BLOCK)+int(threadIdx.x/THREADS_VOXEL);
  int idVOX_inBlock =  threadIdx.x/THREADS_VOXEL;
  int idSubVOX= threadIdx.x%THREADS_VOXEL;
  bool leader = (idSubVOX==0);  // Some steps are performed by only one thread of the warp
    
  ////////// DYNAMIC SHARED MEMORY ///////////
  extern  __shared__ double shared[];				//Size:
  double* pcf = (double*) shared;    				//VOXELS_BLOCK 
  double* ncf = (double*) &pcf[VOXELS_BLOCK];			//VOXELS_BLOCK
  double* lambda = (double*) &ncf[VOXELS_BLOCK];		//VOXELS_BLOCK
  double* olambda = (double*) &lambda[VOXELS_BLOCK];		//VOXELS_BLOCK
  
  T* CFP = (T*) &olambda[VOXELS_BLOCK];				//nmeas*CMP_Tsize
  T* params = (T*) &CFP[nmeas*CFP_Tsize]; 			//NPARAMS*VOXELS_BLOCK
  T* params_transf = (T*) &params[NPARAMS*VOXELS_BLOCK];         //NPARAMS*VOXELS_BLOCK
  T* Gradient = (T*) &params_transf[NPARAMS*VOXELS_BLOCK];     	//NPARAMS*VOXELS_BLOCK
  T* Hessian = (T*) &Gradient[NPARAMS*VOXELS_BLOCK];		//NPARAMS*NPARAMS*VOXELS_BLOCK
  T* step = (T*) &Hessian[NPARAMS*NPARAMS*VOXELS_BLOCK];	//NPARAMS*VOXELS_BLOCK
  
  int* success = (int*) &step[NPARAMS*VOXELS_BLOCK];		//VOXELS_BLOCK
  int* end = (int*) &success[VOXELS_BLOCK];			//VOXELS_BLOCK
  ////////////////////////////////////////////
  
  /// Copy common fixed model parameters to Shared Memory ///
  if(threadIdx.x==0){ // only one thread of the whole block. Common to all voxels
    for(int i=0;i<nmeas*CFP_Tsize;i++){
      CFP[i]=CFP_global[i];
    }
  }
  ///////////////////////////////////////////////////////////

  ///////// each voxel/warp of the block points to its data///////////
  meas = &meas[idVOX*nmeas]; //Global memory
  FixP = &FixP[idVOX*FixP_Tsize]; // Global memory
  pcf = &pcf[idVOX_inBlock];
  ncf = &ncf[idVOX_inBlock];
  lambda = &lambda[idVOX_inBlock];
  olambda = &olambda[idVOX_inBlock];
  params = (T*)&params[idVOX_inBlock*NPARAMS];
  params_transf = &params_transf[idVOX_inBlock*NPARAMS];
  Gradient = &Gradient[idVOX_inBlock*NPARAMS];
  Hessian = &Hessian[idVOX_inBlock*NPARAMS*NPARAMS];
  step = &step[idVOX_inBlock*NPARAMS];
  success = &success[idVOX_inBlock];
  end = &end[idVOX_inBlock];

  int iter=0;
    
  /// Ititialise shared values of each voxel: only the leader///
  if(leader){ 
    *end=false;
    *success=true;
    *lambda=0.1;
    *olambda= 0.0;    
    *ncf=0.0;
    #pragma unroll
    for(int i=0;i<NPARAMS;i++){
      params[i]=parameters[idVOX*NPARAMS+i];
    }
    if(DEBUG){
      if(idVOX==debugVOX){
	printf("\n ----- Levenberg_Marquardt GPU algorithm: voxel %i -----\n",idVOX);
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
	for(int i=0;i< FixP_Tsize;i++){
	  printf("%f, ",FixP[i]);
	}
	printf("\n--------------------------------------------------------\n",idVOX);  
      }
    }
    invtransformAll(params,params_transf); //calculate params_transf  
  }
  // __threadfence_block();
  __syncthreads();
  ///////////////////////////////////////////

  Cost_Function<T,DEBUG>(idSubVOX,nmeas,CFP_Tsize,meas,params,CFP,FixP,pcf,debugVOX);
  if(DEBUG){
    if(idVOX==debugVOX&&leader){
      printf("--------------------------------------------------------\n");  
    }
  }
    
  while (!( (*success) && iter++>=nmax_iters)){
    //if success we don't increase niter (first condition is true)
    //function cost has been decreased, we have advanced.

    if(DEBUG){
      if(idVOX==debugVOX&&leader){
	printf("---------------------- Iteration %i ---------------------\n",iter);
      }
    }

    if(*success){
      Calculate_Gradient<T,DEBUG>(idSubVOX,nmeas,CFP_Tsize,meas,params,params_transf,CFP,FixP,Gradient,debugVOX);
      Calculate_Hessian(idSubVOX,nmeas,CFP_Tsize,meas,params,params_transf,CFP,FixP,Hessian);
    }
        
    if(leader){
      #pragma unroll
      for (int i=0; i<NPARAMS; i++){
	if(MARQUARDT)
	  Hessian[(i*NPARAMS)+i]=((1+(*lambda))/(1+(*olambda)))*Hessian[i*NPARAMS+i];	//Levenberg-Marquardt
	else
	  Hessian[(i*NPARAMS)+i]+=(*lambda)-(*olambda);	//Levenberg
      }
    }

    //__threadfence_block(); // Hessian & Gradient updated
    __syncthreads();    
    LUsolver(idSubVOX,Hessian,Gradient,step);
    // Hessian & Gradient updated // LU solution updated by thread(NPARAMS)
    __syncthreads(); 

    if(leader){
      #pragma unroll
      for(int p=0;p<NPARAMS;p++){
	step[p]=params_transf[p]-step[p];
      }
      transformAll(params,step);

      if(DEBUG){
	if(idVOX==debugVOX){
	  for(int i=0;i<NPARAMS;i++){
	    printf("Gradient[%i]: %f   ",i,Gradient[i]);
	  }
	  printf("\n");
	  for(int i=0;i<NPARAMS;i++){
	    printf("Proposed[%i]: %f   ",i,params[i]);
	  }
	  printf("\n");
	}
      }
    }
    
    //__threadfence_block();  // params updated
    __syncthreads(); 
    Cost_Function<T,DEBUG>(idSubVOX,nmeas,CFP_Tsize,meas,params,CFP,FixP,ncf,debugVOX);
    //__threadfence_block(); // Leader may be faster an update params
    __syncthreads(); 

    if(leader){
      if ( *success = ((*ncf) < (*pcf))){ 
	*olambda = 0.0;
        #pragma unroll
	for(int i=0;i<NPARAMS;i++){
	  params_transf[i]=step[i];
	}
	*lambda=(*lambda)/10.0;
	
	if (zero_cf_diff_conv(pcf,ncf)){
	  *end=true;
	}
	*pcf=*ncf;
      }else{
	*olambda=*lambda;
	*lambda=(*lambda)*10.0;
	if(*lambda > LTOL){ 
	 *end=true;
	}
	// undo step in parameters
	transformAll(params,params_transf);
      }
      if(DEBUG){
	if(idVOX==debugVOX){
	  printf("--------------------------------------------------------\n");  
	}
      }
    }
    //__threadfence_block(); // end,sucess updated
    __syncthreads(); 

    if(*end) break;		
  }
 
  if(leader){
    // Change parameters if needed: FixConstraints()
    FixConstraintsLM(NPARAMS,params);
    // save parameters in global
    #pragma unroll
    for(int i=0;i<NPARAMS;i++){
      parameters[idVOX*NPARAMS+i]=params[i];
    }
    if(DEBUG){
      if(idVOX==debugVOX){
	for(int i=0;i<NPARAMS;i++){
	  printf("Final Parameter[%i]: %f\n",i,params[i]);
	}
      }
    }
  }
}
  
  
  template <typename T>
  Levenberg_Marquardt<T>::Levenberg_Marquardt(vector<int> bou_types, 
					      vector<T> bou_min, 
					      vector<T> bou_max,
					      vector<int> fix)
  {
    cudimotOptions& opts = cudimotOptions::getInstance();
    max_iterations=opts.iterLevMar.value();
    Marquardt=true;
    if(opts.no_Marquardt.value()) Marquardt=false;
    DEBUG=false;
    if(opts.debug.set()){
      DEBUG=true;
      debugVOX= atoi(opts.debug.value().data());
    }

    // Set bounds
    bound_types_host = new int[NPARAMS];
    bounds_min_host = new float[NPARAMS];
    bounds_max_host = new float[NPARAMS];
    
    for(int p=0;p<NPARAMS;p++){
      bound_types_host[p]=bou_types[p];
      bounds_min_host[p]=bou_min[p];
      bounds_max_host[p]=bou_max[p];
    }

    fixed_host = new int[NPARAMS];
    for(int p=0;p<NPARAMS;p++){
       fixed_host[p]=fix[p];
    }

    cudaMemcpyToSymbol(LMbound_types,bound_types_host,NPARAMS*sizeof(int));
    cudaMemcpyToSymbol(LMbounds_min,bounds_min_host,NPARAMS*sizeof(float));
    cudaMemcpyToSymbol(LMbounds_max,bounds_max_host,NPARAMS*sizeof(float));
    cudaMemcpyToSymbol(LMfixed,fixed_host,NPARAMS*sizeof(int));
    sync_check("Levenberg_Marquardt: Setting Bounds - Fixed");
  }
  
  template <typename T>
  void Levenberg_Marquardt<T>::run(
				   int nvox, int nmeas,
				   int CFP_size, int FixP_size,
				   T* meas,
				   T* params,
				   T* CFP, T* FixP) 
  {
    
    long int amount_shared_mem = 0;
    amount_shared_mem += 4*VOXELS_BLOCK*sizeof(double); // Levenberg parameters
    amount_shared_mem += (nmeas*CFP_size)*sizeof(T); // CFP
    amount_shared_mem += (NPARAMS*VOXELS_BLOCK)*sizeof(T); // Parameters
    amount_shared_mem += (NPARAMS*VOXELS_BLOCK)*sizeof(T); // Params_transf
    amount_shared_mem += (NPARAMS*VOXELS_BLOCK)*sizeof(T); // Gradient
    amount_shared_mem += (NPARAMS*NPARAMS*VOXELS_BLOCK)*sizeof(T); // Hessian
    amount_shared_mem += 2*VOXELS_BLOCK*sizeof(int); // Levenberg parameters
    amount_shared_mem += (NPARAMS*VOXELS_BLOCK)*sizeof(T); // step
    
    cout << "Shared Memory used in Levenberg-Marquardt kernel: " << amount_shared_mem << endl;
    
    int threads_block = VOXELS_BLOCK * THREADS_VOXEL;
    int nblocks=(nvox/VOXELS_BLOCK);
    if(nvox%VOXELS_BLOCK) nblocks++;
    
    if(!DEBUG){
      if(Marquardt){
	levenberg_kernel<T,true,false><<<nblocks,threads_block,amount_shared_mem>>>(nmeas,CFP_size,FixP_size,meas,params,CFP,FixP,max_iterations,debugVOX);
      }else{
	levenberg_kernel<T,false,false><<<nblocks,threads_block,amount_shared_mem>>>(nmeas,CFP_size,FixP_size,meas,params,CFP,FixP,max_iterations,debugVOX);
      }
    }else{
      if(Marquardt){
	levenberg_kernel<T,true,true><<<nblocks,threads_block,amount_shared_mem>>>(nmeas,CFP_size,FixP_size,meas,params,CFP,FixP,max_iterations,debugVOX);
      }else{
	levenberg_kernel<T,false,true><<<nblocks,threads_block,amount_shared_mem>>>(nmeas,CFP_size,FixP_size,meas,params,CFP,FixP,max_iterations,debugVOX);
      }
    }
    sync_check("Levenberg_Marquardt Kernel");
  }
  
  // Explicit Instantiations of the template
  template class Levenberg_Marquardt<float>;
  template class Levenberg_Marquardt<double>;
}
