/* GridSearch.cu

   Moises Hernandez-Fernandez - FMRIB Image Analysis Group
   
   Copyright (C) 2005 University of Oxford */

#include "GridSearch.h"
#include "functions_gpu.h"
#include "modelparameters.h"
#include "macro_numerical.h"
#include "modelfunctions.h"

using namespace std;

namespace Cudimot{

#define VOXELS_BLOCK 8
#define THREADS_VOXEL 32 // Multiple of 32: Threads collaborating to compute a voxel. Do not change this, otherwise Synchronization will be needed

  __constant__ int gridParams [NPARAMS]; // may not use all, but max is NPARAMS
  __constant__ int GSbound_types [NPARAMS];
  __constant__ float GSbounds_min [NPARAMS];
  __constant__ float GSbounds_max [NPARAMS];

  template <typename T>
  __device__ inline bool checkBounds(T* params)
  {
    #pragma unroll
    for(int p=0;p<NPARAMS;p++){
      if(GSbound_types[p]==BMIN){
	// Bounded with only min
	if (params[p] < GSbounds_min[p]) 
	  return false;
      }else if(GSbound_types[p]==BMAX){
	// Bounded with only max
	if (params[p] > GSbounds_max[p])
	   return false;
      }else if(GSbound_types[p]==BMINMAX){
	// Bounded with min & max
	if (params[p] < GSbounds_min[p])
	  return false;
	else if (params[p] > GSbounds_max[p])
	  return false;
      }
    }
    return true;
  }

  template <typename T, bool DEBUG>
  __device__ inline void Cost_Function(int idSubVOX,
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
  __global__ void gridSearch_kernel(
				    int nGridParams,
				    int gridCombs,
				    int nmeas, // nmeasurements
				    int CFP_Tsize, // common fixed params: size*M-measurements
				    int FixP_Tsize, // fixed params: size*Nvoxels 
				    T* meas, // measurements
				    T* grid, // values to try
				    T* parameters, // model parameters 
				    T* CFP_global, // common fixed model parameters
				    T* FixP, // fixed model parameters
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
    T* CFP = (T*) &ncf[VOXELS_BLOCK];		       		//nmeas*CMP_Tsize
    T* params = (T*) &CFP[nmeas*CFP_Tsize]; 			//NPARAMS*VOXELS_BLOCK
    T* trial_params = (T*) &params[NPARAMS*VOXELS_BLOCK];       //NPARAMS*VOXELS_BLOCK
    int* boundsTest = (int*) &trial_params[NPARAMS*VOXELS_BLOCK]; //VOXELS_BLOCK
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
    params = (T*)&params[idVOX_inBlock*NPARAMS];
    trial_params = (T*)&trial_params[idVOX_inBlock*NPARAMS];
    boundsTest = (int*)&boundsTest[idVOX_inBlock];
    
    /// Ititialise shared values of each voxel: only the leader///
    if(leader){
      #pragma unroll
      for(int i=0;i<NPARAMS;i++){
	params[i]=parameters[idVOX*NPARAMS+i];
      }
      if(DEBUG){
	if(idVOX==debugVOX){
	  printf("\n ----- GridSearch GPU algorithm: voxel %i -----\n",idVOX);
	  for(int i=0;i<NPARAMS;i++){
	    printf("Initial Parameter[%i]: %f\n",i,params[i]);
	  }
	  for(int i=0;i<CFP_Tsize;i++){
	    printf("Common Fixed Params[%i]: ",i);
	    for(int j=0;j<nmeas;j++){
	      printf("%f ",CFP_global[j*CFP_Tsize+i]);
	    }
	    printf("\n");
	  }
	  printf("Fix Parameters: ");
	  for(int i=0;i< FixP_Tsize;i++){
	    printf("%f, ",FixP[i]);
	  }
	  printf("\n--------------------------------------------------------\n");  
	}
      }
    }
    // __threadfence_block();
    __syncthreads();
    ///////////////////////////////////////////
    
    //Cost_Function<T,DEBUG>(idSubVOX,nmeas,CFP_Tsize,meas,params,CFP,FixP,pcf,debugVOX);
    if(leader){
      *pcf=9e20;
    }
    if(DEBUG){
      if(idVOX==debugVOX&&leader){
	printf("--------------------------------------------------------\n");  
      }
    }
    
    for(int comb=0;comb<gridCombs;comb++){
      if(leader){
        #pragma unroll
	for(int i=0; i<NPARAMS; i++){
	  trial_params[i]= params[i];
	}
	for(int i=0;i<nGridParams;i++){
	  trial_params[gridParams[i]]= grid[(comb*nGridParams)+i];
	}
	*boundsTest=checkBounds(trial_params);
      }
    
      if(DEBUG){
	if(idVOX==debugVOX&&leader){
	  printf("---------------------- Combination %i ---------------------\n",comb);
	  printf("Parameters: "); 
	  for(int i=0;i<NPARAMS;i++) printf("%f ",trial_params[i]);
	  printf("\nBoundsTest: %i\n",*boundsTest);
	}
      }
      
      //__threadfence_block();
      __syncthreads();    
      
      if(*boundsTest){
	Cost_Function<T,DEBUG>(idSubVOX,nmeas,CFP_Tsize,meas,trial_params,CFP,FixP,ncf,debugVOX);
	//__threadfence_block(); // Leader may be faster an update params
	__syncthreads(); 
      
	if(leader){
	  if ((*ncf) < (*pcf)){ 
            #pragma unroll
	    for(int i=0;i<NPARAMS;i++){
	      params[i]=trial_params[i];
	    }
	    *pcf=*ncf;
	  }
	}
      }
      if(DEBUG){
	if(leader&&idVOX==debugVOX){
	  printf("--------------------------------------------------------\n");  
	}
      }
      //__threadfence_block();
      __syncthreads(); 
    }
    
    if(leader){
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
  GridSearch<T>::GridSearch(int nGP, vector<int> gP, int gC, T* grid_host,
			    vector<int> bou_types, vector<T> bou_min, 
			    vector<T> bou_max)
  {
    cudimotOptions& opts = cudimotOptions::getInstance();
    if(opts.gridSearch.value()!=""){
      nGridParams=nGP;
      gridParams_host=new int[NPARAMS];
      for(int i=0;i<NPARAMS;i++) gridParams_host[i]=0;
      for(int i=0;i<gP.size();i++) gridParams_host[i]=gP[i];
      cudaMemcpyToSymbol(gridParams,gridParams_host,NPARAMS*sizeof(int));
      gridCombs = gC;
      cudaMalloc((void**)&grid_gpu,gridCombs*nGridParams*sizeof(T));
      cudaMemcpy(grid_gpu,grid_host,gridCombs*nGridParams*sizeof(T),cudaMemcpyHostToDevice);
      
      sync_check("GridSearch: Copying Grid to GPU");

      // Set bounds
      bound_types_host = new int[NPARAMS];
      bounds_min_host = new float[NPARAMS];
      bounds_max_host = new float[NPARAMS];
      for(int p=0;p<NPARAMS;p++){
	bound_types_host[p]=bou_types[p];
	bounds_min_host[p]=bou_min[p];
	bounds_max_host[p]=bou_max[p];
      }
      cudaMemcpyToSymbol(GSbound_types,bound_types_host,NPARAMS*sizeof(int));
      cudaMemcpyToSymbol(GSbounds_min,bounds_min_host,NPARAMS*sizeof(float));
      cudaMemcpyToSymbol(GSbounds_max,bounds_max_host,NPARAMS*sizeof(float));
      sync_check("GridSearch: Setting Bounds");
   
      DEBUG=false;
      if(opts.debug.set()){
	DEBUG=true;
	debugVOX= atoi(opts.debug.value().data());
      }
    }
  }
  
  template <typename T>
  void GridSearch<T>::run(int nvox, int nmeas,
			  int CFP_size, int FixP_size,
			  T* meas, T* params,
			  T* CFP, T* FixP) 
  {
  
    long int amount_shared_mem = 0;
    amount_shared_mem += 2*VOXELS_BLOCK*sizeof(double); // cost function
    amount_shared_mem += (nmeas*CFP_size)*sizeof(T); // CFP
    amount_shared_mem += (NPARAMS*VOXELS_BLOCK)*sizeof(T); // Parameters
    amount_shared_mem += (NPARAMS*VOXELS_BLOCK)*sizeof(T); // trial_Params
    amount_shared_mem += VOXELS_BLOCK*sizeof(int); // boundsTest
        
    cout << "Shared Memory used in GridSearch kernel: " << amount_shared_mem << endl;
    
    int threads_block = VOXELS_BLOCK * THREADS_VOXEL;
    int nblocks=(nvox/VOXELS_BLOCK);
    if(nvox%VOXELS_BLOCK) nblocks++;
    
    if(!DEBUG){
      gridSearch_kernel<T,false><<<nblocks,threads_block,amount_shared_mem>>>(nGridParams,gridCombs,nmeas,CFP_size,FixP_size,meas,grid_gpu,params,CFP,FixP,debugVOX);
    }else{
      gridSearch_kernel<T,true><<<nblocks,threads_block,amount_shared_mem>>>(nGridParams,gridCombs,nmeas,CFP_size,FixP_size,meas,grid_gpu,params,CFP,FixP,debugVOX);
    }
    sync_check("GridSearch Kernel");
  }
  
  // Explicit Instantiations of the template
  template class GridSearch<float>;
  template class GridSearch<double>;
}

