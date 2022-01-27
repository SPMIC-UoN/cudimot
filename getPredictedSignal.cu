/* getPredictedSignal.cu

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

// Method to calculate the model predicted signal on the GPU of a group of voxels given the model parameters.

#include "getPredictedSignal.h"
#include "functions_gpu.h"
#include "modelparameters.h"
#include "macro_numerical.h"
#include "modelfunctions.h"

namespace Cudimot{

#define VOXELS_BLOCK 8
#define THREADS_VOXEL 32 // Multiple of 32: Threads collaborating to compute a voxel. Do not change this, otherwise Synchronization will be needed and shuffles cannot be used
  
  template <typename T>
  __global__ void getPredictedSignal_kernel(
					    int nvox, // nvoxels
					    int nmeas, // nmeasurements
					    int nsamples,
					    int CFP_Tsize, //size*M-measurements
					    int FixP_Tsize, // fixed params: size*N-voxels
					    T* samples, // samples of estimated parameters 
					    T* CFP_global, // common fixed model parameters
					    T* FixP, // fixed model parameters
					    T* PredictedSignal)
  {
    // 1 block of threads process several voxels
    // Each warp processes 1 voxel
    int idVOX= (blockIdx.x*VOXELS_BLOCK)+int(threadIdx.x/THREADS_VOXEL);
    int idVOX_inBlock =  threadIdx.x/THREADS_VOXEL;
    int idSubVOX= threadIdx.x%THREADS_VOXEL;
    bool leader = (idSubVOX==0);  // Some steps are performed by only one thread of the warp
    
    ////////// DYNAMIC SHARED MEMORY ///////////
    extern __shared__ double shared[];	      	//Size:
    T* CFP = (T*)shared;		 	//nmeas*CMP_Tsize
    T* meanSamples = (T*) &CFP[nmeas*CFP_Tsize];//NPARAMS*VOXELS_BLOCK
    ////////////////////////////////////////////
    
    /// Copy common fixed model parameters to Shared Memory ///
    if(threadIdx.x==0){ // only one thread of the whole block. Common to all voxels
      for(int i=0;i<nmeas*CFP_Tsize;i++){
	CFP[i]=CFP_global[i];
      }
    }
    __syncthreads();
    ///////////////////////////////////////////////////////////
    
    ///////// each voxel/warp of the block points to its data///////////
    samples = &samples[idVOX*NPARAMS*nsamples]; // Global
    meanSamples = &meanSamples[idVOX_inBlock*NPARAMS];
    PredictedSignal = &PredictedSignal[idVOX*nmeas]; //Global
    FixP = &FixP[idVOX*FixP_Tsize]; // Global memory
    ////////////////////////////////////////////////////////////////////

    /// Ititialise shared values of each voxel: only the leader///
    if(leader){
      if(nsamples>1){
	for(int par=0;par<NPARAMS;par++){
	  T value=0;
	  for(int samp=0;samp<nsamples;samp++){
	    value+= samples[par*nsamples+samp];
	  }
	  meanSamples[par]=value/nsamples;
	}
      }else{
        #pragma unroll
	for(int par=0;par<NPARAMS;par++){
	  meanSamples[par]=samples[par];
	}
      }
    }
    __syncthreads();

    int idMeasurement=idSubVOX;
    int nmeas2compute = nmeas/THREADS_VOXEL;
    if (idSubVOX<(nmeas%THREADS_VOXEL)) nmeas2compute++;
    
    for(int iter=0;iter<nmeas2compute;iter++){
      T* myCFP = &CFP[idMeasurement*CFP_Tsize];
      T pred=Predicted_Signal(NPARAMS,meanSamples,myCFP,FixP);
      PredictedSignal[idMeasurement]=pred;
      idMeasurement+=THREADS_VOXEL;
    }
  }
  
  
  template <typename T>
  getPredictedSignal<T>::getPredictedSignal(){}
  
  template <typename T>
  void getPredictedSignal<T>::run(
				  int nvox, int nmeas, int nsamples,
				  int CFP_size, int FixP_size,
				  T* samples, T* CFP, T* FixP,
				  T* PredictedSignal) 
  {
    
    long int amount_shared_mem = 0;
    amount_shared_mem += (nmeas*CFP_size)*sizeof(T); // CFP
    amount_shared_mem += (NPARAMS*VOXELS_BLOCK)*sizeof(T); //mean_samples
    
    cout << "Shared Memory used in PredictedSignal kernel: " << amount_shared_mem << endl;
    
    int threads_block = VOXELS_BLOCK * THREADS_VOXEL;
    int nblocks=(nvox/VOXELS_BLOCK);
    if(nvox%VOXELS_BLOCK) nblocks++;
    
    getPredictedSignal_kernel<T><<<nblocks,threads_block,amount_shared_mem>>>(nvox,nmeas,nsamples,CFP_size,FixP_size,samples,CFP,FixP,PredictedSignal);
    sync_check("getPredictedSignal Kernel");
  }
  
  // Explicit Instantiations of the template
  template class getPredictedSignal<float>;
  template class getPredictedSignal<double>;
}
