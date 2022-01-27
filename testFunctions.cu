/* testFunctions.cu

   A method to evaluate the model predicted signal and the derivatives given all the parameters (fixed and non-fixed). It evaluates the functions for just one set of parameters (does no accept volumes)

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


#include <iostream>
#include "checkcudacalls.h"
#include "functions_gpu.h"
#include "modelparameters.h"
#include "macro_numerical.h"
#include "modelfunctions.h"

using namespace std;

template <typename T>
__global__ void testFunctions_kernel(T* params,
				     T* CFP,
				     T* FixP,
				     int CFP_Tsize,
				     int FixP_Tsize)
 {
   T pred;
   pred=Predicted_Signal(NPARAMS,params,CFP,FixP);

   T myderivatives[NPARAMS];
   Partial_Derivatives(NPARAMS,params,CFP,FixP,myderivatives);
  
   for(int i=0;i<NPARAMS;i++){
     printf("Parameter[%i]: %f\n",i,params[i]);
   }
   printf("\n");
   for(int i=0;i<CFP_Tsize;i++){
     printf("CFP[%i]: %f\n",i,CFP[i]);
   }
   printf("\n");
   for(int i=0;i<FixP_Tsize;i++){
     printf("FixP[%i]: %f\n",i,FixP[i]);
   }
   printf("\n");

   printf("Predicted_Signal: %f\n\n",pred);
   for(int i=0; i<NPARAMS; i++){
     printf("Derivative[%i]: %f\n",i,myderivatives[i]);
   }
   printf("-------------------------------\n",pred);
   for(int i=0; i<NPARAMS; i++){
     printf("Numerical Derivative[%i]: %f\n",i,numerical(i,params,CFP,FixP));
   }
 }



int main(int argc, char** argv){

  int total_args = NPARAMS;
  int CFP_Tsize=0;
  int FixP_Tsize=0;
  for(int i=0;i<NCFP;i++){
    total_args+=MODEL::CFP_size[i];
    CFP_Tsize+=MODEL::CFP_size[i];
  }
  for(int i=0;i<NFIXP;i++){
    total_args+=MODEL::FixP_size[i];
    FixP_Tsize+=MODEL::FixP_size[i];
  }

  if(argc<=total_args){
    cout << "CUDIMOT" << endl;
    cout << "Usage:"<< endl;
    cout << "\t" << argv[0] << " [Model_Parameters] [CFP] [FIXP]" << endl;
    cout << endl << "This model has:" << endl;
    cout << "\t" << NPARAMS << " Model_Parameters" << endl;
    cout << "\t" << NCFP << " CFP";
    for(int i=0;i<NCFP;i++){
      cout << "[" << MODEL::CFP_size[i] << "]";
    }
    cout << endl;
    cout << "\t" << NFIXP << " FIXP";
    for(int i=0;i<NFIXP;i++){
      cout << "[" << MODEL::FixP_size[i] << "]";
    }
    cout << endl << endl;
    exit(-1);
  }
  
  cout << "--- Test Model Functions ---" <<endl;
  
  MyType* params = new MyType[NPARAMS];
  MyType* CFP = new MyType[CFP_Tsize];
  MyType* FixP = new MyType[FixP_Tsize];
  
  int par=1;
  for(int i=0;i<NPARAMS;i++){
    params[i]=atof(argv[par]);
    par++;
  }
  for(int i=0;i<CFP_Tsize;i++){
    CFP[i]=atof(argv[par]);
    par++;
  }
  for(int i=0;i<FixP_Tsize;i++){
    FixP[i]=atof(argv[par]);
    par++;
  }

  MyType* params_gpu;
  MyType* CFP_gpu;
  MyType* FixP_gpu;
  cudaMalloc((void**)&params_gpu,NPARAMS*sizeof(MyType));
  cudaMalloc((void**)&CFP_gpu,CFP_Tsize*sizeof(MyType));
  cudaMalloc((void**)&FixP_gpu,FixP_Tsize*sizeof(MyType));

  cudaMemcpy(params_gpu,params,NPARAMS*sizeof(MyType),cudaMemcpyHostToDevice);
  cudaMemcpy(CFP_gpu,CFP,CFP_Tsize*sizeof(MyType),cudaMemcpyHostToDevice);
  cudaMemcpy(FixP_gpu,FixP,FixP_Tsize*sizeof(MyType),cudaMemcpyHostToDevice);
 
  sync_check("Copying Parameters to GPU\n");

  testFunctions_kernel<<<1,1>>>(params_gpu,CFP_gpu,FixP_gpu,CFP_Tsize,FixP_Tsize);
  sync_check("test_Functions Kernel\n");

  return 1;

}
