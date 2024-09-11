/* testFunctions.cu

   A method to evaluate the model predicted signal and the derivatives given all the parameters (fixed and non-fixed). It evaluates the functions for just one set of parameters (does no accept volumes)

   Moises Hernandez-Fernandez - FMRIB Image Analysis Group
   
   Copyright (C) 2005 University of Oxford */

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
   printf("-------------------------------\n");
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
