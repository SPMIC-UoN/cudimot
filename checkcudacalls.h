#ifndef CUDIMOT_CHEKCUDACALLS_H_INCLUDED
#define CUDIMOT_CHECKCUDACALLS_H_INCLUDED

/*  checkcudacalls.h
    
    Moises Hernandez-Fernandez FMRIB Image Analysis Group
    
    Copyright (C) 2005 University of Oxford  */

#include <fstream>

#define safecall(call) do{\
    cudaError_t err=call;			\
    if (cudaSuccess != err){						\
      fprintf(stderr,"cuda error at %s:%d. %s\n",__FILE__,__LINE__,cudaGetErrorString(err)); \
    }									\
  }while(0)

#define sync_check(message) do{;		\
    safecall(cudaDeviceSynchronize());		\
    cudaError_t error = cudaGetLastError();	\
    if (cudaSuccess != error){						\
      fprintf(stderr,"ERROR: %s: %s\n",message,cudaGetErrorString(error));	\
      exit(-1);								\
    }									\
  }while(0)

#define checkCuda(error) do{;						\
    if (cudaSuccess != error){						\
      fprintf(stderr,"CUDA Runtime Error: %s\n",cudaGetErrorString(error));	\
      exit(-1);								\
    }									\
  }while(0)

#endif

