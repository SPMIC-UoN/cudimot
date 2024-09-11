/*  init_gpu.cu
    
    Moises Hernandez-Fernandez FMRIB Image Analysis Group
    
    Copyright (C) 2005 University of Oxford  */

#include "checkcudacalls.h"
#include <fstream>

void init_gpu(){
  int *q;
  cudaMalloc((void **)&q, sizeof(int));
  cudaFree(q);
  sync_check("init_gpu");
  
  int device;
  cudaGetDevice(&device);
  printf ("\n...................Allocated GPU %d...................\n", device);
  cudaDeviceSetCacheConfig(cudaFuncCachePreferShared);
  sync_check("init_gpu");
} 

