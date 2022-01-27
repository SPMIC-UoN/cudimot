//http://docs.nvidia.com/cuda/cuda-math-api/group__CUDA__MATH__SINGLE.html#group__CUDA__MATH__SINGLE
//http://docs.nvidia.com/cuda/cuda-c-programming-guide/#mathematical-functions-appendix

#define FUNC __device__ inline

FUNC float log_gpu(float x){return logf(x);}
FUNC double log_gpu(double x){return log(x);}
FUNC float exp_gpu(float x){return expf(x);}
FUNC double exp_gpu(double x){return exp(x);}
FUNC float fabs_gpu(float x){return fabsf(x);}
FUNC double fabs_gpu(double x){return fabs(x);}
FUNC float sin_gpu(float x){return sinf(x);}
FUNC double sin_gpu(double x){return sin(x);}
FUNC float cos_gpu(float x){return cosf(x);}
FUNC double cos_gpu(double x){return cos(x);}
FUNC float atan_gpu(float x){return atanf(x);}
FUNC double atan_gpu(double x){return atan(x);}
FUNC float min_gpu(float x1, float x2){return fminf(x1,x2);}
FUNC double min_gpu(double x1, double x2){return fmin(x1,x2);}
FUNC float max_gpu(float x1, float x2){return fmaxf(x1,x2);}
FUNC double max_gpu(double x1, double x2){return fmax(x1,x2);}

FUNC double tan_gpu(float x){return tan(x);}
FUNC double tan_gpu(double x){return tan(x);}
// tan() always in double precision (single 4 range error)
FUNC double sqrt_gpu(float x){return sqrt(x);}
FUNC double sqrt_gpu(double x){return sqrt(x);}
//sqrt() always in double precision (single 3 range error)

FUNC double atan2_gpu(float x,float y){return atan2(x,y);}
FUNC double atan2_gpu(double x,double y){return atan2(x,y);}
//atan2() always in double precision (single 3 range error)

FUNC double acos_gpu(float x){return acos(x);}
FUNC double acos_gpu(double x){return acos(x);}
//acos() always in double precision (single 3 range error)

FUNC double pow_gpu(float x1, float x2){return pow(x1,x2);}
FUNC double pow_gpu(double x1, double x2){return pow(x1,x2);}
//pow() always in double precision (single 8 range error)



// shfl_down function for double precision
FUNC double shfl_down(double x, uint s){
  int lo, hi;
  asm volatile( "mov.b64 {%0,%1}, %2;" : "=r"(lo), "=r"(hi) : "d"(x) );
  lo = __shfl_down(lo,s);
  hi = __shfl_down(hi,s);
  asm volatile( "mov.b64 %0, {%1,%2};" : "=d"(x) : "r"(lo), "r"(hi) );
  return x;
}

// shfl_down function for single precision
FUNC float shfl_down(float x, uint s){
  return __shfl_down(x,s);
} 

// shfl function for double precision
FUNC double shfl(double x, int s) {
  int lo, hi;
  asm volatile( "mov.b64 {%0,%1}, %2;" : "=r"(lo), "=r"(hi) : "d"(x) );
  lo = __shfl(lo,s);
  hi = __shfl(hi,s);
  asm volatile( "mov.b64 %0, {%1,%2};" : "=d"(x) : "r"(lo), "r"(hi) );
  return x;
}

// shfl function for single precision
FUNC float shfl(float x, uint s){
  return __shfl(x,s);
} 

// shfl function for int
FUNC int shfl(int x, int s){
  return __shfl(x,s);
} 
