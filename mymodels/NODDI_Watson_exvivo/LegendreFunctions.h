#define MACRO template <typename T> __device__ inline

#define MN 7

//legendre gaussian integrals, order 6
MACRO void legendreGaussianIntegral(T x, T* L)
{
  if(x>(T)0.05){
    //exact
    T I[MN];
    T sqrtx = sqrt_gpu(x);
    I[0] = sqrt_gpu((T)MPI)*erff(sqrtx)/sqrtx;
    T dx = (T)1.0/x;
    T emx = -exp_gpu(-x);
    for (int i=1;i<MN;i++){
      I[i] = emx + ((i+1)-(T)1.5)*I[i-1];
      I[i] = I[i]*dx;
    }
    
    L[0] = I[0];
	
    L[1] = -(T)0.5*I[0] + (T)1.5*I[1];
	
    L[2] = (T)0.375*I[0] - (T)3.75*I[1] + (T)4.375*I[2];
   
    L[3] = -(T)0.3125*I[0] + (T)6.5625*I[1] - (T)19.6875*I[2] + (T)14.4375*I[3];
	
    L[4] = (T)0.2734375*I[0] - (T)9.84375*I[1] + (T)54.140625*I[2] - (T)93.84375*I[3] + (T)50.2734375*I[4];
    
    L[5] = -((T)63.0/(T)256.0)*I[0] + ((T)3465.0/(T)256.0)*I[1] - ((T)30030.0/(T)256.0)*I[2] + ((T)90090.0/(T)256.0)*I[3] - ((T)109395.0/(T)256.0)*I[4] + ((T)46189.0/(T)256.0)*I[5];
    
    L[6] = ((T)231.0/(T)1024.0)*I[0] - ((T)18018.0/(T)1024.0)*I[1] + ((T)225225.0/(T)1024.0)*I[2] - ((T)1021020.0/(T)1024.0)*I[3] + ((T)2078505.0/(T)1024.0)*I[4] - ((T)1939938.0/(T)1024.0)*I[5] + ((T)676039.0/(T)1024.0)*I[6];
   
  }else{
    // x<=0.05, approx
    // Computing the legendre gaussian integrals for small x
    T x2=x*x;
    T x3=x2*x;
    T x4=x3*x;
    T x5=x4*x;
    T x6=x5*x;
        
    L[0] = (T)2.0 - (T)2.0*x/(T)3.0 + x2/(T)5.0 - x3/(T)21.0 + x4/(T)108.0;
    
    L[1] = -(T)4.0*x/(T)15.0 + (T)4.0*x2/(T)35.0 - (T)2.0*x3/(T)63.0 + (T)2.0*x4/(T)297.0;
    
    L[2] = (T)8.0*x2/(T)315.0 - (T)8.0*x3/(T)693.0 + (T)4.0*x4/(T)1287.0;
	
    L[3] = -(T)16.0*x3/(T)9009.0 + (T)16.0*x4/(T)19305.0;
    
    L[4] = (T)32.0*x4/(T)328185.0;
    
    L[5] = -(T)64.0*x5/(T)14549535.0;
    
    L[6] = (T)128.0*x6/(T)760543875.0;
  }
}


MACRO T legendre_m0(int l, T x)
{
// From numerical recipes
// Computes the associated Legendre polynomial P m0-l. 
// 0 ≤ m ≤ l

  T pmm=(T)1.0;
  if(l==0){
    return pmm;
  }else{
    T pmmp1=x;
    T pll=(T)0.0;
    if (l==1){
      return pmmp1;
    }else{
      for(int ll=2;ll<=l;ll++){
	pll=(x*((T)2.0*ll-(T)1.0)*pmmp1-(ll-(T)1.0)*pmm)/((T)ll);
	pmm=pmmp1;
	pmmp1=pll;
      }
      return pll;
    }
  }
}
    
MACRO void computeSH_values(T xv, T* SH)
{
  SH[0] = sqrt_gpu(((T)1.0-(T)0.75)/(T)MPI);
  SH[0]*= legendre_m0(0,xv); // l=0
  
  SH[1] = sqrt_gpu(((T)2.0-(T)0.75)/(T)MPI);
  SH[1]*= legendre_m0(2,xv); // l=2

  SH[2] = sqrt_gpu(((T)3.0-(T)0.75)/(T)MPI);
  SH[2]*= legendre_m0(4,xv); // l=4

  SH[3] = sqrt_gpu(((T)4.0-(T)0.75)/(T)MPI);
  SH[3]*= legendre_m0(6,xv); // l=6

  SH[4] = sqrt_gpu(((T)5.0-(T)0.75)/(T)MPI);
  SH[4]*= legendre_m0(8,xv); // l=8

  SH[5] = sqrt_gpu(((T)6.0-(T)0.75)/(T)MPI);
  SH[5]*= legendre_m0(10,xv); // l=10

  SH[6] = sqrt_gpu(((T)7.0-(T)0.75)/(T)MPI);
  SH[6]*= legendre_m0(12,xv); // l=12
}
