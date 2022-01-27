#define MACRO template <typename T> __device__ inline

#define MPI 3.14159265358979323846
#define M_SQRT_PI 1.772453850905516
#define NMAX 6
#define H 0.4
#define A1 (2.0/3.0)
#define A2 0.4
#define A3 (2.0/7.0)

#include "diffusivities.h"

MACRO T dawson(T x){  // From Numerical Recipes in C
//Returns Dawson’s integral for any real x

  int i,n0;
  T d1,d2,e1,e2,sum,x2,xp,xx,ans;
  T c[NMAX+1];
  int init=0; // Flag is 0 if we need to initialize, else 1
  if (init==0){
    init=1;
    for(i=0;i<NMAX;i++){
      c[i]=((T)2.0*(i+1)-(T)1.0)*(T)H;
      c[i]=exp_gpu(-c[i]*c[i]);
    }
  }
  if(fabs_gpu(x) < (T)0.2){  // Use series expansion.
    x2=x*x;
    ans=x*((T)1.0-(T)A1*x2*((T)1.0-(T)A2*x2*((T)1.0-(T)A3*x2)));
  }else{  // Use sampling theorem representation.
    xx=fabs_gpu(x);
    n0=2*(int)((T)0.5*xx/(T)H+(T)0.5);
    xp=xx-n0*(T)H;
    e1=exp_gpu((T)2.0*xp*(T)H);
    e2=e1*e1;
    d1=n0+1;
    d2=d1-(T)2.0;
    sum=0.0;
    for (i=0;i<NMAX;i++,d1+=(T)2.0,d2-=(T)2.0,e1*=e2){
      sum += c[i]*(e1/d1+(T)1.0/(d2*e1));
    }
    
    //ans=(T)0.5641895835*SIGN(exp_gpu(-xp*xp),x)*sum; 
    ans=exp_gpu(-xp*xp);
    if(ans<0) ans=-ans; // magnitude
    if(x<0) ans=-ans; // * sign(x)
    ans=(T)0.5641895835*ans*sum; //Constant is 1/sqrt(pi)
    
  }
  return ans;
}

MACRO void WatsonHinderedDiffusionCoeff(//input
					T kappa,
					// dPar is a constant
					T dPerp,
					//output
					T& dPar_equivalent,
					T& dPerp_equivalent)
{
  T dParMdPerp = Dparallel - dPerp;
  if (kappa < (T)1e-5){
    T dParP2dPerp = Dparallel + (T)2.0*dPerp;
    T k2 = kappa*kappa;

    dPar_equivalent = dParP2dPerp/(T)3.0 
      + (T)4.0*dParMdPerp*kappa/(T)45.0 
      + (T)8.0*dParMdPerp*k2/(T)945.0;
    
    dPerp_equivalent = dParP2dPerp/(T)3.0 
      - (T)2.0*dParMdPerp*kappa/(T)45.0 
      - (T)4.0*dParMdPerp*k2/(T)945.0;

  }else{
    T sk;
    T dawsonf;
    //if(kapppa<0){
    // never. kappa is always >=0
    //}else{
    sk = sqrt_gpu(kappa);
    //dawsonf= (T)0.5*exp_gpu(-kappa)*sqrt_gpu(M_PI)*erfi_sk;
    dawsonf=dawson(sk);
    
    T factor = sk/dawsonf;
    
    dPar_equivalent = (-dParMdPerp+(T)2.0*dPerp*kappa+dParMdPerp*factor) / ((T)2.0*kappa);
    dPerp_equivalent = (dParMdPerp+(T)2.0*(Dparallel+dPerp)*kappa-dParMdPerp*factor) / ((T)4.0*kappa);
  }
}

MACRO void WatsonSHCoeff(T k, T* coeff)
{
  //commputes the spherical harmonic (SH) coefficients of the Watson's distribution with the concentration parameter k (kappa) and 12th order
  
  coeff[0]= (T)2.0*sqrt_gpu((T)MPI);
  
  if(k>30){
    //very large kappa
    T lnkd = log_gpu(k) - log_gpu((T)30.0);
    T lnkd2 = lnkd*lnkd;
    T lnkd3 = lnkd2*lnkd;
    T lnkd4 = lnkd3*lnkd;
    T lnkd5 = lnkd4*lnkd;
    T lnkd6 = lnkd5*lnkd;
    coeff[1] = (T)7.52308 + (T)0.411538*lnkd - (T)0.214588*lnkd2 + (T)0.0784091*lnkd3 - (T)0.023981*lnkd4 + (T)0.00731537*lnkd5 - (T)0.0026467*lnkd6;
    coeff[2] = (T)8.93718 + (T)1.62147*lnkd - (T)0.733421*lnkd2 + (T)0.191568*lnkd3 - (T)0.0202906*lnkd4 - (T)0.00779095*lnkd5 + (T)0.00574847*lnkd6;
    coeff[3] = (T)8.87905 + (T)3.35689*lnkd - (T)1.15935*lnkd2 + (T)0.0673053*lnkd3 + (T)0.121857*lnkd4 - (T)0.066642*lnkd5 + (T)0.0180215*lnkd6;
    coeff[4] = (T)7.84352 + (T)5.03178*lnkd - (T)1.0193*lnkd2 - (T)0.426362*lnkd3 + (T)0.328816*lnkd4 - (T)0.0688176*lnkd5 - (T)0.0229398*lnkd6;
    coeff[5] = (T)6.30113 + (T)6.09914*lnkd - (T)0.16088*lnkd2 - (T)1.05578*lnkd3 + (T)0.338069*lnkd4 + (T)0.0937157*lnkd5 - (T)0.106935*lnkd6;
    coeff[6] = (T)4.65678 + (T)6.30069*lnkd + (T)1.13754*lnkd2 - (T)1.38393*lnkd3 - (T)0.0134758*lnkd4 + (T)0.331686*lnkd5 - (T)0.105954*lnkd6;

  }else if(k>0.1){
    //exact

    T sk = sqrt_gpu(k);
    T sk2 = sk*k;
    T sk3 = sk2*k;
    T sk4 = sk3*k;
    T sk5 = sk4*k;
    T sk6 = sk5*k;
    //T sk7 = sk6*k;

    T k2 = k*k;
    T k3 = k2*k;
    T k4 = k3*k;
    T k5 = k4*k;
    T k6 = k5*k;
    //T k7 = k6*k;

    T dawsonk = dawson(sk);
    T erfik = ((T)2.0/(T)M_SQRT_PI)*exp_gpu(sk*sk)*dawsonk;
    //erfi(x) = 2/sqrt(pi) * exp(x*x) * dawson(x)
        
    T ierfik = (T)1.0/erfik;
    T ek = exp_gpu(k);
    //T dawsonk = (T)0.5*sqrt_gpu((T)MPI)*erfik/ek;

    coeff[1] = (T)3.0*sk - ((T)3.0 + (T)2.0*k)*dawsonk;
    coeff[1] = sqrt_gpu((T)5.0)*coeff[1]*ek;
    coeff[1] = coeff[1]*ierfik/k;
    
    coeff[2] = ((T)105.0 + (T)60.0*k + (T)12.0*k2)*dawsonk;
    coeff[2] = coeff[2] -(T)105.0*sk + (T)10.0*sk2;
    coeff[2] = (T)0.375*coeff[2]*ek/k2;
    coeff[2] = coeff[2]*ierfik;
    
    coeff[3] = -(T)3465.0 - (T)1890.0*k - (T)420.0*k2 - (T)40.0*k3;
    coeff[3] = coeff[3]*dawsonk;
    coeff[3] = coeff[3] + (T)3465.0*sk - (T)420.0*sk2 + (T)84.0*sk3;
    coeff[3] = coeff[3]*sqrt_gpu((T)13.0*(T)MPI)/(T)64.0/k3;
    coeff[3] = coeff[3]/dawsonk;
    
    coeff[4] = (T)675675.0 + (T)360360.0*k + (T)83160.0*k2 + (T)10080.0*k3 + (T)560.0*k4;
    coeff[4] = coeff[4]*dawsonk;
    coeff[4] = coeff[4] - (T)675675.0*sk + (T)90090.0*sk2 - (T)23100.0*sk3 + (T)744.0*sk4;
    coeff[4] = sqrt((T)17.0)*coeff[4]*ek;
    coeff[4] = coeff[4]/(T)512.0/k4;
    coeff[4] = coeff[4]*ierfik;

    coeff[5] = -(T)43648605.0 - (T)22972950.0*k - (T)5405400.0*k2 - (T)720720.0*k3 - (T)55440.0*k4 - (T)2016.0*k5;
    coeff[5] = coeff[5]*dawsonk;
    coeff[5] = coeff[5] + (T)43648605.0*sk - (T)6126120.0*sk2 + (T)1729728.0*sk3 - (T)82368.0*sk4 + (T)5104.0*sk5;
    coeff[5] = sqrt((T)21.0*(T)MPI)*coeff[5]/(T)4096.0/k5;
    coeff[5] = coeff[5]/dawsonk;
        
    coeff[6] = (T)7027425405.0 + (T)3666482820.0*k + (T)872972100.0*k2 + (T)122522400.0*k3  + (T)10810800.0*k4 + (T)576576.0*k5 + (T)14784.0*k6;
    coeff[6] = coeff[6]*dawsonk;
    coeff[6] = coeff[6] - (T)7027425405.0*sk + (T)1018467450.0*sk2 - (T)302630328.0*sk3 + (T)17153136.0*sk4 - (T)1553552.0*sk5 + (T)25376.0*sk6;
    coeff[6] = (T)5.0*coeff[6]*ek;
    coeff[6] = coeff[6]/(T)16384.0/k6;
    coeff[6] = coeff[6]*ierfik;
    
  }else{
    // K<=0.1 approx
    
    T k2 = k*k;
    T k3 = k2*k;
    T k4 = k3*k;
    T k5 = k4*k;
    T k6 = k5*k;
    //T k7 = k6*k;
    
    coeff[1] = (T)4.0/(T)3.0*k + (T)8.0/(T)63.0*k2;
    coeff[1] = coeff[1]*sqrt_gpu((T)MPI/(T)5.0);
    coeff[2] = (T)8.0/(T)21.0*k2 + (T)32.0/(T)693.0*k3;
    coeff[2] = coeff[2]*(sqrt_gpu((T)MPI)*(T)0.2);
    coeff[3] = (T)16.0/(T)693.0*k3 + (T)32.0/(T)10395.0*k4;
    coeff[3] = coeff[3]*sqrt_gpu((T)MPI/(T)13.0);
    coeff[4] = (T)32.0/(T)19305.0*k4;
    coeff[4] = coeff[4]*sqrt_gpu((T)MPI/(T)17.0);
    coeff[5] = (T)64.0*sqrt_gpu((T)MPI/(T)21.0)*k5/(T)692835.0;
    coeff[6] = (T)128.0*sqrt_gpu((T)MPI)*k6/(T)152108775.0;
  }
}



