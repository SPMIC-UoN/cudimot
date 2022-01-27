#define MACRO template <typename T> __device__ inline
#define NUMERICAL(idpar) numerical(idpar,P,CFP,FixP);
//Used for numerical differentiation
#define TINY 1e-5    
#define MAXSTEP 0.2
#define MINSTEP 1e-6
#define nonzerosign(a) ((a)<0?-1:1)

MACRO T numerical(
		     int idpar, 
		     T* P, 	
		     T* CFP, 
		     T* FixP)
{
  T scale=(T)1.0;
  T step;
  step = (T)TINY*nonzerosign(P[idpar])*fabs_gpu(P[idpar])*scale;
  step = min_gpu(max_gpu(fabs_gpu(step),(T)MINSTEP),(T)MAXSTEP);
  T save_param=P[idpar];
  P[idpar]+=step*(T)0.5;
  T signalA=Predicted_Signal(NPARAMS,P,CFP,FixP);
  P[idpar]-=2*(step*(T)0.5);
  T signalB=Predicted_Signal(NPARAMS,P,CFP,FixP);
  P[idpar]=save_param;
  return (signalA-signalB)/step;
}
