#define MACRO template <typename T> __device__ inline

#define gyro 26753.803038


// SSFP !!! 
// In Ball & 2 Sticks model the parameters P are:
// d   P[0]
// f1  P[1]
// th1 P[2]
// ph1 P[3]
// f2  P[4]
// th2 P[5]
// ph2 P[6] 
// S0  P[7]

// CFP[0:2] bvecs 
// CFP[3]   bvals
// CFP[4]   TRs
// CFP[5]   diffGradAmps
// CFP[6]   diffGradDurs
// CFP[7]   flipAngle
// CFP[8]   Noisefloor estimate

// FixP[0] T1
// FixP[1] T2
// FixP[2] B1

MACRO T SSFP_Signal(T adc, T* P, T* CFP, T* FixP, T E1, T E2, T sa, T ca)
{
  T S0 = P[7];
  T qval = (T)gyro*CFP[5]*CFP[6]/(T)10.0;
  // Iso
  T A1    = exp_gpu( -qval*qval*CFP[4]*adc );
  T A2    = exp_gpu( -qval*qval*CFP[6]*adc );
  T A2_03 = exp_gpu( -qval*qval*CFP[6]*adc/(T)3.0 );

  T s     = E2*A1/A2_03/A2_03/A2_03/A2_03*((T)1.0-E1*ca)+E2/A2_03*(ca-E1);
  T r     = (T)1.0 - E1*ca+E2*E2*A1*A2_03*(ca-E1);
  T K     = ((T)1.0-E1*A1*ca-E2*E2*A1*A1/A2_03/A2_03*(E1*A1-ca))/(E2*A1/A2_03/A2_03/A2_03/A2_03*((T)1.0+ca)*((T)1.0-E1*A1));

  T F1    = K - sqrt_gpu(K*K-A2*A2);
  T Mminus_top = -((T)1.0-E1)*E2/A2_03/A2_03*(F1-E2*A1*A2_03*A2_03)*sa;
  T Mminus_bottom=r-F1*s;
  T signal = sqrt_gpu(fabs_gpu(S0*Mminus_top/Mminus_bottom)*fabs_gpu(S0*Mminus_top/Mminus_bottom)+CFP[8]*CFP[8]);


  return signal;
}

MACRO T Predicted_Signal(
			 int npar, 	// Number of Parameters to estimate
			 T* P, 		// Estimated parameters
			 T* CFP, 	// Fixed Parameters common to all the voxels
			 T* FixP) 	// Fixed Parameters for each voxel
{
 
  // Relaxation terms
  T E1 = exp_gpu(-CFP[4]/ ((T)1e-3*FixP[0]));
  T E2 = exp_gpu(-CFP[4]/ ((T)1e-3*FixP[1]));

  // Flipping
  T sa = sin_gpu(CFP[7]*FixP[2]*M_PI/(T)180.0);
  T ca = cos_gpu(CFP[7]*FixP[2]*M_PI/(T)180.0);

  // Diffusion terms  
  T adc_iso= P[0]; 	       
  T xv_1 = CFP[0]*sin_gpu(P[2])*cos_gpu(P[3])   // (bvec(1)*sinth1*cosph1
	+ CFP[1]*sin_gpu(P[2])*sin_gpu(P[3]) 	// + bvec(2)*sinth1*sinph1
	+ CFP[2]*cos_gpu(P[2]);			// + bvec(3)*costh1)

  T xv_2 = CFP[0]*sin_gpu(P[5])*cos_gpu(P[6])	// (bvec(1)*sinth2*cosph2
	 + CFP[1]*sin_gpu(P[5])*sin_gpu(P[6]) 	// + bvec(2)*sinth2*sinph2
	 + CFP[2]*cos_gpu(P[5]);		// + bvec(3)*costh2)		
  T adc_aniso1= P[0]*xv_1*xv_1;	// d*(pow(xv1,2))
  T adc_aniso2= P[0]*xv_2*xv_2;	// d*(pow(xv2,2))

  // SSFP Stuff
  T isoterm=SSFP_Signal(adc_iso,P,CFP,FixP,E1,E2,sa,ca);
  T anisoterm_1=SSFP_Signal(adc_aniso1,P,CFP,FixP,E1,E2,sa,ca);
  T anisoterm_2=SSFP_Signal(adc_aniso2,P,CFP,FixP,E1,E2,sa,ca);

  T pred_signal= ((1-P[1]-P[4])*isoterm + P[1]*anisoterm_1 + P[4]*anisoterm_2);

  return pred_signal;
}

// Constraints checked during MCMC (if MCMC is used)
MACRO bool ConstraintsMCMC(
		       int npar, // Number of Parameters to estimate
		       T* P) 	 // Estimated parameters
{
  if(P[1]<P[4] | (P[1]+P[4])>1){
    //if f1<f2 reject 
    //if f1+f2>1 reject
    return false;
  }
  return true;
}


// Partial derivatives respect each model parameter
MACRO void Partial_Derivatives(
			       int npar, 	// Number of Parameters to estimate
			       T* P, 		// Estimated parameters
			       T* CFP, 		// Fixed Parameters common to all the voxels
			       T* FixP, 	// Fixed Parameters for each voxel
			       T* derivatives)  // Derivative respect each model estimated parameter
{
  
}

// Constraints run after LevenbergMarquardt (if LevenbergMarquardt is used)
MACRO void FixConstraintsLM(	
		       int npar, // Number of Parameters to estimate
		       T* P) 	 // Estimated parameters
{
 // if(P[2]<P[5]){	
 //    // if f1 < f2 switch sticks
 //    T tmp=P[2];
 //    P[2]=P[5];
 //    P[5]=tmp;
 //    tmp=P[3];
 //    P[3]=P[6];
 //    P[6]=tmp;
 //    tmp=P[4];
 //    P[4]=P[7];
 //    P[7]=tmp;
 //  }
}

MACRO T custom_priors(
       int id_p,   // the number of parameter in the model (starts at 0)
			 T* P, 		// Estimated parameters
       int nmeas, // Number of measurements per voxel
			 T* CFP, 	// Fixed Parameters common to all the voxels for all measurements !!
			 T* FixP) 	// Fixed Parameters for each voxel
{
	return (T)0.0;
}
