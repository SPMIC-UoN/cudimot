// Ball & 2 Sticks Multiexponential:

// S0: 0
// d: 1
// d_std: 2
// f1: 3
// th1: 4 
// ph1: 5 

// CFP[0:2] are bvecs 
// CFP[3] are bvals

MACRO T Predicted_Signal(
			 int npar, // Number of Parameters to estimate
			 T* P, // Estimated parameters
			 T* CFP, // Fixed Parameters common to all the voxels
			 T* FixP) // Fixed Parameters for each voxel
{
  T isoterm;
  isoterm = exp_gpu( ((P[1]*P[1])/(P[2]*P[2])) * log_gpu(P[1]/(P[1]+P[2]*P[2]*CFP[3])) );
   // exp( d*d/d_std*d_std * log(d/(d+d_std*d_std*bval)) )
   
  T xv_1 = CFP[0]*sin_gpu(P[4])*cos_gpu(P[5])	// (bvec(1)*sinth1*cosph1
	+ CFP[1]*sin_gpu(P[4])*sin_gpu(P[5]) 	// + bvec(2)*sinth1*sinph1
	+ CFP[2]*cos_gpu(P[4]);			// + bvec(3)*costh1)

  T anisoterm_1;
  anisoterm_1 = exp_gpu( ((P[1]*P[1])/(P[2]*P[2])) * log_gpu(P[1]/(P[1]+P[2]*P[2]*CFP[3]*xv_1*xv_1)) );	
  // exp( d*d/d_std*d_std * log(d/(d+d_std*d_std*bval*xv_1*xv_1)) )
  
  T pred_signal = P[0]* (((T)1.0-P[3])*isoterm + P[3]*anisoterm_1);
  // S0*((1-f1-f2)*isoterm + f1*anisoterm_1 + f2*anisoterm_2)
 
  return pred_signal;
}

// Constraints checked during MCMC (if MCMC is used)
MACRO bool ConstraintsMCMC(
		       int npar, // Number of Parameters to estimate
		       T* P) // Estimated parameters
{
  return true;
}

// Partial derivatives respect each model parameter
MACRO void Partial_Derivatives(
			       int npar, // Number of Parameters to estimate
			       T* P, // Estimated parametersß
			       T* CFP, // Fixed Parameters common to all the voxels
			       T* FixP, // Fixed Parameters for each voxel
			       T* derivatives) // Derivative respect each model estimated parameter
{
  T isoterm = exp_gpu( ((P[1]*P[1])/(P[2]*P[2])) * log_gpu(P[1]/(P[1]+P[2]*P[2]*CFP[3])) );
  // exp( d*d/d_std*d_std * log(d/(d+d_std*d_std*bval)) )
  T denlogtermIso = (P[1] + P[2]*P[2]*CFP[3]);
  // d + dstd * dstd * bval

  T xv_1 = CFP[0]*sin_gpu(P[4])*cos_gpu(P[5])	// (bvec(1)*sinth1*cosph1
	+ CFP[1]*sin_gpu(P[4])*sin_gpu(P[5]) 	// + bvec(2)*sinth1*sinph1
	+ CFP[2]*cos_gpu(P[4]);			// + bvec(3)*costh1)
		
  T anisoterm_1 = exp_gpu( ((P[1]*P[1])/(P[2]*P[2])) * log_gpu(P[1]/(P[1]+P[2]*P[2]*CFP[3]*xv_1*xv_1)) );
  // exp( d*d/d_std*d_std * log(d/(d+d_std*d_std*bval*xv_1*xv_1)) )
  T denlogtermAniso1 = (P[1] + P[2]*P[2]*CFP[3]*xv_1*xv_1);
  // d + dstd * dstd * bval * xv_1*xv_1


  // d/dS0
  derivatives[0]= (((T)1.0-P[3])*isoterm + P[3]*anisoterm_1); 
  // ((1-f1-f2)*isoterm + f1*anisoterm_1 + f2*anisoterm_2)				

  // d/dd
  // from isoterm
  derivatives[1]= ((T)1.0-P[3]) * (isoterm *( (((T)2.0*P[1])/(P[1]*P[1]))*log_gpu(P[1]/denlogtermIso) + ((P[1]*P[1])/(P[2]*P[2]))*(denlogtermIso/P[1])*((denlogtermIso-P[1])/(denlogtermIso*denlogtermIso)) ) );  
  // from anisoterm1
  derivatives[1]+= P[3] * (anisoterm_1 *( (((T)2.0*P[1])/(P[1]*P[1]))*log_gpu(P[1]/denlogtermAniso1) + ((P[1]*P[1])/(P[2]*P[2]))*(denlogtermAniso1/P[1])*((denlogtermAniso1-P[1])/(denlogtermAniso1*denlogtermAniso1)) ) );    
  // * S0
  derivatives[1]*=P[0];

  // d/dd_std
  // from isoterm
  derivatives[2]= ((T)1.0-P[3]) * (isoterm *( (((T)-2.0*P[1]*P[1])/(P[1]*P[1]*P[1]))*log_gpu(P[1]/denlogtermIso) + ((P[1]*P[1])/(P[2]*P[2]))*(denlogtermIso/P[1])*(((T)-2.0*P[2]*CFP[3]*P[1])/(denlogtermIso*denlogtermIso)) ) ); 
  // from anisoterm1
  derivatives[2]+= P[3] * (anisoterm_1 *( (((T)-2.0*P[1]*P[1])/(P[1]*P[1]*P[1]))*log_gpu(P[1]/denlogtermAniso1) + ((P[1]*P[1])/(P[2]*P[2]))*(denlogtermAniso1/P[1])*(((T)-2.0*P[2]*CFP[3]*xv_1*xv_1*P[1])/(denlogtermAniso1*denlogtermAniso1)) ) ); 
  // * S0
  derivatives[2]*=P[0];

  // d/df1
  derivatives[3]= P[0]*(-isoterm + anisoterm_1);
  //S0*(-isoterm + anisoterm_1) 

  // d/dth1
  //d xv_1/dth1
  T xv_1_d= cos_gpu(P[4])*(cos_gpu(P[5])*CFP[0] + sin_gpu(P[5])*CFP[1]) -sin_gpu(P[4])*CFP[2];
  derivatives[4]= P[0]*P[3]*anisoterm_1*((P[1]*P[1])/(P[2]*P[2]))*(denlogtermAniso1/P[1])* 
    (((T)-2.0*P[1]*P[2]*P[2]*CFP[3]*xv_1*xv_1_d)/(denlogtermAniso1*denlogtermAniso1));

  // d/dph1
  //d xv_1/dph1
  xv_1_d= sin_gpu(P[4])*(-sin_gpu(P[5])*CFP[0] + cos_gpu(P[5])*CFP[1]);
  derivatives[5]= P[0]*P[3]*anisoterm_1*((P[1]*P[1])/(P[2]*P[2]))*(denlogtermAniso1/P[1])* 
    (((T)-2.0*P[1]*P[2]*P[2]*CFP[3]*xv_1*xv_1_d)/(denlogtermAniso1*denlogtermAniso1));
}

// Constraints run after LevenbergMarquardt (if LevenbergMarquardt is used)
MACRO void FixConstraintsLM(	
		       int npar, // Number of Parameters to estimate
		       T* P) // Estimated parameters
{
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

