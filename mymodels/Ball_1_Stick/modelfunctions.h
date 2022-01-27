// In Ball & Sticks (1 stick) model the parameters P are:
// P[0]: S0
// P[1]: d
// P[2]: f
// P[3]: th
// P[4]: ph

// CFP[0:2] are bvecs 
// CFP[3] are bvals

MACRO T Predicted_Signal(
			 int npar, 	// Number of Parameters to estimate
			 T* P, 		// Estimated parameters
			 T* CFP, 	// Fixed Parameters common to all the voxels
			 T* FixP) 	// Fixed Parameters for each voxel
{
  T isoterm= exp_gpu(-P[1]*CFP[3]); 	// exp(-d*bval)

  T xv = CFP[0]*sin_gpu(P[3])*cos_gpu(P[4])	// (bvec(1)*sinth*cosph
	+ CFP[1]*sin_gpu(P[3])*sin_gpu(P[4])   	// + bvec(2)*sinth*sinph
	+ CFP[2]*cos_gpu(P[3]);			// + bvec(3)*costh)
		
  T anisoterm= exp_gpu(-P[1]*CFP[3]*xv*xv);	// exp(-d*bval* (pow(xv,2))

  T pred_signal= P[0]* (((T)1.0-P[2])*isoterm + P[2]*anisoterm);
  // S0*((1-f)*isoterm + f*anisoterm)
 
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
			       T* P, // Estimated parameters, use P*
			       T* CFP, // Fixed Parameters common to all the voxels
			       T* FixP, // Fixed Parameters for each voxel
			       T* derivatives) // Derivative respect each model estimated parameter
{
  T isoterm= exp_gpu(-P[1]*CFP[3]); 	// exp(-d*bval)

  T xv = CFP[0]*sin_gpu(P[3])*cos_gpu(P[4])	// (bvec(1)*sinth*cosph
	+ CFP[1]*sin_gpu(P[3])*sin_gpu(P[4])	// + bvec(2)*sinth*sinph
	+ CFP[2]*cos_gpu(P[3]);			// + bvec(3)*costh)
		
  T anisoterm= exp_gpu(-P[1]*CFP[3]*xv*xv);	// exp(-d*bval* (pow(xv,2))

  // df/dS0
  derivatives[0]= (((T)1.0-P[2])*isoterm + P[2]*anisoterm); // (1-f)*isoterm + f*anisoterm

  // df/dd
  derivatives[1]= -(P[0]*CFP[3])* (((T)1.0-P[2])*isoterm + P[2]*xv*xv*anisoterm);   
  //-S0*bval* ((1-f)*isoterm + (pow(xv,2))*f*anisoterm)

  // df/df
  derivatives[2]= P[0]*(-isoterm + anisoterm);
  //S0*(-isoterm + anisoterm) 

  // df/dth
  T xv1= cos_gpu(P[3])*(cos_gpu(P[4])*CFP[0] + sin_gpu(P[4])*CFP[1]) -sin_gpu(P[3])*CFP[2];
  //d xv/dth
  derivatives[3]= -(T)2.0*P[0]*P[1]*P[2]*CFP[3]*xv*xv1*anisoterm;
  // -2 * S0 * d * f * bval * xv * xv1 * anisoterm

  // df/dph
  xv1= sin_gpu(P[3])*(-sin_gpu(P[4])*CFP[0] + cos_gpu(P[4])*CFP[1]);
  //d xv/dth
  derivatives[4]= -(T)2.0*P[0]*P[1]*P[2]*CFP[3]*xv*xv1*anisoterm;
  // -2 * S0 * d * f * bval * xv * xv1 * anisoterm
}

// Constraints run after LevenbergMarquardt (if Levenberg-Marquardt is used)
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
