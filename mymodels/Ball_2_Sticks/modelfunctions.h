
// In Ball & 2 Sticks model the parameters P are:
// S0 P[0], d P[1], f1 P[2], th1 P[3], ph1 P[4], f2 P[5], th2 P[6], ph2  P[7] 
// CFP[0:2] are bvecs 
// CFP[3] are bvals

MACRO T Predicted_Signal(
			 int npar, 	// Number of Parameters to estimate
			 T* P, 		// Estimated parameters
			 T* CFP, 	// Fixed Parameters common to all the voxels
			 T* FixP) 	// Fixed Parameters for each voxel
{
  T isoterm= exp_gpu(-P[1]*CFP[3]); 		// exp(-d*bval)

  T xv_1 = CFP[0]*sin_gpu(P[3])*cos_gpu(P[4])   // (bvec(1)*sinth1*cosph1
	+ CFP[1]*sin_gpu(P[3])*sin_gpu(P[4]) 	// + bvec(2)*sinth1*sinph1
	+ CFP[2]*cos_gpu(P[3]);			// + bvec(3)*costh1)

  T xv_2 = CFP[0]*sin_gpu(P[6])*cos_gpu(P[7])	// (bvec(1)*sinth2*cosph2
	+ CFP[1]*sin_gpu(P[6])*sin_gpu(P[7]) 	// + bvec(2)*sinth2*sinph2
	+ CFP[2]*cos_gpu(P[6]);			// + bvec(3)*costh2)
		
  T anisoterm_1= exp_gpu(-P[1]*CFP[3]*xv_1*xv_1);	// exp(-d*bval* (pow(xv1,2))
  T anisoterm_2= exp_gpu(-P[1]*CFP[3]*xv_2*xv_2);	// exp(-d*bval* (pow(xv2,2))

  T pred_signal= P[0]* ((1-P[2]-P[5])*isoterm + P[2]*anisoterm_1 + P[5]*anisoterm_2);
  // S0*((1-f1-f2)*isoterm + f1*anisoterm_1 + f2*anisoterm_2)
  return pred_signal;
}

// Constraints checked during MCMC (if MCMC is used)
MACRO bool ConstraintsMCMC(
		       int npar, // Number of Parameters to estimate
		       T* P) 	 // Estimated parameters
{
  if(P[2]<P[5]){
    //if f1<f2 reject sample
    return false;
  }
  return true;
}


// Partial derivatives respect each model parameter
MACRO void Partial_Derivatives(
			       int npar, 	// Number of Parameters to estimate
			       T* P, 		// Estimated parametersß
			       T* CFP, 		// Fixed Parameters common to all the voxels
			       T* FixP, 	// Fixed Parameters for each voxel
			       T* derivatives)  // Derivative respect each model estimated parameter
{
  T isoterm= exp_gpu(-P[1]*CFP[3]); 		// exp(-d*bval)

  T xv_1 = CFP[0]*sin_gpu(P[3])*cos_gpu(P[4])	// (bvec(1)*sinth1*cosph1
	+ CFP[1]*sin_gpu(P[3])*sin_gpu(P[4]) 	// + bvec(2)*sinth1*sinph1
	+ CFP[2]*cos_gpu(P[3]);			// + bvec(3)*costh1)

  T xv_2 = CFP[0]*sin_gpu(P[6])*cos_gpu(P[7])	// (bvec(1)*sinth2*cosph2
	+ CFP[1]*sin_gpu(P[6])*sin_gpu(P[7]) 	// + bvec(2)*sinth2*sinph2
	+ CFP[2]*cos_gpu(P[6]);			// + bvec(3)*costh2)
		
  T anisoterm_1= exp_gpu(-P[1]*CFP[3]*xv_1*xv_1);	// exp(-d*bval* (pow(xv1,2))
  T anisoterm_2= exp_gpu(-P[1]*CFP[3]*xv_2*xv_2);	// exp(-d*bval* (pow(xv2,2))

  // d/dS0
  derivatives[0]= ((1-P[2]-P[5])*isoterm + P[2]*anisoterm_1 + P[5]*anisoterm_2); 
  // ((1-f1-f2)*isoterm + f1*anisoterm_1 + f2*anisoterm_2)				

  // d/dd
  derivatives[1]= -(P[0]*CFP[3])* ((1-P[2]-P[5])*isoterm + P[2]*xv_1*xv_1*anisoterm_1 + P[5]*xv_2*xv_2*anisoterm_2);   
  //-S0*bval* ((1-f1-f2)*isoterm + f1*(pow(xv1,2))*anisoterm_1  f2*(pow(xv2,2))*anisoterm_2)

  // d/df1
  derivatives[2]= P[0]*(-isoterm + anisoterm_1);
  //S0*(-isoterm + anisoterm_1) 

  // d/dth1
  //d xv_1/dth1
  T xv_1_d= cos_gpu(P[3])*(cos_gpu(P[4])*CFP[0] + sin_gpu(P[4])*CFP[1]) -sin_gpu(P[3])*CFP[2];
  derivatives[3]= -2*P[0]*P[1]*P[2]*CFP[3]*xv_1_d*xv_1*anisoterm_1;
  // -2 * S0 * d * f1 * bval * xv_1_d * xv_1 * anisoterm_1

  // d/dph1
  //d xv_1/dph1
  xv_1_d= sin_gpu(P[3])*(-sin_gpu(P[4])*CFP[0] + cos_gpu(P[4])*CFP[1]);
  derivatives[4]= -2*P[0]*P[1]*P[2]*CFP[3]*xv_1_d*xv_1*anisoterm_1;
  // -2 * S0 * d * f1 * bval * xv_1_d * xv_1 * anisoterm_1

   // df/df2
  derivatives[5]= P[0]*(-isoterm + anisoterm_2);
  //S0*(-isoterm + anisoterm_2) 

  // d/dth2
  //d xv_2/dth2
  T xv_2_d= cos_gpu(P[6])*(cos_gpu(P[7])*CFP[0] + sin_gpu(P[7])*CFP[1]) -sin_gpu(P[6])*CFP[2];
  derivatives[6]= -2*P[0]*P[1]*P[5]*CFP[3]*xv_2_d*xv_2*anisoterm_2;
  // -2 * S0 * d * f2 * bval * xv_2_d * xv_2 * anisoterm_2

  // d/dph2
  //d xv_2/dph2
  xv_1_d= sin_gpu(P[6])*(-sin_gpu(P[7])*CFP[0] + cos_gpu(P[7])*CFP[1]);
  derivatives[7]= -2*P[0]*P[1]*P[5]*CFP[3]*xv_2_d*xv_2*anisoterm_2;
  // -2 * S0 * d * f2 * bval * xv_2_d * xv_2 * anisoterm_2
}

// Constraints run after LevenbergMarquardt (if LevenbergMarquardt is used)
MACRO void FixConstraintsLM(	
		       int npar, // Number of Parameters to estimate
		       T* P) 	 // Estimated parameters
{
 if(P[2]<P[5]){	
    // if f1 < f2 switch sticks
    T tmp=P[2];
    P[2]=P[5];
    P[5]=tmp;
    tmp=P[3];
    P[3]=P[6];
    P[6]=tmp;
    tmp=P[4];
    P[4]=P[7];
    P[7]=tmp;
  }
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

