#define MACRO template <typename T> __device__ inline

// In DTI model the parameters P are:
// S0: 0
// Dxx: 1
// Dxy: 2
// Dxz: 3 
// Dyy: 4 
// Dyz: 5
// Dzz: 6
// CFP[0:2] are bvecs 
// CFP[3] are bvals

#define EXAMPLE_CONSTANT 1.0

MACRO T Example_subroutine(T input){
  return input+0.0;
}

MACRO T Predicted_Signal(
			 int npar, // Number of Parameters to estimate
			 T* P, // Estimated parameters
			 T* CFP, // Fixed Parameters common to all the voxels
			 T* FixP) // Fixed Parameters for each voxel
{
  T pred_signal= P[0]*exp(-CFP[3]*pow_gpu(CFP[0],2)*P[1] 
			  - CFP[3]*pow_gpu(CFP[1],2)*P[4] 
			  - CFP[3]*powf(CFP[2],2)*P[6] 
			  - 2*CFP[3]*CFP[0]*CFP[1]*P[2] 
			  - 2*CFP[3]*CFP[0]*CFP[2]*P[3] 
			  - 2*CFP[3]*CFP[1]*CFP[2]*P[5]);
  
  return Example_subroutine<T>(pred_signal)*(T)EXAMPLE_CONSTANT;
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
  // df/dS0
  derivatives[0]= exp(-CFP[3]*pow_gpu(CFP[0],2)*P[1] 
		      - CFP[3]*pow_gpu(CFP[1],2)*P[4] 
		      - CFP[3]*pow_gpu(CFP[2],2)*P[6] 
		      - 2*CFP[3]*CFP[0]*CFP[1]*P[2] 
		      - 2*CFP[3]*CFP[0]*CFP[2]*P[3] 
		      - 2*CFP[3]*CFP[1]*CFP[2]*P[5]);
  
  // df/dDxx
  derivatives[1]= P[0]*exp(-CFP[3]*pow_gpu(CFP[0],2)*P[1] 
			   - CFP[3]*pow_gpu(CFP[1],2)*P[4] 
			   - CFP[3]*pow_gpu(CFP[2],2)*P[6] 
			   - 2*CFP[3]*CFP[0]*CFP[1]*P[2] 
			   - 2*CFP[3]*CFP[0]*CFP[2]*P[3] 
			   - 2*CFP[3]*CFP[1]*CFP[2]*P[5]) 
    * (-CFP[3]*pow_gpu(CFP[0],2));
  
  
  // df/dDxy
  derivatives[2]= P[0]*exp(-CFP[3]*pow_gpu(CFP[0],2)*P[1] 
			   - CFP[3]*pow_gpu(CFP[1],2)*P[4] 
			   - CFP[3]*pow_gpu(CFP[2],2)*P[6] 
			   - 2*CFP[3]*CFP[0]*CFP[1]*P[2] 
			   - 2*CFP[3]*CFP[0]*CFP[2]*P[3] 
			   - 2*CFP[3]*CFP[1]*CFP[2]*P[5]) 
	* (-2*CFP[3]*CFP[0]*CFP[1]);
  
  // df/dDxz
  derivatives[3]= P[0]*exp(-CFP[3]*pow_gpu(CFP[0],2)*P[1] 
			   - CFP[3]*pow_gpu(CFP[1],2)*P[4] 
			   - CFP[3]*pow_gpu(CFP[2],2)*P[6] 
			   - 2*CFP[3]*CFP[0]*CFP[1]*P[2] 
			   - 2*CFP[3]*CFP[0]*CFP[2]*P[3] 
			   - 2*CFP[3]*CFP[1]*CFP[2]*P[5]) 
    * (-2*CFP[3]*CFP[0]*CFP[2]);
  
  // df/dDyy
  derivatives[4]= P[0]*exp(-CFP[3]*pow_gpu(CFP[0],2)*P[1] 
			   - CFP[3]*pow_gpu(CFP[1],2)*P[4] 
			   - CFP[3]*pow_gpu(CFP[2],2)*P[6] 
			   - 2*CFP[3]*CFP[0]*CFP[1]*P[2] 
			   - 2*CFP[3]*CFP[0]*CFP[2]*P[3] 
			   - 2*CFP[3]*CFP[1]*CFP[2]*P[5]) 
    * (-CFP[3]*pow_gpu(CFP[1],2));
  
  // df/dDyz
  derivatives[5]= P[0]*exp(-CFP[3]*pow_gpu(CFP[0],2)*P[1] 
			   - CFP[3]*pow_gpu(CFP[1],2)*P[4] 
			   - CFP[3]*pow_gpu(CFP[2],2)*P[6] 
			   - 2*CFP[3]*CFP[0]*CFP[1]*P[2] 
			   - 2*CFP[3]*CFP[0]*CFP[2]*P[3] 
			   - 2*CFP[3]*CFP[1]*CFP[2]*P[5]) 
    * (-2*CFP[3]*CFP[1]*CFP[2]);
  
  // df/dDzz
  derivatives[6]= P[0]*exp(-CFP[3]*pow_gpu(CFP[0],2)*P[1] 
			   - CFP[3]*pow_gpu(CFP[1],2)*P[4] 
			   - CFP[3]*pow_gpu(CFP[2],2)*P[6] 
			   - 2*CFP[3]*CFP[0]*CFP[1]*P[2] 
			   - 2*CFP[3]*CFP[0]*CFP[2]*P[3] 
			   - 2*CFP[3]*CFP[1]*CFP[2]*P[5]) 
    * (-CFP[3]*pow_gpu(CFP[2],2));
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
