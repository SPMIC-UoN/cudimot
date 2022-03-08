/*******************************************************************
 * Autogenerated by cudimot_gui
 *
 * Any changes made to this file will be overwritten if it is opened
 * in cudimot_gui!
 *******************************************************************/
// Parameters in 
// P[0]: fiso
// P[1]: fintra
// P[2]: kappa
// P[3]: beta
// P[4]: th
// P[5]: ph
// P[6]: psi

// CFP[0:2] are bvecs 
// CFP[3] are bvals

// FixP[0] is S0

MACRO T Predicted_Signal(
                         int npar, // Number of Parameters to estimate
                         T* P,  // Estimated parameters
                         T* CFP, // Fixed Parameters common to all the voxels
                         T* FixP) // Fixed Parameters for each voxel
{
    // Add code to calculate the model predicted signal here
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
}


// Constraints run after LevenbergMarquardt (if LevenbergMarquardt is used)
MACRO void FixConstraintsLM(	
	int npar,        // Number of Parameters to estimate
	T* P)            // Estimated parameters
{
}


MACRO T custom_priors(
       int id_p,   // the number of parameter in the model (starts at 0)
			 T* P, 		// Estimated parameters
       int nmeas, // Number of measurements per voxel
			 T* CFP, 	// Fixed Parameters common to all the voxels for all measurements !!
			 T* FixP) 	// Fixed Parameters for each voxel
{
}
