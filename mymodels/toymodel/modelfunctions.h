// 5 functions to be implemented:
// - The model predicted signal
// - Constraints during MCMC (Optional)
// - Partial derivatives for Levenberg-Marquardt (Optional)
// - Constraints after Levenberg-Marquardt (Optional)
// - custom priors function (Optional)

// P[0]: a
// P[1]: b
// CFP[0]: Xs
// f(x) = a * exp(-b*x)

MACRO T Predicted_Signal(
			 int npar, 	// Number of Parameters to estimate
			 T* P, 		// Estimated parameters
			 T* CFP, 	// Fixed Parameters common to all the voxels
			 T* FixP) 	// Fixed Parameters for each voxel
{
  return P[0]*exp_gpu(-P[1]*CFP[0]);
}

// Constraints checked during MCMC (if MCMC is used)
MACRO bool ConstraintsMCMC(
			   int npar, // Number of Parameters to estimate
			   T* P) // Estimated parameters
{
  // You can add any constraints (but remember that you can also specify bounds in a different file)
  // These are useful for specify relation between parameters. For instance:
  // if (P[3]>(P[4]+3.14)) return false;

  // Do not modify this.
  return true;
}

// Partial derivatives respect each model parameter
// If Levenberg–Marquardt algorithm is used, the value of the partial derivative for each parameter has to be stored in the outpu array derivatives
// You can use Numerical differentiation using the keyword "NUMERICAL"
MACRO void Partial_Derivatives(
			       int npar, // Number of Parameters to estimate
			       T* P, // Estimated parameters, use P*
			       T* CFP, // Fixed Parameters common to all the voxels
			       T* FixP, // Fixed Parameters for each voxel
			       T* derivatives) // Derivative respect each model estimated parameter
{
  derivatives[0]=exp_gpu(-P[1]*CFP[0]);
  derivatives[1]=-P[0]*CFP[0]*exp_gpu(-P[1]*CFP[0]);
  //derivatives[0]=NUMERICAL(0);
  //derivatives[1]=NUMERICAL(1);

}

// Constraints run after LevenbergMarquardt (if Levenberg-Marquardt is used)
MACRO void FixConstraintsLM(	
			    int npar, // Number of Parameters to estimate
			    T* P) // Estimated parameters
{
  // You can specify what to do with some parameters after Levenberg-Marquardt if a constraint is not satisfied.
  // For instance:
  // if(P[2]>1.0]) P[2]=1.0; 
}



// If needed, here you can define custom priors
// this function will be called only for the parameters that have a prior custom() in modelpriors file
// It must return a single value
MACRO T custom_priors(
       int id_p,   // the number of parameter in the model (starts at 0)
			 T* P, 		// Estimated parameters
       int nmeas, // Number of measurements per voxel
			 T* CFP, 	// Fixed Parameters common to all the voxels for all measurements !!
			 T* FixP) 	// Fixed Parameters for each voxel
{
  // Here you can define several functions, using as conditional the number of parameter p
  // For instance: if(p==4) do this; if(p==5) do that;    
  // but prior must be defined custom() to call this function !!
  // In this case CFP contains a vector with Number_Common_Parameters X NumberMeasurements
  // So, user needs to deal with the indexes of the vector
  // For support this indexig, the fucntion receive in nmeas the Number of measurements

  // Example, a gaussian prior with mean 3.0 ans std 0.1:
  // return ((P[id_p]-(T)3.0)*(P[id_p]-(T)3.0)/((T)2.0*((T)0.1*(T)0.1)));
  
  return (T)0.0;
}

