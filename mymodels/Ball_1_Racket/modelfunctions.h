#include "matrices_operations.h"
#include "eigenValues.h"
#include "hyp_Sapprox.h"
#include "hyp_Sapprox_doublePrecision.h"

// Ball & 1 Racket, with only one Diffusion coefficient
// P[0]: d
// P[1]: f aniso
// P[2]: k1
// P[3]: k2
// P[4]: th
// P[5]: ph
// P[6]: psi

// CFP[0:2] are bvecs 
// CFP[3] are bvals

// FixP[0] is S0 

MACRO T Calculate_Anisoterm(T* P, // Estimated parameters
			    T* CFP) // Fixed Parameters common to all the voxels
{
  T denominator = hyp_Sapprox_doublePrecision(-P[2],-P[3], (T)0.0); //hyp_Sapprox(k1,-k2);
 
  T RPsi[9];
  RPsi[0]=cos_gpu(P[6]); 	RPsi[1]=sin_gpu(P[6]); 	RPsi[2]=(T)0.0;
  RPsi[3]=-sin_gpu(P[6]); 	RPsi[4]=cos_gpu(P[6]); 	RPsi[5]=(T)0.0; 
  RPsi[6]=(T)0.0; 		RPsi[7]=(T)0.0;		RPsi[8]=(T)1.0;
  
  T RTheta[9];
  RTheta[0]=cos_gpu(P[4]);	RTheta[1]=(T)0.0;	RTheta[2]=-sin_gpu(P[4]);
  RTheta[3]=(T)0.0;		RTheta[4]=(T)1.0;	RTheta[5]=(T)0.0;
  RTheta[6]=sin_gpu(P[4]);	RTheta[7]=(T)0.0;	RTheta[8]=cos_gpu(P[4]);
  
  T RPhi[9];
  RPhi[0]=cos_gpu(P[5]);	RPhi[1]=sin_gpu(P[5]);	RPhi[2]=(T)0.0;
  RPhi[3]=-sin_gpu(P[5]);	RPhi[4]=cos_gpu(P[5]);	RPhi[5]=(T)0.0;
  RPhi[6]=(T)0.0;		RPhi[7]=(T)0.0;		RPhi[8]=(T)1.0;
  
  T R_tmp[9]; // R_tmp = RPsi * RTheta
  multiply_square_matrices(3,RPsi,RTheta,R_tmp);
  T R[9]; // R = RPsi * RTheta * RPhi
  multiply_square_matrices(3,R_tmp,RPhi,R);
  
  T Bdiag[9]; // Bdiag [-k1 -k2 0]
  Bdiag[0]=-P[2]; 		Bdiag[1]=(T)0.0;	Bdiag[2]=(T)0.0;
  Bdiag[3]=(T)0.0;		Bdiag[4]=-P[3];		Bdiag[5]=(T)0.0;
  Bdiag[6]=(T)0.0;		Bdiag[7]=(T)0.0;	Bdiag[8]=(T)0.0;
  
  T R_t[9]; // R_inv = R'
  transpose_square_matrix(3,R,R_t);
  T B_tmp[9]; // B_tmp = R' * Bdiag
  multiply_square_matrices(3,R_t,Bdiag,B_tmp);
  T B[9]; // B = R' * Bdiag * R 
  multiply_square_matrices(3,B_tmp,R,B);
  
  B[0]=B[0]-CFP[0]*CFP[0]*P[0]*CFP[3];
  B[1]=B[1]-CFP[0]*CFP[1]*P[0]*CFP[3];
  B[2]=B[2]-CFP[0]*CFP[2]*P[0]*CFP[3];
  B[3]=B[3]-CFP[1]*CFP[0]*P[0]*CFP[3];
  B[4]=B[4]-CFP[1]*CFP[1]*P[0]*CFP[3];
  B[5]=B[5]-CFP[1]*CFP[2]*P[0]*CFP[3];
  B[6]=B[6]-CFP[2]*CFP[0]*P[0]*CFP[3];
  B[7]=B[7]-CFP[2]*CFP[1]*P[0]*CFP[3];
  B[8]=B[8]-CFP[2]*CFP[2]*P[0]*CFP[3];
      
  T E[3];
  getEigenValues_symm_3x3(B,E); // A must be symmetric
	
  T numerator = hyp_Sapprox_doublePrecision(E[0],E[1],E[2]);
  
  return (numerator/denominator);
}
MACRO T Predicted_Signal(
			 int npar, // Number of Parameters to estimate
			 T* P, // Estimated parameters
			 T* CFP, // Fixed Parameters common to all the voxels
			 T* FixP) // Fixed Parameters for each voxel
{
  T isoterm = exp_gpu(-P[0]*CFP[3]); // exp(-d*bvals)
  T anisoterm = Calculate_Anisoterm(P,CFP);
  T pred_signal= FixP[0]*(((T)1.0-P[1])*isoterm + P[1]*anisoterm);
 
  return pred_signal;
}

// Constraints checked during MCMC (if MCMC is used)
MACRO bool ConstraintsMCMC(
			   int npar, // Number of Parameters to estimate
			   T* P) // Estimated parameters
{
  if(P[2]<P[3]) return false; // K1 > k2
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
  T isoterm = exp_gpu(-P[0]*CFP[3]); // exp(-d*bvals)
  T anisoterm = Calculate_Anisoterm(P,CFP);
  
   // d
  derivatives[0]= NUMERICAL(0)

  // df/dfaniso
  derivatives[1]= FixP[0]*(-isoterm + anisoterm); // analytic
  //S0*(-isoterm + anisoterm) 

  // k1,k2,th,ph,psi
  // numerical differentiation
  derivatives[2]= NUMERICAL(2); 
  derivatives[3]= NUMERICAL(3);
  derivatives[4]= NUMERICAL(4);
  derivatives[5]= NUMERICAL(5);
  derivatives[6]= NUMERICAL(6);
}

// Constraints run after LevenbergMarquardt (if LevenbergMarquardt is used)
MACRO void FixConstraintsLM(	
			    int npar, // Number of Parameters to estimate
			    T* P) // Estimated parameters
{
  if(P[2]<P[3]){
    // k1 < k2
    T temp = P[2];
    P[2]=P[3];
    P[3]=temp;
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

