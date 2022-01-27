#include "diffusivities.h"
#include "hyp_Sapprox.h"
#include "hyp_Sapprox_doublePrecision.h"
#include "eigenValues.h"

// Parameters in NODDI-Bingham with 2 Fans
// P[0]: fiso
// P[1]: fFan2
// P[2]: fintra1
// P[3]: kappa1
// P[4]: beta1
// P[5]: th1
// P[6]: ph1
// P[7]: psi1
// P[8]: fintra2
// P[9]: kappa2
// P[10]: beta2
// P[11]: th2
// P[12]: ph2
// P[13]: psi2

// CFP[0:2] are bvecs 
// CFP[3] are bvals

// FixP[0] is S0

MACRO T Calculate_ExtraCellular1(T* P,
				 T* CFP,
				 T &xv,
				 T &gradXLocal,
				 T &gradYLocal)
{
  // Calculate Signal from ExtraCellular compartment //
  
  //dPerp = dPar*(1-fintra);
  T Dperpendicular = Dparallel*(1-P[2]);

  //BinghamHinderedDiffusionCoeff
  T dPar; // diffusivity along the dominant orientation
  T dPerp1; // diffusivity along the fanning direction
  T dPerp2; // diffusivity along the final direction
  //[dPar dPerp1 dPerp2] = BinghamHinderedDiffusionCoeff(x(1), x(2), x(3), x(4));
  T dParMdPerp = Dparallel - Dperpendicular;
  T delta = (T)0.00005;
  T nc =  hyp_Sapprox_doublePrecision(P[3],P[4],(T)0.0);
  //computed using the derivative relation with finite difference
  T diff =  hyp_Sapprox_doublePrecision(P[3]+delta,P[4],(T)0.0) - hyp_Sapprox_doublePrecision(P[3]-delta,P[4],(T)0.0);
  diff = diff/((T)2.0*delta);
  dPar = Dperpendicular + dParMdPerp*diff/nc;

  diff = hyp_Sapprox_doublePrecision(P[3],P[4]+delta,(T)0.0) - hyp_Sapprox_doublePrecision(P[3],P[4]-delta,(T)0.0);
  diff = diff/((T)2.0*delta);
  dPerp1 = Dperpendicular + dParMdPerp*diff/nc;
 
  //computed using the trace constancy
  dPerp2 = Dparallel + 2*Dperpendicular - dPar - dPerp1;
  ////////

  T cosThetaSq = xv*xv;
  T sinThetaSq = 1-(cosThetaSq);
  T phi= atan2_gpu(gradYLocal, gradXLocal);
  T cosPhi = cos_gpu(phi);
  T cosPhiSq = cosPhi*cosPhi;
  T sinPhiSq = 1-cosPhiSq;
  
  T ExtraCell = exp_gpu(-CFP[3]*(dPar*cosThetaSq + sinThetaSq*(dPerp1*sinPhiSq + dPerp2*cosPhiSq)));
      
  return ExtraCell;
}

MACRO T Calculate_ExtraCellular2(T* P,
				 T* CFP,
				 T &xv,
				 T &gradXLocal,
				 T &gradYLocal)
{
  // Calculate Signal from ExtraCellular compartment //
  
  //dPerp = dPar*(1-fintra);
  T Dperpendicular = Dparallel*(1-P[8]);

  //BinghamHinderedDiffusionCoeff
  T dPar; // diffusivity along the dominant orientation
  T dPerp1; // diffusivity along the fanning direction
  T dPerp2; // diffusivity along the final direction
  //[dPar dPerp1 dPerp2] = BinghamHinderedDiffusionCoeff(x(1), x(2), x(3), x(4));
  T dParMdPerp = Dparallel - Dperpendicular;
  T delta = (T)0.00005;
  T nc =  hyp_Sapprox_doublePrecision(P[9],P[10],(T)0.0);
  //computed using the derivative relation with finite difference
  T diff =  hyp_Sapprox_doublePrecision(P[9]+delta,P[10],(T)0.0) - hyp_Sapprox_doublePrecision(P[9]-delta,P[10],(T)0.0);
  diff = diff/((T)2.0*delta);
  dPar = Dperpendicular + dParMdPerp*diff/nc;

  diff = hyp_Sapprox_doublePrecision(P[9],P[10]+delta,(T)0.0) - hyp_Sapprox_doublePrecision(P[9],P[10]-delta,(T)0.0);
  diff = diff/((T)2.0*delta);
  dPerp1 = Dperpendicular + dParMdPerp*diff/nc;
 
  //computed using the trace constancy
  dPerp2 = Dparallel + 2*Dperpendicular - dPar - dPerp1;
  ////////

  T cosThetaSq = xv*xv;
  T sinThetaSq = 1-(cosThetaSq);
  T phi= atan2_gpu(gradYLocal, gradXLocal);
  T cosPhi = cos_gpu(phi);
  T cosPhiSq = cosPhi*cosPhi;
  T sinPhiSq = 1-cosPhiSq;
  
  T ExtraCell = exp_gpu(-CFP[3]*(dPar*cosThetaSq + sinThetaSq*(dPerp1*sinPhiSq + dPerp2*cosPhiSq)));
      
  return ExtraCell;
}


MACRO T Calculate_IntraCellular1(T* P,
				 T* CFP,
				 T& xv,
				 T &gradXLocal,
				 T &gradYLocal)
{
  // Calculate Signal from IntraCellular compartment //

  T parComp =-CFP[3]*Dparallel; // Parallel component: -bval * dintra
  // Radius is 0, so no Perpendicular Component

  T BinghamNC = 1/hyp_Sapprox_doublePrecision(P[3],P[4],(T)0.0);

  T Matrix[9];
  // matrix = [0 0 0; 0 beta 0; 0 0 kappa];
  // projG = [grad_dirs(i,:)*normaldir grad_dirs(i,:)*fanningdir grad_dirs(i,:)*fibredir]';
  // matrix = matrix + Lpmp(i)*projG*projG'; 
  Matrix[0] = gradXLocal*gradXLocal*parComp;
  Matrix[1] = gradXLocal*gradYLocal*parComp;
  Matrix[2] = gradXLocal*xv*parComp;
  Matrix[3] = Matrix[1];
  Matrix[4] = P[4] + gradYLocal*gradYLocal*parComp;
  Matrix[5] = gradYLocal*xv*parComp;
  Matrix[6] = Matrix[2];
  Matrix[7] = Matrix[5];
  Matrix[8] = P[3] + xv*xv*parComp;

  T eigenvalues[3];
  getEigenValues_symm_3x3(Matrix,eigenvalues);
  
  T IntraCell = (hyp_Sapprox_doublePrecision(eigenvalues[0],eigenvalues[1],eigenvalues[2])*BinghamNC);
  
  return IntraCell;
}

MACRO T Calculate_IntraCellular2(T* P,
				 T* CFP,
				 T& xv,
				 T &gradXLocal,
				 T &gradYLocal)
{
  // Calculate Signal from IntraCellular compartment //

  T parComp =-CFP[3]*Dparallel; // Parallel component: -bval * dintra
  // Radius is 0, so no Perpendicular Component

  T BinghamNC = 1/hyp_Sapprox_doublePrecision(P[9],P[10],(T)0.0);

  T Matrix[9];
  // matrix = [0 0 0; 0 beta 0; 0 0 kappa];
  // projG = [grad_dirs(i,:)*normaldir grad_dirs(i,:)*fanningdir grad_dirs(i,:)*fibredir]';
  // matrix = matrix + Lpmp(i)*projG*projG'; 
  Matrix[0] = gradXLocal*gradXLocal*parComp;
  Matrix[1] = gradXLocal*gradYLocal*parComp;
  Matrix[2] = gradXLocal*xv*parComp;
  Matrix[3] = Matrix[1];
  Matrix[4] = P[10] + gradYLocal*gradYLocal*parComp;
  Matrix[5] = gradYLocal*xv*parComp;
  Matrix[6] = Matrix[2];
  Matrix[7] = Matrix[5];
  Matrix[8] = P[9] + xv*xv*parComp;

  T eigenvalues[3];
  getEigenValues_symm_3x3(Matrix,eigenvalues);
  
  T IntraCell = (hyp_Sapprox_doublePrecision(eigenvalues[0],eigenvalues[1],eigenvalues[2])*BinghamNC);
  
  return IntraCell;
}

MACRO void get_GradLocals1(
			   T* P,
			   T* CFP,
			   T* fibredir,
			   T* gradsLocal)
{  
  // gradYLocal = grad_dirs*fanningdir;
  // fanningdir = [ -cosPsi.*sinPhi - cosTheta.*cosPhi.*sinPsi cosPhi.*cosPsi - cosTheta.*sinPhi.*sinPsi sinTheta.*sinPsi]';
  T fanningdir[3];
  fanningdir[0] = -cos_gpu(P[7])*sin_gpu(P[6]) - cos_gpu(P[5])*cos(P[6])*sin_gpu(P[7]);
  fanningdir[1] = cos_gpu(P[6])*cos_gpu(P[7]) - cos_gpu(P[5])*sin_gpu(P[6])*sin_gpu(P[7]);
  fanningdir[2] = sin_gpu(P[5])*sin_gpu(P[7]);

  gradsLocal[1] =  CFP[0]*fanningdir[0] +  CFP[1]*fanningdir[1] + CFP[2]*fanningdir[2];
  
  T crossProduct_fanfib[3];
  crossProduct_fanfib[0] = fanningdir[1]*fibredir[2] - fanningdir[2]*fibredir[1];
  crossProduct_fanfib[1] = fanningdir[2]*fibredir[0] - fanningdir[0]*fibredir[2];
  crossProduct_fanfib[2] = fanningdir[0]*fibredir[1] - fanningdir[1]*fibredir[0];
  
  // gradXLocal = grad_dirs*cross(fanningdir,fibredir);
  gradsLocal[0] =  CFP[0]*crossProduct_fanfib[0] +  CFP[1]*crossProduct_fanfib[1] + CFP[2]*crossProduct_fanfib[2];
}

MACRO void get_GradLocals2(
			   T* P,
			   T* CFP,
			   T* fibredir,
			   T* gradsLocal)
{  
  // gradYLocal = grad_dirs*fanningdir;
  // fanningdir = [ -cosPsi.*sinPhi - cosTheta.*cosPhi.*sinPsi cosPhi.*cosPsi - cosTheta.*sinPhi.*sinPsi sinTheta.*sinPsi]';
  T fanningdir[3];
  fanningdir[0] = -cos_gpu(P[13])*sin_gpu(P[12]) - cos_gpu(P[11])*cos(P[12])*sin_gpu(P[13]);
  fanningdir[1] = cos_gpu(P[12])*cos_gpu(P[13]) - cos_gpu(P[11])*sin_gpu(P[12])*sin_gpu(P[13]);
  fanningdir[2] = sin_gpu(P[11])*sin_gpu(P[13]);

  gradsLocal[1] =  CFP[0]*fanningdir[0] +  CFP[1]*fanningdir[1] + CFP[2]*fanningdir[2];
  
  T crossProduct_fanfib[3];
  crossProduct_fanfib[0] = fanningdir[1]*fibredir[2] - fanningdir[2]*fibredir[1];
  crossProduct_fanfib[1] = fanningdir[2]*fibredir[0] - fanningdir[0]*fibredir[2];
  crossProduct_fanfib[2] = fanningdir[0]*fibredir[1] - fanningdir[1]*fibredir[0];
  
  // gradXLocal = grad_dirs*cross(fanningdir,fibredir);
  gradsLocal[0] =  CFP[0]*crossProduct_fanfib[0] +  CFP[1]*crossProduct_fanfib[1] + CFP[2]*crossProduct_fanfib[2];
}

MACRO T Predicted_Signal(
			 int npar, // Number of Parameters to estimate
			 T* P, 	// Estimated parameters
			 T* CFP, // Fixed Parameters common to all the voxels
			 T* FixP) // Fixed Parameters for each voxel
{
  T isoterm = exp_gpu(-CFP[3]*Diso); 	// exp(-bval*d)

  // cosTheta
  T fibredir[3];
  fibredir[0] = sin_gpu(P[5])*cos_gpu(P[6]);
  fibredir[1] = sin_gpu(P[5])*sin_gpu(P[6]);
  fibredir[2] = cos_gpu(P[5]);
  T xv = CFP[0]*fibredir[0] + CFP[1]*fibredir[1] + CFP[2]*fibredir[2];
  T gradsLocal[2];
  get_GradLocals1(P,CFP,fibredir,gradsLocal);
  T ExtraCell = Calculate_ExtraCellular1(P,CFP,xv,gradsLocal[0],gradsLocal[1]);
  T IntraCell = Calculate_IntraCellular1(P,CFP,xv,gradsLocal[0],gradsLocal[1]);
  T Fib1 = (1-P[2])*ExtraCell + P[2]*IntraCell;

  fibredir[0] = sin_gpu(P[11])*cos_gpu(P[12]);
  fibredir[1] = sin_gpu(P[11])*sin_gpu(P[12]);
  fibredir[2] = cos_gpu(P[11]);
  xv = CFP[0]*fibredir[0] + CFP[1]*fibredir[1] + CFP[2]*fibredir[2];
  get_GradLocals2(P,CFP,fibredir,gradsLocal);
  ExtraCell = Calculate_ExtraCellular2(P,CFP,xv,gradsLocal[0],gradsLocal[1]);
  IntraCell = Calculate_IntraCellular2(P,CFP,xv,gradsLocal[0],gradsLocal[1]);
  T Fib2 = (1-P[8])*ExtraCell + P[8]*IntraCell;
    		
  T anisoterm = (1-P[1])*Fib1 + P[1]*Fib2;
      
  T pred_signal=FixP[0]*((1-P[0])*anisoterm + P[0]*isoterm);
    
  return pred_signal;
}

// Constraints checked during MCMC (if MCMC is used)
MACRO bool ConstraintsMCMC(
		       int npar, // Number of Parameters to estimate
		       T* P) // Estimated parameters
{
  if((P[3]-0.1)<P[4]){
    // kappa1 < beta1
    return false;
  }

  if((P[9]-0.1)<P[10]){
    // kappa2 < beta2
    return false;
  }

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
  derivatives[0]=NUMERICAL(0);
  derivatives[1]=NUMERICAL(1);
  derivatives[2]=NUMERICAL(2);
  derivatives[3]=NUMERICAL(3);
  derivatives[4]=NUMERICAL(4);
  derivatives[5]=NUMERICAL(5);
  derivatives[6]=NUMERICAL(6);
  derivatives[7]=NUMERICAL(7);
  derivatives[8]=NUMERICAL(8);
  derivatives[9]=NUMERICAL(9);
  derivatives[10]=NUMERICAL(10);
  derivatives[11]=NUMERICAL(11);
  derivatives[12]=NUMERICAL(12);
  derivatives[13]=NUMERICAL(13);
}

// Constraints run after LevenbergMarquardt (if LevenbergMarquardt is used)
MACRO void FixConstraintsLM(	
			    int npar, // Number of Parameters to estimate
			    T* P) // Estimated parameters
{
  if(P[3]<P[4]){
    // kappa2 < beta2
    P[4]=P[3]-0.5;
  }
  if(P[9]<P[10]){
    // kappa2 < beta2
    P[10]=P[9]-0.5;
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

