#ifndef MODEL_H_INCLUDED
#define MODEL_H_INCLUDED

///// Edit this /////
typedef double MyType;
#define NPARAMS 14 
// fiso, 
// fFan2 (out of all Fans)
// fintra1, kappa1, beta1, th1, ph1, psi1
// fintra2, kappa2, beta2, th2, ph2, psi2
#define NCFP 2 // bvecs, bvals
#define NFIXP 1 // S0
/////////////////////

///// Do not edit this /////
struct MODEL
{
  static int CFP_size[NCFP];
  static int FixP_size[NFIXP];
};
////////////////////////////

#endif

