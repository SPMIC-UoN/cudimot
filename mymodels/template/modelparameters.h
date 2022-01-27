#ifndef MODEL_H_INCLUDED
#define MODEL_H_INCLUDED

///// Edit this /////
typedef double MyType;
// NPARAMS specifies the number of parameters to be estimated, for instance in Ball and 1 Stick we want to estimate 5 parameters: S0, d, f, th ph
#define NPARAMS 5

// NCFP specifies the number of fixed parameters that are common to all the voxels, for instance 2 parameters: bvecs, bvals
#define NCFP 2 

// NFIXP specifies the number of fixed parameters that are different for each voxel (such as S0, T1, ...), in this case 0
#define NFIXP 0
/////////////////////

///// Do not edit this /////
struct MODEL
{
	static int CFP_size[NCFP];
	static int FixP_size[NFIXP];
};
////////////////////////////

#endif

