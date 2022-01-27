#ifndef MODEL_H_INCLUDED
#define MODEL_H_INCLUDED

///// Edit this /////
typedef double MyType;
#define NPARAMS 5	// S0, d, f, th ph
#define NCFP 2		// bvecs, bvals
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

