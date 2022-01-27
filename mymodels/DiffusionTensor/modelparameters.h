#ifndef MODEL_H_INCLUDED
#define MODEL_H_INCLUDED

///// Edit this /////
typedef double MyType;
#define NPARAMS 7 // S0,Dxx,Dxy,Dxz,Dyy,Dyz,Dzz
#define NCFP 2
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

