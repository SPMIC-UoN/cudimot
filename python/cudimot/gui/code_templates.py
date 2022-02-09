"""
CUDIMOT GUI - Templates for autogenerated code
"""

HEADER = """\
/*******************************************************************
 * Autogenerated by cudimot_gui
 *
 * Any changes made to this file will be overwritten if it is opened
 * in cudimot_gui!
 *******************************************************************/
"""

MODELPARAMETERS_H = """\
#ifndef MODEL_H_INCLUDED
#define MODEL_H_INCLUDED

// Use {precision} precision
typedef {dtype} MyType;

// Define numbers of parameters

// Variable parameters
#define NPARAMS {n_param} // {param_names}

// Common fixed parameters
#define NCFP {n_cfp} // {cfp_names}

// Voxelwise fixed parameters
#define NFIXP {n_vfp} // {vfp_names}

// Model data structure
struct MODEL
{{
  static int CFP_size[NCFP];
  static int FixP_size[NFIXP];
}};
#endif
"""

MODELPARAMETERS_CC = """\
#include "modelparameters.h"

int MODEL::CFP_size[] = {{{cfp_sizes}}};
int MODEL::FixP_size[] = {{{vfp_sizes}}};
"""

MODELFUNCTIONS_H = """\
// Parameters in {name}
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

{code_fwdmodel}

{code_constraints_mcmc}

{code_derivatives_lm}

{code_constraints_lm}

{code_custom_priors}
"""