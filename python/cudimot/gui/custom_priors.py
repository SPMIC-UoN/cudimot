"""
CUDIMOT GUI: Custom priors code
"""
from .code import UserCode

HELP="""\
With this function you can define your own priors for any parameter. If you define the prior as
custom() (see next section), this function will be called. You can use voxelwise information 
using FixP (Fixed Parameters).

You can write a different custom prior for each model parameter, but you must use the id_p 
(Id of parameter) argument and if-else-if conditions to differentiate between actions to take,
since any prior defined as custom() will call this function:
"""

BOILERPLATE = """\
MACRO T custom_priors(
       int id_p,   // the number of parameter in the model (starts at 0)
			 T* P, 		// Estimated parameters
       int nmeas, // Number of measurements per voxel
			 T* CFP, 	// Fixed Parameters common to all the voxels for all measurements !!
			 T* FixP) 	// Fixed Parameters for each voxel
{
}
"""

class CustomPriors(UserCode):
    """
    Tab page containing entry for custom priors code
    """

    def __init__(self, app, parent, idx, n):
        UserCode.__init__(self, app, parent, "Custom priors", idx, n, name="custom_priors", function_names=["custom_priors"], help=HELP, boilerplate=BOILERPLATE)
