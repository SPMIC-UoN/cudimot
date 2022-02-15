"""
CUDIMOT GUI: Forward model code
"""
from .code import UserCode

HELP = """\
Here you can specify the code to calculate your model, i.e. 
given the values of the parameters return the model prediction"""

BOILERPLATE = """\
MACRO T Predicted_Signal(
                         int npar, // Number of Parameters to estimate
                         T* P,  // Estimated parameters
                         T* CFP, // Fixed Parameters common to all the voxels
                         T* FixP) // Fixed Parameters for each voxel
{
    // Add code to calculate the model predicted signal here
}
"""

class ForwardModel(UserCode):
    """
    Tab page containing entry for forward model evaluation code
    """

    def __init__(self, app, parent, idx, n):
        UserCode.__init__(self, app, parent, "Forward model", idx, n, name="fwdmodel", function_names=["Predicted_Signal"], help=HELP, boilerplate=BOILERPLATE)
