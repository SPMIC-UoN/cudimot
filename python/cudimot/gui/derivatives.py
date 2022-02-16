"""
CUDIMOT GUI: LM derivatives code
"""
from .code import UserCode

HELP="""\
From the model predicted signal (Not from the cost function).

It returns an array of values, one for each parameter of the model. If you are not going to use
Levenberg-Marquardt algorithm for fitting your model (--runLevMar=false), this function can be
empty (but the function must be declared).
"""

BOILERPLATE = """\
// Partial derivatives respect each model parameter
MACRO void Partial_Derivatives(
			       int npar, // Number of Parameters to estimate
			       T* P, // Estimated parameters, use P*
			       T* CFP, // Fixed Parameters common to all the voxels
			       T* FixP, // Fixed Parameters for each voxel
			       T* derivatives) // Derivative respect each model estimated parameter
{
}
"""

class DerivativesLM(UserCode):
    """
    Tab page containing entry for LM derivatives code
    """

    def __init__(self, app, parent, idx, n):
        UserCode.__init__(self, app, parent, "Partial derivatives for LM", idx, n, name="derivatives_lm", function_name="Partial_Derivatives", help=HELP, boilerplate=BOILERPLATE)

