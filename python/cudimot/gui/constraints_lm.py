"""
CUDIMOT GUI: LM constraints code
"""
from .code import UserCode

HELP="""\
You can specify what to do with some parameters after Levenberg-Marquardt if a constraint
is not satisfied.
"""

BOILERPLATE = """\
// Constraints run after LevenbergMarquardt (if LevenbergMarquardt is used)
MACRO void FixConstraintsLM(	
	int npar,        // Number of Parameters to estimate
	T* P)            // Estimated parameters
{
}
"""

class ConstraintsLM(UserCode):
    """
    Tab page containing entry for LM constraints code
    """

    def __init__(self, app, parent, idx, n):
        UserCode.__init__(self, app, parent, "LM constraints", idx, n, name="constraints_lm", help=HELP, boilerplate=BOILERPLATE)
