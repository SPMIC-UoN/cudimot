"""
CUDIMOT GUI: MCMC constraints code
"""
from .code import UserCode

HELP="""\
Return false if a constraint is not satisfied during MCMC and the sample will be rejected.
By default, it should return true.
"""

BOILERPLATE = """\
// Constraints checked during MCMC (if MCMC is used)
MACRO bool ConstraintsMCMC(
		       int npar, // Number of Parameters to estimate
		       T* P) // Estimated parameters
{
  return true;
}
"""

class ConstraintsMCMC(UserCode):
    """
    Tab page containing entry for MCMC constraints code
    """

    def __init__(self, app, parent, idx, n):
        UserCode.__init__(self, app, parent, "MCMC constraints", idx, n, name="constraints_mcmc", function_names=["ConstraintsMCMC"], help=HELP, boilerplate=BOILERPLATE)
