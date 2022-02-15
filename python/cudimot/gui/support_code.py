"""
CUDIMOT GUI: LM constraints code
"""
from .code import UserCode

HELP="""\
Additional code required by your other functions
"""

BOILERPLATE = """\
"""

class SupportCode(UserCode):
    """
    Tab page containing entry for supporting code
    """

    def __init__(self, app, parent, idx, n):
        UserCode.__init__(self, app, parent, "Support code", idx, n, name="support", function_names=["FixConstraintsLM"], invert=True, help=HELP, boilerplate=BOILERPLATE)
