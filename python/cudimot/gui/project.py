"""
CUDIMOT GUI: Main project options
"""
import os

import wx
import wx.grid

from . import OptionError
from .widgets import TabPage, NumberChooser, NumberList

class ProjectOptions(TabPage):
    """
    Tab page containing basic project options, e.g. path to project, name of model etc
    """

    def __init__(self, app, parent, idx, n):
        TabPage.__init__(self, app, parent, "Basic options", idx, n, name="input")
        self._nvols = -1 # Cached to improve responsiveness

        self.section("Basic configuration")
        self.projdir = self.text("Project Directory", style=wx.TE_READONLY)
        self.model_name = self.text("Model name", size=(300, -1))
        self.precision = self.choice("Floating point precision", ["single", "double"], initial=1)

        self.sizer.AddGrowableCol(2, 1)
        self.SetSizer(self.sizer)
        self.add_next_prev_btn()
        self.Layout()

    def reset(self, projdir):
        self.projdir.SetValue(projdir)
        self.model_name.SetValue(os.path.basename(projdir))

    def config(self):
        precision = self.precision.GetString(self.precision.GetSelection())
        return {
            "name"        : self.model_name.GetValue(),
            "projdir"     : self.projdir.GetValue(),
            "precision"   : precision,
            "dtype"       : "double" if precision == "double" else "float",
        }

    def load(self, projdir):
        self.projdir.SetValue(projdir)
        precision = self.config_from_line_regex("precision", 
                                                os.path.join(projdir, "modelparameters.h"), 
                                                "typedef\s+(\w+)\s+MyType\s*;")
        if precision in ("single", "double"):
            self.precision.SetSelection(self.precision.FindString(precision))
        elif precision:
            print(f"Unrecogniazed precision: {precision}")
        else:
            print(f"Precision not found")

        modelname = self.config_from_line_regex("modelname",
                                                os.path.join(projdir, "Pipeline_*.sh"), 
                                                "modelname\s*=\s*(\w+)")
        if not modelname:
            # Fall back on project directory basename
            modelname = os.path.basename(projdir)
        self.model_name.SetValue(modelname)




