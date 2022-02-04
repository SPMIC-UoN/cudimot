"""
CUDIMOT GUI: Main project options
"""
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
        self.projdir = self.file_picker("Project Directory", pick_dir=True)
        self.model_name = self.text("Model name", size=(300, -1))
        self.precision = self.choice("Floating point precision", ["single", "double"], initial=1)

        self.sizer.AddGrowableCol(2, 1)
        self.SetSizer(self.sizer)
        self.add_next_prev_btn()
        self.Layout()

    def config(self):
        return {
            "name"        : self.model_name.GetValue(),
            "projdir"     : self.projdir.GetPath(),
            "precision"     : self.precision.GetString(self.precision.GetSelection()),
        }
