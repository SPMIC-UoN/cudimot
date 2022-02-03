"""
CUDIMOT GUI: Main project options
"""
import wx
import wx.grid

from . import get_nvols, OptionError
from .widgets import TabPage, NumberChooser, NumberList

class ProjectOptions(TabPage):
    """
    Tab page containing basic project options, e.g. path to project, name of model etc
    """

    def __init__(self, app, parent, idx, n):
        TabPage.__init__(self, app, parent, "Basic model options", idx, n, name="input")
        self._nvols = -1 # Cached to improve responsiveness

        self.section("Basic configuration")
        self.model_name = self.text("Model name", size=(300, -1))
        self.projdir = self.file_picker("Model Directory", pick_dir=True)
        self.precision = self.choice("Floating point precision", ["single", "double"], initial=1)

        self.sizer.AddGrowableCol(2, 1)
        self.SetSizer(self.sizer)
        self.add_next_prev_btn()
        self.Layout()

    def options(self):
        options = {
            "name"        : self.model_name.GetText(),
            "projdir"     : self.projdir.GetPath(),
            "precision"     : self.precision.GetValue(),
        }
        return options

    def check_options(self, options):
        pass

    def option_changed(self, options, key, value):
        pass
