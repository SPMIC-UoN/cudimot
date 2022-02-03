"""
CUDIMOT GUI: Forward model code
"""
import wx
import wx.grid

from . import get_nvols, OptionError
from .widgets import TabPage, NumberChooser, NumberList

TEXT = """\
Here you can specify the code to calculate your model, i.e. 
given the values of the parameters return the model prediction"""

class ForwardModelCode(TabPage):
    """
    Tab page containing entry for forward model evaluation code
    """

    def __init__(self, app, parent, idx, n):
        TabPage.__init__(self, app, parent, "Forward model", idx, n, name="fwdmodel")

        self.section("Forward model calculation")
        self.pack("", wx.StaticText(self, label=TEXT))
        self.fwdmodel_code = wx.TextCtrl(self, style=wx.TE_MULTILINE, size=(-1, 200))
        self.pack("", self.fwdmodel_code)
        self.sizer.AddGrowableCol(0, 1)
        self.SetSizer(self.sizer)
        self.add_next_prev_btn()
        self.Layout()

    def options(self):
        options = {
        }
        return options

    def check_options(self, options):
        pass

    def option_changed(self, options, key, value):
        pass
