"""
CUDIMOT GUI: Partial derivatives code
"""
import wx
import wx.grid

from . import get_nvols, OptionError
from .widgets import TabPage, NumberChooser, NumberList

class DerivativesCode(TabPage):
    """
    Tab page containing entry for forward model evaluation code
    """

    def __init__(self, app, parent, idx, n):
        TabPage.__init__(self, app, parent, "Derivatives", idx, n, name="derivatives")

        self.section("Partial derivatives for Levenberg-Marquardt")

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
