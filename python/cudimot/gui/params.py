"""
CUDIMOT GUI: Model parameters setup
"""
import wx
import wx.grid

from . import get_nvols, OptionError
from .widgets import TabPage, ParameterList, NumberChooser, NumberList

class ModelParameters(TabPage):
    """
    Tab page containing definitions of model parameters
    """

    def __init__(self, app, parent, idx, n):
        TabPage.__init__(self, app, parent, "Model parameters", idx, n, name="params")

        self.section("Variable parameters")
        self.ptables = []
        self.ptables.append(ParameterList(self))
        self.pack("", self.ptables[0], span=3)
        self.section("Common fixed parameters")
        self.ptables.append(ParameterList(self))
        self.pack("", self.ptables[1], span=3)
        self.section("Voxelwise fixed parameters")
        self.ptables.append(ParameterList(self))
        self.pack("", self.ptables[2], span=3)
        self.add_next_prev_btn()

        self.sizer.AddGrowableCol(2, 1)
        self.SetSizer(self.sizer)
        self.Fit()
        self.Layout()
        self.Bind(wx.EVT_SIZE, self._resize_tables)

    def _resize_tables(self, event):
        width = self.GetClientSize()[0] - 10
        for table in self.ptables:
            table.ResizeCols(width)
        event.Skip()

    def options(self):
        options = {
        }
        return options

    def check_options(self, options):
        pass

    def option_changed(self, options, key, value):
        pass
