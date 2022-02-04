"""
CUDIMOT GUI: Tab page for user-entered code

FIXME: Warn if user tries to modify boilerplate
FIXME: Reset boilerplate only not user code
"""
import wx
import wx.grid

from . import OptionError
from .widgets import TabPage

class UserCode(TabPage):
    """
    Tab page containing entry for forward model evaluation code
    """

    def __init__(self, app, parent, title, idx, n, name="", help="", boilerplate=""):
        TabPage.__init__(self, app, parent, title, idx, n, name=name)

        self.section(title)
        self.pack("", wx.StaticText(self, label=help))
        self.code = wx.TextCtrl(self, style=wx.TE_MULTILINE, size=(-1, 300))
        self.code.SetValue(boilerplate)
        self.pack("", self.code, span=3, expand=True)
        self.reset_btn = wx.Button(self, label="Reset")
        self.reset_btn.Bind(wx.EVT_BUTTON, self._reset)
        self.pack("", self.reset_btn)

        self.sizer.AddGrowableCol(2, 1)
        self.SetSizer(self.sizer)
        self.add_next_prev_btn()
        self.Layout()

    def config(self):
        return {
            "%s_code" % self.name : self.code.GetValue()
        }

    def _reset(self):
        self.code.SetValue(BOILERPLATE)
