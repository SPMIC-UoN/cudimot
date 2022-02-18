"""
CUDIMOT GUI: Settings dialog
"""
import os

import wx

class SettingsDialog(wx.Dialog):
    def __init__(self, parent, settings):
        wx.Dialog.__init__(self, parent, title="Settings")
        grid = wx.GridBagSizer(5, 5)

        grid.Add(wx.StaticText(self, label="FSL directory"), pos=(0, 0), flag=wx.ALIGN_CENTRE_VERTICAL | wx.EXPAND | wx.LEFT | wx.ALL, border=5)
        self._fsldir = wx.TextCtrl(self, size=(300, -1))
        self._fsldir.SetValue(settings.get("fsldir", ""))
        self._fsldir.Bind(wx.EVT_TEXT, self._changed)
        grid.Add(self._fsldir, pos=(0, 1), flag=wx.ALIGN_CENTRE_VERTICAL | wx.EXPAND | wx.LEFT | wx.ALL, border=5)

        grid.Add(wx.StaticText(self, label="CUDIMOT executables"), pos=(1, 0), flag=wx.ALIGN_CENTRE_VERTICAL | wx.EXPAND | wx.LEFT | wx.ALL, border=5)
        self._bindir = wx.TextCtrl(self)
        self._bindir.SetValue(settings.get("bindir", ""))
        self._bindir.Bind(wx.EVT_TEXT, self._changed)
        grid.Add(self._bindir, pos=(1, 1), flag=wx.ALIGN_CENTRE_VERTICAL | wx.EXPAND | wx.LEFT | wx.ALL, border=5)

        grid.Add(wx.StaticText(self, label="CUDIMOT source code"), pos=(2, 0), flag=wx.ALIGN_CENTRE_VERTICAL | wx.EXPAND | wx.LEFT | wx.ALL, border=5)
        self._srcdir = wx.TextCtrl(self)
        self._srcdir.SetValue(settings.get("srcdir", ""))
        self._srcdir.Bind(wx.EVT_TEXT, self._changed)
        grid.Add(self._srcdir, pos=(2, 1), flag=wx.ALIGN_CENTRE_VERTICAL | wx.EXPAND | wx.LEFT | wx.ALL, border=5)

        self._warning = wx.StaticText(self)
        grid.Add(self._warning, pos=(3, 0), flag=wx.ALIGN_CENTRE_VERTICAL | wx.EXPAND | wx.LEFT | wx.ALL, border=5, span=(1, 2))
        grid.AddGrowableCol(1, 1)

        button_sizer = self.CreateButtonSizer(wx.OK | wx.CANCEL)

        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(grid, 1, wx.ALL|wx.EXPAND)
        mainSizer.Add(button_sizer, 0, wx.ALL|wx.EXPAND, 5)
        self.SetSizer(mainSizer)
        self._changed()
        self.Fit()

    def config(self):
        return {
            "fsldir" : self._fsldir.GetValue(),
            "bindir" : self._bindir.GetValue(),
            "srcdir" : self._srcdir.GetValue(),
        }

    def _changed(self, evt=None):
        config = self.config()
        if not os.path.exists(os.path.join(config["fsldir"], "bin", "fslmaths")):
            self._warning.SetLabel("WARNING: FSLDIR does not seem to contain an FSL installation")
            return
        
        if not os.path.exists(os.path.join(config["bindir"], "Run_dtifit.sh")):
            self._warning.SetLabel("WARNING: Executable dir does not seem to contain CUDIMOT executables")
            return

        if not os.path.exists(os.path.join(config["srcdir"], "cudimot.cc")):
            self._warning.SetLabel("WARNING: Source dir does not seem to contain CUDIMOT source code")
            return
