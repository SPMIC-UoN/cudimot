#!/usr/bin/env python
"""
CUDIMOT GUI

The purpose of this GUI is to offer a friendlier way for users to 
create and compile their own models using the CUDIMOT framework
"""
import sys
import os
import traceback

import wx
import wx.grid
from wx.lib.pubsub import pub

from .project import ProjectOptions
from .params import ModelParameters
from .fwdmodel import ForwardModelCode

class CudimotGui(wx.Frame):
    """
    Main GUI window
    """

    def __init__(self):
        # Initialize main window title, icon etc and vertical box sizer
        wx.Frame.__init__(self, None, title="CUDIMOT", style=wx.DEFAULT_FRAME_STYLE)
        icon_fname = os.path.join(os.path.abspath(os.path.dirname(__file__)), "cudimot.png")
        self.SetIcon(wx.Icon(icon_fname))
        main_panel = wx.Panel(self)
        main_vsizer = wx.BoxSizer(wx.VERTICAL)

        # Add title banner
        banner = wx.Panel(main_panel, size=(-1, 80))
        banner.SetBackgroundColour((54, 122, 157))
        banner_fname = os.path.join(os.path.abspath(os.path.dirname(__file__)), "banner.png")
        wx.StaticBitmap(banner, -1, wx.Bitmap(banner_fname, wx.BITMAP_TYPE_ANY))
        main_vsizer.Add(banner, 0, wx.EXPAND)

        # Main GUI is in tab format
        notebook = wx.Notebook(main_panel, id=wx.ID_ANY, style=wx.BK_DEFAULT)
        main_vsizer.Add(notebook, 1, wx.EXPAND | wx.ALL, 5)

        # Add pages to the tab box
        tabs = [ProjectOptions, ModelParameters, ForwardModelCode]
        for idx, cls in enumerate(tabs):
            tab = cls(self, notebook, idx, len(tabs))
            notebook.AddPage(tab, tab.title)

        # Save/Create/Compile buttons at bottom of window
        line = wx.StaticLine(main_panel, style=wx.LI_HORIZONTAL)
        main_vsizer.Add(line, 0, wx.EXPAND)

        bottom_panel = wx.Panel(main_panel)
        bottom_sizer = wx.BoxSizer(wx.HORIZONTAL)
        bottom_panel.SetSizer(bottom_sizer)
        
        self.save_btn = wx.Button(bottom_panel, label="Save project")
        bottom_sizer.Add(self.save_btn, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)
        self.build_btn = wx.Button(bottom_panel, label="Create code")
        bottom_sizer.Add(self.build_btn, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)
        self.compile_btn = wx.Button(bottom_panel, label="Compile code")
        bottom_sizer.Add(self.compile_btn, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)
        self.status = wx.StaticText(bottom_panel, label="")
        self.status.SetFont(wx.Font(12, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        bottom_sizer.Add(self.status, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)
        #self.run_btn.Bind(wx.EVT_BUTTON, self._do_run)
        main_vsizer.Add(bottom_panel, 0, wx.EXPAND)

        #self.SetMinSize(self.GetSize())
        #self.SetMaxSize(self.GetSize())
        main_panel.SetSizerAndFit(main_vsizer)
        self.Fit()

def main():
    """
    GUI entry point
    """
    app = wx.App(redirect=False)
    top = CudimotGui()
    top.Show()
    app.MainLoop()

if __name__ == '__main__':
    main()