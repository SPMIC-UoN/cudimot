#!/usr/bin/env python
"""
CUDIMOT GUI

The purpose of this GUI is to offer a friendlier way for users to 
create and compile their own models using the CUDIMOT framework
"""
import sys
import os
import traceback
import yaml

import wx
import wx.grid
from wx.lib.pubsub import pub

from .project import ProjectOptions
from .params import ModelParameters
from .fwdmodel import ForwardModel
from .constraints_mcmc import ConstraintsMCMC
from .derivatives import DerivativesLM
from .constraints_lm import ConstraintsLM
from .custom_priors import CustomPriors
from .support_code import SupportCode
from .settings import SettingsDialog
import cudimot.gui.code_templates as templates

class CudimotGui(wx.Frame):
    """
    Main GUI window
    """

    def __init__(self):
        # Initialize main window title, icon etc and vertical box sizer
        wx.Frame.__init__(self, None, title="CUDIMOT", style=wx.DEFAULT_FRAME_STYLE)
        self.projdir = ""
        self.fsldir = os.environ.get("FSLDIR", "")
        self.settings = {
            "fsldir" : self.fsldir,
            "bindir" : os.path.join(self.fsldir, "bin"),
            "srcdir" : os.path.join(self.fsldir, "src", "cudimot"),
        }
        icon_fname = os.path.join(os.path.abspath(os.path.dirname(__file__)), "cudimot.png")
        self.SetIcon(wx.Icon(icon_fname))
        main_panel = wx.Panel(self)
        main_vsizer = wx.BoxSizer(wx.VERTICAL)

        # Add title banner
        banner = wx.Panel(main_panel, size=(-1, 80))
        banner.SetBackgroundColour((0, 0, 0))
        banner_fname = os.path.join(os.path.abspath(os.path.dirname(__file__)), "banner.png")
        wx.StaticBitmap(banner, -1, wx.Bitmap(banner_fname, wx.BITMAP_TYPE_ANY))
        main_vsizer.Add(banner, 0, wx.EXPAND)

        # Main GUI is in tab format
        self.notebook = wx.Notebook(main_panel, id=wx.ID_ANY, style=wx.BK_DEFAULT)
        main_vsizer.Add(self.notebook, 1, wx.EXPAND | wx.ALL, 5)

        # Add pages to the tab box
        tabs = [ProjectOptions, ModelParameters, ForwardModel, ConstraintsMCMC, DerivativesLM, ConstraintsLM, CustomPriors, SupportCode]
        for idx, cls in enumerate(tabs):
            tab = cls(self, self.notebook, idx, len(tabs))
            self.notebook.AddPage(tab, tab.title)
            tab.Enable(False)
        self.notebook.SendSizeEvent()

        # Save/Create/Compile buttons at bottom of window
        line = wx.StaticLine(main_panel, style=wx.LI_HORIZONTAL)
        main_vsizer.Add(line, 0, wx.EXPAND)

        bottom_panel = wx.Panel(main_panel)
        bottom_sizer = wx.BoxSizer(wx.HORIZONTAL)
        bottom_panel.SetSizer(bottom_sizer)

        self.new_btn = wx.Button(bottom_panel, label="New")
        bottom_sizer.Add(self.new_btn, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)
        self.new_btn.Bind(wx.EVT_BUTTON, self._new)
        self.open_btn = wx.Button(bottom_panel, label="Open")
        bottom_sizer.Add(self.open_btn, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)
        self.open_btn.Bind(wx.EVT_BUTTON, self._open)
        self.save_btn = wx.Button(bottom_panel, label="Save")
        self.save_btn.Enable(False)
        self.save_btn.Bind(wx.EVT_BUTTON, self._save)
        bottom_sizer.Add(self.save_btn, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)
        self.save_as_btn = wx.Button(bottom_panel, label="Save As")
        self.save_as_btn.Enable(False)
        self.save_as_btn.Bind(wx.EVT_BUTTON, self._save_as)
        bottom_sizer.Add(self.save_as_btn, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)
        self.compile_btn = wx.Button(bottom_panel, label="Build")
        self.compile_btn.Enable(False)
        #self.compile_btn.Bind(wx.EVT_BUTTON, self._build)
        bottom_sizer.Add(self.compile_btn, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)
        self.status = wx.StaticText(bottom_panel, label="")
        self.status.SetFont(wx.Font(12, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        bottom_sizer.Add(self.status, 1, wx.EXPAND | wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)
        self.settings_btn = wx.Button(bottom_panel, label="Settings")
        self.settings_btn.Enable(True)
        self.settings_btn.Bind(wx.EVT_BUTTON, self._settings)
        bottom_sizer.Add(self.settings_btn, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)
        main_vsizer.Add(bottom_panel, 0)

        #self.SetMinSize(self.GetSize())
        #self.SetMaxSize(self.GetSize())
        main_panel.SetSizerAndFit(main_vsizer)
        self.Fit()

    def _settings(self, evt=None):
        with SettingsDialog(self, self.settings) as dialog:
            if dialog.ShowModal() == wx.ID_CANCEL:
                return
            self.settings = dialog.config()
            print(self.settings)

    def _new(self, evt=None):
        with wx.DirDialog(self, "Select project folder") as fileDialog:
            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return

            projdir = fileDialog.GetPath()
            if os.path.exists(projdir) and os.listdir(projdir):
                with wx.MessageDialog(self, "Directory is not empty - Overwrite?", caption="Confirm",
                                        style=wx.OK|wx.CANCEL) as confirm_dialog:
                    if confirm_dialog.ShowModal() == wx.ID_CANCEL:
                        return

            self.projdir = projdir
            for idx in range(self.notebook.PageCount):
                page = self.notebook.GetPage(idx)
                page.Enable(True)
                page.reset(self.projdir)
            self.save_btn.Enable(True)
            self.save_as_btn.Enable(True)

    def _open(self, evt=None):
        with wx.DirDialog(self, "Open project folder",
                       style=wx.DD_DIR_MUST_EXIST) as fileDialog:
            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return

            self.projdir = fileDialog.GetPath()
            for idx in range(self.notebook.PageCount):
                page = self.notebook.GetPage(idx)
                page.Enable(True)
                page.load(self.projdir)
            self.save_btn.Enable(True)
            self.save_as_btn.Enable(True)

    def _save_as(self, _event=None):
        """
        Save project in a new location
        """
        with wx.DirDialog(self, "Select save folder") as fileDialog:
            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return

            projdir = fileDialog.GetPath()
            if os.path.exists(projdir) and os.listdir(projdir):
                with wx.MessageDialog(self, "Directory is not empty - Overwrite?", caption="Confirm",
                                      style=wx.OK|wx.CANCEL) as confirm_dialog:
                    if confirm_dialog.ShowModal() == wx.ID_CANCEL:
                        return

            self.projdir = projdir
            # FIXME hack
            self.notebook.GetPage(0).projdir.SetValue(self.projdir)
            self._save()

    def _save(self, _event=None):
        """
        Save project
        
        We save the project configuration in a YAML file and also autogenerate the
        code from it
        """
        config = dict(self.settings)
        for idx in range(self.notebook.PageCount):
            config.update(self.notebook.GetPage(idx).config())
        print(config)

        if not os.path.exists(self.projdir):
            os.makedirs(self.projdir)
        elif not os.path.isdir(self.projdir):
            raise ValueError(f"{self.projdir} already exists and is not a directory")

        # YAML config - FIXME do we need this or does it just muddy the waters
        with open(os.path.join(self.projdir, "cudimot_config.yml"), "w") as f:
            f.write(yaml.dump(config))

        # Create modelparameters.h
        with open(os.path.join(self.projdir, "modelparameters.h"), "w") as f:
            f.write(templates.HEADER)
            f.write(templates.MODELPARAMETERS_H.format(**config))

        # Create modelparameters.cc
        with open(os.path.join(self.projdir, "modelparameters.cc"), "w") as f:
            f.write(templates.HEADER)
            f.write(templates.MODELPARAMETERS_CC.format(**config))

        # Create modelfunctions.h
        with open(os.path.join(self.projdir, "modelfunctions.h"), "w") as f:
            f.write(templates.MODELFUNCTIONS_H.format(**config))

        # Create model info
        with open(os.path.join(self.projdir, f"{config['name']}.info"), "w") as f:
            f.write(templates.MODELINFO.format(**config))

        # Create modelpriors
        bounds_spec, priors_spec = "", ""
        for idx, param in enumerate(config["params"]):
            bounds_spec += f"bounds[{idx}]=({param['lbound']}, {param['ubound']})\n"
            if param["prior"] != "None":
                priors_spec += f"priors[{idx}]={param['prior']}\n"
        config["bounds_spec"] = bounds_spec
        config["priors_spec"] = priors_spec
        with open(os.path.join(self.projdir, "modelpriors"), "w") as f:
            f.write(templates.HEADER)
            f.write(templates.MODELPRIORS.format(**config))

        # Create support files
        for fname, code in config["support_files"].items():
            with open(os.path.join(self.projdir, fname), "w") as f:
                f.write(code)

        # Create Makefile
        with open(os.path.join(self.projdir, "Makefile"), "w") as f:
            f.write(templates.MAKEFILE.format(**config))

        # Create wrapper scripts
        # FIXME this is the 'generic script' also finish script?
        with open(os.path.join(self.projdir, f"cudimot_{config['name']}.sh"), "w") as f:
            f.write(templates.WRAPPER_SCRIPT.format(**config))

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
