"""
CUDIMOT GUI: Tab page for user-entered code

FIXME: Warn if user tries to modify boilerplate
FIXME: Reset boilerplate only not user code
"""
import wx
import wx.grid
import wx.stc

from . import OptionError
from .widgets import TabPage

class UserCode(TabPage):
    """
    Tab page containing entry for forward model evaluation code
    """

    def __init__(self, app, parent, title, idx, n, name="", help="", function_names=[], boilerplate="", invert=False):
        TabPage.__init__(self, app, parent, title, idx, n, name=name)
        self.function_names = function_names
        self.boilerplate=boilerplate
        self.invert = invert

        self.section(title)
        self.pack("", wx.StaticText(self, label=help))

        self.code = wx.stc.StyledTextCtrl(self, size=(-1, 300))
        #self.code.StyleClearAll()
        self.code.SetLexer(wx.stc.STC_LEX_CPP)
        font = wx.Font(10, wx.FONTFAMILY_TELETYPE, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL)
        face = font.GetFaceName()
        size = font.GetPointSize()
        self.code.StyleSetSpec(wx.stc.STC_STYLE_DEFAULT, "face:%s,size:%d" % (face, size))
        self.code.SetWrapMode (wx.stc.STC_WRAP_WORD) # other choice is wxSCI_WRAP_NONE
        self.code.StyleSetForeground (wx.stc.STC_C_STRING,            wx.Colour(150,0,0))
        self.code.StyleSetForeground (wx.stc.STC_C_PREPROCESSOR,      wx.Colour(165,105,0))
        self.code.StyleSetForeground (wx.stc.STC_C_IDENTIFIER,        wx.Colour(40,0,60))
        self.code.StyleSetForeground (wx.stc.STC_C_NUMBER,            wx.Colour(0,150,0))
        self.code.StyleSetForeground (wx.stc.STC_C_CHARACTER,         wx.Colour(150,0,0))
        self.code.StyleSetForeground (wx.stc.STC_C_WORD,              wx.Colour(0,0,150))
        self.code.StyleSetForeground (wx.stc.STC_C_WORD2,             wx.Colour(0,150,0))
        self.code.StyleSetForeground (wx.stc.STC_C_COMMENT,           wx.Colour(150,150,150))
        self.code.StyleSetForeground (wx.stc.STC_C_COMMENTLINE,       wx.Colour(150,150,150))
        self.code.StyleSetForeground (wx.stc.STC_C_COMMENTDOC,        wx.Colour(150,150,150))
        self.code.StyleSetForeground (wx.stc.STC_C_COMMENTDOCKEYWORD, wx.Colour(0,0,200))
        self.code.StyleSetForeground (wx.stc.STC_C_COMMENTDOCKEYWORDERROR, wx.Colour(0,0,200))
        self.code.StyleSetBold(wx.stc.STC_C_WORD, True)
        self.code.StyleSetBold(wx.stc.STC_C_WORD2, True)
        self.code.StyleSetBold(wx.stc.STC_C_COMMENTDOCKEYWORD, True)
        self.code.SetText(boilerplate)
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
            "code_%s" % self.name : self.code.GetValue()
        }

    def _reset(self):
        self.code.SetValue(self.boilerplate)

    def load(self, projdir):
        import os
        code = self._function_code(os.path.join(projdir, "modelfunctions.h"), [f"^.*{name}\s*\(" for name in self.function_names])
        #print(code)
        if not code:
            print(f"No code found")
            code = self.boilerplate

        self.code.SetValue(code)

    def _function_code(self, fname, func_start_regexes):
        import re
        patterns = [re.compile(r) for r in func_start_regexes]
        code = ""
        in_function = False
        open_brace, close_brace = 0, 0
        try:
            with open(fname, "r") as f:
                for line in f:
                    if not in_function:
                        for pattern in patterns:
                            if pattern.match(line):
                                in_function = True
                                break
                    if in_function:
                        if not self.invert:
                            code += line
                        #print(code)
                        open_brace += line.count("{")
                        close_brace += line.count("}")
                        if open_brace > 0 and open_brace == close_brace:
                            in_function = False
                    elif self.invert:
                        code += line

        except FileNotFoundError:
            print(f"Couldn't find expected file: {self.fname}")
            return None
        return code
        