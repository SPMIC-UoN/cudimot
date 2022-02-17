"""
CUDIMOT GUI: Tab page for user-entered code

FIXME: Warn if user tries to modify boilerplate
FIXME: Reset boilerplate only not user code
"""
import os
import re

import wx
import wx.grid
import wx.stc

from . import OptionError
from .widgets import TabPage

def comment_status(line, currently_in_comment):
    """:return: line_is_comment, now_in_comment"""
    # FIXME not getting C style comments
    line = line.strip()
    cpp = line.startswith("//")
    cstart = line.startswith("/*")
    cend = line.endswith("*/")
    if currently_in_comment and not cend:
        now_in_comment = True
    elif not currently_in_comment and cstart:
        now_in_comment = True
    else:
        now_in_comment = False

    if currently_in_comment:
        line_is_comment = True
    elif cpp or cstart:
        line_is_comment = True
    else:
        line_is_comment = False
       
    return line_is_comment, now_in_comment

class CodeEditor(wx.stc.StyledTextCtrl):
    def __init__(self, parent):
        wx.stc.StyledTextCtrl.__init__(self, parent, size=(-1, 300))
        #self.StyleClearAll()
        self.SetLexer(wx.stc.STC_LEX_CPP)
        font = wx.Font(10, wx.FONTFAMILY_TELETYPE, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL)
        face = font.GetFaceName()
        size = font.GetPointSize()
        self.StyleSetSpec(wx.stc.STC_STYLE_DEFAULT, "face:%s,size:%d" % (face, size))
        self.SetWrapMode (wx.stc.STC_WRAP_WORD) # other choice is wxSCI_WRAP_NONE
        self.StyleSetForeground (wx.stc.STC_C_STRING,            wx.Colour(150,0,0))
        self.StyleSetForeground (wx.stc.STC_C_PREPROCESSOR,      wx.Colour(165,105,0))
        self.StyleSetForeground (wx.stc.STC_C_IDENTIFIER,        wx.Colour(40,0,60))
        self.StyleSetForeground (wx.stc.STC_C_NUMBER,            wx.Colour(0,150,0))
        self.StyleSetForeground (wx.stc.STC_C_CHARACTER,         wx.Colour(150,0,0))
        self.StyleSetForeground (wx.stc.STC_C_WORD,              wx.Colour(0,0,150))
        self.StyleSetForeground (wx.stc.STC_C_WORD2,             wx.Colour(0,150,0))
        self.StyleSetForeground (wx.stc.STC_C_COMMENT,           wx.Colour(150,150,150))
        self.StyleSetForeground (wx.stc.STC_C_COMMENTLINE,       wx.Colour(150,150,150))
        self.StyleSetForeground (wx.stc.STC_C_COMMENTDOC,        wx.Colour(150,150,150))
        self.StyleSetForeground (wx.stc.STC_C_COMMENTDOCKEYWORD, wx.Colour(0,0,200))
        self.StyleSetForeground (wx.stc.STC_C_COMMENTDOCKEYWORDERROR, wx.Colour(0,0,200))
        self.StyleSetBold(wx.stc.STC_C_WORD, True)
        self.StyleSetBold(wx.stc.STC_C_WORD2, True)
        self.StyleSetBold(wx.stc.STC_C_COMMENTDOCKEYWORD, True)

class UserCode(TabPage):
    """
    Tab page containing entry for forward model evaluation code
    """

    def __init__(self, app, parent, title, idx, n, name="", help="", function_name="", boilerplate=""):
        TabPage.__init__(self, app, parent, title, idx, n, name=name)
        self.function_name = function_name
        self.boilerplate=boilerplate

        self.section(title)
        self.pack("", wx.StaticText(self, label=help), span=3)
        self.code = CodeEditor(self)
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

    def reset(self, projdir):
        self.code.SetValue(self.boilerplate)

    def _reset(self, evt=None):
        self.reset("")

    def load(self, projdir):
        code = self._function_code(os.path.join(projdir, "modelfunctions.h"), )
        #print(code)
        if not code:
            print(f"No code found")
            code = self.boilerplate

        self.code.SetValue(code)

    def _function_code(self, fname):
        """
        Extract function code (and pre-function comments) from a file
        """
        pattern = re.compile(f"^.*{self.function_name}\s*\(")
        code, pre_func_comments = "", ""
        in_function, in_comment = False, False
        open_brace, close_brace = 0, 0
        try:
            with open(fname, "r") as f:
                for line in f:
                    if not in_function and pattern.match(line):
                        in_function = True
                    if not in_function:
                        line_is_comment, in_comment = comment_status(line, in_comment)
                        if line_is_comment:
                            pre_func_comments += line
                        elif line.strip():
                            pre_func_comments = ""
                    else:
                        code += pre_func_comments
                        pre_func_comments = ""
                        code += line
                        #print(code)
                        open_brace += line.count("{")
                        close_brace += line.count("}")
                        if open_brace > 0 and open_brace == close_brace:
                            open_brace, close_brace = 0, 0
                            in_function = False

        except FileNotFoundError:
            print(f"Couldn't find expected file: {self.fname}")
            return None
        return code
