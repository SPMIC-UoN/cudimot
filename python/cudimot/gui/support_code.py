"""
CUDIMOT GUI: LM constraints code
"""
import os
import re

import wx.stc

from .code import CodeEditor, comment_status
from .widgets import TabPage

HELP="""\
Here you can add additional code required by your other functions. If you want to
split code into other files, simply #include a filename and an empty file will
be created for you to put the code in.
"""

class SupportCode(TabPage):
    """
    Tab page containing entry for supporting code
    """

    def __init__(self, app, parent, idx, n):
        TabPage.__init__(self, app, parent, "Support code", idx, n, name="support")
        self.skip_functions = ["FixConstraintsLM", "ConstraintsMCMC", "custom_priors", "Partial_Derivatives", "Predicted_Signal"]
        self.pattern = re.compile('#include\s+(?:\"|<)([\w\-\.\d]+)(?:\"|>)')
        self.includes = set()
        self.include_editors = {}

        self.section(self.title)
        self.pack("", wx.StaticText(self, label=HELP))

        self.notebook = wx.Notebook(self, style=wx.BK_DEFAULT)
        self.pack("", self.notebook, span=3, expand=True)

        self.code = CodeEditor(self.notebook)
        self.notebook.AddPage(self.code, "Support code")
        self.code.Bind(wx.stc.EVT_STC_MODIFIED, self._check_includes)

        self.sizer.AddGrowableCol(2, 1)
        self.SetSizer(self.sizer)
        self.add_next_prev_btn()
        self.Layout()
        
    def config(self):
        return {
            "code_%s" % self.name : self.code.GetValue(),
            "support_files" : {f : e.GetValue() for f, e in self.include_editors.items()},
        }

    def load(self, projdir):
        code = self._support_code(os.path.join(projdir, "modelfunctions.h"))
        if not code:
            print(f"No code found")
            code = ""

        self.code.SetValue(code)

    def _support_code(self, fname):
        patterns = [re.compile(f"^.*{name}\s*\(") for name in self.skip_functions]
        code, pre_func_comments = "", ""
        in_skipped_function, in_comment = False, False
        open_brace, close_brace = 0, 0
        try:
            with open(fname, "r") as f:
                for line in f:
                    if not in_skipped_function:
                        for pattern in patterns:
                            if pattern.match(line):
                                in_skipped_function = True
                                break
                    if not in_skipped_function:
                        line_is_comment, in_comment = comment_status(line, in_comment)
                        if line_is_comment:
                            pre_func_comments += line
                            line = ""
                        elif line.strip():
                            code += pre_func_comments
                            pre_func_comments = ""
                    if in_skipped_function:
                        pre_func_comments = ""
                        open_brace += line.count("{")
                        close_brace += line.count("}")
                        if open_brace > 0 and open_brace == close_brace:
                            open_brace, close_brace = 0, 0
                            in_skipped_function = False
                    else:
                        code += line

        except FileNotFoundError:
            print(f"Couldn't find expected file: {self.fname}")
            return None
        return code
        
    def _check_includes(self, evt=None):
        # This could be made more efficient by tracking the line changed but
        # seems to be performant enough for now
        includes = set()
        text = self.code.GetValue()
        for line in text.splitlines():
            match = self.pattern.match(line)
            if match:
                includes.add(match.group(1))
        for fname in includes - self.includes:
            self._create_include(fname)

        for fname in self.includes - includes:
            self._remove_include(fname)

        self.includes = includes

    def _create_include(self, fname):
        self.include_editors[fname] = CodeEditor(self.notebook)
        self.notebook.AddPage(self.include_editors[fname], fname)
        fpath = os.path.join(self.app.projdir, fname)
        if os.path.exists(fpath):
            with open(fpath, "r") as f:
                self.include_editors[fname].SetValue(f.read())

    def _remove_include(self, fname):
        for idx in range(self.notebook.GetPageCount()):
            if self.notebook.GetPageText(idx) == fname:
                self.notebook.DeletePage(idx)
                self.notebook.SendSizeEvent()
                break
        del self.include_editors[fname]
        # FIXME when to remove file? On save?
