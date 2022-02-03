"""
BASIL GUI for Oxford ASL - Base classes for pages in the tab notebook

Copyright (C) 2020 University of Oxford
"""
import os

import numpy as np

import wx
import wx.grid
import wx.adv

from . import OptionComponent

class TabPage(wx.Panel, OptionComponent):
    """
    Shared methods used by the various tab pages in the GUI
    """
    def __init__(self, app, notebook, title, tab_idx, num_tabs, name=None, **kwargs):
        OptionComponent.__init__(self, app)
        wx.Panel.__init__(self, parent=notebook, id=wx.ID_ANY, **kwargs)
        self.notebook = notebook
        self.tab_idx = tab_idx
        self.num_tabs = num_tabs
        self.sizer = wx.GridBagSizer(vgap=5, hgap=5)
        self.row = 0
        self.title = title
        if name is None:
            self.name = title.lower()
        else:
            self.name = name

    def add_next_prev_btn(self):
        """
        Add next/previous buttons
        """
        if self.tab_idx < self.num_tabs-1:
            next_btn = wx.Button(self, label="Next", id=wx.ID_FORWARD)
            next_btn.Bind(wx.EVT_BUTTON, self._next)
        else:
            next_btn = wx.StaticText(self, label="")

        if self.tab_idx > 0:
            prev_btn = wx.Button(self, label="Previous", id=wx.ID_BACKWARD)
            prev_btn.Bind(wx.EVT_BUTTON, self._prev)
        else:
            prev_btn = wx.StaticText(self, label="")

        self.pack(" ")
        self.sizer.AddGrowableRow(self.row-1, 1)
        self.sizer.Add(prev_btn, pos=(self.row, 0), border=10, flag=wx.ALIGN_CENTRE_VERTICAL | wx.ALIGN_LEFT)
        self.sizer.Add(wx.StaticText(self, label=""), pos=(self.row, 1), border=10, flag=wx.ALIGN_CENTRE_VERTICAL | wx.ALIGN_LEFT)
        #self.sizer.Add(wx.StaticText(self, label=""), pos=(self.row, 2), border=10, flag=wx.ALIGN_CENTRE_VERTICAL | wx.ALIGN_LEFT)
        self.sizer.Add(next_btn, pos=(self.row, 2), border=10, flag=wx.ALIGN_CENTRE_VERTICAL | wx.ALIGN_RIGHT)

    def _next(self, _evt):
        self.notebook.SetSelection(self.tab_idx+1)

    def _prev(self, _evt):
        self.notebook.SetSelection(self.tab_idx-1)

    def pack(self, label, *widgets, **kwargs):
        """
        Add a horizontal line to the tab with a label and series of widgets

        If label is empty, first widget is used instead (usually to provide a checkbox)
        """
        col = 0
        border = kwargs.get("border", 10)
        font = self.GetFont()
        if "font_size" in kwargs:
            font.SetPointSize(kwargs["font_size"])
        if kwargs.get("bold", False):
            font.SetWeight(wx.BOLD)

        if label:
            text = wx.StaticText(self, label=label)
            text.SetFont(font)
            self.sizer.Add(text, pos=(self.row, col), border=border, flag=wx.ALIGN_CENTRE_VERTICAL | wx.LEFT)
            col += 1
        else:
            text = None

        for widget in widgets:
            widget.label = text
            if hasattr(widget, "span"):
                span = (1, widget.span)
            else:
                span = (1, kwargs.get("span", 1))
            widget.SetFont(font)
            widget.Enable(col == 0 or kwargs.get("enable", True))
            self.sizer.Add(widget, pos=(self.row, col), border=border, flag=wx.ALIGN_CENTRE_VERTICAL | wx.EXPAND | wx.LEFT, span=span)
            col += span[1]
        self.row += 1

    def file_picker(self, label, pick_dir=False, handler=None, optional=False, initial_on=False, pack=True, **kwargs):
        """
        Add a file picker to the tab
        """
        if not handler:
            handler = self.state_changed
        if pick_dir:
            picker = wx.DirPickerCtrl(self, style=wx.DIRP_USE_TEXTCTRL)
            picker.Bind(wx.EVT_DIRPICKER_CHANGED, handler)
        else:
            picker = wx.FilePickerCtrl(self)
            picker.Bind(wx.EVT_FILEPICKER_CHANGED, handler)
        picker.span = 2
        if optional:
            checkbox = wx.CheckBox(self, label=label)
            checkbox.SetValue(initial_on)
            checkbox.Bind(wx.EVT_CHECKBOX, self._checkbox_toggle_cb(checkbox, picker))
            picker.checkbox = checkbox
            if pack:
                self.pack("", checkbox, picker, enable=initial_on, **kwargs)
        elif pack:
            self.pack(label, picker, **kwargs)

        return picker

    def set_picker_dir(self, filepicker, fpath):
        """
        Set the initial directory for a file or dir picker

        :param fpath: Path to file or directory
        """
        dirname = os.path.dirname(fpath)
        try:
            filepicker.SetInitialDirectory(dirname)
        except AttributeError:
            # WX version dependent - so try alternate name
            filepicker.SetPath(dirname)

    def _checkbox_toggle_cb(self, checkbox, widget):
        def _toggled(_event):
            widget.Enable(checkbox.IsChecked())
            self.state_changed(_event)
        return _toggled

    def choice(self, label, choices, initial=0, optional=False, initial_on=False, handler=None, pack=True, **kwargs):
        """
        Add a widget to choose from a fixed set of options
        """
        if not handler:
            handler = self.state_changed
        choice = wx.Choice(self, choices=choices)
        choice.SetSelection(initial)
        choice.Bind(wx.EVT_CHOICE, handler)
        if optional:
            checkbox = wx.CheckBox(self, label=label)
            checkbox.SetValue(initial_on)
            checkbox.Bind(wx.EVT_CHECKBOX, self.state_changed)
            choice.checkbox = checkbox
            if pack:
                self.pack("", checkbox, choice, enable=initial_on, **kwargs)
        elif pack:
            self.pack(label, choice, **kwargs)
        return choice

    def text(self, label, handler=None, **kwargs):
        """
        Add a widget to input text
        """
        if not handler:
            handler = self.state_changed
        text = wx.TextCtrl(self, -1, **kwargs)
        text.Bind(wx.EVT_TEXT, handler)
        text.span = 2
        self.pack(label, text, **kwargs)
        return text

    def number(self, label, handler=None, **kwargs):
        """
        Add a widget to choose a floating point number
        """
        if not handler:
            handler = self.state_changed
        num = NumberChooser(self, changed_handler=handler, **kwargs)
        num.span = 2
        self.pack(label, num, **kwargs)
        return num

    def integer(self, label, handler=None, pack=True, **kwargs):
        """
        Add a widget to choose an integer
        """
        if not handler:
            handler = self.state_changed
        spin = wx.SpinCtrl(self, **kwargs)
        spin.SetValue(kwargs.get("initial", 0))
        spin.Bind(wx.EVT_SPINCTRL, handler)
        if pack:
            self.pack(label, spin)
        return spin

    def checkbox(self, label, *extra_widgets, initial=False, handler=None, **kwargs):
        """
        Add a simple on/off option
        """
        checkbox = wx.CheckBox(self, label=label)
        checkbox.span = kwargs.get("span", 2)
        checkbox.SetValue(initial)
        if handler:
            checkbox.Bind(wx.EVT_CHECKBOX, handler)
        else:
            checkbox.Bind(wx.EVT_CHECKBOX, self.state_changed)
        self.pack("", checkbox, *extra_widgets, **kwargs)
        return checkbox

    def section(self, label):
        """
        Add a section heading
        """
        self.pack(label, bold=True)

class NumberChooser(wx.Panel):
    """
    Widget for choosing a floating point number
    """

    def __init__(self, parent, label=None, minval=0, maxval=1, initial=0.5, step=0.1, digits=2, changed_handler=None):
        super(NumberChooser, self).__init__(parent)
        self.minval, self.orig_minval, self.maxval, self.orig_maxval = minval, minval, maxval, maxval
        self.handler = changed_handler
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        if label is not None:
            self.label = wx.StaticText(self, label=label)
            self.hbox.Add(self.label, proportion=0, flag=wx.ALIGN_CENTRE_VERTICAL)
        # Set a very large maximum as we want to let the user override the default range
        #self.spin = wx.SpinCtrl(self, minval=0, maxval=100000, initial=initial)
        #self.spin.Bind(wx.EVT_SPINCTRL, self._spin_changed)
        self.spin = wx.SpinCtrlDouble(self, min=0, max=100000, inc=step, initial=initial)
        self.spin.SetDigits(digits)
        self.spin.Bind(wx.EVT_SPINCTRLDOUBLE, self._spin_changed)
        self.slider = wx.Slider(self, value=initial, minValue=0, maxValue=100)
        self.slider.SetValue(100*(initial-self.minval)/(self.maxval-self.minval))
        self.slider.Bind(wx.EVT_SLIDER, self._slider_changed)
        self.hbox.Add(self.slider, proportion=1, flag=wx.EXPAND)
        self.hbox.Add(self.spin, proportion=0, flag=wx.EXPAND)
        self.SetSizer(self.hbox)

    def GetValue(self):
        """
        Get the currently selected number
        """
        return self.spin.GetValue()

    def SetValue(self, val):
        """
        Set the selected number
        """
        self.spin.SetValue(val)
        self.slider.SetValue(100*(val-self.minval)/(self.maxval-self.minval))

    def _slider_changed(self, event):
        slider_pos = event.GetInt()
        val = self.minval + (self.maxval-self.minval)*float(slider_pos)/100
        self.spin.SetValue(val)
        if self.handler:
            self.handler(event)
        event.Skip()

    def _spin_changed(self, event):
        """ If user sets the spin outside the current range, update the slider range
        to match. However if they go back inside the current range, revert to this for
        the slider"""
        val = event.GetValue()
        if val < self.minval:
            self.minval = val
        elif val > self.orig_minval:
            self.minval = self.orig_minval
        if val > self.maxval:
            self.maxval = val
        elif val < self.orig_maxval:
            self.maxval = self.maxval
        self.slider.SetValue(100*(val-self.minval)/(self.maxval-self.minval))
        if self.handler:
            self.handler()
        event.Skip()

class ParameterList(wx.Panel):
    """
    Widget for specifying a list of parameters
    """

    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        self.nparams = 0

        self.sizer = wx.BoxSizer(wx.VERTICAL)

        self.grid = wx.grid.Grid(self)
        self.grid.CreateGrid(0, 2)
        self.grid.SetRowLabelSize(0)
        self.grid.SetColLabelValue(0, "Parameter name")
        self.grid.SetColLabelValue(1, "Parameter size")
        font = self.GetFont()
        font.SetWeight(wx.NORMAL)
        self.grid.SetLabelFont(font)
        #self.SetColLabelSize(0)
        self.grid.SetSelectionMode(self.grid.SelectRows)
        self.AddParameter()
        self.grid.Fit()
        #self.grid.Bind(wx.EVT_SIZE, self._on_size)
        self.sizer.Add(self.grid, 1, wx.EXPAND)

        self.button_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.add = wx.Button(self, label="Add")
        self.button_sizer.Add(self.add)
        self.add.Bind(wx.EVT_BUTTON, self.AddParameter)
        self.delete = wx.Button(self, label="Delete")
        self.button_sizer.Add(self.delete)
        self.delete.Bind(wx.EVT_BUTTON, self.DeleteSelectedParameters)
        
        self.sizer.Add(self.button_sizer)

        #self.SetNumValues(size)
        #self.Bind(wx.EVT_SIZE, self._on_size)
        self.SetSizer(self.sizer)
        self.Fit()

    def AddParameter(self, event=None):
        self.grid.AppendRows()
        self.nparams += 1
        self.grid.SetCellValue(self.nparams-1, 0, "param%i" % self.nparams)
        self.grid.SetCellValue(self.nparams-1, 1, "1")

    def DeleteSelectedParameters(self, event=None):
        for row in self.grid.GetSelectedRows():
            self.grid.DeleteRows(row)
            self.nparams -= 1

    def GetParameters(self):
        """
        :return: Sequence of parameters in the list
        """
        params = []
        for p_idx in range(self.nparams):
            try:
                params.append({
                    "name" : self.grid.GetCellValue(p_idx, 0),
                    "size" : int(self.grid.GetCellValue(p_idx, 1)),
                })
            except ValueError:
                pass # FIXME how to handle invalid sizes
        return params

    def ResizeCols(self, width):
        num_width = 100
        name_width = max(100, width - num_width)
        self.grid.SetColSize(0, name_width)
        self.grid.SetColSize(1, num_width)

class NumberList(wx.grid.Grid):
    """
    Widget for specifying a list of numbers
    """

    def __init__(self, parent, size, default=1.8):
        super(NumberList, self).__init__(parent, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, 0)
        self.size = 0
        self.default = default
        self.CreateGrid(1, 0)
        self.SetRowLabelSize(0)
        self.SetColLabelSize(0)
        self.SetNumValues(size)
        self.Bind(wx.EVT_SIZE, self._on_size)

    def GetValues(self):
        """
        :return: Sequence of values in the list
        """
        vals = []
        for c in range(self.size):
            try:
                vals.append(float(self.GetCellValue(0, c)))
            except ValueError:
                vals.append(-999)
        return vals

    def SetNumValues(self, size, default=None):
        """
        Set the size of the number list

        :param size: Number of items in list
        :param default: Default value to use for newly created columns. If not specified
                        uses default defined in constructor
        """
        if default is None:
            if self.size == 0:
                default = self.default
            else:
                default = self.GetCellValue(0, self.size-1)
        if size > self.size:
            self.AppendCols(size - self.size)
            for col in range(self.size, size):
                self.SetCellValue(0, col, str(default))
        elif size < self.size:
            self.DeleteCols(size, self.size-size)
        self.size = size
        self._resize_cols()

    def _resize_cols(self):
        if self.size == 0:
            return

        width, _height = self.GetClientSize()
        col_width = (width - 5) / self.size
        for i in range(self.size):
            self.SetColSize(i, col_width)

    def _on_size(self, event):
        self._resize_cols()
        event.Skip()
