# -*- coding: UTF-8 -*-
#
# generated by wxGlade 1.0.5 on Wed May 17 19:37:45 2023
#

import wx

# begin wxGlade: dependencies
# end wxGlade

# begin wxGlade: extracode
# end wxGlade


class MyDialog(wx.Dialog):
    def __init__(self, *args, **kwds):
        # begin wxGlade: MyDialog.__init__
        kwds["style"] = kwds.get("style", 0) | wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER
        wx.Dialog.__init__(self, *args, **kwds)
        self.SetTitle("Manage Experiments")

        sizer_4 = wx.BoxSizer(wx.VERTICAL)

        grid_sizer_1 = wx.FlexGridSizer(3, 2, 10, 0)
        sizer_4.Add(grid_sizer_1, 1, wx.ALL | wx.EXPAND, 10)

        grid_sizer_1.Add((20, 20), 0, 0, 0)

        sizer_2 = wx.BoxSizer(wx.HORIZONTAL)
        grid_sizer_1.Add(sizer_2, 1, wx.EXPAND, 0)

        sizer_2_copy = wx.BoxSizer(wx.HORIZONTAL)
        sizer_2.Add(sizer_2_copy, 1, wx.EXPAND, 0)

        self.label_2 = wx.StaticText(self, wx.ID_ANY, "Isotope:")
        sizer_2_copy.Add(self.label_2, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)

        self.ComboIsotope = wx.ComboBox(self, wx.ID_ANY, choices=[], style=wx.CB_DROPDOWN | wx.CB_READONLY)
        sizer_2_copy.Add(self.ComboIsotope, 0, 0, 0)

        sizer_1 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_2.Add(sizer_1, 1, wx.EXPAND, 0)

        self.label_1 = wx.StaticText(self, wx.ID_ANY, "B0:")
        sizer_1.Add(self.label_1, 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT | wx.RIGHT, 5)

        self.ComboB0 = wx.ComboBox(self, wx.ID_ANY, choices=[], style=wx.CB_DROPDOWN | wx.CB_READONLY)
        sizer_1.Add(self.ComboB0, 0, 0, 0)

        sizer_5 = wx.BoxSizer(wx.VERTICAL)
        grid_sizer_1.Add(sizer_5, 1, wx.BOTTOM | wx.EXPAND | wx.LEFT | wx.RIGHT, 10)

        self.ButtonView = wx.Button(self, wx.ID_ANY, "&View...")
        sizer_5.Add(self.ButtonView, 0, wx.TOP, 5)

        self.ButtonDelete = wx.Button(self, wx.ID_DELETE, "")
        sizer_5.Add(self.ButtonDelete, 0, wx.TOP, 30)

        self.ListExperiments = wx.ListCtrl(self, wx.ID_ANY, style=wx.BORDER_SUNKEN | wx.LC_REPORT)
        grid_sizer_1.Add(self.ListExperiments, 1, wx.EXPAND, 0)

        grid_sizer_1.Add((20, 20), 0, 0, 0)

        sizer_6 = wx.BoxSizer(wx.HORIZONTAL)
        grid_sizer_1.Add(sizer_6, 1, wx.BOTTOM | wx.EXPAND, 10)

        self.ButtonImport = wx.Button(self, wx.ID_ANY, "&Import...")
        sizer_6.Add(self.ButtonImport, 0, 0, 0)

        self.ButtonExport = wx.Button(self, wx.ID_ANY, "E&xport...")
        sizer_6.Add(self.ButtonExport, 0, wx.LEFT, 10)

        sizer_6.Add((20, 20), 1, wx.EXPAND, 0)

        self.ButtonClose = wx.Button(self, wx.ID_CLOSE, "")
        sizer_6.Add(self.ButtonClose, 0, 0, 0)

        grid_sizer_1.AddGrowableRow(1)
        grid_sizer_1.AddGrowableCol(1)

        self.SetSizer(sizer_4)
        sizer_4.Fit(self)

        self.Layout()

        self.Bind(wx.EVT_COMBOBOX, self.on_isotope, self.ComboIsotope)
        self.Bind(wx.EVT_COMBOBOX, self.on_b0, self.ComboB0)
        self.Bind(wx.EVT_BUTTON, self.on_view, self.ButtonView)
        self.Bind(wx.EVT_BUTTON, self.on_delete, self.ButtonDelete)
        self.Bind(wx.EVT_BUTTON, self.on_import, self.ButtonImport)
        self.Bind(wx.EVT_BUTTON, self.on_export, self.ButtonExport)
        self.Bind(wx.EVT_BUTTON, self.on_close, self.ButtonClose)
        # end wxGlade

    def on_isotope(self, event):  # wxGlade: MyDialog.<event_handler>
        print("Event handler 'on_isotope' not implemented!")
        event.Skip()

    def on_b0(self, event):  # wxGlade: MyDialog.<event_handler>
        print("Event handler 'on_b0' not implemented!")
        event.Skip()

    def on_view(self, event):  # wxGlade: MyDialog.<event_handler>
        print("Event handler 'on_view' not implemented!")
        event.Skip()

    def on_delete(self, event):  # wxGlade: MyDialog.<event_handler>
        print("Event handler 'on_delete' not implemented!")
        event.Skip()

    def on_import(self, event):  # wxGlade: MyDialog.<event_handler>
        print("Event handler 'on_import' not implemented!")
        event.Skip()

    def on_export(self, event):  # wxGlade: MyDialog.<event_handler>
        print("Event handler 'on_export' not implemented!")
        event.Skip()

    def on_close(self, event):  # wxGlade: MyDialog.<event_handler>
        print("Event handler 'on_close' not implemented!")
        event.Skip()

# end of class MyDialog
