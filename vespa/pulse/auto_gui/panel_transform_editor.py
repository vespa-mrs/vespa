# -*- coding: UTF-8 -*-
#
# generated by wxGlade 1.0.5 on Mon May 15 22:17:33 2023
#

import wx

# begin wxGlade: dependencies
# end wxGlade

# begin wxGlade: extracode
# end wxGlade


class FrameTransformEditor(wx.Frame):
    def __init__(self, *args, **kwds):
        # begin wxGlade: FrameTransformEditor.__init__
        kwds["style"] = kwds.get("style", 0) | wx.DEFAULT_FRAME_STYLE
        wx.Frame.__init__(self, *args, **kwds)
        self.SetTitle("frame_transform_core")

        sizer_1 = wx.BoxSizer(wx.VERTICAL)

        self.PanelTransformKernelEditor = wx.Panel(self, wx.ID_ANY)
        sizer_1.Add(self.PanelTransformKernelEditor, 1, wx.EXPAND, 0)

        sizer_2 = wx.BoxSizer(wx.VERTICAL)

        sizer_global_parameters = wx.StaticBoxSizer(wx.StaticBox(self.PanelTransformKernelEditor, wx.ID_ANY, "Global Parameters"), wx.VERTICAL)
        sizer_2.Add(sizer_global_parameters, 0, wx.ALL | wx.EXPAND, 6)

        grid_sizer_8 = wx.FlexGridSizer(6, 4, 4, 4)
        sizer_global_parameters.Add(grid_sizer_8, 1, wx.EXPAND, 0)

        self.CheckFile1 = wx.CheckBox(self.PanelTransformKernelEditor, wx.ID_ANY, " Hide")
        self.CheckFile1.SetValue(1)
        grid_sizer_8.Add(self.CheckFile1, 0, wx.ALIGN_CENTER_VERTICAL, 0)

        label_1 = wx.StaticText(self.PanelTransformKernelEditor, wx.ID_ANY, " File1     Top Label: ", style=wx.ALIGN_RIGHT)
        grid_sizer_8.Add(label_1, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT, 0)

        self.TextFile1Label = wx.TextCtrl(self.PanelTransformKernelEditor, wx.ID_ANY, "")
        grid_sizer_8.Add(self.TextFile1Label, 0, wx.ALIGN_CENTER_VERTICAL | wx.EXPAND, 0)

        grid_sizer_8.Add((20, 20), 0, 0, 0)

        self.CheckFile2 = wx.CheckBox(self.PanelTransformKernelEditor, wx.ID_ANY, " Hide")
        self.CheckFile2.SetValue(1)
        grid_sizer_8.Add(self.CheckFile2, 0, wx.ALIGN_CENTER_VERTICAL, 0)

        label_6 = wx.StaticText(self.PanelTransformKernelEditor, wx.ID_ANY, " File2     Top Label: ", style=wx.ALIGN_RIGHT)
        grid_sizer_8.Add(label_6, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT, 0)

        self.TextFile2Label = wx.TextCtrl(self.PanelTransformKernelEditor, wx.ID_ANY, "")
        grid_sizer_8.Add(self.TextFile2Label, 0, wx.ALIGN_CENTER_VERTICAL | wx.EXPAND, 0)

        grid_sizer_8.Add((20, 20), 0, 0, 0)

        self.CheckTipAngle = wx.CheckBox(self.PanelTransformKernelEditor, wx.ID_ANY, " Hide")
        grid_sizer_8.Add(self.CheckTipAngle, 0, wx.ALIGN_CENTER_VERTICAL, 0)

        label_8 = wx.StaticText(self.PanelTransformKernelEditor, wx.ID_ANY, "Tip Angle [degrees]: ", style=wx.ALIGN_RIGHT)
        label_8.SetMinSize((102, -1))
        grid_sizer_8.Add(label_8, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT, 0)

        self.TextTipAngle = wx.TextCtrl(self.PanelTransformKernelEditor, wx.ID_ANY, "")
        grid_sizer_8.Add(self.TextTipAngle, 0, wx.ALIGN_CENTER_VERTICAL | wx.EXPAND, 0)

        label_3 = wx.StaticText(self.PanelTransformKernelEditor, wx.ID_ANY, "default")
        grid_sizer_8.Add(label_3, 0, wx.ALIGN_CENTER_VERTICAL, 0)

        self.CheckTimeSteps = wx.CheckBox(self.PanelTransformKernelEditor, wx.ID_ANY, " Hide")
        grid_sizer_8.Add(self.CheckTimeSteps, 0, wx.ALIGN_CENTER_VERTICAL, 0)

        label_9 = wx.StaticText(self.PanelTransformKernelEditor, wx.ID_ANY, "Time Steps [increments]: ", style=wx.ALIGN_RIGHT)
        label_9.SetMinSize((123, -1))
        grid_sizer_8.Add(label_9, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT, 0)

        self.TextTimeSteps = wx.TextCtrl(self.PanelTransformKernelEditor, wx.ID_ANY, "")
        grid_sizer_8.Add(self.TextTimeSteps, 0, wx.ALIGN_CENTER_VERTICAL | wx.EXPAND, 0)

        label_4 = wx.StaticText(self.PanelTransformKernelEditor, wx.ID_ANY, "default")
        grid_sizer_8.Add(label_4, 0, wx.ALIGN_CENTER_VERTICAL, 0)

        self.CheckDuration = wx.CheckBox(self.PanelTransformKernelEditor, wx.ID_ANY, " Hide")
        grid_sizer_8.Add(self.CheckDuration, 0, wx.ALIGN_CENTER_VERTICAL, 0)

        label_10 = wx.StaticText(self.PanelTransformKernelEditor, wx.ID_ANY, "Duration [milliseconds]: ", style=wx.ALIGN_RIGHT)
        label_10.SetMinSize((115, -1))
        grid_sizer_8.Add(label_10, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT, 0)

        self.TextDuration = wx.TextCtrl(self.PanelTransformKernelEditor, wx.ID_ANY, "")
        grid_sizer_8.Add(self.TextDuration, 0, wx.ALIGN_CENTER_VERTICAL | wx.EXPAND, 0)

        label_5 = wx.StaticText(self.PanelTransformKernelEditor, wx.ID_ANY, "default")
        grid_sizer_8.Add(label_5, 0, wx.ALIGN_CENTER_VERTICAL, 0)

        self.CheckBandwidth = wx.CheckBox(self.PanelTransformKernelEditor, wx.ID_ANY, " Hide")
        grid_sizer_8.Add(self.CheckBandwidth, 0, wx.ALIGN_CENTER_VERTICAL, 0)

        label_11 = wx.StaticText(self.PanelTransformKernelEditor, wx.ID_ANY, "Bandwidth [kilohertz]: ", style=wx.ALIGN_RIGHT)
        label_11.SetMinSize((109, -1))
        grid_sizer_8.Add(label_11, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT, 0)

        self.TextBandwidth = wx.TextCtrl(self.PanelTransformKernelEditor, wx.ID_ANY, "")
        grid_sizer_8.Add(self.TextBandwidth, 0, wx.ALIGN_CENTER_VERTICAL | wx.EXPAND, 0)

        label_12 = wx.StaticText(self.PanelTransformKernelEditor, wx.ID_ANY, "default")
        grid_sizer_8.Add(label_12, 0, wx.ALIGN_CENTER_VERTICAL, 0)

        sizer_user_parameters = wx.StaticBoxSizer(wx.StaticBox(self.PanelTransformKernelEditor, wx.ID_ANY, "User Defined Parameters"), wx.VERTICAL)
        sizer_2.Add(sizer_user_parameters, 1, wx.ALL | wx.EXPAND, 4)

        sizer_6 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_user_parameters.Add(sizer_6, 0, wx.BOTTOM | wx.EXPAND, 10)

        self.ButtonAddUserParameter = wx.Button(self.PanelTransformKernelEditor, wx.ID_ANY, "Add")
        sizer_6.Add(self.ButtonAddUserParameter, 0, wx.RIGHT, 5)

        self.ButtonRemoveUserParameter = wx.Button(self.PanelTransformKernelEditor, wx.ID_ANY, "Remove Selected")
        sizer_6.Add(self.ButtonRemoveUserParameter, 0, wx.LEFT, 5)

        grid_sizer_2 = wx.FlexGridSizer(1, 6, 5, 0)
        sizer_user_parameters.Add(grid_sizer_2, 1, wx.EXPAND, 0)

        self.InfoLabelParametersPlaceholder = wx.StaticText(self.PanelTransformKernelEditor, wx.ID_ANY, "LabelParametersPlaceholder")
        grid_sizer_2.Add(self.InfoLabelParametersPlaceholder, 0, 0, 0)

        sizer_comments = wx.StaticBoxSizer(wx.StaticBox(self.PanelTransformKernelEditor, wx.ID_ANY, "Comments"), wx.HORIZONTAL)
        sizer_2.Add(sizer_comments, 1, wx.ALL | wx.EXPAND, 4)

        self.TextComment = wx.TextCtrl(self.PanelTransformKernelEditor, wx.ID_ANY, "", style=wx.TE_MULTILINE)
        sizer_comments.Add(self.TextComment, 1, wx.EXPAND, 0)

        grid_sizer_2.AddGrowableCol(3)
        grid_sizer_2.AddGrowableCol(5)

        grid_sizer_8.AddGrowableCol(2)

        self.PanelTransformKernelEditor.SetSizer(sizer_2)

        self.SetSizer(sizer_1)
        sizer_1.Fit(self)

        self.Layout()

        self.Bind(wx.EVT_CHECKBOX, self.on_hide_file1, self.CheckFile1)
        self.Bind(wx.EVT_CHECKBOX, self.on_hide_file2, self.CheckFile2)
        self.Bind(wx.EVT_CHECKBOX, self.on_hide_tip_angle, self.CheckTipAngle)
        self.Bind(wx.EVT_CHECKBOX, self.on_hide_time_steps, self.CheckTimeSteps)
        self.Bind(wx.EVT_CHECKBOX, self.on_hide_duration, self.CheckDuration)
        self.Bind(wx.EVT_CHECKBOX, self.on_hide_bandwidth, self.CheckBandwidth)
        self.Bind(wx.EVT_BUTTON, self.on_add_user_parameter, self.ButtonAddUserParameter)
        self.Bind(wx.EVT_BUTTON, self.on_remove_user_parameter, self.ButtonRemoveUserParameter)
        # end wxGlade

    def on_hide_file1(self, event):  # wxGlade: FrameTransformEditor.<event_handler>
        print("Event handler 'on_hide_file1' not implemented!")
        event.Skip()

    def on_hide_file2(self, event):  # wxGlade: FrameTransformEditor.<event_handler>
        print("Event handler 'on_hide_file2' not implemented!")
        event.Skip()

    def on_hide_tip_angle(self, event):  # wxGlade: FrameTransformEditor.<event_handler>
        print("Event handler 'on_hide_tip_angle' not implemented!")
        event.Skip()

    def on_hide_time_steps(self, event):  # wxGlade: FrameTransformEditor.<event_handler>
        print("Event handler 'on_hide_time_steps' not implemented!")
        event.Skip()

    def on_hide_duration(self, event):  # wxGlade: FrameTransformEditor.<event_handler>
        print("Event handler 'on_hide_duration' not implemented!")
        event.Skip()

    def on_hide_bandwidth(self, event):  # wxGlade: FrameTransformEditor.<event_handler>
        print("Event handler 'on_hide_bandwidth' not implemented!")
        event.Skip()

    def on_add_user_parameter(self, event):  # wxGlade: FrameTransformEditor.<event_handler>
        print("Event handler 'on_add_user_parameter' not implemented!")
        event.Skip()

    def on_remove_user_parameter(self, event):  # wxGlade: FrameTransformEditor.<event_handler>
        print("Event handler 'on_remove_user_parameter' not implemented!")
        event.Skip()

# end of class FrameTransformEditor
