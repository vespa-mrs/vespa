# -*- coding: UTF-8 -*-
#
# generated by wxGlade 1.0.5 on Wed May 17 12:28:20 2023
#

import wx

# begin wxGlade: dependencies
from vespa.common.wx_gravy.widgets.floatspin_multiplier.floatspin_multiplier_base import FloatSpinMultiplier
from wx.lib.agw.floatspin import FloatSpin, EVT_FLOATSPIN, FS_LEFT, FS_RIGHT, FS_CENTRE, FS_READONLY
import wx.grid
# end wxGlade

# begin wxGlade: extracode
# end wxGlade


class DatasimUI(wx.Panel):
    def __init__(self, *args, **kwds):
        # begin wxGlade: DatasimUI.__init__
        kwds["style"] = kwds.get("style", 0) | wx.TAB_TRAVERSAL
        wx.Panel.__init__(self, *args, **kwds)

        self.SizerSplitterWindow = wx.BoxSizer(wx.VERTICAL)

        self.SplitterWindow = wx.SplitterWindow(self, wx.ID_ANY, style=wx.SP_3D | wx.SP_BORDER)
        self.SplitterWindow.SetMinimumPaneSize(20)
        self.SizerSplitterWindow.Add(self.SplitterWindow, 1, wx.EXPAND, 0)

        self.PanelDatasimControl = wx.Panel(self.SplitterWindow, wx.ID_ANY)

        sizer_39 = wx.BoxSizer(wx.VERTICAL)

        sizer_7 = wx.BoxSizer(wx.VERTICAL)
        sizer_39.Add(sizer_7, 0, wx.ALL | wx.EXPAND, 10)

        sizer_30 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_7.Add(sizer_30, 0, wx.ALL | wx.EXPAND, 4)

        label_11 = wx.StaticText(self.PanelDatasimControl, wx.ID_ANY, "Source:")
        sizer_30.Add(label_11, 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT | wx.RIGHT, 4)

        self.TextData = wx.TextCtrl(self.PanelDatasimControl, wx.ID_ANY, "", style=wx.TE_READONLY)
        sizer_30.Add(self.TextData, 1, wx.ALIGN_CENTER_VERTICAL, 0)

        label_4 = wx.StaticText(self.PanelDatasimControl, wx.ID_ANY, " Scale:")
        sizer_30.Add(label_4, 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT | wx.RIGHT, 4)

        self.FloatScale = FloatSpinMultiplier(self.PanelDatasimControl, wx.ID_ANY, value=0.0, digits=5, min_val=0.0, max_val=100.0, multiplier=1.1, agwStyle=FS_LEFT, style=wx.SP_ARROW_KEYS | wx.SP_WRAP | wx.TE_PROCESS_ENTER)
        sizer_30.Add(self.FloatScale, 0, 0, 0)

        sizer_8 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_7.Add(sizer_8, 0, wx.ALL | wx.EXPAND, 4)

        label_1 = wx.StaticText(self.PanelDatasimControl, wx.ID_ANY, "Loop1:")
        sizer_8.Add(label_1, 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT, 6)

        self.SpinIndex1 = wx.SpinCtrl(self.PanelDatasimControl, wx.ID_ANY, "1", min=0, max=100, style=0)
        sizer_8.Add(self.SpinIndex1, 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT, 3)

        label_100 = wx.StaticText(self.PanelDatasimControl, wx.ID_ANY, "Loop2:")
        sizer_8.Add(label_100, 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT, 6)

        self.SpinIndex2 = wx.SpinCtrl(self.PanelDatasimControl, wx.ID_ANY, "1", min=0, max=100, style=0)
        sizer_8.Add(self.SpinIndex2, 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT, 3)

        label_101 = wx.StaticText(self.PanelDatasimControl, wx.ID_ANY, "Loop3:")
        sizer_8.Add(label_101, 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT, 6)

        self.SpinIndex3 = wx.SpinCtrl(self.PanelDatasimControl, wx.ID_ANY, "1", min=0, max=100, style=0)
        sizer_8.Add(self.SpinIndex3, 0, wx.LEFT, 3)

        self.DatasimNotebook = wx.Notebook(self.PanelDatasimControl, wx.ID_ANY, style=wx.NB_BOTTOM)
        sizer_39.Add(self.DatasimNotebook, 1, wx.EXPAND, 0)

        self.PanelSpectralSettings = wx.Panel(self.DatasimNotebook, wx.ID_ANY)
        self.DatasimNotebook.AddPage(self.PanelSpectralSettings, "Spectral Settings")

        sizer_2_copy_1 = wx.BoxSizer(wx.VERTICAL)

        sizer_3 = wx.StaticBoxSizer(wx.StaticBox(self.PanelSpectralSettings, wx.ID_ANY, "Metabolite Spectral Settings"), wx.VERTICAL)
        sizer_2_copy_1.Add(sizer_3, 0, wx.EXPAND | wx.TOP, 8)

        self.ButtonSetSpectralResolution = wx.Button(self.PanelSpectralSettings, wx.ID_ANY, "Set metabolite spectral resolution and basis PPM range")
        sizer_3.Add(self.ButtonSetSpectralResolution, 0, wx.ALL, 8)

        grid_sizer_2 = wx.FlexGridSizer(4, 4, 4, 6)
        sizer_3.Add(grid_sizer_2, 0, wx.EXPAND, 0)

        labelTa = wx.StaticText(self.PanelSpectralSettings, wx.ID_ANY, "Ta [sec]:")
        grid_sizer_2.Add(labelTa, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT, 3)

        self.FloatTa = FloatSpin(self.PanelSpectralSettings, wx.ID_ANY, value=0.0, digits=3, min_val=0.0, max_val=1000000.0, increment=1.0, agwStyle=FS_LEFT, style=0)
        grid_sizer_2.Add(self.FloatTa, 0, wx.ALIGN_CENTER_VERTICAL, 0)

        label_17 = wx.StaticText(self.PanelSpectralSettings, wx.ID_ANY, "(for LW display only)")
        label_17.SetToolTip("Use individual metabolite Ta parameters to set the Voigt 'T2' value. This widget is only for calculating estimated LW for display below.")
        grid_sizer_2.Add(label_17, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT, 3)

        grid_sizer_2.Add((20, 20), 0, 0, 0)

        labelTb = wx.StaticText(self.PanelSpectralSettings, wx.ID_ANY, "Tb [sec]:")
        grid_sizer_2.Add(labelTb, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT, 3)

        self.FloatTb = FloatSpin(self.PanelSpectralSettings, wx.ID_ANY, value=0.0, digits=3, min_val=0.0, max_val=1000000.0, increment=1.0, agwStyle=FS_LEFT, style=0)
        grid_sizer_2.Add(self.FloatTb, 0, wx.ALIGN_CENTER_VERTICAL, 0)

        label_1_copy = wx.StaticText(self.PanelSpectralSettings, wx.ID_ANY, "Line Width [Hz]:")
        grid_sizer_2.Add(label_1_copy, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT, 3)

        self.TextLinewidth = wx.StaticText(self.PanelSpectralSettings, wx.ID_ANY, "")
        grid_sizer_2.Add(self.TextLinewidth, 0, wx.ALIGN_CENTER_VERTICAL, 0)

        label_5 = wx.StaticText(self.PanelSpectralSettings, wx.ID_ANY, "Phase 0 [deg]:")
        grid_sizer_2.Add(label_5, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT, 3)

        self.FloatPhase0 = FloatSpin(self.PanelSpectralSettings, wx.ID_ANY, value=0.0, digits=3, min_val=0.0, max_val=100.0, increment=1.0, agwStyle=FS_LEFT, style=0)
        grid_sizer_2.Add(self.FloatPhase0, 0, wx.ALIGN_CENTER_VERTICAL, 0)

        label_3 = wx.StaticText(self.PanelSpectralSettings, wx.ID_ANY, "B0 shift [Hz]:")
        grid_sizer_2.Add(label_3, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT, 3)

        self.FloatB0Shift = FloatSpin(self.PanelSpectralSettings, wx.ID_ANY, value=0.0, digits=3, min_val=0.0, max_val=100.0, increment=1.0, agwStyle=FS_LEFT, style=0)
        grid_sizer_2.Add(self.FloatB0Shift, 0, wx.ALIGN_CENTER_VERTICAL, 0)

        label_10 = wx.StaticText(self.PanelSpectralSettings, wx.ID_ANY, "Left Shift [pts]:")
        grid_sizer_2.Add(label_10, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT, 0)

        self.SpinLeftShift = wx.SpinCtrl(self.PanelSpectralSettings, wx.ID_ANY, "0", min=0, max=100, style=0)
        grid_sizer_2.Add(self.SpinLeftShift, 0, wx.ALIGN_CENTER_VERTICAL, 0)

        grid_sizer_2.Add((20, 20), 0, 0, 0)

        grid_sizer_2.Add((20, 20), 0, 0, 0)

        self.panel_1 = wx.Panel(self.PanelSpectralSettings, wx.ID_ANY)
        sizer_3.Add(self.panel_1, 0, wx.BOTTOM | wx.EXPAND | wx.TOP, 6)

        sizer_4 = wx.StaticBoxSizer(wx.StaticBox(self.panel_1, wx.ID_ANY, "Applied only on the plot, not on output"), wx.HORIZONTAL)

        label_8 = wx.StaticText(self.panel_1, wx.ID_ANY, "Phase 1 [deg]: ")
        sizer_4.Add(label_8, 0, wx.ALIGN_CENTER_VERTICAL, 0)

        self.FloatPhase1 = FloatSpin(self.panel_1, wx.ID_ANY, value=0.0, digits=3, min_val=0.0, max_val=100.0, increment=1.0, agwStyle=FS_LEFT, style=0)
        sizer_4.Add(self.FloatPhase1, 0, wx.ALIGN_CENTER_VERTICAL, 0)

        label_9 = wx.StaticText(self.panel_1, wx.ID_ANY, " Phase 1 Pivot [ppm]: ")
        sizer_4.Add(label_9, 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT, 8)

        self.FloatPhase1Pivot = FloatSpin(self.panel_1, wx.ID_ANY, value=0.0, digits=3, min_val=0.0, max_val=100.0, increment=1.0, agwStyle=FS_LEFT, style=0)
        sizer_4.Add(self.FloatPhase1Pivot, 0, wx.ALIGN_CENTER_VERTICAL, 0)

        sizer_12 = wx.StaticBoxSizer(wx.StaticBox(self.PanelSpectralSettings, wx.ID_ANY, "Datasim Noise Settings"), wx.VERTICAL)
        sizer_2_copy_1.Add(sizer_12, 0, wx.EXPAND | wx.TOP, 8)

        self.CheckDisplayNoiseInPlot = wx.CheckBox(self.PanelSpectralSettings, wx.ID_ANY, "  Display noise in plot")
        sizer_12.Add(self.CheckDisplayNoiseInPlot, 0, wx.BOTTOM | wx.TOP, 8)

        grid_sizer_4 = wx.FlexGridSizer(4, 4, 4, 4)
        sizer_12.Add(grid_sizer_4, 1, wx.EXPAND, 0)

        label_2 = wx.StaticText(self.PanelSpectralSettings, wx.ID_ANY, "Noise Ref Peak Area [spins]:")
        grid_sizer_4.Add(label_2, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT, 0)

        self.FloatRefPeakArea = FloatSpin(self.PanelSpectralSettings, wx.ID_ANY, value=0.0, digits=3, min_val=0.0, max_val=100.0, increment=1.0, agwStyle=FS_LEFT, style=0)
        grid_sizer_4.Add(self.FloatRefPeakArea, 0, 0, 0)

        label_6 = wx.StaticText(self.PanelSpectralSettings, wx.ID_ANY, "Effective SNR: ")
        grid_sizer_4.Add(label_6, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT | wx.LEFT, 8)

        self.TextEffectiveSnr = wx.StaticText(self.PanelSpectralSettings, wx.ID_ANY, "")
        grid_sizer_4.Add(self.TextEffectiveSnr, 0, wx.ALIGN_CENTER_VERTICAL, 0)

        label_23 = wx.StaticText(self.PanelSpectralSettings, wx.ID_ANY, "Noise Ref Peak Ta [sec]:")
        grid_sizer_4.Add(label_23, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT, 3)

        self.FloatRefPeakTa = FloatSpin(self.PanelSpectralSettings, wx.ID_ANY, value=0.0, digits=3, min_val=0.0, max_val=100.0, increment=1.0, agwStyle=FS_LEFT, style=0)
        grid_sizer_4.Add(self.FloatRefPeakTa, 0, wx.ALIGN_CENTER_VERTICAL, 0)

        label_36 = wx.StaticText(self.PanelSpectralSettings, wx.ID_ANY, "Effective Linewidth: ")
        grid_sizer_4.Add(label_36, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT | wx.LEFT, 8)

        self.TextEffectiveLinewidth = wx.StaticText(self.PanelSpectralSettings, wx.ID_ANY, "")
        grid_sizer_4.Add(self.TextEffectiveLinewidth, 0, wx.ALIGN_CENTER_VERTICAL, 0)

        label_23_copy = wx.StaticText(self.PanelSpectralSettings, wx.ID_ANY, "Noise Ref Peak Tb [sec]:")
        grid_sizer_4.Add(label_23_copy, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT, 3)

        self.FloatRefPeakTb = FloatSpin(self.PanelSpectralSettings, wx.ID_ANY, value=0.0, digits=3, min_val=0.0, max_val=100.0, increment=1.0, agwStyle=FS_LEFT, style=0)
        grid_sizer_4.Add(self.FloatRefPeakTb, 0, wx.ALIGN_CENTER_VERTICAL, 0)

        grid_sizer_4.Add((20, 20), 0, 0, 0)

        grid_sizer_4.Add((20, 20), 0, 0, 0)

        label_22 = wx.StaticText(self.PanelSpectralSettings, wx.ID_ANY, "Noise RMS Multiplier [0-1000]:")
        grid_sizer_4.Add(label_22, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT, 3)

        self.FloatNoiseRmsMultiplier = FloatSpin(self.PanelSpectralSettings, wx.ID_ANY, value=0.0, digits=3, min_val=0.0, max_val=100.0, increment=1.0, agwStyle=FS_LEFT, style=0)
        grid_sizer_4.Add(self.FloatNoiseRmsMultiplier, 0, 0, 0)

        grid_sizer_4.Add((20, 20), 0, 0, 0)

        grid_sizer_4.Add((20, 20), 0, 0, 0)

        sizer_6 = wx.StaticBoxSizer(wx.StaticBox(self.PanelSpectralSettings, wx.ID_ANY, "Monte Carlo Dataset Settings"), wx.VERTICAL)
        sizer_2_copy_1.Add(sizer_6, 0, wx.EXPAND | wx.TOP, 8)

        grid_sizer_3 = wx.FlexGridSizer(1, 2, 4, 4)
        sizer_6.Add(grid_sizer_3, 1, wx.EXPAND, 0)

        label_7 = wx.StaticText(self.PanelSpectralSettings, wx.ID_ANY, "Number of voxels in dataset [eg. 10 or 1000]:")
        grid_sizer_3.Add(label_7, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT, 3)

        self.SpinMonteCarloVoxels = wx.SpinCtrl(self.PanelSpectralSettings, wx.ID_ANY, "10", min=1, max=100, style=0)
        grid_sizer_3.Add(self.SpinMonteCarloVoxels, 0, wx.ALIGN_CENTER_VERTICAL, 0)

        sizer_13 = wx.StaticBoxSizer(wx.StaticBox(self.PanelSpectralSettings, wx.ID_ANY, "Datasim Creation Comments"), wx.VERTICAL)
        sizer_2_copy_1.Add(sizer_13, 1, wx.EXPAND | wx.TOP, 8)

        self.TextComment = wx.TextCtrl(self.PanelSpectralSettings, wx.ID_ANY, "", style=wx.TE_BESTWRAP | wx.TE_MULTILINE)
        self.TextComment.SetToolTip("Add description of simulated data creation here.")
        sizer_13.Add(self.TextComment, 1, wx.EXPAND, 0)

        self.PanelMetaboliteSignals = wx.Panel(self.DatasimNotebook, wx.ID_ANY)
        self.DatasimNotebook.AddPage(self.PanelMetaboliteSignals, "Metabolite Signals")

        sizer_24 = wx.BoxSizer(wx.VERTICAL)

        sizer_15 = wx.BoxSizer(wx.VERTICAL)
        sizer_24.Add(sizer_15, 1, wx.EXPAND, 0)

        self.PanelMetabolite = wx.Panel(self.PanelMetaboliteSignals, wx.ID_ANY)
        sizer_15.Add(self.PanelMetabolite, 1, wx.EXPAND, 0)

        sizer_14 = wx.BoxSizer(wx.VERTICAL)

        sizer_16 = wx.StaticBoxSizer(wx.StaticBox(self.PanelMetabolite, wx.ID_ANY, "Metabolite Signals"), wx.VERTICAL)
        sizer_14.Add(sizer_16, 1, wx.EXPAND | wx.TOP, 8)

        self.PanelMetaboliteLines = wx.ScrolledWindow(self.PanelMetabolite, wx.ID_ANY, style=wx.BORDER_SIMPLE | wx.TAB_TRAVERSAL)
        self.PanelMetaboliteLines.SetBackgroundColour(wx.Colour(255, 255, 255))
        self.PanelMetaboliteLines.SetScrollRate(10, 10)
        sizer_16.Add(self.PanelMetaboliteLines, 1, wx.EXPAND, 0)

        MetaboliteListGridSizer = wx.FlexGridSizer(1, 3, 5, 10)

        self.LabelMetabolites = wx.StaticText(self.PanelMetaboliteLines, wx.ID_ANY, "placeholder")
        MetaboliteListGridSizer.Add(self.LabelMetabolites, 0, wx.ALL, 2)

        MetaboliteListGridSizer.Add((0, 0), 0, 0, 0)

        MetaboliteListGridSizer.Add((0, 0), 0, 0, 0)

        sizer_5 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_16.Add(sizer_5, 0, wx.ALL | wx.EXPAND, 2)

        self.ButtonSelectAll = wx.Button(self.PanelMetabolite, wx.ID_ANY, "Select All")
        sizer_5.Add(self.ButtonSelectAll, 0, wx.RIGHT, 5)

        self.ButtonSelectNone = wx.Button(self.PanelMetabolite, wx.ID_ANY, "Select None")
        sizer_5.Add(self.ButtonSelectNone, 0, wx.RIGHT, 2)

        self.PanelMacromoleculeSignals = wx.Panel(self.DatasimNotebook, wx.ID_ANY)
        self.DatasimNotebook.AddPage(self.PanelMacromoleculeSignals, "Macromolecule Signals")

        sizer_25 = wx.BoxSizer(wx.VERTICAL)

        sizer_10 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_25.Add(sizer_10, 0, wx.EXPAND, 0)

        self.RadioMmolLineshape = wx.RadioBox(self.PanelMacromoleculeSignals, wx.ID_ANY, "Macromolecule Lineshape", choices=["Lorentzian", "Gaussian"], majorDimension=1, style=wx.RA_SPECIFY_ROWS)
        self.RadioMmolLineshape.SetSelection(0)
        sizer_10.Add(self.RadioMmolLineshape, 0, wx.ALL | wx.EXPAND, 1)

        sizer_17 = wx.StaticBoxSizer(wx.StaticBox(self.PanelMacromoleculeSignals, wx.ID_ANY, "Group Scale"), wx.HORIZONTAL)
        sizer_10.Add(sizer_17, 0, wx.EXPAND, 1)

        self.FloatMmolGroupScale = FloatSpin(self.PanelMacromoleculeSignals, wx.ID_ANY, value=1.0, digits=3, min_val=0.0, max_val=100000.0, increment=0.1, agwStyle=FS_LEFT)
        sizer_17.Add(self.FloatMmolGroupScale, 0, wx.EXPAND, 0)

        sizer_10.Add((20, 20), 1, 0, 0)

        sizer_11 = wx.StaticBoxSizer(wx.StaticBox(self.PanelMacromoleculeSignals, wx.ID_ANY, "Macromolecule Presets"), wx.VERTICAL)
        sizer_10.Add(sizer_11, 1, wx.EXPAND, 0)

        self.ChoicePresetMmol = wx.Choice(self.PanelMacromoleculeSignals, wx.ID_ANY, choices=["choice 1"])
        self.ChoicePresetMmol.SetSelection(0)
        sizer_11.Add(self.ChoicePresetMmol, 0, wx.ALL | wx.EXPAND, 1)

        sizer_19 = wx.StaticBoxSizer(wx.StaticBox(self.PanelMacromoleculeSignals, wx.ID_ANY, "Macromolecule Signals"), wx.VERTICAL)
        sizer_25.Add(sizer_19, 1, wx.EXPAND | wx.TOP, 8)

        self.ScrolledWindowMmol = wx.ScrolledWindow(self.PanelMacromoleculeSignals, wx.ID_ANY, style=wx.TAB_TRAVERSAL)
        self.ScrolledWindowMmol.SetScrollRate(10, 10)
        sizer_19.Add(self.ScrolledWindowMmol, 1, wx.BOTTOM | wx.EXPAND | wx.TOP, 4)

        sizer_2 = wx.BoxSizer(wx.VERTICAL)

        self.GridMmol = wx.grid.Grid(self.ScrolledWindowMmol, wx.ID_ANY, size=(1, 1))
        self.GridMmol.CreateGrid(0, 0)
        sizer_2.Add(self.GridMmol, 1, wx.EXPAND, 0)

        sizer_28 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_19.Add(sizer_28, 0, wx.ALL | wx.EXPAND, 2)

        self.ButtonMmolSelectAll = wx.Button(self.PanelMacromoleculeSignals, wx.ID_ANY, "Select All")
        sizer_28.Add(self.ButtonMmolSelectAll, 0, wx.RIGHT, 5)

        self.ButtonMmolSelectNone = wx.Button(self.PanelMacromoleculeSignals, wx.ID_ANY, "Select None")
        sizer_28.Add(self.ButtonMmolSelectNone, 0, wx.RIGHT, 5)

        self.ButtonMmolDeleteLine = wx.Button(self.PanelMacromoleculeSignals, wx.ID_ANY, "Delete Selected")
        sizer_28.Add(self.ButtonMmolDeleteLine, 0, wx.RIGHT, 5)

        self.ButtonMmolAddLine = wx.Button(self.PanelMacromoleculeSignals, wx.ID_ANY, "Add Line")
        sizer_28.Add(self.ButtonMmolAddLine, 0, wx.RIGHT, 0)

        sizer_28.Add((20, 20), 1, 0, 0)

        self.ButtonMmolImportHlsvd = wx.Button(self.PanelMacromoleculeSignals, wx.ID_ANY, "Import HLSVD File")
        sizer_28.Add(self.ButtonMmolImportHlsvd, 0, 0, 0)

        self.PanelBaselineSignals = wx.Panel(self.DatasimNotebook, wx.ID_ANY)
        self.DatasimNotebook.AddPage(self.PanelBaselineSignals, "Baseline Signals")

        sizer_38 = wx.BoxSizer(wx.VERTICAL)

        sizer_18 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_38.Add(sizer_18, 0, wx.EXPAND, 0)

        self.RadioBaseLineshape = wx.RadioBox(self.PanelBaselineSignals, wx.ID_ANY, "Baseline Lineshape", choices=["Lorentzian", "Gaussian"], majorDimension=1, style=wx.RA_SPECIFY_ROWS)
        self.RadioBaseLineshape.SetSelection(0)
        sizer_18.Add(self.RadioBaseLineshape, 0, wx.ALL | wx.EXPAND, 1)

        sizer_20 = wx.StaticBoxSizer(wx.StaticBox(self.PanelBaselineSignals, wx.ID_ANY, "Group Scale"), wx.HORIZONTAL)
        sizer_18.Add(sizer_20, 0, wx.EXPAND, 1)

        self.FloatBaseGroupScale = FloatSpin(self.PanelBaselineSignals, wx.ID_ANY, value=1.0, digits=3, min_val=0.0, max_val=100000.0, increment=0.1, agwStyle=FS_LEFT)
        sizer_20.Add(self.FloatBaseGroupScale, 0, wx.EXPAND, 0)

        sizer_18.Add((20, 20), 1, 0, 0)

        sizer_26 = wx.StaticBoxSizer(wx.StaticBox(self.PanelBaselineSignals, wx.ID_ANY, "Other Baseline Signals"), wx.VERTICAL)
        sizer_38.Add(sizer_26, 1, wx.EXPAND | wx.TOP, 8)

        self.ScrolledWindowBase = wx.ScrolledWindow(self.PanelBaselineSignals, wx.ID_ANY, style=wx.TAB_TRAVERSAL)
        self.ScrolledWindowBase.SetScrollRate(10, 10)
        sizer_26.Add(self.ScrolledWindowBase, 1, wx.BOTTOM | wx.EXPAND | wx.TOP, 4)

        sizer_9 = wx.BoxSizer(wx.VERTICAL)

        self.GridBase = wx.grid.Grid(self.ScrolledWindowBase, wx.ID_ANY, size=(1, 1))
        self.GridBase.CreateGrid(0, 0)
        sizer_9.Add(self.GridBase, 1, wx.EXPAND, 0)

        sizer_27 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_26.Add(sizer_27, 0, wx.ALL | wx.EXPAND, 2)

        self.ButtonBaseSelectAll = wx.Button(self.PanelBaselineSignals, wx.ID_ANY, "Select All")
        sizer_27.Add(self.ButtonBaseSelectAll, 0, wx.RIGHT, 5)

        self.ButtonBaseSelectNone = wx.Button(self.PanelBaselineSignals, wx.ID_ANY, "Select None")
        sizer_27.Add(self.ButtonBaseSelectNone, 0, wx.RIGHT, 5)

        self.ButtonBaseDeleteLine = wx.Button(self.PanelBaselineSignals, wx.ID_ANY, "Delete Selected")
        sizer_27.Add(self.ButtonBaseDeleteLine, 0, wx.RIGHT, 5)

        self.ButtonBaseAddLine = wx.Button(self.PanelBaselineSignals, wx.ID_ANY, "Add Line")
        sizer_27.Add(self.ButtonBaseAddLine, 0, wx.RIGHT, 0)

        sizer_27.Add((20, 20), 1, 0, 0)

        self.ButtonBaseImportHlsvd = wx.Button(self.PanelBaselineSignals, wx.ID_ANY, "Import HLSVD File")
        sizer_27.Add(self.ButtonBaseImportHlsvd, 0, 0, 0)

        self.PanelDatasimPlot = wx.Panel(self.SplitterWindow, wx.ID_ANY)

        self.ScrolledWindowBase.SetSizer(sizer_9)

        self.PanelBaselineSignals.SetSizer(sizer_38)

        self.ScrolledWindowMmol.SetSizer(sizer_2)

        self.PanelMacromoleculeSignals.SetSizer(sizer_25)

        self.PanelMetaboliteLines.SetSizer(MetaboliteListGridSizer)

        self.PanelMetabolite.SetSizer(sizer_14)

        self.PanelMetaboliteSignals.SetSizer(sizer_24)

        self.panel_1.SetSizer(sizer_4)

        self.PanelSpectralSettings.SetSizer(sizer_2_copy_1)

        self.PanelDatasimControl.SetSizer(sizer_39)

        self.SplitterWindow.SplitVertically(self.PanelDatasimControl, self.PanelDatasimPlot)

        self.SetSizer(self.SizerSplitterWindow)

        self.Layout()

        self.Bind( EVT_FLOATSPIN, self.on_scale, self.FloatScale)
        self.Bind(wx.EVT_SPINCTRL, self.on_index1, self.SpinIndex1)
        self.Bind(wx.EVT_SPINCTRL, self.on_index2, self.SpinIndex2)
        self.Bind(wx.EVT_SPINCTRL, self.on_index3, self.SpinIndex3)
        self.Bind(wx.EVT_BUTTON, self.on_set_spectral_resolution, self.ButtonSetSpectralResolution)
        self.Bind( EVT_FLOATSPIN, self.on_ta, self.FloatTa)
        self.Bind( EVT_FLOATSPIN, self.on_tb, self.FloatTb)
        self.Bind( EVT_FLOATSPIN, self.on_phase0, self.FloatPhase0)
        self.Bind( EVT_FLOATSPIN, self.on_b0_shift, self.FloatB0Shift)
        self.Bind(wx.EVT_SPINCTRL, self.on_left_shift, self.SpinLeftShift)
        self.Bind( EVT_FLOATSPIN, self.on_phase1, self.FloatPhase1)
        self.Bind( EVT_FLOATSPIN, self.on_phase_1_pivot, self.FloatPhase1Pivot)
        self.Bind(wx.EVT_CHECKBOX, self.on_display_noise_in_plot, self.CheckDisplayNoiseInPlot)
        self.Bind( EVT_FLOATSPIN, self.on_ref_peak_area, self.FloatRefPeakArea)
        self.Bind( EVT_FLOATSPIN, self.on_ref_peak_ta, self.FloatRefPeakTa)
        self.Bind( EVT_FLOATSPIN, self.on_ref_peak_tb, self.FloatRefPeakTb)
        self.Bind( EVT_FLOATSPIN, self.on_noise_rms_multiplier, self.FloatNoiseRmsMultiplier)
        self.Bind(wx.EVT_SPINCTRL, self.on_matrix_size, self.SpinMonteCarloVoxels)
        self.Bind(wx.EVT_TEXT, self.on_comment, self.TextComment)
        self.Bind(wx.EVT_TEXT_ENTER, self.on_comment, self.TextComment)
        self.Bind(wx.EVT_BUTTON, self.on_select_all, self.ButtonSelectAll)
        self.Bind(wx.EVT_BUTTON, self.on_select_none, self.ButtonSelectNone)
        self.Bind(wx.EVT_RADIOBOX, self.on_mmol_lineshape, self.RadioMmolLineshape)
        self.Bind( EVT_FLOATSPIN, self.on_mmol_group_scale, self.FloatMmolGroupScale)
        self.Bind(wx.EVT_CHOICE, self.on_preset_mmol, self.ChoicePresetMmol)
        self.Bind(wx.EVT_BUTTON, self.on_mmol_select_all, self.ButtonMmolSelectAll)
        self.Bind(wx.EVT_BUTTON, self.on_mmol_select_none, self.ButtonMmolSelectNone)
        self.Bind(wx.EVT_BUTTON, self.on_mmol_delete_selected, self.ButtonMmolDeleteLine)
        self.Bind(wx.EVT_BUTTON, self.on_mmol_add_line, self.ButtonMmolAddLine)
        self.Bind(wx.EVT_BUTTON, self.on_mmol_import_hlsvd_file, self.ButtonMmolImportHlsvd)
        self.Bind(wx.EVT_RADIOBOX, self.on_base_lineshape, self.RadioBaseLineshape)
        self.Bind( EVT_FLOATSPIN, self.on_base_group_scale, self.FloatBaseGroupScale)
        self.Bind(wx.EVT_BUTTON, self.on_base_select_all, self.ButtonBaseSelectAll)
        self.Bind(wx.EVT_BUTTON, self.on_base_select_none, self.ButtonBaseSelectNone)
        self.Bind(wx.EVT_BUTTON, self.on_base_delete_selected, self.ButtonBaseDeleteLine)
        self.Bind(wx.EVT_BUTTON, self.on_base_add_line, self.ButtonBaseAddLine)
        self.Bind(wx.EVT_BUTTON, self.on_base_import_hlsvd_file, self.ButtonBaseImportHlsvd)
        # end wxGlade

    def on_scale(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_scale' not implemented!")
        event.Skip()

    def on_index1(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_index1' not implemented!")
        event.Skip()

    def on_index2(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_index2' not implemented!")
        event.Skip()

    def on_index3(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_index3' not implemented!")
        event.Skip()

    def on_set_spectral_resolution(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_set_spectral_resolution' not implemented!")
        event.Skip()

    def on_ta(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_ta' not implemented!")
        event.Skip()

    def on_tb(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_tb' not implemented!")
        event.Skip()

    def on_phase0(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_phase0' not implemented!")
        event.Skip()

    def on_b0_shift(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_b0_shift' not implemented!")
        event.Skip()

    def on_left_shift(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_left_shift' not implemented!")
        event.Skip()

    def on_phase1(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_phase1' not implemented!")
        event.Skip()

    def on_phase_1_pivot(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_phase_1_pivot' not implemented!")
        event.Skip()

    def on_display_noise_in_plot(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_display_noise_in_plot' not implemented!")
        event.Skip()

    def on_ref_peak_area(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_ref_peak_area' not implemented!")
        event.Skip()

    def on_ref_peak_ta(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_ref_peak_ta' not implemented!")
        event.Skip()

    def on_ref_peak_tb(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_ref_peak_tb' not implemented!")
        event.Skip()

    def on_noise_rms_multiplier(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_noise_rms_multiplier' not implemented!")
        event.Skip()

    def on_matrix_size(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_matrix_size' not implemented!")
        event.Skip()

    def on_comment(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_comment' not implemented!")
        event.Skip()

    def on_select_all(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_select_all' not implemented!")
        event.Skip()

    def on_select_none(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_select_none' not implemented!")
        event.Skip()

    def on_mmol_lineshape(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_mmol_lineshape' not implemented!")
        event.Skip()

    def on_mmol_group_scale(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_mmol_group_scale' not implemented!")
        event.Skip()

    def on_preset_mmol(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_preset_mmol' not implemented!")
        event.Skip()

    def on_mmol_select_all(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_mmol_select_all' not implemented!")
        event.Skip()

    def on_mmol_select_none(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_mmol_select_none' not implemented!")
        event.Skip()

    def on_mmol_delete_selected(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_mmol_delete_selected' not implemented!")
        event.Skip()

    def on_mmol_add_line(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_mmol_add_line' not implemented!")
        event.Skip()

    def on_mmol_import_hlsvd_file(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_mmol_import_hlsvd_file' not implemented!")
        event.Skip()

    def on_base_lineshape(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_base_lineshape' not implemented!")
        event.Skip()

    def on_base_group_scale(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_base_group_scale' not implemented!")
        event.Skip()

    def on_base_select_all(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_base_select_all' not implemented!")
        event.Skip()

    def on_base_select_none(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_base_select_none' not implemented!")
        event.Skip()

    def on_base_delete_selected(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_base_delete_selected' not implemented!")
        event.Skip()

    def on_base_add_line(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_base_add_line' not implemented!")
        event.Skip()

    def on_base_import_hlsvd_file(self, event):  # wxGlade: DatasimUI.<event_handler>
        print("Event handler 'on_base_import_hlsvd_file' not implemented!")
        event.Skip()

# end of class DatasimUI

class MyFrame(wx.Frame):
    def __init__(self, *args, **kwds):
        # begin wxGlade: MyFrame.__init__
        kwds["style"] = kwds.get("style", 0) | wx.DEFAULT_FRAME_STYLE
        wx.Frame.__init__(self, *args, **kwds)
        self.SetSize((973, 844))
        self.SetTitle("frame_1")

        sizer_1 = wx.BoxSizer(wx.VERTICAL)

        self.DatasimUI = DatasimUI(self, wx.ID_ANY)
        sizer_1.Add(self.DatasimUI, 1, wx.EXPAND, 0)

        self.SetSizer(sizer_1)

        self.Layout()
        # end wxGlade

# end of class MyFrame
