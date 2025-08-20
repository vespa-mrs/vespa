# Python modules
import os
import sys
import importlib
import xml.etree.cElementTree as ElementTree
from xml.etree.cElementTree import Element

# 3rd party modules
import wx
import wx.grid as gridlib
import numpy as np
from pubsub import pub as pubsub
from wx.lib.mixins.listctrl import CheckListCtrlMixin, ColumnSorterMixin

# Our modules
import vespa.analysis.mrs_dataset as mrs_dataset
import vespa.analysis.tab_base as tab_base
import vespa.analysis.constants as constants
import vespa.analysis.util_menu as util_menu
import vespa.analysis.prefs as prefs_module
import vespa.analysis.util_analysis_config as util_analysis_config
import vespa.analysis.dialog_dataset_browser as dialog_dataset_browser
import vespa.analysis.auto_gui.spectral as spectral
import vespa.analysis.functors.funct_ecc as funct_ecc
import vespa.analysis.functors.funct_water_filter as funct_watfilt

from vespa.analysis.plot_panel_spectral import PlotPanelSpectral
from vespa.analysis.plot_panel_svd_filter import PlotPanelSvdFilter

import vespa.common.wx_gravy.util as wx_util
import vespa.common.constants as common_constants
import vespa.common.wx_gravy.common_dialogs as common_dialogs
import vespa.common.util.misc as util_misc
import vespa.common.util.ppm as util_ppm
import vespa.common.util.fileio as util_fileio
import vespa.common.util.export as util_export
import vespa.common.util.xml_ as util_xml
import vespa.common.util.time_ as util_time
import vespa.common.mrs_data_raw as mrs_data_raw

from vespa.common.constants import DEGREES_TO_RADIANS


#------------------------------------------------------------------------------
HLSVD_MAX_SINGULAR_VALUES = 256
HLSVD_MIN_DATA_POINTS = 128
_HLSVD_RESULTS_DISPLAY_SIZE = 6


class CheckListCtrl(wx.ListCtrl, ColumnSorterMixin):
    def __init__(self, _inner_notebook, tab):
        style = wx.LC_REPORT | wx.LC_HRULES | wx.LC_VRULES
        wx.ListCtrl.__init__(self, _inner_notebook, -1, style=style)
        ColumnSorterMixin.__init__(self, _HLSVD_RESULTS_DISPLAY_SIZE)
        self.itemDataMap = {}
        self._tab_dataset = _inner_notebook
        self.tab = tab
        self.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnItemActivated)
        self.Bind(wx.EVT_LIST_ITEM_CHECKED, self.OnCheckItem)
        self.Bind(wx.EVT_LIST_ITEM_UNCHECKED, self.OnCheckItem)

    def GetListCtrl(self):
        return self

    def OnItemActivated(self, evt):
        item_index = evt.GetIndex()
        if self.IsItemChecked(item_index):
            self.CheckItem(item_index, False)
        else:
            self.CheckItem(item_index, True)

    # this is called by the base class when an item is checked/unchecked
    def OnCheckItem(self, evt):
        self.tab.on_check_item(self, evt.Index, evt.Label, self.IsItemChecked(evt.Index))


#------------------------------------------------------------------------------

class CustomDataTable(gridlib.GridTableBase):

    def __init__(self):
        gridlib.GridTableBase.__init__(self)

        self.colLabels = ['','Rank', 'PPM', 'Freq', 'Damping', 'Phase', 'Area']

        self.dataTypes = [gridlib.GRID_VALUE_BOOL,
                          gridlib.GRID_VALUE_NUMBER,
                          gridlib.GRID_VALUE_FLOAT + ':6,2',
                          gridlib.GRID_VALUE_FLOAT + ':6,2',
                          gridlib.GRID_VALUE_FLOAT + ':6,2',
                          gridlib.GRID_VALUE_FLOAT + ':6,2',
                          gridlib.GRID_VALUE_FLOAT + ':6,2', ]

        # self.colWidths = [60,60,80,80,80,80,80]
        #
        # self.list_svd_results.SetColumnWidth(0, 60)
        # self.list_svd_results.SetColumnWidth(1, 60)
        # self.list_svd_results.SetColumnWidth(2, 80)
        # self.list_svd_results.SetColumnWidth(3, 80)
        # self.list_svd_results.SetColumnWidth(4, 80)
        # self.list_svd_results.SetColumnWidth(5, 80)
        # self.list_svd_results.SetColumnWidth(6, 80)

        self.data = [
            [0, 1, 1.11, 1.12, 1.13, 1.14, 1.15],
            [0, 2, 1.21, 1.22, 1.23, 1.24, 1.25],
            [0, 3, 1.31, 1.32, 1.33, 1.34, 1.35]

            ]


    #--------------------------------------------------
    # required methods for the wxPyGridTableBase interface

    def GetNumberRows(self):
        return len(self.data) + 1

    def GetNumberCols(self):
        return len(self.data[0])

    def IsEmptyCell(self, row, col):
        try:
            return not self.data[row][col]
        except IndexError:
            return True

    # Get/Set values in the table.  The Python version of these
    # methods can handle any data-type, (as long as the Editor and
    # Renderer understands the type too,) not just strings as in the
    # C++ version.
    def GetValue(self, row, col):
        try:
            return self.data[row][col]
        except IndexError:
            return ''

    def SetValue(self, row, col, value):
        def innerSetValue(row, col, value):
            try:
                self.data[row][col] = value
            except IndexError:
                # add a new row
                self.data.append([''] * self.GetNumberCols())
                innerSetValue(row, col, value)

                # tell the grid we've added a row
                msg = gridlib.GridTableMessage(self,            # The table
                        gridlib.GRIDTABLE_NOTIFY_ROWS_APPENDED, # what we did to it
                        1   )                                   # how many

                self.GetView().ProcessTableMessage(msg)
        innerSetValue(row, col, value)

    #--------------------------------------------------
    # Some optional methods

    # Called when the grid needs to display labels
    def GetColLabelValue(self, col):
        return self.colLabels[col]

    # Called to determine the kind of editor/renderer to use by
    # default, doesn't necessarily have to be the same type used
    # natively by the editor/renderer if they know how to convert.
    def GetTypeName(self, row, col):
        return self.dataTypes[col]

    # Called to determine how the data can be fetched and stored by the
    # editor and renderer.  This allows you to enforce some type-safety
    # in the grid.
    def CanGetValueAs(self, row, col, typeName):
        colType = self.dataTypes[col].split(':')[0]
        if typeName == colType:
            return True
        else:
            return False

    def CanSetValueAs(self, row, col, typeName):
        return self.CanGetValueAs(row, col, typeName)


#---------------------------------------------------------------------------

class CustTableGrid(gridlib.Grid):

    def __init__(self, parent):

        gridlib.Grid.__init__(self, parent, -1)
        table = CustomDataTable()

        # The second parameter means that the grid is to take ownership of the
        # table and will destroy it when done.  Otherwise you would need to keep
        # a reference to it and call it's Destroy method later.
        self.SetTable(table, True)
        self.SetRowLabelSize(0)
        self.SetMargins(0,0)
        self.AutoSizeColumns(False)

        self.Bind(gridlib.EVT_GRID_CELL_LEFT_DCLICK, self.OnLeftDClick)


    # I do this because I don't like the default behaviour of not starting the
    # cell editor on double clicks, but only a second click.
    def OnLeftDClick(self, evt):
        if self.CanEnableCellControl():
            self.EnableCellEditControl()

#------------------------------------------------------------------------------

def _configure_combo(control, choices, selection=''):
        lines = list(choices.values())
        control.SetItems(lines)
        if selection in lines:
            control.SetStringSelection(selection)
        else:
            control.SetStringSelection(lines[0])

def _paired_event(obj_min, obj_max):
        val_min = obj_min.GetValue()
        val_max = obj_max.GetValue()
        pmin = min(val_min, val_max)
        pmax = max(val_min, val_max)
        obj_min.SetValue(pmin)
        obj_max.SetValue(pmax)
        return pmin, pmax    


#------------------------------------------------------------------------------

def derived_dataset(ds, dsB, view, _tab_dataset, top, mode):

    if mode in ['plotc_a_plus_b', 'plotc_a_minus_b']:

        # Create new tab label from current - check for uniqueness in Tab label list
        tab_name = _tab_dataset.indexAB[0]
        label_base = tab_name.split('.')
        label_new = label_base[0]+'.Sum' if mode == 'plotc_a_plus_b' else label_base[0]+'.Dif'

        tab_labels = list(top.datasets.keys())
        if label_new in tab_labels:
            msg = 'Warning (derived_dataset) - New Tab/Dataset has same label as existing Tab, adding Suffix.'
            r = common_dialogs.message(msg, style=common_dialogs.I_OK)
            for s in range(2,10):
                ss = str(s)
                label_new = label_base[0]+'.Sum'+ss if mode == 'plotc_a_plus_b' else label_base[0]+'.Dif'+ss
                if not label_new in tab_labels:
                    break

        # Grab the A and B data, specify the derived data

        stamp = util_time.now(util_time.ISO_TIMESTAMP_FORMAT).split('T')
        if mode == 'plotc_a_plus_b':
            labl_func = 'PlotA + PlotB'
            dat = view.get_data_phased(0)+view.get_data_phased(1)
        elif mode == 'plotc_a_minus_b':
            labl_func = 'PlotA - PlotB'
            dat = view.get_data_phased(0)-view.get_data_phased(1)
        dat = dat.ravel()
        dat = np.fft.ifft(np.fft.fftshift(dat))*len(dat)
        dat.shape = 1, 1, 1, len(dat)

        # Create informative header - re. provenance of this spectrum

        lines = ['Derived Dataset from Analysis Spectral General Tab - PlotC ('+labl_func+')']
        lines.append('------------------------------------------------------------------')
        lines.append('An IFFT() was performed on the Plot C spectrum and the resultant FID written into')
        lines.append('this Dataset object. All other parameters were taken from the parent Dataset')

        lines.append(' ')
        lines.append('Creation_date            - '+stamp[0])
        lines.append('Creation_time            - '+stamp[1])
        lines.append(' ')
        lines.append('Filename Dataset A - '+str(ds.blocks['raw'].data_sources[0]))
        lines.append('UUID Dataset A - '+str(ds.id))
        lines.append('Filename Dataset B - '+str(dsB.blocks['raw'].data_sources[0]))
        lines.append('UUID Dataset B - '+str(dsB.id))
        lines = "\n".join(lines)
        if (sys.platform == "win32"):
            lines = lines.replace("\n", "\r\n")

        # Create raw data object
        met = mrs_data_raw.DataRawEdit()
        met.data_sources = ['analysis_spectral_plotc_derived_dataset_edit_sum_dif']  # todo bjs - set on/off filenames
        met.sw = ds.sw
        met.frequency = ds.frequency
        met.resppm = ds.resppm
        met.echopeak = ds.echopeak
        met.nucleus = ds.nucleus
        met.seqte = ds.seqte
        met.seqtr = ds.seqtr
        met.voxel_dimensions = ds.blocks['raw'].voxel_dimensions
        met.headers = [lines]
        met.data = dat
        met.transform = ds.blocks['raw'].transform

        # Convert Raw to Dataset - from util_file_import.get_datasets() line 356-ish
        block_classes = {}
        block_classes["raw"] = getattr(importlib.import_module('vespa.analysis.block_raw_edit'), 'BlockRawEdit')
        derived = mrs_dataset.dataset_from_raw(met, block_classes, ds.zero_fill_multiplier)

        # set up associated datasets in original and derived objects
        # - get existing values first from original
        # - overwrite with derived values

        derived.blocks['raw'].data_off_id = ds.id
        derived.blocks['raw'].data_off = ds
        derived.blocks['raw'].data_on_id = dsB.id
        derived.blocks['raw'].data_on = dsB
        derived.blocks['raw'].data_sum_id = ds.blocks['raw'].data_sum_id
        derived.blocks['raw'].data_sum = ds.blocks['raw'].data_sum
        derived.blocks['raw'].data_dif_id = ds.blocks['raw'].data_dif_id
        derived.blocks['raw'].data_dif = ds.blocks['raw'].data_dif
        derived.blocks['raw'].data_sum_indiv_id = ds.blocks['raw'].data_sum_indiv_id
        derived.blocks['raw'].data_sum_indiv = ds.blocks['raw'].data_sum_indiv
        derived.blocks['raw'].data_dif_indiv_id = ds.blocks['raw'].data_dif_indiv_id
        derived.blocks['raw'].data_dif_indiv = ds.blocks['raw'].data_dif_indiv

        if mode == 'plotc_a_plus_b':
            derived.blocks['raw'].data_sum_id = derived.id
            derived.blocks['raw'].data_sum = derived
            for ds in list(top.datasets.values()):
                raw = ds.blocks['raw']
                if hasattr(raw, 'data_sum_id'): raw.data_sum_id = derived.id
                if hasattr(raw, 'data_sum'): raw.data_sum = derived
        elif mode == 'plotc_a_minus_b':
            derived.blocks['raw'].data_dif_id = derived.id
            derived.blocks['raw'].data_dif = derived
            for ds in list(top.datasets.values()):
                raw = ds.blocks['raw']
                if hasattr(raw, 'data_dif_id'): raw.data_dif_id = derived.id
                if hasattr(raw, 'data_dif'): raw.data_dif = derived

        # insert into Tab Notebook
        _tab_dataset._outer_notebook.add_dataset_tab(datasets=[derived, ], force_name=label_new)






#------------------------------------------------------------------------------
#
#  Tab SPECTRAL
#
#------------------------------------------------------------------------------

class TabSpectral(tab_base.Tab, spectral.PanelSpectralUI):
    
    # self-identify tab to notebook, value does not matter, its presence is sufficient.
    IS_SPECTRAL = True

    def __init__(self, tab_dataset, top, block):

        spectral.PanelSpectralUI.__init__(self, tab_dataset.NotebookDataset)
        tab_base.Tab.__init__(self, tab_dataset, top, prefs_module.PrefsSpectral)

        self.top      = top               # application frame
        self.block    = block             # processing object

        self.cursor_span_picks_lines = False

        # the button at bottom of tab can be set to a number of user defined
        # purposes, namely:
        #
        # 1. Automatic Phasing of data
        #
        # 2. Output of current area value to a text file
        #    - each time it is hit a file select dialog comes up
        #    - default file is last filename selected
        #    - if file exists, value is appended, else file created
        #    - both the dataset name, x,y,z voxel numbers, and area values
        #      are saved in a tab delineated text line

        self.user_button_area_fname = ''

        # _svd_scale_intialized performs the same role as
        # tab_base.Tab._scale_intialized (q.v.)
        self._svd_scale_intialized = False

        # Plotting is disabled during some of init. That's because the plot
        # isn't ready to plot, but the population of some controls
        # (e.g. spin controls on the water filter panel) fires their
        # respective change event which triggers a call to plot(). This
        # is a Windows-specific bug.
        # See http://trac.wxwidgets.org/ticket/14583
        # In any case, skipping some calls to plot() will speed things up. =)
        self._plotting_enabled = False
        self.plot_results = None

        # self.plot_C_function refers to either None or a callable object
        # (i.e. a function) that takes two params which are typically the
        # data from plots A & B.
        # The plot_C_map relates the various options to the
        # appropriate function. The functions currently used are so simple
        # that I implement them with anonymous lambda functions
        self.plot_C_map =   {
            util_menu.ViewIdsSpectral.PLOT_C_FUNCTION_NONE      : (lambda a, b: (a+b)*0.0),
            util_menu.ViewIdsSpectral.PLOT_C_FUNCTION_A_MINUS_B : (lambda a, b: a - b),
            util_menu.ViewIdsSpectral.PLOT_C_FUNCTION_B_MINUS_A : (lambda a, b: b - a),
            util_menu.ViewIdsSpectral.PLOT_C_FUNCTION_A_PLUS_B  : (lambda a, b: a + b),
                            }
        # The default plot 3 function is set here
        self.plot_C_function = self.plot_C_map[util_menu.ViewIdsSpectral.PLOT_C_FUNCTION_A_MINUS_B]
        self.plot_C_final = None

        self.initialize_controls()
        self.populate_controls()

        self._plotting_enabled = True

        #------------------------------------------------------------
        # Setup the canvas
        self.process_and_plot()

        #------------------------------------------------------------
        # PubSub subscriptions
        pubsub.subscribe(self.on_check_datasetb_status, "check_datasetb_status")
        pubsub.subscribe(self.on_dataset_keys_change, "dataset_keys_change")

        #------------------------------------------------------------
        # Set window events
        self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, self.on_tab_changed, self.NotebookSpectral)

        # If the sash position isn't recorded in the INI file, we use the
        # arbitrary-ish value of 400.
        if not self._prefs.sash_position_main:
            self._prefs.sash_position_main = 400
        if not self._prefs.sash_position_svd:
            self._prefs.sash_position_svd = 400

        # Under OS X, wx sets the sash position to 10 (why 10?) *after*
        # this method is done. So setting the sash position here does no
        # good. We use wx.CallAfter() to (a) set the sash position and
        # (b) fake an EVT_SPLITTER_SASH_POS_CHANGED.
        wx.CallAfter(self.SplitterWindow.SetSashPosition,
                     self._prefs.sash_position_main, True)
        wx.CallAfter(self.SplitterWindowSvd.SetSashPosition,
                     self._prefs.sash_position_svd, True)
        wx.CallAfter(self.on_splitter)


    @property
    def datasetB(self):
        return self.top.datasets.get(self._tab_dataset.indexAB[1], None)

    @property
    def tabB_spectral(self):
        label = self._tab_dataset.indexAB[1]
        if not label:
            return None
        tab = self._tab_dataset._outer_notebook.get_tab_by_label(label)
        if not tab:
            return None
        tab = tab.get_tab("spectral")
        return tab

    @property
    def do_sync(self):
        sync = self._tab_dataset.sync
        indexB = self._tab_dataset.indexAB[1]
        if indexB and sync:
            return True
        else:
            return False

    @property
    def svd_tab_active(self):
        """Returns True if HLSVD is the active tab on the spectral notebook,
        False otherwise."""
        tab = self.NotebookSpectral.GetPage(self.NotebookSpectral.GetSelection())
        return (tab.Id == self.PanelSvdFilter.Id)


    #=======================================================
    #
    #           GUI Setup Handlers
    #
    #=======================================================

    def initialize_controls(self):
        """
        This methods goes through the widgets and sets up certain sizes
        and constraints for those widgets. This method does not set the
        value of any widget except insofar that it is outside a min/max
        range as those are being set up.

        Use populate_controls() to set the values of the widgets from
        a data object.

        """
        ds = self.dataset
        dim0, dim1, dim2, dim3 = ds.spectral_dims
        sw      = ds.sw
        maxppm  = ds.pts2ppm(0)
        minppm  = ds.pts2ppm(dim0-1)
        ppmlim  = (minppm, maxppm)
        
        # Here I set up the spin controls. Some are standard wxSpinCtrls and
        # some are floatspins (from wx.lib.agw.floatspin or wx_common).
        # For some of our variables, the default increment of +/-1 makes
        # sense, but that's not true for all of them. For example, the phase 1
        # spin control changes in increments of 2.
        # These increments are arbitrary and expressed only once in the code
        # below. If you don't like them, feel free to change them.

        wx_util.configure_spin(self.FloatWidth,       70, 3, 0.5, (0.0, 100.0))
        wx_util.configure_spin(self.FloatFrequency,   70, 3, 0.5, (-1e4, 1e4))
        wx_util.configure_spin(self.FloatAmplitude,   70, 3, None, (0,1e12))
        self.FloatAmplitude.multiplier = 1.1
        wx_util.configure_spin(self.FloatPhase0,      70, 3, 1.0, (-360.0, 360.0))
        wx_util.configure_spin(self.FloatPhase1,      70, 3, 10.0, (-1e5, 1e5))
        wx_util.configure_spin(self.FloatPhase1Pivot, 70, 3, 0.5, (-1000.0, 1000.0))
        wx_util.configure_spin(self.FloatDcOffset,    70, 3, 0.25, (-1e5, 1e5))
        wx_util.configure_spin(self.SpinLeftShift,    70)

        # set up combo selections
        _configure_combo(self.ComboApodization,      constants.Apodization.choices)

        # set up the user function button initial setting
        if self._prefs.user_button_phasing:
            self.ButtonUserFunction.SetLabel('Do Automatic Phasing')
        elif self._prefs.user_button_area:
            self.ButtonUserFunction.SetLabel('Output Area Value')

        # Water Filter widget constraints
        #
        # - note. by setting the default value first, we avoid having an event
        #     sent that sets the block attribute to the min value because the
        #     widget's original value was outside min/max range when they're set

        self.SpinFirLength.SetValue(funct_watfilt.FIR_LENGTH_DEFAULT)
        wx_util.configure_spin(self.SpinFirLength, 60, None, None,
                               (funct_watfilt.FIR_LENGTH_MIN,
                                funct_watfilt.FIR_LENGTH_MAX))

        self.FloatFirWidth.SetValue(funct_watfilt.FIR_HALF_WIDTH_DEFAULT)
        wx_util.configure_spin(self.FloatFirWidth, 60, None,
                                funct_watfilt.FIR_HALF_WIDTH_STEP,
                               (funct_watfilt.FIR_HALF_WIDTH_MIN,
                                funct_watfilt.FIR_HALF_WIDTH_MAX))

        self.FloatFirRipple.SetValue(funct_watfilt.FIR_RIPPLE_DEFAULT)
        wx_util.configure_spin(self.FloatFirRipple, 60, None,
                                funct_watfilt.FIR_RIPPLE_STEP,
                               (funct_watfilt.FIR_RIPPLE_MIN,
                                funct_watfilt.FIR_RIPPLE_MAX))

        self.SpinFirExtrapValue.SetValue(funct_watfilt.FIR_EXTRAPOLATION_POINTS_DEFAULT)
        wx_util.configure_spin(self.SpinFirExtrapValue, 60, None, None,
                               (funct_watfilt.FIR_EXTRAPOLATION_POINTS_MIN,
                                funct_watfilt.FIR_EXTRAPOLATION_POINTS_MAX))

        self.SpinHamLength.SetValue(funct_watfilt.HAM_LENGTH_DEFAULT)
        wx_util.configure_spin(self.SpinHamLength, 60, None, None,
                               (funct_watfilt.HAM_LENGTH_MIN,
                                funct_watfilt.HAM_LENGTH_MAX))

        self.SpinHamExtrapValue.SetValue(funct_watfilt.HAM_EXTRAPOLATION_POINTS_DEFAULT)
        wx_util.configure_spin(self.SpinHamExtrapValue, 60, None, None,
                               (funct_watfilt.HAM_EXTRAPOLATION_POINTS_MIN,
                                funct_watfilt.HAM_EXTRAPOLATION_POINTS_MAX))

        #-------------------------------------------------------------
        # Dataset View setup
        #-------------------------------------------------------------

        self.view = PlotPanelSpectral(self.PanelViewSpectral,
                                      self,
                                      self._tab_dataset,
                                      naxes=3,
                                      reversex=True,
                                      zoom='span',
                                      reference=True,
                                      middle=True,
                                      zoom_button=1,
                                      middle_button=3,
                                      refs_button=2,
                                      do_zoom_select_event=True,
                                      do_zoom_motion_event=True,
                                      do_refs_select_event=True,
                                      do_refs_motion_event=True,
                                      do_middle_select_event=True,
                                      do_middle_motion_event=True,
                                      do_scroll_event=True,
                                      props_zoom=dict(alpha=0.2, facecolor='yellow'),
                                      props_cursor=dict(alpha=0.2, facecolor='gray'),
                                      xscale_bump=0.0,
                                      yscale_bump=0.05,
                                      data = [],
                                      prefs=self._prefs,
                                      dataset=self.dataset,
                                      )

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.view, 1, wx.LEFT | wx.TOP | wx.EXPAND)
        self.PanelViewSpectral.SetSizer(sizer)
        self.view.Fit()
        self.view.change_naxes(1)

        #------------------------------------------------------------
        # SVD tab settings
        #------------------------------------------------------------

        self.view_svd = PlotPanelSvdFilter( self.PanelViewSvd,
                                            self,
                                            self._tab_dataset,
                                            naxes=3,
                                            reversex=True,
                                            zoom='span',
                                            reference=True,
                                            middle=True,
                                            zoom_button=1,
                                            middle_button=3,
                                            refs_button=2,
                                            do_zoom_select_event=True,
                                            do_zoom_motion_event=True,
                                            do_refs_select_event=True,
                                            do_refs_motion_event=True,
                                            do_middle_select_event=True,
                                            do_middle_motion_event=True,
                                            do_scroll_event=True,
                                            props_zoom=dict(alpha=0.2, facecolor='yellow'),
                                            props_cursor=dict(alpha=0.2, facecolor='gray'),
                                            xscale_bump=0.0,
                                            yscale_bump=0.05,
                                            data = [],
                                            prefs=self._prefs,
                                            dataset=self.dataset,
                                            )

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.view_svd, 1, wx.LEFT | wx.TOP | wx.EXPAND)
        self.PanelViewSvd.SetSizer(sizer)
        self.view_svd.Fit()


        #------------------------------------------------------------
        # Results Control

        self.list_svd_results = CheckListCtrl(self.PanelResults, self)
        sizer = wx.BoxSizer()
        sizer.Add(self.list_svd_results, 1, wx.EXPAND)
        self.PanelResults.SetSizer(sizer)
        self.Fit()

        self.list_svd_results.InsertColumn(0, "Rank", wx.LIST_FORMAT_CENTER)
        self.list_svd_results.InsertColumn(1, "PPM", wx.LIST_FORMAT_CENTER)
        self.list_svd_results.InsertColumn(2, "Freq", wx.LIST_FORMAT_CENTER)
        self.list_svd_results.InsertColumn(3, "Damping", wx.LIST_FORMAT_CENTER)
        self.list_svd_results.InsertColumn(4, "Phase", wx.LIST_FORMAT_CENTER)
        self.list_svd_results.InsertColumn(5, "Area", wx.LIST_FORMAT_CENTER)

        self.list_svd_results.SetColumnWidth(0, 60)
        self.list_svd_results.SetColumnWidth(1, 60)
        self.list_svd_results.SetColumnWidth(2, 80)
        self.list_svd_results.SetColumnWidth(3, 80)
        self.list_svd_results.SetColumnWidth(4, 80)
        self.list_svd_results.SetColumnWidth(5, 80)

        self.list_svd_results.EnableCheckBoxes(enable=True)

        # Manual selection is the default
        self.RadioSvdManual.SetValue(True)


        #------------------------------------------------------------
        # Algorithm Controls

        # Here's our safe hard limit for # of data points.
        min_data_points = HLSVD_MIN_DATA_POINTS + 1

        # min_data_points is now in safe territory. Since everything in this
        # code is expressed in powers of 2, we'll keep the same style with
        # min_data_points. Plus it gives us a little margin of error.
        i = 1
        while 2**i < min_data_points:
            i += 1
        min_data_points = 2**i

        # Max # of data points can't be more than the # of points in the data.
        max_data_points = ds.spectral_dims[0]
        self.SliderDataPoints.SetRange(min_data_points, max_data_points)

        # Range is set; now set value.
        n_data_points = min(funct_watfilt.SVD_N_DATA_POINTS, max_data_points)
        self.SliderDataPoints.SetValue(n_data_points)

        # Singular values
        self.SliderSingularValues.SetRange(1, int(HLSVD_MAX_SINGULAR_VALUES))
        self.SliderSingularValues.SetValue(funct_watfilt.SVD_N_SINGULAR_VALUES)

        # There's a lot of ways to change sliders -- dragging the thumb,
        # clicking on either side of the thumb, using the arrow keys to
        # move one tick at a time, and hitting home/end. Fortunately all
        # platforms cook these down into a simple "the value changed" event.
        # Unfortunately it has different names depending on the platform.
        if "__WXMAC__" in wx.PlatformInfo:
            event = wx.EVT_SCROLL_THUMBRELEASE
        else:
            event = wx.EVT_SCROLL_CHANGED

        for slider in (self.SliderDataPoints, self.SliderSingularValues):
                self.Bind(event, self.on_slider_changed, slider)

        wx_util.configure_spin(self.FloatSvdThreshold, 70, 2, 0.1,
                               (constants.SvdThreshold.MIN, constants.SvdThreshold.MAX))

        _configure_combo(self.ComboSvdThresholdUnit, constants.SvdThresholdUnit.choices)

        wx_util.configure_spin(self.FloatSvdExcludeLipidStart, 50, 2, 0.5, ppmlim)
        wx_util.configure_spin(self.FloatSvdExcludeLipidEnd,   50, 2, 0.5, ppmlim)


    def populate_controls(self, preset=False):
        """
        Populates the widgets with relevant values from the data object.
        It's meant to be called when a new data object is loaded.

        This function trusts that the data object it is given doesn't violate
        any rules. Whatever is in the data object gets slapped into the
        controls, no questions asked.

        This function is, however, smart enough to enable and disable
        other widgets depending on settings.

        """
        ds = self.dataset
        dim0, dim1, dim2, dim3 = ds.spectral_dims
        voxel   = self._tab_dataset.voxel
        maxppm  = ds.pts2ppm(0)
        minppm  = ds.pts2ppm(dim0-1)

        #------------------------------------------------------------
        # ECC Filters

        # (Re)populate the ECC list
        self.ComboEcc.Clear()
        for item in funct_ecc.ECC_MENU_ITEMS:
            self.ComboEcc.Append(item, "")
        if self.block.set.ecc_method == '':
            self.block.set.ecc_method = 'None'
        self.ComboEcc.SetStringSelection(self.block.set.ecc_method)
        if self.block.set.ecc_method == 'None':
            self.PanelEccBrowse.Hide()
        else:
            self.PanelEccBrowse.Show()
            self.TextEccFilename.SetValue(self.block.set.ecc_filename)

        #------------------------------------------------------------
        # Water Filter settings

        # (Re)populate the Water Filter list

        self.ComboWater.Clear()
        for item in funct_watfilt.WATFILT_MENU_ITEMS:
            self.ComboWater.Append(item, "")
        if self.block.set.water_filter_method == '':
            self.block.set.water_filter_method = 'None'
        self.ComboWater.SetStringSelection(self.block.set.water_filter_method)

        # set all the other filter panel widgets to attribute settings
        self.SpinFirLength.SetValue(self.block.set.fir_length)
        self.FloatFirWidth.SetValue(self.block.set.fir_half_width)
        self.FloatFirRipple.SetValue(self.block.set.fir_ripple)
        self.SpinFirExtrapValue.SetValue(self.block.set.fir_extrapolation_point_count)
        self.ComboFirExtrapMethod.Clear()
        for label in funct_watfilt.FIR_EXTRAPOLATION_ALL:
            self.ComboFirExtrapMethod.Append(label, "")
            if self.block.set.fir_extrapolation_method == label:
                # This is the active filter, so I select it in the list
                self.ComboFirExtrapMethod.SetSelection(self.ComboFirExtrapMethod.GetCount() - 1)
        self.SpinFirExtrapValue.Enable(self.block.set.fir_extrapolation_method == 'AR Model')

        self.SpinHamLength.SetValue(self.block.set.ham_length)
        self.SpinHamExtrapValue.SetValue(self.block.set.ham_extrapolation_point_count)
        self.ComboHamExtrapMethod.Clear()
        for label in funct_watfilt.HAM_EXTRAPOLATION_ALL:
            self.ComboHamExtrapMethod.Append(label, "")
            if self.block.set.ham_extrapolation_method == label:
                # This is the active filter, so I select it in the list
                self.ComboHamExtrapMethod.SetSelection(self.ComboHamExtrapMethod.GetCount() - 1)
        self.SpinHamExtrapValue.Enable(self.block.set.ham_extrapolation_method == 'AR Model')

        choice = self.ComboWater.GetStringSelection()
        if choice == 'None' or choice == 'SVD - water filter':
            self.PanelWaterFir.Hide()
            self.PanelWaterHamming.Hide()
        elif choice == 'FIR - water filter':
            self.PanelWaterFir.Show()
            self.PanelWaterHamming.Hide()
        elif choice == 'Hamming - water filter':
            self.PanelWaterFir.Hide()
            self.PanelWaterHamming.Show()

        #------------------------------------------------------------
        # Spectral Control setup

        self.ComboDataB.SetStringSelection(str(self._tab_dataset.indexAB[1]))
        self.CheckSync.Disable()

        self.CheckFlip.SetValue(self.block.set.flip)
        self.CheckFFT.SetValue(self.block.set.fft)
        self.CheckChop.SetValue(self.block.set.chop)

        # Apodization width is disabled if there's no method chosen
        apodize = constants.Apodization.choices[self.block.set.apodization]
        self.ComboApodization.SetStringSelection(apodize)
        self.FloatWidth.SetValue(self.block.set.apodization_width)
        self.FloatWidth.Enable(bool(self.block.set.apodization))
        self.ComboZeroFill.SetStringSelection(str(int(self.block.set.zero_fill_multiplier)))
        self.FloatFrequency.SetValue(ds.get_frequency_shift(voxel))
        self.CheckFreqLock.SetValue(self.block.frequency_shift_lock)
        self.FloatAmplitude.SetValue(self.block.set.amplitude)
        self.FloatPhase0.SetValue(ds.get_phase_0(voxel))
        self.CheckPhaseLock.SetValue(self.block.phase_lock)
        self.FloatPhase1.SetValue(ds.get_phase_1(voxel))
        self.CheckZeroPhase1.SetValue(self.block.phase_1_lock_at_zero)
        self.FloatPhase1Pivot.SetValue(self.block.set.phase_1_pivot)
        self.FloatDcOffset.SetValue(self.block.set.dc_offset)
        self.SpinLeftShift.SetValue(self.block.set.left_shift_value)
        self.SpinLeftShift.SetRange(0, ds.raw_dims[0])
        self.CheckCorrectPhase1.SetValue(self.block.set.left_shift_correct)

        #------------------------------------------------------------
        # SVD tab settings
        #------------------------------------------------------------

        # SVD Algorithm Controls
        self.SliderDataPoints.SetValue(self.block.get_data_point_count(voxel))
        self.SliderSingularValues.SetValue(self.block.get_signal_singular_value_count(voxel))

        # SVD Results Controls
        self.svd_checklist_update()

        self.RadioSvdApplyThreshold.SetValue(self.block.set.svd_apply_threshold)
        self.FloatSvdThreshold.SetValue(self.block.set.svd_threshold)
        if self.block.set.svd_apply_threshold:
            self.FloatSvdThreshold.Enable()
        else:
            self.FloatSvdThreshold.Disable()

        item = constants.SvdThresholdUnit.choices[self.block.set.svd_threshold_unit]
        self.ComboSvdThresholdUnit.SetStringSelection(item)
        if item == 'PPM':
            # threshold value used to be in Hz and allowed to range +/- 200
            #  units can now be PPM, so check if we are outside min/max ppm
            val = self.block.set.svd_threshold
            if val < minppm:
                self.block.set.svd_threshold = minppm
            elif val > maxppm:
                self.block.set.svd_threshold = maxppm
                self.FloatSvdThreshold.SetValue(self.block.set.svd_threshold)
        
        self.CheckSvdExcludeLipid.SetValue(self.block.set.svd_exclude_lipid)
        self.FloatSvdExcludeLipidStart.SetValue(self.block.set.svd_exclude_lipid_start)    
        self.FloatSvdExcludeLipidEnd.SetValue(self.block.set.svd_exclude_lipid_end)
        if self.block.set.svd_exclude_lipid:
            self.FloatSvdExcludeLipidStart.Enable()
            self.FloatSvdExcludeLipidEnd.Enable()
        else:
            self.FloatSvdExcludeLipidStart.Disable()
            self.FloatSvdExcludeLipidEnd.Disable()
            
            



    #=======================================================
    #
    #           Global and Menu Event Handlers
    #
    #=======================================================

    def on_destroy(self, event):
        tab_base.Tab.on_destroy(self, event)
        pubsub.subscribe(self.on_check_datasetb_status, "check_datasetb_status")
        pubsub.subscribe(self.on_dataset_keys_change, "dataset_keys_change")

    def on_activation(self):
        tab_base.Tab.on_activation(self)

        # these BlockSpectral object values may be changed by other tabs, so
        # update their widget values on activation of this tab
        voxel      = self._tab_dataset.voxel
        freq_shift = self.dataset.get_frequency_shift(voxel)
        phase_0    = self.dataset.get_phase_0(voxel)
        phase_1    = self.dataset.get_phase_1(voxel)
        self.FloatFrequency.SetValue(freq_shift)
        self.FloatPhase0.SetValue(phase_0)
        self.FloatPhase1.SetValue(phase_1)

        # Force re-process when switching from a PreProcess tab in case it has
        # changed. It will also trigger other times, but ...
        self.process_and_plot()


    def on_menu_view_option(self, event):
        event_id = event.GetId()

        if self._prefs.handle_event(event_id):

            if event_id in (util_menu.ViewIdsSpectral.ZERO_LINE_SHOW,
                            util_menu.ViewIdsSpectral.ZERO_LINE_TOP,
                            util_menu.ViewIdsSpectral.ZERO_LINE_MIDDLE,
                            util_menu.ViewIdsSpectral.ZERO_LINE_BOTTOM,
                            util_menu.ViewIdsSpectral.XAXIS_SHOW,
                           ):
                self.view.update_axes()
                self.view_svd.update_axes()
                self.view.canvas.draw()
                self.view_svd.canvas.draw()


            elif event_id in (util_menu.ViewIdsSpectral.DATA_TYPE_REAL,
                              util_menu.ViewIdsSpectral.DATA_TYPE_IMAGINARY,
                              util_menu.ViewIdsSpectral.DATA_TYPE_MAGNITUDE,
                              util_menu.ViewIdsSpectral.DATA_TYPE_SUMMED,
                              util_menu.ViewIdsSpectral.XAXIS_PPM,
                              util_menu.ViewIdsSpectral.XAXIS_HERTZ,
                           ):
                if event_id == util_menu.ViewIdsSpectral.DATA_TYPE_REAL:
                    self.view.set_data_type_real()
                    self.view_svd.set_data_type_real()
                elif event_id == util_menu.ViewIdsSpectral.DATA_TYPE_IMAGINARY:
                    self.view.set_data_type_imaginary()
                    self.view_svd.set_data_type_imaginary()
                elif event_id == util_menu.ViewIdsSpectral.DATA_TYPE_MAGNITUDE:
                    self.view.set_data_type_magnitude()
                    self.view_svd.set_data_type_magnitude()
                elif event_id == util_menu.ViewIdsSpectral.DATA_TYPE_SUMMED:
                    # no effect on main view plots
                    self.view_svd.set_data_type_summed(index=[1])


                self.view.update(no_draw=True)
                self.view.set_phase_0(0.0, no_draw=True)
                self.view_svd.update(no_draw=True)
                self.view_svd.set_phase_0(0.0, no_draw=True)
                self.set_plot_c()
                self.view.canvas.draw()
                self.view_svd.canvas.draw()

            elif event_id in (util_menu.ViewIdsSpectral.AREA_CALC_PLOT_A,
                              util_menu.ViewIdsSpectral.AREA_CALC_PLOT_B,
                              util_menu.ViewIdsSpectral.AREA_CALC_PLOT_C,
                             ):
                area, rms = self.view.calculate_area()
                if self._prefs.area_calc_plot_a:
                    index = 0
                    labl = 'A'
                elif  self._prefs.area_calc_plot_b:
                    index = 1
                    labl = 'B'
                elif  self._prefs.area_calc_plot_c:
                    index = 2
                    labl = 'C'
                self.top.statusbar.SetStatusText(self.build_area_text(area[index], rms[index], plot_label=labl), 3)

            elif event_id in (util_menu.ViewIdsSpectral.PLOT_C_FUNCTION_NONE,
                              util_menu.ViewIdsSpectral.PLOT_C_FUNCTION_A_MINUS_B,
                              util_menu.ViewIdsSpectral.PLOT_C_FUNCTION_B_MINUS_A,
                              util_menu.ViewIdsSpectral.PLOT_C_FUNCTION_A_PLUS_B,
                             ):
                if event_id == util_menu.ViewIdsSpectral.PLOT_C_FUNCTION_NONE:
                    if len(self.view.axes) == 3:
                        if self.datasetB:
                            self.view.change_naxes(2)
                        else:
                            self.view.change_naxes(1)
                else:
                    if len(self.view.axes) != 3:
                        self.view.change_naxes(3)

                self.plot_C_function = self.plot_C_map[event_id]

                self.set_plot_c()
                self.view.canvas.draw()

            elif event_id in (util_menu.ViewIdsSpectral.USER_BUTTON_PHASING,
                              util_menu.ViewIdsSpectral.USER_BUTTON_AREA,
                             ):
                if event_id == util_menu.ViewIdsSpectral.USER_BUTTON_PHASING:
                    label = 'Do Automatic Phasing'
                elif event_id == util_menu.ViewIdsSpectral.USER_BUTTON_AREA:
                    label = 'Output Area Value'
                self.ButtonUserFunction.SetLabel(label)


    def on_menu_view_derived(self, event):
        event_id = event.GetId()

        if event_id == util_menu.ViewIdsSpectral.PLOTA_FID_TO_BINARY:
            filename = common_dialogs.save_as("Plot A FID to Binary", "BIN files (*.bin)|*.bin", default_filename='analysis_plota_fid_to_binary.bin')
            if filename:
                data = self.view.get_data(0)
                data = np.fft.ifft(np.fft.fftshift(data)) * len(data)
                util_fileio.dump_iterable(filename, data)
        elif event_id == util_menu.ViewIdsSpectral.PLOTA_FID_TO_ASCII:
            filename = common_dialogs.save_as("Plot A FID to ASCII", "TXT files (*.txt)|*.txt", default_filename='analysis_plota_fid_to_text.txt')
            if filename:
                data = self.view.get_data(0)
                data = np.fft.ifft(np.fft.fftshift(data)) * len(data)
                data.tofile(filename, sep=' ')
        elif event_id == util_menu.ViewIdsSpectral.PLOTA_SPECTRUM_TO_BINARY:
            filename = common_dialogs.save_as("Plot A Spectrum to Binary", "BIN files (*.bin)|*.bin", default_filename='analysis_plota_spectrum_to_binary.bin')
            if filename:
                data = self.view.get_data(0)
                util_fileio.dump_iterable(filename, data)
        elif event_id == util_menu.ViewIdsSpectral.PLOTA_SPECTRUM_TO_ASCII:
            filename = common_dialogs.save_as("Plot A Spectrum to ASCII", "TXT files (*.txt)|*.txt", default_filename='analysis_plota_spectrum_to_text.txt')
            if filename:
                data = self.view.get_data(0)
                data.tofile(filename, sep=' ')
        elif event_id == util_menu.ViewIdsSpectral.PLOTA_FID_TO_VIFF:
            self.dump_to_viff('plota_fid')

        if event_id == util_menu.ViewIdsSpectral.PLOTC_FID_TO_BINARY:
            filename = common_dialogs.save_as("Plot C FID to Binary", "BIN files (*.bin)|*.bin", default_filename='analysis_plotc_fid_to_binary.bin')
            if filename:
                data = self.view.get_data(2)
                data = np.fft.ifft(np.fft.fftshift(data)) * len(data)
                util_fileio.dump_iterable(filename, data)
        elif event_id == util_menu.ViewIdsSpectral.PLOTC_FID_TO_ASCII:
            filename = common_dialogs.save_as("Plot C FID to ASCII", "TXT files (*.txt)|*.txt", default_filename='analysis_plotc_fid_to_text.txt')
            if filename:
                data = self.view.get_data(2)
                data = np.fft.ifft(np.fft.fftshift(data)) * len(data)
                data.tofile(filename, sep=' ')
        elif event_id == util_menu.ViewIdsSpectral.PLOTC_SPECTRUM_TO_BINARY:
            filename = common_dialogs.save_as("Plot C Spectrum to Binary", "BIN files (*.bin)|*.bin", default_filename='analysis_plotc_spectrum_to_binary.bin')
            if filename:
                data = self.view.get_data(2)
                util_fileio.dump_iterable(filename, data)
        elif event_id == util_menu.ViewIdsSpectral.PLOTC_SPECTRUM_TO_ASCII:
            filename = common_dialogs.save_as("Plot C Spectrum to ASCII", "TXT files (*.txt)|*.txt", default_filename='analysis_plotc_spectrum_to_text.txt')
            if filename:
                data = self.view.get_data(2)
                data.tofile(filename, sep=' ')
        elif event_id == util_menu.ViewIdsSpectral.PLOTC_FID_TO_VIFF:
            self.dump_to_viff('plotc_fid')

        elif event_id == util_menu.ViewIdsSpectral.SVDA_SPECTRUM_TO_VIFF:
            self.dump_to_viff('svda_spectrum')
        elif event_id == util_menu.ViewIdsSpectral.SVDB_SPECTRUM_CHECKED_SUM_TO_VIFF:
            self.dump_to_viff('svdb_spectrum_checked_sum')
        elif event_id == util_menu.ViewIdsSpectral.SVDB_FIDS_CHECKED_SUM_TO_VIFF:
            self.dump_to_viff('svdb_fids_checked_sum')
        elif event_id == util_menu.ViewIdsSpectral.SVDC_SPECTRUM_TO_VIFF:
            self.dump_to_viff('svdc_spectrum')
        elif event_id == util_menu.ViewIdsSpectral.PLOTC_A_PLUS_B_TO_NEW_TAB:
            if self.datasetB is None:
                msg = 'Warning (Derived Dataset PlotC): Please select Dataset for Plot B'
                r = common_dialogs.message(msg, style=common_dialogs.I_OK)
                return
            derived_dataset(self.dataset, self.datasetB, self.view, self._tab_dataset, self.top, 'plotc_a_plus_b')
        elif event_id == util_menu.ViewIdsSpectral.PLOTC_A_MINUS_B_TO_NEW_TAB:
            if self.datasetB is None:
                msg = 'Warning (Derived Dataset PlotC): Please select Dataset for Plot B'
                r = common_dialogs.message(msg, style=common_dialogs.I_OK)
                return
            derived_dataset(self.dataset, self.datasetB, self.view, self._tab_dataset, self.top, 'plotc_a_minus_b')



    def on_menu_view_output(self, event):
        event_id = event.GetId()
        
        formats = { util_menu.ViewIdsSpectral.VIEW_TO_PNG : "PNG",
                    util_menu.ViewIdsSpectral.VIEW_TO_SVG : "SVG",
                    util_menu.ViewIdsSpectral.VIEW_TO_EPS : "EPS",
                    util_menu.ViewIdsSpectral.VIEW_TO_PDF : "PDF",
                  }

        if event_id in formats:
            format = formats[event_id]
            lformat = format.lower()
            filter_ = "%s files (*.%s)|*.%s" % (format, lformat, lformat)

            if self.svd_tab_active:
                figure = self.view_svd.figure
            else:
                figure = self.view.figure

            filename = common_dialogs.save_as("", filter_)

            if filename:
                msg = ""
                try:
                    figure.savefig( filename,
                                    dpi=300,
                                    facecolor='w',
                                    edgecolor='w',
                                    orientation='portrait',
                                    #papertype='letter',
                                    format=None,
                                    transparent=False)
                except IOError:
                    msg = """I can't write the file "%s".""" % filename

                if msg:
                    common_dialogs.message(msg, style=common_dialogs.E_OK)

        elif event_id == util_menu.ViewIdsSpectral.SVD_TABLE_VALUES_TO_VIFF:
            self.export_hlsvd_table_to_xml()
        elif event_id == util_menu.ViewIdsSpectral.SVD_TABLE_VALUES_TO_CSV:
            self.export_hlsvd_table_to_csv()



    def on_dataset_keys_change(self, keys):
        """ Pubsub notifications handler. """
        # The list of available datasets has changed so we update the
        # list of datasets available for plot B. We make note of the
        # currently selected item so we can re-select it if it's
        # still in the list.
        _keys = [item for item in keys]
        do_update = False
        if self._tab_dataset.indexAB[0] in _keys:
            # need this if because for some reason a deleted tab is still being
            # reported to by PubSub but 'DatasetX' isn't in the _keys list.
            #
            # next lines removes current dataset from plot B list
            _keys.remove(self._tab_dataset.indexAB[0])
            labels = ["None"] + _keys
            # save current value to see if we can re-select it
            valueB = self.ComboDataB.GetStringSelection()
            self.ComboDataB.SetItems(labels)
            indexB = self.ComboDataB.FindString(valueB)
            if (indexB == wx.NOT_FOUND) or (valueB == "None"):
                # Ooops, the previously selected item is gone.
                indexB = 0
                if valueB != "None":
                    do_update = True
                    self._tab_dataset.sync = False
                    self.CheckSync.SetValue(False)
                    self.CheckSync.Disable()
                    self.view.change_naxes(1)
                valueB = None
    
            self.ComboDataB.SetSelection(indexB)
            self._tab_dataset.indexAB[1] = valueB
    
            if do_update:
                self.process_and_plot()


    def on_check_datasetb_status(self, label):
        """ Pubsub notifications handler. """
        # This gets called after any change to a dataset that requires
        # it to be re-processed. Any other dataset using the referenced
        # dataset will have to update its PlotB
        # 
        # first if checks if the main Dataset tab is still in existence. Had
        # bug where this pubsub was called for a tab being closed.
        if self._tab_dataset.indexAB[0] in list(self.top.datasets.keys()):
            if label == self._tab_dataset.indexAB[1]:
                self.plot()


    def dump_to_viff(self, key):
        """
        This method is used to 'dump' spectra from the View canvas, in it's
        visible form (ie. apod, b0 zf, etc. applied) into a VIFF Raw Data file.
        Typically, this is done by taking the spectrum and performing a naive
        ifft() and fftshift() to the spectrum to get FID data to put into the
        VIFF Raw Data file.  In some cases, the data is already in time domain
        FID form, so nothing more is done to it.

        A limited amount of header data and comment are created to document the
        data output. See below for the 'key' values that control what data and
        comments are output.

        Key Values
        -----------------
        'svd_data'
            - data from Plot A in Spectral-SVD tab, this is the processed form
               for the raw MRS data that the HLSVD algorith acts on.
            - has b0 shift, apod, zf, fft, flip, dc, ampl, NO: ph0/ph1
        'svd_peaks_checked_sum'
            - data from Plot B, summed spectrum of HLSVD results that are
               checked in the results table
            - these are processed similar to the data in Plot A for equivalence
            - has b0 shift, apod, zf, fft, flip, dc, ampl, NO: ph0/ph1
        'svd_peaks_checked_diff'
            - data from Plot C, difference between Plot A - Plot B
            - inherently processed similar to the data in Plot A for equivalence
            - has b0 shift, apod, zf, fft, flip, dc, ampl, NO: ph0/ph1
        'svd_fids_checked_sum'
            - Calculated HLSVD FIDs
            - b0 shift is applied to these FIDs before they are output
            - NO apod, zf, fft, flip, ph0/ph1, dc, ampl applied
        'plotc_fid'
            - Derived spectrum, but without FFT or ph0/1 or chop, or flip

        """

        msg   = ''
        voxel = self._tab_dataset.voxel
        ds    = self.dataset
        if ds:
            block    = ds.blocks["spectral"]
            results  = block.chain.run([voxel], entry='viff')
            defname  = "analysis_export_hlsvd_data.xml"
            label    = "Save HLSVDPro Data to VIFF XML Raw Filename"
            filename = common_dialogs.save_as(label, "XML files (*.xml)|*.xml", default_filename=defname)
            if filename:

                stamp = util_time.now(util_time.ISO_TIMESTAMP_FORMAT).split('T')
            
                if key == 'svda_spectrum':
                    lines = ['Export from Analysis Spectral SVD Tab - PlotA Spectrum Data']
                    lines.append('------------------------------------------------------------------')
                    lines.append('An IFFT() was performed on the Plot A spectrum and the resultant FID written into')
                    lines.append('this VIFF XML Raw Data file.')
                elif key == 'svdb_spectrum_checked_sum':
                    lines = ['Export from Analysis SVD Spectral Tab - PlotB Checked SVD Summed Spectrum']
                    lines.append('------------------------------------------------------------------')
                    lines.append('An IFFT() was performed on the summed SVD peak spectrum and the resultant FID written')
                    lines.append('into this VIFF XML Raw Data file.')
                elif key == 'svdb_fids_checked_sum':
                    lines = ['Export from Analysis Spectral SVD Tab - PlotB Checked SVD Summed FID Data']
                    lines.append('------------------------------------------------------------------')
                    lines.append('The time data for the summed SVD FIDs was written into this VIFF XML Raw Data file.')
                    lines.append('You can process it the same way your would any other time domain FID data.')
                elif key == 'svdc_spectrum':
                    lines = ['Export from Analysis Spectral SVD Tab - PlotC (PlotA - PlotB) Spectrum Data']
                    lines.append('------------------------------------------------------------------')
                    lines.append('An IFFT() was performed on the Plot C spectrum and the resultant FID written into')
                    lines.append('this VIFF XML Raw Data file.')

                elif key == 'plota_fid':
                    lines = ['Export from Analysis Spectral General Tab - PlotA FID Data']
                    lines.append('------------------------------------------------------------------')
                    lines.append('An IFFT() was performed on the Plot A spectrum and the resultant FID written into')
                    lines.append('this VIFF XML Raw Data file.')

                elif key == 'plotc_fid':
                    if self.plot_C_function == self.plot_C_map[util_menu.ViewIdsSpectral.PLOT_C_FUNCTION_NONE]:
                        labl_func = 'None'
                    elif self.plot_C_function == self.plot_C_map[util_menu.ViewIdsSpectral.PLOT_C_FUNCTION_A_MINUS_B]:
                        labl_func = 'PlotA - PlotB'
                    elif self.plot_C_function == self.plot_C_map[util_menu.ViewIdsSpectral.PLOT_C_FUNCTION_B_MINUS_A]:
                        labl_func = 'PlotB - PlotA'
                    elif self.plot_C_function == self.plot_C_map[util_menu.ViewIdsSpectral.PLOT_C_FUNCTION_A_PLUS_B]:
                        labl_func = 'PlotA + PlotB'

                    lines = ['Export from Analysis Spectral General Tab - PlotC ('+labl_func+') FID Data']
                    lines.append('------------------------------------------------------------------')
                    lines.append('An IFFT() was performed on the Plot C spectrum and the resultant FID written into')
                    lines.append('this VIFF XML Raw Data file.')

                lines.append(' ')
                lines.append('Creation_date            - '+stamp[0])
                lines.append('Creation_time            - '+stamp[1])
                lines.append(' ')
                lines = "\n".join(lines)
                if (sys.platform == "win32"):
                    lines = lines.replace("\n", "\r\n")

                if key == 'svda_spectrum':
                    dat = results['svd_data'].copy()
                    dat = np.fft.ifft(np.fft.fftshift(dat)) * len(dat)
                elif key == 'svdb_spectrum_checked_sum':
                    dat = results['svd_peaks_checked_sum'].copy()
                    dat = np.fft.ifft(np.fft.fftshift(dat)) * len(dat)
                elif key == 'svdb_fids_checked_sum':
                    dat = results[key].copy()
                elif key =='svd_spectrum_checked_diff':
                    dat = results['svd_data'].copy() - results['svd_peaks_checked_sum'].copy()
                    dat = np.fft.ifft(np.fft.fftshift(dat)) * len(dat)
                elif key == 'plota_fid':
                    dat = self.view.get_data(0)
                    dat = np.fft.ifft(np.fft.fftshift(dat)) * len(dat)
                elif key == 'plotc_fid':
                    dat = self.view.get_data(2)
                    dat = np.fft.ifft(np.fft.fftshift(dat)) * len(dat)

                met = mrs_data_raw.DataRaw()
                met.data_sources = ['analysis_spectral_general_tab_plotc_data_export']
                met.headers      = [lines]
                met.sw           = ds.sw
                met.frequency    = ds.frequency
                met.resppm       = ds.resppm
                met.data         = dat
                try:
                    util_export.export(filename, [met], None, lines, False)
                except IOError:
                    msg = """I can't write the file "%s".""" % filename
            
                if msg:
                    common_dialogs.message(msg, style=common_dialogs.E_OK)
                    return


    def export_hlsvd_table_to_xml(self):
        """
        Converts the HLSVD results values seen in the Spectral SVD tab into an
        Vespa VIFF XML file that can be read into Datasim's Macromolecule tab
        to recreate those lines/peaks.  The comment in the XML file has a text
        print out of the table. The XML nodes contain separate numpy encoded
        arrays of ppm, area, phase, damp, and flag values.

        Note. The 'Freq' column values have been converted into PPMs, and the
        'Damping' column values have also been converted into PPM widths, both
        to preserve field strength independence of the model.

        In this method, we also clip the damping values to be a min of -0.001
        (as opposed to the occasional positive values HLSVD serves up), since
        these appear as very narrow lines, and Datasim did not know what to do
        with a positive growing exponential lineshape value.

        Note. Datasim MMol peak widths are described as (+) values, while HLSVD
        dampings are typically (-). We convert HLSVD damping values to positive
        numbers with units of [sec] (by dividing them by 1000.0) and then
        convert that number to a relative PPM width value.

        """

        defname = "analysis_export_hlsvd_table.xml"
        label = "Save HLSVDPro Table Values to XML Filename"
        fname = common_dialogs.save_as(label, "XML files (*.xml)|*.xml", default_filename=defname)
        if fname:

            voxel   = self._tab_dataset.voxel
            ndp     = self.block.get_data_point_count(voxel)
            nssv    = self.block.get_signal_singular_value_count(voxel)
            svd_out = self.block.get_svd_output(voxel)

            freqs = svd_out.frequencies
            ampls = svd_out.amplitudes
            damps = svd_out.damping_factors
            phase = svd_out.phases
            flags = svd_out.in_model

            # convert freqs->ppms and damps->width_ppm for field independence
            ppms = self.dataset.resppm - (freqs*1000.0/self.dataset.frequency)
            damps_ppm = (1.0/(np.pi*(-damps/1000.0))) / self.dataset.frequency
            # clip damps_ppm to not have (-) numbers
            damps_ppm = np.clip(damps_ppm, 0.001, None)   # min based on Datasim min PPM linewidth in ppms

            comment = []
            comment.append('Exported HLSVDPro water filter results from Vespa-Analysis Spectral Tab')
            comment.append(' Source:  '+self.dataset.blocks["raw"].data_source)
            comment.append(' Voxel (x,y,z):  %d, %d, %d' % (voxel[0],voxel[1],voxel[2]))
            comment.append(' Data points fitted:  %d,   Number of singular values: %d' % (ndp, nssv))
            comment.append(' ')
            comment.append('HLSVDPro Table values - just FYI ')
            comment.append(' Freq [ppm]\t Freq [kHz]\t Area\t Ph0 [deg]\t Damp [ms]\t Damp [ppm]\t Flag')
            for ppm, fre, amp, pha, damp, damp_ppm, flag in zip(ppms, freqs, ampls, phase, damps, damps_ppm, flags):
                txt =  '{0:.2f}\t {1:.7f}\t {2:.7f}\t {3:.3f}\t {4:.7f}\t {5:.5f}\t {6}\t'.format(ppm, fre, amp, pha, damp, damp_ppm, flag)
                comment.append(txt)
            comment = "\n".join(comment)


            LOCAL_ENCODING = "npy zlib base64"
            e = Element("analysis_export_hlsvd", {"id": 'none', "version": '1.0.0'})
            e.append(util_xml.numpy_array_to_element(ppms,  'ppms',   encoding=LOCAL_ENCODING))
            e.append(util_xml.numpy_array_to_element(ampls, 'ampls',  encoding=LOCAL_ENCODING))
            e.append(util_xml.numpy_array_to_element(phase, 'phases', encoding=LOCAL_ENCODING))
            e.append(util_xml.numpy_array_to_element(flags, 'flags',  encoding=LOCAL_ENCODING))
            e.append(util_xml.numpy_array_to_element(damps_ppm, 'damps_ppm',  encoding=LOCAL_ENCODING))

            STANDARD_COMMENT = """This XML file is in Vespa Interchange File Format (VIFF). \nIt was created with Vespa version %s. """ % util_misc.get_vespa_version()

            root = ElementTree.Element(common_constants.Export.ROOT_ELEMENT_NAME,
                                       {"version": common_constants.Export.VERSION})

            # We append an XML comment that we hope is informative. This is the same
            # for every exported file. Not to be confused with user comment.
            root.append(ElementTree.Comment(STANDARD_COMMENT))
            util_xml.TextSubElement(root, "timestamp", util_time.now().isoformat())
            util_xml.TextSubElement(root, "comment", comment)
            root.append(e)
            util_xml.indent(root)   # Prettify the XML
            tree = ElementTree.ElementTree(root)
            tree.write(fname, "utf-8")


    def export_hlsvd_table_to_csv(self):

        defname = "analysis_export_hlsvd_table.csv"
        label = "Save HLSVDPro Table Values to CSV Filename"
        fname = common_dialogs.save_as(label, "CSV files (*.csv)|*.csv", default_filename=defname)
        if fname:

            voxel = self._tab_dataset.voxel
            ndp = self.block.get_data_point_count(voxel)
            nssv = self.block.get_signal_singular_value_count(voxel)
            svd_out = self.block.get_svd_output(voxel)

            freqs = svd_out.frequencies
            ampls = svd_out.amplitudes
            damps = svd_out.damping_factors
            phase = svd_out.phases
            flags = svd_out.in_model

            # convert freqs->ppms and damps->width_ppm for field independence
            ppms = self.dataset.resppm - (freqs * 1000.0 / self.dataset.frequency)
            damps_ppm = (1.0 / (np.pi * (-damps / 1000.0))) / self.dataset.frequency
            # clip damps_ppm to not have (-) numbers
            damps_ppm = np.clip(damps_ppm, 0.001, None)  # min based on Datasim min PPM linewidth in ppms

            comment = []
            comment.append('Exported HLSVDPro water filter results from Vespa-Analysis Spectral Tab')
            comment.append(' Source:  ' + self.dataset.blocks["raw"].data_source)
            comment.append(' Voxel (x y z):  %d  %d  %d' % (voxel[0], voxel[1], voxel[2]))
            comment.append(' Data points fitted: %d    Number of singular values: %d' % (ndp, nssv))
            comment.append(' ')
            comment.append('HLSVDPro Table Values ')
            comment.append(' Freq [ppm], Freq [kHz], Area, Ph0 [deg], Damp [ms], Damp [ppm], Flag')
            for ppm, fre, amp, pha, damp, damp_ppm, flag in zip(ppms, freqs, ampls, phase, damps, damps_ppm, flags):
                txt = '{0:.2f}, {1:.7f}, {2:.7f}, {3:.3f}, {4:.7f}, {5:.5f}, {6}'.format(ppm, fre, amp, pha, damp, damp_ppm, flag)
                comment.append(txt)
            comment = "\n".join(comment)

            with open(fname, 'w') as f:
                f.write(comment)


    #=======================================================
    #
    #           Widget Event Handlers
    #
    #=======================================================

    def on_tab_changed(self, event):
        voxel = self._tab_dataset.voxel
        ph0 = self.dataset.get_phase_0(voxel)
        ph1 = self.dataset.get_phase_1(voxel)
        # refresh the plot in the current sub-tab
        if self.svd_tab_active:
            view = self.view_svd
            view.set_phase_0(ph0, absolute=True)
            view.set_phase_1(ph1, absolute=True)
        else:
            view = self.view
            view.set_phase_0(ph0, index=[0], absolute=True)
            view.set_phase_1(ph1, index=[0], absolute=True)
            if self.tabB_spectral:
                ph0b = self.datasetB.get_phase_0(voxel)
                ph1b = self.datasetB.get_phase_1(voxel)
                view.set_phase_0(ph0b, index=[1], absolute=True)
                view.set_phase_1(ph1b, index=[1], absolute=True)

        self._tab_dataset.FloatScale.SetValue(view.vertical_scale)


    def on_splitter(self, event=None):
        # This is sometimes called programmatically, in which case event is None
        self._prefs.sash_position_main = self.SplitterWindow.GetSashPosition()
        self._prefs.sash_position_svd = self.SplitterWindowSvd.GetSashPosition()


    def on_combo_dataB(self, event):
        val = event.GetEventObject().GetStringSelection()
        if val == "" or val == "None":
            self._tab_dataset.indexAB[1] = None
            self._tab_dataset.sync = False
            self.CheckSync.SetValue(False)
            self.CheckSync.Disable()
            self.view.change_naxes(1)
        else:
            self._tab_dataset.indexAB[1] = val
            if self._prefs.plot_c_function_none:
                self.view.change_naxes(2)
            else:
                self.view.change_naxes(3)
            self.CheckSync.Enable()
            self.process_and_plot('all', 1)


    def on_flip(self, event):
        # flip respects the sync A/B setting
        value = event.GetEventObject().GetValue()
        nb = self.top.notebook_datasets
        poll_labels = [self._tab_dataset.indexAB[0]]
        if self.do_sync:
            poll_labels = [self._tab_dataset.indexAB[0],self._tab_dataset.indexAB[1]]
        nb.global_poll_sync_event(poll_labels, value, event='flip')

    def set_flip(self, value):
        self.block.set.flip = value
        self.CheckFlip.SetValue(value)
        self.process_and_plot()


    def on_fft(self, event):
        # FFT is always the same across datasets regardless of sync A/B
        value = event.GetEventObject().GetValue()
        self.block.set.fft = value
        if self.do_sync:
            self.datasetB.blocks["spectral"].set.fft = value
            self.tabB_spectral.CheckFFT.SetValue(value)
            self.tabB_spectral.process_and_plot()
        self.process_and_plot()


    def on_chop(self, event):
        # chop respects the sync A/B setting
        value = event.GetEventObject().GetValue()
        nb = self.top.notebook_datasets
        poll_labels = [self._tab_dataset.indexAB[0]]
        if self.do_sync:
            poll_labels = [self._tab_dataset.indexAB[0],self._tab_dataset.indexAB[1]]
        nb.global_poll_sync_event(poll_labels, value, event='chop')


    def on_sync(self, event):
        """
        This check should only be turned on if there is a dataset selected for
        plot B. Otherwise this is always turned off.
        """
        value = event.GetEventObject().GetValue()
        self._tab_dataset.sync = value


    def on_ecc_method(self, event=None):
        # ecc method ignores the sync A/B setting
        # This event handler is sometimes called programmatically, in which
        # case event is None. Don't rely on its existence!
        label = event.GetEventObject().GetStringSelection()
        self.block.set.ecc_method = label

        self.top.Freeze()
        if label == 'None':
            self.PanelEccBrowse.Hide()
        else:
            self.PanelEccBrowse.Show()

        self.top.Layout()
        self.PanelSpectral.Layout()
        self.top.Thaw()

        self.process_and_plot()


    def on_ecc_browse(self, event):
        # Allows the user to select an ECC dataset.
        dialog = dialog_dataset_browser.DialogDatasetBrowser(self.top.datasets)
        dialog.ShowModal()
        ecc_dataset = dialog.dataset
        dialog.Destroy()

        if ecc_dataset:
            self.dataset.set_associated_dataset_ecc(ecc_dataset)
            self.TextEccFilename.SetValue(ecc_dataset.blocks["raw"].data_source)
            self.process_and_plot()


    def on_water_method(self, event=None):
        # water filter type ignores the sync A/B setting
        # This event handler is sometimes called programmatically, in which
        # case event is None. Don't rely on its existence!
        label = event.GetEventObject().GetStringSelection()
        self.block.set.water_filter_method = label

        self.top.Freeze()
        if label == 'None' or label == 'SVD - water filter':
            self.PanelWaterFir.Hide()
            self.PanelWaterHamming.Hide()
        elif label == 'FIR - water filter':
            self.PanelWaterFir.Show()
            self.PanelWaterHamming.Hide()
        elif label == 'Hamming - water filter':
            self.PanelWaterFir.Hide()
            self.PanelWaterHamming.Show()

        self.top.Layout()
        self.PanelSpectral.Layout()
        self.top.Thaw()
        self.process_and_plot()


    def on_fir_length(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.fir_length = value
        self.process_and_plot()

    def on_fir_width(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.fir_half_width = value
        self.process_and_plot()

    def on_fir_ripple(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.fir_ripple = value
        self.process_and_plot()

    def on_fir_extrap_method(self, event):
        value = event.GetEventObject().GetStringSelection()
        self.block.set.fir_extrapolation_method = value
        flag = self.block.set.fir_extrapolation_method == 'AR Model'
        self.SpinFirExtrapValue.Enable(flag)
        self.process_and_plot()

    def on_fir_extrap_value(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.fir_extrapolation_point_count = value
        self.process_and_plot()


    def on_ham_length(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.ham_length = value
        self.process_and_plot()

    def on_ham_extrap_method(self, event):
        value = event.GetEventObject().GetStringSelection()
        self.block.set.ham_extrapolation_method = value
        flag = self.block.set.ham_extrapolation_method == 'AR Model'
        self.SpinHamExtrapValue.Enable(flag)
        self.process_and_plot()

    def on_ham_extrap_value(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.ham_extrapolation_point_count = value
        self.process_and_plot()


    def on_frequency_shift_lock(self, event):
        # frequency shift lock respects the sync A/B setting
        value = event.GetEventObject().GetValue()
        nb = self.top.notebook_datasets
        poll_labels = [self._tab_dataset.indexAB[0]]
        if self.do_sync:
            poll_labels = [self._tab_dataset.indexAB[0],self._tab_dataset.indexAB[1]]
        nb.global_poll_sync_event(poll_labels, value, event='frequency_shift_lock')

    def set_frequency_shift_lock(self, value):
        self.block.frequency_shift_lock = value
        self.CheckFreqLock.SetValue(value)


    def on_phase_lock(self, event):
        # phase 0 lock respects the sync A/B setting
        value = event.GetEventObject().GetValue()
        nb = self.top.notebook_datasets
        poll_labels = [self._tab_dataset.indexAB[0]]
        if self.do_sync:
            poll_labels = [self._tab_dataset.indexAB[0],self._tab_dataset.indexAB[1]]
        nb.global_poll_sync_event(poll_labels, value, event='phase_lock')

    def set_phase_lock(self, value):
        self.block.phase_lock = value
        self.CheckPhaseLock.SetValue(value)


    def on_phase1_zero(self, event):
        # phase 1 zero respects the sync A/B setting
        value = event.GetEventObject().GetValue()
        nb = self.top.notebook_datasets
        poll_labels = [self._tab_dataset.indexAB[0]]
        if self.do_sync:
            poll_labels = [self._tab_dataset.indexAB[0],self._tab_dataset.indexAB[1]]
        nb.global_poll_sync_event(poll_labels, value, event='phase1_zero')

    def set_phase1_zero(self, value):
        self.block.phase_1_lock_at_zero = value
        voxel = self._tab_dataset.voxel
        self.dataset.set_phase_1(0.0, voxel)
        self.FloatPhase1.SetValue(0.0)
        self.CheckZeroPhase1.SetValue(value)
        self.process_and_plot( )


    def on_left_shift_correction(self, event):
        # left shift correction respects the sync A/B setting
        value = event.GetEventObject().GetValue()
        nb = self.top.notebook_datasets
        poll_labels = [self._tab_dataset.indexAB[0]]
        if self.do_sync:
            poll_labels = [self._tab_dataset.indexAB[0],self._tab_dataset.indexAB[1]]
        nb.global_poll_sync_event(poll_labels, value, event='left_shift_correction')

    def set_left_shift_correction(self, value):
        self.block.set.left_shift_correct = value
        self.CheckCorrectPhase1.SetValue(value)
        self.process_and_plot()


    def on_zero_fill(self, event):
        # zero fill multiplier is always synched across spectra
        zf_value = int(event.GetEventObject().GetStringSelection())

        # reset ALL dataset then ALL gui stuff as needed
        self._tab_dataset._outer_notebook.global_block_zerofill_update(zf_value)
        self._tab_dataset._outer_notebook.global_tab_zerofill_update(zf_value)


    def on_apodization_method(self, event):
        # apodization type ignores the synch A/B setting
        index = event.GetEventObject().GetSelection()
        apodization = list(constants.Apodization.choices.keys())[index]
        self.block.set.apodization = apodization

        # Enable the value textbox if a method is selected
        self.FloatWidth.Enable(bool(apodization))
        self.process_and_plot()


    def on_apodization_value(self, event):
        # apodization width ignores the sync A/B setting
        value = event.GetEventObject().GetValue()
        self.block.set.apodization_width = value
        self.process_and_plot()


    def on_b0_shift(self, event):
        # frequency shift respects the sync A/B setting
        value = event.GetEventObject().GetValue()
        voxel = self._tab_dataset.voxel
        orig = self.dataset.get_frequency_shift(voxel)
        # we use the notebook level method to deal with this change because it
        # covers all the actions that need to be taken for manual changes
        poll_labels = [self._tab_dataset.indexAB[0]]
        if self.do_sync:
            poll_labels = [self._tab_dataset.indexAB[0],self._tab_dataset.indexAB[1]]
        self.top.notebook_datasets.global_poll_frequency_shift(poll_labels, value-orig, voxel)


    def on_amplitude(self, event):
        # amplitude multiplier ignores the sync A/B setting
        value = event.GetEventObject().GetValue()
        self.block.set.amplitude = value
        self.process_and_plot()

    def on_phase0(self, event):
        # phase 0 respects the sync A/B setting
        value = event.GetEventObject().GetValue()
        voxel = self._tab_dataset.voxel
        orig = self.dataset.get_phase_0(voxel)
        # we use the notebook level method to deal with this change because it
        # covers all the actions that need to be taken for manual changes
        poll_labels = [self._tab_dataset.indexAB[0]]
        if self.do_sync:
            poll_labels = [self._tab_dataset.indexAB[0],self._tab_dataset.indexAB[1]]
        self.top.notebook_datasets.global_poll_phase(poll_labels, value-orig, voxel, do_zero=True)

    def on_phase1(self, event):
        # phase 1 respects the sync A/B setting
        value = event.GetEventObject().GetValue()
        voxel = self._tab_dataset.voxel
        orig = self.dataset.get_phase_1(voxel)
        # we use the notebook level method to deal with this change because it
        # covers all the actions that need to be taken for manual changes
        poll_labels = [self._tab_dataset.indexAB[0]]
        if self.do_sync:
            poll_labels = [self._tab_dataset.indexAB[0],self._tab_dataset.indexAB[1]]
        self.top.notebook_datasets.global_poll_phase(poll_labels, value-orig, voxel, do_zero=False)


    def on_phase1_pivot(self, event):
        # phase 1 pivot respects the sync A/B setting
        value = event.GetEventObject().GetValue()
        nb = self.top.notebook_datasets
        poll_labels = [self._tab_dataset.indexAB[0]]
        if self.do_sync:
            poll_labels = [self._tab_dataset.indexAB[0],self._tab_dataset.indexAB[1]]
        nb.global_poll_sync_event(poll_labels, value, event='phase1_pivot')

    def set_phase1_pivot(self, value):
        self.block.set.phase_1_pivot = value
        self.FloatPhase1Pivot.SetValue(value)
        self.process_and_plot()


    def on_dc_offset(self, event):
        # DC offset respects the sync A/B setting
        value = event.GetEventObject().GetValue()
        nb = self.top.notebook_datasets
        poll_labels = [self._tab_dataset.indexAB[0]]
        if self.do_sync:
            poll_labels = [self._tab_dataset.indexAB[0],self._tab_dataset.indexAB[1]]
        nb.global_poll_sync_event(poll_labels, value, event='dc_offset')

    def set_dc_offset(self, value):
        self.block.set.dc_offset = value
        self.FloatDcOffset.SetValue(value)
        self.process_and_plot()


    def on_left_shift_value(self, event):
        # left shift points respects the sync A/B setting
        value = event.GetEventObject().GetValue()
        nb = self.top.notebook_datasets
        poll_labels = [self._tab_dataset.indexAB[0]]
        if self.do_sync:
            poll_labels = [self._tab_dataset.indexAB[0],self._tab_dataset.indexAB[1]]
        nb.global_poll_sync_event(poll_labels, value, event='left_shift')

    def set_left_shift(self, value):
        self.block.set.left_shift_value = value
        self.SpinLeftShift.SetValue(value)
        self.process_and_plot()


    def on_user_function(self, event):
        label = event.GetEventObject().GetLabel()

        if label == 'Do Automatic Phasing':

            freq = self.plot_results['freq']
            phase = self.dataset.automatic_phasing_max_real_freq(freq)
            voxel = self._tab_dataset.voxel
            self.dataset.set_phase_0(phase, voxel)
            self.FloatPhase0.SetValue(phase)
            self.process_and_plot()

        elif label == 'Output Area Value':

            ini_name = "spectral_output_area_as_csv"
            default_path_name = util_analysis_config.get_path(ini_name)
            default_path, default_fname = os.path.split(default_path_name)
            filetype_filter = "CSV (*.csv)|*.csv"

            filename = common_dialogs.save_as(default_path=default_path,
                                              filetype_filter=filetype_filter,
                                              default_filename=default_fname)
            if filename:
                # Create output header and results strings, check element count.
                # If the file exists, check that the element count is the same in
                # in the last line as for this results line. If it is, just write
                # out the results string. If different length, output both the
                # header and results strings.

                voxel = self._tab_dataset.voxel
                dfname = self.dataset.dataset_filename
                if not dfname:
                    dfname = self.dataset.blocks['raw'].data_sources[voxel[0]]
                dpath, dfname = os.path.split(dfname)
                area, rms = self.view.calculate_area()
                refs      = self.view.ref_locations_ppm()
                if   self._prefs.area_calc_plot_a:
                    labl = 'A'
                    area = area[0]
                    rms  = rms[0]
                elif self._prefs.area_calc_plot_b:
                    labl = 'B'
                    area = area[1]
                    rms  = rms[1]
                elif self._prefs.area_calc_plot_c:
                    labl = 'C'
                    area = area[2]
                    rms  = rms[2]

                xstr = str(voxel[0])
                ystr = str(voxel[1])
                zstr = str(voxel[2])
                sppm = "%.4f" % refs[0]
                eppm = "%.4f" % refs[1]

                hdr = ['Dataset File', 'x', 'y', 'z', 'Plot', 'Start PPM', 'End PPM', 'Area', 'RMS', 'Path']
                val = [dfname, xstr, ystr, zstr, labl, sppm, eppm, str(area), str(rms), dpath]

                nhdr = len(hdr)
                val = ",".join(val)
                hdr = ",".join(hdr)
                val += "\n"
                hdr += "\n"

                hdr_flag = True
                if os.path.isfile(filename):
                    with open(filename, 'r+') as f:
                        data = f.readlines()
                        if len(data)>1:
                            last = data[-1]
                            nlast = len(last.split(','))
                            if nlast == nhdr:
                                hdr_flag = False

                with open(filename, 'a') as f:
                    if hdr_flag:
                        f.write(hdr)
                    f.write(val)

                # We saved results, so we write the path to the INI file.
                util_analysis_config.set_path(ini_name, filename)



    # start SVD tab event handlers --------------------------------------------

    def on_reset_all(self, event):

        msg = "This will set parameters to default values and erase \nall results. Are you sure you want to continue?"
        if wx.MessageBox(msg, "Reset All Voxels", wx.YES_NO, self) == wx.YES:
            # set all results to 0, turn off all lines. Note. when we replot
            # the current voxel, this will calculate a result for that voxel
            self.block.set_dims(self.dataset)
            self.on_voxel_change(self._tab_dataset.voxel)
            self.process_and_plot()

    def on_slider_changed(self, event):
        # One of the sliders was changed. Each change requires re-applying
        # the HLSVD algorithm.
        # We allow the control to update itself before performing the HLSVD.
        # If we don't, then there can be a noticeable & confusing pause
        # between interacting with the control and seeing it actually change.
        wx.CallAfter(self._apply_hlsvd)

    def on_check_item(self, listctrl, index, item_text, flag):
        # Clicking on a check box in the table causes both the Threshold and
        # Exclude Lipid automatic calculations to be turned off. If you are
        # clicking you obviously want a 'manual' mode. No other boxes are set
        # or unset when the Threshold or Exclude Lipid boxes are unchecked,
        # but those algorithms don't run when their flag are off.

        voxel = self._tab_dataset.voxel
        # because the user can sort the list using any column, we need to
        # get the rank value for the item at index that is causing the event,
        # this is the actual index of the line in the block
        block_index = int(item_text)-1

        svd_output = self.block.get_svd_output(voxel)
        svd_output.in_model[block_index] = flag

        self.block.set.svd_apply_threshold = False
        if self.RadioSvdApplyThreshold.GetValue():
            self.RadioSvdManual.SetValue(True)
        self.FloatSvdThreshold.Disable()
        self.ComboSvdThresholdUnit.Disable()

        self.block.set.svd_exclude_lipid = False
        self.CheckSvdExcludeLipid.SetValue(False)
        self.FloatSvdExcludeLipidStart.Disable()
        self.FloatSvdExcludeLipidEnd.Disable()

        self.process_and_plot()

    def on_svd_manual(self, event):
        self.cursor_span_picks_lines = False
        self.block.set.svd_apply_threshold = False
        self.FloatSvdThreshold.Disable()
        self.ComboSvdThresholdUnit.Disable()

    def on_svd_cursor_span_picks_lines(self, event):
        self.cursor_span_picks_lines = True
        self.block.set.svd_apply_threshold = False
        self.FloatSvdThreshold.Disable()
        self.ComboSvdThresholdUnit.Disable()

    def on_svd_apply_threshold(self, event):
        self.cursor_span_picks_lines = False
        self.block.set.svd_apply_threshold = True
        self.FloatSvdThreshold.Enable()
        self.ComboSvdThresholdUnit.Enable()
        self.process_and_plot()

    def on_svd_threshold(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.svd_threshold = value
        if self.block.set.svd_threshold_unit == 'PPM':
            dataset = self.dataset
            dim0, dim1, dim2, dim3 = dataset.spectral_dims
            maxppm  = dataset.pts2ppm(0)
            minppm  = dataset.pts2ppm(dim0-1)
            # threshold value used to be in Hz and allowed to range +/- 200
            #  units can now be PPM, so check if we are outside min/max ppm
            val = self.block.set.svd_threshold
            if val < minppm:
                self.block.set.svd_threshold = minppm
                self.FloatSvdThreshold.SetValue(self.block.set.svd_threshold)
            elif val > maxppm:
                self.block.set.svd_threshold = maxppm
                self.FloatSvdThreshold.SetValue(self.block.set.svd_threshold)
        self.process_and_plot()

    def on_svd_threshold_unit(self, event):
        index = event.GetEventObject().GetSelection()
        item = list(constants.SvdThresholdUnit.choices.keys())[index]
        self.block.set.svd_threshold_unit = item

        if self.block.set.svd_threshold_unit == 'PPM':
            dataset = self.dataset
            dim0, dim1, dim2, dim3 = dataset.spectral_dims
            maxppm  = dataset.pts2ppm(0)
            minppm  = dataset.pts2ppm(dim0-1)
            # threshold value used to be in Hz and allowed to range +/- 200
            #  units can now be PPM, so check if we are outside min/max ppm
            if self.block.set.svd_threshold < minppm:
                self.block.set.svd_threshold = minppm
                self.FloatSvdThreshold.SetValue(self.block.set.svd_threshold)
            elif self.block.set.svd_threshold > maxppm:
                self.block.set.svd_threshold = maxppm
                self.FloatSvdThreshold.SetValue(self.block.set.svd_threshold)
        self.process_and_plot()

    def on_svd_exclude_lipid(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.svd_exclude_lipid = value
        if value:
            self.FloatSvdExcludeLipidStart.Enable()
            self.FloatSvdExcludeLipidEnd.Enable()
        else:
            self.FloatSvdExcludeLipidStart.Disable()
            self.FloatSvdExcludeLipidEnd.Disable()
        self.process_and_plot()

    def on_svd_exclude_lipid_start(self, event):
        # Note. min=End and max=Start because dealing with PPM range
        min, max = _paired_event(self.FloatSvdExcludeLipidEnd,
                                 self.FloatSvdExcludeLipidStart)
        self.block.set.svd_exclude_lipid_start = max 
        self.block.set.svd_exclude_lipid_end   = min
        self.process_and_plot()

    def on_svd_exclude_lipid_end(self, event):
        # Note. min=End and max=Start because dealing with PPM range
        min, max = _paired_event(self.FloatSvdExcludeLipidEnd,
                                 self.FloatSvdExcludeLipidStart)
        self.block.set.svd_exclude_lipid_start = max 
        self.block.set.svd_exclude_lipid_end   = min
        self.process_and_plot()

    def on_all_on(self, event):
        # Spectral water filter HLSVD threshold value will take precedence
        # over this setting if it is on. The change in the in_model array takes
        # place in the process() call to run the chain.
        voxel = self._tab_dataset.voxel
        svd_output = self.block.get_svd_output(voxel)
        svd_output.in_model.fill(True)

        self.block.set.svd_apply_threshold = False
        self.FloatSvdThreshold.Disable()
        if self.RadioSvdApplyThreshold.GetValue():
            self.RadioSvdManual.SetValue(True)

        self.block.set.svd_exclude_lipid = False
        self.CheckSvdExcludeLipid.SetValue(False)
        self.FloatSvdExcludeLipidStart.Disable()
        self.FloatSvdExcludeLipidEnd.Disable()

#        self.set_check_boxes()

        self.process_and_plot()


    def on_all_off(self, event):
        # Spectral water filter HLSVD threshold value will take precedence
        # over this setting if it is on. The change in the in_model array takes
        # place in the process() call to run the chain.
        voxel = self._tab_dataset.voxel
        svd_output = self.block.get_svd_output(voxel)
        svd_output.in_model.fill(False)

        self.block.set.svd_apply_threshold = False
        self.FloatSvdThreshold.Disable()
        if self.RadioSvdApplyThreshold.GetValue():
            self.RadioSvdManual.SetValue(True)

        self.block.set.svd_exclude_lipid = False
        self.CheckSvdExcludeLipid.SetValue(False)
        self.FloatSvdExcludeLipidStart.Disable()
        self.FloatSvdExcludeLipidEnd.Disable()

#        self.set_check_boxes()

        self.process_and_plot()









    #=======================================================
    #
    #           Public Methods
    #
    #=======================================================

    def process_and_plot(self, entry='all',
                               dataset_to_process=(0, 1),
                               no_draw=False):
        """
        The process(), plot() and process_and_plot() methods are standard in
        all processing tabs. They are called to update the data in the plot
        results dictionary, the plot_panel in the View side of the tab or both.

        """
        tab_base.Tab.process_and_plot(self, entry)

        if self._plotting_enabled:
            self.process(entry, dataset_to_process)
            self.plot(no_draw=no_draw)
            self.plot_svd(no_draw=no_draw)



    def process(self, entry='all', dataset_to_process=(0,1)):
        """
        Data processing results are stored into the Block inside the Chain,
        but the View results are returned as a dictionary from the Chain.run()
        method. The plot routine takes its inputs from this dictionary.

        The dataset_to_process param can be a single integer or a tuple/list
        in the range (0, 1). It defaults to (0, 1). Since the Spectral tab has
        the option to compare a "datasetB" to the "dataset" in the tab, and
        datasetB can also be synched to the tab's dataset, we may need to
        process both dataset and datasetB for each call to process(). The
        parameter "dataset_to_process" can be set to process either dataset or
        datasetB or both by setting it to a tuple of (0,) or (1,) or (0,1)
        respectively
        """
        tab_base.Tab.process(self, entry)

        if self._plotting_enabled:

            if not util_misc.is_iterable(dataset_to_process):
                # Make it a tuple
                dataset_to_process = (dataset_to_process,)

            # update plot results arrays if required
            voxel = self._tab_dataset.voxel

            for i in dataset_to_process:
                dataset = self.dataset if i==0 else self.datasetB
                tab = self if i==0 else self.tabB_spectral
                if dataset:
                    block = dataset.blocks["spectral"]
                    do_fit = block.get_do_fit(voxel)
                    tab.plot_results = block.chain.run([voxel], entry=entry)

                    if dataset == self.dataset:
                        # refresh the hlsvd sub-tab on the active dataset tab
                        if do_fit:
                            # we changed results, now need to update results widget
                            self.svd_checklist_update()
                        else:
                            # same results, but check boxes may have changed
                            self.set_check_boxes()
                            
                        if not self.top.notebook_datasets.suspend_pubsub:
                            # need this check for global zerofill update
                            pubsub.sendMessage("check_datasetb_status",
                                               label=self._tab_dataset.indexAB[0])


    def plot_svd(self, no_draw=False):
        """
        The set_data() method sets data into the plot_panel_spectrum object
        in the plot in the right panel.

        """
        if self._plotting_enabled:
            if not self._tab_dataset.indexAB[0]:
                return

            voxel = self._tab_dataset.voxel
            results = self.plot_results

            data1 = {'data' : results['svd_data'],
                     'line_color_real'      : self._prefs.line_color_real,
                     'line_color_imaginary' : self._prefs.line_color_imaginary,
                     'line_color_magnitude' : self._prefs.line_color_magnitude }

            data2 = {'data' : results['svd_peaks_checked'],
                     'line_color_real'      : self._prefs.line_color_svd,
                     'line_color_imaginary' : self._prefs.line_color_svd,
                     'line_color_magnitude' : self._prefs.line_color_svd }

            data3 = {'data' : results['svd_data'] - results['svd_peaks_checked_sum'],
                     'line_color_real'      : self._prefs.line_color_real,
                     'line_color_imaginary' : self._prefs.line_color_imaginary,
                     'line_color_magnitude' : self._prefs.line_color_magnitude }

            data = [[data1], [data2], [data3]]
            self.view_svd.set_data(data)
            self.view_svd.update(no_draw=True, set_scale=not self._svd_scale_intialized)

            if not self._svd_scale_intialized:
                self._svd_scale_intialized = True

            ph0 = self.dataset.get_phase_0(voxel)
            ph1 = self.dataset.get_phase_1(voxel)
            self.view_svd.set_phase_0(ph0, absolute=True, no_draw=True)
            self.view_svd.set_phase_1(ph1, absolute=True)



    def plot(self, no_draw=False):
        """
        The set_data() method sets data into the plot_panel_spectrum object
        in the plot in the right panel.

        """
        tab_base.Tab.plot(self)

        if self._plotting_enabled:

            voxel   = self._tab_dataset.voxel
            results = self.plot_results
            data1   = results['freq']
            ph0_1   = self.dataset.get_phase_0(voxel)
            ph1_1   = self.dataset.get_phase_1(voxel)

            if self.datasetB:
                resultsB = self.tabB_spectral.plot_results
                data2    = resultsB['freq']
                ph0_2    = self.datasetB.get_phase_0(voxel)
                ph1_2    = self.datasetB.get_phase_1(voxel)
            else:
                data2 = np.zeros_like(data1)
                ph0_2 = 0.0
                ph1_2 = 0.0

            data3 = np.zeros_like(data1)

            # these data will use default line colors in view  data1 == data2
            data = [[data1], [data2], [data3]]
            self.view.set_data(data)
            self.view.update(no_draw=True, set_scale=not self._scale_intialized)

            if not self._scale_intialized:
                self._scale_intialized = True

            # we take this opportunity to ensure that our phase values reflect
            # the values in the block.
            self.view.set_phase_0(ph0_1, absolute=True, no_draw=True, index=[0])
            self.view.set_phase_1(ph1_1, absolute=True, no_draw=True, index=[0])
            self.view.set_phase_0(ph0_2, absolute=True, no_draw=True, index=[1])
            self.view.set_phase_1(ph1_2, absolute=True, no_draw=True, index=[1])

            self.set_plot_c()

            self.view.canvas.draw()

            # Calculate the new area after phasing
            area, rms = self.view.calculate_area()
            if self._prefs.area_calc_plot_a:
                index = 0
                labl = 'A'
            elif  self._prefs.area_calc_plot_b:
                index = 1
                labl = 'B'
            elif  self._prefs.area_calc_plot_c:
                index = 2
                labl = 'C'
            self.top.statusbar.SetStatusText(self.build_area_text(area[index], rms[index], plot_label=labl), 3)


    def set_check_boxes(self):
        """
        This method only refreshes the checkboxes in the current checklist
        widget.

        NB. Bind/Unbind prevents process_plot from being called for each CheckItem()

        """
        self.list_svd_results.Unbind(wx.EVT_LIST_ITEM_CHECKED, handler=self.list_svd_results.OnCheckItem)
        self.list_svd_results.Unbind(wx.EVT_LIST_ITEM_UNCHECKED, handler=self.list_svd_results.OnCheckItem)

        voxel   = self._tab_dataset.voxel
        num     = self.list_svd_results.GetItemCount()
        svd_output = self.block.get_svd_output(voxel)
        for i in range(num):
            if i < self.list_svd_results.GetItemCount():
                index = int(self.list_svd_results.GetItemText(i))-1
                self.list_svd_results.CheckItem(i, svd_output.in_model[index])

        self.list_svd_results.Bind(wx.EVT_LIST_ITEM_CHECKED, self.list_svd_results.OnCheckItem)
        self.list_svd_results.Bind(wx.EVT_LIST_ITEM_UNCHECKED, self.list_svd_results.OnCheckItem)


    def set_chop(self, value):
        self.block.set.chop = value
        self.CheckChop.SetValue(value)
        self.process_and_plot()


    def set_plot_c(self):
        data1 = self.view.all_axes[0].lines[0].get_ydata()
        data2 = self.view.all_axes[1].lines[0].get_ydata()
        data3 = self.plot_C_function(data1, data2)
        self.view.all_axes[2].lines[0].set_ydata(data3)

        view_data1 = self.view.get_data(0)
        view_data2 = self.view.get_data(1)
        view_data3 = self.plot_C_function(view_data1, view_data2)
        self.view.set_data_direct(view_data3, 2)

    def svd_checklist_update(self):
        """
        This method totally rebuilds the result set in the checklist widget.

        Take the hlsvd results for the current voxel and set them into the
        checklist widget. If the spectral tab hlsvd water filter is on and
        has checked the "Apply Threshold" box, then we modify (and save) the
        index values according to the threshold.

        """
        voxel = self._tab_dataset.voxel

        svd_output = self.block.get_svd_output(voxel)

        amp = svd_output.amplitudes
        dam = svd_output.damping_factors
        fre = svd_output.frequencies
        pha = svd_output.phases
        in_model = svd_output.in_model
        nsvd = len(svd_output)

        ppm = self.dataset.resppm - (fre*1000.0/self.dataset.frequency)

        # Update the list_svd_results widget
        if sum(amp) == 0:
            in_model.fill(False)

        res = {}
        self.list_svd_results.DeleteAllItems()
        for i in range(nsvd):
#            res[i] = (i + 1, ppm[i], fre[i]*1000, dam[i], pha[i], amp[i])
#            index = self.list_svd_results.InsertItem(sys.maxsize, ' '+str(res[i][0])+' ')    # bjs_ccx
            res[i] = (i + 1, float(ppm[i]), float(fre[i]*1000), float(dam[i]), float(pha[i]), float(amp[i]))
            index = self.list_svd_results.InsertItem(i,' '+str(res[i][0])+' ')    
            self.list_svd_results.SetItemImage(index, in_model[i])
            self.list_svd_results.SetItem(index, 1, '%.2f'%(res[i][1])) # bjs_ccx
            self.list_svd_results.SetItem(index, 2, '%.1f'%(res[i][2])) # bjs_ccx
            self.list_svd_results.SetItem(index, 3, '%.1f'%(res[i][3])) # bjs_ccx
            self.list_svd_results.SetItem(index, 4, '%.1f'%(res[i][4])) # bjs_ccx
            self.list_svd_results.SetItem(index, 5, '%.1f'%(res[i][5])) # bjs_ccx
            self.list_svd_results.SetItemData(index, i)

        self.list_svd_results.itemDataMap = res


    def on_voxel_change(self, voxel):
        # this just updates widgets that vary based on the voxel number
        # selection. We do not update plot here because that is only done
        # for the active tab in the inner notebook.
        freq_shift = self.dataset.get_frequency_shift(voxel)
        phase_0    = self.dataset.get_phase_0(voxel)
        phase_1    = self.dataset.get_phase_1(voxel)
        self.FloatFrequency.SetValue(freq_shift)
        self.FloatPhase0.SetValue(phase_0)
        self.FloatPhase1.SetValue(phase_1)

        self.SliderDataPoints.SetValue(self.block.get_data_point_count(voxel))
        self.SliderSingularValues.SetValue(self.block.get_signal_singular_value_count(voxel))

        self.svd_checklist_update()







    #=======================================================
    #
    #           Internal Helper Functions
    #
    #=======================================================

    def _apply_hlsvd(self):
        # Updates the plot in response to changes in the HLSVD inputs.
        # This exists just so that we can call it via wx.CallAfter().
        voxel = self._tab_dataset.voxel

        n_data_points = self.SliderDataPoints.GetValue()
        n_singular_values = self.SliderSingularValues.GetValue()
        
        self.block.set.svd_last_n_data_points = n_data_points
        self.block.set.svd_last_n_singular_values = n_singular_values

        self.block.set_data_point_count(n_data_points, voxel)
        self.block.set_signal_singular_value_count(n_singular_values, voxel)
        self.block.set_do_fit(True, voxel)
        self.process_and_plot(no_draw=True)
        self.set_check_boxes()







