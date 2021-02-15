# Python modules
import copy

# 3rd party modules
import wx
import numpy as np

# Our modules
import vespa.analysis.mrs_dataset as mrs_dataset
import vespa.analysis.constants as constants
import vespa.analysis.prefs as prefs_module
import vespa.common.util.misc as util_misc
import vespa.common.wx_gravy.util as wx_util
import vespa.common.wx_gravy.common_dialogs as common_dialogs
import vespa.common.util.ppm as util_ppm

#from vespa.common.wx_gravy.widgets.floatspin_multiplier.floatspin_multiplier_base import FloatSpinMultiplier
from wx.lib.agw.floatspin import FloatSpin, EVT_FLOATSPIN, FS_LEFT, FS_RIGHT, FS_CENTRE, FS_READONLY

from vespa.analysis.plot_panel_user_prior import PlotPanelUserPrior
from vespa.analysis.auto_gui.dialog_user_prior import MyDialog


# _NEW_LINE_DEFAULTS are the default values used when the user adds a line.
_NEW_LINE_DEFAULTS = [True, 4.7, 1.0, 0.0, 5.0, 1.0, 1.0, 1.0, 1.0] 


def _paired_event(obj_min, obj_max):
        val_min = obj_min.GetValue()
        val_max = obj_max.GetValue()
        pmin = min(val_min, val_max)
        pmax = max(val_min, val_max)
        obj_min.SetValue(pmin)
        obj_max.SetValue(pmax)
        return pmin, pmax


# PlotPanel expects a Prefs object, which are intimately tied with View menu.
# This modal dialog uses a cobbled together Prefs object with hardcoded values.

class _UserPriorPrefs(object):
    def __init__(self):
        self.bgcolor = "white"
        self.line_color_imaginary = "red"
        self.line_color_magnitude = "purple"
        self.line_color_real = "black"
        self.zero_line_color = "goldenrod"
        self.zero_line_style = "solid"
        self.line_width = 1.0
        self.line_color_individual = "green"
        self.line_color_summed = "black"
        self.zero_line_bottom = False
        self.zero_line_middle = True
        self.zero_line_top = False
        self.zero_line_show = True
        self.xaxis_ppm = True
        self.xaxis_hertz = False
        self.xaxis_show = True

    @property
    def _ini_section_name(self):
        return "spectral_prefs"



class DialogUserPrior(MyDialog):
    """
    Displays metabolite prior info that can be used to create a simulated
    spectrum and allows the user to edit it. Returns True if the user hits OK
    and False if user hits CANCEL. Changes are written to an internal object
    and must be copied from the dialog after it returns.

    """
    def __init__(self, parent, dataset, prior, show_ph1=True):

        # self.SizerSplitterWindow.Fit(self)         # bugfix wxGlade 0.9.6 to 1.0.0

        super().__init__(parent)

        self.dataset = dataset
        self.prior = copy.deepcopy(prior)       # do not change original!

        # cache the result of self._update_spectrum().
        self._prefs = _UserPriorPrefs()

        # We add buttons dynamically to get right order under OS X, GTK, Windows, etc.
        self.ButtonOpen, self.ButtonCancel = \
            wx_util.add_ok_cancel(self, self.LabelOKCancelPlaceholder, self.on_ok, self.on_cancel)

        self.SetSize((1000, 600))
        self.Layout()
        self.Center()

        # Plotting is disabled during some of init to avoid circular event firings
        self._plotting_enabled = False
        self._initialize_controls()
        self._populate_controls()
        self._plotting_enabled = True

        self.show_ph1(show_ph1)

        # Init data scale to a sane value after passing to the view.
        self._scale_intialized = False
        self._plot()

        # Sash position is hardcoded here, not part of Prefs.
        wx.CallAfter(self.WindowSplitter.SetSashPosition, 500, True)


    #################     Event Handlers   ##################

    def on_ok(self, event):
        self.EndModal(wx.ID_OK)

    def on_cancel(self, event):
        self.EndModal(wx.ID_CANCEL)


    #################     Internal Helpers    ##################
    def _initialize_controls(self):
        """ 
        Set up sizes and constraints for widgets. Does not set values.

        """
        ds = self.dataset
        dim0, dim1, dim2, dim3 = ds.spectral_dims
        ppmlim  = (ds.pts2ppm(dim0-1), ds.pts2ppm(0))
        dim0lim = (0, dim0 - 1)

        # Configure spin controls size, # of digits displayed, increment and min/max
        wx_util.configure_spin(self.FloatAutoB0RangeStart,     70, 3, 0.1, ppmlim)
        wx_util.configure_spin(self.FloatAutoB0RangeEnd,       70, 3, 0.1, ppmlim)
        wx_util.configure_spin(self.FloatAutoPhase0RangeStart, 70, 3, 0.1, ppmlim)
        wx_util.configure_spin(self.FloatAutoPhase0RangeEnd,   70, 3, 0.1, ppmlim)
        wx_util.configure_spin(self.FloatAutoPhase1RangeStart, 70, 3, 0.1, ppmlim)
        wx_util.configure_spin(self.FloatAutoPhase1RangeEnd,   70, 3, 0.1, ppmlim)
        wx_util.configure_spin(self.FloatAutoPhase1Pivot,      70, 3, 0.1, ppmlim)

        self.view = PlotPanelUserPrior( self.PanelView,
                                        naxes=2,
                                        reversex=True,
                                        zoom='span',
                                        reference=True,
                                        middle=True,
                                        do_zoom_select_event=True,
                                        do_zoom_motion_event=True,
                                        do_refs_select_event=True,
                                        do_refs_motion_event=True,
                                        do_middle_select_event=False,
                                        do_middle_motion_event=False,
                                        do_scroll_event=True,
                                        props_zoom=dict(alpha=0.2, facecolor='yellow'),
                                        props_cursor=dict(alpha=0.2, facecolor='gray'),
                                        xscale_bump=0.0,
                                        yscale_bump=0.05,
                                        data = [],
                                        prefs=self._prefs,
                                        dataset=self.dataset   )

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.view, 1, wx.LEFT | wx.TOP | wx.EXPAND)
        self.PanelView.SetSizer(sizer)
        self.view.Fit()    

        #------------------------------------------------------------
        # set up dynamic list
        
        self.LinesGridSizer = self.LabelPriorLinesPlaceholder.GetContainingSizer()
        parent = self.LabelPriorLinesPlaceholder.GetParent()
        self.LabelPriorLinesPlaceholder.Destroy()
       
        # Add headings to the first row of the grid sizer.
        self.LinesGridSizer.Clear()
        self.LinesGridSizer.SetRows(1)
        headings = (None, "Peak Center\n[ppm]", "Peak Area", "Peak Phase\n[deg]", "Linewidth\n[Hz]")
        
        for heading in headings:
            if heading:
                label = wx.StaticText(parent, label=heading, style=wx.ALIGN_CENTRE)
                self.LinesGridSizer.Add(label, 0, wx.EXPAND|wx.ALIGN_CENTER_VERTICAL)
            else:
                self.LinesGridSizer.AddSpacer( 1 )

        self.dynlist1 = DynamicList1(self.PanelLines, self.PanelPrior, self.LinesGridSizer,
                                     self.dataset, self.prior, self.on_dynamic_user_prior_list)




    def _populate_controls(self):
        """ 
        Populates the widgets with relevant values from the data object. 
        It's meant to be called when a new data object is loaded.
        
        This function trusts that the data object it is given doesn't violate
        any rules. Whatever is in the data object gets slapped into the 
        controls, no questions asked. 
        
        This function is, however, smart enough to enable and disable 
        other widgets depending on settings.

        """
        self.FloatAutoB0RangeStart.SetValue(    self.prior.auto_b0_range_start)
        self.FloatAutoB0RangeEnd.SetValue(      self.prior.auto_b0_range_end)
        self.FloatAutoPhase0RangeStart.SetValue(self.prior.auto_phase0_range_start)
        self.FloatAutoPhase0RangeEnd.SetValue(  self.prior.auto_phase0_range_end)
        self.FloatAutoPhase1RangeStart.SetValue(self.prior.auto_phase1_range_start)
        self.FloatAutoPhase1RangeEnd.SetValue(  self.prior.auto_phase1_range_end)
        self.FloatAutoPhase1Pivot.SetValue(     self.prior.auto_phase1_pivot)

        # - set default values into DynList
        # - push defaults to prior object
        self.dynlist1.set_new_values()
        self.prior.basis.set_values(self.dynlist1.lines)
        

    def _plot(self):

        if self._plotting_enabled:

            # Get auto_prior spectral arrays (indiv and summed)
            lines    = self.prior.basis.get_spectrum_all(self.dataset)
            line_sum = self.prior.basis.get_spectrum_sum(self.dataset)
            
            data1 = {'data' : lines, 
                     'line_color_real'      : self._prefs.line_color_individual,
                     'line_color_imaginary' : self._prefs.line_color_individual,
                     'line_color_magnitude' : self._prefs.line_color_individual }
            
            data2 = {'data' : line_sum,
                     'line_color_real'      : self._prefs.line_color_summed,
                     'line_color_imaginary' : self._prefs.line_color_summed,
                     'line_color_magnitude' : self._prefs.line_color_summed }
    
            data = [[data1], [data2]]
            self.view.set_data(data)
            self.view.update(set_scale=not self._scale_intialized)
    
            if not self._scale_intialized:
                self._scale_intialized = True
    
            # Calculate the new area after phasing
            area, rms = self.view.calculate_area()
            area = area[1]
            rms = rms[1]


    def show_ph1(self, flag):
        self.LabelAutoPhase1RangeStart.Show(flag)
        self.LabelAutoPhase1RangeEnd.Show(flag)
        self.LabelPhase1Pivot.Show(flag)
        self.FloatAutoPhase1RangeStart.Show(flag)
        self.FloatAutoPhase1RangeEnd.Show(flag)
        self.FloatAutoPhase1Pivot.Show(flag)


    ##############   Event Handlers   ################

    ########  Handlers for line controls ########

    def on_add_line(self, event):
        # This fires on_dynamic_user_prior_list() so nothing else needed here.
        self.dynlist1.add_row(_NEW_LINE_DEFAULTS, update=True)

    def on_delete_line(self, event):
        # This fires on_dynamic_user_prior_list() so nothing else needed here.
        self.dynlist1.remove_checked_rows()

    def on_restore_defaults(self, event): 
        self.prior.basis.reset_to_default()
        self.dynlist1.set_new_values()
        self.on_dynamic_user_prior_list()

    def on_all_on(self, event):
        self.dynlist1.select_all()

    def on_all_off(self, event):
        self.dynlist1.deselect_all()

    def on_dynamic_user_prior_list(self, event=None):
        # This is a multi-use event handler. It is called from another actual
        # event that occurs in the dynamic list class. But, it can also be
        # invoked programmatically as needed. In the latter case, event is None.

        lines = self.dynlist1.lines
        self.prior.basis.set_values(lines)
        self._plot()

    def reset_to_default(self):
        pass

    ####### Handlers for algorithm params controls #######

    def on_auto_b0_range_start(self, event):
        min_, max_ = _paired_event(self.FloatAutoB0RangeStart, self.FloatAutoB0RangeEnd)
        self.prior.auto_b0_range_start = min_
        self.prior.auto_b0_range_end = max_
        self._plot()

    def on_auto_b0_range_end(self, event):
        min_, max_ = _paired_event(self.FloatAutoB0RangeStart, self.FloatAutoB0RangeEnd)
        self.prior.auto_b0_range_start = min_
        self.prior.auto_b0_range_end = max_
        self._plot()

    def on_auto_phase0_range_start(self, event):
        min_, max_ = _paired_event(self.FloatAutoPhase0RangeStart, self.FloatAutoPhase0RangeEnd)
        self.prior.auto_phase0_range_start = min_
        self.prior.auto_phase0_range_end = max_
        self._plot()

    def on_auto_phase0_range_end(self, event): 
        min_, max_ = _paired_event(self.FloatAutoPhase0RangeStart, self.FloatAutoPhase0RangeEnd)
        self.prior.auto_phase0_range_start = min_
        self.prior.auto_phase0_range_end = max_
        self._plot()

    def on_auto_phase1_range_start(self, event): 
        min_, max_ = _paired_event(self.FloatAutoPhase1RangeStart, self.FloatAutoPhase1RangeEnd)
        self.prior.auto_phase1_range_start = min_
        self.prior.auto_phase1_range_end = max_
        self._plot()

    def on_auto_phase1_range_end(self, event): 
        min_, max_ = _paired_event(self.FloatAutoPhase1RangeStart, self.FloatAutoPhase1RangeEnd)
        self.prior.auto_phase1_range_start = min_
        self.prior.auto_phase1_range_end = max_
        self._plot()

    def on_auto_phase1_pivot(self, event): 
        val = event.GetEventObject().GetValue()
        self.prior.auto_phase1_pivot = val
        self._plot()


#---------------------------------------------------------------------------------------------------

class DynamicList1(object):

    def __init__(self, PanelLines, PanelPrior, GridSizer, dataset, prior, external_event_handler):
        
        self.dataset = dataset
        self.prior = prior
        self._PanelLines = PanelLines
        self._PanelPrior = PanelPrior
        self._GridSizer = GridSizer
        self.external_event_handler = external_event_handler
        self._list_lines = []


    @property
    def lines(self):
        return [self._get_line_values(line) for line in self._list_lines]


    def set_new_values(self, new_rows=None):
        if not new_rows:
            new_rows = self.prior.basis.get_rows()
        self.remove_all_rows()
        for row_vals in new_rows:
            self.add_row(row_vals)


    def add_row(self, row_vals, update=False):
        """ Adds a row to the end of the list.  """

        ppmlim = (self.dataset.pts2ppm(self.dataset.spectral_dims[0] - 1), self.dataset.pts2ppm(0))
        self._GridSizer.SetRows(self._GridSizer.GetRows() + 1)

        # create widgets to go into the line
        checkbox    = wx.CheckBox(self._PanelLines)
        value_ppm   = FloatSpin(self._PanelLines, agwStyle=FS_LEFT)
        value_area  = FloatSpin(self._PanelLines, agwStyle=FS_LEFT)
        value_phase = FloatSpin(self._PanelLines, agwStyle=FS_LEFT)
        value_lwhz  = FloatSpin(self._PanelLines, agwStyle=FS_LEFT)
        
        # keep a copy of panel and widgets to access later
        line = { "check"        : checkbox, 
                 "value_ppm"    : value_ppm, 
                 "value_area"   : value_area, 
                 "value_phase"  : value_phase,
                 "value_lwhz"   : value_lwhz,  }

        # Add the controls to the grid sizer, configure and set values
        self._GridSizer.Add(line["check"], 0, wx.ALIGN_CENTER_VERTICAL)
        for key in ("value_ppm", "value_area", "value_phase", "value_lwhz"):
            self._GridSizer.Add(line[key], 0, wx.EXPAND)

        wx_util.configure_spin(value_ppm,  70, 2, 0.05, ppmlim)
        wx_util.configure_spin(value_area, 70, 3, 0.1, (0.001,100000.0))
        wx_util.configure_spin(value_phase,70, 1, 5.0, (-360,360))
        wx_util.configure_spin(value_lwhz, 70, 2, 1.0, (0.001,10000.0))

        checkbox.SetValue(row_vals[0])
        value_ppm.SetValue(row_vals[1])
        value_area.SetValue(row_vals[2])
        value_phase.SetValue(row_vals[3])
        value_lwhz.SetValue(row_vals[4])

        self._list_lines.append(line)

        # only need to update if a metabolite is added/removed from basis
        self._PanelLines.Bind(wx.EVT_CHECKBOX, self.event_handler, checkbox)
        self._PanelLines.Bind(EVT_FLOATSPIN, self.event_handler, value_ppm)
        self._PanelLines.Bind(EVT_FLOATSPIN, self.event_handler, value_area)
        self._PanelLines.Bind(EVT_FLOATSPIN, self.event_handler, value_phase)
        self._PanelLines.Bind(EVT_FLOATSPIN, self.event_handler, value_lwhz)
        self._PanelPrior.Layout()

        if update: self.event_handler()
        
        
    def remove_checked_rows(self):

        checklist = []
        for i, line in enumerate(self._list_lines):
            if line["check"].GetValue():
                checklist.append(i)
        
        # remove in reverse order to not invalidate indices earlier in list
        checklist.reverse()
        for i in checklist:
            for item in list(self._list_lines[i].values()):
                if hasattr(item, "Destroy"):
                    item.Destroy()              # It's a wx control
            del self._list_lines[i]
            
        # Reduce the # of rows in the grid sizer
        rows = self._GridSizer.GetRows()
        self._GridSizer.SetRows(rows - len(checklist))
        self._GridSizer.Layout()
        self._PanelPrior.Layout()        
        self.event_handler()

    def remove_all_rows(self):
        self.select_all()
        self.remove_checked_rows()

    def select_all(self):
        for line in self._list_lines:
            line["check"].SetValue(True)
        self.event_handler()

    def deselect_all(self):
        for line in self._list_lines:
            line["check"].SetValue(False)
        self.event_handler()

    def event_handler(self, event=None):
        self.external_event_handler(event)


    #######   Private Methods   #########

    def _get_line_values(self, line):
        # Returns a dict containing  values of the controls in the line.
        return { "check"        : line["check"].GetValue(),
                 "ppm"          : line["value_ppm"].GetValue(),
                 "area"         : line["value_area"].GetValue(),
                 "phase"        : line["value_phase"].GetValue(),
                 "lw"           : line["value_lwhz"].GetValue(),
                 "ppm_lim"      : 0.0,
                 "area_lim"     : 0.0,
                 "phase_lim"    : 0.0,
                 "lw_lim"       : 0.0
               }






## Test code  ##########################################################

class MyForm(wx.Frame):
 
    def __init__(self, dataset):

        self.dataset = dataset
        self.prior   = dataset.user_prior

        wx.Frame.__init__(self, None, wx.ID_ANY, "Testing the Dialog")
        panel = wx.Panel(self, wx.ID_ANY)
        buttn = wx.Button(panel, label="Show Dialog")
        buttn.Bind(wx.EVT_BUTTON, self.on_dialog)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(buttn, 0, wx.ALL|wx.CENTER, 5)
        panel.SetSizer(sizer)

        self.statusbar = self.CreateStatusBar(4, 0)
        self.statusbar.SetStatusText("Ready")
 
    def on_dialog(self, event):

        dialog = DialogUserPrior(self, self.dataset, self.dataset.user_prior, show_ph1=True)
        if dialog.ShowModal() == wx.ID_OK:
            self.dataset.user_prior = dialog.prior

        dialog.Destroy()

 
#----------------------------------------------------------------------

if __name__ == "__main__":
    
    dataset = mrs_dataset.Dataset()
    app = wx.App(False)
    frame = MyForm(dataset)
    frame.Show()
    app.MainLoop()