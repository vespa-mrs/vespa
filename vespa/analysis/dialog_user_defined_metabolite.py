# Python modules
import re
import copy

# 3rd party modules
import wx
import numpy as np

# Our modules
import vespa.analysis.mrs_dataset as mrs_dataset
import vespa.common.util.misc as util_misc
import vespa.common.mrs_prior_metabolite as mrs_prior_metabolite
import vespa.common.wx_gravy.util as wx_util
import vespa.common.wx_gravy.common_dialogs as common_dialogs
import vespa.analysis.auto_gui.dialog_user_defined_metabolite as dialog_user_defined_metabolite

from vespa.common.util.generic_spectral import create_spectrum
from vespa.analysis.plot_panel_user_defined_metabolite import PlotPanelUserDefinedMetabolite
from wx.lib.agw.floatspin import FloatSpin, EVT_FLOATSPIN, FS_LEFT, FS_RIGHT, FS_CENTRE, FS_READONLY


# default values used when the user adds a line.
_NEW_LINE_DEFAULTS = [False, 2.0, 1.0, 0.0] 
_NAME_PATTERN = r'[^\.a-zA-Z0-9_-]'


def _configure_combo(control, choices, selection=''):
        if isinstance(choices,dict):
            lines = [choices[key] for key in list(choices.keys())]
        elif isinstance(choices,list):
            lines = [str(item) for item in choices]
        else:
            return
            
        control.Clear()        
        control.AppendItems( lines )
        if selection in lines:
            control.SetStringSelection(selection)
        else: 
            control.SetStringSelection(lines[0])




class DialogUserDefinedMetabolite(dialog_user_defined_metabolite.MyDialog):
    """
    Displays user prior info and allows the user to edit it. Returns True
    if the user changed settings, False otherwise. Changes are written to the
    dataset by this dialog, so the return code is informative only. The caller
    doesn't need to react to it (except perhaps to update the GUI).

    """
    def __init__(self, parent, dataset, metabolite=None, unique_names=[]):

        dialog_user_defined_metabolite.MyDialog.__init__(self, parent)

        self.SizerSplitterWindow.Fit(self)  # bugfix wxGlade 0.9.6 to 1.0.0

        if metabolite is None:
            # set a default, one line, metabolite at water ppm
            metabolite = mrs_prior_metabolite.PriorMetabolite()
            metabolite.dims[0] = util_misc.uuid()
            metabolite.spins  = 1
            metabolite.ppms   = [4.7,]
            metabolite.areas  = [1.0,]
            metabolite.phases = [0.0,]

        self.dataset     = dataset
        self.metabolite  = metabolite
        self.unique_names = unique_names

        # copy original metabolite in case user hits cancel
        self.original   = copy.deepcopy(metabolite)

        # cache result of update_spectrum().
        self.spectrum_all = None
        self.spectrum_sum = None

        self.linewidth = 5.0        # display line broadening value

        # Plot panel expects a Prefs object, but we are in modal dialog,
        # so we cobble together a fake Prefs object with hardcoded values.
        class FakePrefs(object):
            def __init__(self):
                self.foreground_color = "black"
                self.bgcolor = "white"
                self.zero_line_show = True
                self.zero_line_top = False
                self.zero_line_middle = True
                self.zero_line_bottom = False
                self.xaxis_show = True
                self.xaxis_ppm = True
                self.xaxis_hertz = False
                self.xaxis_sec = False
                self.data_type_real = True
                self.data_type_imaginary = False
                self.data_type_magnitude = False
                self.zero_line_color = "goldenrod"
                self.zero_line_style = "solid"
                self.line_color_real = "black"
                self.line_color_imaginary = "red"
                self.line_color_magnitude = "purple"
                self.line_width = 1.0
                self.line_color_individual = "green"
                self.line_color_summed = "black"
                self.plot_view = "all"

        self.prefs = FakePrefs()

        # Open & Cancel buttons added dynamically so that they're in the
        # right order under OS X, GTK, Windows, etc.
        self.ButtonOpen, self.ButtonCancel = \
            wx_util.add_ok_cancel(self, self.LabelOKCancelPlaceholder, self.on_ok, self.on_cancel)

        self.SetSize( (750, 600) )
        self.Layout()
        self.Center()

        # disable plot during init to avoid change events that trigger call
        # to plot(). This is a Windows-specific bug.
        self.plotting_enabled = False
        
        self.initialize_controls()
        self.populate_controls()

        self.plotting_enabled = True

        # one time flag to show plot scale initialize after data passed to view.
        self.scale_intialized = False
        self.update_spectrum()

        # Set up the canvas
        self.plot()

        # Sash position normally in self.prefs; here just hardcoded
        wx.CallAfter(self.SplitterWindow.SetSashPosition, 450, True)




    #################     Internal Helpers    ##################
    def initialize_controls(self):
        """ 
        This methods goes through the widgets and sets up certain sizes
        and constraints for those widgets. This method does not set the 
        value of any widget except insofar that it is outside a min/max
        range as those are being set up. 
        
        Use populate_controls() to set the values of the widgets from
        a data object.

        """
        # calculate a few useful values
        
        dim0, dim1, dim2, dim3 = self.dataset.spectral_dims
        sw      = self.dataset.sw
        maxppm  = self.dataset.pts2ppm(0)
        minppm  = self.dataset.pts2ppm(dim0-1)
        ppmlim  = (minppm, maxppm)
        dim0lim = (0, dim0 - 1)
        
        
        # The many spin controls on various tabs need configuration of 
        # their size, # of digits displayed, increment and min/max. 
        wx_util.configure_spin(self.FloatLinewidth, 70, 1, 1.0, [0.1,500.0])

        _configure_combo(self.ComboZeroLineLocation, ['Top','Middle','Bottom'])


        self.view = PlotPanelUserDefinedMetabolite( self.PanelView, 
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
                                                    yscale_bump=0.0,
                                                    data = [],
                                                    prefs=self.prefs, 
                                                    dataset=self.dataset,
                                                    )
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.view, 1, wx.LEFT | wx.TOP | wx.EXPAND)
        self.PanelView.SetSizer(sizer)
        self.view.Fit()    


        #------------------------------------------------------------
        # set up dynamic list
        
        # The list grid sizer is marked so we can find it at run-time
        self.GridSizerMetabolite = self.LabelMetaboliteLinesPlaceholder.GetContainingSizer()
        parent = self.LabelMetaboliteLinesPlaceholder.GetParent()
        self.LabelMetaboliteLinesPlaceholder.Destroy()
       
        # Add headings to the first row of the grid sizer.
        self.GridSizerMetabolite.Clear()
        self.GridSizerMetabolite.SetRows(1)
        headings = (None, "Location\n[ppm]", "Peak Area\nRatio", "Peak Phase\n[deg]")
        
        for heading in headings:
            if heading:
                label = wx.StaticText(parent, label=heading, style=wx.ALIGN_CENTRE)
                self.GridSizerMetabolite.Add(label, 0, wx.EXPAND|wx.ALIGN_CENTER_VERTICAL)
            else:
                self.GridSizerMetabolite.AddSpacer( 1 )

        self.dynamic_list = DynamicList1( self.PanelMetabolite,
                                           self.PanelUserDefinedMetabolite,
                                           self.GridSizerMetabolite,
                                           self.dataset,
                                           self.metabolite,
                                           self.on_dynamic_list)


    def populate_controls(self):
        """ 
        Populates the widgets with relevant values from the data object. 
        It's meant to be called when a new data object is loaded.
        
        This function trusts that the data object it is given doesn't violate
        any rules. Whatever is in the data object gets slapped into the 
        controls, no questions asked. 
        
        This function is, however, smart enough to enable and disable 
        other widgets depending on settings.
        """
        dataset    = self.dataset
        metabolite = self.metabolite

        #------------------------------------------------------------
        # set up Basic widgets
        user_prior = self.dataset.user_prior

        self.TextMetaboliteName.SetValue( metabolite.dims[0] )
        self.FloatLinewidth.SetValue( self.linewidth )
        self.ComboZeroLineLocation.SetStringSelection( 'Middle' )

        # we do not match old settings to any new ones in this widget
        self.dynamic_list.set_new_values(metabolite, add_check=False)
        self.update_spectrum()
        

    def plot(self):

        if self.plotting_enabled:

            data1 = {'data' : self.spectrum_all, 
                     'line_color_real'      : self.prefs.line_color_individual,
                     'line_color_imaginary' : self.prefs.line_color_individual,
                     'line_color_magnitude' : self.prefs.line_color_individual }
            
            data2 = {'data' : self.spectrum_sum,
                     'line_color_real'      : self.prefs.line_color_summed,
                     'line_color_imaginary' : self.prefs.line_color_summed,
                     'line_color_magnitude' : self.prefs.line_color_summed }
    
            data = [[data1], [data2]]
            self.view.set_data(data)
            self.view.update(set_scale=not self.scale_intialized)
    
            if not self.scale_intialized:
                self.scale_intialized = True
    
            # Calculate the new area after phasing
            area, rms = self.view.calculate_area()
            area = area[1]
            rms = rms[1]


    def update_spectrum(self):
        """
        Calculates the list of spectral lines that are 
        described in the widget and returns them. It also updates the
        value returned by this object's summed property.

        There are two cases where this method will return a single line of 
        zeros, first there are no lines in the widget, second there are lines
        but none have been checked to include.
        """
        dataset = self.dataset
        lines   = self.dynamic_list.lines
        areas   = [item['value_area']  for item in lines]
        ppms    = [item['value_ppm']   for item in lines]
        phases  = [item['value_phase'] for item in lines]
        spectrum_all = []

        if areas:
            for i in range(len(areas)):
                apod = 1.0 / (self.linewidth * 0.34 * np.pi)
                spectrum = create_spectrum( areas[i],ppms[i],phases[i],dataset,10000.0,apod)
                spectrum_all.append(np.array(spectrum))
            spectrum_all = np.array(spectrum_all)
        else:
            zfmult   = dataset.zero_fill_multiplier
            acqdim0  = dataset.raw_dims[0]
            spectrum_all = np.zeros((acqdim0 * zfmult), complex)
            
        if len(spectrum_all.shape)==1:
            spectrum_all.shape = 1, spectrum_all.shape[-1]

        # Update the sum.
        self.spectrum_all = spectrum_all 
        self.spectrum_sum = np.sum(spectrum_all, axis=0)
            
        return
        



##############################    Event Handlers

    def on_ok(self, event):
        """ copy widget values into the PriorMetabolite object passed in """
        mname  = self.TextMetaboliteName.GetValue()

        msg = ''
        if mname in self.unique_names:
            msg = "Metabolite name is not unique, do you want to fix it?\n" \
                  "    If YES, you will return to the dialog.\n" \
                  "    If NO, then your settings will be lost."
        
        if re.search(_NAME_PATTERN, mname):
            "name contains other than a-zA-Z0-9"
            msg = "Metabolite name has invalid characters, do you want to fix it?\n" \
                  " (Please only use alpha-numeric, underscore and hyphens.)\n" \
                  "    If YES, you will return to the dialog.\n" \
                  "    If NO, then your settings will be lost."

        if msg:
            res = common_dialogs.message(msg, "Error Saving Metabolite",wx.ICON_ERROR|wx.YES|wx.NO)
            if res==wx.YES:
                return
            else:
                self.metabolite = self.original
                self.EndModal(wx.ID_CANCEL)
                return
        
        lines  = self.dynamic_list.lines
        areas  = [item['value_area']  for item in lines]
        ppms   = [item['value_ppm']   for item in lines]
        phases = [item['value_phase'] for item in lines]
        
        self.metabolite.dims[0] = mname
        self.metabolite.spins   = len(areas)
        self.metabolite.areas   = areas
        self.metabolite.ppms    = ppms
        self.metabolite.phases  = phases
        
        self.EndModal(wx.ID_OK)


    def on_cancel(self, event):
        # Restore original settings before exiting.
        self.metabolite = self.original
        self.EndModal(wx.ID_CANCEL)


########  Handlers for line controls
    def on_metabolite_name(self, event):
        pass
    
    
    def on_add_line(self, event):
        self.dynamic_list.add_row(_NEW_LINE_DEFAULTS, update=True)
        # The call above fires on_dynamic_list() so I don't need
        # to do anything else here.


    def on_delete_line(self, event): 
        self.dynamic_list.remove_checked_rows()
        # The call above fires on_dynamic_list() so I don't need
        # to do anything else here.


    def on_dynamic_list(self, event=None):
        # This is a fudged event called from the actual event that occurs 
        # inside the dynamic list class but its also invoked programmatically
        # (i.e. by our code, not by wx). In the latter case, event is None.

        self.update_spectrum()
        self.plot()


####### Handlers for preview control params

    def on_linewidth(self, event): 
        val = event.GetEventObject().GetValue()
        self.linewidth = val
        self.update_spectrum()
        self.plot()


    def on_zero_line_location(self, event): 
        val = event.GetEventObject().GetValue()
        self.prefs.zero_line_top = False
        self.prefs.zero_line_middle = False
        self.prefs.zero_line_bottom = False
        if val == 'Top':
            self.prefs.zero_line_top = True
        elif val == 'Middle':
            self.prefs.zero_line_middle = True
        else:
            self.prefs.zero_line_bottom = True
        self.plot()



#------------------------------------------------------------------------------

class DynamicList1(object):

    def __init__(self, PanelParent, PanelTop, GridSizer, dataset, metabolite, external_event_handler):
        
        self.dataset     = dataset
        self.metabolite  = metabolite
        self.PanelParent = PanelParent
        self.PanelTop    = PanelTop
        self.GridSizer = GridSizer
        
        self.external_event_handler = external_event_handler
        self.list_lines = []


    @property
    def lines(self):
        return [self.get_line_values(line) for line in self.list_lines]

        
    def set_new_values(self, previous=None, add_check=None):
        if not previous:
            previous = self.metabolite

        self.select_all()
        self.remove_checked_rows()

        for row in previous.get_rows():
            if add_check is not None:
                if add_check is not True:
                    row.insert(0,False)
                else:
                    row.insert(0,True)
                
            self.add_row(row)

    
    def add_row(self, row_vals, update=False):
        '''
        Adds a row to the end of the list. 
        
        '''
        maxppm = self.dataset.pts2ppm(0)
        minppm = self.dataset.pts2ppm(self.dataset.spectral_dims[0]-1)

        self.GridSizer.SetRows(self.GridSizer.GetRows() + 1)

        # create widgets to go into the line
        list_line = {}

        checkbox     = wx.CheckBox(self.PanelParent)
        value_ppm    = FloatSpin(self.PanelParent, agwStyle=FS_CENTRE)
        value_area   = FloatSpin(self.PanelParent, agwStyle=FS_CENTRE)
        value_phase  = FloatSpin(self.PanelParent, agwStyle=FS_CENTRE)
        
        # keep a copy of panel and widgets to access later
        line = { "check"        : checkbox, 
                 "value_ppm"    : value_ppm, 
                 "value_area"   : value_area, 
                 "value_phase"  : value_phase,
               }

        # Add the controls to the grid sizer
        self.GridSizer.Add(line["check"], 0, wx.ALIGN_CENTER_VERTICAL)
        for key in ("value_ppm", "value_area", "value_phase"):
            self.GridSizer.Add(line[key], 0, wx.EXPAND, wx.ALIGN_CENTER)

        # Configure the controls I just created

        # Note. On these Spin and FloatSpin widgets, if the value you want to
        #    set is outside the wxGlade standard range, you should make the 
        #    call to reset the range first and then set the value you want.

        wx_util.configure_spin(value_ppm,  70, 2, 0.05,(minppm,maxppm))
        wx_util.configure_spin(value_area, 70, 3, 0.1, (0.001,100000.0))
        wx_util.configure_spin(value_phase,70, 1, 5.0, (-360,360))

        checkbox.SetValue(row_vals[0])
        value_ppm.SetValue(row_vals[1])
        value_area.SetValue(row_vals[2])
        value_phase.SetValue(row_vals[3])

        self.list_lines.append(line)

        # only need to update if a metabolite is added/removed from basis
        self.PanelParent.Bind(wx.EVT_CHECKBOX, self.event_handler, checkbox)
        self.PanelParent.Bind(EVT_FLOATSPIN,   self.event_handler, value_ppm)
        self.PanelParent.Bind(EVT_FLOATSPIN,   self.event_handler, value_area)
        self.PanelParent.Bind(EVT_FLOATSPIN,   self.event_handler, value_phase)

        self.PanelTop.Layout()  

        if update:
            self.event_handler()
        
        
    def remove_checked_rows(self):

        checklist = []
        for i, line in enumerate(self.list_lines):
            # gather indices of all checked boxes
            if line["check"].GetValue():
                checklist.append(i)
        
        # remove in reverse order so we don't invalidate later
        # indices by removing the ones preceding them in the list
        checklist.reverse()
        
        for i in checklist:
            # Each line is a dict of controls + mixture info
            for item in list(self.list_lines[i].values()):
                if hasattr(item, "Destroy"):
                    item.Destroy()              # It's a wx control
                
            del self.list_lines[i]
            
        # Reduce the # of rows in the grid sizer
        rows = self.GridSizer.GetRows()
        self.GridSizer.SetRows(rows - len(checklist))
        self.GridSizer.Layout()
        self.PanelTop.Layout()        
        self.event_handler()


    def select_all(self):
        for line in self.list_lines:
            line["check"].SetValue(True)


    def deselect_all(self):
        for line in self.list_lines:
            line["check"].SetValue(False)
            
            
    def event_handler(self, event=None):
        self.external_event_handler(event)


    def get_line_values(self, line):
        # Returns a dict containing  values of the controls in the line.
        return { "check"        : line["check"].GetValue(),
                 "value_ppm"    : line["value_ppm"].GetValue(),
                 "value_area"   : line["value_area"].GetValue(),
                 "value_phase"  : line["value_phase"].GetValue()
               }



## Test Code ######################################################################

class MyForm(wx.Frame):
 
    def __init__(self, dataset, metabolite):

        self.dataset = dataset
        self.metabolite = metabolite

        wx.Frame.__init__(self, None, wx.ID_ANY, "Testing the Dialog")

        self.statusbar = self.CreateStatusBar(4, 0)
        self.statusbar.SetStatusText("Ready")

        panel = wx.Panel(self, wx.ID_ANY)
        buttn = wx.Button(panel, label="Show Dialog")
        buttn.Bind(wx.EVT_BUTTON, self.on_dialog)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(buttn, 0, wx.ALL|wx.CENTER, 5)
        panel.SetSizer(sizer)

    def on_dialog(self, event):
        """
        Show the Dialog 
        """
        dialog = DialogUserDefinedMetabolite(self, self.dataset, self.metabolite, unique_names=['bob','bob3'])
        code = dialog.ShowModal()
        
        if code == wx.ID_OK:
            print("UserDefinedMetabolite  User hit OK")
        elif code == wx.ID_CANCEL:
            print("UserDefinedMetabolite  User Cancelled")

        print(dialog.metabolite.name)
        print(dialog.metabolite.spins)
        print(dialog.metabolite.ppms)
        print(dialog.metabolite.areas)
        print(dialog.metabolite.phases)
        
        dialog.Destroy()
 

 
#----------------------------------------------------------------------
# Run the program

if __name__ == "__main__":
    
    dataset = mrs_dataset.Dataset()
    metabolite = mrs_prior_metabolite.PriorMetabolite()
    metabolite.dims[0] = 'default'
    metabolite.spins  = 2
    metabolite.ppms   = [4.7,2.0,]
    metabolite.areas  = [3.0,1.0,]
    metabolite.phases = [0.0,0.0,]
    metabolite = None
    
    app = wx.App(False)
    frame = MyForm(dataset, metabolite)
    frame.Show()
    app.MainLoop()
