# Python modules

import os
import struct
from datetime import datetime


# 3rd party modules
import wx
import wx.grid as gridlib
import numpy as np

# Our modules
import vespa.datasim.prefs as prefs
import vespa.datasim.util_menu as util_menu
import vespa.datasim.util_datasim as util_datasim
import vespa.datasim.util_datasim_config as util_datasim_config
import vespa.datasim.util_import as util_import
import vespa.datasim.plot_panel_datasim as plot_panel_datasim
import vespa.datasim.auto_gui.datasim as datasim_ui
import vespa.datasim.mrs_datasim as mrs_datasim

from vespa.datasim.plot_panel_datasim import PlotPanelDatasim
from vespa.datasim.dynamic_metabolite_list import DynamicMetaboliteList
from vespa.datasim.dialog_datasim_resolution import DialogDatasimResolution

import vespa.common.util.ppm as util_ppm
import vespa.common.util.xml_ as util_xml
import vespa.common.util.time_ as util_time
import vespa.common.util.math_ as util_math
import vespa.common.util.export as util_export
import vespa.common.util.fileio as util_fileio
import vespa.common.wx_gravy.util as wx_util
import vespa.common.wx_gravy.common_dialogs as common_dialogs
import vespa.common.mrs_data_raw as mrs_data_raw
import vespa.common.configobj as configobj

from vespa.common.constants import DEGREES_TO_RADIANS
from vespa.datasim.util_datasim import calc_lw


#------------------------------------------------------------------------------

def _configure_combo(control, choices, selection=''):
    lines = list(choices.values())
    control.SetItems(lines)
    if selection in lines:
        control.SetStringSelection(selection)
    else:
        control.SetStringSelection(lines[0])


def _chop(data):
    return data * ((((np.arange(len(data)) + 1) % 2) * 2) - 1)


    



class TabDatasim(datasim_ui.DatasimUI):
    
    def __init__(self, outer_notebook, top, datasim=None):

        datasim_ui.DatasimUI.__init__(self, outer_notebook)
        
        self.top       = top
        self.parent    = outer_notebook
        self.datasim   = datasim

        self._prefs = prefs.PrefsMain()
        
        self.dataymax               = 1.0          # used for zoom out
        self.vertical_scale         = 1.0
        
        # values used in plot and export routines, filled in process()
        self.display_noise_in_plot  = True      
        self.metabolites_all        = None      
        self.metabolites_sum        = None
        self.metabolites_time_all   = None
        self.metabolites_time_sum   = None
        self.macromolecule_all      = None
        self.macromolecule_sum      = None
        self.macromolecule_time_all = None
        self.macromolecule_time_sum = None
        self.baseline_all           = None
        self.baseline_sum           = None
        self.baseline_time_all      = None
        self.baseline_time_sum      = None
        self.noise_freq             = None
        self.noise_time             = None
        self.last_export_filename   = ''
        
        self.plotting_enabled = False
        self.initialize_controls()
        self.populate_controls()
        self.plotting_enabled = True
        
        self.process_and_plot(set_scale=True)

        self.Bind(wx.EVT_WINDOW_DESTROY, self.on_destroy, self)
        self.Bind(wx.EVT_SPLITTER_SASH_POS_CHANGED, self.on_splitter, self.SplitterWindow)

        # Under OS X, wx sets the sash position to 10 (why 10?) *after*
        # this method is done. So setting the sash position here does no
        # good. We use wx.CallAfter() to (a) set the sash position and
        # (b) fake an EVT_SPLITTER_SASH_POS_CHANGED.
        wx.CallAfter(self.SplitterWindow.SetSashPosition,  self._prefs.sash_position, True)


    @property
    def view_mode(self):
        return (3 if self._prefs.plot_view_all else 1)


    ##### GUI Setup Handlers ##################################################

    def initialize_controls(self):
        """ 
        Set up sizes and bounds for widgets, but not values.
        - populate_controls() sets the values from data object.

        """
        ds = self.datasim
        dim0, dim1, dim2, dim3 = ds.dims
        sw      = ds.sw
        hpp     = sw / dim0
        resppm  = ds.resppm
        freq    = ds.frequency
        ppm_min = ds.pts2ppm(dim0-1)
        ppm_max = ds.pts2ppm(0)
        
        # Configure size, # of digits displayed, increment and min/max.

        wx_util.configure_spin(self.FloatScale, 70, 4, None, (0.0001, 1000))
        self.FloatScale.multiplier = 1.1
        
        wx_util.configure_spin(self.SpinIndex1, 60, None, None, (1,ds.loop_dims[1]))
        wx_util.configure_spin(self.SpinIndex2, 60, None, None, (1,ds.loop_dims[2]))
        wx_util.configure_spin(self.SpinIndex3, 60, None, None, (1,ds.loop_dims[3]))

        wx_util.configure_spin(self.FloatTa,           70, 3, 0.01, (0.001, 1000000.0))
        wx_util.configure_spin(self.FloatTb,           70, 3, 0.01, (0.001, 1000000.0))
        wx_util.configure_spin(self.FloatPhase0,       70, 2, 5, (-360,360))
        wx_util.configure_spin(self.FloatPhase1,       70, 2, 100, (-10000,10000))
        wx_util.configure_spin(self.FloatPhase1Pivot,  70, 2, 0.5, (ppm_min,ppm_max))
        wx_util.configure_spin(self.FloatB0Shift,      70, 2, 5, (-10000,10000))
        wx_util.configure_spin(self.SpinLeftShift,     70, None, None, (0,dim0-1))
        wx_util.configure_spin(self.FloatRefPeakArea,  70, 3, 0.5, (0.0,10000.0))
        wx_util.configure_spin(self.FloatRefPeakTa,    70, 3, 0.01, (0,1000.0))
        wx_util.configure_spin(self.FloatRefPeakTb,    70, 3, 0.01, (0,1000.0))
        wx_util.configure_spin(self.FloatNoiseRmsMultiplier, 70, 2, 5, (0.0000001,1000000.0))
        wx_util.configure_spin(self.SpinMonteCarloVoxels, 70, None, None, (1,10000))

        wx_util.configure_spin(self.FloatMmolGroupScale, 90, 10, 0.1, (0.0000000001,1000000.0))
        lines = ['None', 'Birch_MRM_2017', 'Wilson_MRM_2021']
        self.ChoicePresetMmol.SetItems(lines)
        self.ChoicePresetMmol.SetStringSelection(lines[0])

        wx_util.configure_spin(self.FloatBaseGroupScale, 90, 10, 0.1, (0.0000000001,1000000.0))


    def populate_controls(self):
        """ 
        Set values from the data object into widgets.
        - no range/value checking done in here
        
        """
        ds   = self.datasim
        loop = ds.loop
        hpp  = ds.hpp

        #----------------------------------------------------------------------
        # Dataset View setup 

        self.view = PlotPanelDatasim(self.PanelDatasimPlot,
                                     self,
                                     self.parent,
                                     naxes=3,
                                     reversex=True,
                                     zoom='span',
                                     reference=True,
                                     middle=True,
                                     do_zoom_select_event=True,
                                     do_zoom_motion_event=True,
                                     do_refs_select_event=True,
                                     do_refs_motion_event=True,
                                     do_middle_select_event=True,
                                     do_middle_motion_event=True,
                                     do_scroll_event=True,
                                     props_zoom=dict(alpha=0.2, facecolor='yellow'),
                                     props_cursor=dict(alpha=0.3, facecolor='gray'),
                                     xscale_bump=0.0,
                                     yscale_bump=0.05,
                                     data=[],
                                     prefs=self._prefs,
                                     dataset=ds)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.view, 1, wx.LEFT | wx.TOP | wx.EXPAND)
        self.PanelDatasimPlot.SetSizer(sizer)
        self.view.Fit() 

        #----------------------------------------------------------------------
        # Global controls            

        self.TextData.SetValue(ds.experiment.name)
        self.FloatScale.SetValue(self.vertical_scale)
        self.SpinIndex1.SetValue(loop[0]+1)
        self.SpinIndex2.SetValue(loop[1]+1)
        self.SpinIndex3.SetValue(loop[2]+1)

        #----------------------------------------------------------------------
        # Spectral settings controls            

        self.FloatTa.SetValue(ds.ta)
        self.FloatTb.SetValue(ds.tb)
        self.FloatPhase0.SetValue(ds.phase0)
        self.FloatPhase1.SetValue(ds.phase1)
        self.FloatPhase1Pivot.SetValue(ds.phase_1_pivot)
        self.FloatB0Shift.SetValue(ds.b0shift)
        self.FloatRefPeakArea.SetValue(ds.noise_ref_peak_area)
        self.FloatRefPeakTa.SetValue(ds.noise_ref_peak_ta)
        self.FloatRefPeakTb.SetValue(ds.noise_ref_peak_tb)
        self.FloatNoiseRmsMultiplier.SetValue(ds.noise_rms_multiplier)
        self.SpinMonteCarloVoxels.SetValue(ds.montecarlo_voxels)

        self.CheckDisplayNoiseInPlot.SetValue(self.display_noise_in_plot)
        self.TextLinewidth.SetLabel("%.3f" % (calc_lw(ds.ta,ds.tb),))
        self.TextEffectiveLinewidth.SetLabel("%.3f" % (calc_lw(ds.noise_ref_peak_ta,ds.noise_ref_peak_tb),))
        self.TextComment.SetValue(ds.comment)

        #----------------------------------------------------------------------
        # Metabolite signal controls            
        # - set up Metabolite dynamic list

        # The list grid sizer is marked so we can find it at run-time
        self.MetaboliteGridSizer = self.LabelMetabolites.GetContainingSizer()
        grid_parent = self.LabelMetabolites.GetParent()
        self.LabelMetabolites.Destroy()
        
        # By setting the row count to 0, we signal to wx that we will let it 
        # decide how many rows the grid sizer needs.
        self.MetaboliteGridSizer.Clear()
        self.MetaboliteGridSizer.SetRows(0)

        # Add headings to the first row of the grid sizer.
        headings = ( "Metabolites", "Area Scale\nFactor", "T2 (Ta) Decay\n[sec]")
        
        for heading in headings:
            if heading:
                label = wx.StaticText(grid_parent, label=heading, style=wx.ALIGN_CENTRE)
                self.MetaboliteGridSizer.Add(label, 0, wx.EXPAND|wx.ALIGN_CENTER_VERTICAL|wx.TOP, 8)
            else:
                self.MetaboliteGridSizer.AddSpacer(1)
                
        self.metab_list = DynamicMetaboliteList(self.PanelMetaboliteLines, self,
                                                self.MetaboliteGridSizer, ds)

        # ----------------------------------------------------------------------
        # Macromolecule Grid
        self.grid_set_values(grid=self.GridMmol)
        self.GridMmol.Bind(gridlib.EVT_GRID_LABEL_LEFT_CLICK, self.on_label_left_clicked)
        val = 0 if ds.mmol_lineshape=='lorentzian' else 1
        self.RadioMmolLineshape.SetSelection(val)
        self.FloatMmolGroupScale.SetValue(ds.mmol_group_scale)

        #----------------------------------------------------------------------
        # Baseline Grid
        self.grid_set_values(grid=self.GridBase)
        self.GridBase.Bind(gridlib.EVT_GRID_LABEL_LEFT_CLICK, self.on_label_left_clicked)
        val = 0 if ds.base_lineshape=='lorentzian' else 1
        self.RadioBaseLineshape.SetSelection(val)
        self.FloatBaseGroupScale.SetValue(ds.base_group_scale)


    ##### Menu & Notebook Event Handlers ######################################################

    def on_label_left_clicked(self, event):
        ds = self.datasim
        obj, col = event.GetEventObject(), event.GetCol()
        ascend = True
        if col == obj.GetSortingColumn():
            if obj.IsSortOrderAscending():
                ascend = False
        obj.SetSortingColumn(col, ascending=ascend)

        r = np.array(self.grid_get_values(grid=obj))
        if ascend:
            r = r[:,r[col].argsort()]
        else:
            r = r[:,r[col].argsort()[::-1]]

        if obj == self.GridMmol:
            ds.mmol_flags  = [item==1.0 for item in r[0]]
            ds.mmol_ppms   = [item for item in r[1]]
            ds.mmol_areas  = [item for item in r[2]]
            ds.mmol_phases = [item for item in r[3]]
            ds.mmol_widths = [item for item in r[4]]
        elif obj == self.GridBase:
            ds.base_flags  = [item==1.0 for item in r[0]]
            ds.base_ppms   = [item for item in r[1]]
            ds.base_areas  = [item for item in r[2]]
            ds.base_phases = [item for item in r[3]]
            ds.base_widths = [item for item in r[4]]

        self.grid_set_values(grid=obj)

        #print('on_label_left_clicked')


    def on_activation(self):
        # This is a faux event handler. wx doesn't call it directly. It's 
        # a notification from my parent (the dataset notebook) to let
        # me know that this tab has become the current one.
        
        # Force the View menu to match the current plot options.
        util_menu.bar.set_menu_from_state(self._prefs.menu_state)
       

    def on_destroy(self, event):
        self._prefs.save()


    def on_menu_view_option(self, event):
        event_id = event.GetId()

        if self._prefs.handle_event(event_id):
            if event_id in (util_menu.ViewIds.ZERO_LINE_SHOW,
                            util_menu.ViewIds.ZERO_LINE_TOP,
                            util_menu.ViewIds.ZERO_LINE_MIDDLE,
                            util_menu.ViewIds.ZERO_LINE_BOTTOM,
                            util_menu.ViewIds.XAXIS_SHOW,
                           ):
                self.view.update_axes()
                self.view.canvas.draw()

            if event_id in (util_menu.ViewIds.DATA_TYPE_REAL,
                            util_menu.ViewIds.DATA_TYPE_IMAGINARY,
                            util_menu.ViewIds.DATA_TYPE_MAGNITUDE,
                            util_menu.ViewIds.DATA_TYPE_SUMMED,
                            util_menu.ViewIds.XAXIS_PPM,
                            util_menu.ViewIds.XAXIS_HERTZ,
                           ):
                if event_id == util_menu.ViewIds.DATA_TYPE_REAL:
                    self.view.set_data_type_real()
                elif event_id == util_menu.ViewIds.DATA_TYPE_IMAGINARY:
                    self.view.set_data_type_imaginary()
                elif event_id == util_menu.ViewIds.DATA_TYPE_MAGNITUDE:
                    self.view.set_data_type_magnitude()
                elif event_id == util_menu.ViewIds.DATA_TYPE_SUMMED:
                    self.view.set_data_type_summed(index=[1,2])

                self.view.update(no_draw=True)
                self.view.set_phase_0(0.0, no_draw=True)
                self.view.canvas.draw()

            if event_id in (util_menu.ViewIds.PLOT_VIEW_FINAL,
                            util_menu.ViewIds.PLOT_VIEW_ALL,
                           ):
                if event_id == util_menu.ViewIds.PLOT_VIEW_FINAL:
                    self.view.change_naxes(1)
                else:
                    self.view.change_naxes(3)
            

    def on_menu_view_output(self, event):

        event_id = event.GetId()

        formats = { util_menu.ViewIds.VIEW_TO_PNG : "PNG",
                    util_menu.ViewIds.VIEW_TO_SVG : "SVG", 
                    util_menu.ViewIds.VIEW_TO_EPS : "EPS", 
                    util_menu.ViewIds.VIEW_TO_PDF : "PDF", 
                  }

        if event_id in formats:
            format   = formats[event_id]
            lformat  = format.lower()
            filter_  = "%s files (*.%s)|*.%s" % (format, lformat, lformat)
            figure   = self.view.figure
            filename = common_dialogs.save_as("", filter_)

            if filename:
                msg = ""
                try:
                    figure.savefig( filename,
                                    dpi=300, 
                                    facecolor='w', 
                                    edgecolor='w',
                                    orientation='portrait', 
                                    papertype='letter', 
                                    format=None,
                                    transparent=False)
                except IOError:
                    msg = """I can't write the file "%s".""" % filename
                
                if msg:
                    common_dialogs.message(msg, style=common_dialogs.E_OK)


        if event_id in (util_menu.ViewIds.EXPERIMENT_TO_TEXT,):
            experiment = self.datasim.experiment
            lines = str(experiment)
            lines += "\n\nVespa Simulation Experment Results\n" + "-" * 75 + "\n\n"
            lines += "\n".join([simulation.summary() for simulation in experiment.simulations])
            wx_util.display_text_as_file(lines)


    ########## Widget Event Handlers ####################    
    
    def on_scale(self, event):
        view = self.view
        scale = self.FloatScale.GetValue()
        if scale > view.vertical_scale:
            view.set_vertical_scale(1.0, scale_mult=1.1)
        else:
            view.set_vertical_scale(-1.0, scale_mult=1.1)
        self.FloatScale.SetValue(view.vertical_scale)

    def on_index1(self, event):
        self.datasim.loop[0] = self.SpinIndex1.GetValue()-1
        self.process_and_plot()

    def on_index2(self, event):
        self.datasim.loop[1] = self.SpinIndex2.GetValue()-1
        self.process_and_plot()

    def on_index3(self, event):
        self.datasim.loop[2] = self.SpinIndex3.GetValue()-1
        self.process_and_plot()

    def on_set_spectral_resolution(self, event):
        ds = self.datasim
        pts       = ds.dims[0]
        sw        = ds.sw
        frequency = ds.frequency
        ppm_start = ds.mets_ppm_start
        ppm_end   = ds.mets_ppm_end

        dialog = DialogDatasimResolution(self, sw, pts, frequency, ppm_start, ppm_end)
        if dialog.ShowModal() == wx.ID_OK:
            
            ds.sw             = dialog.sweep_width
            ds.dims[0]        = int(dialog.points)
            ds.frequency      = dialog.frequency
            ds.mets_ppm_start = dialog.ppm_range_start
            ds.mets_ppm_end   = dialog.ppm_range_end

            # reset mmol tab max/min ppm range and Width (ppm, hz, damp) values
            self.grid_check_ppms(grid=self.GridMmol)
            self.grid_check_ppms(grid=self.GridBase)

            wx.SetCursor(wx.HOURGLASS_CURSOR)
            ds.calculate_basis()
            self.mmol_update_basis()
            self.base_update_basis()
            self.process_and_plot()
            wx.SetCursor(wx.NullCursor)
            
        dialog.Destroy()


    def on_splitter(self, event):
        self._prefs.sash_position = self.SplitterWindow.GetSashPosition()


    def on_ta(self, event):
        self.datasim.ta = self.FloatTa.GetValue()
        self.TextLinewidth.SetLabel("%.3f" % (calc_lw(self.datasim.ta,self.datasim.tb),))
        self.process_and_plot()

    def on_tb(self, event):
        self.datasim.tb = self.FloatTb.GetValue()
        self.TextLinewidth.SetLabel("%.3f" % (calc_lw(self.datasim.ta,self.datasim.tb),))
        self.process_and_plot()

    def on_phase0(self, event):
        # Sets the 0th order phase for display
        val = self.FloatPhase0.GetValue()
        self.datasim.phase0 = val
        self.view.set_phase_0(val, absolute=True)
        self.view.canvas.draw()        

    def on_phase1(self, event):
        # Sets the 1th order phase for display
        val = self.FloatPhase1.GetValue()
        self.datasim.phase1 = val
        self.view.set_phase_1(val, absolute=True)
        self.view.canvas.draw()        

    def on_phase_1_pivot(self, event):
        # Sets the 1th order phase pivot location for display
        val = self.FloatPhase1Pivot.GetValue()
        self.datasim.phase_1_pivot = val
        # next line just refreshes the plot with new pivot
        self.view.set_phase_1(0.0, absolute=False)
        self.view.canvas.draw()        

    def on_b0_shift(self, event):
        # Sets the B0 shift in Hz for display
        self.datasim.b0shift = self.FloatB0Shift.GetValue()
        self.process_and_plot()

    def on_left_shift(self, event):
        # Sets the left shift in pts for display
        self.datasim.left_shift = self.SpinLeftShift.GetValue()
        self.process_and_plot()

    def on_display_noise_in_plot(self, event):    
        self.display_noise_in_plot = self.CheckDisplayNoiseInPlot.GetValue()
        self.process_and_plot()
        self.TextEffectiveSnr.SetLabel("%.3f" % (self.snr_actual),)

    def on_ref_peak_area(self, event):
        self.datasim.noise_ref_peak_area = self.FloatRefPeakArea.GetValue()
        self.process_and_plot()
        self.TextEffectiveSnr.SetLabel("%.3f" % (self.snr_actual),)

    def on_ref_peak_ta(self, event):
        self.datasim.noise_ref_peak_ta = self.FloatRefPeakTa.GetValue()
        self.process_and_plot()
        self.TextEffectiveSnr.SetLabel("%.3f" % (self.snr_actual),)
        val = calc_lw(self.datasim.noise_ref_peak_ta,self.datasim.noise_ref_peak_tb)
        self.TextEffectiveLinewidth.SetLabel("%.3f" % (val,))

    def on_ref_peak_tb(self, event):
        self.datasim.noise_ref_peak_tb = self.FloatRefPeakTb.GetValue()
        self.process_and_plot()
        self.TextEffectiveSnr.SetLabel("%.3f" % (self.snr_actual,))
        val = calc_lw(self.datasim.noise_ref_peak_ta,self.datasim.noise_ref_peak_tb)
        self.TextEffectiveLinewidth.SetLabel("%.3f" % (val,))
    
    def on_noise_rms_multiplier(self, event):
        self.datasim.noise_rms_multiplier = self.FloatNoiseRmsMultiplier.GetValue() 
        self.process_and_plot()
        self.TextEffectiveSnr.SetLabel("%.3f" % (self.snr_actual,))
        
    def on_matrix_size(self, event):
        self.datasim.montecarlo_voxels = self.SpinMonteCarloVoxels.GetValue()

    def on_comment(self, event):
        self.datasim.comment = event.GetEventObject().GetValue()


    # --- Metabolite events -----------------

    def on_dynamic_metabolite_list(self, event):
        # fudged event called from actual event inside the dynamic list class.
        self.TextLinewidth.SetLabel("%.3f" % (calc_lw(self.datasim.ta,self.datasim.tb),))
        self.process_and_plot()

    def on_select_all(self, event):
        self.metab_list.select_all()
        self.process_and_plot()

    def on_select_none(self, event):
        self.metab_list.deselect_all()
        self.process_and_plot()


    # --- Macromolecule events -----------------

    def on_mmol_lineshape(self, event):
        obj = event.GetEventObject()
        self.datasim.mmol_lineshape = 'gaussian'
        if obj.GetSelection() == 0:
            self.datasim.mmol_lineshape = 'lorentzian'
        self.mmol_update_basis()
        self.process_and_plot()

    def on_mmol_group_scale(self, event):
        self.datasim.mmol_group_scale = self.FloatMmolGroupScale.GetValue()
        self.process_and_plot()

    def on_mmol_select_all(self, event):
        for i in range(self.GridMmol.GetNumberRows()):
            self.GridMmol.SetCellValue(i, 0, '1')
        self.GridMmol.ForceRefresh()
        self.mmol_update_basis()
        self.process_and_plot()

    def on_mmol_select_none(self, event):
        for i in range(self.GridMmol.GetNumberRows()):
            self.GridMmol.SetCellValue(i, 0, '')
        self.GridMmol.ForceRefresh()
        self.mmol_update_basis()
        self.process_and_plot()

    def on_mmol_delete_selected(self, event):
        nrow = self.GridMmol.GetNumberRows()
        self.GridMmol.BeginBatch()
        remove = [i for i in range(nrow) if self.GridMmol.GetCellValue(i, 0) == '1']
        for i in remove[::-1]:
            self.GridMmol.DeleteRows(i, 1)
        self.GridMmol.EndBatch()
        self.Layout()
        self.Refresh()
        self.mmol_update_basis()
        self.process_and_plot()

    def on_mmol_add_line(self, event):
        pg = (self.datasim.mmol_lineshape=='gaussian')
        val_ppm  = self.datasim.hz2ppm(20.0, rel=True)
        val_damp = self.datasim.width_ppm2damp(val_ppm, pure_gauss=pg)
        new_line = ['', '4.7', '1.0', '0.0', str(val_ppm), '20.0', str(val_damp)]
        if self.GridMmol.AppendRows():
            nrow = self.GridMmol.GetNumberRows() - 1
            for i,val in enumerate(new_line):
                self.GridMmol.SetCellValue(nrow, i, val)
                self.GridMmol.SetCellAlignment(wx.ALIGN_CENTER, nrow, i)
        self.mmol_update_basis()
        self.process_and_plot()

    def on_mmol_cell_changed(self, event):
        ds = self.datasim
        lim = [ [0,500],
                [ds.pts2ppm(ds.dims[0]-1),ds.pts2ppm(0)],   # ppm
                [0.00000001,100000.0],                      # area scale
                [-360.0,360.0],                             # phase deg
                [0.00000001,1000000.0],                     # width ppm
                [0.00000001,1000000.0],                     # width hz
                [0.00000001,1000000.0]]                     # damp sec
        lshape = (ds.mmol_lineshape=='gaussian')
        self.grid_cell_changed(event, lshape, lim)

    def on_mmol_import_hlsvd_file(self, event):
        self.grid_import_hlsvd_file(grid=self.GridMmol)

    def on_preset_mmol(self, event):
        ds = self.datasim
        if event.String == 'None':
            ppms   = mrs_datasim.DEFAULT_MMOL_PPMS
            areas  = mrs_datasim.DEFAULT_MMOL_AREAS
            flags  = mrs_datasim.DEFAULT_MMOL_FLAGS
            phases = mrs_datasim.DEFAULT_MMOL_PHASES
            widths = mrs_datasim.DEFAULT_MMOL_WIDTHS
            lshape = 'lorentzian'
        elif event.String == 'Birch_MRM_2017':
            ppms   = [0.91,1.21,1.38,1.63,2.01,2.09,2.25,2.61,2.96,3.11,3.67,3.80,3.96]
            areas  = [0.72,0.28,0.38,0.05,0.45,0.36,0.36,0.04,0.20,0.11,0.64,0.07,1.00]
            flags  = [True] * len(areas)
            phases = [0.0] * len(areas)
            widths = [21.2,19.16,15.9,7.5,29.03,20.53,17.89,5.3,14.02,17.89,33.52,11.85,37.48] # Hz
            widths = [ds.hz2ppm(item, rel=True) for item in widths]
            lshape = 'gaussian'
        elif event.String == 'Wilson_MRM_2021':
            ppms   = [1.28,1.28,0.89,0.91,2.04,2.25,2.80,2.08,2.25,1.95,3.00,1.21,1.43,1.67]
            areas  = [2.0,2.0,3.0,3.0,1.33,0.67,0.87,1.33,0.33,0.33,0.4,2.0,2.0,2.0]
            flags  = [True] * len(areas)
            phases = [0.0] * len(areas)
            widths = [0.15,0.089,0.14,0.14,0.15,0.15,0.20,0.15,0.20,0.15,0.20,0.15,0.17,0.15] # Hz
            lshape = 'lorentzian'

        ds.mmol_flags = flags
        ds.mmol_ppms  = ppms
        ds.mmol_areas = areas
        ds.mmol_phases = phases
        ds.mmol_widths = widths
        ds.mmol_lineshape = lshape
        if lshape == 'lorentzian':
            self.RadioMmolLineshape.SetSelection(0)
        else:
            self.RadioMmolLineshape.SetSelection(1)
        self.grid_set_values(grid=self.GridMmol)
        ds.macromolecule_basis = ds.calculate_signals(ppms, areas, phases, widths, lshape)
        self.process_and_plot()

    # --- Baseline events -----------------

    def on_base_lineshape(self, event):
        obj = event.GetEventObject()
        self.datasim.base_lineshape = 'gaussian'
        if obj.GetSelection() == 0:
            self.datasim.base_lineshape = 'lorentzian'
        self.base_update_basis()
        self.process_and_plot()

    def on_base_group_scale(self, event):
        self.datasim.base_group_scale = event.GetEventObject().GetValue()
        self.process_and_plot()

    def on_base_select_all(self, event):
        for i in range(self.GridBase.GetNumberRows()):
            self.GridBase.SetCellValue(i, 0, '1')
        self.GridBase.ForceRefresh()
        self.base_update_basis()
        self.process_and_plot()

    def on_base_select_none(self, event):
        for i in range(self.GridBase.GetNumberRows()):
            self.GridBase.SetCellValue(i, 0, '')
        self.GridBase.ForceRefresh()
        self.base_update_basis()
        self.process_and_plot()

    def on_base_delete_selected(self, event):
        nrow = self.GridBase.GetNumberRows()
        self.GridBase.BeginBatch()
        remove = [i for i in range(nrow) if self.GridBase.GetCellValue(i, 0) == '1']
        for i in remove[::-1]:
            self.GridBase.DeleteRows(i, 1)
        self.GridBase.EndBatch()
        self.Layout()
        self.Refresh()
        self.base_update_basis()
        self.process_and_plot()

    def on_base_add_line(self, event):
        pg = (self.datasim.base_lineshape=='gaussian')
        val_ppm  = self.datasim.hz2ppm(45.0, rel=True)
        val_damp = self.datasim.width_ppm2damp(val_ppm, pure_gauss=pg)
        new_line = ['', '4.7', '1.0', '0.0', str(val_ppm), '45', str(val_damp)]
        if self.GridBase.AppendRows():
            nrow = self.GridBase.GetNumberRows() - 1
            for i,val in enumerate(new_line):
                self.GridBase.SetCellValue(nrow, i, val)
                self.GridBase.SetCellAlignment(wx.ALIGN_CENTER, nrow, i)
        self.base_update_basis()
        self.process_and_plot()

    def on_base_cell_changed(self, event):
        ds = self.datasim
        lim = [ [0,500],
                [ds.pts2ppm(ds.dims[0]-1),ds.pts2ppm(0)],   # ppm
                [0.00000001,100000.0],                      # area scale
                [-360.0,360.0],                             # phase deg
                [0.00000001,1000000.0],                     # width ppm
                [0.00000001,1000000.0],                     # width hz
                [0.00000001,1000000.0]]                     # damp sec
        lshape = (ds.base_lineshape=='gaussian')
        self.grid_cell_changed(event, lshape, lim)


    def on_base_import_hlsvd_file(self, event):
        self.grid_import_hlsvd_file(grid=self.GridBase)



    ##### Internal helper functions  ##########################################
    
    def grid_check_ppms(self, grid=None):
        """
        PPM range can change if spectral resolution set to new values. This
        method checks current values in Mmol/Base Grid widget and sets to edge
        values if outside new range. Background set to red if changed.

        """
        if grid is None:
            grid = self.GridMmol
        ds = self.datasim
        ppm_max, ppm_min = ds.pts2ppm(0), ds.pts2ppm(ds.dims[0]-1)

        for row in range(grid.GetNumberRows()):
            val = grid.GetCellValue(row,1)
            if float(val) < ppm_min:
                grid.SetCellValue(row, 1, str(ppm_min))
                grid.SetCellTextColour(row, 1, wx.RED)
            elif float(val) > ppm_max:
                grid.SetCellValue(row, 1, str(ppm_max))
                grid.SetCellTextColour(row, 1, wx.RED)
            else:
                grid.SetCellTextColour(row, 1, wx.BLACK)


    def grid_get_values(self, grid=None):
        """ return current values in the grid object provided """
        if grid is None:
            grid = self.GridMmol
        if grid==self.GridMmol:
            pg = (self.datasim.mmol_lineshape == 'gaussian')
        else:
            pg = (self.datasim.base_lineshape == 'gaussian')
        nrow = grid.GetNumberRows()
        self.grid_check_ppms(grid=grid)
        flags  = [grid.GetCellValue(i,0)=='1'   for i in range(nrow)]
        ppms   = [float(grid.GetCellValue(i,1)) for i in range(nrow)]
        areas  = [float(grid.GetCellValue(i,2)) for i in range(nrow)]
        phases = [float(grid.GetCellValue(i,3)) for i in range(nrow)]
        widths = [float(grid.GetCellValue(i,4)) for i in range(nrow)]
        widhz  = [self.datasim.ppm2hz(item,rel=True) for item in widths]
        damps  = [self.datasim.width_ppm2damp(item, pure_gauss=pg) for item in widths]

        return flags, ppms, areas, phases, widths, widhz, damps


    def grid_set_values(self, grid=None):
        """
        Sets current Mmol or Base values in Datasim object into the grid object
        provided. Generally used when Datasim object has been changed/reset.

        """
        if grid is self.GridMmol:
            lines = self.datasim.mmol_to_lines
        elif grid is self.GridBase:
            lines = self.datasim.base_to_lines

        ncol, nrow = grid.GetNumberCols(), grid.GetNumberRows()
        if ncol > 0: grid.DeleteCols(pos=0, numCols=ncol)
        if nrow > 0: grid.DeleteRows(pos=0, numRows=nrow)
        grid.ClearGrid()
        grid.AppendCols(7)
        grid.AppendRows(len(lines))

        items = ("Peaks", "  PPM  ", " Area Factor ", "Phase [deg]", " Width [ppm] ", " Width [hz] ", " Damp [sec] ")
        for col, labl in enumerate(items):
            grid.SetColLabelValue(col, labl)

        grid.SetColFormatBool(0)
        for i,val in [(1,3),(2,8),(3,2),(4,8),(5,8),(6,8)]:
            grid.SetColFormatFloat(i, precision=val)

        grid.BeginBatch()
        for row, line in enumerate(lines):
            grid.SetCellValue( row, 0, '')
            grid.SetCellEditor(row, 0, gridlib.GridCellBoolEditor())
            grid.SetCellRenderer(row, 0, gridlib.GridCellBoolRenderer())
            for col, item in enumerate(line):
                if col == 0:
                    val = '1' if item else ''
                    grid.SetCellValue(row, 0, val)
                    grid.SetCellEditor(row, 0, gridlib.GridCellBoolEditor())
                    grid.SetCellRenderer(row, 0, gridlib.GridCellBoolRenderer())
                elif col in [1,2,3,4,5]:
                    grid.SetCellValue(row, col, str(float(item)))
                    grid.SetCellEditor(row, col, gridlib.GridCellFloatEditor())
                elif col == 6:
                    grid.SetCellValue(row, col, str(float(item)))
                    grid.SetCellEditor(row, col, gridlib.GridCellFloatEditor())
                    grid.SetReadOnly(row, col, isReadOnly=True)
                # bjs grid.SetCellAlignment(wx.ALIGN_CENTER, row, col)
                grid.SetCellAlignment(row, col, wx.ALIGN_CENTER, wx.ALIGN_CENTER)
        grid.EndBatch()

        grid.SetColLabelAlignment(wx.ALIGN_CENTER, wx.ALIGN_CENTER)
        grid.HideRowLabels()
        grid.SetMargins(0, 0)
        grid.AutoSizeColumns(True)

        if grid is self.GridMmol:
            grid.Bind(gridlib.EVT_GRID_CELL_CHANGED, self.on_mmol_cell_changed)
        elif grid is self.GridBase:
            grid.Bind(gridlib.EVT_GRID_CELL_CHANGED, self.on_base_cell_changed)
        self.grid_check_ppms(grid=grid)


    def grid_cell_changed(self, event, lshape, lim):
        """
        Check if cell value is between limits, if outside, set to min/max value
        and set font color to red. Cols 4/5/6 are inter-related, need to update
        depending on which changed

        """
        ds  = self.datasim
        obj = event.GetEventObject()
        row, col = event.GetRow(), event.GetCol()
        value = obj.GetCellValue(row,col)
        pg = (ds.mmol_lineshape=='gaussian') if obj==self.GridMmol else (ds.base_lineshape=='gaussian')
        if col > 0:
            value = float(value)
            if value < lim[col][0]:
                value = lim[col][0]
                obj.SetCellTextColour(row, col, wx.RED)
            elif value > lim[col][1]:
                value = lim[col][1]
                obj.SetCellTextColour(row, col, wx.RED)
            else:
                self.GridMmol.SetCellTextColour(row, col, wx.BLACK)
            obj.SetCellValue(row, col, str(value))

        if col==4:          # width - ppm
            obj.SetCellValue(row, 5, str(ds.width_ppm2hz(value)))
            obj.SetCellValue(row, 6, str(ds.width_ppm2damp(value, pure_gauss=pg)))
        elif col==5:        # width - hz
            val_ppm = ds.width_hz2ppm(value)
            obj.SetCellValue(row, 4, str(val_ppm))
            obj.SetCellValue(row, 6, str(ds.width_ppm2damp(val_ppm, pure_gauss=pg)))
        elif col==6:        # damp - sec
            val_ppm = ds.width_damp2ppm(value)
            obj.SetCellValue(row, 4, str(val_ppm))
            obj.SetCellValue(row, 5, str(ds.width_ppm2hz(val_ppm)))

        if obj==self.GridMmol:
            self.mmol_update_basis()
        else:
            self.base_update_basis()
        self.process_and_plot()


    def grid_import_hlsvd_file(self, grid):
        ds = self.datasim
        label = "Select HLSVD file with 'analysis_export_hlsvd' Node"
        fname = common_dialogs.pickfile(label, "XML files (*.xml)|*.xml")
        if fname:
            msg = ""
            try:
                importer = util_import.AnalysisHlsvdImporter(fname)
            except IOError:
                msg = """I can't read the file "%s".""" % fname
            except SyntaxError:
                msg = """The file "%s" isn't valid Vespa Interchange File Format.""" % fname

            if msg:
                common_dialogs.message(msg, "TabDatasim - Import Analysis Hlsvd Results File", common_dialogs.E_OK)
            else:
                # Time to rock and roll!
                wx.BeginBusyCursor()
                objs = importer.go()
                wx.EndBusyCursor()

                if objs == []:
                    msg = """No HLSVD results found in: "%s".""" % fname
                    common_dialogs.message(msg, "TabDatasim - Import Analysis Hlsvd Results File", common_dialogs.E_OK)
                    return
                else:
                    source = objs[0]

                item = source.find("ppms")      # If ppms present, so are others.
                if item is not None:

                    ppms   = util_xml.element_to_numpy_array(item)
                    areas  = util_xml.element_to_numpy_array(source.find("ampls"))
                    phases = util_xml.element_to_numpy_array(source.find("phases"))
                    flags  = util_xml.element_to_numpy_array(source.find("flags"))
                    widths = util_xml.element_to_numpy_array(source.find("damps_ppm"))
                    lshape = 'lorentzian'

                    if grid==self.GridMmol:
                        ds.mmol_ppms   = ppms
                        ds.mmol_areas  = areas
                        ds.mmol_phases = phases
                        ds.mmol_flags  = flags
                        ds.mmol_widths = widths
                        ds.mmol_lineshape = lshape
                        self.RadioMmolLineshape.SetSelection(0)
                        self.grid_set_values(grid=grid)
                        ds.macromolecule_basis = ds.calculate_signals(ppms, areas, phases, widths, lshape)
                    elif grid==self.GridBase:
                        ds.base_ppms   = ppms
                        ds.base_areas  = areas
                        ds.base_phases = phases
                        ds.base_flags  = flags
                        ds.base_widths = widths
                        ds.base_lineshape = lshape
                        self.RadioBaseLineshape.SetSelection(0)
                        self.grid_set_values(grid=grid)
                        ds.baseline_basis = ds.calculate_signals(ppms, areas, phases, widths, lshape)

                    self.process_and_plot()


    def mmol_update_basis(self):
        """ update datasim Mmol values from Grid and recreate basis array """
        ds = self.datasim
        r = self.grid_get_values(grid=self.GridMmol)
        flags, ppms, areas, phases, widths, widths_hz, widths_damp = r
        ds.mmol_flags = flags
        ds.mmol_ppms = ppms
        ds.mmol_areas = areas
        ds.mmol_phases = phases
        ds.mmol_widths = widths
        lshape = ds.mmol_lineshape
        r = ds.calculate_signals(ppms, areas, phases, widths, lshape)
        ds.macromolecule_basis = r


    def base_update_basis(self):
        """ update datasim Base values from Grid and recreate basis array """
        ds = self.datasim
        r = self.grid_get_values(grid=self.GridBase)
        flags, ppms, areas, phases, widths, widths_hz, widths_damp = r
        ds.base_flags = flags
        ds.base_ppms = ppms
        ds.base_areas = areas
        ds.base_phases = phases
        ds.base_widths = widths
        lshape = ds.base_lineshape
        r = ds.calculate_signals(ppms, areas, phases, widths, lshape)
        ds.baseline_basis = r


    def rms(self, data):
        return np.sqrt(np.sum( (data - np.mean(data))**2 )/(len(data)-1.0) )


    def export_spectrum_to_viff(self, monte_carlo_flag=False, comment=''):

        if comment == '':
            comment = 'Exported from Vespa-Datasim using Export Spectrum to VIFF'

        default_path = util_datasim_config.get_last_export_path()
        
        filename = common_dialogs.save_as("Save As XML/VIFF (Vespa Interchange File Format)",
                                          "VIFF/XML files (*.xml)|*.xml",
                                          default_path, self.last_export_filename)
        if filename:
            ds = self.datasim
            sw = ds.sw
            if not monte_carlo_flag:
                dims0, dims1 = ds.dims[0], 1
            else:
                dims0, dims1 = ds.dims[0], ds.montecarlo_voxels

            val   = 1j * ds.phase0 * DEGREES_TO_RADIANS
            phase = np.exp(val)
            
            if not monte_carlo_flag:
                data  = self.metabolites_time_sum.copy() * phase
                base  = self.baseline_time_sum.copy() * phase
                mmol  = self.macromolecule_time_sum.copy() * phase
                noise = self.noise_time.copy() * phase
                final = (self.metabolites_time_sum + self.macromolecule_time_sum + self.baseline_time_sum + self.noise_time).copy() * phase
                snra  = self.snr_actual
                results = [data, mmol, base, noise, final]
                for item in results:
                    item.shape = 1,1,1,item.shape[-1]
            else:
                data  = np.zeros((1, 1, dims1, dims0), complex)
                mmol  = np.zeros((1, 1, dims1, dims0), complex)
                base  = np.zeros((1, 1, dims1, dims0), complex)
                final = np.zeros((1, 1, dims1, dims0), complex)
                noise = np.zeros((1, 1, dims1, dims0), complex)
                snra  = np.zeros((dims1,), float)            # actual measured SNR for each voxel

                save_flag = self.display_noise_in_plot
                self.display_noise_in_plot = True

                for i in range(dims1):
                    self.process()
                    data[0,0,i,:]  = self.metabolites_time_sum.copy() * phase
                    mmol[0,0,i,:]  = self.macromolecule_time_sum.copy()* phase
                    base[0,0,i,:]  = self.baseline_time_sum.copy() * phase
                    noise[0,0,i,:] = self.noise_time.copy() * phase
                    final[0,0,i,:] = (self.metabolites_time_sum + self.macromolecule_time_sum + self.baseline_time_sum + self.noise_time).copy() * phase
                    snra[i]        = self.snr_actual

                results = [data, mmol, base, noise, final]
                self.display_noise_in_plot = save_flag

            #--------------------------------------------------------
            # get all spectral line data as strings for header

            unames = []
            flags, areas, decays = self.metab_list.get_metabolite_settings()
            indx = np.where(flags)[0]
            if len(indx) > 0:
                unames = np.array(self.datasim.names)[indx]
                decays = np.array(decays)[indx]
                areas = np.array(areas)[indx]
            unames = ", ".join(unames)
            areas = ", ".join([str(val) for val in areas])
            decays = ", ".join([str(val) for val in decays])

            r = self.grid_get_values(grid=self.GridMmol)
            flags, m_ppms, m_areas, m_phases, m_widths, m_widths_hz, m_widths_damp = r
            indx = np.where(flags)[0]
            if len(indx) > 0:
                m_ppms = np.array(m_ppms)[indx]
                m_areas = np.array(m_areas)[indx]
                m_phases = np.array(m_phases)[indx]
                m_widths = np.array(m_widths)[indx]
                m_widths_hz = np.array(m_widths_hz)[indx]
                m_widths_damp = np.array(m_widths_damp)[indx]
            m_ppms = ", ".join([str(val) for val in m_ppms])
            m_areas = ", ".join([str(val) for val in m_areas])
            m_phases = ", ".join([str(val) for val in m_phases])
            m_widths = ", ".join([str(val) for val in m_widths])
            m_widths_hz = ", ".join([str(val) for val in m_widths_hz])
            m_widths_damp = ", ".join([str(val) for val in m_widths_damp])

            r = self.grid_get_values(grid=self.GridBase)
            flags, b_ppms, b_areas, b_phases, b_widths, b_widths_hz, b_widths_damp = r
            indx = np.where(flags)[0]
            if len(indx) > 0:
                b_ppms = np.array(b_ppms)[indx]
                b_areas = np.array(b_areas)[indx]
                b_phases = np.array(b_phases)[indx]
                b_widths = np.array(b_widths)[indx]
                b_widths_hz = np.array(b_widths)[indx]
                b_widths_damp = np.array(b_widths)[indx]
            b_ppms = ", ".join([str(val) for val in b_ppms])
            b_areas = ", ".join([str(val) for val in b_areas])
            b_phases = ", ".join([str(val) for val in b_phases])
            b_widths = ", ".join([str(val) for val in b_widths])
            b_widths_hz = ", ".join([str(val) for val in b_widths_hz])
            b_widths_damp = ", ".join([str(val) for val in b_widths_damp])

            stamp = util_time.now(util_time.ISO_TIMESTAMP_FORMAT).split('T')
        
            lines     = ['Vespa-DataSim - An Application for Simulated MRS Data']
            lines.append('-----------------------------------------------------')
            lines.append('The following information is a summary of the settings')
            lines.append('used to generate the enclosed MRS data.')
            lines.append(' ')
            lines.append('Creation_date                 - '+stamp[0])
            lines.append('Creation_time                 - '+stamp[1])
            lines.append(' ')
            lines.append('Original Export Filename      - '+filename)
            lines.append('User Comment                  - '+ds.comment)
            lines.append(' ')
            lines.append('Prior Experiment Name         - '+ds.experiment.name)
            lines.append('Prior Experiment ID           - '+ds.experiment.id)
            lines.append('Metabolite Names              - '+unames)
            lines.append('Metabolite Area Values        - '+areas)
            lines.append('Metabolite Ta Decays          - '+decays)
            lines.append('Metabolite Tb Decay           - '+str(ds.tb))
            lines.append('Metabolite Range Start        - '+str(ds.mets_ppm_start))
            lines.append('Metabolite Range End          - '+str(ds.mets_ppm_end))
            lines.append('Macromolecule Lineshape       - '+ds.mmol_lineshape)
            lines.append('Macromolecule Group Scale     - '+str(ds.mmol_group_scale))
            lines.append('Macromolecule PPM Values      - '+m_ppms)
            lines.append('Macromolecule Area Values     - '+m_areas)
            lines.append('Macromolecule Phase Values    - '+m_phases)
            lines.append('Macromolecule Peak Widths PPM - '+m_widths)
            lines.append('Macromolecule Peak Widths Hz  - '+m_widths_hz)
            lines.append('Macromolecule Peak Tb Decays  - '+m_widths_damp)
            lines.append('Baseline Lineshape            - '+ds.base_lineshape)
            lines.append('Baseline Group Scale          - '+str(ds.base_group_scale))
            lines.append('Baseline PPM Values           - '+b_ppms)
            lines.append('Baseline Area Values          - '+b_areas)
            lines.append('Baseline Phase Values         - '+b_phases)
            lines.append('Baseline Peak Widths PPM      - '+b_widths)
            lines.append('Baseline Peak Widths Hz       - '+b_widths_hz)
            lines.append('Baseline Peak Tb Decays       - '+b_widths_damp)
            lines.append('Noise RMS Multiplier          - '+str(ds.noise_rms_multiplier / 100.0))
            lines.append('Noise Ref Peak Area           - '+str(ds.noise_ref_peak_area))
            lines.append('Noise Ref Peak Ta             - '+str(ds.noise_ref_peak_ta))
            lines.append('Noise Ref Peak Tb             - '+str(ds.noise_ref_peak_tb))
            lines.append('Measured SNR (actual)         - '+str(self.snr_actual))
            lines.append('Spectrum Phase0               - '+str(ds.phase0))
            lines = "\n".join(lines)
            # if (wx.Platform == "__WXMSW__"):
            #     lines = lines.replace("\n", "\r\n")

            datasets = [] 
            msg = ''   
            filename, _ = os.path.splitext(filename)
            out_filenames = [ ]
            for extension in ('metabolites', 'macromolecules', 'baselines', 'noise', 'summed'):
                out_filenames.append("%s_%s.xml" % (filename, extension))

            for i,item in enumerate(results):
                raw = mrs_data_raw.DataRaw() 
                raw.data_sources = [out_filenames[i]]
                raw.headers = [lines]
                raw.sw = ds.sw
                raw.frequency = ds.frequency
                raw.resppm = ds.resppm
                raw.data = item.copy()
                filename = out_filenames[i]
                try:
                    util_export.export(filename, [raw], None, comment, False)
                except IOError:
                    msg = """I can't write the file "%s".""" % filename
                else:
                    path, _ = os.path.split(filename)
                    util_datasim_config.set_last_export_path(path)
                
                if msg:
                    common_dialogs.message(msg, style=common_dialogs.E_OK)
    
    
    def export_monte_carlo_to_viff(self):
        comment = 'Exported from Vespa-Datasim using Export Monte Carlo to VIFF'
        self.export_spectrum_to_viff(monte_carlo_flag=True)
        return


    def export_spectrum_to_siemens_rda(self):
        """ 
        Generates Siemens *.RDA formatted fake data
        
        Five files are output, four data and one provenance. The user selects a 
        base filename (fbase) and then extension strings are added to fbase to
        indicate which data is in the file
        
        1. fbase_metabolites - time domain summed metabolite signals, no noise
        2. fbase_baseline    - time domain summed baseline signals, no noise
        3. fbase_noise       - the time domain added noise in the summed data
        4. fbase_summed      - the sum of the three above file, final spectrum
        
        5. fbase_provenance - text file, info about Vespa-Datasim setting used
                              to create the above files. This is created because
                              Siemens RDA has limited ability to store provenance
                              information within its header 
        """
        datasim = self.datasim

        filter = "Pick Filename for Siemens RDA file (*.rda)|*.rda"
        filename = common_dialogs.save_as(filetype_filter=filter)
        
        if filename: 

            hdr   = self.make_rda_header()
            phase = np.exp(1j * datasim.phase0 * DEGREES_TO_RADIANS)
            data  = self.metabolites_time_sum.copy() * phase
            mmol  = self.macromolecule_time_sum.copy() * phase
            base  = self.baseline_time_sum.copy() * phase
            noise = self.noise_time.copy() * phase
            final = (self.metabolites_time_sum + self.macromolecule_time_sum + self.baseline_time_sum + self.noise_time).copy() * phase
            snra  = self.snr_actual
                
            results = [data, mmol, base, noise, final]
            
            #-------------------------------------------
            
            filename, _ = os.path.splitext(filename)
            rda_filenames = [ ]
            for extension in ('metabolites', 'macromolecules', 'baselines', 'noise', 'summed'):
                rda_filenames.append("%s_%s.rda" % (filename, extension))
            fname_provenance = "%s_provenance.txt" % filename

            # write the four numpy results arrays to file
            for i, item in enumerate(results):
                try:
                    # write the header to a file
                    with open(rda_filenames[i], 'w') as f:
                        for line in hdr:
                            f.write("%s\n" % line)
                    
                    # append data after header text
                    format = 'd'
                    item = np.complex128(item)
                    item = item.ravel().tolist()
                    item = util_fileio.expand_complexes(item)
                    item = struct.pack(format * len(item), *item)
                    
                    f = open(rda_filenames[i], 'ab')
                    f.write(item)
                    f.close()
                        
                except IOError:
                    msg = "I was unable to write the following file:\n"
                    msg += rda_filenames[i] + '\n'
                    common_dialogs.message(msg, style=common_dialogs.E_OK)   
                
                # write the Datasim provenance to text file    
                provenance = str(datasim)
                f = open(fname_provenance,'w')
                f.write(provenance)
                f.close()

    def make_rda_header(self):
        
        dt = datetime.now()
        _date = '{:%Y%m%d}'.format(dt)
        _time = '{:%Y%m%d}'.format(dt)

        sw    = self.datasim.sw
        dim0  = self.datasim.dims[0]
        freq  = self.datasim.frequency
        nuc   = '1H'
        field = '0.1'
        if freq > 60 and freq < 68:
            field = '1.5'
        elif freq > 80 and freq < 100:
            field = '2.0'
        elif freq > 120 and freq < 130:
            field = '3.0'
        elif freq > 170 and freq < 190:
            field = '4.0'
        elif freq > 200 and freq < 225:
            field = '4.7'
        elif freq > 290 and freq < 335:
            field = '7.0'
        dwell = 1000000.0 / sw

        hdr = []
        hdr.append(">>> Begin of header <<<")
        hdr.append("PatientName: Vespa-Datasim Simulated Data")
        hdr.append("PatientID: Vespa-Datasim")
        hdr.append("PatientSex: O")
        hdr.append("PatientBirthDate: 19190101")
        hdr.append("StudyDate: {0}".format(_date))
        hdr.append("StudyTime: 100000.000000")
        hdr.append("StudyDescription: Vespa-Datasim Simulated Data")
        hdr.append("PatientAge: 050Y")
        hdr.append("PatientWeight: 60.0000")
        hdr.append("SeriesDate: {0}".format(_date))
        hdr.append("SeriesTime: 100000.000000")
        hdr.append("SeriesDescription: Vespa-Datasim Simulated Data")
        hdr.append("ProtocolName: Vespa-Datasim Simulated Data")
        hdr.append("PatientPosition: HFS")
        hdr.append("SeriesNumber: 1")
        hdr.append("InstitutionName: Vespa-Datasim")
        hdr.append("StationName: Vespa-Datasim")
        hdr.append("ModelName: Vespa-Datasim")
        hdr.append("DeviceSerialNumber: 76543")
        hdr.append("SoftwareVersion[0]: Vespa-Datasim")
        hdr.append("InstanceDate: {0}".format(_date))
        hdr.append("InstanceTime: 100000.000000")
        hdr.append("InstanceNumber: 1")
        hdr.append("InstanceComments: Vespa-Datasim Simulated Data")
        hdr.append("AcquisitionNumber: 1")
        hdr.append("SequenceName: {0}".format(self.datasim.experiment.name))
        hdr.append("SequenceDescription: {0}".format(self.datasim.experiment.id))
        hdr.append("TR: 1000.000000")
        hdr.append("TE: 30.000000")
        hdr.append("TM: 0.000000")
        hdr.append("TI: 0.000000")
        hdr.append("DwellTime: {0}".format(str(dwell)))  # sw in Hz, dwell in usec
        hdr.append("EchoNumber: 0")
        hdr.append("NumberOfAverages: 1.000000")
        hdr.append("MRFrequency: {0}".format(str(freq)))            # in MHz
        hdr.append("Nucleus: {0}".format(nuc))                      # string
        hdr.append("MagneticFieldStrength: {0}".format(str(field))) # should be float 1.5 or 3.0 or 7.0 etc.
        hdr.append("NumOfPhaseEncodingSteps: 1")
        hdr.append("FlipAngle: 90.000000")
        hdr.append("VectorSize: {0}".format(dim0))                  # should be int
        hdr.append("CSIMatrixSize[0]: 1")
        hdr.append("CSIMatrixSize[1]: 1")
        hdr.append("CSIMatrixSize[2]: 1")
        hdr.append("CSIMatrixSizeOfScan[0]: 1")
        hdr.append("CSIMatrixSizeOfScan[1]: 1")
        hdr.append("CSIMatrixSizeOfScan[2]: 1")
        hdr.append("CSIGridShift[0]: 0")
        hdr.append("CSIGridShift[1]: 0")
        hdr.append("CSIGridShift[2]: 0")
        hdr.append("HammingFilter: Off")
        hdr.append("FrequencyCorrection: NO")
        hdr.append("TransmitCoil: TxRx_Head")
        hdr.append("TransmitRefAmplitude[1H]: 100.000000")
        hdr.append("SliceThickness: 20.000000")
        hdr.append("PositionVector[0]: 0.000000")
        hdr.append("PositionVector[1]: 0.000000")
        hdr.append("PositionVector[2]: 0.000000")
        hdr.append("RowVector[0]: -1.000000")
        hdr.append("RowVector[1]: 0.000000")
        hdr.append("RowVector[2]: 0.000000")
        hdr.append("ColumnVector[0]: 0.000000")
        hdr.append("ColumnVector[1]: 1.000000")
        hdr.append("ColumnVector[2]: 0.000000")
        hdr.append("VOIPositionSag: 0.000000")
        hdr.append("VOIPositionCor: 0.000000")
        hdr.append("VOIPositionTra: 0.000000")
        hdr.append("VOIThickness: 20.000000")
        hdr.append("VOIPhaseFOV: 20.000000")
        hdr.append("VOIReadoutFOV: 20.000000")
        hdr.append("VOINormalSag: 0.000000")
        hdr.append("VOINormalCor: 0.000000")
        hdr.append("VOINormalTra: 1.000000")
        hdr.append("VOIRotationInPlane: 0.000000")
        hdr.append("FoVHeight: 20.000000")
        hdr.append("FoVWidth: 20.000000")
        hdr.append("FoV3D: 20.000000")
        hdr.append("PercentOfRectFoV: 1.000000")
        hdr.append("NumberOfRows: 1")
        hdr.append("NumberOfColumns: 1")
        hdr.append("NumberOf3DParts: 1")
        hdr.append("PixelSpacingRow: 20.000000")
        hdr.append("PixelSpacingCol: 20.000000")
        hdr.append("PixelSpacing3D: 20.000000")
        hdr.append(">>> End of header <<<")
        
        return hdr

    
    def process_and_plot(self, set_scale=False):
        self.process()
        self.plot(set_scale=set_scale)


    def process(self):
        """ 
        Converts spectral (metabolite, macromolecule, baseline and noise) signals
        from ideal time basis functions to spectral domain peaks by applying a
        line shape envelope, scaling and phasing and then applying the FFT.
        
        """
        if not self.plotting_enabled: 
            return

        ds   = self.datasim
        loop = ds.loop
        dims = ds.dims
        xx   = np.arange(dims[0]) / ds.sw
        b0   = np.exp(1j * ds.b0shift * 2.0 * np.pi * xx)       # in Radians

        #----------------------------------------------------------------------
        # METABOLITE Section
        #
        # The time signal in both metabs and baseline have _chop() applied to
        # shift the water signal to the zero frequency location. This is so
        # that the hlsvd algorithm (in Analysis) works appropriately. However,
        # this does require the Analysis Spectral tab to default to 'Chop' ON.

        flags, scales, decays = self.metab_list.get_metabolite_settings()

        ds.mets_flags = flags
        ds.mets_scales = scales
        ds.mets_decays = decays

        # get the proper loop and metabolite names
        metab_time     = np.zeros((1,dims[0]), np.complex128)
        metab_freq     = np.zeros((1,dims[0]), np.complex128)
        metab_freq_sum = np.zeros((1,dims[0]), np.complex128)
        metab_time_sum = np.zeros((1,dims[0]), np.complex128)

        indx  = np.where(flags)[0]
        if len(indx) > 0:
            metab_time = ds.basis[loop[2],loop[1],loop[0],indx,:].copy()
            metab_freq = np.zeros((len(indx),dims[0]), np.complex128)
            arr_ta = np.array(decays)[indx]
            area   = np.array(scales)[indx]
            tb     = -(xx/ds.tb)**2

            for i,time in enumerate(metab_time):
                ta   = -xx/arr_ta[i]
                time = time * area[i] * util_math.safe_exp(ta+tb) * b0

                # apply left_shift early, as per Analysis
                if ds.left_shift:
                    time = np.roll(time, -ds.left_shift)
                    time[-ds.left_shift:] = time[0]*0.0
                
                metab_time[i,:] = _chop(time.copy())  # consistent with Analysis
                time[0] *= 0.5
                metab_freq[i,:] = np.fft.fft(time) / dims[0]     

                if len(indx) > 1:
                    metab_freq_sum = np.sum(metab_freq, axis=0)
                    metab_time_sum = np.sum(metab_time, axis=0)
                else:
                    metab_freq_sum = metab_freq[0,:]
                    metab_time_sum = metab_time[0,:]

        self.metabolites_all      = metab_freq
        self.metabolites_sum      = metab_freq_sum
        self.metabolites_time_all = metab_time
        self.metabolites_time_sum = metab_time_sum

        #----------------------------------------------------------------------
        # MACROMOLECULE Section

        mmol_time     = np.zeros((1,dims[0]), np.complex128)
        mmol_freq     = np.zeros((1,dims[0]), np.complex128)
        mmol_time_sum = np.zeros((1,dims[0]), np.complex128)
        mmol_freq_sum = np.zeros((1,dims[0]), np.complex128)

        if ds.macromolecule_basis is not None:
            r = self.grid_get_values(grid=self.GridMmol)
            flags, ppms, areas, phases, widths, _, _ = r

            ds.mmol_flags = flags
            ds.mmol_ppms  = ppms
            ds.mmol_areas = areas
            ds.mmol_phases = phases
            ds.mmol_widths = widths
            
            indx = np.where(flags)[0]
            if len(indx) > 0:
                mmol_time = ds.macromolecule_basis[indx,:].copy() * ds.mmol_group_scale
                mmol_freq = np.zeros((len(indx),dims[0]), np.complex128)
                for i, time in enumerate(mmol_time):
                    time = time.copy() * b0
                    mmol_time[i,:] = _chop(mmol_time[i,:]) # consistent with Analysis
                    time[0] *= 0.5
                    mmol_freq[i,:] = np.fft.fft(time) / dims[0]

                if len(indx) > 1:
                    mmol_freq_sum = np.sum(mmol_freq, axis=0)
                    mmol_time_sum = np.sum(mmol_time, axis=0)
                else:
                    mmol_freq_sum = mmol_freq[0,:]
                    mmol_time_sum = mmol_time[0,:]

        self.macromolecule_all      = mmol_freq
        self.macromolecule_sum      = mmol_freq_sum
        self.macromolecule_time_all = mmol_time
        self.macromolecule_time_sum = mmol_time_sum

        #----------------------------------------------------------------------
        # BASELINE Section

        base_time     = np.zeros((1,dims[0]), np.complex128)
        base_freq     = np.zeros((1,dims[0]), np.complex128)
        base_time_sum = np.zeros((1,dims[0]), np.complex128)
        base_freq_sum = np.zeros((1,dims[0]), np.complex128)

        if ds.baseline_basis is not None:
            r = self.grid_get_values(grid=self.GridBase)
            flags, ppms, areas, phases, widths, _, _ = r

            ds.base_flags = flags
            ds.base_ppms  = ppms
            ds.base_areas = areas
            ds.base_phases = phases
            ds.base_widths = widths

            indx = np.where(flags)[0]
            if len(indx) > 0:
                base_time = ds.baseline_basis[indx, :].copy() * ds.base_group_scale
                base_freq = np.zeros((len(indx), dims[0]), np.complex128)
                for i, time in enumerate(base_time):
                    time = time.copy() * b0
                    base_time[i, :] = _chop(base_time[i, :])  # consistent with Analysis
                    time[0] *= 0.5
                    base_freq[i, :] = np.fft.fft(time) / dims[0]

                if len(indx) > 1:
                    base_freq_sum = np.sum(base_freq, axis=0)
                    base_time_sum = np.sum(base_time, axis=0)
                else:
                    base_freq_sum = base_freq[0, :]
                    base_time_sum = base_time[0, :]

        self.baseline_all      = base_freq
        self.baseline_sum      = base_freq_sum
        self.baseline_time_all = base_time
        self.baseline_time_sum = base_time_sum

        #----------------------------------------------------------------------
        # NOISE Section
            
        self.noise_time = np.zeros(dims[0], np.complex128)
        self.noise_freq = np.zeros(dims[0], np.complex128)
        self.snr_actual = 0.0
        
        if self.display_noise_in_plot:
            ref_peak_area    = ds.noise_ref_peak_area
            ref_peak_ta      = ds.noise_ref_peak_ta
            ref_peak_tb      = ds.noise_ref_peak_tb
            noise_multiplier = ds.noise_rms_multiplier / 100.0

            ta = -xx/ref_peak_ta
            tb = -(xx/ref_peak_tb)**2
            nexpo = ref_peak_area * util_math.safe_exp(ta+tb)
            nexpo[0] *= 0.5
            noise_time = np.random.randn(dims[0]) * noise_multiplier
            rms_noise_time = self.rms(noise_time.real)

            noise_freq     =  np.fft.fft(noise_time) / dims[0]
            noise_time     = _chop(noise_time)                          # consistent with Analysis
            refpeak_height = (np.fft.fft(nexpo) / dims[0]).real[0]
            rms_noise_freq = self.rms(noise_freq.real) 
            snr_actual     = refpeak_height / rms_noise_freq

            self.noise_freq = noise_freq
            self.noise_time = noise_time
            self.snr_actual = snr_actual

    
    def plot(self, is_replot=False, set_scale=False):
        # The parameter, is_replot, is used to determine if the plot is
        # the initial plot or a replot.  The variable will fix the axis
        # xscale to the full width of the spectrum for the initial plot,
        # and the scale will remain constant if it is a replot.  For
        # example, the xscale should remain constant during phasing, but
        # should reset if a new spectrum is loaded.

        if not self.plotting_enabled: 
            return
        
        if self.metabolites_all is None:
            return

        data1 = {'data' : (self.metabolites_sum + self.macromolecule_sum + self.baseline_sum + self.noise_freq),
                 'line_color_real'      : self._prefs.line_color_real,
                 'line_color_imaginary' : self._prefs.line_color_imaginary,
                 'line_color_magnitude' : self._prefs.line_color_magnitude }

        data2 = {'data' : self.metabolites_all, 
                 'line_color_real'      : self._prefs.line_color_metabolite,
                 'line_color_imaginary' : self._prefs.line_color_metabolite,
                 'line_color_magnitude' : self._prefs.line_color_metabolite }

        data3 = {'data' : np.concatenate([self.baseline_all,self.macromolecule_all],axis=0),
                 'line_color_real'      : self._prefs.line_color_baseline,
                 'line_color_imaginary' : self._prefs.line_color_baseline,
                 'line_color_magnitude' : self._prefs.line_color_baseline }

        data = [[data1], [data2], [data3]]
        self.view.set_data(data)
        self.view.update(no_draw=True, set_scale=set_scale)
        self.view.set_phase_0(0.0)


