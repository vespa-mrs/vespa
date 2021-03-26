# Python modules
import warnings

# 3rd party modules
import wx
import numpy as np
import matplotlib.cm as cm
from pubsub import pub as pubsub

# Our modules
import vespa.analysis.tab_base as tab_base
import vespa.analysis.prefs as prefs_module
import vespa.analysis.util_menu as util_menu
import vespa.analysis.block_raw_edit_fidsum as block_raw_edit_fidsum
import vespa.analysis.auto_gui.fidsum as fidsum
import vespa.analysis.functors.funct_fidsum_coil_combine as coil_combine
import vespa.analysis.functors.funct_fidsum_exclude as exclude_methods
import vespa.analysis.functors.funct_fidsum_correction as correction_methods
import vespa.analysis.dialog_dataset_browser as dialog_dataset_browser
import vespa.common.wx_gravy.util as wx_util
import vespa.common.wx_gravy.common_dialogs as common_dialogs

from vespa.analysis.dialog_user_prior import DialogUserPrior
from vespa.analysis.plot_panel_prep_fidsum import   PlotPanelPrepFidsum, \
                                                    PlotPanelPrepFidsumSeries, \
                                                    PlotPanelPrepFidsumImage


EXCLUDE_DISPLAY_CHOICES = ["FID Abs First Point",
                           "Peak Shift [Hz]",
                           "Peak Phase [deg]"
                          ]

TARGET_CHOICES = ["Average all", "Avg first 4", "Avg first 10%", "Avg first 25%", \
                  "Avg first 50%", "Avg mid 10%", "Avg mid 30%"]    # not included yet ... , "User define"]

#------------------------------------------------------------------------------
# GUI Helper methods

def _configure_combo(control, lines, selection=''):
    control.Clear()
    control.Set(lines)
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
#
#  Tab PREP FIDSUM
#
#------------------------------------------------------------------------------

class TabPrepFidsum(tab_base.Tab, fidsum.PanelPrepFidsumUI):

    # self-identify tab to notebook, value does not matter, its presence is sufficient.
    IS_PREP_FIDSUM = True
    
    def __init__(self, tab_dataset, top, block):
        
        fidsum.PanelPrepFidsumUI.__init__(self, tab_dataset.NotebookDataset)
        
        tab_base.Tab.__init__(self, tab_dataset, top, prefs_module.PrefsPrepFidsum)
        
        self.top      = top
        self.block    = block

        # Disable plotting during init because some controls fire change events
        # as they are inited triggering a call to plot().
        self._plotting_enabled = False
        self.plot_results = None

        self.waterfall_range_start = 3.99       # [ppm]
        self.waterfall_range_end   = 1.01

        self.fid_index = 0
        self.initialize_controls()
        self.populate_controls()

        self._plotting_enabled = True
        self._update_freq_raw  = True

        # Setup the canvas
        self.process_and_plot()

        # PubSub subscriptions
        pubsub.subscribe(self.on_push_prep_fidsum_results, "push_prep_fidsum_results")


        # Set sash position to INI file value or arbitrary-ish value of 400.
        if not self._prefs.sash_position:
            self._prefs.sash_position = 400

        # Under OS X, wx sets the sash position to 10 (why 10?) *after*
        # this method is done. So setting the sash position here does no
        # good. We use wx.CallAfter() to (a) set the sash position and
        # (b) fake an EVT_SPLITTER_SASH_POS_CHANGED.
        wx.CallAfter(self.SplitterWindow.SetSashPosition, self._prefs.sash_position, True)
        wx.CallAfter(self.on_splitter)


    #=======================================================
    #
    #           GUI Setup Handlers 
    #
    #=======================================================

    def initialize_controls(self):
        """ 
        One-time setup for controls: size, range, number of decimal places.
        Default values done in populate_controls().
        
        """
        dataset = self.dataset
        dim0, dim1, dim2, dim3 = dataset.spectral_dims
        nrep, ncoil, nfid, npts = dataset.raw_shape
        sw      = dataset.sw
        maxppm  = dataset.pts2ppm(0, acq=True)
        minppm  = dataset.pts2ppm(dim0-1, acq=True)
        ppmlim  = (minppm, maxppm)


        wx_util.configure_spin(self.FloatExclusionAutoFidaBadThreshold,70, 2, 0.5, (0, 100))
        wx_util.configure_spin(self.SpinExclusionAutoFidaWorstN,       70, None, None,(0, 100))

        wx_util.configure_spin(self.FloatVespaReferencePeakCenter,     70, 3, 0.1, (-30.0, 30.0))
        wx_util.configure_spin(self.FloatVespaPeakSearchWidth,         70, 3, 0.1, (0.0, 100.0))
        wx_util.configure_spin(self.FloatVespaPhase0RangeStart,        70, 2, 0.1, ppmlim)
        wx_util.configure_spin(self.FloatVespaPhase0RangeEnd,          70, 2, 0.1, ppmlim)

        wx_util.configure_spin(self.FloatSuspectInitialGuessFreqShift, 70, 2, 0.5, (-1000.0, 1000.0))
        wx_util.configure_spin(self.FloatSuspectInitialGuessPhase,     70, 2, 0.5, (-360.0, 360.0))
        wx_util.configure_spin(self.FloatSuspectOptimizeRangeStart,    70, 2, 0.1, ppmlim)
        wx_util.configure_spin(self.FloatSuspectOptimizeRangeEnd,      70, 2, 0.1, ppmlim)

        wx_util.configure_spin(self.FloatRatsInitialGuessFreqShift,    70, 2, 0.5, (-1000.0, 1000.0))
        wx_util.configure_spin(self.FloatRatsInitialGuessPhase,        70, 2, 0.5, (-360.0, 360.0))
        wx_util.configure_spin(self.FloatRatsOptimizeRangeStart,       70, 2, 0.1, ppmlim)
        wx_util.configure_spin(self.FloatRatsOptimizeRangeEnd,         70, 2, 0.1, ppmlim)
        wx_util.configure_spin(self.SpinRatsBaselineOrder,             70, None, None, (0, 10))

        wx_util.configure_spin(self.SpinFidIndex, 60)
        wx_util.configure_spin(self.FloatCurrentPeakShiftValue,   70, 3, 0.5, (-1000.0, 1000.0))
        wx_util.configure_spin(self.FloatCurrentPhase0Value,      70, 3, 0.5, (-1000.0, 1000.0))
        wx_util.configure_spin(self.FloatWaterfallRangeStart,     70, 3, 0.5, ppmlim)
        wx_util.configure_spin(self.FloatWaterfallRangeEnd,       70, 3, 0.5, ppmlim)

        wx_util.configure_spin(self.FloatGlobalPhase0,            70, 1,  1.0, (-360.0, 360.0))
        wx_util.configure_spin(self.FloatGlobalPhase1,            70, 1, 50.0, (-50000.0, 50000.0))
        wx_util.configure_spin(self.SpinGlobalLeftShift,          70, None, None, (0, 50))
        wx_util.configure_spin(self.FloatGlobalGaussApodize,      70, 1, 1.0, (0.0, 1000.0))

        # Set this Combo here as it is part of the Tab display options, not the Block
        _configure_combo(self.ComboDataExclusionPlotDisplay, EXCLUDE_DISPLAY_CHOICES,
                         selection=EXCLUDE_DISPLAY_CHOICES[0])

        # self.spintxt = self.FloatGlobalGaussApodize.Children[0]
        # self.spintxt.Bind(wx.EVT_CHAR_HOOK, self.on_spin_txt_char_hook)

        #-------------------------------------------------------------
        # Raw Fidsum View setup 

        self.view = PlotPanelPrepFidsum(self.PanelViewPrepFidsum,
                                        self,
                                        self._tab_dataset,
                                        naxes=2,
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
        self.PanelViewPrepFidsum.SetSizer(sizer)
        self.view.Fit()

        self.view_series = PlotPanelPrepFidsumSeries(self.PanelViewPrepFidsumSeries,
                                                     self,
                                                     self._tab_dataset,
                                                     naxes=1,
                                                     reversex=False,
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
                                                     do_middle_select_event=False,
                                                     do_middle_motion_event=False,
                                                     do_middle_press_event=True,
                                                     do_scroll_event=True,
                                                     props_zoom=dict(alpha=0.2, facecolor='yellow'),
                                                     props_cursor=dict(alpha=0.2, facecolor='gray'),
                                                     xscale_bump=0.0,
                                                     yscale_bump=0.05,
                                                     data = [],
                                                     prefs=self._prefs,
                                                     )

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.view_series, 1, wx.LEFT | wx.TOP | wx.EXPAND)
        self.PanelViewPrepFidsumSeries.SetSizer(sizer)
        self.view_series.Fit()

        xstr, xend = (-0.05, 0.05) if nfid <= 1 else (0, nfid-1)
        self.view_series.axes[0].set_xlim(xstr, xend)

        self.view_image = PlotPanelPrepFidsumImage(self.PanelViewPrepFidsumImage,
                                                   self,
                                                   naxes=2,
                                                   reversex=False,
                                                   zoom='box',
                                                   reference=False,
                                                   middle=True,
                                                   zoom_button=1,
                                                   middle_button=3,
                                                   do_zoom_select_event=False,
                                                   do_zoom_motion_event=True,
                                                   do_refs_select_event=False,
                                                   do_refs_motion_event=False,
                                                   do_middle_press_event=True,
                                                   props_zoom=dict(alpha=0.2, facecolor='yellow'),
                                                   props_cursor=dict(alpha=0.1, facecolor='gray'))
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.view_image, 1, wx.LEFT | wx.TOP | wx.EXPAND)
        self.PanelViewPrepFidsumImage.SetSizer(sizer)
        self.view_image.Fit()

        self.view_image.figure.set_facecolor('white')
        self.view_image.axes[0].set_facecolor('gray')
        self.view_image.axes[0].set_xlim(xstr, xend)
        self.view_image.axes[0].set_ylim(self.waterfall_range_start, self.waterfall_range_end)
        self.view_image.axes[0].axes.set_xticklabels([])
#        self.view_image.axes[0].axes.set_yticklabels([])
        self.view_image.axes[1].set_facecolor('gray')
        self.view_image.axes[1].set_xlim(xstr, xend)
        self.view_image.axes[1].set_ylim(self.waterfall_range_start, self.waterfall_range_end)
#        self.view_image.axes[1].axes.set_xticklabels([])
#        self.view_image.axes[1].axes.set_yticklabels([])


        # width, height = self.SpinFidLeftShift.GetSize()
        #
        # if "__WXMAC__" in wx.PlatformInfo:
        #     # Under OS X this spin control needs a poke to paint itself
        #     # properly. Without this code, the control doesn't appear on screen
        #     # even though the wx Inspector reports that it has an appropriate
        #     # size & position. Go figger.
        #     wx.CallAfter(self.SpinFidLeftShift.SetSize, (width + 1, height))

  
    def populate_controls(self, preset=False):
        """ 
        Populate controls with values from the dataset['prep'] block object.
        Called when a new data object is loaded.

        """
        block = self.block
        set   = self.block.set

        raw = self._tab_dataset.dataset.get_source_data('prep')

        # These Combos set in here they are part of the Block

        if set.coil_combine_method == '': set.coil_combine_method = coil_combine.COILCOMBINE_MENU_ITEMS[0]
        if set.exclude_method      == '': set.exclude_method = exclude_methods.EXCLUDE_MENU_ITEMS[0]
        if set.correction_method   == '': set.correction_method = correction_methods.CORRECTION_MENU_ITEMS[0]

        _configure_combo(self.ComboCoilCombineMethod, 
                         coil_combine.COILCOMBINE_MENU_ITEMS, 
                         selection=set.coil_combine_method)

        _configure_combo(self.ComboExclusionMethod, 
                         exclude_methods.EXCLUDE_MENU_ITEMS, 
                         selection=set.exclude_method)

        _configure_combo(self.ComboCorrectionMethod, 
                         correction_methods.CORRECTION_MENU_ITEMS, 
                         selection=set.correction_method)

        _configure_combo(self.ComboTargetSpectrumVespa,
                         TARGET_CHOICES,
                         selection=set.vespa_target_method)

        _configure_combo(self.ComboTargetSpectrumSuspectSpectral,
                         TARGET_CHOICES,
                         selection=set.suspect_target_method)

        _configure_combo(self.ComboTargetSpectrumSuspectRats,
                         TARGET_CHOICES,
                         selection=set.rats_target_method)


        method = set.coil_combine_method
        if method != 'External Dataset' and method != 'External Dataset with Offset':
            self.PanelCoilCombineExternalDataset.Hide()
        else:
            self.PanelCoilCombineExternalDataset.Show()
            self.TextCoilCombineExternalDataset.SetValue(self.block.set.coil_combine_external_filename)

        method = set.exclude_method
        if method == 'Manual':
            self.PanelExclusionManualOptions.Show()
            self.PanelExclusionFidaRmBad.Hide()
            self.PanelExclusionFidaWorstN.Hide()
        elif method == 'Remove Bad Averages (fid-a)':
            self.PanelExclusionManualOptions.Hide()
            self.PanelExclusionFidaRmBad.Show()
            self.PanelExclusionFidaWorstN.Hide()
        elif method == 'Remove N Worst Averages (fid-a)':
            self.PanelExclusionManualOptions.Hide()
            self.PanelExclusionFidaRmBad.Hide()
            self.PanelExclusionFidaWorstN.Show()

        method = set.correction_method
        if method == 'Manual':
            self.PanelCorrectionsVespa.Hide()
            self.PanelCorrectionsVespaCorrelate.Hide()
            self.PanelCorrectionsSuspectSpectralRegistration.Hide()
            self.PanelCorrectionsSuspectRats.Hide()
        elif method == 'Optimized Search (vespa)':
            self.PanelCorrectionsVespa.Show()
            self.PanelCorrectionsVespaCorrelate.Hide()
            self.PanelCorrectionsSuspectSpectralRegistration.Hide()
            self.PanelCorrectionsSuspectRats.Hide()
        elif method == 'Correlation (vespa)':
            self.PanelCorrectionsVespa.Hide()
            self.PanelCorrectionsVespaCorrelate.Show()
            self.PanelCorrectionsSuspectSpectralRegistration.Hide()
            self.PanelCorrectionsSuspectRats.Hide()
        elif method == 'Spectral Registration (suspect)':
            self.PanelCorrectionsVespa.Hide()
            self.PanelCorrectionsVespaCorrelate.Hide()
            self.PanelCorrectionsSuspectSpectralRegistration.Show()
            self.PanelCorrectionsSuspectRats.Hide()
        elif method == 'RATS (suspect)':
            self.PanelCorrectionsVespa.Hide()
            self.PanelCorrectionsVespaCorrelate.Hide()
            self.PanelCorrectionsSuspectSpectralRegistration.Hide()
            self.PanelCorrectionsSuspectRats.Show()


        # only turn on combine method widget if data has indiv coil FIDs in it        
        if raw.shape[1] > 1:
            self.PanelCoilCombination.Enable(True)
        else:
            self.PanelCoilCombination.Enable(False)
        self.CheckApplyNoiseWhitening.SetValue(set.coil_combine_noise_whitening)

        self.CheckApplyDataExclusion.SetValue(set.apply_data_exclusion)
        self.CheckExclusionInputAdjust.SetValue(set.exclusion_input_adjust)
        self.FloatExclusionAutoFidaBadThreshold.SetValue(set.fida_bad_threshold)
        self.SpinExclusionAutoFidaWorstN.SetValue(set.fida_n_worst)
        
        self.CheckCorrectionInputAdjust.SetValue(set.correction_input_adjust)
        self.FloatVespaReferencePeakCenter.SetValue(set.vespa_reference_peak_center)
        self.FloatVespaPeakSearchWidth.SetValue(set.vespa_peak_search_width)
        self.FloatVespaPhase0RangeStart.SetValue(set.vespa_phase0_range_start)
        self.FloatVespaPhase0RangeEnd.SetValue(set.vespa_phase0_range_end)
        
        self.FloatSuspectInitialGuessFreqShift.SetValue(set.suspect_initial_guess_freq)
        self.FloatSuspectInitialGuessPhase.SetValue(set.suspect_initial_guess_phase)
        self.FloatSuspectOptimizeRangeStart.SetValue(set.suspect_optimization_range_start)
        self.FloatSuspectOptimizeRangeEnd.SetValue(set.suspect_optimization_range_end)
        
        self.FloatRatsInitialGuessFreqShift.SetValue(set.rats_initial_guess_freq)
        self.FloatRatsInitialGuessPhase.SetValue(set.rats_initial_guess_phase)
        self.FloatRatsOptimizeRangeStart.SetValue(set.rats_optimization_range_start)
        self.FloatRatsOptimizeRangeEnd.SetValue(set.rats_optimization_range_end)
        self.SpinRatsBaselineOrder.SetValue(set.rats_baseline_order)

        self.SpinFidIndex.SetValue(self.fid_index)
        self.SpinFidIndex.SetRange(0, raw.shape[-2]-1)
        self.FloatCurrentPeakShiftValue.SetValue(self.block.frequency_shift[self.fid_index])
        self.FloatCurrentPhase0Value.SetValue(self.block.phase_0[self.fid_index])
        self.FloatWaterfallRangeStart.SetValue(self.waterfall_range_start)
        self.FloatWaterfallRangeEnd.SetValue(self.waterfall_range_end)

        self.SpinGlobalLeftShift.SetValue(set.global_left_shift)
        self.FloatGlobalPhase0.SetValue(set.global_phase0)
        self.FloatGlobalPhase1.SetValue(set.global_phase1)
        self.FloatGlobalGaussApodize.SetValue(set.global_gaussian_apodization)
        self.CheckChop.SetValue(set.chop_data)
        self.CheckZeroGlobalPhase1.SetValue(set.zero_phase1)
        self.CheckApplyPeakShift.SetValue(set.apply_peak_shift)
        self.CheckApplyPhase0.SetValue(set.apply_phase0)


    
    #=======================================================
    #
    #           Global and Menu Event Handlers 
    #
    #=======================================================

    def on_activation(self):
        super().on_activation()
        # plot is necessary here to resize img0/img1 waterfall plots
        # after a dataset loads tried a bunch of other locations, and
        # CallLater/CallAfter, but only this works ... annoying.
        self.plot()

    # def on_spin_txt_char_hook(self, event):
    #     event.Skip()  # allow the event to propagate further
    #     print("Event handler 'on_spin_txt_char_hook'")
    #     code = event.GetKeyCode()
    #     if code == wx.WXK_RETURN:
    #         bob = 10
    #         self.on_global_gauss_apodize(event)


    def on_menu_view_option(self, event):
        event_id = event.GetId()

        if self._prefs.handle_event(event_id):

            if event_id in (util_menu.ViewIdsPrepFidsum.ZERO_LINE_SHOW,
                            util_menu.ViewIdsPrepFidsum.ZERO_LINE_TOP,
                            util_menu.ViewIdsPrepFidsum.ZERO_LINE_MIDDLE,
                            util_menu.ViewIdsPrepFidsum.ZERO_LINE_BOTTOM,
                            util_menu.ViewIdsPrepFidsum.XAXIS_SHOW,
                          ):
                self.view.update_axes()
                self.view.canvas.draw()

            if event_id in (util_menu.ViewIdsPrepFidsum.ZERO_LINE_PLOT_SHOW,
                            util_menu.ViewIdsPrepFidsum.ZERO_LINE_PLOT_TOP,
                            util_menu.ViewIdsPrepFidsum.ZERO_LINE_PLOT_MIDDLE,
                            util_menu.ViewIdsPrepFidsum.ZERO_LINE_PLOT_BOTTOM,
                            util_menu.ViewIdsPrepFidsum.XAXIS_SHOW,
                          ):
                self.view_series.update_axes()
                self.view_series.canvas.draw()

            # note. these need to come before next
            if event_id == util_menu.ViewIdsPrepFidsum.DATA_TYPE_REAL:
                self.view.set_data_type_real()
            if event_id == util_menu.ViewIdsPrepFidsum.DATA_TYPE_IMAGINARY:
                self.view.set_data_type_imaginary()
            if event_id == util_menu.ViewIdsPrepFidsum.DATA_TYPE_MAGNITUDE:
                self.view.set_data_type_magnitude()

            if event_id in (util_menu.ViewIdsPrepFidsum.DATA_TYPE_REAL,
                            util_menu.ViewIdsPrepFidsum.DATA_TYPE_IMAGINARY,
                            util_menu.ViewIdsPrepFidsum.DATA_TYPE_MAGNITUDE,
                            util_menu.ViewIdsPrepFidsum.XAXIS_PPM,
                            util_menu.ViewIdsPrepFidsum.XAXIS_HERTZ,
                          ):
                self.view.update()
                self.view.canvas.draw()

            if event_id in (util_menu.ViewIdsPrepFidsum.AREA_CALC_PLOT_A,
                            util_menu.ViewIdsPrepFidsum.AREA_CALC_PLOT_B,
                           ):
                area, rms = self.view.calculate_area()
                if self._prefs.area_calc_plot_a:
                    index = 0
                else:
                    index = 1
                self.top.statusbar.SetStatusText(self.build_area_text(area[index], rms[index]), 3)
                

    def on_menu_view_output(self, event):
        event_id = event.GetId()

        formats = { util_menu.ViewIdsPrepFidsum.VIEW_TO_PNG : "PNG",
                    util_menu.ViewIdsPrepFidsum.VIEW_TO_SVG : "SVG", 
                    util_menu.ViewIdsPrepFidsum.VIEW_TO_EPS : "EPS", 
                    util_menu.ViewIdsPrepFidsum.VIEW_TO_PDF : "PDF", 
                  }

        if event_id in formats:
            format = formats[event_id]
            lformat = format.lower()
            filter_ = "%s files (*.%s)|*.%s" % (format, lformat, lformat)
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
                                    papertype='letter', 
                                    format=None,
                                    transparent=False)
                except IOError:
                    msg = """I can't write the file "%s".""" % filename
                
                if msg:
                    common_dialogs.message(msg, style=common_dialogs.E_OK)


    #=======================================================
    #
    #           PubSub Events
    #
    #=======================================================

    def on_push_prep_fidsum_results(self, uids=None, shifts=None, phase0=None, phase1=None):
        """ 
        Pubsub notifications handler. 
        
        Calling tab sends list of all associated dataset ids, and new freq
        and phase values.  If self.dataset.id in the uids list, we update the
        results arrays, the widgets, and re-process and plot spectra
        
        """
        if self.dataset.id in uids:
        
            block = self.dataset.blocks['prep']
    
            block.frequency_shift    = shifts.copy()
            block.phase_0            = phase0.copy()
            block.set.global_phase1  = phase1

            self.FloatCurrentPeakShiftValue.SetValue(block.frequency_shift[self.fid_index])
            self.FloatCurrentPhase0Value.SetValue(block.phase_0[self.fid_index])
            self.FloatGlobalPhase1.SetValue(block.set.global_phase1)
        
            self.process_and_plot()


    #=======================================================
    #
    #           Widget Event Handlers  
    #
    #=======================================================

    def on_coil_combine_method(self, event):
        value = event.GetEventObject().GetStringSelection()

        if value in ['External Dataset','External Dataset with Offset']:
            if self.block.set.coil_combine_external_dataset is None:
                if not self._coil_combine_external_dataset_browse():
                    event.GetEventObject().SetStringSelection(self.block.set.coil_combine_method)
                    return

        self.block.set.coil_combine_method = value
        
        self.top.Freeze()
        if value in ['External Dataset','External Dataset with Offset']:
            self.PanelCoilCombineExternalDataset.Show()
        else:
            self.PanelCoilCombineExternalDataset.Hide()
        self.top.Layout()
        self.PanelPrepFidsum.Layout()
        self.top.Thaw()
        
        self.process_and_plot()


    def on_apply_noise_whitening(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.coil_combine_noise_whitening = value
        self.process_and_plot()

    def on_coil_combine_external_dataset_browse(self, event):
            if self._coil_combine_external_dataset_browse():
                self.process_and_plot()

    def _coil_combine_external_dataset_browse(self):
        # Allows the user to select an ECC dataset.
        dialog = dialog_dataset_browser.DialogDatasetBrowser(self.top.datasets)
        dialog.ShowModal()
        the_dataset = dialog.dataset
        dialog.Destroy()
        if the_dataset:
            self.dataset.set_associated_dataset_combine(the_dataset)
            self.TextCoilCombineExternalDataset.SetValue(the_dataset.blocks["raw"].data_source)
            return True
        else:
            return False
        
    def on_apply_data_exclusion(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.apply_data_exclusion = value
        self.process_and_plot()

    def on_data_exclusion_plot_display(self, event):
        self.plot()

    def on_exclusion_method(self, event):
        value = event.GetEventObject().GetStringSelection()

        self.block.set.exclude_method = value

        self.top.Freeze()
        if value == 'Manual':
            self.PanelExclusionManualOptions.Show()
            self.PanelExclusionFidaRmBad.Hide()
            self.PanelExclusionFidaWorstN.Hide()
        elif value == 'Remove Bad Averages (fid-a)':
            self.PanelExclusionManualOptions.Hide()
            self.PanelExclusionFidaRmBad.Show()
            self.PanelExclusionFidaWorstN.Hide()
        elif value == 'Remove N Worst Averages (fid-a)':
            self.PanelExclusionManualOptions.Hide()
            self.PanelExclusionFidaRmBad.Hide()
            self.PanelExclusionFidaWorstN.Show()
        self.top.Layout()
        self.PanelPrepFidsum.Layout()
        self.top.Thaw()
        self.process_and_plot()

    def on_exclusions_input_adjust(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.exclusion_input_adjust = value
        self.process_and_plot()

    def on_data_exclusion_indices(self, event):
        pass

    def on_toggle_current_index(self, event):
        self.block.toggle_exclude_index(self.dataset, self.fid_index)
        val = ','.join(str(x) for x in self.block.exclude_indices)
        self.TextDataExclusion.SetValue(val)
        self.process_and_plot()
    
    def on_clear_indices(self, event):
        self.block.toggle_exclude_index(self.dataset, None)     # none indicates a reset
        self.TextDataExclusion.SetValue('')
        self.process_and_plot()

    def on_exclusion_auto_fida_bad_threshold(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.fida_bad_threshold = value
        self.process_and_plot()

    def on_exclusion_auto_fida_bad_domain(self, event):
        val = event.GetEventObject().GetStringSelection()
        value = 't' if val == 'Time' else 'f'
        self.block.set.fida_bad_domain = value
        self.process_and_plot()

    def on_exclusion_auto_fida_worst_n(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.fida_n_worst = value
        self.process_and_plot()

    def on_correction_method(self, event):
        value = event.GetEventObject().GetStringSelection()

        self.block.set.correction_method = value

        self.top.Freeze()
        if value == 'Manual':
            self.PanelCorrectionsVespa.Hide()
            self.PanelCorrectionsVespaCorrelate.Hide()
            self.PanelCorrectionsSuspectSpectralRegistration.Hide()
            self.PanelCorrectionsSuspectRats.Hide()
        elif value == 'Optimized Search (vespa)':
            self.PanelCorrectionsVespa.Show()
            self.PanelCorrectionsVespaCorrelate.Hide()
            self.PanelCorrectionsSuspectSpectralRegistration.Hide()
            self.PanelCorrectionsSuspectRats.Hide()
        elif value == 'Correlation (vespa)':
            self.PanelCorrectionsVespa.Hide()
            self.PanelCorrectionsVespaCorrelate.Show()
            self.PanelCorrectionsSuspectSpectralRegistration.Hide()
            self.PanelCorrectionsSuspectRats.Hide()
        elif value == 'Spectral Registration (suspect)':
            self.PanelCorrectionsVespa.Hide()
            self.PanelCorrectionsVespaCorrelate.Hide()
            self.PanelCorrectionsSuspectSpectralRegistration.Show()
            self.PanelCorrectionsSuspectRats.Hide()
        elif value == 'RATS (suspect)':
            self.PanelCorrectionsVespa.Hide()
            self.PanelCorrectionsVespaCorrelate.Hide()
            self.PanelCorrectionsSuspectSpectralRegistration.Hide()
            self.PanelCorrectionsSuspectRats.Show()
        self.top.Layout()
        self.PanelPrepFidsum.Layout()
        self.top.Thaw()

        if self.block.set.auto_correct:
            self.process_and_plot()

    def on_calculate_corrections(self, event):
        self.process_and_plot()

    def on_auto_calculate(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.auto_correct = value
        if value:
            self.process_and_plot()

    def on_correction_input_adjust(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.correction_input_adjust = value
        if self.block.set.auto_correct:
            self.process_and_plot()

    def on_vespa_reference_peak_center(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.vespa_reference_peak_center = value
        if self.block.set.auto_correct:
            self.process_and_plot()

    def on_vespa_peak_search_width(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.vespa_peak_search_width = value
        if self.block.set.auto_correct:
            self.process_and_plot()

    def on_vespa_phase0_range_start(self, event):
        # Note. min=End and max=Start because dealing with PPM range
        min, max = _paired_event(self.FloatVespaPhase0RangeEnd,
                                 self.FloatVespaPhase0RangeStart)
        self.block.set.vespa_phase0_range_start = max
        self.block.set.vespa_phase0_range_end   = min
        if self.block.set.auto_correct:
            self.process_and_plot()

    def on_vespa_phase0_range_end(self, event):
        # Note. min=End and max=Start because dealing with PPM range
        min, max = _paired_event(self.FloatVespaPhase0RangeEnd,
                                 self.FloatVespaPhase0RangeStart)
        self.block.set.vespa_phase0_range_start = max
        self.block.set.vespa_phase0_range_end   = min
        if self.block.set.auto_correct:
            self.process_and_plot()

    def on_vespa_target_spectrum(self, event):
        value = event.GetEventObject().GetStringSelection()
        self.block.set.vespa_target_method = value
        if self.block.set.auto_correct:
            self.process_and_plot()

    def on_vespa_correlate_prior(self, event):
        prior = self.block.set.vespa_preprocess_prior
        dialog = DialogUserPrior(self, self.dataset, prior, show_ph1=False)
        if dialog.ShowModal():
            # Changes were made
            self.block.set.vespa_preprocess_prior = dialog.prior
            self.block.set.vespa_preprocess_prior.basis.update(self.dataset, zfmult=4)
            if self.block.set.auto_correct:
                self.process_and_plot()

    def on_suspect_initial_guess_freq(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.suspect_initial_guess_freq = value
        if self.block.set.auto_correct:
            self.process_and_plot()

    def on_suspect_initial_guess_phase(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.suspect_initial_guess_phase = value
        if self.block.set.auto_correct:
            self.process_and_plot()

    def on_suspect_optimize_range_end(self, event):
        # Note. min=End and max=Start because dealing with PPM range
        min, max = _paired_event(self.FloatSuspectOptimizeRangeEnd,
                                 self.FloatSuspectOptimizeRangeStart)
        self.block.set.suspect_optimization_range_start = max
        self.block.set.suspect_optimization_range_end   = min
        if self.block.set.auto_correct:
            self.process_and_plot()

    def on_suspect_optimize_range_start(self, event):
        # Note. min=End and max=Start because dealing with PPM range
        min, max = _paired_event(self.FloatSuspectOptimizeRangeEnd,
                                 self.FloatSuspectOptimizeRangeStart)
        self.block.set.suspect_optimization_range_start = max
        self.block.set.suspect_optimization_range_end   = min
        if self.block.set.auto_correct:
            self.process_and_plot()

    def on_suspect_target_spectrum(self, event):
        value = event.GetEventObject().GetStringSelection()
        self.block.set.suspect_target_method = value
        if self.block.set.auto_correct:
            self.process_and_plot()

    def on_rats_initial_guess_freq(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.rats_initial_guess_freq = value
        if self.block.set.auto_correct:
            self.process_and_plot()

    def on_rats_initial_guess_phase(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.rats_initial_guess_phase = value
        if self.block.set.auto_correct:
            self.process_and_plot()

    def on_rats_optimize_range_start(self, event):
        # Note. min=End and max=Start because dealing with PPM range
        min, max = _paired_event(self.FloatRatsOptimizeRangeEnd,
                                 self.FloatRatsOptimizeRangeStart)
        self.block.set.rats_optimization_range_start = max
        self.block.set.rats_optimization_range_end   = min
        if self.block.set.auto_correct:
            self.process_and_plot()

    def on_rats_optimize_range_end(self, event):
        # Note. min=End and max=Start because dealing with PPM range
        min, max = _paired_event(self.FloatRatsOptimizeRangeEnd,
                                 self.FloatRatsOptimizeRangeStart)
        self.block.set.rats_optimization_range_start = max
        self.block.set.rats_optimization_range_end   = min
        if self.block.set.auto_correct:
            self.process_and_plot()

    def on_rats_baseline_order(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.rats_baseline_order = value
        if self.block.set.auto_correct:
            self.process_and_plot()

    def on_rats_target_spectrum(self, event):
        value = event.GetEventObject().GetStringSelection()
        self.block.set.rats_target_method = value
        if self.block.set.auto_correct:
            self.process_and_plot()

    def on_fid_index(self, event):
        value = event.GetEventObject().GetValue()
        self.fid_index = value # - 1
        self.process_and_plot(entry='adjust')

    def on_current_peak_shift_value(self, event):
        value = event.GetEventObject().GetValue()
        self.block.frequency_shift[self.fid_index] = value
        self.process_and_plot(entry='correct_adjust')
        
    def on_current_phase0_value(self, event):
        value = event.GetEventObject().GetValue()
        self.block.phase_0[self.fid_index] = value
        self.process_and_plot(entry='correct_adjust')

    def on_waterfall_range_start(self, event):
        # Note. min=End and max=Start because dealing with PPM range
        min, max = _paired_event(self.FloatWaterfallRangeEnd,
                                 self.FloatWaterfallRangeStart)
        self.waterfall_range_start = max
        self.waterfall_range_end = min
        self.process_and_plot(entry='adjust')

    def on_waterfall_range_end(self, event):
        # Note. min=End and max=Start because dealing with PPM range
        min, max = _paired_event(self.FloatWaterfallRangeEnd,
                                 self.FloatWaterfallRangeStart)
        self.waterfall_range_start = max
        self.waterfall_range_end = min
        self.process_and_plot(entry='adjust')

    def on_global_left_shift(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.global_left_shift = value
        self._update_freq_raw = True
        set = self.block.set
        if set.exclusion_input_adjust or set.correction_input_adjust:
            entry = 'all'
        else:
            entry = 'correct_adjust'
        self.process_and_plot(entry=entry)

    def on_global_phase0(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.global_phase0 = value
        self._update_freq_raw = True
        set = self.block.set
        if set.exclusion_input_adjust or set.correction_input_adjust:
            entry = 'all'
        else:
            entry = 'correct_adjust'
        self.process_and_plot(entry=entry)

    def on_global_phase1(self, event):
        value = event.GetEventObject().GetValue()
        entry = 'correct_adjust'
        set = self.block.set
        if not set.zero_phase1:
            set.global_phase1 = value
            self._update_freq_raw = True
            if set.exclusion_input_adjust or set.correction_input_adjust:
                entry = 'all'
        else:
            self.block.set.global_phase1 = 0.0
            self.FloatGlobalPhase1.SetValue(0.0)
        self.process_and_plot(entry=entry)

    def on_global_gauss_apodize(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.global_gaussian_apodization = value
        self._update_freq_raw = True
        set = self.block.set
        if set.exclusion_input_adjust or set.correction_input_adjust:
            entry = 'all'
        else:
            entry = 'correct_adjust'
        self.process_and_plot(entry=entry)

    def on_chop(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.chop_data = value
        self._update_freq_raw = True
        self.process_and_plot()

    def on_zero_global_phase1(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.zero_phase1 = value
        self._update_freq_raw = True
        if value:
            self.block.set.global_phase1 = 0.0
            self.FloatGlobalPhase1.SetValue(0.0)
        self.process_and_plot(entry='correct_adjust')

    def on_apply_peak_shift(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.apply_peak_shift = value
        self.process_and_plot(entry='correct_adjust')

    def on_reset_peak_shift(self, event): 
        self.block.frequency_shift *= 0
        self.FloatCurrentPeakShiftValue.SetValue(self.block.frequency_shift[self.fid_index])
        self.process_and_plot(entry='correct_adjust')

    def on_apply_phase0(self, event):
        value = event.GetEventObject().GetValue()
        self.block.set.apply_phase0 = value
        self.process_and_plot(entry='correct_adjust')

    def on_reset_phase0(self, event):
        self.block.phase_0 *= 0
        self.FloatCurrentPhase0Value.SetValue(self.block.phase_0[self.fid_index])
        self.process_and_plot(entry='correct_adjust')

    def on_push_results(self, event):
    
        raw = self.dataset.blocks['raw']
        if isinstance(raw, block_raw_edit_fidsum.BlockRawEditFidsum):
            shifts = self.block.frequency_shift
            phase0 = self.block.phase_0
            phase1 = self.block.set.global_phase1
            
            uids = [raw.data_on.id, raw.data_off.id, raw.data_sum.id, raw.data_dif.id]

            pubsub.sendMessage("push_prep_fidsum_results", 
                               uids=uids, 
                               shifts=shifts, 
                               phase0=phase0, 
                               phase1=phase1)


    def on_splitter(self, event=None):
        # This is sometimes called programmatically, in which case event is None
        self._prefs.sash_position = self.SplitterWindow.GetSashPosition()

    def on_splitter_leftright(self, event):
        pass
        #wx.CallLater(200, self.plot)

    def on_splitter_topbottom(self, event):
        pass
        #wx.CallLater(200, self.plot)

    #=======================================================
    #
    #           Public Methods
    #
    #=======================================================

    def process_and_plot(self, entry='all', plot_entry=''):
        """ process(), plot() and process_and_plot() are standard in all tabs """
        tab_base.Tab.process_and_plot(self, entry)
        self.process(entry=entry)
        self.plot(plot_entry=plot_entry)


    def process(self, entry='all'):
        """
        The Chain stores processing results into the Block automatically in
        chain.run(). View results are returned as a dictionary from chain.run().
        The plot() method takes its inputs from this dictionary.
        
        """
        tab_base.Tab.process(self, entry)

        voxel = self.fid_index
        self.plot_results = self.block.chain.run([voxel], entry=entry, freq_raw=self._update_freq_raw)

        # update exclusion 'TextList' and 'PlotMarks'
        nexclude = len(self.block.exclude_indices)
        nfids = self.block.chain.raw.shape[2] - nexclude
        val = '' if nexclude == 0 else ','.join(str(x) for x in self.block.exclude_indices)
        self.TextDataExclusion.SetValue(val)
        msg = 'Indices Excluded from Result    -    FIDs Remaining = %s,  Excluded = %s' % (nfids, nexclude)
        self.LabelRemainingFidCount.SetLabel(msg)

        # update correction 'ResultsControls'
        self.FloatCurrentPeakShiftValue.SetValue(self.block.frequency_shift[self.fid_index])
        self.FloatCurrentPhase0Value.SetValue(self.block.phase_0[self.fid_index])

        # reset control flags
        self._update_freq_raw = False

        
    def plot(self, plot_entry=''):
        """
        set_data() puts results into the PlotPanel object in the right panel.

        """
        if self._plotting_enabled:
            tab_base.Tab.plot(self)

            results = self.plot_results

            data1 = np.squeeze(results['freq_current'])
            data2 = np.squeeze(results['freq_summed'])
            
            data = [[data1], [data2]]
            self.view.set_data(data)
            self.view.update(set_scale=not self._scale_intialized)
            
            select = self.ComboDataExclusionPlotDisplay.GetStringSelection()
            if select == EXCLUDE_DISPLAY_CHOICES[0]:
                data3  = np.abs(self.block.chain.raw[0,0,:,0])
                do_scale = (self.view_series.xlabel != 'fid[0] Magn')
                self.view_series.xlabel = 'fid[0] Magn'
            elif select == EXCLUDE_DISPLAY_CHOICES[1]:
                data3  = self.block.frequency_shift.copy()
                do_scale = (self.view_series.xlabel != 'B0 Shift [Hz]')
                self.view_series.xlabel = 'B0 Shift [Hz]'
            elif select == EXCLUDE_DISPLAY_CHOICES[2]:
                data3  = self.block.phase_0.copy()
                do_scale = (self.view_series.xlabel != 'Phase 0 [deg]')
                self.view_series.xlabel = 'Phase 0 [deg]'

            exclude = self.block.exclude_indices
            data3 = {'data':data3, 'markevery':exclude, 'markevery_color':'red'}

            data_series = [[data3],]
            self.view_series.set_data(data_series)
            self.view_series.update(set_scale=do_scale)

            if plot_entry != 'no_image':

                # do view_image here
                istr = np.int(self.dataset.ppm2pts(self.waterfall_range_start, acq=True))
                iend = np.int(self.dataset.ppm2pts(self.waterfall_range_end, acq=True))

                img0 = np.squeeze(results['freq_adjusted'])
                if len(img0.shape) == 1:
                    img0.shape = 1,img0.shape[0]

                img0 = img0[:,istr:iend]
                img0 = np.transpose(img0.real)
                navg = img0.shape[1]
                vmax0 = img0.max() - 0.3*(img0.max() - img0.min())
                vmin0 = img0.min() + 0.1*(img0.max() - img0.min())

                img1 = np.squeeze(results['freq_raw'])
                if len(img1.shape) == 1:
                    img1.shape = 1,img1.shape[0]

                img1 = img1[:,istr:iend]
                img1 = np.transpose(img1.real)
                vmax1 = img1.max() - 0.3*(img1.max() - img1.min())
                vmin1 = img1.min() + 0.1*(img1.max() - img1.min())

                navg = img1.shape[1]
                xstr, xend = 0, navg - 1
                if navg <= 1:
                    xstr, xend = -0.05,  0.05

                txt = 'FIDs Before (bottom) and After (top) Alignment'
                self.view_image.figure.suptitle(txt, fontsize='medium')
                ax0, ax1 = self.view_image.axes[0:2]
                x0,   x1 = ax0.get_xlim(), ax1.get_xlim()
                ax0.clear()
                ax1.clear()
                yrange = [self.waterfall_range_end, self.waterfall_range_start]

                ax0.get_xaxis().set_ticks([])

                ax0.imshow(img0, cmap=cm.jet, origin='upper',
                           aspect='auto', vmin=vmin0, vmax=vmax0,
                           extent=[xstr, xend, self.waterfall_range_end, self.waterfall_range_start])
                ax1.imshow(img1, cmap=cm.jet, origin='upper',
                           aspect='auto', vmin=vmin1, vmax=vmax1,
                           extent=[xstr, xend, self.waterfall_range_end, self.waterfall_range_start])

                ax0.set_xlim(x0)
                ax0.set_ylim(yrange)
                ax1.set_xlim(x1)
                ax1.set_ylim(yrange)

                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    self.view_image.canvas.draw()

            if not self._scale_intialized:
                self._scale_intialized = True

            # Calculate the new area after phasing
            area, rms = self.view.calculate_area()
            index = (0 if self._prefs.area_calc_plot_a else 1)
            area = area[index]
            rms = rms[index]

            self.top.statusbar.SetStatusText(self.build_area_text(area, rms), 3)


    #=======================================================
    #
    #           Internal Helper Functions  
    #
    #=======================================================
    

    def _display_header_text(self):
    
        index = self.fid_index
        header = self.dataset.blocks["raw"].headers[index]
        lines = "\nCurrent Header, FID index = "+str(index)+"\n" + "-" * 75 + "\n\n"
        lines += str(header)
        wx_util.display_text_as_file(lines)

