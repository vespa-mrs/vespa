# Python modules

# 3rd party modules
import wx

# Our modules
import vespa.analysis.tab_base as tab_base
import vespa.analysis.prefs as prefs_module
import vespa.analysis.constants as constants
import vespa.analysis.util_menu as util_menu
import vespa.analysis.auto_gui.timeseries as timeseries
import vespa.analysis.functors.funct_fidsum_coil_combine as funct_coil_combine
import vespa.common.wx_gravy.util as wx_util
import vespa.common.wx_gravy.common_dialogs as common_dialogs

from vespa.analysis.plot_panel_prep_timeseries import PlotPanelPrepTimeseries


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
#
#  Tab RAW FIDSUM 
#
#------------------------------------------------------------------------------

class TabPrepTimeseries(tab_base.Tab, timeseries.PanelPrepTimeseriesUI):
    
    # self-identify tab to notebook, value does not matter, its presence is sufficient.
    IS_PREP_TIMESERIES = True
    
    def __init__(self, tab_dataset, top, block):
        
        timeseries.PanelPrepTimeseriesUI.__init__(self, tab_dataset.NotebookDataset)
        
        tab_base.Tab.__init__(self, tab_dataset, top, prefs_module.PrefsPrepFidsum)
        
        self.top      = top               # application frame
        self.block    = block             # processing object

        # Plotting is disabled during some of init. That's because the plot
        # isn't ready to plot, but the population of some controls 
        # (e.g. spin controls on the water filter panel) fires their 
        # respective change event which triggers a call to plot(). This
        # appears to happen only under Windows; it might be a Windows-specific
        # bug.
        # In any case, skipping some calls to plot() will speed things up. =)
        self._plotting_enabled = False
        self.plot_results = None
        
        self.fid_index      = 0
        self.initialize_controls()
        self.populate_controls()

        self._plotting_enabled = True 

        #------------------------------------------------------------
        # Setup the canvas 
        self.process_and_plot()

        #------------------------------------------------------------
        # PubSub subscriptions
#        pubsub.subscribe(self.on_push_prep_fidsum_results, "push_prep_fidsum_results")


        # If the sash position isn't recorded in the INI file, we use the
        # arbitrary-ish value of 400.
        if not self._prefs.sash_position:
            self._prefs.sash_position = 400

        # Under OS X, wx sets the sash position to 10 (why 10?) *after*
        # this method is done. So setting the sash position here does no
        # good. We use wx.CallAfter() to (a) set the sash position and
        # (b) fake an EVT_SPLITTER_SASH_POS_CHANGED.
        wx.CallAfter(self.SplitterWindow.SetSashPosition, 
                     self._prefs.sash_position, True)
        wx.CallAfter(self.on_splitter)


    #=======================================================
    #
    #           GUI Setup Handlers 
    #
    #=======================================================

    def initialize_controls(self):
        """ 
        Initializes the controls to be the right size or have the right
        range or number of decimal places. It typically does not set the
        default value (that's for populate_controls method to do). This
        method does the one-time setup bits.
        
        """
        dataset = self.dataset
        dim0, dim1, dim2, dim3 = dataset.spectral_dims
        sw      = dataset.sw
        maxppm  = dataset.pts2ppm(0)
        minppm  = dataset.pts2ppm(dim0-1)
        ppmlim  = (minppm, maxppm)
        
        
        wx_util.configure_spin(self.SpinFidIndex, 60)
        wx_util.configure_spin(self.SpinFidsToAverage, 60)

        wx_util.configure_spin(self.FloatGaussianApodization, 70, 2, constants.PrepFidsum.STEP_APODIZE,
                       (constants.PrepFidsum.MIN_APODIZE, constants.PrepFidsum.MAX_APODIZE))

        wx_util.configure_spin(self.SpinFidLeftShift, 70, None, None,
                       (constants.PrepFidsum.MIN_LEFT, constants.PrepFidsum.MAX_LEFT))

        wx_util.configure_spin(self.FloatPeakShiftValue, 70, 3, constants.PrepFidsum.STEP_PEAK,
                       (constants.PrepFidsum.MIN_PEAK,constants.PrepFidsum.MAX_PEAK))

        wx_util.configure_spin(self.FloatPhase0, 70, 3, constants.PrepFidsum.STEP_PHASE,
                       (constants.PrepFidsum.MIN_PHASE,constants.PrepFidsum.MAX_PHASE))

        wx_util.configure_spin(self.FloatGlobalPhase1, 70, 3, constants.PrepFidsum.STEP_PHASE1,
                       (constants.PrepFidsum.MIN_PHASE1,constants.PrepFidsum.MAX_PHASE1))

        wx_util.configure_spin(self.FloatGlobalPhase0, 70, 3, constants.PrepFidsum.STEP_PHASE,
                       (constants.PrepFidsum.MIN_PHASE,constants.PrepFidsum.MAX_PHASE))

        wx_util.configure_spin(self.FloatReferencePeakCenter, 70, 3, constants.PrepFidsum.STEP_CENTER,
                       (constants.PrepFidsum.MIN_CENTER,constants.PrepFidsum.MAX_CENTER))

        wx_util.configure_spin(self.FloatPeakSearchWidth, 70, 3, constants.PrepFidsum.STEP_WIDTH,
                       (constants.PrepFidsum.MIN_WIDTH,constants.PrepFidsum.MAX_WIDTH))

        wx_util.configure_spin(self.FloatPhase0RangeStart, 70, 2, 0.25, ppmlim)
        wx_util.configure_spin(self.FloatPhase0RangeEnd,   70, 2, 0.25, ppmlim)


        #-------------------------------------------------------------
        # Raw Fidsum View setup 

        self.view = PlotPanelPrepTimeseries(self.PanelViewPrepFidsum, 
                                            self,
                                            self._tab_dataset,
                                            naxes=2,
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

        width, height = self.SpinFidLeftShift.GetSize()

        if "__WXMAC__" in wx.PlatformInfo:
            # Under OS X this spin control needs a poke to paint itself 
            # properly. Without this code, the control doesn't appear on screen
            # even though the wx Inspector reports that it has an appropriate
            # size & position. Go figger.
            wx.CallAfter(self.SpinFidLeftShift.SetSize, (width + 1, height))

  
    def populate_controls(self, preset=False):
        """ 
        Populates the raw data tab with values from the dataset.raw 
        object. It's meant to be called when a new data object is loaded.

        """
        block = self.block

        raw = self._tab_dataset.dataset.get_source_data('prep')

        self.ComboCoilCombineMethod.Clear()
        for item in funct_coil_combine.COILCOMBINE_MENU_ITEMS:
            self.ComboCoilCombineMethod.Append(item, "")
        if self.block.set.coil_combine_method == '':
            self.block.set.coil_combine_method = 'None'
        self.ComboCoilCombineMethod.SetStringSelection(self.block.set.coil_combine_method)

        # only turn on combine method widget if data has indiv coil FIDs in it        
        if raw.shape[1] > 1:
            self.PanelCoilCombination.Enable(True)
        else:
            self.PanelCoilCombination.Enable(False)

        imax = int(raw.shape[-2]/self.block.set.fids_to_average)        # depends on if averaging fids or not

        self.SpinFidIndex.SetValue(self.fid_index+1)
        self.SpinFidIndex.SetRange(1, imax)

        self.SpinFidsToAverage.SetValue(self.block.set.fids_to_average)
        self.SpinFidsToAverage.SetRange(1, raw.shape[-2])

        self.FloatGaussianApodization.SetValue(self.block.set.gaussian_apodization)
        self.SpinFidLeftShift.SetValue(self.block.set.fid_left_shift)
        
        self.FloatPeakShiftValue.SetValue(self.block.frequency_shift[self.fid_index])
        self.FloatPhase0.SetValue(self.block.phase_0[self.fid_index])
        self.FloatGlobalPhase1.SetValue(self.block.set.global_phase1)
        self.FloatGlobalPhase0.SetValue(self.block.set.global_phase0)
        
        self.CheckApplyPeakShift.SetValue(self.block.set.apply_peak_shift)
        self.FloatReferencePeakCenter.SetValue(self.block.set.reference_peak_center)
        self.FloatPeakSearchWidth.SetValue(self.block.set.peak_search_width)
        self.CheckApplyPhase0.SetValue(self.block.set.apply_phase0)
        self.FloatPhase0RangeStart.SetValue(self.block.set.phase0_range_start)
        self.FloatPhase0RangeEnd.SetValue(self.block.set.phase0_range_end)
        

    
    #=======================================================
    #
    #           Global and Menu Event Handlers 
    #
    #=======================================================

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

#    def on_push_prep_fidsum_results(self, uids=None, shifts=None, phase0=None, phase1=None):
#        """ 
#        Pubsub notifications handler. 
#        
#        Calling tab sends list of all associated dataset ids, and new freq
#        and phase values.  If self.dataset.id in the uids list, we update the
#        results arrays, the widgets, and re-process and plot spectra
#        
#        """
#        if self.dataset.id in uids:
#        
#            block = self.dataset.blocks['prep']
#    
#            block.frequency_shift    = shifts.copy()
#            block.phase_0            = phase0.copy()
#            block.set.global_phase1  = phase1
#
#            self.FloatPeakShiftValue.SetValue(block.frequency_shift[self.fid_index])
#            self.FloatPhase0.SetValue(block.phase_0[self.fid_index])
#            self.FloatPhase1.SetValue(block.set.global_phase1)
#        
#            self.process_and_plot()


    #=======================================================
    #
    #           Widget Event Handlers  
    #
    #=======================================================

    def on_coil_combine_method(self, event):
        value = event.GetEventObject().GetStringSelection()
        self.block.set.coil_combine_method = value
        self.process_and_plot()

    def on_fids_to_average(self, event):
        value = event.GetEventObject().GetValue()
        if value != self.block.set.fids_to_average:
            self.block.set.fids_to_average = value
            #self.block._reset_dimensional_data(self.dataset)
            self.dataset.update_for_preprocess_change()

            nfids = self.block.data.shape[-2]
            val = self.SpinFidIndex.GetValue()
            if val >= nfids:
                self.SpinFidIndex.SetValue(nfids) 
                self.fid_index = nfids - 1
                shft = self.block.frequency_shift[self.fid_index]
                phas = self.block.phase_0[self.fid_index]
                self.FloatPeakShiftValue.SetValue(shft)        
                self.FloatPhase0.SetValue(phas)
            self.SpinFidIndex.SetRange(1, nfids)          
            self._tab_dataset.set_voxel_range(nfids)
            
            self.process(entry='all', do_calculate=False)

            self.plot()            

            
        
    def on_splitter(self, event=None):
        # This is sometimes called programmatically, in which case event is None
        self._prefs.sash_position = self.SplitterWindow.GetSashPosition()

    def on_fid_index(self, event): 
        value = event.GetEventObject().GetValue()
        self.fid_index = value - 1
        shft = self.block.frequency_shift[self.fid_index]
        phas = self.block.phase_0[self.fid_index]
        self.FloatPeakShiftValue.SetValue(shft)        
        self.FloatPhase0.SetValue(phas)
        self.FloatGlobalPhase0.SetValue(self.block.set.global_phase0)
        self.FloatGlobalPhase1.SetValue(self.block.set.global_phase1)
        self.process_and_plot()

    def on_gaussian_apodization(self, event): 
        value = event.GetEventObject().GetValue()
        self.block.set.gaussian_apodization = value
        self.process_and_plot()

    def on_fid_left_shift(self, event): 
        value = event.GetEventObject().GetValue()
        self.block.set.fid_left_shift = value
        self.process_and_plot()

    def on_peak_shift_value(self, event): 
        value = event.GetEventObject().GetValue()
        self.block.frequency_shift[self.fid_index] = value
        self.process_and_plot()
        
    def on_phase0(self, event): 
        value = event.GetEventObject().GetValue()
        self.block.phase_0[self.fid_index] = value
        self.process_and_plot()

    def on_global_phase1(self, event): 
        value = event.GetEventObject().GetValue()
        self.block.set.global_phase1 = value
        self.process_and_plot()

    def on_global_phase0(self, event): 
        value = event.GetEventObject().GetValue()
        self.block.set.global_phase0 = value
        self.process_and_plot()
        
    def on_apply_peak_shift(self, event): 
        value = event.GetEventObject().GetValue()
        self.block.set.apply_peak_shift = value

    def on_reset_peak_shift(self, event): 
        self.block.frequency_shift *= 0
        self.FloatPeakShiftValue.SetValue(self.block.frequency_shift[self.fid_index])
        self.process_and_plot()

    def on_reference_peak_center(self, event): 
        value = event.GetEventObject().GetValue()
        self.block.set.reference_peak_center = value

    def on_peak_search_width(self, event): 
        value = event.GetEventObject().GetValue()
        self.block.set.peak_search_width = value

    def on_apply_phase0(self, event): 
        value = event.GetEventObject().GetValue()
        self.block.set.apply_phase0 = value

    def on_reset_phase0(self, event): 
        self.block.phase_0 *= 0
        self.FloatPhase0.SetValue(self.block.phase_0[self.fid_index])
        self.process_and_plot()

    def on_phase0_range_start(self, event): 
        # Note. min=End and max=Start because dealing with PPM range
        min, max = _paired_event(self.FloatPhase0RangeEnd,
                                 self.FloatPhase0RangeStart)
        self.block.set.phase0_range_start = max 
        self.block.set.phase0_range_end   = min

    def on_phase0_range_end(self, event): 
        # Note. min=End and max=Start because dealing with PPM range
        min, max = _paired_event(self.FloatPhase0RangeEnd,
                                 self.FloatPhase0RangeStart)
        self.block.set.phase0_range_start = max
        self.block.set.phase0_range_end   = min

    def on_calculate_corrections(self, event): 
        self.process_and_plot(do_calculate=True)    
        shft = self.block.frequency_shift[self.fid_index]
        phas = self.block.phase_0[self.fid_index]
        self.FloatPeakShiftValue.SetValue(shft)        
        self.FloatPhase0.SetValue(phas)

    def on_push_results(self, event):
        pass
#        raw = self.dataset.blocks['raw']
#        if isinstance(raw, block_raw_edit_fidsum.BlockRawEditFidsum):
#            shifts = self.block.frequency_shift
#            phase0 = self.block.phase_0
#            phase1 = self.block.set.global_phase1
#            
#            uids = [raw.data_on.id, raw.data_off.id, raw.data_sum.id, raw.data_dif.id]
#
#            pubsub.sendMessage("push_prep_fidsum_results", 
#                               uids=uids, 
#                               shifts=shifts, 
#                               phase0=phase0, 
#                               phase1=phase1)


    #=======================================================
    #
    #           Public Methods
    #
    #=======================================================

    def process_and_plot(self, entry='all', do_calculate=False): 
        """
        The process(), plot() and process_and_plot() methods are standard in
        all processing tabs. They are called to update the data in the plot
        results dictionary, the plot_panel in the View side of the tab or both.

        """
        tab_base.Tab.process_and_plot(self, entry)

        self.process(entry=entry, do_calculate=do_calculate)
        self.plot()


    def process(self, entry='all', do_calculate=False):
        """
        Data processing results are stored into the Block inside the Chain,
        but the View results are returned as a dictionary from the Chain.run()
        method. The plot routine takes its inputs from this dictionary.
        
        """
        tab_base.Tab.process(self, entry)

        voxel = self.fid_index
        self.plot_results = self.block.chain.run([voxel], 
                                                 entry=entry, 
                                                 do_calculate=do_calculate)
        
        
    def plot(self):
        """
        The set_data() method sets data into the plot_panel_spectrum object
        in the plot in the right panel. 

        """
        if self._plotting_enabled:
            tab_base.Tab.plot(self)

            results = self.plot_results

            data1 = results['freq_current']
            data2 = results['freq_all']

            data = [[data1], [data2]]
            self.view.set_data(data)
            self.view.update(set_scale=not self._scale_intialized)

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

