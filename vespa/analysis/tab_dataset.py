# Python modules
import math
from collections import OrderedDict

# 3rd party modules
import wx
#import wx.aui as aui
import wx.lib.agw.aui as aui        # NB. wx.aui version throws odd wxWidgets exception on Close/Exit  ?? Not anymore in wxPython 4.0.6 ??

# Our modules
import vespa.analysis.util_menu as util_menu
import vespa.analysis.block_raw as block_raw
import vespa.analysis.block_raw_probep as block_raw_probep
import vespa.analysis.block_raw_edit as block_raw_edit
import vespa.analysis.block_raw_edit_fidsum as block_raw_edit_fidsum
import vespa.analysis.block_prep_fidsum as block_prep_fidsum
import vespa.analysis.block_prep_timeseries as block_prep_timeseries
import vespa.analysis.block_prep_wbnaa as block_prep_wbnaa
import vespa.analysis.block_prep_edit_fidsum as block_prep_edit_fidsum
import vespa.analysis.block_spectral as block_spectral
import vespa.analysis.block_fit_voigt as block_fit_voigt
import vespa.analysis.block_fit_giso as block_fit_giso
import vespa.analysis.block_quant_watref as block_quant_watref
import vespa.analysis.constants as constants
import vespa.analysis.tab_raw as tab_raw
import vespa.analysis.tab_prep_fidsum as tab_prep_fidsum
import vespa.analysis.tab_prep_timeseries as tab_prep_timeseries
import vespa.analysis.tab_prep_wbnaa as tab_prep_wbnaa
import vespa.analysis.tab_spectral as tab_spectral
import vespa.analysis.tab_voigt as tab_voigt
import vespa.analysis.tab_giso as tab_giso
import vespa.analysis.tab_watref as tab_watref
import vespa.analysis.dialog_user_metabolite_info as dialog_user_metabolite_info
import vespa.analysis.auto_gui.dataset as dataset_module
import vespa.common.wx_gravy.notebooks as notebooks

from vespa.analysis.dialog_user_prior import DialogUserPrior

PI = math.pi


#------------------------------------------------------------------------------
#
#  Tab DATASET
#
#------------------------------------------------------------------------------

class TabDataset(dataset_module.DatasetUI):

    def __init__(self, _outer_notebook, top, name):

        dataset_module.DatasetUI.__init__(self, _outer_notebook)

        self.top              = top
        self._outer_notebook  = _outer_notebook
        self.indexAB          = [name, None]  # indices for displayed datasets

        # _tabs is akin to Dataset.blocks. One tab per block. A slot can
        # contain None, no tab in that slot, or a Tab object.
        self._tabs = OrderedDict([["raw",None], ["prep",None], ["spectral",None], ["fit",None], ["quant",None]])

        # Determine whether certain options synchronized between datasets
        self.sync = False

        # Need a flag to let me know if we don't need to update
        self.turn_off_update = False

        # Create dataset notebook here (not wxGlade) to use our custom class
        sizer = self.LabelNotebookPlaceholder.GetContainingSizer()
        self.LabelNotebookPlaceholder.Destroy()
        style     = wx.BORDER_NONE
        agw_style = aui.AUI_NB_SCROLL_BUTTONS | aui.AUI_NB_BOTTOM          
        self.NotebookDataset = notebooks.VespaAuiNotebook(self, style, agw_style)
        sizer.Add(self.NotebookDataset, 1, wx.EXPAND, 0)

        self.Bind(aui.EVT_AUINOTEBOOK_PAGE_CHANGED, self.on_tab_changed, self.NotebookDataset)

        #------------------------------------------------------------
        # Init code - Dataset

        self.populate_controls()

        self.Bind(wx.EVT_WINDOW_DESTROY, self.on_destroy, self)


    @property
    def dataset(self):
        return self.top.datasets[self.indexAB[0]]

    @property
    def datasetB(self):
        return self.top.datasets[self.indexAB[1]]

    @property
    def voxel(self):
        return [self.SpinX.GetValue() - 1,
                self.SpinY.GetValue() - 1,
                self.SpinZ.GetValue() - 1]



    #=======================================================
    #
    #           GUI Setup Handlers
    #
    #=======================================================

    def populate_controls(self, preset=False):
        """
        Populates dataset tab with pages depending on the data object. This is
        called when a new data object is loaded.

        """
        dataset = self.dataset

        #-------------------------------------------------------------
        # Create notebook tabs based on dataset objects

        for block in dataset.blocks.values():
            if block.is_identity:
                # Identity blocks don't need a tab.
                pass
            else:
                # Match block/tab types. Block order determines GUI tab order

                if isinstance(block, block_raw.BlockRaw) or \
                   isinstance(block, block_raw_probep.BlockRawProbep) or \
                   isinstance(block, block_raw_edit.BlockRawEdit) or \
                   isinstance(block, block_raw_edit_fidsum.BlockRawEditFidsum):
                    tab = tab_raw.TabRaw(self, self.top, block)
                    self._tabs["raw"] = tab
                    self.NotebookDataset.AddPage(tab, "Raw", False)

                elif isinstance(block, block_prep_fidsum.BlockPrepFidsum) or \
                     isinstance(block, block_prep_edit_fidsum.BlockPrepEditFidsum):
                    tab = tab_prep_fidsum.TabPrepFidsum(self, self.top, block)
                    self._tabs["prep"] = tab
                    self.NotebookDataset.AddPage(tab, "Preprocess", False)

                elif isinstance(block, block_prep_timeseries.BlockPrepTimeseries):
                    tab = tab_prep_timeseries.TabPrepTimeseries(self, self.top, block)
                    self._tabs["prep"] = tab
                    self.NotebookDataset.AddPage(tab, "Preprocess", False)

                elif isinstance(block, block_prep_wbnaa.BlockPrepWbnaa):
                    tab = tab_prep_wbnaa.TabPrepWbnaa(self, self.top, block)
                    self._tabs["prep"] = tab
                    self.NotebookDataset.AddPage(tab, "Preprocess", False)

                elif isinstance(block, block_spectral.BlockSpectral):
                    tab = tab_spectral.TabSpectral(self, self.top, block)
                    self._tabs["spectral"] = tab
                    self.NotebookDataset.AddPage(tab, "Spectral", False)

                elif isinstance(block, block_fit_voigt.BlockFitVoigt):
                    tab = tab_voigt.TabVoigt(self, self.top, block)
                    self._tabs["fit"] = tab
                    self.NotebookDataset.AddPage(tab, "Fitting", False)

                elif isinstance(block, block_fit_giso.BlockFitGiso):
                    tab = tab_giso.TabGiso(self, self.top, block)
                    self._tabs["fit"] = tab
                    self.NotebookDataset.AddPage(tab, "Fitting", False)

                elif isinstance(block, block_quant_watref.BlockQuantWatref):
                    tab = tab_watref.TabWatref(self, self.top, block)
                    self._tabs["quant"] = tab
                    self.NotebookDataset.AddPage(tab, "Quant", False)

        #-------------------------------------------------------------
        # View Controls
        self.turn_off_update = True

        self.SpinX.SetRange(1, dataset.spectral_dims[1])
        self.SpinY.SetRange(1, dataset.spectral_dims[2])
        self.SpinZ.SetRange(1, dataset.spectral_dims[3])
        self.SpinX.SetValue(1)
        self.SpinY.SetValue(1)
        self.SpinZ.SetValue(1)

        self.turn_off_update = False

        self.on_voxel_change()

        if self._tabs["spectral"] is not None:
            self.FloatScale.SetRange(0.0, 1000000000.0)
            self.NotebookDataset.activate_tab(self._tabs["spectral"])
            self.FloatScale.multiplier = 1.1



    #=======================================================
    #
    #           Global and Menu Event Handlers
    #
    #=======================================================

    def on_add_voigt_tab(self, event):
        """
        Adds a Voigt tab if there isn't already a fitting tab, or activates
        the existing fitting tab if there is one.
        """
        # If there's already a voigt tab, just activate it.
        if self._tabs["fit"]:
            self.NotebookDataset.activate_tab(self._tabs["fit"])
        else:
            # Create a new voigt block and add a tab to go with it.
            block = self.dataset.add_voigt()

            # The fitting tab follows spectral
            tab = self._tabs["spectral"]
            index = self.NotebookDataset.get_tab_index(tab) + 1

            self.top.Freeze()

            # OK, now we know the insertion point. Create the new tab and
            # insert it.
            tab = tab_voigt.TabVoigt(self, self.top, block)
            self._tabs["fit"] = tab

            self.NotebookDataset.InsertPage(index, tab, 'Fitting', True)
            
            self.top.Thaw()

    def on_add_giso_tab(self, event):
        """
        Adds a Giso tab if there isn't already a fitting tab, or activates
        the existing fitting tab if there is one.
        """
        # If there's already a fit tab, just activate it.
        if self._tabs["fit"]:
            self.NotebookDataset.activate_tab(self._tabs["fit"])
        else:
            # Create a new giso block and add a tab to go with it.
            block = self.dataset.add_giso()

            # The fitting tab follows spectral
            tab = self._tabs["spectral"]
            index = self.NotebookDataset.get_tab_index(tab) + 1

            # OK, now we know the insertion point. Create the new tab and
            # insert it.
            tab = tab_giso.TabGiso(self, self.top, block)
            self._tabs["fit"] = tab

            self.NotebookDataset.InsertPage(index, tab, 'Fitting', True)            


    def on_add_watref_tab(self, event):
        """
        Adds a Water Reference Quantitation tab if there isn't already a quant
        tab, or activates the existing quant tab if there is one.
        """
        # If there's already a tab, just activate it.
        if self._tabs["quant"]:
            self.NotebookDataset.activate_tab(self._tabs["quant"])
        else:
            if self._tabs["fit"]:
                # Create a new quant block and add a tab to go with it.
                block = self.dataset.add_watref()
    
                # The quant tab follows fitting
                tab = self._tabs["fit"]
                index = self.NotebookDataset.get_tab_index(tab) + 1
    
                # We know the insertion point. Create and insert new tab.
                tab = tab_watref.TabWatref(self, self.top, block)
                self._tabs["quant"] = tab
    
                self.NotebookDataset.InsertPage(index, tab, 'Quant', True)


    def on_destroy(self, event):
        if self.indexAB[0] in self.top.datasets:
            del(self.top.datasets[self.indexAB[0]])


    def on_menu_view_option(self, event):
        self.NotebookDataset.active_tab.on_menu_view_option(event)

    def on_menu_view_output(self, event):
        self.NotebookDataset.active_tab.on_menu_view_output(event)

    def on_menu_view_results(self, event):
        self.NotebookDataset.active_tab.on_menu_view_results(event)

    def on_menu_view_debug(self, event):
        self.NotebookDataset.active_tab.on_menu_view_debug(event)

    def on_menu_plot_x(self, event):
        self.NotebookDataset.active_tab.on_menu_plot_x(event)


    def on_tab_changed(self, event=None):
        # This can be called programmatically, in which case event is None.
        tab = self.NotebookDataset.active_tab

        # If the voigt or giso tab is activated, we present fit-specific menus
        if tab.is_raw:
            type_ = util_menu.bar.TYPE_RAW
        elif tab.is_prep_fidsum:
            type_ = util_menu.bar.TYPE_PREP_FIDSUM
        elif tab.is_prep_timeseries:
            type_ = util_menu.bar.TYPE_PREP_FIDSUM
        elif tab.is_spectral:
            type_ = util_menu.bar.TYPE_SPECTRAL
        elif tab.is_voigt:
            type_ = util_menu.bar.TYPE_VOIGT
        elif tab.is_giso:
            type_ = util_menu.bar.TYPE_GISO
        elif tab.is_watref:
            type_ = util_menu.bar.TYPE_WATREF
        else:
            # I don't think this case ever happens.
            type_ = util_menu.bar.TYPE_RAW

        util_menu.bar.show_menus(type_)

        # Let the tab know it's active
        tab.on_activation()


#    def on_user_prior_orig(self, event):
#        dialog = dialog_user_prior.DialogUserPrior(self, self.dataset)
#        if dialog.ShowModal():
#            # Changes were made
#            self.NotebookDataset.active_tab.process_and_plot()
#        #else:
#            # User hit cancel, prior info didn't change.


    def on_user_prior(self, event):
        prior = self.dataset.user_prior
        dialog = DialogUserPrior(self, self.dataset, prior)
        if dialog.ShowModal():
            # Changes were made
            self.dataset.user_prior = dialog.prior
            self.dataset.user_prior.basis.update(self.dataset)        
            self.NotebookDataset.active_tab.process_and_plot()


    def on_user_metabolite_info(self, event):
        dialog = dialog_user_metabolite_info.DialogUserMetaboliteInfo(self, self.dataset)
        if dialog.ShowModal():
            # Changes were made
            self.NotebookDataset.active_tab.process_and_plot()
        #else:
            # User hit cancel, prior info didn't change.


    def on_preset_loaded(self):

        self.top.Freeze()
        wx.BeginBusyCursor()


        self.set_voxel_range(self.dataset.raw_dims[1])

        # reload preset values into each tab
        for item in ['raw', 'prep', 'spectral', 'fit', 'quant']:
            tab = self.get_tab(item)
            if tab:
                self.top.Freeze()
                tab.block = tab.dataset.blocks[item]
                tab.populate_controls(preset=True)
                if item == 'spectral':
                    # this refreshes the ECC or WaterFilter panels as needed
                    self.top.Layout()
                    tab.PanelSpectral.Layout()
                self.top.Thaw()
            else:
                if item == 'fit':
                    # no 'fit' tab exists, make one if needed
                    if not self.dataset.blocks['fit'].is_identity:
                        # The fitting tab follows spectral, find index
                        tab   = self._tabs["spectral"]
                        index = self.NotebookDataset.get_tab_index(tab) + 1

                        # Create the new tab and insert it
                        if isinstance(self.dataset.blocks['fit'],block_fit_voigt.BlockFitVoigt):
                            tab = tab_voigt.TabVoigt(self, self.top, self.dataset.blocks['fit'])
                        elif isinstance(self.dataset.blocks['fit'],block_fit_giso.BlockFitGiso):
                            tab = tab_giso.TabGiso(self, self.top, self.dataset.blocks['fit'])
                        
                        self._tabs["fit"] = tab
                        self.NotebookDataset.InsertPage(index, tab, 'Fitting', True)
                if item == 'quant':
                    # no 'quant' tab exists, make one if needed
                    if not self.dataset.blocks['quant'].is_identity:
                        # The fitting tab follows fitting, find index
                        tab   = self._tabs["fit"]
                        index = self.NotebookDataset.get_tab_index(tab) + 1

                        # Create the new tab and insert it
                        tab = tab_watref.TabWatref(self, self.top, self.dataset.blocks['quant'])
                        self._tabs["quant"] = tab
                        self.NotebookDataset.InsertPage(index, tab, 'Quant', True)

        self.top.notebook_datasets.Layout()
        wx.EndBusyCursor()
        self.top.Thaw()

#         for tab in self._tabs.itervalues():
#             if tab is not None:
#                 tab.process()

        # call these individually so we can turn keywords on/off
        
        for item in ['raw', 'prep', 'spectral', 'fit', 'quant']:
            tab = self.get_tab(item)
            if tab is not None:
                if item == 'raw':
                    tab.process()
                elif item == 'prep':
                    tab.process()
                elif item == 'spectral':
                    tab.process()
                    tab.svd_checklist_update()
                elif item == 'fit':
                    tab.process(entry='voxel_change')
                elif item == 'quant':
                    tab.process()
                    


        self.NotebookDataset.active_tab.plot()


    def batch_process_all(self, statusbar=None):
        """
        We are trying to do all this 'behind the scenes' in that the user does
        not see Tabs/Widgets changing as we go through the voxels. So we will
        not call the Tab.on_voxel_change() for any of the tabs.
        
        For now, we are ASSUMING that this is only called from the fit tab.
        
        """
        self.dataset.batch_process_all(statusbar=statusbar)


    def set_voxel_range(self, xmax):
        
        xval = self.voxel[0]
        
        if xval > xmax:
            self.SpinX.SetValue(xmax)
            self.on_voxel_change()
            
        self.SpinX.SetRange(1,xmax)
            
        
         
               
            

    #=======================================================
    #
    #           Widget Event Handlers
    #
    #=======================================================

    def on_voxel_change(self, event=None):
        # This event can be called programmatically, in which case event is None.
        raw = self.dataset.blocks["raw"]
        self.TextData.ChangeValue(raw.get_data_source(self.voxel))

        # Let everyone know that the current voxel changed.
        for tab in self.NotebookDataset.tabs:
            tab.on_voxel_change(self.voxel)

            if tab.is_spectral and not self.turn_off_update:
                tab.process()

        # not all pages need process() updated
        active_tab = self.NotebookDataset.active_tab

        if active_tab is not None:
    
            if active_tab.is_spectral:
                active_tab.plot()
                active_tab.plot_svd()
            elif active_tab.is_voigt:
                active_tab.process_and_plot(entry='plot_refresh')
                active_tab.view.canvas.draw()
            elif active_tab.is_giso:
                active_tab.process_and_plot(entry='voxel_change')
                active_tab.view.canvas.draw()
            elif active_tab.is_watref:
                active_tab.process_and_plot()
                #active_tab.view.canvas.draw()


    def on_scale(self, event):
        """
        This is a FloatMultiplier widget. The spin button arrow events are
        caught separate from the text field events. The spin buttons multiply
        up or divide down by the amount in the self.multiplier attribute.

        Bottom line ... by the time this event is called, the value in the
        text field is changed to some new amount.  So, now we just set an
        absolute scale value within the plot_panel_spectrum View object.
        """
        tab  = self.NotebookDataset.active_tab
        if tab.view:
            view = tab.view
            if tab.is_spectral and tab.svd_tab_active:
                view = tab.view_svd

            scale = self.FloatScale.GetValue()
            view.set_vertical_scale(1.0, abs_scale=scale)
            self.FloatScale.SetValue(view.vertical_scale)



    #=======================================================
    #
    #           Public methods
    #
    #=======================================================

    def get_tab(self, slot_name):
        """Given a slot name ('raw', 'prep', 'spectral', or 'fit'), returns
        the tab in that slot (which may be None).
        """
        return self._tabs[slot_name]


    def get_frequency_shift(self, voxel):
        '''
        Phase0, phase 1 and frequency shift are all parameters that affect the
        data in the spectral tab, however, they can also be changed in the
        fitting (all three) and svd (phase0/1 only) tabs using either widgets
        or mouse canvas events. In the end, these GUI interactions are all
        changing the same variables located in the block_spectral object.

        Because these can be changed by "between tabs" actions, I've located
        these methods at the dataset_tab level so that one tab does not ever
        talk directly to another tab, but just to a parent (or grandparent).

        '''
        return self._tabs["spectral"].block.get_frequency_shift(voxel)


    def set_frequency_shift(self, delta, voxel, auto_calc=False, entry='all'):
        '''
        Phase0, phase 1 and frequency shift are all parameters that affect the
        data in the spectral tab, however, they can also be changed in the
        fitting (all three) and svd (phase0/1 only) tabs using either widgets
        or mouse canvas events. In the end, these GUI interactions are all
        changing the same variables located in the block_spectral object.

        Because these can be changed by "between tabs" actions, I've located
        these methods at the notebook_dataset level so that one datasest does
        not ever talk directly to another tab, but just to a parent (or
        grandparent).

        '''
        spectral_tab = self._tabs["spectral"]
        fitting_tab  = self._tabs["fit"]

        b0shift = spectral_tab.block.get_frequency_shift(voxel)
        b0shift = b0shift + delta
        spectral_tab.block.set_frequency_shift(b0shift,voxel)
        spectral_tab.FloatFrequency.SetValue(b0shift)
        spectral_tab.plot_results = spectral_tab.block.chain.run([voxel]) #, entry=entry)
        spectral_tab.set_check_boxes()

        if fitting_tab:
            fitting_tab.FloatB0ShiftValue.SetValue(b0shift)
            if not auto_calc:
                key = constants.FitInitialB0ShiftMethod.MANUAL
                val = constants.FitInitialB0ShiftMethod.choices[key]
                fitting_tab.ComboInitialB0ShiftMethod.SetStringSelection( val )
                fitting_tab.block.set.initial_b0_shift_method = key
            fitting_tab.plot_results = fitting_tab.block.chain.run([voxel]) #, entry=entry)

        if spectral_tab:
            #spectral_tab.process_and_plot()    # was getting second 'process' call here, trying just a plot ... bjs
            spectral_tab.plot()
        if fitting_tab:
            fitting_tab.plot()


    def set_phase_0(self, delta, voxel, auto_calc=False):
        '''
        This method only updates block values and widget settings, not view
        display. That is done in the set_xxx_x_view() method.
        '''
        spectral_tab = self._tabs["spectral"]
        fitting_tab  = self._tabs["fit"]

        phase_0 = spectral_tab.block.get_phase_0(voxel)
        phase_0 = (phase_0  + delta) % 360
        spectral_tab.block.set_phase_0(phase_0,voxel)
        spectral_tab.FloatPhase0.SetValue(phase_0)

        if fitting_tab:
            fitting_tab.FloatInitialPhase0Value.SetValue(phase_0)
            if not auto_calc:
                key = constants.FitInitialPhaseMethod.MANUAL
                val = constants.FitInitialPhaseMethod.choices[key]
                fitting_tab.ComboInitialPhaseMethod.SetStringSelection( val )
                fitting_tab.block.set.initial_phase_method = key


    def set_phase_0_view(self, voxel):

        tab = self.NotebookDataset.active_tab
        if tab.is_spectral:
            # update phase in PlotA
            phase0   = tab.block.get_phase_0(voxel)
            if tab.svd_tab_active:
                tab.view_svd.set_phase_0(phase0, absolute=True)
            else:
                tab.view.set_phase_0(phase0, index=[0], absolute=True, no_draw=True)
                # update phase in PlotB
                if self.indexAB[1]:
                    tab_datasetb = self._outer_notebook.get_tab_by_label(self.indexAB[1])
                    tabb    = tab_datasetb.get_tab("spectral")
                    phase0b = tabb.block.get_phase_0(voxel)
                    tab.view.set_phase_0(phase0b, index=[1], absolute=True, no_draw=True)
                # update PlotC and Draw it
                tab.set_plot_c()
                tab.view.canvas.draw()

        elif tab.is_voigt:
            spectral_tab = self._tabs["spectral"]
            phase0 = spectral_tab.block.get_phase_0(voxel)
            tab.view.set_phase_0(phase0, absolute=True)

        elif tab.is_giso:
            spectral_tab = self._tabs["spectral"]
            phase0 = spectral_tab.block.get_phase_0(voxel)
            tab.view.set_phase_0(phase0, absolute=True)



    def set_phase_1(self, delta, voxel, auto_calc=False):
        '''
        Phase0, phase 1 and frequency shift are all parameters that affect the
        data in the spectral tab, however, they can also be changed in the
        fitting (all three) and svd (phase0/1 only) tabs using either widgets
        or mouse canvas events. In the end, these GUI interactions are all
        changing the same variables located in the block_spectral object.

        Because these can be changed by "between tabs" actions, I've located
        these methods at the dataset_tab level so that one tab does not ever
        talk directly to another tab, but just to a parent (or grandparent).

        I am not using pubsub here, because the Spectral tab needs to be
        updated before the Voigt tab, and blocks before tabs, and tabs before
        plots.  So, it's a function call for us to do it all, no asynchronous
        communications here, sigh.

        '''
        spectral_tab = self._tabs["spectral"]
        fitting_tab  = self._tabs["fit"]

        # check if phase1 is locked at zero
        if not spectral_tab.block.phase_1_lock_at_zero:
            phase_1 = spectral_tab.block.get_phase_1(voxel)
            phase_1 = phase_1  + delta
        else:
            phase_1 = 0.0

        spectral_tab.block.set_phase_1(phase_1,voxel)
        spectral_tab.FloatPhase1.SetValue(phase_1)
        if fitting_tab:
            fitting_tab.FloatInitialPhase1Value.SetValue(phase_1)
            if not auto_calc:
                key = constants.FitInitialPhaseMethod.MANUAL
                val = constants.FitInitialPhaseMethod.choices[key]
                fitting_tab.ComboInitialPhaseMethod.SetStringSelection( val )
                fitting_tab.block.set.initial_phase_method = key


    def set_phase_1_view(self, voxel):

        tab = self.NotebookDataset.active_tab
        if tab.is_spectral:
            # update phase in PlotA
            phase1   = tab.block.get_phase_1(voxel)
            if tab.svd_tab_active:
                tab.view_svd.set_phase_1(phase1, absolute=True)
            else:
                tab.view.set_phase_1(phase1, index=[0], absolute=True, no_draw=True)
                # update phase in PlotB
                if self.indexAB[1]:
                    tab_datasetb = self._outer_notebook.get_tab_by_label(self.indexAB[1])
                    tabb    = tab_datasetb.get_tab("spectral")
                    phase1b = tabb.block.get_phase_1(voxel)
                    tab.view.set_phase_1(phase1b, index=[1], absolute=True, no_draw=True)
                # update PlotC and Draw it
                tab.set_plot_c()
                tab.view.canvas.draw()

        elif tab.is_voigt:
            phase1 = self._tabs["spectral"].block.get_phase_1(voxel)
            tab.view.set_phase_1(phase1, absolute=True)

        elif tab.is_giso:
            phase1 = self._tabs["spectral"].block.get_phase_1(voxel)
            tab.view.set_phase_1(phase1, absolute=True)

    #=======================================================
    #
    #           Internal helpers
    #
    #=======================================================





