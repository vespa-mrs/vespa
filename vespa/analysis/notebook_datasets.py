# Python modules
import os.path

# 3rd party modules
import wx
import wx.html
#import wx.aui as aui
import wx.lib.agw.aui as aui        # NB. wx.aui version throws odd wxWidgets exception on Close/Exit  ?? Not anymore in wxPython 4.0.6 ??
from pubsub import pub as pubsub

# Our modules
import vespa.analysis.util_menu as util_menu
import vespa.analysis.tab_dataset as tab_dataset
import vespa.analysis.block_prep_fidsum as block_prep_fidsum
import vespa.analysis.block_raw_probep as block_raw_probep
import vespa.analysis.block_prep_wbnaa as block_prep_wbnaa
import vespa.analysis.block_raw_cmrr_slaser as block_raw_cmrr_slaser
import vespa.analysis.block_raw_edit_fidsum as block_raw_edit_fidsum
import vespa.analysis.block_prep_timeseries as block_prep_timeseries

import vespa.common.util.misc as util_misc
import vespa.common.wx_gravy.util as wx_util
import vespa.common.wx_gravy.notebooks as vespa_notebooks


class NotebookDatasets(vespa_notebooks.VespaAuiNotebook):
    # I need the path to the welcome page image which is in vespa/common.
    _path = util_misc.get_vespa_install_directory()
    _path = os.path.join(_path, "common", "resources", "analysis_welcome.png")

    WELCOME_TAB_TEXT = """
    <html><body>
    <h1>Welcome to Vespa - Analysis</h1>
    <img src="%s" alt="Time-Freq Plots" />
    <p><b>Currently there are no datasets loaded.</b></p>
    <p>You can use the File menu to load data files.</p>
    </body></html>
    """ % _path
    # I tidy up my namespace by deleting this temporary variable.
    del _path


    def __init__(self, top, parent):

        vespa_notebooks.VespaAuiNotebook.__init__(self, parent)

        self.top    = top
        self.parent = parent
        self.count  = 0

        self.suspend_pubsub = False

        self.associated_close = []

        self.show_welcome_tab()

        self.Bind(aui.EVT_AUINOTEBOOK_PAGE_CLOSE,    self.on_tab_close)
        self.Bind(aui.EVT_AUINOTEBOOK_PAGE_CLOSED,   self.on_tab_closed)

        # Apparently this was a bug waiting to happen because I started to get
        # a lot of 'can not access bitmap' in the plot_panels in Spectral Tab.
        # Never fully worked out why, but we have a tab_changed event at the
        # TabDataset level, so we likely did not need this one anyway.
        #
        # self.Bind(aui.EVT_AUINOTEBOOK_PAGE_CHANGED,  self.on_tab_changed)


    #=======================================================
    #
    #           Global and Menu Event Handlers
    #
    #=======================================================

    def on_menu_view_option(self, event):
        if self.active_tab:
            self.active_tab.on_menu_view_option(event)

    def on_menu_view_output(self, event):
        if self.active_tab:
            self.active_tab.on_menu_view_output(event)

    def on_menu_view_results(self, event):
        if self.active_tab:
            self.active_tab.on_menu_view_results(event)

    def on_menu_view_debug(self, event):
        if self.active_tab:
            self.active_tab.on_menu_view_debug(event)

    def on_menu_plot_x(self, event):
        if self.active_tab:
            self.active_tab.on_menu_plot_x(event)

    def on_add_voigt_tab(self, event):
        if self.active_tab:
            self.active_tab.on_add_voigt_tab(event)

    def on_add_giso_tab(self, event):
        if self.active_tab:
            self.active_tab.on_add_giso_tab(event)

    def on_add_watref_tab(self, event):
        if self.active_tab:
            self.active_tab.on_add_watref_tab(event)

    def on_tab_changed(self, event):
        wx.GetApp().vespa.update_title()

        if self.active_tab:
            self.active_tab.on_tab_changed()


    def on_tab_close(self, event):
        """
        This is a two step event. Here we give the user a chance to cancel
        the Close action. If user selects to continue, then the on_tab_closed()
        event will also fire.  event.GetPageCound()

        """
        assoc_classes = (block_raw_edit_fidsum.BlockRawEditFidsum,
                         block_raw_cmrr_slaser.BlockRawCmrrSlaser,
                         block_raw_probep.BlockRawProbep)
        
        
        raw = self.active_tab.dataset.blocks['raw']
        if isinstance(raw, assoc_classes):
            if not self.associated_close:
                if not self.top.close_all:
                    # This is not a menu close_all event ... give them a chance to back out of it
                    msg = "Are you sure you want to close this dataset and all associated datasets?"
                    if wx.MessageBox(msg, "Close Dataset", wx.YES_NO, self) != wx.YES:
                        event.Veto()
                    else:
                        # check if there are any other tabs to close
                        tabs = self.global_poll_associated_tabs(self.active_tab.dataset)
                        self.associated_close = []
                        for tab in tabs:
                            self.associated_close = [tab,]
                            self.activate_tab(tab=tab)
                            self.close_active_dataset()
                            self.associated_close = []
                else:
                    # This is a menu close_all event ...
                    # check if there are any other tabs to close
                    tabs = self.global_poll_associated_tabs(self.active_tab.dataset)
                    self.associated_close = []
                    for tab in tabs:
                        self.associated_close = [tab,]
                        self.activate_tab(tab=tab)
                        self.close_active_dataset()
                        self.associated_close = []
                        
#        elif isinstance(raw, block_raw_cmrr_slaser.BlockRawCmrrSlaser):
#            # "CMRR sLaser" blocks have multiple datasets that are all opened
#            #  from the same twix file. We want to close them all at the same time.
#            if not self.associated_close:
#                if not self.top.close_all:
#                    # This is not a menu close_all event ...
#                    # give them a chance to back out of it
#                    msg = "Are you sure you want to close this dataset and all associated datasets?"
#                    if wx.MessageBox(msg, "Close Dataset", wx.YES_NO, self) != wx.YES:
#                        event.Veto()
#                    else:
#                        # check if there are any other tabs to close
#                        tabs = self.global_poll_associated_tabs(self.active_tab.dataset)
#                        self.associated_close = []
#                        for tab in tabs:
#                            self.associated_close = [tab,]
#                            self.activate_tab(tab=tab)
#                            self.close_active_dataset()
#                            self.associated_close = []
#                else:
#                    # This is a menu close_all event ...
#                    # check if there are any other tabs to close
#                    tabs = self.global_poll_associated_tabs(self.active_tab.dataset)
#                    self.associated_close = []
#                    for tab in tabs:
#                        self.associated_close = [tab,]
#                        self.activate_tab(tab=tab)
#                        self.close_active_dataset()
#                        self.associated_close = []                    
        else:
            if not self.top.close_all:
                # Not a close_all menu event ...
                # This is a single file, no need to populate associated_close list
                msg = "Are you sure you want to close this dataset?"
                if wx.MessageBox(msg, "Close Dataset", wx.YES_NO, self) != wx.YES:
                    event.Veto()


    def on_tab_closed(self, event):
        """
        At this point the page is already closed and the dataset removed from
        the datasets list. We need to pubsub that the keys have changed.

        """
        if self.associated_close == []:

            if self.GetPageCount() == 0:
                util_menu.bar.show_menus(util_menu.AnalysisMenuBar.TYPE_START)
                self.show_welcome_tab()
            else:
                # let the world know the keys have changed
                keys = sorted(self.top.datasets.keys())
                pubsub.sendMessage("dataset_keys_change", keys=keys)
        # else:
        #     there are associated tabs that need to be closed

    def on_user_prior(self, event):
        if self.active_tab:
            self.active_tab.on_user_prior(event)


    def on_user_metabolite_info(self, event):
        if self.active_tab:
            self.active_tab.on_user_metabolite_info(event)


    def on_preset_loaded(self):
        self.active_tab.on_preset_loaded()


    #=======================================================
    #
    #           Public methods shown below
    #             in alphabetical order
    #
    #=======================================================

    def add_dataset_tab(self, datasets=None):
        # If the welcome tab is open, close it.
        if self.is_welcome_tab_open:
            self.DeletePage(0)

        names, count = self.custom_tab_names(datasets, self.count)
        self.count = count

        for i,dataset in enumerate(datasets):

            name = names[i]

            # register Dataset object at top level
            self.top.datasets[name] = dataset

            # create new notebook tab with process controls
            dtab = tab_dataset.TabDataset(self, self.top, name)
            self.AddPage(dtab, name, True)

        # let the world know the keys have changed
        keys = sorted(self.top.datasets.keys())
        pubsub.sendMessage("dataset_keys_change", keys=keys)


    def close_active_dataset(self):
        if self.active_tab:
            wx_util.send_close_to_active_tab(self)


    def get_tab_by_label(self, label):
        for i in range(self.GetPageCount()):
            if label == self.GetPageText(i):
                return self.GetPage(i)

        return None


    def global_poll_associated_tabs(self, dataset):
        """
        Here is the list of object classes that we polling for associated tabs:
        
          block_raw_edit_fidsum.BlockRawEditFidsum
          block_raw_cmrr_slaser.BlockRawCmrrSlaser  
        
        
        """
        assoc_classes = (block_raw_edit_fidsum.BlockRawEditFidsum,
                         block_raw_cmrr_slaser.BlockRawCmrrSlaser,
                         block_raw_probep.BlockRawProbep)
        tabs = []
        raw = dataset.blocks['raw']

        # if isinstance(raw, block_raw_edit_fidsum.BlockRawEditFidsum):
        #     ids = [raw.data_on.id, raw.data_off.id, raw.data_sum.id, raw.data_dif.id]
        #     for label in list(self.top.datasets.keys()):
        #         tab = self.get_tab_by_label(label)
        #         if tab.dataset.id in ids and tab.dataset.id != dataset.id:
        #             tabs.append(tab)
        #
        # elif isinstance(raw, block_raw_probep.BlockRawProbep):
        #     ids = [raw.data_on.id, raw.data_off.id, raw.data_sum.id, raw.data_dif.id]
        #     for label in list(self.top.datasets.keys()):
        #         tab = self.get_tab_by_label(label)
        #         if tab.dataset.id in ids and tab.dataset.id != dataset.id:
        #             tabs.append(tab)
        #
        # elif isinstance(raw, block_raw_cmrr_slaser.BlockRawCmrrSlaser):
        #     ids = [raw.data_coil_combine, raw.data_ecc1, raw.data_water1, raw.data_metab64, raw.data_ecc2, raw.data_water2]
        #     ids = [item.id for item in ids if item is not None]
        #     for label in list(self.top.datasets.keys()):
        #         tab = self.get_tab_by_label(label)
        #         if tab.dataset.id in ids and tab.dataset.id != dataset.id:
        #             tabs.append(tab)

        if isinstance(raw, assoc_classes):
            ids = raw.get_associated_datasets()
            ids = [item.id for item in ids if item is not None]
            for label in list(self.top.datasets.keys()):
                tab = self.get_tab_by_label(label)
                if tab.dataset.id in ids and tab.dataset.id != dataset.id:
                    tabs.append(tab)
        else:
            return tabs

        
        return tabs


    def global_block_zerofill_update(self, zf_mult):
        """
        We have a rule that only datasets with the same dimensionality are
        allowed to be loaded into the notebook at the same time. This is
        checked when the raw data is  loaded, but also has to be checked
        when the zerofill widget is changed on any tab. This ensures that
        multiple sets of data can be rationally compared in plots.

        """
        for tab in self.tabs:
            # Reset dataset results arrays in blocks and chains.
            tab.dataset.update_for_zerofill_change(zf_mult)


    def global_tab_zerofill_update(self, zf_mult):
        """
        As opposed to the global_block_zerofill_update() method, here
        we need to do all data processing (with no plotting) FIRST just
        in case there are tabs where sync is on, all datasets have to be
        the same size to plot correctly.

        Thus, this code is implemented at the Notebook level rather than
        at each TabDataset in the self.tabs list.

        """
        self.suspend_pubsub = True
        for tab in self.tabs:
            # update GUI stuff

            # tab_preprocess = tab.get_tab("prep")      # bjs - not sure if needed for ZF change
            # if tab_preprocess:
            #     tab_preprocess.process()

            # do processing in spectral and fitting tabs as needed
            tab_spectral = tab.get_tab("spectral")
            tab_spectral.ComboZeroFill.SetStringSelection(str(zf_mult))
            tab_spectral.process()

            tab_fit = tab.get_tab("fit")
            if tab_fit:
                tab_fit.process()

            tab_quant = tab.get_tab("quant")
            if tab_quant:
                tab_quant.process()

        self.suspend_pubsub = False

        # now do plotting
        for tab in self.tabs:
            block_tab = tab.NotebookDataset.active_tab
            if block_tab.view:
                block_tab.plot()


    def global_poll_phase(self, poll_labels, delta, voxel, do_zero=True):
        """
        Phase 0 and phase 1 are parameters that can affect one or more
        view panels in tabs within a dataset tab. They can be changed using
        either widgets or mouse canvas events. For a given dataset, these
        events only change one phase 0 or phase 1 variable located in the
        block_spectral object.

        These values can also be changed "between dataset tabs" due to the
        PlotB, option on the Spectral tab when combined with the Sync check
        box. So, I've located these methods at the notebook level so that one
        tab does not ever talk directly to another tab, but just to a parent
        (or grandparent).

        The "poll_labels" list whether one or more datasets are involved
        from the calling event. This is typically only more than one label
        when the event is from the Spectral tab.

        """
        # determine which dataset tabs are affected by this event
        poll_labels, view_labels = self._get_poll_and_view_labels(poll_labels)
        # update values in dataset tabs
        for label in poll_labels:
            tab = self.get_tab_by_label(label)
            if do_zero:
                tab.set_phase_0(delta, voxel)
            else:
                tab.set_phase_1(delta, voxel)

        # refresh views in tabs where values have changed
        for label in view_labels:
            tab = self.get_tab_by_label(label)
            if do_zero:
                tab.set_phase_0_view(voxel)
            else:
                tab.set_phase_1_view(voxel)


    def global_poll_frequency_shift(self, poll_labels, delta, voxel):
        """
        Frequency shift changes can affect one or more view panels in tabs
        within a dataset tab. For a given dataset, this event really only
        changes one variable located in the block_spectral object.

        This value can also be changed "between dataset tabs" due to the
        PlotB, option on the Spectral tab when combined with the Sync check
        box. So, I've located this methods at the notebook level so that one
        tab does not ever talk directly to another tab, but just to a parent
        (or grandparent).

        The "poll_labels" list whether one or more datasets are involved
        from the calling event. This is typically only more than one label
        when the event is from the Spectral tab.

        """
        poll_labels, view_labels = self._get_poll_and_view_labels(poll_labels)

        # update phase values in dataset tabs
        for label in poll_labels:
            tab = self.get_tab_by_label(label)
            tab.set_frequency_shift(delta, voxel)


    def global_poll_sync_event(self, poll_labels, value, voxel=None, event=''):
        """
        Various spectral parameters can affect one or more view panels in tabs
        within a dataset tab. For a given dataset, this event really only
        changes one variable located in the block_spectral object.

        These value can also be changed "between dataset tabs" due to the
        PlotB, option on the Spectral tab when combined with the Sync check
        box. So, I've located this methods at the notebook level so that one
        tab does not ever talk directly to another tab, but just to a parent
        (or grandparent).

        The "poll_labels" list whether one or more datasets are involved
        from the calling event. This is typically only more than one label
        when the event is from the Spectral tab.

        """
        if event:

            poll_labels, _ = self._get_poll_and_view_labels(poll_labels)

            # update phase values in dataset tabs
            for label in poll_labels:
                tab_dataset  = self.get_tab_by_label(label)
                tab_spectral = tab_dataset.get_tab("spectral")
                if event == 'flip':
                    tab_spectral.set_flip(value)
                elif event == 'chop':
                    tab_spectral.set_chop(value)
                elif event == 'frequency_shift_lock':
                    tab_spectral.set_frequency_shift_lock(value)
                elif event == 'phase_lock':
                    tab_spectral.set_phase_lock(value)
                elif event == 'phase1_zero':
                    tab_spectral.set_phase1_zero(value)
                elif event == 'left_shift_correction':
                    tab_spectral.set_left_shift_correction(value)
                elif event == 'phase1_pivot':
                    tab_spectral.set_phase1_pivot(value)
                elif event == 'dc_offset':
                    tab_spectral.set_dc_offset(value)
                elif event == 'left_shift':
                    tab_spectral.set_left_shift(value)



    def _get_poll_and_view_labels(self, poll_labels):
        """
        This is a helper function used by the global_poll_xxx() methods.

        We poll all dataset tabs to determine which are affected by the event
        triggered by the calling dataset(s) listed in the poll_labels list as
        it is sent into the method.

        If a dataset is in the list, and it has a PlotB selected, the PlotB
        dataset is added to the list.

        If a dataset is synched to a PlotB that is in the list, then that
        dataset is added to the list.

        We also poll tab labels for ANY dataset that is involved in this event
        whether in PlotA or PlotB, because we need to refresh those views.

        Unique values for both the poll_list and view_list are returned.

        """
        view_labels = list(poll_labels)
        for label in list(self.top.datasets.keys()):
            tab_dataset  = self.get_tab_by_label(label)
            tab_spectral = tab_dataset.get_tab("spectral")
            if label in poll_labels:
                poll_labels.append(tab_dataset.indexAB[0])
                if tab_spectral.do_sync:
                    poll_labels.append(tab_dataset.indexAB[1])
                    view_labels.append(tab_dataset.indexAB[1])
                view_labels.append(tab_dataset.indexAB[0])
            if tab_dataset.indexAB[1] in poll_labels:
                poll_labels.append(tab_dataset.indexAB[1])
                if tab_spectral.do_sync:
                    poll_labels.append(tab_dataset.indexAB[0])
                view_labels.append(tab_dataset.indexAB[0])
                view_labels.append(tab_dataset.indexAB[1])
        poll_labels = list(set(poll_labels))
        view_labels = list(set(view_labels))
        return poll_labels, view_labels


    #=======================================================
    #
    #           Internal methods start here
    #
    #=======================================================



    ###############    Overloaded event handlers re. wx.lib.aui.auibook.py 
    ###############    
    
    
    def OnTabButton(self, event):
        """
        Handles the ``EVT_AUINOTEBOOK_BUTTON`` event for :class:`AuiNotebook`.

        :param `event`: a :class:`AuiNotebookEvent` event to be processed.
        
        This is overloaded here because the original calls Destroy() through
        a CallAfter() method. This messes up a number of code steps downstream
        that need to know how many tabs remain in the notebook (e.g. putting
        the StartTab back in place if the last Experiment tab is closed).
        
        """

        tabs = event.GetEventObject()
        button_id = event.GetInt()

        if button_id == aui.AUI_BUTTON_CLOSE:

            selection = event.GetSelection()

            if selection == -1:

                # if the close button is to the right, use the active
                # page selection to determine which page to close
                selection = tabs.GetActivePage()

            if selection == -1 or not tabs.GetEnabled(selection):
                return

            if selection != -1:

                close_wnd = tabs.GetWindowFromIdx(selection)

                if close_wnd.GetName() == "__fake__page__":
                    # This is a notebook preview
                    previous_active, page_status = close_wnd.__previousStatus
                    for page, status in zip(tabs.GetPages(), page_status):
                        page.enabled = status

                    main_idx = self._tabs.GetIdxFromWindow(close_wnd)
                    self.DeletePage(main_idx)

                    if previous_active >= 0:
                        tabs.SetActivePage(previous_active)
                        page_count = tabs.GetPageCount()
                        selection = -1

                        for page in range(page_count):
                            # remove the page from the source tabs
                            page_info = tabs.GetPage(page)
                            if page_info.active:
                                selection = page
                                break

                        tabs.DoShowHide()
                        self.DoSizing()
                        tabs.Refresh()

                        if selection >= 0:
                            wx.CallAfter(tabs.MakeTabVisible, selection, self)

                    # Don't fire the event
                    return

                # ask owner if it's ok to close the tab
                e = aui.AuiNotebookEvent(aui.wxEVT_COMMAND_AUINOTEBOOK_PAGE_CLOSE, self.GetId())
                idx = self._tabs.GetIdxFromWindow(close_wnd)
                e.SetSelection(idx)
                e.SetOldSelection(event.GetSelection())
                e.SetEventObject(self)
                self.GetEventHandler().ProcessEvent(e)
                if not e.IsAllowed():
                    return

                if repr(close_wnd.__class__).find("AuiMDIChildFrame") >= 0:
                    close_wnd.Close()

                else:
                    main_idx = self._tabs.GetIdxFromWindow(close_wnd)
                    self.DeletePage(main_idx)
                    #wx.CallAfter(self.DeletePage, main_idx)

                # notify owner that the tab has been closed
                e2 = aui.AuiNotebookEvent(aui.wxEVT_COMMAND_AUINOTEBOOK_PAGE_CLOSED, self.GetId())
                e2.SetSelection(idx)
                e2.SetEventObject(self)
                self.GetEventHandler().ProcessEvent(e2)

                if self.GetPageCount() == 0:
                    mgr = self.GetAuiManager()
                    win = mgr.GetManagedWindow()
                    win.SendSizeEvent()


    def custom_tab_names(self, datasets, count):
        """ Based on 'raw' Block class, create more informative Tab label"""

        dataset = datasets[0] if util_misc.is_iterable(datasets) else datasets

        if isinstance(dataset.blocks['raw'], block_raw_probep.BlockRawProbep):
            count += 1
            names = []
            for dataset in datasets:
                fname = dataset.blocks['raw'].data_sources[0]
                if ".water" in fname:
                    names.append("Probe%d.Water" % count)
                elif ".metab" in fname:
                    names.append("Probe%d.Metab" % count)
                else:
                    names.append("Probe%d.Unknown" % count)

        elif isinstance(dataset.blocks['raw'], block_raw_edit_fidsum.BlockRawEditFidsum):
            base = "Edit%d." % count
            names = [base+"On", base+"Off", base+"Sum", base+"Dif"]
            count += 1

        elif isinstance(dataset.blocks['raw'], block_raw_cmrr_slaser.BlockRawCmrrSlaser):
            count += 1
            iother = 1
            names = []
            for dataset in datasets:
                fname = dataset.blocks['raw'].data_sources[0]
                if ".combine" in fname:
                    names.append("CMRR%d." % count + "Coil")
                elif ".ecc1" in fname:
                    names.append("CMRR%d." % count + "Ecc1")
                elif ".water1" in fname:
                    names.append("CMRR%d." % count + "Water1")
                elif ".metab64" in fname:
                    names.append("CMRR%d." % count + "Metab64")
                elif ".ecc2" in fname:
                    names.append("CMRR%d." % count + "Ecc2")
                elif ".water2" in fname:
                    names.append("CMRR%d." % count + "Water2")
                else:
                    names.append("CMRR%d." % count + "Other" + str(iother))  # typically for MMol Basis dataset, but unknowns too.
                    iother += 1

        elif isinstance(dataset.blocks['prep'], block_prep_fidsum.BlockPrepFidsum):
            names = []
            for dataset in datasets:
                count += 1
                names.append("Fidsum%d" % count)

        elif isinstance(dataset.blocks['prep'], block_prep_timeseries.BlockPrepTimeseries):
            names = []
            for dataset in datasets:
                count += 1
                names.append("Timeseries%d" % count)

        elif isinstance(dataset.blocks['prep'], block_prep_wbnaa.BlockPrepWbnaa):
            count += 1
            names = ["Wbnaa %d" % count,]

        else:
            names = []
            for dataset in datasets:
                count += 1
                names.append("Dataset%d" % count)

        return names, count