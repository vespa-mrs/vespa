#!/usr/bin/env python

# Python modules
import os
import sys
import webbrowser
import struct

# 3rd party modules
import wx
import wx.adv as wx_adv
import wx.lib.agw.aui as aui        # NB. wx.aui version throws odd wxWidgets exception on Close/Exit  ?? Not anymore in wxPython 4.0.6 ??
import numpy as np
#import wx.aui as aui

# Our modules
import vespa.analysis.mrs_dataset as mrs_dataset
import vespa.analysis.notebook_datasets as notebook_datasets
import vespa.analysis.util_menu as util_menu
import vespa.analysis.util_analysis_config as util_analysis_config
import vespa.analysis.util_import as util_import
import vespa.analysis.util_file_import as util_file_import

import vespa.common.constants as common_constants
import vespa.common.util.init as util_init
import vespa.common.util.misc as util_misc
import vespa.common.util.export as util_export
import vespa.common.wx_gravy.common_dialogs as common_dialogs
import vespa.common.wx_gravy.util as wx_util
import vespa.common.images as images
import vespa.common.dialog_export as dialog_export


_MSG_PRESET_MISMATCH = """
One or more of the selected VIFF files is an Analysis Preset file. These can not be opened as datasets.
"""
 
_MSG_OPEN_ATTRIBUTE_MISMATCH = """
The dimensions and/or sweep width of the currently open datasets differ from those of the file(s) you asked to open.
 
You can open these files, but first you have to close all currently open datasets.
"""
 
_MSG_OPEN_ZEROFILL_MISMATCH = """
The zerofill factor of the currently open datasets differ from those of the file you asked to open.
 
You can open this file, but the zero fill factor of the open datasets needs to be changed.
"""
 
_MSG_NO_DATASETS_FOUND = """The file "%s" doesn't contain any datasets."""


#----------------------------------------------------------------------


class Main(wx.Frame):

    def __init__(self, db, position, size):

        self._left,  self._top    = position
        self._width, self._height = size

        style = wx.CAPTION | wx.CLOSE_BOX | wx.MINIMIZE_BOX |           \
                wx.MAXIMIZE_BOX | wx.SYSTEM_MENU | wx.RESIZE_BORDER |   \
                wx.CLIP_CHILDREN

        wx.Frame.__init__(self, None, wx.ID_ANY, "Analysis",
                          (self._left, self._top),
                          (self._width, self._height), style)

        self.db = db
        
        # flags for global control ----------------------------------
        
        self.close_all = False

        # GUI Creation ----------------------------------------------

        self.datasets = {}

        self._mgr = aui.AuiManager()
        self._mgr.SetManagedWindow(self)

        self.SetIcon(images.icon4_128_analysis_monogram.GetIcon())

        self.statusbar = self.CreateStatusBar(4, 0)
        self.statusbar.SetStatusText("Ready")

        # I make the status bar and the update_title method globally available
        # because multiple places in the app want to use them.
        wx.GetApp().vespa.statusbar = self.statusbar
        wx.GetApp().vespa.update_title = self.update_title

        # set up default and user import data classes
        #
        # TODO - FIXME - bjs this is a hack for 2to3 convert ... rethink in future
        if (sys.version_info > (3, 0)):
            # Python 3 code in this block
            fname = os.path.join(util_misc.get_data_dir(), "analysis_import_menu_additions_py3.ini")
        else:
            fname = os.path.join(util_misc.get_data_dir(), "analysis_import_menu_additions.ini")
        classes, full_cfg, msg = util_file_import.set_import_data_classes(filename=fname)
        if msg:
            # some class was not found ... but we continue
            common_dialogs.message(msg)
        self._import_data_classes = classes

        bar = util_menu.AnalysisMenuBar(self, full_cfg)
        self.SetMenuBar(bar)
        util_menu.bar = bar
 
        self.build_panes()
        self.bind_events()


    ##############                                    ############
    ##############     Internal helpers are below     ############
    ##############       in alphabetical order        ############
    ##############                                    ############

    def bind_events(self):
        self.Bind(wx.EVT_MENU,  self.on_menu_exit, id=wx.ID_EXIT)
        self.Bind(wx.EVT_CLOSE, self.on_self_close)
        self.Bind(wx.EVT_SIZE,  self.on_self_coordinate_change)
        self.Bind(wx.EVT_MOVE,  self.on_self_coordinate_change)


    def build_panes(self):
        # Create & show these panes

        self.notebook_datasets  = notebook_datasets.NotebookDatasets(self, self)

        # Center pane
        self._mgr.AddPane(self.notebook_datasets, aui.AuiPaneInfo().CenterPane())


        # "commit" all changes made to AuiManager
        self._mgr.Update()



    ##############                                    ############
    ##############      Event handlers are below      ############
    ##############       in alphabetical order        ############
    ##############                                    ############

    def on_menu_exit(self, event):
        self.Close(False)


    def on_self_close(self, event):

        if not self.notebook_datasets.is_welcome_tab_open:
            msg = "Are you sure you want to exit Analysis?"
            if wx.MessageBox(msg, "Quit Analysis", wx.YES_NO, self) != wx.YES:
                event.Veto()
                return

        # Save my coordinates
        config = util_analysis_config.Config()

        config.set_window_coordinates("main", self._left, self._top,
                                      self._width, self._height)

        config.set_window_maximized("main", self.IsMaximized())
        config.write()

        self._mgr.UnInit()      # needed to avoid wx._core.wxAssertionError: C++ assertion "GetEventHandler() == this" failed at ..\..\src\common\wincmn.cpp 
        self.Destroy()


    def on_self_coordinate_change(self, event):
        # This is invoked for move & size events
        if self.IsMaximized() or self.IsIconized():
            # Bah, forget about this. Recording coordinates doesn't make sense
            # when the window is maximized or minimized. This is only a
            # concern on Windows; GTK and OS X don't produce move or size
            # events when a window is minimized or maximized.
            pass
        else:
            if event.GetEventType() == wx.wxEVT_MOVE:
                self._left, self._top = self.GetPosition()
            else:
                # This is a size event
                self._width, self._height = self.GetSize()



    ##############                                    ############
    ##############       Menu handlers are below      ############
    ##############       in the order they appear     ############
    ##############             on the menu            ############
    ##############                                    ############

    ######  File menu  ######

    def on_open_viff(self, event):
        self.Freeze()
        wx.BeginBusyCursor()
        self._open_viff_dataset_file()
        wx.EndBusyCursor()
        self.Thaw()
        self.notebook_datasets.Layout()


    def on_import_mrs_data_raw(self, event):
        datasets = self._import_viff_raw_file()
        if datasets:
            # although this format only allows us to import single files, it
            # still returns the one dataset in a list to mesh with the
            # functionality in add_dataset_tab() to deal with multiple files
            self.Freeze()
            wx.BeginBusyCursor()
            self.notebook_datasets.add_dataset_tab(datasets)
            wx.EndBusyCursor()
            self.Thaw()
            self.notebook_datasets.Layout()


    def on_import_user_item(self, event):
        # This is triggered for all user-defined menu items created via
        # analysis_menu_additions.ini. They hang off of File/Import.

        # Get the id of the menu item that triggered this event
        id_ = event.GetId()

        # Find the reader class & INI file name associated with that id.
        section_name, ini_name = util_menu.bar.get_user_menu_item_info(id_)

        reader, ini_name = self._import_data_classes[section_name]

        self._import_file(reader(), ini_name)


    def on_save_viff(self, event, save_as=False):
        # This event is also called programmatically by on_save_as_viff().
        if self.notebook_datasets.active_tab:
            dataset = self.notebook_datasets.active_tab.dataset

            filename = dataset.dataset_filename
            if filename and (not save_as):
                # This dataset already has a filename which means it's already
                # associated with a VIFF file. We don't bug the user for a
                # filename, we just save it.
                pass
            else:
                # Prompt the user for the save filename & location
                if filename:
                    # Prompt using the existing VIFF filename
                    path, filename = os.path.split(filename)
                else:
                    # Construct a filename from the raw filename.
                    raw = dataset.blocks["raw"]
                    filename = raw.data_source

                    path, filename = os.path.split(filename)
                    filename = os.path.splitext(filename)[0] + ".xml"

                filename = common_dialogs.save_as("Save As XML/VIFF (Vespa Interchange Format File)",
                                                  "VIFF/XML files (*.xml)|*.xml",
                                                  path, filename)

            if filename:
                dataset.dataset_filename = filename

                wx.BeginBusyCursor()
                self._save_viff(dataset)
                wx.EndBusyCursor()


    def on_save_as_viff(self, event):
        self.on_save_viff(event, True)


    def on_close_dataset(self, event):
        if self.notebook_datasets:
            self.notebook_datasets.close_active_dataset()

    def on_close_all(self, event):
        msg = "This will close all open datasets with no opportunity to save results, continue?"
        if wx.MessageBox(msg, "Close All Datasets", wx.YES_NO, self) != wx.YES:
            event.Veto()
        else:
            self.close_all = True
            while self.datasets:
                self.notebook_datasets.close_active_dataset()
            self.close_all = False


    def on_load_preset_from_file(self, event):

        tab = self.notebook_datasets.active_tab

        if tab:

            dataset = tab.dataset
            file_type = common_constants.MrsFileTypes.VIFF
            ini_name  = "load_preset"       # "open_viff"
            default_path = util_analysis_config.get_path(ini_name)

            # Note that we only allow people to open a single VIFF file as
            # opposed to DICOM & VASF where we allow them to open multiple.
            filetype_filter="Spectra Preset (*.xml,*.xml.gz,*.viff,*.vif)|*.xml;*.xml.gz;*.viff;*.vif"
            filename = common_dialogs.pickfile(filetype_filter=filetype_filter,
                                                multiple=False,
                                                default_path=default_path)
            if filename:
                msg = ""
                try:
                    importer = util_import.DatasetImporter(filename)
                except IOError:
                    msg = """I can't read the preset file "%s".""" % filename
                except SyntaxError:
                    msg = """The preset file "%s" isn't valid Vespa Interchange File Format.""" % filename

                if msg:
                    common_dialogs.message(msg, "Analysis - Open Preset File",  common_dialogs.E_OK)
                else:
                    # Time to rock and roll!
                    wx.BeginBusyCursor()
                    presets = importer.go()
                    wx.EndBusyCursor()

                    preset = presets[0]

                    wx.BeginBusyCursor()

                    # update dataset object with preset blocks and chains
                    dataset.apply_preset(presets[0], voxel=tab.voxel)

                    # refresh chain in each block for new preset values
                    chain_outputs = {}
                    chain_outputs['raw']      = dataset.blocks['raw'].chain.run([(0,0,0)], entry='all')
                    chain_outputs['prep']     = dataset.blocks['prep'].chain.run([(0,0,0)], entry='all', freq_raw=True)
                    chain_outputs['spectral'] = dataset.blocks['spectral'].chain.run([(0, 0, 0)], entry='all')
                    chain_outputs['fit']      = dataset.blocks['fit'].chain.run([(0, 0, 0)], entry='initial_only')
                    chain_outputs['spectral'] = dataset.blocks['spectral'].chain.run([(0, 0, 0)], entry='all')
                    chain_outputs['fit']      = dataset.blocks['fit'].chain.run([(0, 0, 0)], entry='initial_only')
                    chain_outputs['quant']    = dataset.blocks['quant'].chain.run([(0, 0, 0)], entry='initial_only')

                    # Get the active dataset tab to update itself
                    self.notebook_datasets.on_preset_loaded()

                    # Check dimensionality on ALL other loaded datasets
                    preset_zf = dataset.zero_fill_multiplier
                    self.notebook_datasets.global_block_zerofill_update(preset_zf)

                    # Check gui on ALL tabs for zero fill consistency
                    self.notebook_datasets.global_tab_zerofill_update(preset_zf)

                    path, _ = os.path.split(filename)
                    util_analysis_config.set_path(ini_name, path)

                    wx.EndBusyCursor()


    def on_load_preset_from_tab(self, event):
        print("Not yet Implemented - on_load_preset_from_tab")


    def on_save_preset_to_file(self, event):

        tab = self.notebook_datasets.active_tab
        if tab:
            # get the dataset object if it exists, save settings to preset file
            dataset  = tab.dataset
            filename = 'vespa_analysis_preset.xml'
            dialog = dialog_export.DialogExport(self, False, filename=filename)
            dialog.ShowModal()

            if dialog.export:
                filename = dialog.filename
                comment  = dialog.comment
                compress = dialog.compress

                wx.BeginBusyCursor()
                # Many of these objects can deflate themselves into a
                # more compact form when being deflated as presets, so
                # we temporarily set a flag on them to let them know
                # how to behave during deflate(). This hack turns out
                # to be easier than passing a flag to the deflate()
                # calls.
                dataset.set_behave_as_preset(True)

                try:
                    util_export.export(filename, [dataset], None, comment, compress)
                except IOError as xxx_todo_changeme:
                    (error_number, error_string) = xxx_todo_changeme.args
                    msg = """Exporting to "%s" failed. The operating system message is below --\n\n""" % filename
                    msg += error_string
                    common_dialogs.message(msg, "Vespa Export", common_dialogs.E_OK)

                wx.EndBusyCursor()

                # Here we undo the behave_as_preset hack from above.
                dataset.set_behave_as_preset(False)


    def on_close_window(self, event):
        self.Destroy()


    ######  Processing menu  ######

    def on_add_voigt_tab(self, event):
        self.notebook_datasets.on_add_voigt_tab(event)

    def on_add_giso_tab(self, event):
        self.notebook_datasets.on_add_giso_tab(event)

    def on_add_watref_tab(self, event):
        self.notebook_datasets.on_add_watref_tab(event)

    def on_user_prior(self, event):
        self.notebook_datasets.on_user_prior(event)

    def on_user_metabolite_info(self, event):
        self.notebook_datasets.on_user_metabolite_info(event)

    def on_show_inspection_tool(self, event):
        wx_util.show_wx_inspector(self)


    ######  View  menu  ######

    # View options affect plots in a dataset - we pass these events on to the
    # specific dataset notebook tab that needs to react to them.

    def on_menu_view_option(self, event):
        self.notebook_datasets.on_menu_view_option(event)

    def on_menu_view_output(self, event):
        self.notebook_datasets.on_menu_view_output(event)

    def on_menu_view_results(self, event):
        self.notebook_datasets.on_menu_view_results(event)

    def on_menu_view_debug(self, event):
        self.notebook_datasets.on_menu_view_debug(event)

    def on_menu_plot_x(self, event):
        self.notebook_datasets.on_menu_plot_x(event)


    ######  Help menu  ######

    def on_user_manual(self, event):
        # DEPRECATED
        path = util_misc.get_vespa_install_directory()
        path = os.path.join(path, "docs", "analysis_user_manual.pdf")
        webbrowser.open(path, 1)

    def on_analysis_online_user_manual(self, event):
        webbrowser.open("https://vespa-mrs.github.io/vespa.io/user_manuals/analysis_user_manual.html", 1)

    def on_vespa_help_online(self, event):
        webbrowser.open("https://vespa-mrs.github.io/vespa.io/", 1)

    def on_about(self, event):

        pyver = str(sys.version_info.major) + '.' + str(sys.version_info.minor)
        bit = str(8 * struct.calcsize('P')) + '-bit Python ' + pyver
        info = wx_adv.AboutDialogInfo()
        info.SetVersion(util_misc.get_vespa_version())
        info.SetCopyright("Copyright 2010, Brian J Soher. All rights reserved.")
        info.SetDescription("Analysis is an advanced MRS spectral processing and analysis environment. Running on "+bit)
        info.SetWebSite("https://github.com/vespa-mrs/vespa")
        wx_adv.AboutBox(info)



    ##############
    ##############   Public  functions  alphabetized  below
    ##############

    def update_title(self):
        """Updates the main window title to reflect the current dataset."""
        name = ""

        # Create an appropriate name for whatever is selected.
        tab = self.notebook_datasets.active_tab
        if tab and tab.dataset.dataset_filename:
            name = " - " + tab.dataset.dataset_filename
        #else:
            # If there's no active tab or the dataset_filename isn't set, we
            # don't add anything in the titlebar.
            # At present, opening a VIFF file always populates dataset_filename.
            # Importing does not.

        self.SetTitle("Analysis" + name)


    ##############
    ##############   Internal  helper  functions  alphabetized  below
    ##############

    def _import_file(self, reader, ini_name):
        datasets = [ ]

        default_path = util_analysis_config.get_path(ini_name)

        tab = self.notebook_datasets.active_tab

        if tab:     
            # there's an open dataset, let's be sure to get the first one
            # - this was issue for loading MMol basis into Fidsum metab data Notebook
            # - current solution is a hack ... make sure first dataset loaded is a 
            #    Fidsum, too, and the MMol basis data is in a middle tab.
            
            if self.notebook_datasets.GetPageCount() > 1:
                tab = self.notebook_datasets.GetPage(0)

        open_dataset = (tab.dataset if tab else None)

        msg = ""
        wx.BeginBusyCursor()

        if reader.pickfile(default_path):
            
            datasets, msg = util_file_import.get_datasets(reader, open_dataset=open_dataset)

            if msg:
                common_dialogs.message(msg, "Analysis - Open File")
            wx.EndBusyCursor()

            if datasets:
                # Successful open, write path to INI file, add dataset to Notebook
                
                path, _ = os.path.split(reader.filenames[0])
                util_analysis_config.set_path(ini_name, path)
    
                self.Freeze()
                wx.BeginBusyCursor()
                self.notebook_datasets.add_dataset_tab(datasets)
                wx.EndBusyCursor()
                self.Thaw()
                self.notebook_datasets.Layout()
                self.Layout()
        else:
            wx.EndBusyCursor()


    def _import_viff_raw_file(self):
        """
        VIFF is Vespa XML format - NB. this is NOT a processed Analysis dataset
        but only raw mrs data, so we treat it differently from ordinary VIFF.

        """
        dataset = []

        ini_name = "import_viff_raw"
        default_path = util_analysis_config.get_path(ini_name)

        filetype_filter="Raw Spectra (*.xml,*.xml.gz,*.viff,*.vif)|*.xml;*.xml.gz;*.viff;*.vif"
        filename = common_dialogs.pickfile(filetype_filter=filetype_filter,
                                            multiple=False, default_path=default_path)
        if filename:
            msg = ""
            try:
                importer = util_import.DataRawImporter(filename)
            except IOError:
                msg = """I can't read the file "%s".""" % filename
            except SyntaxError:
                msg = """The file "%s" isn't valid Vespa Interchange File Format.""" % filename

            if msg:
                common_dialogs.message(msg, "Analysis - Open File",
                                       common_dialogs.E_OK)
            else:
                # Time to rock and roll!
                wx.BeginBusyCursor()
                raws = importer.go()
                wx.EndBusyCursor()

                if raws:
                    # VIFF only has a single top level element, so keep only first element.
                    raw = raws[0]

                    if self.datasets:
                        # At least one dataset already open. New dataset attributes must match.
                        open_dataset = list(self.datasets.values())[0]
                        if (raw.dims[0] == open_dataset.raw_dims[0]) and \
                           (np.round(raw.sw,1) == np.round(open_dataset.sw,1)):
                            # All is well!
                            pass
                        else:
                            # The dimensions don't match. We can't open these files.
                            common_dialogs.message(_MSG_OPEN_ATTRIBUTE_MISMATCH,
                                                   "Analysis - Dimension Mismatch")
                            return
                    else:
                        open_dataset = None

                    # enforce multi-dataset condition: zero fill
                    if open_dataset:
                        zero_fill_multiplier = open_dataset.zero_fill_multiplier
                    else:
                        zero_fill_multiplier = 0

                    # set up user-defined changes to blocks that will be created
                    block_classes = {}
                    if raw.block_prep_flavor in list(mrs_dataset._XML_TAG_TO_SLOT_CLASS_MAP.keys()):
                        block_classes['prep'] = mrs_dataset._XML_TAG_TO_SLOT_CLASS_MAP[raw.block_prep_flavor][1]

                    dataset = mrs_dataset.dataset_from_raw(raw, block_classes, zero_fill_multiplier)

                    path, _ = os.path.split(filename)
                    util_analysis_config.set_path(ini_name, path)

                    dataset = [dataset,]
                else:
                    msg = """No Vespa raw data found in that VIFF file."""
                    common_dialogs.message(msg, "Analysis - Import Data", common_dialogs.E_OK)


            return dataset


    def _open_viff_dataset_file(self):
        """
        VIFF - Vespa Interchange File Format - in this case, this IS a file
        in the Vespa XML format for data saved from or opened up into the
        Analysis application. This is the only format that is actually 'opened'
        by Analysis, all other formats are considered 'imports'.

        Note that we only allow a single VIFF file to be opened as opposed
        to DICOM, VASF or other imported formats where we allow multiple files
        to be opened which are then concatenated into one Dataset object.

        If the open is successful (and the dimensions match to any existing
        data), the Dataset is opened into a new dataset tab.

        The Dataset object is returned (or None if the user doesn't choose a
        file), along with a list of the filenames opened.

        """
        file_type = common_constants.MrsFileTypes.VIFF
        ini_name = "open_viff"
        default_path = util_analysis_config.get_path(ini_name)

        # Note that we only allow a single VIFF file to be opened.
        filetype_filter="Spectra (*.xml,*.xml.gz,*.viff,*.vif)|*.xml;*.xml.gz;*.viff;*.vif"
        filenames = common_dialogs.pickfile(filetype_filter=filetype_filter,
                                            multiple=True,
                                            default_path=default_path)
        if filenames:
            datasets = []
            for filename in filenames:
                msg = ""
                try:
                    importer = util_import.DatasetImporter(filename)
                except IOError:
                    msg = """I can't read the file "%s".""" % filename
                except SyntaxError:
                    msg = """The file "%s" isn't valid Vespa Interchange File Format.""" % filename
    
                if msg:
                    common_dialogs.message(msg, "Analysis - Open File",
                                           common_dialogs.E_OK)
                else:
                    # Time to rock and roll!
                    wx.BeginBusyCursor()
                    dsets = importer.go()
                    wx.EndBusyCursor()
    
                    for dataset in dsets:
                        # check to ensure that none of the selected files is
                        # actually an Analysis Preset file
                        if dataset.behave_as_preset:
                            # No data in Preset file, can't load
                            common_dialogs.message(_MSG_PRESET_MISMATCH,
                                                   "Analysis - Preset Filetype Mismatch")
                            return
                    
                    for item in dsets:
                        datasets.append(item)


            if datasets:

                for dataset in datasets:
                    if self.datasets:
                        # There are one or more datasets already open. The
                        # attributes open_dataset.raw_dims of the
                        # currently open dataset(s) must match those of the
                        # dataset(s) that we're trying to open.
                        # To compare, we grab one of the currently open
                        # datasets. It doesn't matter which one since they
                        # all have matching attributes.
                        #
                        # Note. Dimensionality rules also apply to zerofill

                        open_dataset = list(self.datasets.values())[0]
                        if (dataset.raw_dims[0] == open_dataset.raw_dims[0]) and \
                           (np.round(dataset.sw,2) == np.round(open_dataset.sw,2)):
                            # All is well!
                            pass
                        else:
                            # The dimensions don't match. We can't open these files.
                            common_dialogs.message(_MSG_OPEN_ATTRIBUTE_MISMATCH, "Analysis - Dimension Mismatch")
                            return

                        open_dataset = list(self.datasets.values())[0]
                        if (dataset.spectral_dims == open_dataset.spectral_dims):
                            # All is well!
                            pass
                        else:
                            # The zerofill factors don't match. We can't open these files.
                            common_dialogs.message(_MSG_OPEN_ZEROFILL_MISMATCH, "Analysis - Dimension Mismatch")
                            return

                for dataset in datasets:
                    dataset.set_associated_datasets(datasets)
                    if dataset.id == datasets[-1].id:
                        dataset.dataset_filename = filename
                        # dataset.filename is an attribute set only at run-time
                        # to maintain the name of the VIFF file that was read in
                        # rather than deriving a filename from the raw data
                        # filenames with *.xml appended. But we need to set this
                        # filename only for the primary dataset, not the associated
                        # datasets. Associated datasets will default back to their
                        # raw filenames if we go to save them for any reason
                    else:
                        dataset.dataset_filename = ''


                dtabs = None
                if datasets:
                    dtabs = self.notebook_datasets.add_dataset_tab(datasets)

                path, _ = os.path.split(filenames[0])
                util_analysis_config.set_path(ini_name, path)

                return dtabs

            else:
                if not datasets:
                    common_dialogs.message(_MSG_NO_DATASETS_FOUND % filename, "Analysis - Open VIFF")


    def _save_viff(self, dataset):
        msg = ""
        filename = dataset.dataset_filename
        try:
            util_export.export(filename, [dataset], None, None, False)
            path, _ = os.path.split(filename)
            util_analysis_config.set_path("save_viff", path)
        except IOError:
            msg = """I can't write the file "%s".""" % filename

        if msg:
            common_dialogs.message(msg, style=common_dialogs.E_OK)
        else:
            # dataset.filename is set only at run-time to hold the name of the
            # VIFF file read in rather than being derived from the raw data
            # filename with *.xml appended. It is set here to indicate the
            # current name the dataset has been saved to VIFF format as.
            dataset.dataset_filename = filename

        self.update_title()




#--------------------------------------------------------------------

if __name__ == "__main__":

    # Having a function for app init is handy for profiling with cProfile

    app, db_path = util_init.init_app("Analysis")

    if util_misc.get_platform() == "windows":
        # code to allow Windows to set Windows taskbar icon correctly
        # - without this, Win was using default icon for python.exe
        # - bugfix from, https://stackoverflow.com/questions/1551605/how-to-set-applications-taskbar-icon-in-windows-7/1552105#1552105

        import ctypes
        # myappid is an arbitrary string - but different for each Vespa App
        myappid = u'mycompany.myproduct.subproduct.analysis' 
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)    


    import vespa.analysis.util_db as util_db

    db = util_db.Database(db_path, True)

    # My settings are in simulation.ini
    config = util_analysis_config.Config()

    position, size = config.get_window_coordinates("main")
    frame = Main(db, position, size)
    app.SetTopWindow(frame)
    frame.Show()

#    import cProfile
    app.MainLoop()
#    cProfile.run('app.MainLoop()')

#    to time this call user
#    >python -m cProfile -s cumulative bloch_lib_hargreaves.py
