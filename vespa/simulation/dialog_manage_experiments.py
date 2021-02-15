# Python modules

import os

# 3rd party modules
import wx
import numpy

# Our modules
import vespa.simulation.auto_gui.manage_experiments as gui_manage_experiments
import vespa.common.util.time_ as util_time
import vespa.common.util.import_ as util_import
import vespa.common.util.export as util_export
import vespa.common.wx_gravy.common_dialogs as common_dialogs
import vespa.common.constants as common_constants
import vespa.common.wx_gravy.util as wx_util
import vespa.common.dialog_export as dialog_export

#------------------------------------------------------------------------------
# Note. GUI Architecture/Style
#
# Many of the GUI components in Vespa are designed using the WxGlade 
# application to speed up development times. The GUI components are designed
# interactively and users can preview the resultant window/panel/dialog, but
# while event functions can be specified, only stub functions with those 
# names are created. The WxGlade files (with *.wxg extensions) are stored in
# the 'wxglade' subdirectory. The ouput of their code generation are stored
# in the 'auto_gui' subdirectory. 
#
# To used these GUI classes, each one is inherited into a unique 'vespa'
# class, where program specific initialization and other functionality are
# written.  Also, the original stub functions for widget event handlers are 
# overloaded to provide program specific event handling.
#------------------------------------------------------------------------------


class DialogManageExperiments(gui_manage_experiments.MyDialog):
    """Displays the dialog for experiment managment (view, delete, import
    export, etc.). The parent param is for the parent window and may be
    None. The db param must be a Database instance.
    """
    def __init__(self, parent, db):
        if not parent:
            parent = wx.GetApp().GetTopWindow()

        gui_manage_experiments.MyDialog.__init__(self, parent)

        self.db = db
        
        # self.experiments is a list of experiment preview objects; there's a 
        # 1:1 correspondence between items in this list and items in the
        # dialog's listbox.
        self.experiments = ( )

        # Windows & OS X make comboboxes no larger than the longest item
        # of content, and for these very short isotope strings, that's fine.
        # GTK, however, uses some large-ish default size in this context and
        # that has the combbox overwriting the out-of-service checkbox
        # to the right. This code fixes that problem.
        if wx.Platform == "__WXGTK__":
            size = (100, -1)
            self.ComboIsotope.SetMinSize(size)
            self.ComboIsotope.SetSize(size)
            self.ComboB0.SetMinSize(size)
            self.ComboB0.SetSize(size)

        # self.isotope is the currently selected isotope by which to filter.
        # It's None (==> no filtering criteria) if "any" is selected in
        # the combobox.
        self.isotope = None

        # self.b0 is the currently selected B0 by which to filter.
        # It's None (==> no filtering criteria) if "any" is selected in
        # the combobox.
        self.b0 = None

        self.ComboB0.Append("any")

        b0s = [str(b0) for b0 in self.db.fetch_b0_bins()]
        self.ComboB0.AppendItems(b0s)
        
        self.ComboB0.SetSelection(0)

        self.ListExperiments.InsertColumn(0, "Name")
        self.ListExperiments.InsertColumn(1, "Public", wx.LIST_FORMAT_CENTER)

        self.populate_isotopes()
        self.populate_experiments()

        self.SetSize( (500, 400) )

        self.Layout()

        self.Center()

        # Set focus on the list, otherwise it'll be on the isotope combobox.
        self.ListExperiments.SetFocus()

        # Grab key events so I can look for Ctrl/Cmd+A ==> select all
        self.ListExperiments.Bind(wx.EVT_KEY_DOWN, self.on_key_down)

        self.set_button_status()


    ##################    EVENT HANDLERS     ##################
    
    def on_key_down(self, event):
        if wx_util.is_select_all(event):
            wx_util.select_list_ctrl_items(self.ListExperiments, 
                                   list(range(self.ListExperiments.GetItemCount())))
            
            self.on_experiment_selection_change()
            
        event.Skip()                    


    def on_view(self, event):
        experiments = self.get_selected_experiments()

        if len(experiments):
            self.view_experiment(experiments[0])


    def on_delete(self, event):
        trouble = self.get_selected_open_intersection()
        
        if trouble:
            msg = "The following %d experiment(s) are currently open. "    \
                  "Please close them before deleting them.\n\n%s"
            names = [experiment.name for experiment in trouble]
            msg = msg % (len(names), ", ".join(names))
                    
            common_dialogs.message(msg, "Delete Experiment(s)",
                                   common_dialogs.X_OK)
        else:
            # no trouble
            experiments = self.get_selected_experiments()

            if experiments:
                names = [experiment.name for experiment in experiments]
                msg = "Are you sure you want to delete the following %d "   \
                      "experiment(s)?\n\n%s"
                msg = msg % (len(names), ", ".join(names))
                if wx.YES == common_dialogs.message(msg,
                                                    "Delete Experiment(s)",
                                                    common_dialogs.Q_YES_NO):
                    # See ya!
                    ids = [experiment.id for experiment in experiments]
                    self.db.delete_experiments(ids)
                    self.populate_experiments()
                    self.ListExperiments.SetFocus()
            else:
                common_dialogs.message("Please select one or more experiments to delete.")


    def on_import(self, event):
        filename = common_dialogs.pickfile("Select Import File")
        
        if filename:
            msg = ""
            try:
                importer = util_import.ExperimentImporter(filename, self.db)
            except IOError:
                msg = """I can't read the file "%s".""" % filename
            except SyntaxError:
                msg = """The file "%s" isn't a valid Simulation export file.""" % filename
                
            if msg:
                common_dialogs.message(msg, style=common_dialogs.E_OK)
            else:
                # Time to rock and roll!
                wx.BeginBusyCursor()
                importer.go()
                wx.EndBusyCursor()

                self.populate_experiments()

                self.ListExperiments.SetFocus()
                
                msg = """Simulation imported %d of %d experiment(s).

Would you like to see the import log?"""            
                msg = msg % (len(importer.imported), importer.found_count)
                if wx.YES == common_dialogs.message(msg, "Import Complete", 
                                                    common_dialogs.Q_YES_NO):
                    wx_util.display_file(importer.log_filename)


    def on_export(self, event):
        trouble = self.get_selected_open_intersection()
        
        if trouble:
            msg = "The following %d experiment(s) are currently open. "    \
                  "Please close them before exportig them.\n\n%s"
            names = [experiment.name for experiment in trouble]
            msg = msg % (len(names), ", ".join(names))
                    
            common_dialogs.message(msg, "Export Experiment(s)",
                                   common_dialogs.X_OK)
        else:
            # no trouble
            experiments = self.get_selected_experiments()
        
            if experiments:
                # Remember that these are lightweight experiment preview 
                # objects and not fully fledged experiment objects. I have to 
                # load proper copies of the experiments before invoking the 
                # exporter.
                experiments = [self.db.fetch_experiment(e.id) for e in experiments]
            
                # It's optional to to export the experiment's simulations; see
                # whether or not the user wants to do so.
                if len(experiments) == 1:
                    msg = "When this experiment is"
                else:
                    msg = "When these %d experiments are" % len(experiments)
                
                msg += " exported, do you want to include the experiment's" \
                       " results in addition to the experiment's definition?" 
                if wx.NO == common_dialogs.message(msg, "Export Experiments", 
                                                   common_dialogs.Q_YES_NO):
                    # Delete the experiments' results before passing them to
                    # the exporter.
                    for experiment in experiments:
                        empty = numpy.array([])
                        for simulation in experiment.simulations:
                            simulation.started = None
                            simulation.completed = None
                            simulation.ppms = empty
                            simulation.areas = empty
                            simulation.phases = empty
            
                dialog = dialog_export.DialogExport(self)
                dialog.ShowModal()

                if dialog.export:
                    filename = dialog.filename
                    comment = dialog.comment
                    compress = dialog.compress
                    
                    wx.BeginBusyCursor()
                    try:
                        util_export.export(filename, experiments, self.db,
                                           comment, compress)
                    except IOError as xxx_todo_changeme:
                        (error_number, error_string) = xxx_todo_changeme.args
                        msg = """Exporting to "%s" failed. The operating system message is below --\n\n""" % filename
                        msg += error_string
                        common_dialogs.message(msg, "Vespa Export", 
                                               common_dialogs.E_OK)
                
                    # Reload list in case any pulse_projects are newly public.
                    self.populate_experiments()

                    wx.EndBusyCursor()
                #else:
                    # The user cancelled the export dialog.

                dialog.Destroy()

                self.ListExperiments.SetFocus()
            else:
                common_dialogs.message("Please select one or more experiments to export.")


    def on_close(self, event):
        self.Close()


    def on_isotope(self, event):
        self.isotope = None if (event.GetSelection() == 0) else event.GetString()

        self.populate_experiments()


    def on_b0(self, event):
        self.b0 = None if (event.GetSelection() == 0) else event.GetString()

        self.populate_experiments()


    def on_experiment_selection_change(self, event=None):
        # This function is sometimes called internally (rather than by wx)
        # in which case event is None. 
        self.set_button_status()


    def on_experiment_activated(self, event):
        if self.ButtonView.IsEnabled():
            experiments = self.get_selected_experiments()

            if len(experiments):
                self.view_experiment(experiments[0])

    ##################    Internal helper functions     ##################
    ##################      in alphabetical order       ##################
    
    def get_selected_experiments(self):
        """Returns a (possibly empty) list of the currently selected
        experiments.
        """
        indices = wx_util.get_selected_item_indices(self.ListExperiments)
        
        return [experiment for i, experiment in enumerate(self.experiments)
                                             if (i in indices)]
        
        
    def get_selected_open_intersection(self):
        """Returns the intersection between the set of currently selected
        projects and the projects open in tabs. This is useful because some
        operations (e.g. delete) can't be performed when this set is not 
        empty.
        """
        open_experiments = wx.GetApp().vespa.get_open_experiments()
        
        selected_experiments = self.get_selected_experiments()
        
        # I can't use set() to construct the intersection because these might
        # be (and probably will be) different objects representing the same
        # experiment.
        open_ids = [experiment.id for experiment in open_experiments]
        
        return [experiment for experiment in selected_experiments
                           if experiment.id in open_ids]


    def populate_experiments(self):
        """Clears & repopulates the list of experiments. Selected items are
        reselected, if possible. 
        
        For purposes of reselecting, items are identified by id rather than
        by list position. This is useful e.g. when editing because an item
        that gets renamed will probably move to a very different place in 
        the list."""
        # Build a list of the ids of the currently selected items. I'll
        # use this after re-populating the list.
        selected = [experiment.id for experiment 
                                  in self.get_selected_experiments()]

        self.ListExperiments.DeleteAllItems()

        self.experiments = \
                    self.db.fetch_experiment_previews(self.b0, self.isotope)

        # Populate the listbox
        frozen_color = wx.Colour(*common_constants.FROZEN_COLOR)
        for i, experiment in enumerate(self.experiments):
            self.ListExperiments.InsertItem(i, experiment.name)

            # Mark the public column if necessary
            public = "x" if experiment.is_public else " "
            self.ListExperiments.SetItem(i, 1, public)

            if experiment.is_public:
                # Frozen!
                self.ListExperiments.SetItemBackgroundColour(i, frozen_color)

        self.ListExperiments.SetColumnWidth(0, wx.LIST_AUTOSIZE)
        self.ListExperiments.SetColumnWidth(1, wx.LIST_AUTOSIZE_USEHEADER)

        # Reselect all the items that were selected before, if possible.
        if self.experiments:
            if selected:
                self.select_ids(selected)

            if not self.get_selected_experiments():            
                # Nothing is selected, so I select the first item.
                wx_util.select_list_ctrl_items(self.ListExperiments, 0)

        self.on_experiment_selection_change()


    def populate_isotopes(self):
        """Clears and populates the list of isotopes"""
        isotopes = ["any"]

        isotopes += self.db.fetch_isotopes()

        self.ComboIsotope.SetItems(isotopes)

        if self.isotope:
            index = self.ComboIsotope.FindString(self.isotope)

            if index == wx.NOT_FOUND:
                index = 0
        else:
            index = 0

        self.ComboIsotope.SetSelection(index)

        self.isotope = self.ComboIsotope.GetStringSelection()

        if self.isotope == "any":
            self.isotope = None

    
    def select_ids(self, ids):
        """Attempts to select the list indices that correspond to the
        uuids in the ids param. Any non-existent ids are ignored.
        """
        indices = [ ]
        
        indices = [i for i, experiment in enumerate(self.experiments)
                                       if experiment.id in ids]
        
        wx_util.select_list_ctrl_items(self.ListExperiments, indices)
        
        if indices:
            self.ListExperiments.EnsureVisible(indices[0])
        
                
    def set_button_status(self):
        """Enables/disables buttons based on what is selected in the list"""
        experiments = self.get_selected_experiments()
        
        enable = (len(experiments) == 1)
        self.ButtonView.Enable(enable)

        enable = bool(len(experiments))
        self.ButtonDelete.Enable(enable)
        self.ButtonExport.Enable(enable)


    def view_experiment(self, experiment):            
    
        lines = str(experiment)
        lines += "\n\nSimulation Results\n" + ("-" * 75) + "\n"
 
        lines += "  For a complete view of the experiment and all \n"+ \
                 "  Simulation results, please open the experiment \n"+ \
                 "  and use the View->Output->Text menu item."
        
        wx_util.display_text_as_file(lines)
