# Python modules

import os.path
import copy

# 3rd party modules
import wx

# Our modules
import vespa.simulation.pane_metabolite_info as pane_metabolite_info
import vespa.simulation.auto_gui.manage_metabolites as gui_manage_metabolites
import vespa.common.util.time_ as util_time
import vespa.common.util.import_ as util_import
import vespa.common.util.export as util_export
import vespa.common.mrs_metabolite as mrs_metabolite
import vespa.common.wx_gravy.common_dialogs as common_dialogs
import vespa.common.wx_gravy.util as wx_util
import vespa.common.dialog_export as dialog_export
import vespa.common.constants as common_constants

METABOLITE_IN_USE_MESSAGE = "Sorry, no metabolites have been deleted \
because one or more of the metabolites you selected are used in an \
experiment.\n\nIn-use metabolites may not be deleted."

MODIFICATIONS_LIMITED_MESSAGE = "Sorry, metabolites may not be edited \
nor deleted while you are editing an experiment."


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


class DialogManageMetabolites(gui_manage_metabolites.MyDialog):
    """Displays the dialog for metab managment (new, edit, delete, import
    export, etc.). The parent param is for the parent window and may be
    None. The db param must be a Database instance.
    """
    def __init__(self, parent, db):
        if not parent:
            parent = wx.GetApp().GetTopWindow()

        gui_manage_metabolites.MyDialog.__init__(self, parent)

        self.db = db

        # self.metabolites is a list of metabolite objects; there's a 1:1
        # correspondence between items in this list and metabs in the
        # dialog's listbox.
        self.metabolites = ( )

        # Windows & OS X make comboboxes no larger than the longest item
        # of content, and for these very short isotope strings, that's fine.
        # GTK, however, uses some large-ish default size in this context and
        # that has the combbox overwriting the out-of-service checkbox
        # to the right. This code fixes that problem.
        if wx.Platform == "__WXGTK__":
            size = (100, -1)
            self.ComboIsotope.SetMinSize(size)
            self.ComboIsotope.SetSize(size)

        self.show_deactivated = True

        self.CheckboxShowDeactivated.SetValue(self.show_deactivated)

        self.ListMetabolites.InsertColumn(0, "Name")
        self.ListMetabolites.InsertColumn(1, "Public", wx.LIST_FORMAT_CENTER)
        self.ListMetabolites.InsertColumn(2, "Use Count", wx.LIST_FORMAT_CENTER)

        # self.isotope is the currently selected isotope by which to filter.
        # It's None (==> no filtering criteria) if "all" is selected in
        # the combobox.
        self.isotope = None

        self.populate_isotopes()
        self.populate_metabolites()

        self.SetSize( (600, 400) )
        
        # The Deactivate button tends to be bigger than all the others 
        # due to its long label. Make them consistent.
        size = self.ButtonDeactivate.GetSize()
        for button in (self.ButtonNew, self.ButtonEdit, self.ButtonView,
                       self.ButtonClone, self.ButtonDelete, self.ButtonImport,
                       self.ButtonExport, self.ButtonClose):
            button.SetMinSize(size)

        self.Layout()

        self.Center()

        # Set focus on the list, otherwise it'll be on the isotope combobox.
        self.ListMetabolites.SetFocus()

        # Grab key events so I can look for Ctrl/Cmd+A ==> select all
        self.ListMetabolites.Bind(wx.EVT_KEY_DOWN, self.on_key_down)

        self.set_button_status()


    ##################    EVENT HANDLERS     ##################
    
    def on_key_down(self, event):
        if wx_util.is_select_all(event):
            wx_util.select_list_ctrl_items(self.ListMetabolites, 
                                   list(range(self.ListMetabolites.GetItemCount())))
                
            self.on_metabolite_selection_change()

        event.Skip()


    def on_new(self, event):
        dialog = pane_metabolite_info.DialogMetaboliteInfo(self, self.db)
        if dialog.ShowModal() == wx.ID_OK:
            self.db.insert_metabolite(dialog.pane.metabolite)
            self.populate_metabolites()
            
            # Deselect all items in the list.
            wx_util.select_list_ctrl_items(self.ListMetabolites, select=False)
            # Select the new item.
            self.select_ids( (dialog.pane.metabolite.id, ) )

        dialog.Destroy()

        self.ListMetabolites.SetFocus()


    def on_edit(self, event):
        # Even if multiple metabolites are selected, we only edit the first
        metabolites = self.get_selected_metabolites()
    
        if len(metabolites):
            self.edit_metabolite(metabolites[0])


    def on_view(self, event):
        metabolites = self.get_selected_metabolites()

        if len(metabolites):
            metabolite = metabolites[0]

            dialog = \
                pane_metabolite_info.DialogMetaboliteInfo(self, self.db, 
                                                          metabolite, True)
            dialog.ShowModal()

            dialog.Destroy()

            self.ListMetabolites.SetFocus()


    def on_deactivate(self, event):
        metabolites = self.get_selected_metabolites()

        if len(metabolites):
            metabolite = metabolites[0]
            
            if metabolite.deactivated:
                # currently inactive, user wants it active
                comment = "Activated"
                metabolite.deactivated = None
            else:
                # currently active, user wants it out of service
                comment = "Deactivated"
                metabolite.deactivated = util_time.now()
            
            comment += " %s\n" % util_time.now(util_time.DISPLAY_TIMESTAMP_FORMAT)
            
            if metabolite.comment:
                metabolite.comment += "\n" + comment
            else:
                metabolite.comment = comment
                        
            self.db.update_metabolite(metabolite)
            self.populate_metabolites()
            self.ListMetabolites.SetFocus()


    def on_clone(self, event):
        metabolites = self.get_selected_metabolites()

        for metabolite in metabolites:
            # Create the clone
            original = metabolite
            metabolite = metabolite.clone()

            # Generate a unique name for it.
            metabolite.name = self.db.find_unique_name(metabolite, "clone")

            # Append a comment marking the clone
            comment = "Cloned %s from %s (%s)\n" % \
                        (util_time.now(util_time.DISPLAY_TIMESTAMP_FORMAT),
                         original.id, original.name)
            if metabolite.comment:
                metabolite.comment += "\n" + comment
            else:
                metabolite.comment = comment

            self.db.insert_metabolite(metabolite)

        if len(metabolites):
            self.populate_metabolites()
            self.ListMetabolites.SetFocus()


    def on_delete(self, event):
        msg = ""
        
        # We don't allow the user to delete metabs when an experiment is
        # being edited. It makes the experiment GUI code too complicated.
        experiments = wx.GetApp().vespa.get_open_experiments()
        experiments = [experiment for experiment in experiments 
                                                 if not experiment.is_public]
                          
        if experiments:
            msg = MODIFICATIONS_LIMITED_MESSAGE
            
        if not msg:
            metabolites = self.get_selected_metabolites()

            if metabolites:
                in_use = [bool(metabolite.experiment_names) for metabolite 
                                                            in metabolites]
                if any(in_use):
                    msg = METABOLITE_IN_USE_MESSAGE
            else:
                msg = "Please select one or more metabolites to delete."
                
        if msg:
            common_dialogs.message(msg, "Delete Metabolite(s)")
        else:
            # All is well, deletion can proceed. Confirm with user.
            names = [metabolite.name for metabolite in metabolites]
            msg = "Are you sure you want to delete the following %d metabolite(s)?\n\n%s"
            msg = msg % (len(names), ", ".join(names))
            if wx.YES == common_dialogs.message(msg,
                                                "Delete Metabolite(s)",
                                                common_dialogs.Q_YES_NO):
                # See ya!
                ids = [metabolite.id for metabolite in metabolites]
                self.db.delete_metabolites(ids)
                self.populate_metabolites()


    def on_import(self, event):
        filename = common_dialogs.pickfile("Select Import File")
        
        if filename:
            msg = ""
            try:
                importer = util_import.MetaboliteImporter(filename, self.db)
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

                self.populate_metabolites()
                wx.EndBusyCursor()

                self.ListMetabolites.SetFocus()
                
                msg = """Simulation imported %d of %d metabolite(s).

Would you like to see the import log?"""            
                msg = msg % (len(importer.imported), importer.found_count)
                if wx.YES == common_dialogs.message(msg, "Import Complete", 
                                                    common_dialogs.Q_YES_NO):
                    wx_util.display_file(importer.log_filename)


    def on_export(self, event):
        metabolites = self.get_selected_metabolites()
        
        if metabolites:
            dialog = dialog_export.DialogExport(self)
            dialog.ShowModal()

            export = dialog.export
            filename = dialog.filename
            comment = dialog.comment
            compress = dialog.compress
            
            dialog.Destroy()
            
            if export:
                wx.BeginBusyCursor()
                try:
                    util_export.export(filename, metabolites, self.db,
                                       comment, compress)
                except IOError as xxx_todo_changeme:
                    (error_number, error_string) = xxx_todo_changeme.args
                    msg = """Exporting to "%s" failed. The operating system message is below --\n\n""" % filename
                    msg += error_string
                    common_dialogs.message(msg, "Vespa Export", common_dialogs.E_OK)

                wx.EndBusyCursor()

                # Populate metabs in case any are newly frozen
                self.populate_metabolites()
            #else:
                # The user cancelled the export dialog.

            self.ListMetabolites.SetFocus()
        else:
            common_dialogs.message("Please select one or more metabolites to export.")


    def on_close(self, event):
        self.Close()


    def on_show_deactivated(self, event):
        self.show_deactivated = event.IsChecked()

        self.populate_metabolites()


    def on_isotope(self, event):
        self.isotope = None if (event.GetSelection() == 0) else event.GetString()

        self.populate_metabolites()


    def on_metabolite_selection_change(self, event=None):
        # This function is sometimes called internally (rather than by wx)
        # in which case event is None. 
        self.set_button_status()


    def on_metabolite_activated(self, event):
        if self.ButtonEdit.IsEnabled():
            self.edit_metabolite(self.metabolites[event.Index])


    ##################    Internal helper functions     ##################
    ##################      in alphabetical order       ##################

    def edit_metabolite(self, metabolite):
        """shows a dialog for editing the metabolite"""
        
        # We don't allow the user to edit metabs when an experiment is
        # being edited. It makes the experiment GUI code too complicated.
        experiments = wx.GetApp().vespa.get_open_experiments()
        experiments = [experiment for experiment in experiments 
                                                 if not experiment.is_public]
                                                 
        if experiments:
            common_dialogs.message(MODIFICATIONS_LIMITED_MESSAGE)
        else:
            # In case the editing dialog decides to fiddle with the metab, I
            # give it a copy of mine so it won't mess up my pristine version.
            metabolite = copy.deepcopy(metabolite)
        
            dialog = pane_metabolite_info.DialogMetaboliteInfo(self, self.db, 
                                                               metabolite)
            if dialog.ShowModal() == wx.ID_OK:
                self.db.update_metabolite(metabolite)
                self.populate_metabolites()

            dialog.Destroy()

        self.ListMetabolites.SetFocus()


    def get_selected_metabolites(self):
        """Returns a (possibly empty) list of the currently selected
        metabolites.
        """        
        indices = wx_util.get_selected_item_indices(self.ListMetabolites)
        
        return [metabolite for i, metabolite in enumerate(self.metabolites)
                                             if (i in indices)]
        

    def populate_isotopes(self):
        """Clears and populates the list of isotopes"""
        isotopes = ["all"]

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

        if self.isotope == "all":
            self.isotope = None


    def populate_metabolites(self):
        """Clears & repopulates the list of metabs. Selected items are
        reselected, if possible. 
        
        For purposes of reselecting, items are identified by id rather than
        by list position. This is useful e.g. when editing because an item
        that gets renamed will probably move to a very different place in 
        the list.
        """
        # Build a list of the ids of the currently selected items. I'll
        # use this after re-populating the list.
        selected = [metabolite.id for metabolite 
                                  in self.get_selected_metabolites()]

        self.ListMetabolites.DeleteAllItems()

        self.metabolites = self.db.fetch_metabolites(self.isotope,
                                                     self.show_deactivated,
                                                     True)
            
        # Populate the listbox
        frozen_color = wx.Colour(*common_constants.FROZEN_COLOR)
        for i, metabolite in enumerate(self.metabolites):
            # If the metab is out of service, that gets appended to the name
            name = metabolite.name

            if metabolite.deactivated:
                name += " (not active)"
                
            self.ListMetabolites.InsertItem(i, name)

            # Mark the public column if necessary
            public = "x" if metabolite.is_public else " "
            self.ListMetabolites.SetItem(i, 1, public)

            # I only show a number for the usage count if it's non-zero.
            # It's easier to read that way.
            usage_count = " "
            if len(metabolite.experiment_names):
                usage_count = str(len(metabolite.experiment_names))
            self.ListMetabolites.SetItem(i, 2, usage_count)

            if metabolite.is_frozen:
                # Frozen!
                self.ListMetabolites.SetItemBackgroundColour(i, frozen_color)

        self.ListMetabolites.SetColumnWidth(0, wx.LIST_AUTOSIZE)
        self.ListMetabolites.SetColumnWidth(1, wx.LIST_AUTOSIZE_USEHEADER)
        self.ListMetabolites.SetColumnWidth(2, wx.LIST_AUTOSIZE_USEHEADER)

        # Reselect all the items that were selected before, if possible.
        if self.metabolites:
            if selected:
                self.select_ids(selected)

            if not self.get_selected_metabolites():            
                # Nothing is selected, so I select the first item.
                wx_util.select_list_ctrl_items(self.ListMetabolites, 0)

        self.on_metabolite_selection_change()


    def select_ids(self, ids):
        """Attempts to select the list indices that correspond to the metab
        uuids in the ids param. Any non-existent ids are ignored.
        """
        indices = [ ]
        
        indices = [i for i, metabolite in enumerate(self.metabolites)
                                       if metabolite.id in ids]
        
        wx_util.select_list_ctrl_items(self.ListMetabolites, indices)
        
        if indices:
            self.ListMetabolites.EnsureVisible(indices[0])
        
                
    def set_button_status(self):
        """Enables/disables buttons based on what is selected in the list"""
        metabolites = self.get_selected_metabolites()

        enable = (len(metabolites) == 1)
        self.ButtonEdit.Enable(enable)
        self.ButtonView.Enable(enable)
        self.ButtonDeactivate.Enable(enable)

        enable = bool(len(metabolites))
        self.ButtonClone.Enable(enable)
        self.ButtonDelete.Enable(enable)
        self.ButtonExport.Enable(enable)

