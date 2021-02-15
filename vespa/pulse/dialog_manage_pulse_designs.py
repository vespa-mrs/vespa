# Python modules

import os

# 3rd party modules
import wx

# Our modules
import vespa.pulse.dialog_view_pulse_design         as dialog_view_pulse_design
import vespa.pulse.auto_gui.manage_pulse_designs    as manage_pulse_designs
import vespa.common.constants                       as common_constants
import vespa.common.dialog_export           as dialog_export
import vespa.common.util.time_              as util_time
import vespa.common.util.import_            as util_import
import vespa.common.util.export             as util_export
import vespa.common.wx_gravy.util           as common_wx_util
import vespa.common.wx_gravy.common_dialogs as common_dialogs



class DialogManagePulseDesigns(manage_pulse_designs.MyDialog):
    """
    Displays the dialog for pulse_design management (view, delete, import
    export, etc.). The parent param is for the parent window and may be
    None. The db param must be a Database instance.
    """
    def __init__(self, parent, db):
        if not parent:
            parent = wx.GetApp().GetTopWindow()

        manage_pulse_designs.MyDialog.__init__(self, parent)

        self.db = db
        
        # self.pulse_designs is a list of pulse_design preview objects; 
        # there's a 1:1 correspondence between items in this list and 
        # items in the dialog's listbox.
        self.pulse_designs = ( )

        self.ListPulseDesigns.InsertColumn(0, "Name")
        self.ListPulseDesigns.InsertColumn(1, "Public", wx.LIST_FORMAT_CENTER)
        self.ListPulseDesigns.InsertColumn(2, "Use Count", wx.LIST_FORMAT_CENTER)

        self.populate_pulse_designs()

        self.SetSize( (500, 400) )
        self.Layout()
        self.Center()

        # Grab key events so I can look for Ctrl/Cmd+A ==> select all
        self.ListPulseDesigns.Bind(wx.EVT_KEY_DOWN, self.on_key_down)

        self.set_button_status()


    ##################    EVENT HANDLERS     ##################
    
    def on_key_down(self, event):
        if common_wx_util.is_select_all(event):
            common_wx_util.select_list_ctrl_items(self.ListPulseDesigns)

            self.on_selection_changed()
            
        event.Skip()                    


    def on_view(self, event):
        pulse_designs = self.get_selected_pulse_designs()

        if len(pulse_designs):
            self.view_pulse_design(pulse_designs[0])


    def on_clone(self, event):
        pulse_designs = self.get_selected_pulse_designs()

        for pulse_design in pulse_designs:
            # Since this is just an pulse_design preview, I have to fetch the 
            # full pulse_design in order to clone it.
            original = self.db.fetch_pulse_design(pulse_design.id)

            # Create the clone
            new_pulse_design = original.clone()
            
            new_pulse_design.name = self.db.find_unique_name(new_pulse_design, "clone")
            
            # Append a comment marking the clone
            comment = "Cloned %s from %s (%s)\n" % \
                        (util_time.now(util_time.DISPLAY_TIMESTAMP_FORMAT),
                         original.id, original.name)
            if new_pulse_design.comment:
                new_pulse_design.comment += "\n" + comment
            else:
                new_pulse_design.comment = comment

            self.db.insert_pulse_design(new_pulse_design)

        if len(pulse_designs):
            self.populate_pulse_designs()


    def on_delete(self, event):
        trouble = self.get_selected_open_intersection()
        
        if trouble:
            msg = "The following %d design(s) are currently open. "    \
                  "Please close them before deleting them.\n\n%s"
            names = [design.name for design in trouble]
            msg = msg % (len(names), ", ".join(names))
                    
            common_dialogs.message(msg, "Delete Pulse Design(s)",
                                   common_dialogs.X_OK)
        else:
            pulse_designs = self.get_selected_pulse_designs()
        
            if pulse_designs:
                names = [pulse_design.name for pulse_design in pulse_designs]
                msg = "Are you sure you want to delete the following %d pulse_design(s)?\n\n%s"
                msg = msg % (len(names), ", ".join(names))
                if wx.YES == common_dialogs.message(msg,
                                                    "Delete pulse_design(s)",
                                                    common_dialogs.Q_YES_NO):
                    # Very important test here --
                    # We must not allow users to delete pulse designs that
                    # are referred to by pulse seqs. Since pulse seqs are 
                    # manipulated outside of this app (in Simulation), 
                    # changes to them can happen at any time. That's why we
                    # postpone this check until the last possible moment
                    # before deletion.
                    # Since the database isn't locked in between this check
                    # and the actual deletion, there's a window in which 
                    # a reference could be created before the deletion. 
                    # However this window is milliseconds wide and can't be
                    # exploited by a user clicking stuff with the mouse.
                    in_use = False
                    for pulse_design in pulse_designs:
                        if self.db.fetch_pulse_design_referrers(pulse_design.id):
                            in_use = True
                            break
                            
                    if in_use:
                        msg = "One or more of these designs is in use "     \
                              "(referenced by a pulse sequence).\n\n"       \
                              "Pulse designs that are in use may not "      \
                              "be deleted."
                        common_dialogs.message(msg, "Delete Pulse Design(s)")
                    else:
                        # Hasta la vista!
                        ids = [pulse_design.id for pulse_design in pulse_designs]
                        self.db.delete_pulse_designs(ids)
                        self.populate_pulse_designs()
            else:
                common_dialogs.message("Please select one or more pulse designs to delete.")


    def on_import(self, event):
        filename = common_dialogs.pickfile("Select Import File")
        
        if filename:
            msg = ""
            try:
                importer = util_import.PulseDesignImporter(filename, self.db)
            except IOError:
                msg = """I can't read the file "%s".""" % filename
            except SyntaxError:
                msg = """The file "%s" isn't a valid Vespa export file.""" % filename
                
            if msg:
                common_dialogs.message(msg, style=common_dialogs.E_OK)
            else:
                # Time to rock and roll!
                wx.BeginBusyCursor()
                importer.go()
                wx.EndBusyCursor()

                self.populate_pulse_designs()
                
                msg = """Pulse imported %d of %d pulse_design(s).
                         Would you like to see the import log?"""            
                msg = msg % (len(importer.imported), importer.found_count)
                if wx.YES == common_dialogs.message(msg, "Import Complete", 
                                                    common_dialogs.Q_YES_NO):
                    common_wx_util.display_file(importer.log_filename)


    def on_export(self, event):
        trouble = self.get_selected_open_intersection()
        
        if trouble:
            msg = "The following %d design(s) are currently open. "    \
                  "Please close them before exporting them.\n\n%s"
            names = [design.name for design in trouble]
            msg = msg % (len(names), ", ".join(names))
                    
            common_dialogs.message(msg, "Export Pulse Design(s)",
                                   common_dialogs.X_OK)
        else:
            pulse_designs = self.get_selected_pulse_designs()
        
            if pulse_designs:
                dialog = dialog_export.DialogExport(self)
                dialog.ShowModal()

                if dialog.export:
                    filename = dialog.filename
                    comment  = dialog.comment
                    compress = dialog.compress
                    
                    wx.BeginBusyCursor()
                    # Remember that these are lightweight pulse_design preview 
                    # objects and not fully fledged pulse_design objects. I have to 
                    # load proper copies of the pulse_designs before invoking the 
                    # exporter.
                    pulse_designs = [self.db.fetch_pulse_design(pp.id) for pp in pulse_designs]
                    try:
                        util_export.export(filename, pulse_designs, self.db, comment, compress)
                    except IOError as xxx_todo_changeme:
                        (error_number, error_string) = xxx_todo_changeme.args
                        msg = """Exporting to "%s" failed. The operating system message is below --\n\n""" % filename
                        msg += error_string
                        common_dialogs.message(msg, "Vespa Export", 
                                               common_dialogs.E_OK)
                
                    # Reload list in case any pulse_designs are newly public.
                    self.populate_pulse_designs()

                    wx.EndBusyCursor()
                #else:
                    # The user cancelled the export dialog.
                
                dialog.Destroy()

                self.ListPulseDesigns.SetFocus()
            else:
                common_dialogs.message("Please select one or more pulse designs to export.")


    def on_close(self, event):
        self.Close()


    def on_selection_changed(self, event=None):
        # This function is sometimes called internally (rather than by wx)
        # in which case event is None. 
        self.set_button_status()


    def on_pulse_design_activated(self, event):
        if self.ButtonView.IsEnabled():
            pulse_designs = self.get_selected_pulse_designs()

            if len(pulse_designs):
                self.view_pulse_design(pulse_designs[0])


    ##################    Internal helper functions     ##################
    
    def get_selected_open_intersection(self):
        """
        Returns the intersection between the set of currently selected
        designs and the designs open in tabs. This is useful because some
        operations (e.g. delete) can't be performed when this set is not 
        empty.
        
        """
        open_pulse_designs = wx.GetApp().vespa.get_open_designs()
        
        selected_pulse_designs = self.get_selected_pulse_designs()
        
        # I can't use set() to construct the intersection because these might
        # be (and probably will be) different objects representing the same
        # pulse design.
        open_ids = [pulse_design.id for pulse_design in open_pulse_designs]
        
        return [pulse_design for pulse_design in selected_pulse_designs
                              if pulse_design.id in open_ids]
        
    
    def get_selected_pulse_designs(self):
        """
        Returns a (possibly empty) list of the currently selected pulse_designs
        
        """
        indices = common_wx_util.get_selected_item_indices(self.ListPulseDesigns)
        
        return [pulse_design for i, pulse_design in enumerate(self.pulse_designs)
                                                   if (i in indices)]
        

    def populate_pulse_designs(self):
        """Clears & repopulates the list of pulse_designs"""
        # Build a list of the ids of the currently selected items. I'll
        # use this after re-populating the list.
        selected = [pulse_design.id for pulse_design 
                                     in  self.get_selected_pulse_designs()]

        self.ListPulseDesigns.DeleteAllItems()

        self.pulse_designs = self.db.fetch_pulse_design_previews()

        # Populate the listbox
        frozen_color = wx.Colour(*common_constants.FROZEN_COLOR)
        for i, pulse_design in enumerate(self.pulse_designs):
            self.ListPulseDesigns.InsertItem(i, pulse_design.name)

            # Mark the public column if necessary
            public = "x" if pulse_design.is_public else " "
            self.ListPulseDesigns.SetItem(i, 1, public)

            # Display referrer count if non-zero
            referrers = len(pulse_design.referrers)
            referrers = (str(referrers) if referrers else "")
            self.ListPulseDesigns.SetItem(i, 2, referrers)

            if pulse_design.is_frozen:
                # Frozen!
                self.ListPulseDesigns.SetItemBackgroundColour(i, frozen_color)

        self.ListPulseDesigns.SetColumnWidth(0, wx.LIST_AUTOSIZE)
        self.ListPulseDesigns.SetColumnWidth(1, wx.LIST_AUTOSIZE_USEHEADER)
        self.ListPulseDesigns.SetColumnWidth(2, wx.LIST_AUTOSIZE_USEHEADER)

        # Reselect all the items that were selected before, if possible.
        if self.pulse_designs:
            if selected:
                self.select_ids(selected)

            if not self.get_selected_pulse_designs():            
                # Nothing is selected, so I select the first item.
                common_wx_util.select_list_ctrl_items(self.ListPulseDesigns, 0)

        self.ListPulseDesigns.SetFocus()
        self.on_selection_changed()

    
    def select_ids(self, ids):
        """
        Attempts to select the list indices that correspond to the ids
        in the param. Any non-existent ids are ignored.
        """
        indices = [ ]
        
        indices = [i for i, pulse_design in enumerate(self.pulse_designs)
                                         if pulse_design.id in ids]
        
        common_wx_util.select_list_ctrl_items(self.ListPulseDesigns, indices)
        
        if indices:
            self.ListPulseDesigns.EnsureVisible(indices[0])
            
            
    def set_button_status(self):
        """Enables/disables buttons based on what is selected in the list"""
        pulse_designs = self.get_selected_pulse_designs()
        
        enable = (len(pulse_designs) == 1)
        self.ButtonView.Enable(enable)

        enable = bool(len(pulse_designs))
        self.ButtonClone.Enable(enable)
        self.ButtonDelete.Enable(enable)
        self.ButtonExport.Enable(enable)
        
        # Can't delete pulse designs that are in use
        if any([bool(pulse_design.referrers) for pulse_design 
                                              in pulse_designs]):
            self.ButtonDelete.Disable()


    def view_pulse_design(self, pulse_design):
        # Since this is just an pulse_design preview, I have to fetch the 
        # full pulse_design in order to display it properly.
        pulse_design = self.db.fetch_pulse_design(pulse_design.id)

        dialog = dialog_view_pulse_design.DialogViewPulseDesign(self, pulse_design)
        
        dialog.ShowModal()
        dialog.Destroy()
