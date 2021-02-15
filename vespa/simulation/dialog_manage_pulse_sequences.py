# Python modules

import os
import copy

# 3rd party modules
import wx

# Our modules
import vespa.simulation.auto_gui.manage_pulse_sequences as gui_manage_pulse_sequences
import vespa.simulation.dialog_pulse_sequence_info as dialog_pulse_sequence_info
import vespa.simulation.dialog_pulse_sequence_editor as dialog_pulse_sequence_editor
import vespa.common.util.db as util_db
import vespa.common.util.time_ as util_time
import vespa.common.util.import_ as util_import
import vespa.common.util.export as util_export
import vespa.common.mrs_pulse_sequence as mrs_pulse_sequence
import vespa.common.wx_gravy.util as wx_util
import vespa.common.wx_gravy.common_dialogs as common_dialogs
import vespa.common.constants as common_constants
import vespa.common.dialog_export as dialog_export

PULSE_SEQUENCE_IN_USE_MESSAGE = "Sorry, no pulse sequences have been deleted \
because one or more of the pulse sequences you selected is used in an \
experiment.\n\nIn-use pulse sequences may not be deleted."

MODIFICATIONS_LIMITED_MESSAGE = "Sorry, pulse sequences may not be edited \
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


class DialogManagePulseSequences(gui_manage_pulse_sequences.MyDialog):
    """Displays the dialog for pulse sequence managment (new, edit, delete,
    import export, etc.). The parent param is for the parent window and may
    be None. The db param must be a Database instance.
    """
    def __init__(self, parent, db):
        if not parent:
            parent = wx.GetApp().GetTopWindow()

        gui_manage_pulse_sequences.MyDialog.__init__(self, parent)

        self.db = db

        # self.pulse_sequences is a list of pulse sequence instances; there's
        # a 1:1 correspondence between items in this list and in the
        # dialog's listbox.
        self.pulse_sequences = ( )

        self.ListPulseSequences.InsertColumn(0, "Name")
        self.ListPulseSequences.InsertColumn(1, "Public", wx.LIST_FORMAT_CENTER)
        self.ListPulseSequences.InsertColumn(2, "Use Count", wx.LIST_FORMAT_CENTER)

        self.populate_pulse_sequences()
        

        if self.pulse_sequences:
            # I want to make the list tall enough to display fifteen (or so)
            # items. I calculate the height of an item and use that to figure
            # out how much extra the other parts of the list occupy.
            ITEMS_DISPLAYED = 15

            item_id = self.ListPulseSequences.GetItem(0).GetId()

            rect = self.ListPulseSequences.GetItemRect(item_id)

            item_height = rect.height

            list_width, list_height = self.ListPulseSequences.GetSize()

            extra = list_height % item_height

            list_height = item_height * ITEMS_DISPLAYED
            list_height += extra

            self.ListPulseSequences.SetMinSize( (list_width, list_height) )
        #else:
            # The list is empty; no point in futzing with the size.

        self.Layout()

        self.Fit()
        
        width, height = self.GetSize()
        
        self.SetSize( (max(width, 500), height) )

        self.Center()

        # Set focus on the list, otherwise it'll be somewhere less useful.
        self.ListPulseSequences.SetFocus()

        # Grab key events so I can look for Ctrl/Cmd+A ==> select all
        self.ListPulseSequences.Bind(wx.EVT_KEY_DOWN, self.on_key_down)

        self.set_button_status()


    ##################    EVENT HANDLERS     ##################

    def on_key_down(self, event):
        if wx_util.is_select_all(event):
            wx_util.select_list_ctrl_items(self.ListPulseSequences)
                    
            self.on_pulse_sequence_selection_change()
            
        event.Skip()


    def on_new(self, event):

        dialog = dialog_pulse_sequence_editor.DialogPulseSequenceEditor(self,
                                                                        self.db)
        if dialog.ShowModal() == wx.ID_OK:
            self.db.insert_pulse_sequence(dialog.pulse_sequence)
            
            self.populate_pulse_sequences()
            # Deselect all items in the list.
            wx_util.select_list_ctrl_items(self.ListPulseSequences, select=False)
            # Select the new item.
            self.select_ids( (dialog.pulse_sequence.id, ) )
            
        dialog.Destroy()
        self.ListPulseSequences.SetFocus()


    def on_edit(self, event):
        # Even if multiple pulse sequences are selected, we only edit the first
        pulse_sequences = self.get_selected_pulse_sequences()

        if len(pulse_sequences):
            self.edit_pulse_sequence(pulse_sequences[0])


    def on_view(self, event):
        pulse_sequences = self.get_selected_pulse_sequences()

        if len(pulse_sequences):
            dialog = \
                dialog_pulse_sequence_info.DialogPulseSequenceInfo(self,
                                                                   self.db,
                                                            pulse_sequences[0],
                                                                   True)
            dialog.ShowModal()
            dialog.Destroy()
            self.ListPulseSequences.SetFocus()


    def on_clone(self, event):
        pulse_sequences = self.get_selected_pulse_sequences()

        for pulse_sequence in pulse_sequences:
            # Create the clone
            original = pulse_sequence
            pulse_sequence = pulse_sequence.clone()
            
            # Generate a unique name for it.
            pulse_sequence.name = self.db.find_unique_name(pulse_sequence,
                                                           "clone")

            # Append a comment marking the clone
            comment = "Cloned %s from %s (%s)\n" % \
                        (util_time.now(util_time.DISPLAY_TIMESTAMP_FORMAT),
                         original.id, original.name)
            if pulse_sequence.comment:
                pulse_sequence.comment += "\n" + comment
            else:
                pulse_sequence.comment = comment

            self.db.insert_pulse_sequence(pulse_sequence)

        if len(pulse_sequences):
            self.populate_pulse_sequences()
            self.ListPulseSequences.SetFocus()


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
            pulse_sequences = self.get_selected_pulse_sequences()

            if pulse_sequences:
                in_use = [bool(pulse_sequence.experiment_names) for pulse_sequence 
                                                                in pulse_sequences]
                if any(in_use):
                    msg = PULSE_SEQUENCE_IN_USE_MESSAGE
            else:
                msg = "Please select one or more pulse sequences to delete."
                
        if msg:
            common_dialogs.message(msg, "Delete Pulse Sequence(s)")
        else:
            # All is well, deletion can proceed. Confirm with user.
            names = [pulse_sequence.name for pulse_sequence in pulse_sequences]
            msg = "Are you sure you want to delete the following %d pulse_sequence(s)?\n\n%s"
            msg = msg % (len(names), ", ".join(names))
            if wx.YES == common_dialogs.message(msg,
                                                "Delete Pulse Sequence(s)",
                                                common_dialogs.Q_YES_NO):
                # See ya!
                ids = [pulse_sequence.id for pulse_sequence in pulse_sequences]
                self.db.delete_pulse_sequences(ids)
                self.populate_pulse_sequences()
                self.ListPulseSequences.SetFocus()


    def on_import(self, event):
        filename = common_dialogs.pickfile("Select Import File")
        
        if filename:
            msg = ""
            try:
                importer = util_import.PulseSequenceImporter(filename, self.db)
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

                self.populate_pulse_sequences()
                wx.EndBusyCursor()

                self.ListPulseSequences.SetFocus()
                
                msg = """Simulation imported %d of %d pulse sequence(s).

Would you like to see the import log?"""            
                msg = msg % (len(importer.imported), importer.found_count)
                if wx.YES == common_dialogs.message(msg, "Import Complete", 
                                                    common_dialogs.Q_YES_NO):
                    wx_util.display_file(importer.log_filename)


    def on_export(self, event):
        pulse_sequences = self.get_selected_pulse_sequences()

        if pulse_sequences:
            dialog = dialog_export.DialogExport(self)
            dialog.ShowModal()

            if dialog.export:
                filename = dialog.filename
                comment = dialog.comment
                compress = dialog.compress
            
                wx.BeginBusyCursor()
                try:
                    util_export.export(filename, pulse_sequences, self.db,
                                       comment, compress)
                except IOError as xxx_todo_changeme:
                    (error_number, error_string) = xxx_todo_changeme.args
                    msg = """Exporting to "%s" failed. The operating system message is below --\n\n""" % filename
                    msg += error_string
                    common_dialogs.message(msg, "Vespa Export", 
                                           common_dialogs.E_OK)
        
                # Reload list in case any pulse_projects are newly public.
                self.populate_pulse_sequences()

                wx.EndBusyCursor()
            #else:
                # The user cancelled the export dialog.

            dialog.Destroy()
            
            self.ListPulseSequences.SetFocus()
        else:
            common_dialogs.message("Please select one or more pulse sequences to export.")


    def on_close(self, event):
        self.Close()


    def on_pulse_sequence_selection_change(self, event=None):
        # This function is sometimes called internally (rather than by wx)
        # in which case event is None. 
        self.set_button_status()


    def on_pulse_sequence_activated(self, event):
        if self.ButtonEdit.IsEnabled():
            self.edit_pulse_sequence(self.pulse_sequences[event.Index])


    ##################    Internal helper functions     ##################
    ##################      in alphabetical order       ##################

    def edit_pulse_sequence(self, pulse_sequence):
        """displays a dialog for editing the pulse sequence"""
        # We don't allow the user to edit pulse seqs when an experiment is
        # being edited. It makes the experiment GUI code too complicated.
        experiments = wx.GetApp().vespa.get_open_experiments()
        experiments = [experiment for experiment in experiments 
                                                 if not experiment.is_public]
                                                 
        if experiments:
            common_dialogs.message(MODIFICATIONS_LIMITED_MESSAGE)
        else:
            # Make a copy of the pulse seq so that the edit dialog can make 
            # changes to it without mucking up my copy.
            pulse_sequence = copy.deepcopy(pulse_sequence)
        
            dialog = dialog_pulse_sequence_editor.DialogPulseSequenceEditor(self,
                                                                         self.db,
                                                                  pulse_sequence)
            if dialog.ShowModal() == wx.ID_OK:
                self.db.update_pulse_sequence(pulse_sequence)

                self.populate_pulse_sequences()

            dialog.Destroy()
            
        self.ListPulseSequences.SetFocus()


    def get_selected_pulse_sequences(self):
        """Returns a (possibly empty) list of the currently selected
        pulse sequences.
        """
        indices = wx_util.get_selected_item_indices(self.ListPulseSequences)
        
        return [pulse_sequence for i, pulse_sequence 
                               in enumerate(self.pulse_sequences)
                               if (i in indices)]


    def populate_pulse_sequences(self):
        """Clears & repopulates the list. Selected items are
        reselected, if possible. 
        
        For purposes of reselecting, items are identified by id rather than
        by list position. This is useful e.g. when editing because an item
        that gets renamed will probably move to a very different place in 
        the list.
        """
        # Build a list of the ids of the currently selected items. I'll
        # use this after re-populating the list.
        selected = [pulse_sequence.id for pulse_sequence 
                                      in self.get_selected_pulse_sequences()]

        self.ListPulseSequences.DeleteAllItems()

        self.pulse_sequences = self.db.fetch_pulse_sequences(True)

        # Populate the listbox
        frozen_color = wx.Colour(*common_constants.FROZEN_COLOR)
        
        for i, pulse_sequence in enumerate(self.pulse_sequences):
            self.ListPulseSequences.InsertItem(i, pulse_sequence.name)

            # Mark the public column if necessary
            public = "x" if pulse_sequence.is_public else " "
            self.ListPulseSequences.SetItem(i, 1, public)

            # I only show a number for the usage count if it's non-zero.
            # It's easier to read that way.
            usage_count = " "
            if len(pulse_sequence.experiment_names):
                usage_count = str(len(pulse_sequence.experiment_names))
            self.ListPulseSequences.SetItem(i, 2, usage_count)

            if pulse_sequence.is_frozen:
                # Frozen!
                self.ListPulseSequences.SetItemBackgroundColour(i, frozen_color)

        self.ListPulseSequences.SetColumnWidth(0, wx.LIST_AUTOSIZE)
        self.ListPulseSequences.SetColumnWidth(1, wx.LIST_AUTOSIZE_USEHEADER)
        self.ListPulseSequences.SetColumnWidth(2, wx.LIST_AUTOSIZE_USEHEADER)

        # Reselect all the items that were selected before, if possible.
        if self.pulse_sequences:
            if selected:
                self.select_ids(selected)

            if not self.get_selected_pulse_sequences():            
                # Nothing is selected, so I select the first item.
                wx_util.select_list_ctrl_items(self.ListPulseSequences, 0)

        self.on_pulse_sequence_selection_change()


    def select_ids(self, ids):
        """Attempts to select the list indices that correspond to the
        uuids in the ids param. Any non-existent ids are ignored.
        """
        indices = [ ]
        
        indices = [i for i, pulse_sequence in enumerate(self.pulse_sequences)
                                       if pulse_sequence.id in ids]
        
        wx_util.select_list_ctrl_items(self.ListPulseSequences, indices)
        
        if indices:
            self.ListPulseSequences.EnsureVisible(indices[0])
        
                
    def set_button_status(self):
        """Enables/disables buttons based on what is selected in the list"""
        pulse_sequences = self.get_selected_pulse_sequences()

        enable = (len(pulse_sequences) == 1)
        self.ButtonEdit.Enable(enable)
        self.ButtonView.Enable(enable)

        enable = bool(len(pulse_sequences))
        self.ButtonDelete.Enable(enable)
        self.ButtonClone.Enable(enable)
        self.ButtonExport.Enable(enable)


