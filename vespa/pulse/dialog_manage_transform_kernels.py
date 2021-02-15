# Python modules

import os

# 3rd party modules
import wx

# Our modules
import vespa.pulse.dialog_editor_transform           as dialog_editor_transform
import vespa.pulse.dialog_view_transform_kernel      as dialog_view_transform_kernel
import vespa.pulse.auto_gui.manage_transform_kernels as manage_transform_kernels
import vespa.common.util.time_                          as util_time
import vespa.common.util.import_                        as util_import
import vespa.common.util.export                         as util_export
import vespa.common.wx_gravy.common_dialogs             as common_dialogs
import vespa.common.wx_gravy.util                       as wx_util
import vespa.common.dialog_export                       as dialog_export
import vespa.common.constants                           as common_constants



class DialogManageTransformKernels(manage_transform_kernels.MyDialog):
    """
    Displays the dialog for transform_kernel management (view, delete, import
    export, etc.). The parent param is for the parent window and may be
    None. The db param must be a Database instance.
    """
    def __init__(self, parent, db):
        if not parent:
            parent = wx.GetApp().GetTopWindow()

        manage_transform_kernels.MyDialog.__init__(self, parent)

        self.db = db
        
        # this is a list of transform_kernel preview objects; there's a 
        # 1:1 correspondence between items in this list and items in the
        # dialog's listbox.
        self.transform_kernels = []
        self.all_names = []
        self.all_labels = []

        self.ListTransformKernels.InsertColumn(0, "  Name  ")
        self.ListTransformKernels.InsertColumn(1, "  Menu Label  ")
        self.ListTransformKernels.InsertColumn(2, " Public ",    wx.LIST_FORMAT_CENTER)
        self.ListTransformKernels.InsertColumn(3, " Use Count ", wx.LIST_FORMAT_CENTER)

        self.CheckHideDeprecated.SetValue(True)         # default value

        self.populate_transform_kernels()

        self.SetSize( (600, 500) )
        self.Layout()
        self.Center()

        # Grab key events so I can look for Ctrl/Cmd+A ==> select all
        self.ListTransformKernels.Bind(wx.EVT_KEY_DOWN, self.on_key_down)

        self.set_button_status()


    ##################    EVENT HANDLERS     ##################
    
    def on_key_down(self, event):
        if wx_util.is_select_all(event):
            wx_util.select_list_ctrl_items(self.ListTransformKernels)
            
            self.on_selection_changed()
            
        event.Skip()                    

    def on_new(self, event):

        # this editor uses internal default values for machine_specs and global
        # parameters, so kernel may work differently in main application
                
        dialog = dialog_editor_transform.DialogEditorTransform(self, 
                                                               self.all_names, 
                                                               self.all_labels)
        if dialog.ShowModal() == wx.ID_OK:

            self.db.insert_transform_kernel(dialog.transform_kernel)
            
            # Re-populate list, de-select all items, and select new item
            self.populate_transform_kernels()
            wx_util.select_list_ctrl_items(self.ListTransformKernels, select=False)
            self.select_ids( (dialog.transform_kernel.id, ) )
            
        dialog.Destroy()
        self.ListTransformKernels.SetFocus()
    
    
    def on_edit(self, event):
        transform_kernels = self.get_selected_transform_kernels()

        if len(transform_kernels):
            kernel = transform_kernels[0]
            
            # exclude current kernel from these lists so we don't raise
            # an exception for the name/label it currently uses.
            all_names  = [item for item in self.all_names  if item != kernel.name]
            all_labels = [item for item in self.all_labels if item != kernel.menu_label]

            # this editor uses internal default values for machine_specs and global
            # parameters, so kernel may work differently in main application

            dialog = dialog_editor_transform.DialogEditorTransform(self, 
                                                                   all_names, 
                                                                   all_labels, 
                                                                   kernel)
            if dialog.ShowModal() == wx.ID_OK:
                self.db.replace_transform_kernel(dialog.transform_kernel)
                
                # Re-populate list, de-select all items, and select new item
                self.populate_transform_kernels()
                wx_util.select_list_ctrl_items(self.ListTransformKernels, select=False)
                self.select_ids( (dialog.transform_kernel.id, ) )
                
            dialog.Destroy()
            self.ListTransformKernels.SetFocus()
            

    def on_view(self, event):
        kernels = self.get_selected_transform_kernels()

        if len(kernels):
            self.view_transform_kernel(kernels[0])


    def on_clone(self, event):
        kernels = self.get_selected_transform_kernels()

        for kernel in kernels:
            # Since this is just an transform_kernel preview, I have to fetch the 
            # full transform_kernel in order to clone it.
            original = self.db.fetch_transform_kernel(kernel.id)

            # Create the clone
            new_transform_kernel = original.clone()
            
            new_transform_kernel.name = self.db.find_unique_name(new_transform_kernel, "clone")
            
            extra = new_transform_kernel.name[len(original.name):]
            new_transform_kernel.menu_label += extra
            
            # Append a comment marking the clone
            comment = "Cloned %s from %s (%s)\n" % \
                        (util_time.now(util_time.DISPLAY_TIMESTAMP_FORMAT),
                         original.id, original.name)
            if new_transform_kernel.comment:
                new_transform_kernel.comment += "\n" + comment
            else:
                new_transform_kernel.comment = comment

            self.db.insert_transform_kernel(new_transform_kernel)

        if len(kernels):
            self.populate_transform_kernels()


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
            kernels = self.get_selected_transform_kernels()
        
            if kernels:
                names = [kernel.name for kernel in kernels]
                msg = "Are you sure you want to delete the following %d transform_kernel(s)?\n\n%s"
                msg = msg % (len(names), ", ".join(names))
                if wx.YES == common_dialogs.message(msg,
                                                    "Delete transform_kernel(s)",
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
                    for kernel in kernels:
                        if self.db.fetch_transform_kernel_referrers(kernel.id):
                            in_use = True
                            break
                            
                    if in_use:
                        msg = "One or more of these design is in use "     \
                              "(referenced by a pulse sequence).\n\n"       \
                              "Pulse designs that are in use may not "     \
                              "be deleted."
                        common_dialogs.message(msg, "Delete Pulse Design(s)")
                    else:
                        # Hasta la vista!
                        ids = [kernel.id for kernel in kernels]
                        self.db.delete_transform_kernels(ids)
                        self.populate_transform_kernels()
            else:
                common_dialogs.message("Please select one or more pulse designs to delete.")


    def on_import(self, event):
        filename = common_dialogs.pickfile("Select Import File")
        
        if filename:
            msg = ""
            try:
                importer = util_import.TransformKernelImporter(filename, self.db)
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

                self.populate_transform_kernels()
                
                msg = """Pulse imported %d of %d transform_kernel(s).
                         Would you like to see the import log?"""            
                msg = msg % (len(importer.imported), importer.found_count)
                if wx.YES == common_dialogs.message(msg, "Import Complete", 
                                                    common_dialogs.Q_YES_NO):
                    wx_util.display_file(importer.log_filename)


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
            kernels = self.get_selected_transform_kernels()
        
            if kernels:
                dialog = dialog_export.DialogExport(self)
                dialog.ShowModal()

                if dialog.export:
                    filename = dialog.filename
                    comment  = dialog.comment
                    compress = dialog.compress
                    
                    wx.BeginBusyCursor()
                    # Remember that these are lightweight transform_kernel preview 
                    # objects and not fully fledged transform_kernel objects. I have to 
                    # load proper copies of the transform_kernels before invoking the 
                    # exporter.
                    kernels = [self.db.fetch_transform_kernel(pp.id) for pp in kernels]
                    try:
                        util_export.export(filename, kernels, self.db, comment, compress)
                    except IOError as xxx_todo_changeme:
                        (error_number, error_string) = xxx_todo_changeme.args
                        msg = """Exporting to "%s" failed. The operating system message is below --\n\n""" % filename
                        msg += error_string
                        common_dialogs.message(msg, "Vespa Export",  common_dialogs.E_OK)
                
                    # Reload list in case any transform_kernels are newly public.
                    self.populate_transform_kernels()

                    wx.EndBusyCursor()
                #else:
                    # The user cancelled the export dialog.
                
                dialog.Destroy()

                self.ListTransformKernels.SetFocus()
            else:
                common_dialogs.message("Please select one or more pulse designs to export.")


    def on_close(self, event):
        self.Close()


    def on_selection_changed(self, event=None):
        # This function is sometimes called internally (rather than by wx)
        # in which case event is None. 
        self.set_button_status()


    def on_transform_kernel_activated(self, event):
        if self.ButtonView.IsEnabled():
            kernels = self.get_selected_transform_kernels()

            if len(kernels):
                self.view_transform_kernel(kernels[0])

    def on_hide_deprecated(self, event):
        self.populate_transform_kernels()
    
    

    ##################    Internal helper functions     ##################
    
    def get_selected_open_intersection(self):
        """
        Returns the intersection between the set of currently selected
        designs and the designs open in tabs. This is useful because some
        operations (e.g. delete) can't be performed when this set is not 
        empty.
        """
        open_transform_kernels = wx.GetApp().vespa.get_open_designs()
        
        selected_transform_kernels = self.get_selected_transform_kernels()
        
        # I can't use set() to construct the intersection because these might
        # be (and probably will be) different objects representing the same
        # pulse design.
        open_ids = [transform_kernel.id for transform_kernel in open_transform_kernels]
        
        return [transform_kernel for transform_kernel in selected_transform_kernels
                              if transform_kernel.id in open_ids]
        
    
    def get_selected_transform_kernels(self):
        """
        Returns a (possibly empty) list of the currently selected
        transform_kernels.
        """
        indices = wx_util.get_selected_item_indices(self.ListTransformKernels)
        
        return [transform_kernel for i, transform_kernel in enumerate(self.transform_kernels)
                                                   if (i in indices)]
        

    def populate_transform_kernels(self):
        """
        Clears & re-populates the list of transform_kernels
        
        """
        # get control settings
        hide_deprecated = self.CheckHideDeprecated.GetValue()
        
        # Build a list of the ids of the currently selected items. I'll
        # use this after re-populating the list.
        selected = [kernel.id for kernel in self.get_selected_transform_kernels()]

        self.ListTransformKernels.DeleteAllItems()

        self.transform_kernels = self.db.fetch_transform_kernels()
        
        if hide_deprecated:
            self.transform_kernels = [item for item in self.transform_kernels if not item.deprecated]
        
        self.all_names  = [kernel.name       for kernel in self.transform_kernels]
        self.all_labels = [kernel.menu_label for kernel in self.transform_kernels]

        # Populate the listbox - guided by hide_deprecated checkbox
        frozen_color = wx.Colour(*common_constants.FROZEN_COLOR)
        for i, transform_kernel in enumerate(self.transform_kernels):
            self.ListTransformKernels.InsertItem(i, transform_kernel.name)
            self.ListTransformKernels.SetItem(i, 1, transform_kernel.menu_label)

            # Mark the public column if necessary
            public = "x" if transform_kernel.is_public else " "
            self.ListTransformKernels.SetItem(i, 2, public)

#             # Display referrer count if non-zero
#             referrers = len(transform_kernel.referrers)
#             referrers = (str(referrers) if referrers else "")
#             self.ListTransformKernels.SetItem(i, 3, referrers)

            if transform_kernel.is_frozen:
                # Frozen!
                self.ListTransformKernels.SetItemBackgroundColour(i, frozen_color)

        self.ListTransformKernels.SetColumnWidth(0, wx.LIST_AUTOSIZE)
        self.ListTransformKernels.SetColumnWidth(1, wx.LIST_AUTOSIZE_USEHEADER)
        self.ListTransformKernels.SetColumnWidth(2, wx.LIST_AUTOSIZE_USEHEADER)

        # Reselect all the items that were selected before, if possible.
        if self.transform_kernels:
            if selected:
                self.select_ids(selected)

            if not self.get_selected_transform_kernels():            
                # Nothing is selected, so I select the first item.
                wx_util.select_list_ctrl_items(self.ListTransformKernels, 0)

        self.ListTransformKernels.SetFocus()
        self.on_selection_changed()

    
    def select_ids(self, ids):
        """Attempts to select the list indices that correspond to the ids
        in the param. Any non-existent ids are ignored.
        """
        indices = [ ]
        
        indices = [i for i, transform_kernel in enumerate(self.transform_kernels)
                                          if transform_kernel.id in ids]
        
        wx_util.select_list_ctrl_items(self.ListTransformKernels, indices)
        
        if indices:
            self.ListTransformKernels.EnsureVisible(indices[0])
            
            
    def set_button_status(self):
        """Enables/disables buttons based on what is selected in the list"""
        transform_kernels = self.get_selected_transform_kernels()
        
        enable = (len(transform_kernels) == 1)
        self.ButtonView.Enable(enable)

        enable = bool(len(transform_kernels))
        self.ButtonClone.Enable(enable)
        self.ButtonDelete.Enable(enable)
        self.ButtonExport.Enable(enable)
        
#         # Can't delete pulse designs that are in use
#         if any([bool(transform_kernel.referrers) for transform_kernel 
#                                               in transform_kernels]):
#             self.ButtonDelete.Disable()


    def view_transform_kernel(self, kernel):
        # Since this is just an transform_kernel preview, I have to fetch the 
        # full transform_kernel in order to display it properly.
#        kernel = self.db.fetch_transform_kernel(kernel.id)

        dialog = dialog_view_transform_kernel.DialogViewTransformKernel(self, kernel)
        
        dialog.ShowModal()
        
        dialog.Destroy()
