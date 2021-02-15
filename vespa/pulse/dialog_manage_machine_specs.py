# Python modules

import copy  

# 3rd party modules
import wx

# Our modules
import vespa.pulse.dialog_machine_specs as dialog_machine_specs
import vespa.pulse.auto_gui.manage_machine_specs as manage_machine_specs

import vespa.common.util.import_ as util_import
import vespa.common.wx_gravy.common_dialogs as common_dialogs
import vespa.common.wx_gravy.util as wx_util



class DialogManageMachineSpecs(manage_machine_specs.MyDialog):
    """
    Displays the dialog for machine_setting managment (view, delete, import
    export, etc.). The parent param is for the parent window and may be
    None. The db param must be a Database instance.
    """
    def __init__(self, parent, db):
        if not parent:
            parent = wx.GetApp().GetTopWindow()

        manage_machine_specs.MyDialog.__init__(self, parent)

        self.db = db
        
        # self.templates is a list of machine_specs template objects;
        # there's a 1:1 correspondence between items in this list and items
        #  in the dialog's combobox.
        self.templates = ( )

        self.ListMachineSpecs.InsertColumn(0, "Name")
        self.ListMachineSpecs.InsertColumn(1, "Is Default", wx.LIST_FORMAT_CENTER)

        self.populate_machine_specs()

        self.Layout()
        self.Fit()
        # Sometimes the default width is pretty darn skinny so we make it
        # a little wider.
        width, height = self.GetSize()
        self.SetMinSize( (450, height) )
        
        self.Center()

        # Set focus on the list.
        self.ListMachineSpecs.SetFocus()

        # Grab key events so I can look for Ctrl/Cmd+A ==> select all
        self.ListMachineSpecs.Bind(wx.EVT_KEY_DOWN, self.on_key_down)

        self.set_button_status()



    ##################    EVENT HANDLERS     ##################
    
    def on_key_down(self, event):
        if wx_util.is_select_all(event):
            wx_util.select_list_ctrl_items(self.ListMachineSpecs)
            
            self.on_selection_changed()
            
        event.Skip()                    


    def on_new(self, event):
        dialog = dialog_machine_specs.DialogMachineSpecs(self, self.db, is_template=True)
        if dialog.ShowModal() == wx.ID_OK:
            self.db.insert_machine_specs(dialog.machine_specs)
            self.populate_machine_specs()
            
            # Deselect all items in the list and select new item.
            wx_util.select_list_ctrl_items(self.ListMachineSpecs, select=False)
            index = self.ListMachineSpecs.FindItem(-1, dialog.machine_specs.name)
            self.ListMachineSpecs.Select(index)

        dialog.Destroy()

        self.ListMachineSpecs.SetFocus()


    def on_edit(self, event):
        # Even if multiple templates are selected, we only edit the first
        templates = self.get_selected_items()
        if len(templates):
            template = self.db.fetch_machine_specs(templates[0].id)
            self.edit_template(template)
            

    def on_view(self, event):
        # Even if multiple templates are selected, we only view the first
        templates = self.get_selected_items()
        if len(templates):
            template = self.db.fetch_machine_specs(templates[0].id)
            self.view_machine_setting(template)


    def on_clone(self, event):
        templates = self.get_selected_items()
        for template in templates:
            # Create the copy
            new_template = copy.copy(template)
            
            # Don't clone the default setting!
            new_template.is_default = False
            
            new_template.name = self.db.find_unique_name(new_template, "clone")

            self.db.insert_machine_specs(new_template)

            self.populate_machine_specs()
            
            # Select the newly-added item
            wx_util.select_list_ctrl_items(self.ListMachineSpecs, select=False)
            index = self.ListMachineSpecs.FindItem(-1, new_template.name)
            self.ListMachineSpecs.Select(index)


    def on_default(self, event):
        templates = self.get_selected_items()
        
        # There should be exactly one (but we won't assume)
        if templates:
            template = templates[0]
        
            self.db.mark_machine_specs_template_as_default(template.id)
            self.populate_machine_specs()
        


    def on_delete(self, event):
        msg = ""
        templates = self.get_selected_items()

        if not templates:
            msg = "Please select one or more machine specs templates to delete."
            
        if not msg:
            is_default = [template.is_default for template in templates]
            
            if any(is_default):
                msg = "You may not delete the default template."
            
        if not msg:
            names = [template.name for template in templates]
            msg = "Are you sure you want to delete the following %d machine setting template(s)?\n\n%s"
            
            msg = msg % (len(names), ", ".join(names))
            if wx.YES == common_dialogs.message(msg,
                                                "Delete Machine Setting Template(s)",
                                                common_dialogs.Q_YES_NO):
                # Hasta la vista!
                ids = [template.id for template in templates]
                self.db.delete_machine_specs(ids)
                self.populate_machine_specs()
        else:
            common_dialogs.message(msg)


    def on_close(self, event):
        self.Close()


    def on_selection_changed(self, event=None):
        # This function is sometimes called internally (rather than by wx)
        # in which case event is None.
        self.set_button_status()


    def on_template_activated(self, event):
        # "activated" means double clicked (or the user hit enter).
        if self.ButtonEdit.IsEnabled():
            self.edit_template(self.templates[event.Index])




    ##################    Internal helper functions     ##################
    ##################      in alphabetical order       ##################
    
    def edit_template(self, template):
        """shows a dialog for editing the machine setting template"""
        # In case the editing dialog decides to fiddle with the template, I
        # give it a copy of mine so it won't mess up my pristine version.
        template = copy.deepcopy(template)
        
        dialog = dialog_machine_specs.DialogMachineSpecs(self, 
                                                               self.db,
                                                               template,
                                                               False, True)        
        if dialog.ShowModal() == wx.ID_OK:
            self.db.update_machine_specs(template)
            self.populate_machine_specs()

        dialog.Destroy()
        self.ListMachineSpecs.SetFocus()
        # Select the changed item.
        wx_util.select_list_ctrl_items(self.ListMachineSpecs, select=False)
        index = self.ListMachineSpecs.FindItem(-1, template.name)
        self.ListMachineSpecs.Select(index)

    
    def get_selected_items(self):
        """
        Returns a (possibly empty) list of the currently selected
        machine_specs.
        """
        selected = [ ]
        indices = wx_util.get_selected_item_indices(self.ListMachineSpecs)
        for i, template in enumerate(self.templates):
            if i in indices:
                selected.append(template)
        
        return selected


    def populate_machine_specs(self):
        """Clears & repopulates the list of machine_specs"""
        self.ListMachineSpecs.DeleteAllItems()

        self.templates = self.db.fetch_machine_specs_templates()

        # Populate the listbox
        for i, template in enumerate(self.templates):
            self.ListMachineSpecs.InsertItem(i, template.name)
            default = "x" if template.is_default else " "
            self.ListMachineSpecs.SetItem(i, 1, default)
            
        self.ListMachineSpecs.SetColumnWidth(0, wx.LIST_AUTOSIZE)
        self.ListMachineSpecs.SetColumnWidth(1, wx.LIST_AUTOSIZE_USEHEADER)

        if self.templates:
            # Select the first item. 
            wx_util.select_list_ctrl_items(self.ListMachineSpecs, 0)

        self.on_selection_changed()

        self.ListMachineSpecs.SetFocus()


    def set_button_status(self):
        """Enables/disables buttons based on what is selected in the list"""
        templates = self.get_selected_items()
        
        enable = (len(templates) == 1)
        self.ButtonView.Enable(enable)
        self.ButtonEdit.Enable(enable)
        self.ButtonDefault.Enable(enable)

        enable = bool(len(templates))
        self.ButtonClone.Enable(enable)
        self.ButtonDelete.Enable(enable)


    def view_machine_setting(self, template):            
        dialog = dialog_machine_specs.DialogMachineSpecs(self, self.db, template, True, True)
        dialog.ShowModal()
        dialog.Destroy()
