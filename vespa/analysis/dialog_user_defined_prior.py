# Python modules
import copy
import collections

# 3rd party modules
import wx

# Our modules
import vespa.analysis.mrs_dataset as mrs_dataset
import vespa.common.util.misc as util_misc
import vespa.common.mrs_prior as mrs_prior_module
import vespa.common.mrs_prior_metabolite as mrs_prior_metabolite
import vespa.common.wx_gravy.util as wx_util
import vespa.analysis.auto_gui.dialog_user_defined_prior as dialog_user_defined_prior

from vespa.analysis.dialog_user_defined_metabolite import DialogUserDefinedMetabolite




class DialogUserDefinedPrior(dialog_user_defined_prior.MyDialog):
    """
    Used by GISO model to allow user to define a multi-line, multi-metabolite Prior
    
    """
    def __init__(self, parent, dataset, mrs_prior=None, 
                                        constraints_eq=[], 
                                        constraints_ineq=[]):
        
        dialog_user_defined_prior.MyDialog.__init__(self, parent)

        if mrs_prior is None:
            mrs_prior = mrs_prior_module.Prior()
            mrs_prior.source    = 'DialogUserDefinedPrior'
            mrs_prior.source_id = util_misc.uuid()
            mrs_prior.comment   = 'none'

        self.dataset          = dataset
        self.mrs_prior        = mrs_prior
        self.constraints_eq   = constraints_eq
        self.constraints_ineq = constraints_ineq

        # copy originals in case user hits cancel
        self.original_mrs_prior        = copy.deepcopy(mrs_prior)
        self.original_constraints_eq   = list(constraints_eq)
        self.original_constraints_ineq = list(constraints_ineq)

        # Open & Cancel added dynamically so in order under OS X, GTK, Windows
        r = wx_util.add_ok_cancel(self, self.LabelOKCancelPlaceholder, self.on_ok, self.on_cancel)
        self.ButtonOpen, self.ButtonCancel = r

        self.SetSize( (550, 750) )
        self.Layout()
        self.Center()

        self.initialize_controls()
        self.populate_controls()
    

    #################     Internal Helpers    ##################
    
    def initialize_controls(self):
        """ 
        Set up sizes and constraints for widgets. Does not set widget values
        - See populate_controls() to set values from data object.
        
        """

        #------------------------------------------------------------
        # set up dynamic list - metabolite groups
        
        self.GridSizerGroups = self.LabelPlaceholderGroups.GetContainingSizer()
        parent = self.LabelPlaceholderGroups.GetParent()
        self.LabelPlaceholderGroups.Destroy()
       
        # Add headings to the first row of the grid sizer.
        self.GridSizerGroups.Clear()
        self.GridSizerGroups.SetRows(1)
        headings = (None, ("Unique Group Name",wx.ALIGN_LEFT), ("Number of Lines",wx.ALIGN_CENTRE))
        for heading in headings:
            if heading:
                label = wx.StaticText(parent, label=heading[0], style=heading[1])
                self.GridSizerGroups.Add(label, 0, wx.EXPAND|wx.ALIGN_CENTER_VERTICAL)
            else:
                self.GridSizerGroups.AddSpacer( 1 )

        self.dynamic_list1 = DynamicList1( self.PanelGroups,
                                           self.PanelUserDefinedPrior,
                                           self.GridSizerGroups,
                                           self.mrs_prior,
                                           self.on_dynamic_list1)

        #------------------------------------------------------------
        # set up dynamic list - equality constraints
         
        # The list grid sizer is marked so we can find it at run-time
        self.GridSizerEquality = self.LabelPlaceholderEquality.GetContainingSizer()
        parent = self.LabelPlaceholderEquality.GetParent()
        self.LabelPlaceholderEquality.Destroy()
        
        # Add headings to the first row of the grid sizer.
        self.GridSizerEquality.Clear()
        self.GridSizerEquality.SetRows(1)
        headings = (None, ("     Constraint String",wx.ALIGN_LEFT))
        for heading in headings:
            if heading:
                label = wx.StaticText(parent, label=heading[0], style=heading[1])
                self.GridSizerEquality.Add(label, 0, wx.EXPAND|wx.ALIGN_CENTER_VERTICAL)
            else:
                self.GridSizerEquality.AddSpacer( 1 )

 
        self.dynamic_list2 = DynamicList2( self.PanelEquality,
                                           self.PanelUserDefinedPrior,
                                           self.GridSizerEquality,
                                           self.constraints_eq,
                                           self.on_dynamic_list2)

        #------------------------------------------------------------
        # set up dynamic list - inequality constraints
         
        # The list grid sizer is marked so we can find it at run-time
        self.GridSizerInequality = self.LabelPlaceholderInequality.GetContainingSizer()
        parent = self.LabelPlaceholderInequality.GetParent()
        self.LabelPlaceholderInequality.Destroy()
        
        # Add headings to the first row of the grid sizer.
        self.GridSizerInequality.Clear()
        self.GridSizerInequality.SetRows(1)
        headings = (None, ("     Constraint String",wx.ALIGN_LEFT))
        for heading in headings:
            if heading:
                label = wx.StaticText(parent, label=heading[0], style=heading[1])
                self.GridSizerInequality.Add(label, 0, wx.EXPAND|wx.ALIGN_CENTER_VERTICAL)
            else:
                self.GridSizerInequality.AddSpacer( 1 )
 
        self.dynamic_list3 = DynamicList3( self.PanelInequality,
                                           self.PanelUserDefinedPrior,
                                           self.GridSizerInequality,
                                           self.constraints_ineq,
                                           self.on_dynamic_list3)


    def populate_controls(self):
        """ 
        Populate widgets with values from the data object. Method trusts that
        the data object does not violate ranges/constraints. It does enable
        and disable other widgets depending on settings.

        """
        mrs_prior  = self.mrs_prior
        equality   = self.constraints_eq
        inequality = self.constraints_ineq

        self.TextCtrlComment.SetValue(mrs_prior.comment)

        self.dynamic_list1.set_new_values(mrs_prior,  add_check=False)
        self.dynamic_list2.set_new_values(equality,   add_check=False)
        self.dynamic_list3.set_new_values(inequality, add_check=False)
        


##############################    Event Handlers

    def on_ok(self, event):
        """ 
        Collate MrsPrior and Equality and Inequality constraints, and quit
         - note. mrs_prior is updated as we go
         - constraint strings could be changed any time, thus we collect them here
    
        """
        self.mrs_prior.comment = self.TextCtrlComment.GetValue()
        lines2 = [item['constraint'] for item in self.dynamic_list2.lines]
        lines3 = [item['constraint'] for item in self.dynamic_list3.lines]
        self.constraints_eq   = lines2
        self.constraints_ineq = lines3
        self.EndModal(wx.ID_OK)

    def on_cancel(self, event):
        """ 
        Restore original settings before exiting.
        
        """
        self.mrs_prior        = self.original_mrs_prior
        self.constraints_eq   = self.original_constraints_eq
        self.constraints_ineq = self.original_constraints_ineq
        self.EndModal(wx.ID_CANCEL)
    
    def on_add_group(self, event):
        """
        Create a PriorMetabolite in a separate dialog
        Stuff it into the mrs_prior object back in this dialog
        Add line to the Metabolite Group dynamic list
        
        """
        dialog = DialogUserDefinedMetabolite(self, self.dataset, unique_names=self.mrs_prior.names)
        code = dialog.ShowModal() 
        
        if code == wx.ID_OK:
            val = dialog.metabolite
            self.mrs_prior.metabolites[val.name] = val  # added at end of ordered dict
            self.dynamic_list1.add_row([False,val.name,len(val.areas)], update=True)
            
            # add/remove a group invalidates contraint equations, delete them all
            self.dynamic_list2.select_all()
            self.dynamic_list3.select_all()
            self.dynamic_list2.remove_checked_rows()
            self.dynamic_list3.remove_checked_rows()
            
        self.Layout()

    def on_edit_group(self, event):
        """ 
        Pick first selected line in Metabolite Group list
        Display selected metabolite in UserDefinedMetabolite dialog
         - edit values and hit return/cancel
        Return to this dialog and:
          - remove old metabolite entry from mrs_prior dict
          - insert edited metabolite into mrs_prior dict, may have new name
          - update line for selected metabolite group to new values
          
        The MAIN ISSUE ... all the constraints entered in the lower part of the
        GUI depend on the 'position' of each metabolite in the mrs_prior 
        ordered dict to determine the x[0] or x[10] parameter values that make
        up the constraint. If I edit a metabolite, and especially if I change 
        the name (which is the key in the ordered dict) I have to be sure that
        the 'newly edited' metabolite holds the same position.  I do this by
        creating a whole new ordered dict and copy all the metabs into it except
        the newly edited one which is given back its original position in the 
        new dict.
        
        """
        key = self.dynamic_list1.get_first_checked_group_name()
        if key:
            metabolite   = self.mrs_prior.metabolites[key]
            unique_names = self.mrs_prior.names
            unique_names.remove(key)                # want to be able to save to same name
    
            dialog = DialogUserDefinedMetabolite(self, self.dataset, 
                                                       metabolite=metabolite, 
                                                       unique_names=unique_names)
            code = dialog.ShowModal() 
            
            if code == wx.ID_OK:
                val = dialog.metabolite
                
                new_dict = collections.OrderedDict()
                old_dict = self.mrs_prior.metabolites
                for new_key in list(self.mrs_prior.metabolites.keys()):
                    if new_key == key:
                        new_dict[val.name] = val    # in case name changed
                    else:
                        new_dict[new_key]  = old_dict[new_key]
                self.mrs_prior.metabolites = new_dict
                self.dynamic_list1.update_row(key, [True,val.name,len(val.areas)], update=True)
                self.Layout()

    def on_delete_group(self, event): 
        """
        See if any metabolite group lines checked
        Remove all checked rows, collect metabolite names from each
        Remove metabolite objects from mrs_prior using their name/keys
        
        """
        keys = self.dynamic_list1.remove_checked_rows()
        for key in keys:
            self.mrs_prior.metabolites.pop(key, None)

        # add/remove a group invalidates contraint equations, delete them all
        self.dynamic_list2.select_all()
        self.dynamic_list3.select_all()
        self.dynamic_list2.remove_checked_rows()
        self.dynamic_list3.remove_checked_rows()

        self.Layout()

    def on_dynamic_list1(self, event=None):
        # This is a fudged event called from the actual event that occurs 
        # inside the dynamic list class but its also invoked programmatically
        # (i.e. by our code, not by wx). In the latter case, event is None.
        pass
    
    def on_add_equality(self, event):
        """
        We do not update strings until the dialog quits
        Also, this event triggers on_dynamic_list2() 
        
        """
        self.dynamic_list2.add_row([False,''])
        self.Layout()

    def on_delete_equality(self, event): 
        """
        We do not update strings until the dialog quits
        Also, this event triggers on_dynamic_list2() 
        
        """
        self.dynamic_list2.remove_checked_rows()
        self.Layout()

    def on_dynamic_list2(self, event=None):
        # This is a fudged event called from the actual event that occurs 
        # inside the dynamic list class but its also invoked programmatically
        # (i.e. by our code, not by wx). In the latter case, event is None.
        pass    

    def on_add_inequality(self, event):
        """
        We do not update strings until the dialog quits
        Also, this event triggers on_dynamic_list3() 
        
        """
        self.dynamic_list3.add_row([False,''])
        self.Layout()

    def on_delete_inequality(self, event): 
        """
        We do not update strings until the dialog quits
        Also, this event triggers on_dynamic_list3() so no more needed here
        
        """
        self.dynamic_list3.remove_checked_rows()
        self.Layout()

    def on_dynamic_list3(self, event=None):
        # This is a fudged event called from the actual event that occurs 
        # inside the dynamic list class but its also invoked programmatically
        # (i.e. by our code, not by wx). In the latter case, event is None.
        pass    




#------------------------------------------------------------------------------

class DynamicList1(object):

    def __init__(self, PanelParent, PanelTop, GridSizer, mrs_prior, external_event_handler):
        
        self.mrs_prior   = mrs_prior
        self.PanelParent = PanelParent
        self.PanelTop    = PanelTop
        self.GridSizer   = GridSizer
        
        self.external_event_handler = external_event_handler
        self.list_lines = []

    @property
    def lines(self):
        return [self.get_line_values(line) for line in self.list_lines]
        
    def set_new_values(self, previous=None, add_check=None):
        if not previous:
            previous = self.mrs_prior

        self.select_all()
        self.remove_checked_rows()

        for key in list(previous.metabolites.keys()):
            item = previous.metabolites[key]
            row = []
            if add_check is not None:
                row.append(add_check==True)
                    
            row.append(item.name)
            row.append(len(item.areas))
            self.add_row(row)
    
    def add_row(self, row_vals, update=False):
        '''
        Adds a row to the end of the list. 
        
        '''
        self.GridSizer.SetRows(self.GridSizer.GetRows() + 1)

        # create widgets to go into the line
        list_line = {}

        checkbox     = wx.CheckBox(self.PanelParent)
        value_name   = wx.StaticText(self.PanelParent)
        value_lines  = wx.StaticText(self.PanelParent)
        
        # keep a copy of panel and widgets to access later
        line = { "check"        : checkbox, 
                 "value_name"   : value_name, 
                 "value_lines"  : value_lines
               }

        # Add the controls to the grid sizer
        self.GridSizer.Add(line["check"], 0, wx.ALIGN_CENTER_VERTICAL)
        self.GridSizer.Add(line['value_name'], 0, wx.EXPAND)
        self.GridSizer.Add(line['value_lines'], 0, wx.ALIGN_CENTER, wx.EXPAND)

        # Configure the controls I just created
        checkbox.SetValue(       row_vals[0])
        value_name.SetLabel( str(row_vals[1]))
        value_lines.SetLabel(str(row_vals[2]))

        self.list_lines.append(line)

        # only need to update if a metabolite is added/removed from basis
        self.PanelParent.Bind(wx.EVT_CHECKBOX, self.event_handler, checkbox)

        self.PanelTop.Layout()  

        if update:
            self.event_handler()

    def update_row(self, key, row_vals, update=False):
        '''
        Adds a row to the end of the list. 
        
        '''
        indx = []
        for i, line in enumerate(self.list_lines):
            if line["value_name"].GetLabel() == key:
                indx.append(i)
        if indx:
            self.list_lines[indx[0]]["check"].SetValue(row_vals[0])
            self.list_lines[indx[0]]["value_name"].SetLabel(str(row_vals[1]))
            self.list_lines[indx[0]]["value_lines"].SetLabel(str(row_vals[2]))        

        self.PanelTop.Layout()  

        if update:
            self.event_handler()
                    
    def get_first_checked_group_name(self):
        for i, line in enumerate(self.list_lines):
            if line["check"].GetValue():
                vals = self.get_line_values(line)
                return vals['value_name']
        return ''
        
    def remove_checked_rows(self):
        # gather indices of all checked boxes
        checklist = []
        
        for i, line in enumerate(self.list_lines):
            if line["check"].GetValue():
                checklist.append(i)
        if checklist:
            keys = self._remove_lines_from_checklist(checklist)
            return keys

    def remove_by_name(self, val):
        # gather indices of all checked boxes
        checklist = []
        
        for i, line in enumerate(self.list_lines):
            if line["value_name"].GetLabel() == val:
                checklist.append(i)
        if checklist:
            keys = self._remove_lines_from_checklist(checklist)
    
    def _remove_lines_from_checklist(self, checklist):    
        # get list of names being removed
        keys = []
        for i in checklist:
            keys.append(self.list_lines[i]['value_name'].GetLabel())
        
        # remove in reverse order so we don't invalidate later
        # indices by removing the ones preceding them in the list
        checklist.reverse()
        
        for i in checklist:
            # Each line is a dict of controls + mixture info
            for item in list(self.list_lines[i].values()):
                if hasattr(item, "Destroy"):
                    # It's a wx control
                    item.Destroy()
                
            del self.list_lines[i]
            
        # Reduce the # of rows in the grid sizer
        rows = self.GridSizer.GetRows()
        self.GridSizer.SetRows(rows - len(checklist))
        self.GridSizer.Layout()
        self.PanelTop.Layout()        
        self.event_handler()
        
        return keys

    def select_all(self):
        for line in self.list_lines:
            line["check"].SetValue(True)

    def deselect_all(self):
        for line in self.list_lines:
            line["check"].SetValue(False)
            
    def event_handler(self, event=None):
        self.external_event_handler(event)

    def get_line_values(self, line):
        # Returns a dict containing  values of the controls in the line.
        return { "check"        : line["check"].GetValue(),
                 "value_name"   : line["value_name"].GetLabel(),
                 "value_lines"  : line["value_lines"].GetLabel()
               }


#------------------------------------------------------------------------------

class DynamicList2(object):

    def __init__(self, PanelParent, PanelTop, GridSizer, oplist, external_event_handler):
        
        self.oplist      = oplist
        self.PanelParent = PanelParent
        self.PanelTop    = PanelTop
        self.GridSizer   = GridSizer
        
        self.external_event_handler = external_event_handler
        self.list_lines = []

    @property
    def lines(self):
        return [self.get_line_values(line) for line in self.list_lines]

    def set_new_values(self, previous=None, add_check=None):
        if not previous:
            previous = self.oplist

        self.select_all()
        self.remove_checked_rows()

        for item in previous:
            row = []
            if add_check is not None:
                if add_check is not True:
                    row.append(False)
                else:
                    row.append(True)
                    
            row.append(item)
            self.add_row(row)

    
    def add_row(self, row_vals, update=False):
        '''
        Adds a row to the end of the list. 
        
        '''
        self.GridSizer.SetRows(self.GridSizer.GetRows() + 1)

        # create widgets to go into the line
        list_line = {}

        checkbox   = wx.CheckBox(self.PanelParent)
        constraint = wx.TextCtrl(self.PanelParent, wx.ID_ANY, "")
        
        # keep a copy of panel and widgets to access later
        line = { "check"        : checkbox, 
                 "constraint"   : constraint
               }

        # Add the controls to the grid sizer
        self.GridSizer.Add(line["check"],      0, wx.ALIGN_CENTER_VERTICAL)
        self.GridSizer.Add(line['constraint'], 0, wx.EXPAND, wx.ALIGN_CENTER_VERTICAL)

        # Configure the controls I just created
        checkbox.SetValue(       row_vals[0])
        constraint.SetValue( str(row_vals[1]))

        self.list_lines.append(line)

        # only need to update if a metabolite is added/removed from basis
        self.PanelParent.Bind(wx.EVT_CHECKBOX, self.event_handler, checkbox)

        self.PanelTop.Layout()  

        if update:
            self.event_handler()
        
        
    def remove_checked_rows(self):
        # gather indices of all checked boxes
        checklist = []
        
        for i, line in enumerate(self.list_lines):
            if line["check"].GetValue():
                checklist.append(i)
        
        # remove in reverse order so we don't invalidate later
        # indices by removing the ones preceding them in the list
        checklist.reverse()
        
        for i in checklist:
            # Each line is a dict of controls + mixture info
            for item in list(self.list_lines[i].values()):
                if hasattr(item, "Destroy"):
                    # It's a wx control
                    item.Destroy()
                
            del self.list_lines[i]
            
        # Reduce the # of rows in the grid sizer
        rows = self.GridSizer.GetRows()
        self.GridSizer.SetRows(rows - len(checklist))
        self.GridSizer.Layout()
        self.PanelTop.Layout()        
        self.event_handler()


    def select_all(self):
        for line in self.list_lines:
            line["check"].SetValue(True)


    def deselect_all(self):
        for line in self.list_lines:
            line["check"].SetValue(False)
            
            
    def event_handler(self, event=None):
        self.external_event_handler(event)


    def get_line_values(self, line):
        # Returns a dict containing  values of the controls in the line.
        return { "check"        : line["check"].GetValue(),
                 "constraint"   : line["constraint"].GetValue()
               }


        
#------------------------------------------------------------------------------

class DynamicList3(object):

    def __init__(self, PanelParent, PanelTop, GridSizer, oplist, external_event_handler):
        
        self.oplist      = oplist
        self.PanelParent = PanelParent
        self.PanelTop    = PanelTop
        self.GridSizer   = GridSizer
        
        self.external_event_handler = external_event_handler
        self.list_lines = []


    @property
    def lines(self):
        return [self.get_line_values(line) for line in self.list_lines]

        
    def set_new_values(self, previous=None, add_check=None):
        if not previous:
            previous = self.oplist

        self.select_all()
        self.remove_checked_rows()

        for item in previous:
            row = []
            if add_check is not None:
                if add_check is not True:
                    row.append(False)
                else:
                    row.append(True)
                    
            row.append(item)
            self.add_row(row)

    
    def add_row(self, row_vals, update=False):
        '''
        Adds a row to the end of the list. 
        
        '''
        self.GridSizer.SetRows(self.GridSizer.GetRows() + 1)

        # create widgets to go into the line
        list_line = {}

        checkbox   = wx.CheckBox(self.PanelParent)
        constraint = wx.TextCtrl(self.PanelParent, wx.ID_ANY, "")
        
        # keep a copy of panel and widgets to access later
        line = { "check"        : checkbox, 
                 "constraint"   : constraint
               }

        # Add the controls to the grid sizer
        self.GridSizer.Add(line["check"],      0, wx.ALIGN_CENTER_VERTICAL)
        self.GridSizer.Add(line['constraint'], 0, wx.EXPAND, wx.ALIGN_CENTER_VERTICAL)

        # Configure the controls I just created
        checkbox.SetValue(       row_vals[0])
        constraint.SetValue( str(row_vals[1]))

        self.list_lines.append(line)

        # only need to update if a metabolite is added/removed from basis
        self.PanelParent.Bind(wx.EVT_CHECKBOX, self.event_handler, checkbox)

        self.PanelTop.Layout()  

        if update:
            self.event_handler()
        
        
    def remove_checked_rows(self):
        # gather indices of all checked boxes
        checklist = []
        
        for i, line in enumerate(self.list_lines):
            if line["check"].GetValue():
                checklist.append(i)
        
        # remove in reverse order so we don't invalidate later
        # indices by removing the ones preceding them in the list
        checklist.reverse()
        
        for i in checklist:
            # Each line is a dict of controls + mixture info
            for item in list(self.list_lines[i].values()):
                if hasattr(item, "Destroy"):
                    # It's a wx control
                    item.Destroy()
                
            del self.list_lines[i]
            
        # Reduce the # of rows in the grid sizer
        rows = self.GridSizer.GetRows()
        self.GridSizer.SetRows(rows - len(checklist))
        self.GridSizer.Layout()
        self.PanelTop.Layout()        
        self.event_handler()


    def select_all(self):
        for line in self.list_lines:
            line["check"].SetValue(True)


    def deselect_all(self):
        for line in self.list_lines:
            line["check"].SetValue(False)
            
            
    def event_handler(self, event=None):
        self.external_event_handler(event)


    def get_line_values(self, line):
        # Returns a dict containing  values of the controls in the line.
        return { "check"        : line["check"].GetValue(),
                 "constraint"   : line["constraint"].GetValue()
               }




##  Local test functionality  ###############################################

class MyForm(wx.Frame):
 
    def __init__(self, dataset, mrs_prior=None,constraints_eq=[],constraints_ineq=[]):

        self.dataset          = dataset
        self.mrs_prior        = mrs_prior
        self.constraints_eq   = constraints_eq
        self.constraints_ineq = constraints_ineq

        wx.Frame.__init__(self, None, wx.ID_ANY, "Testing the Dialog")
        panel = wx.Panel(self, wx.ID_ANY)
        buttn = wx.Button(panel, label="Show Dialog")
        buttn.Bind(wx.EVT_BUTTON, self.on_dialog)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(buttn, 0, wx.ALL|wx.CENTER, 5)
        panel.SetSizer(sizer)

        self.statusbar = self.CreateStatusBar(4, 0)
        self.statusbar.SetStatusText("Ready")


    def on_dialog(self, event):

        dialog = DialogUserDefinedPrior(self, self.dataset,
                                        mrs_prior=self.mrs_prior,
                                        constraints_eq=self.constraints_eq,
                                        constraints_ineq=self.constraints_ineq)
        code = dialog.ShowModal()
        if code == wx.ID_OK:
            print("UserDefinedPrior User hit OK")
        elif code == wx.ID_CANCEL:
            print("UserDefinedPrior User Cancelled")

        self.mrs_prior        = dialog.mrs_prior
        self.constraints_eq   = dialog.constraints_eq
        self.constraints_ineq = dialog.constraints_ineq
        
        print(self.mrs_prior)
        print("---------------------------")
        print("Constraints Equality")
        print("---------------------------")
        for i,item in enumerate(self.constraints_eq):
            print(str(i)+' '+item)
        print("---------------------------")
        print("Constraints Inquality")
        print("---------------------------")
        for i,item in enumerate(self.constraints_ineq):
            print(str(i)+' '+item)

        dialog.Destroy()        
 

 
#----------------------------------------------------------------------
# Run the program
if __name__ == "__main__":
    
    dataset = mrs_dataset.Dataset()

    met1 = mrs_prior_metabolite.PriorMetabolite()
    met1.dims[0] = 'met1'
    met1.spins  = 2
    met1.ppms   = [4.7,]
    met1.areas  = [1.0,]
    met1.phases = [0.0,]
    
    met2 = mrs_prior_metabolite.PriorMetabolite()
    met2.dims[0] = 'met2'
    met2.spins  = 2
    met2.ppms   = [4.7,2.0,]
    met2.areas  = [3.0,1.0,]
    met2.phases = [0.0,0.0,]

    met3 = mrs_prior_metabolite.PriorMetabolite()
    met3.dims[0] = 'met3'
    met3.spins  = 2
    met3.ppms   = [4.7,2.2,2.5]
    met3.areas  = [3.0,1.0,0.5]
    met3.phases = [0.0,0.0,0.0]
    
    mrs_prior = mrs_prior_module.Prior()
    mrs_prior.source    = 'Run time test object'
    mrs_prior.source_id = util_misc.uuid()
    mrs_prior.comment   = 'none as of yet'
    mrs_prior.seqte     = 1.0
    mrs_prior.nucleus   = '1H'
    for item in [met1,met2,met3]:
        key = item.dims[0]
        mrs_prior.metabolites[key] = item
    
    constraints_eq   = ['1111111','222222222','33333333']
    constraints_ineq = ['444444444','555555555','666666']
    
    app = wx.App(False)
    frame = MyForm(dataset, mrs_prior=mrs_prior,
                            constraints_eq=constraints_eq,
                            constraints_ineq=constraints_ineq)
    frame.Show()
    app.MainLoop()
