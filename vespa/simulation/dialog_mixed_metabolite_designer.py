# Python modules


# 3rd party modules
import wx

# Our modules
import vespa.simulation.auto_gui.mixed_metabolite_designer as mixed_metabolite_designer
import vespa.common.wx_gravy.common_dialogs as common_dialogs
import vespa.common.wx_gravy.util as wx_util
from wx.lib.agw.floatspin import FloatSpin, EVT_FLOATSPIN, FS_LEFT, FS_RIGHT, FS_CENTRE, FS_READONLY
from vespa.common.wx_gravy.widgets.floatspin_multiplier.floatspin_multiplier_base import FloatSpinMultiplier


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


class DialogMixedMetaboliteDesigner(mixed_metabolite_designer.MyDialog):
    """
    Displays the dialog for defining mixtures of metabolite results for use
    in the DialogMixedMetaboliteOutput dialog. This dialog is accessed when
    the user hits the "Add Metabolite Mixture ..." button. The parent param 
    is for the parent window and may not be None. The names variable can be
    a list of default metbolite names that populate the dynamic list when it 
    is initialized.
    """

    def __init__(self, parent, names, abbr):
        
        if not parent:
            parent = wx.GetApp().GetTopWindow()

        mixed_metabolite_designer.MyDialog.__init__(self, parent)

        self.names       = names
        self.abbr        = abbr
        self.unique_name = None
        self.values      = None
        
        self.SetTitle("Mixed Metabolite Designer")
           
        #------------------------------
        # Initialize widget controls
        
        self.initialize_controls()

        self.Layout()

        self.Fit()
        
        # The dialog will fit itself so that everything is visible. We force
        # this to be the min size.
        self.SetMinSize(self.GetSize())
        
        self.Center()



    ##### Event Handlers ######################################################

    def on_select_all(self, event): 
        self.dynamic_output_list.select_all()

    def on_deselect_all(self, event): 
        self.dynamic_output_list.deselect_all()

    def on_add_metabolite(self, event): 
        self.dynamic_output_list.add_row()
        self.Layout()
        self.Fit()

    def on_remove_selected(self, event): 
        self.dynamic_output_list.remove_checked_rows()
        self.Layout()
        self.Fit()

    def on_ok(self,event):
        
        if self.TextUniqueName.GetValue() in self.abbr:
            msg = "You already have a mixture with this name. Mixture names must be unique."
            common_dialogs.message(msg, None, common_dialogs.I_OK)
        else:
            self.unique_name = self.TextUniqueName.GetValue()
            self.values = self.dynamic_output_list.get_values()
    
            self.EndModal(wx.ID_OK)
            

    ##### Internal helper functions  ##########################################

    def initialize_controls(self):
        placeholder = self.LabelGridPlaceholder
        self.GridSizer = placeholder.GetContainingSizer()
        placeholder.Destroy()
        
        # Add headings to the first row of the grid sizer.
        self.GridSizer.Clear()
        self.GridSizer.SetRows(1)
        headings = (None, None, "Scale")
        
        for heading in headings:
            if heading:
                label = wx.StaticText(self, label=heading, style=wx.ALIGN_CENTRE)
                self.GridSizer.Add(label, 0, wx.EXPAND|wx.ALIGN_CENTER_VERTICAL)
            else:
                self.GridSizer.AddSpacer( 1 )

        
        # Add widgets to dialog that wxGlade could not
        self.dynamic_output_list = DynamicDesignList(self,
                                                     self.GridSizer,
                                                     self.names)

        # We add the OK & Cancel buttons dynamically so that they're in the 
        # right order under OS X, GTK, Windows, etc.
        self.ButtonOk, self.ButtonCancel = \
                    wx_util.add_ok_cancel(self, self.LabelOkCancelPlaceholder,
                                          self.on_ok)

       
###############################################################################
##### Dynamic Design List Class ###############################################

class DynamicDesignList(object):

    def __init__(self, parent, grid_sizer, names):
        
        self.parent = parent
        self.names = names
        # We follow the wx CamelCaps naming convention for this wx object.
        self.GridSizer = grid_sizer
        
        self.list_lines = []

        # default to having two copies of first metabolite in list
        self.add_row()
        self.add_row()

    
    def add_row(self, selected_metabolite=""):
        # Build the individual controls
        checkbox = wx.CheckBox(self.parent)
        metabolites = wx.Choice(self.parent, choices=self.names)
        scale = FloatSpinMultiplier(self.parent, 
                                     style=wx.SP_ARROW_KEYS|wx.SP_WRAP|wx.TE_PROCESS_ENTER, 
                                     agwStyle=FS_LEFT)

        # Keep a reference to the controls so I can use them later
        self.list_lines.append({"checkbox" : checkbox, 
                                "metabolites" : metabolites, 
                                "scale" : scale
                                })

        # Configure the controls I just created
        i = metabolites.FindString(selected_metabolite)
        if i == wx.NOT_FOUND:
            i = 0
        metabolites.SetSelection(i) 
        
        scale.SetValue(1.0)
        scale.multiplier = 1.25
        scale.SetDigits(5)
        scale.SetIncrement(1.0)
        scale.SetRange(0.00001,10000.0)
        scale.SetMinSize((90, -1))
        
        # add controls to the grid sizer
        self.GridSizer.SetRows(self.GridSizer.GetRows() + 1)

        self.GridSizer.Add(checkbox, 0, wx.ALIGN_CENTER_VERTICAL)
        self.GridSizer.Add(metabolites, 0, wx.EXPAND)
        self.GridSizer.Add(scale, 0, wx.EXPAND)
        

        
    def remove_checked_rows(self):
        # gather indices of all checked boxes
        checklist = []
        
        for i, line in enumerate(self.list_lines):
            if line["checkbox"].GetValue():
                checklist.append(i)
        
        # remove in reverse order so we don't invalidate later
        # indices by removing the ones preceding them in the list
        checklist.reverse()
        
        for i in checklist:
            # Each line is a dict of controls
            for control in list(self.list_lines[i].values()):
                control.Destroy()
                
            del self.list_lines[i]
            
        # Reduce the # of rows in the grid sizer
        rows = self.GridSizer.GetRows()
        self.GridSizer.SetRows(rows - len(checklist))
        self.GridSizer.Layout()
        

    def get_values(self):
        return [self.get_line_values(line) for line in self.list_lines]


    def get_line_values(self, line):
        return (line["metabolites"].GetStringSelection(),
                line["scale"].GetValue()
               )
                

    def set_check(self, i, value=False):
        if self.list_lines == []: return
        nlines = len(self.list_lines)        
        if i < 0 or i >= nlines: return
        self.list_lines[i][1].SetValue(value==True)

        
    def select_all(self):
        for line in self.list_lines:
            line["checkbox"].SetValue(True)

        
    def deselect_all(self):
        for line in self.list_lines:
            line["checkbox"].SetValue(False)
   
