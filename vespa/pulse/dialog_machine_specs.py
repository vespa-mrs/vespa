# Python modules

import copy

# 3rd party modules
import wx

# Our modules
import vespa.pulse.auto_gui.machine_specs as machine_specs_module
import vespa.common.util.misc as util_misc
import vespa.common.rfp_machine_specs as rfp_machine_specs
import vespa.common.constants as constants
import vespa.common.wx_gravy.common_dialogs as common_dialogs
import vespa.common.wx_gravy.util as wx_util



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


class DialogMachineSpecs(machine_specs_module.MyDialog):

    def __init__(self, parent, db, machine_specs=None, read_only=False,
                                                          is_template=False):
        
        # This dialog can be invoked in four different states. They are --
        # 1. Edit machine specs
        # 2. New machine specs template
        # 3. Edit machine specs template
        # 4. View machine specs template
        
        # Notably absent from the list are "New machine specs" and "View
        # machine specs". This dialog could be used for those purposes,
        # but RFPulse doesn't offer a path to get to them.
        
        # In state 1 (when is_template is False), the "Copy from template"
        # label & combobox are visible. They're hidden in all other states.
        # In states 2-4 (when is_template is True), the name label & textbox
        # are visible, otherwise they're hidden.
        if not parent:
            parent = wx.GetApp().GetTopWindow()

        machine_specs_module.MyDialog.__init__(self, parent)
        
        self.db                 = db
        self.machine_specs   = machine_specs
        self.read_only = read_only
        self.is_template = is_template
        
        # Fetch the list o' templates 
        self.templates = db.fetch_machine_specs_templates()
        
        # Find the default template. There's guaranteed to be one.
        self.default_template = None
        for template in self.templates:
            if template.is_default:
                self.default_template = template
                break
        
        
        if machine_specs:
            # We're viewing or editing an existing object
            title = "View" if read_only else "Edit"
        else:
            # We're creating a new machine specs template
            self.read_only = False
            # Create a new machine_specs object based on the default
            self.machine_specs = copy.copy(self.default_template)
            self.machine_specs.is_default = False
            self.machine_specs.id = 0
            self.machine_specs.name = ""
            
            title = "New"
            
        title += " Machine Specs"
        if is_template:
            title += " Template"

        self.SetTitle(title)
        
        self.initialize_controls()

        self.Fit()
        self.Center()



    ##### Event Handlers ######################################################
    
    def on_template_selected(self, event):
        # The user selected a template. I populate the dialog with values from
        # that template.
        # This can only happen in case 1 (see __init__() comment)
        template = self.templates[event.GetSelection()]
        self.machine_specs = rfp_machine_specs.specs_from_template(template)
        
        self.populate_controls()

    
    def on_ok(self, event):
        if self.validate_machine_specs():
            self.EndModal(wx.ID_OK)        


    ##### Internal helper functions  ##########################################

    def initialize_controls(self):
        # sets up widgets in dialog for either a new or existing 
        # rfp_machine_specs object.
        #
        # finally if the dialog is read only, the widgets are disabled
        
        self.LabelName.Show(self.is_template)
        self.TextName.Show(self.is_template)
        
        if self.is_template:
            # cases 2-4, hide the panel that contains these controls
            self.PanelCopyFrom.Hide()
        else:
            # Editing machine specs -- case 1 (see __init__ comment)
            # Populate templates combobox
            names = [template.name for template in self.templates]
            
            self.ComboTemplates.SetItems(names)

            # Select the default (there will always be one)
            i = 0
            while not self.templates[i].is_default:
                i += 1
                
            self.ComboTemplates.SetSelection(i)
            
        # Disable controls as appropriate.
        if self.read_only:
            for control in (self.TextName, 
                            self.ComboMachineType, 
                            self.TextFieldStrength,
                            self.TextMaxB1Field, 
                            self.TextZeroPadding, 
                            self.TextMinDwellTime, 
                            self.TextDwellTimeIncrement, 
                            self.TextGradientRasterTime, 
                            self.TextGradientSlewRate,
                            self.TextGradientMaximum, 
                            self.PanelCopyFrom, ):
                control.Disable()

        # We add the OK & Cancel buttons dynamically so that they're in the 
        # right order under OS X, GTK, Windows, etc.
        self.ButtonOk, self.ButtonCancel = \
                    wx_util.add_ok_cancel(self, self.LabelOkCancelPlaceholder,
                                          self.on_ok)

        if self.read_only:
            # I hide OK rather than destroying it; other code references
            # self.ButtonOk and it's easier to assume it will always exist.
            self.ButtonOk.Hide()
            self.ButtonCancel.SetLabel("Done")

        self.populate_controls()


    def populate_controls(self):
        # populates controls with the data from self.machine_specs.
        if self.is_template:
            self.TextName.ChangeValue(self.machine_specs.name)
        self.TextFieldStrength.ChangeValue(str(self.machine_specs.field_strength))
        self.TextMaxB1Field.ChangeValue(str(self.machine_specs.max_b1_field))
        self.TextZeroPadding.ChangeValue(str(self.machine_specs.zero_padding))
        self.TextMinDwellTime.ChangeValue(str(self.machine_specs.min_dwell_time))
        self.TextDwellTimeIncrement.ChangeValue(str(self.machine_specs.dwell_time_increment))
        self.TextGradientRasterTime.ChangeValue(str(self.machine_specs.gradient_raster_time))
        self.TextGradientSlewRate.ChangeValue(str(self.machine_specs.gradient_slew_rate))
        self.TextGradientMaximum.ChangeValue(str(self.machine_specs.gradient_maximum))
        
        # Populate the machine types combobox
        names = constants.MachineType.get_all_display()
        self.ComboMachineType.Clear()        
        self.ComboMachineType.SetItems(names)
                
        machine_type = self.machine_specs.machine_type
        if machine_type in constants.MachineType.ALL:
            machine_type = machine_type['display']
            self.ComboMachineType.SetStringSelection(machine_type)
        else:
            self.ComboMachineType.SetValue(machine_type)
        
        
    def validate_machine_specs(self):
        # Validate controls one by one
        # A non-empty string message indicates failure
        
        msg = ""

        if self.is_template:
            # name
            name = self.TextName.GetValue().strip()        
            if not name:
                msg = "Please enter a name for this template."
                
        # machine type
        if not msg:
            machine_type = self.ComboMachineType.GetValue().strip()
            if not machine_type:
                msg = "Please enter a machine type."
                
            # Translate this to a constant if it's a known type
            if constants.MachineType.get_type_for_value(machine_type, "display"):
                machine_type = constants.MachineType.get_type_for_value(machine_type, "display")

        # field_strength    
        if not msg:
            field_strength = self.TextFieldStrength.GetValue().strip()
            if not util_misc.is_floatable(field_strength):
                msg = """I don't understand the field strength - "%s".""" % field_strength
    
        # max_b1_field    
        if not msg:
            max_b1_field = self.TextMaxB1Field.GetValue().strip()
            if not util_misc.is_floatable(max_b1_field):
                msg = """I don't understand the max b1 field - "%s".""" % max_b1_field

        # zero_padding    
        if not msg:
            zero_padding = self.TextZeroPadding.GetValue().strip()
            if not util_misc.is_intable(zero_padding):
                msg = """I don't understand the zero padding - "%s".""" % zero_padding

        # min_dwell_time    
        if not msg:
            min_dwell_time = self.TextMinDwellTime.GetValue().strip()
            if not util_misc.is_floatable(min_dwell_time):
                msg = """I don't understand the min dwell time - "%s".""" % min_dwell_time

        # dwell_time_increment    
        if not msg:
            dwell_time_increment = self.TextDwellTimeIncrement.GetValue().strip()
            if not util_misc.is_floatable(dwell_time_increment):
                msg = """I don't understand the dwell time increment - "%s".""" % dwell_time_increment

        # gradient_raster_time    
        if not msg:
            gradient_raster_time = self.TextGradientRasterTime.GetValue().strip()
            if not util_misc.is_floatable(gradient_raster_time):
                msg = """I don't understand the gradient raster time - "%s".""" % gradient_raster_time

        # gradient_slew_rate    
        if not msg:
            gradient_slew_rate = self.TextGradientSlewRate.GetValue().strip()
            if not util_misc.is_floatable(gradient_slew_rate):
                msg = """I don't understand the gradient slew rate - "%s".""" % gradient_slew_rate

        # gradient_maximum
        if not msg:
            gradient_maximum = self.TextGradientMaximum.GetValue().strip()
            if not util_misc.is_floatable(gradient_maximum):
                msg = """I don't understand the gradient maximum - "%s".""" % gradient_maximum

        if not msg and self.is_template and not self.machine_specs.id:
            # I still have yet to check whether or not the name is unique.
            # I save it until last because it is a database hit.
            names = [template.name for template in self.templates]
            if name in names:
                msg = "A machine specs template with this name already exists.\n"  \
                      "Please enter a unique name."

        if msg:
            # validation failed
            common_dialogs.message(msg, None, common_dialogs.I_OK)
        else:
            # All is well!
            if self.is_template:
                self.machine_specs.name                 = name
            self.machine_specs.machine_type         = machine_type
            self.machine_specs.field_strength       = float(field_strength)
            self.machine_specs.max_b1_field         = float(max_b1_field)
            self.machine_specs.zero_padding         = int(zero_padding)
            self.machine_specs.min_dwell_time       = float(min_dwell_time)
            self.machine_specs.dwell_time_increment = float(dwell_time_increment)
            self.machine_specs.gradient_raster_time = float(gradient_raster_time)
            self.machine_specs.gradient_slew_rate   = float(gradient_slew_rate)
            self.machine_specs.gradient_maximum     = float(gradient_maximum)
                                
        return not bool(msg)






