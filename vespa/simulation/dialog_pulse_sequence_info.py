# Python modules


# 3rd party modules
import wx


# Our modules
import vespa.simulation.auto_gui.pulse_sequence_info as pulse_sequence_info
import vespa.simulation.dialog_experiment_list as dialog_experiment_list
import vespa.common.mrs_pulse_sequence as mrs_pulse_sequence
import vespa.common.util.time_ as util_time
import vespa.common.util.misc as util_misc
import vespa.common.constants as common_constants
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


class DialogPulseSequenceInfo(pulse_sequence_info.MyDialog):
    def __init__(self, parent, db, pulse_sequence=None, read_only=False):
        if not parent:
            parent = wx.GetApp().GetTopWindow()
            
        pulse_sequence_info.MyDialog.__init__(self, parent)
        
        self.pulse_sequence = pulse_sequence
        
        if self.pulse_sequence:
            self.read_only = read_only
            
            title = "Pulse Sequence " + self.pulse_sequence.name
            
            states = [ ]
            if self.pulse_sequence.is_frozen:
                states.append("frozen")
            else:
                states.append("not frozen")
            if self.pulse_sequence.experiment_names:
                states.append("in use")
            else:
                states.append("not in use")
            if self.pulse_sequence.is_public:
                states.append("public")
            else:
                states.append("private")

            if states:
                title += " (" + ", ".join(states) + ")"
        else:
            # This is a pulse sequence
            # When creating a new pulse sequence, read_only can't be True.
            self.read_only = False
            # Create a new metab object
            title = "New Pulse Sequence"
            self.pulse_sequence = mrs_pulse_sequence.PulseSequence()
            self.pulse_sequence.id = util_misc.uuid()
            # This non-standard attribute must exist for all pulse sequences 
            # in this dialog. 
            self.pulse_sequence.experiment_names = [ ]

        self.SetTitle(title)

        # This db handle is only used to check that a pulse sequence name is 
        # unique. It's not for writing!
        self.db = db

        # I need a reference to this sizer for adding and removing controls
        self.ParameterSizer = self.LabelParametersPlaceholder.GetContainingSizer()
        self.LabelParametersPlaceholder.Destroy()
                
        # I create loop label aliases that are easier to work with than the 
        # names that wxGlade assigns.
        self.loop_labels = [ ]
        for i in range(1, common_constants.RESULTS_SPACE_DIMENSIONS):
            self.loop_labels.append(getattr(self, "TextLoop%dLabel" % i))

        # parameter_controls is a list of lists. The inner lists contain
        # one row each of parameter controls.
        self.parameter_controls = [ ]
        
        self.initialize_controls()

        self.Fit()

        self.Center()



    ##################    EVENT HANDLERS     ##################

    def on_ok(self, event):
        # Validate controls one by one
        msg = ""

        # name
        name = self.TextName.GetValue().strip()        
        if not name:
            msg = "Please enter a name for the pulse sequence."

        # creator
        if not msg:
            creator = self.TextCreator.GetValue().strip()        

        # comment
        if not msg:
            # Note that I don't call .strip() on this field. It's totally
            # freeform and whitespace might be important.
            comment = self.TextComment.GetValue()
            
        # parameters.    
        if not msg:
            # Remember that rows that the user has removed are merely hidden.
            i = 0
            for row in self.parameter_controls:
                _, combobox, _, text_name, _, text_value = row
                if combobox.IsShown():
                    # This row is visible and might contain data.
                    i += 1
                    combobox = bool(combobox.GetStringSelection())
                    text_name = bool(text_name.GetValue())
                    text_value = bool(text_value.GetValue())

                    if (combobox != text_name) or (text_name != text_value):
                        msg = "Please complete parameter %d." % i
                #else:
                    # This row was "removed" by the user

        # sequence code.    
        if not msg:
            # Since this is someone's code, whitespace might be important
            # and I don't want to strip() it. However, I want to detect 
            # the situation where a user has entered only whitespace, so
            # it's appropriate to strip() there.
            sequence_code = self.TextSequenceCode.GetValue()
            
            if not sequence_code.strip():
                msg = "Please enter some sequence code."

        # binning code.    
        if not msg:
            binning_code = self.TextBinningCode.GetValue().strip()
            
            if not binning_code.strip():
                msg = "Please enter some binning code."

        if not msg:
            # I still have yet to check whether or not this name is unique.
            # I save it until last because it is a database hit.
            pulse_sequences = self.db.fetch_pulse_sequences_by_name(name)
            ids = [pulse_sequence.id for pulse_sequence in pulse_sequences \
                               if pulse_sequence.id != self.pulse_sequence.id]

            if ids:
                msg = "A pulse sequence with this name already exists."

        if msg:
            # validation failed
            common_dialogs.message(msg, None, common_dialogs.X_OK)
        else:
            # All is well!
            self.pulse_sequence.name = name
            self.pulse_sequence.creator = creator
            self.pulse_sequence.comment = comment
            self.pulse_sequence.sequence_code = sequence_code
            self.pulse_sequence.binning_code = binning_code

            loop_labels = [label.GetValue().strip() for label in \
                                                            self.loop_labels]
            # Remove blanks
            self.pulse_sequence.loop_labels = \
                                    [label for label in loop_labels if label]

            self.pulse_sequence.user_static_parameters = [ ]
            for row in self.parameter_controls:
                _, combobox, _, text_name, _, text_value = row
                if combobox.IsShown():
                    # This row is visible and might contain data.
                    _type = combobox.GetStringSelection()

                    if _type:
                        parameter = mrs_pulse_sequence.UserStaticParameter()
                        parameter.type = _type
                        parameter.name = text_name.GetValue().strip()
                        parameter.default = text_value.GetValue()
                        
                        self.pulse_sequence.user_static_parameters.append(parameter)
                #else:
                    # This row was "removed" by the user

                self.EndModal(wx.ID_OK)


    def on_add(self, event):
        sizer = self.ParameterSizer

        sizer.SetRows(len(self.parameter_controls) + 1)

        self.add_parameter_row()
        row = self.parameter_controls[-1]

        checkbox, combobox, label_name, text_name, label_value, text_value = row
        sizer.Add(checkbox,    0, wx.ALIGN_CENTER_VERTICAL|wx.RIGHT, 5)
        sizer.Add(combobox,    0, wx.RIGHT, 5)
        sizer.Add(label_name,  0, wx.LEFT|wx.ALIGN_CENTER_VERTICAL, 5)
        sizer.Add(text_name,   1, wx.LEFT|wx.RIGHT|wx.EXPAND, 5)
        sizer.Add(label_value, 0, wx.LEFT|wx.ALIGN_CENTER_VERTICAL, 5)
        sizer.Add(text_value,  1, wx.LEFT|wx.RIGHT|wx.EXPAND, 5)

        sizer.Layout()
        self.Layout()
        self.Fit()
        
        self.correct_tab_order()


    def on_remove(self, event):
        # From the Crude-But-Effective Dept. -- it turns out that it is
        # much easier to hide these controls than it is to remove them.
        # Removing them isn't difficult, but doing so renders inaccurate the
        # indices passed to the controls' validators when they were created.
        at_least_one_row_checked = False
        for i, row in enumerate(self.parameter_controls):
            checkbox, combobox, _, text_name, _, text_value = row
            if checkbox.IsChecked():
                at_least_one_row_checked = True
                # Empty the controls of input
                checkbox.SetValue(False)
                combobox.SetSelection(-1)
                text_name.ChangeValue("")
                text_value.ChangeValue("")
                for control in row:
                    control.Hide()

        if not at_least_one_row_checked:
            common_dialogs.message("Please check the parameter you want to remove.", style=common_dialogs.I_OK)
        else:
            self.ParameterSizer.Layout()
            self.Fit()
            
            
    def on_experiments(self, event):
        dialog = dialog_experiment_list.DialogExperimentList(self, self.pulse_sequence)
        dialog.ShowModal()


    ##################    Internal helper functions     ##################

    def add_parameter_row(self):
        index = len(self.parameter_controls)
        
        # These controls need to be children of the same window as the other
        # controls. I think that's the notebook, but there's no point to me
        # guessing when I can just ask wx.
        parent = self.TextName.GetParent()

        types = ("Double", "Long", "String")
        row = [ ]
        checkbox = wx.CheckBox(parent)
        row.append(checkbox)
        if self.read_only or self.pulse_sequence.is_frozen:
            # In read-only mode, the checkboxes are useless. I hide them
            # rather than destroying them because the grid sizer still
            # expects a certain number of columns per row.
            checkbox.Hide()

        combobox = wx.ComboBox(parent, style=wx.CB_DROPDOWN|wx.CB_READONLY,
                               size=(100, -1), choices=types)
        if self.read_only or self.pulse_sequence.is_frozen:
            combobox.Disable()
        row.append(combobox)

        row.append(wx.StaticText(parent, wx.ID_ANY, "Name:"))

        textbox = wx.TextCtrl(parent)
        if self.read_only or self.pulse_sequence.is_frozen:
            textbox.Disable()
        row.append(textbox)

        row.append(wx.StaticText(parent, wx.ID_ANY, "Default Value:"))

        textbox = wx.TextCtrl(parent)
        if self.read_only or self.pulse_sequence.is_frozen:
            textbox.Disable()
        row.append(textbox)

        self.parameter_controls.append(row)

        return row


    def initialize_controls(self):
        # We add the OK & Cancel buttons dynamically so that they're in the 
        # right order under OS X, GTK, Windows, etc.
        self.ButtonOk, self.ButtonCancel = \
                    wx_util.add_ok_cancel(self, self.LabelOkCancelPlaceholder,
                                          self.on_ok)
        
        self.TextName.ChangeValue(self.pulse_sequence.name)

        self.LabelUuid.SetLabel(self.pulse_sequence.id)
        self.TextCreator.ChangeValue(self.pulse_sequence.creator)
        created = self.pulse_sequence.created.strftime(util_time.DISPLAY_DATE_FORMAT)
        self.LabelCreated.SetLabel(created)                
        self.TextComment.ChangeValue(self.pulse_sequence.comment)

        for i, label in enumerate(self.pulse_sequence.loop_labels):
            self.loop_labels[i].ChangeValue(label)

        # Read-only and frozen are almost the same except frozen permits 
        # editing of name & comment
        if self.read_only or self.pulse_sequence.is_frozen:
            controls = (self.TextName, self.TextCreator,
                        self.TextComment)

            for control in controls:
                control.Disable()

            for label in self.loop_labels:
                label.Disable()
                
        if not self.read_only:
            # Ensure these are enabled
            self.TextName.Enable()
            self.TextComment.Enable()

        # Here I figure out how many rows of parameter controls I need.
        # I create one row of param controls for each param, or one empty
        # row if there are no params & we're in edit mode.
        if self.pulse_sequence.user_static_parameters or self.read_only:
            for parameter in self.pulse_sequence.user_static_parameters:
                self.on_add(None)
            
                row = self.parameter_controls[-1]

                _, combobox, _, text_name, _, text_value = row
                
                combobox.SetStringSelection(parameter.type)
                text_name.ChangeValue(parameter.name)
                text_value.ChangeValue(parameter.default)
        else:
            # Create a single blank row
            self.on_add(None)
            

        # Sequence & Binning Code tabs 
        # Set font to monospace
        wx_util.set_font_to_monospace(self.TextSequenceCode)
        wx_util.set_font_to_monospace(self.TextBinningCode)
        
        code = ""
        if self.pulse_sequence.sequence_code:
            code = self.pulse_sequence.sequence_code
        self.TextSequenceCode.ChangeValue(code)
        if self.read_only or self.pulse_sequence.is_frozen:
            self.TextSequenceCode.SetEditable(False)
        
        code = ""
        if self.pulse_sequence.binning_code:
            code = self.pulse_sequence.binning_code
        self.TextBinningCode.ChangeValue(code)
        if self.read_only or self.pulse_sequence.is_frozen:
            self.TextBinningCode.SetEditable(False)

        if self.read_only:
            self.ButtonAdd.Hide()
            self.ButtonRemove.Hide()
            self.ButtonOk.Hide()
            self.ButtonCancel.SetLabel("Done")
            self.ButtonCancel.SetFocus()

        if self.pulse_sequence.is_frozen:
            self.ButtonAdd.Hide()
            self.ButtonRemove.Hide()
            
        if not self.pulse_sequence.experiment_names:
            # No point in showing this button which will only display an
            # empty list.
            self.ButtonListExperiments.Hide()

        # These text controls have a habit of sometimes assigning themselves
        # a huge width and/or height. Here we ensure that they have a sane size.
        screen_width = wx.SystemSettings_GetMetric(wx.SYS_SCREEN_X)
        screen_height = wx.SystemSettings_GetMetric(wx.SYS_SCREEN_Y)

        width = int(screen_width * 0.5)
        height = int(screen_height * 0.5)

        self.TextSequenceCode.SetSize( (width, height) )
        self.TextBinningCode.SetSize( (width, height) )



    def correct_tab_order(self):
        # wx sets the tab order based on the order in which controls were
        # created. wxGlade writes the control creation statements in an order 
        # that results in logical tabbing. When controls are created after
        # dialog init (as they are in this dialog), that messes up the tab
        # order and it needs to be corrected. 
        
        # Adding to the quirkiness here is the fact that many of the controls
        # in the lower half of this dialog are not always useful (e.g.
        # "show experiment list" when there are no experiments to show). 

        # The correct tab order is Add, Remove, param controls left to right,
        # top to bottom, List Experiments, leftmost OK/Cancel button, 
        # rightmost OK/Cancel button.
        
        controls = [ ]
        for row in self.parameter_controls:
            controls += row
        controls.append(self.ButtonListExperiments)
        
        # The location of OK & Cancel is OS-dependent.
        ok_left, _ = self.ButtonOk.GetPosition()
        cancel_left, _ = self.ButtonCancel.GetPosition()
        
        if (ok_left == -1) and (cancel_left == -1):
            # This happens under GTK during __init__, and on that platform I 
            # know that cancel is on the left.
            cancel_left = ok_left - 1
        
        if ok_left > cancel_left:
            controls.append(self.ButtonCancel)
            controls.append(self.ButtonOk)
        else:
            controls.append(self.ButtonOk)
            controls.append(self.ButtonCancel)
            
        previous = self.ButtonRemove
        for control in controls:
            control.MoveAfterInTabOrder(previous)
            previous = control

