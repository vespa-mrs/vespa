# Python modules
import math

# 3rd party modules
import wx
#import wx.aui as aui
import wx.lib.agw.aui as aui        # NB. wx.aui version throws odd wxWidgets exception on Close/Exit ?? Not anymore in wxPython 4.0.6 ??
import wx.stc as wx_stc
import numpy as np

# Our modules
import vespa.pulse.prefs as prefs
import vespa.pulse.util_pulse_config as rf_config
import vespa.pulse.tab_transform as tab_transform
import vespa.pulse.auto_gui.editor_transform as editor_transform
import vespa.common.rfp_rf_result as rfp_rf_result
import vespa.common.rfp_transform as rfp_transform
import vespa.common.rfp_transform_kernel as rfp_transform_kernel
import vespa.common.styled_text_control as our_stc
import vespa.common.constants as common_constants
import vespa.common.wx_gravy.common_dialogs as common_dialogs
import vespa.common.wx_gravy.util as common_wx_util
import vespa.common.util.misc as util_misc
import vespa.common.util.time_ as util_time
from vespa.common.rfp_machine_specs import MachineSpecs


PI = math.pi

# This is the row of lines that goes in the console to separate the messages
# of one run from another.
SEPARATOR = "-" * 60


def _find_run_button(window):
    # copied from tab_pulse_design
    #
    # Given a window (usually a tab), will walk through the tree of children
    # looking for a button with the text "Run". If it can't find such a 
    # button, it returns None.
    #
    # This function is recursive.
    run_button = None

    subwindows = [ ]
    for kid in window.GetChildren():
        if hasattr(kid, "GetChildren"):
            # This kid has kids of its own that we might need to check
            subwindows.append(kid)

        if hasattr(kid, "GetLabel") and isinstance(kid, wx.Button):
            # Looks like a button...
            if kid.GetLabel() == "Run":
                # It is a button, and it's the one we're looking for.
                run_button = kid
           
    # As long as we haven't found our target yet and there's still kids to
    # traverse, we keep looking.
    while (not bool(run_button)) and bool(subwindows):
        run_button = _find_run_button(subwindows.pop())
        
    return run_button


class LocalPulseDesign(object):

    def __init__(self, machine_specs):

        self.calc_resolution      = common_constants.DEFAULT_CALC_RESOLUTION       #5000 # num_points for Bloch and rfpulse profile calcs.
        self.pulse_bandwidth_type = common_constants.DEFAULT_PULSE_BANDWIDTH_TYPE  #'half_height'   # 'half_height', 'minimum' or 'maximum'
        self.gyromagnetic_nuclei  = common_constants.DEFAULT_GYROMAGNETIC_NUCLEI
        
        self.bloch_range_value    = common_constants.DEFAULT_BLOCH_RANGE_VALUE
        self.bloch_range_units    = common_constants.DEFAULT_BLOCH_RANGE_UNITS
        self.bloch_offset_value   = common_constants.DEFAULT_BLOCH_OFFSET_VALUE
        if machine_specs is None:
            machine_specs = MachineSpecs()
            
        self.machine_specs = machine_specs

    def get_previous_result(self, transform):
        """ for equivalent functionality as tab_pulse_design """
        if transform.type in ['Create Transform',]:
            return None
        else:
            # create previous RF waveform - simple hamming-sinc
            npts      = 128     # int
            duration  = 1.0     # float, msec
            ncycles   = 6
            
            dwell_time = 0.001 * float(duration)/float(npts)
            
            t = (np.arange(npts) - (npts/2.0)) / (npts/2.0)   # time, symmetric about zero
            y = 2 * np.pi * ncycles * t + 0.00001   # 0.0001 to avoid div by zero
            b1 = np.sin(y) / y
            
            # apply a hamming filter
            b1 = b1 * 4 * ncycles * (0.54 + 0.46*np.cos(np.pi*t)) / npts
            
            b1x = np.arange(npts) * dwell_time      # time from zero in sec
            grx = np.arange(npts) * dwell_time      # time from zero in sec
            
            # scale empirically to 90 degree tip
            b1 = 0.1266667 * b1/np.max(b1)
            gr = None  #np.ones(npts)
        
            b1 = np.array(b1)
            if not np.iscomplexobj(b1):
                b1r = b1.real.copy()
                b1i = b1r.copy() * 0
                b1 = b1r + 1j*b1i            

            result = rfp_rf_result.RfResults()
            result.rf_waveform = np.array(b1)
            result.rf_xaxis    = np.array(b1x)
            if gr is not None:
                gr = np.array(gr)
            if grx is not None:
                grx = np.array(grx)
            result.gradient    = gr
            result.grad_xaxis  = grx            
            
            return result


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


class DialogEditorTransform(editor_transform.MyDialog):

    def __init__(self, parent, all_names, all_labels, transform_kernel=None, machine_specs=None):

        config = rf_config.Config()
        position, size = config.get_window_coordinates("dialog_editor_transform")
        self._left, self._top = position
        self._width, self._height = size

        if not parent:
            parent = wx.GetApp().GetTopWindow()

        editor_transform.MyDialog.__init__(self, parent)
        
        self.parent             = parent
        self.all_names          = all_names
        self.all_labels         = all_labels
        self.test_transform     = rfp_transform.Transform()
        self.transform_kernel   = transform_kernel

        self.test_transform.reset_result()
        
        # this dialog is launched from the application rather than from a 
        # given pulse_design, so most times we will use default specs 
        # for these two global objects
        self.pulse_design = LocalPulseDesign(machine_specs)        
        
        if self.transform_kernel:
            
            title = transform_kernel.type + " - " + self.transform_kernel.name

            states = [ ]
            if self.transform_kernel.is_frozen:
                states.append("frozen")
            else:
                states.append("not frozen")
            if self.transform_kernel.is_public:
                states.append("public")
            else:
                states.append("private")

            if states:
                title += " (" + ", ".join(states) + ")"
        else:
            self.transform_kernel = rfp_transform_kernel.TransformKernel()
            self.transform_kernel.id = util_misc.uuid()
            title = "New Create Transform"
            
        self.SetTitle(title)

           
        #------------------------------
        # Setup plot parameters

        self.foreground = "black"
        self.background = "white"
        self._prefs = prefs.PrefsMain()

        #------------------------------
        # Initialize widget controls
        
        self.tab_transform = None       # this is same Tab as in main program
        self.initialize_controls()

        self.test_transform.reset_parameters()
        self.test_transform.reset_result()

        # Set dialog size & position from the INI file.
        self.SetSize( (self._width, self._height) )
        self.SetPosition( (self._left, self._top) )

        # The sash (a.k.a. grabber) on Design needs to be 550 wide to show the entire
        # left-hand pane under OS X. Under Windows & Ubuntu, it needs less.
        sash_position = 550
        self.SplitDesign.SetSashPosition(sash_position, True)

        # The sash (a.k.a. grabber) on Test needs to be 550 wide to show the entire
        # left-hand pane under OS X. Under Windows & Ubuntu, it needs less.
        sash_position = 650
        self.SplitTest.SetSashPosition(sash_position, True)  
        
        #self.Center()

        self.bind_events()


    def bind_events(self):
        self.Bind(wx.EVT_SIZE, self.on_self_coordinate_change)
        self.Bind(wx.EVT_MOVE, self.on_self_coordinate_change)


    def on_self_coordinate_change(self, event):
        # This is invoked for move & size events
        if self.IsMaximized() or self.IsIconized():
            # Bah, forget about this. Recording coordinates doesn't make sense
            # when the window is maximized or minimized. This is only a
            # concern on Windows; GTK and OS X don't produce move or size
            # events when a window is minimized or maximized.
            pass
        else:
            if event.GetEventType() == wx.wxEVT_MOVE:
                self._left, self._top = self.GetPosition()
            else:
                # This is a size event
               self._width, self._height = self.GetSize()
        event.Skip()


    ##### Event Handlers ######################################################

    def on_add_user_parameter(self, event=None): 
        # This event is called both from wx as a proper event and from 
        # our code in which case event will be None. Don't rely on it being
        # present!
        
        # sizer is a GridSizer, thus we add a row to it and fill in the
        # columns with widgets to describe a user parameter
        
        sizer = self.sizer_design_parameters
        sizer.SetRows(len(self.design_parameter_controls) + 1)

        self.design_add_parameter_row()
        row = self.design_parameter_controls[-1]

        checkbox, combobox, label_name, text_name, label_value, text_value, label_var, text_var = row
        sizer.Add(checkbox,    0, wx.ALIGN_CENTER_VERTICAL|wx.RIGHT, 5)
        sizer.Add(combobox,    0, wx.RIGHT, 5)
        sizer.Add(label_name,  0, wx.LEFT|wx.ALIGN_CENTER_VERTICAL, 5)
        sizer.Add(text_name,   1, wx.LEFT|wx.RIGHT|wx.EXPAND, 5)
        sizer.Add(label_value, 0, wx.LEFT|wx.ALIGN_CENTER_VERTICAL, 5)
        sizer.Add(text_value,  1, wx.LEFT|wx.RIGHT|wx.EXPAND, 5)
        sizer.Add(label_var,   0, wx.LEFT|wx.ALIGN_CENTER_VERTICAL, 5)
        sizer.Add(text_var,    1, wx.LEFT|wx.RIGHT|wx.EXPAND, 5)

        self.PanelDesignControls.Layout()
        
        self.correct_tab_order()
        

    def on_remove_user_parameter(self, event):
        # Remove all checked rows from the user user param definition 
        # controls.
        # I loop over the rows of controls but I can't delete them right
        # away because it's not kosher to alter a list while one is
        # iterating over it.
        delete_these_rows = [ ]
        for i, row in enumerate(self.design_parameter_controls):
            checkbox, _, _, _, _, _, _, _ = row
            if checkbox.IsChecked():
                delete_these_rows.append(i)

        # Now I can clean up the list. I do it from the highest index to the
        # lowest. If I did it lowest to highest, then deleting the first
        # item would throw off the indices of the subsequent items.
        delete_these_rows.reverse()

        for row in delete_these_rows:
            for control in self.design_parameter_controls[row]:
                # Been nice knowin' ya!
                control.Destroy()
            del self.design_parameter_controls[row]
            
        if delete_these_rows:
            # Refresh the grid.
            self.sizer_design_parameters.SetRows(len(self.design_parameter_controls))
            self.PanelDesign.Layout()
        else:
            common_dialogs.message("Please check the parameter you want to remove.", style=common_dialogs.I_OK)


    def on_page_changed(self, event): 
        
        index = self.NotebookTools.GetSelection()
        
        # we check the designer to see if the current set of parameters
        # pass the validation step, and if so, we update the tester widget
        
        if index == 1:
            if self.validate_design_transform():
                new_kernel = self.transform_kernel.clone()
                self.tab_transform.update_kernel_controls(new_kernel)
                
            else:
                self.NotebookTools.SetSelection(0)


    def on_results(self, event):
        # pops up a textual description of the test Experiment and results
        # in a native text editor.  This allows a user to keep a copy of a
        # previous test run to compare to subsequent runs.
        self.display_results_text()


    def on_ok(self, event):
        if self.validate_design_transform():
            # Save my coordinates
            config = rf_config.Config()
            config.set_window_coordinates("dialog_editor_transform", self._left, self._top, self._width, self._height)
            config.write()
            # close dialog
            self.EndModal(wx.ID_OK)


    ##### Internal helper functions  ##########################################

    def initialize_controls(self):
        # sets up widgets in designer tab and test tab (including user control 
        # values) from values in either a new or existing transform kernel.
        #
        # finally if the pulse sequence is 'frozen' most of the widgets 
        # are disabled since only the name and comment should be editable in
        # that state

        # clear the labels that act as status bar locations
        self.LabelStatus0.SetLabel('')
        self.LabelStatus1.SetLabel('')
        self.LabelStatus2.SetLabel('')
        self.LabelStatus3.SetLabel('')
        
        # We add the OK & Cancel buttons dynamically so that they're in the 
        # right order under OS X, GTK, Windows, etc.
        self.ButtonOk, self.ButtonCancel = common_wx_util.add_ok_cancel(self, self.LabelOkCancelPlaceholder, self.on_ok)
        
        self.initialize_design_controls()
        self.initialize_test_controls()

        # FIXME bjs - decide which things need to be disabled for FROZEN state

        if self.transform_kernel.is_frozen:
            controls = (self.TextCreator,
                        self.ButtonAddUserParameter, 
                        self.ButtonRemoveUserParameter,
                        self.PanelTest)

            for control in controls:
                control.Disable()
            
            msg = "Welcome to the Pulse Transform Editor.\n" + \
                  "\nThe transform that you are editing is 'frozen' so only "+ \
                  "the name and comment fields can be edited.\n" 
        else:
            msg = "Welcome to the Pulse Transform Editor.\n" + \
                  "\nStatus messages will appear here as you test your code.\n"

        self.TextConsole.Clear()
        self.TextConsole.AppendText(msg)


    def initialize_test_controls(self):

        self.tab_transform = tab_transform.TabTransform(self, self, 
                                                        self.pulse_design,
                                                        self.test_transform, 
                                                        parent=self.PanelTransformBase,
                                                        in_editor=True)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.tab_transform, 1, wx.EXPAND)
        self.PanelTransformBase.SetSizer(sizer)

        self.PanelTest.Layout()


    def initialize_design_controls(self):

        kernel = self.transform_kernel
        
        # a persistent reference to this sizer to add/remove controls
        self.sizer_design_parameters = self.LabelPlaceholder1.GetContainingSizer()
        self.LabelPlaceholder1.Destroy()
        
        # I force this to a minimum size so that it is at least wide enough
        # to be usable. The size I chose is arbitrary. 
        self.TextComment.SetMinSize( (350, -1) )        

        # Global Parameters panel
        
        self.TextName.ChangeValue(kernel.name)
        self.TextMenuLabel.ChangeValue(kernel.menu_label)
        self.LabelUuid.SetLabel(kernel.id)
        self.TextCreator.ChangeValue(kernel.creator)
        created = kernel.created.strftime(util_time.DISPLAY_DATE_FORMAT)
        self.LabelCreated.SetLabel(created)
        
        self.ChoiceType.Clear()
        for item in ['Create Transform','Modify Transform','Visualize Transform']:
            self.ChoiceType.Append(item)
        self.ChoiceType.SetStringSelection(kernel.type)
        
        self.TextComment.ChangeValue(kernel.comment)

        self.CheckFile1.SetValue(kernel.hide_file1)
        self.CheckFile2.SetValue(kernel.hide_file2)
        self.TextFile1Label.ChangeValue(str(kernel.file1_label))
        self.TextFile2Label.ChangeValue(str(kernel.file2_label))

        self.TextTipAngle.ChangeValue(str(kernel.tip_angle))
        self.TextTimeSteps.ChangeValue(str(kernel.time_steps))
        self.TextDuration.ChangeValue(str(kernel.duration))
        self.TextBandwidth.ChangeValue(str(kernel.bandwidth))
        
        self.CheckTimeSteps.SetValue(kernel.hide_time_steps)
        self.CheckDuration.SetValue(kernel.hide_duration)
        self.CheckTipAngle.SetValue(kernel.hide_tip_angle)
        self.CheckBandwidth.SetValue(kernel.hide_bandwidth)

        # User Defined Parameters panel
        #
        # design_parameter_controls is a list of lists. The inner lists contain
        # one row each of parameter controls.
        self.design_parameter_controls = []

        # Here I hide all existing controls and then figure out 
        # how many rows of parameter controls I need for the update.
        # I create one row of param controls for each param, or one empty
        # row if there are no params & we're in edit mode.
        for i, parameter in enumerate(self.design_parameter_controls):
            row = self.design_parameter_controls[i]
            checkbox, combobox, _, text_name, _, text_value = row
            checkbox.SetValue(False)
            combobox.SetSelection(-1)
            text_name.ChangeValue("")
            text_value.ChangeValue("")
            for control in row:
                control.Hide()
        
        if kernel.transform_kernel_controls:
            for control in kernel.transform_kernel_controls:
                if control.variable not in ['time_steps', 'duration', 'tip_angle', 'bandwidth']:
                    self.on_add_user_parameter(None)
                    row = self.design_parameter_controls[-1]
                    _, combobox, _, text_name, _, text_value, _, text_variable = row
                    combobox.SetStringSelection(control.type)
                    text_name.ChangeValue(control.name)
                    text_value.ChangeValue(control.default)
                    text_variable.ChangeValue(control.variable)
        else:
            # Create a single blank row
            self.on_add_user_parameter()

        self.sizer_design_parameters.Layout()
        
        # Algorithm Code Editor panel
        #
        # NB. the End Of Line (EOL) mode below is set to 2 which results in
        # an LF only, the default mode is 0 = CR,LF and mode = 1 is CR only
        self.tab_algorithm_code = our_stc.PythonSTC(self.PanelDesignAlgorithm, -1)
        self.tab_algorithm_code.EmptyUndoBuffer()
        self.tab_algorithm_code.SetEOLMode(2)
        self.tab_algorithm_code.Colourise(0, -1)
        # line numbers in the margin
        self.tab_algorithm_code.SetMarginType(1, wx_stc.STC_MARGIN_NUMBER)
        self.tab_algorithm_code.SetMarginWidth(1, 25)        
        if self.transform_kernel.is_frozen:
            if self.transform_kernel.algorithm_code:
                # need to add code prior to setting 'read only'
                code = self.transform_kernel.algorithm_code
                self.tab_algorithm_code.EmptyUndoBuffer()
                self.tab_algorithm_code.SetText(code)
            else:
                self.tab_algorithm_code.SetText("")
            self.tab_algorithm_code.SetReadOnly(True)
        else:
            self.tab_algorithm_code.SetText("")        

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.tab_algorithm_code, 1, wx.EXPAND)
        self.PanelDesignAlgorithm.SetSizer(sizer)
        
        # put algorithm code into the editor 

        code = ""
        if kernel.algorithm_code:
            code = kernel.algorithm_code
        self.tab_algorithm_code.EmptyUndoBuffer()
        self.tab_algorithm_code.SetText(code)


    ##### Internal Design Tab helper functions  ##########################################

    def design_add_parameter_row(self):
        index = len(self.design_parameter_controls)
        
        # These controls need to be children of the same window as the other
        # controls. I think that's the notebook, but there's no point to me
        # guessing when I can just ask wx.
        parent = self.TextName.GetParent()

        types = ("Double", "Long", "String", "Choice", "Output")
        row = [ ]
        checkbox = wx.CheckBox(parent)
        row.append(checkbox)
        if self.transform_kernel.is_frozen:
            # In read-only mode, the checkboxes are useless. I hide them
            # rather than destroying them because the grid sizer still
            # expects a certain number of columns per row.
            checkbox.Hide()

        combobox = wx.ComboBox(parent, style=wx.CB_DROPDOWN|wx.CB_READONLY,
                               size=(-1, -1), choices=types)
        # Under OS X, wx2.9 selects the first item in the list automatically.
        combobox.SetSelection(-1)
        if self.transform_kernel.is_frozen:
            combobox.Disable()
        row.append(combobox)

        row.append(wx.StaticText(parent, wx.ID_ANY, "Name:"))

        textbox = wx.TextCtrl(parent)
        if self.transform_kernel.is_frozen:
            textbox.Disable()
        row.append(textbox)

        row.append(wx.StaticText(parent, wx.ID_ANY, "Default:"))

        textbox = wx.TextCtrl(parent)
        if self.transform_kernel.is_frozen:
            textbox.Disable()
        row.append(textbox)

        row.append(wx.StaticText(parent, wx.ID_ANY, "Variable:"))

        textbox = wx.TextCtrl(parent)
        if self.transform_kernel.is_frozen:
            textbox.Disable()
        row.append(textbox)
        
        self.design_parameter_controls.append(row)

        return row

        
    def validate_design_transform(self):
        # Validate controls one by one
        # A non-empty string message indicates failure
        
        msg = ""

        # name
        name = self.TextName.GetValue().strip()        
        if not name:
            msg = "Please enter a name for the transform."
        elif name in self.all_names:
            msg = "Name is not unique, please change."

        # menu_label
        menu_label = self.TextMenuLabel.GetValue().strip()        
        if not menu_label:
            msg = "Please enter a menu label for the transform."
        elif menu_label in self.all_labels:
            msg = "Menu Label is not unique, please change."

        # creator (optional)
        if not msg:
            creator = self.TextCreator.GetValue().strip()        
                
        # comment
        if not msg:
            # Note that I don't call .strip() on this field. It's totally
            # freeform and whitespace might be important.
            comment = self.TextComment.GetValue()

        # type
        if not msg:
            type_ = str(self.ChoiceType.GetStringSelection())

        # Hide certain attributes
        hide_file1       = self.CheckFile1.IsChecked()
        hide_file2       = self.CheckFile2.IsChecked()
        hide_time_steps  = self.CheckTimeSteps.IsChecked()
        hide_duration    = self.CheckDuration.IsChecked()
        hide_tip_angle   = self.CheckTipAngle.IsChecked()
        hide_bandwidth   = self.CheckBandwidth.IsChecked()

        if not msg:
            file1       = ''
            file1_label = self.TextFile1Label.GetValue()
        if not msg:
            file2       = ''
            file2_label = self.TextFile2Label.GetValue()

        # time_steps
        time_steps = self.TextTimeSteps.GetValue().strip()        
        if not time_steps and not hide_time_steps:
            msg = "Please enter a default value for the time steps."
        elif not util_misc.is_intable(time_steps):
            msg = "Please enter an integer value for the tip angle."

        # duration
        duration = self.TextDuration.GetValue().strip()        
        if not duration and not hide_duration:
            msg = "Please enter a default value for the duration."
        elif not util_misc.is_floatable(duration):
            msg = "Please enter a float value for the duration."

        # tip_angle
        tip_angle = self.TextTipAngle.GetValue().strip()        
        if not tip_angle and not hide_tip_angle:
            msg = "Please enter a default value for the tip angle."
        elif not util_misc.is_floatable(tip_angle):
            msg = "Please enter a float value for the tip angle."

        # bandwidth
        bandwidth = self.TextBandwidth.GetValue().strip()        
        if not bandwidth and not hide_bandwidth:
            msg = "Please enter a default value for the bandwidth."
        elif not util_misc.is_floatable(bandwidth):
            msg = "Please enter a float value for the bandwidth."

        # parameters - test completeness    
        if not msg:
            for i, row in enumerate(self.design_parameter_controls):
                _, combobox, _, text_name, _, text_value, _, text_variable = row
                combobox = bool(combobox.GetStringSelection())
                test_name = bool(text_name.GetValue())
                test_value = bool(text_value.GetValue())
                test_variable = bool(text_variable.GetValue())

                if (combobox != test_name) or (test_name != test_value) or (test_value != test_variable):
                    msg = "Please complete User parameter %d, '%s'." % (i+1, text_name.GetValue())
                    break

        # parameters - test value types    
        if not msg:
            for i, row in enumerate(self.design_parameter_controls):
                _, combobox, _, text_name, _, text_value, _, text_variable = row
                combo = combobox.GetStringSelection()
                value = text_value.GetValue()

                if combo == "Double" and not util_misc.is_floatable(value):
                    msg = "Please enter float value for parameter %d." % (i + 1)
                    break
                if combo == "Long" and not util_misc.is_intable(value):
                    msg = "Please enter integer value for parameter %d." % (i + 1)
                    break

        # parameters - unique variable names    
        if not msg:
            vars = []
            for i, row in enumerate(self.design_parameter_controls):
                _, combobox, _, text_name, _, text_value, _, text_variable = row
                vars.append(text_variable.GetValue())
            if len(set(vars)) != len(vars):
                msg = "Variable names are not unique, please check and change."
                

        # algorithm code.    
        if not msg:
            # Since this is someone's code, whitespace might be important
            # and I don't want to strip() it. However, I want to detect 
            # the situation where a user has entered only whitespace, so
            # it's appropriate to strip() there.
            algorithm_code = self.tab_algorithm_code.GetText()

            if not algorithm_code.strip():
                msg = "In order to save/test a transform, there must \n"+\
                      "be Python code in the Algorithm Code window. \n\n"+\
                      "Please enter code for the transform."

        if not msg:
            # I still have yet to check whether or not this name is unique.
            # I save it until last because it is a database hit.
            bob = 10
#FIXME bjs - I pass in transform names and menu_labels to avoid this db hit.        
        
#bjs             pulse_sequences = self.db.fetch_pulse_sequences_by_name(name)
#bjs             ids = [transform.id for transform in pulse_sequences \
#bjs                                if transform.id != self.transform.id]
#bjs 
#bjs             if ids:
#bjs                 msg = "A pulse sequence with this name already exists.\n"+\
#bjs                       "Please change the name to a unique value."

        if msg:
            # validation failed
            common_dialogs.message(msg, None, common_dialogs.X_OK)
            return False
        else:
            # All is well!
            self.transform_kernel.type           = type_
            self.transform_kernel.name           = name
            self.transform_kernel.menu_label     = menu_label
            self.transform_kernel.creator        = creator
            self.transform_kernel.comment        = comment
            
            self.transform_kernel.hide_file1     = hide_file1
            self.transform_kernel.hide_file2     = hide_file2
            self.transform_kernel.hide_time_steps = hide_time_steps
            self.transform_kernel.hide_duration  = hide_duration
            self.transform_kernel.hide_tip_angle = hide_tip_angle
            self.transform_kernel.hide_bandwidth = hide_bandwidth
            
            self.transform_kernel.file1_label    = file1_label
            self.transform_kernel.file2_label    = file2_label
            self.transform_kernel.file1          = file1
            self.transform_kernel.file2          = file2
            self.transform_kernel.tip_angle      = tip_angle
            self.transform_kernel.time_steps     = time_steps
            self.transform_kernel.duration       = duration
            self.transform_kernel.bandwidth      = bandwidth
            
            self.transform_kernel.algorithm_code = algorithm_code

            self.transform_kernel.transform_kernel_controls = []
            
            globs = [[hide_time_steps,  'Long', 'Time Steps [int]', time_steps, 'time_steps'],
                     [hide_duration,  'Double', 'Duration [msec]',  duration,   'duration'],
                     [hide_tip_angle, 'Double', 'Tip Angle [deg]',  tip_angle,  'tip_angle'],
                     [hide_bandwidth, 'Double', 'Bandwidth [kHz]',  bandwidth,  'bandwidth']]
            
            for row in globs:
                if not row[0]:
                    control = rfp_transform_kernel.TransformKernelControl()
                    control.type = row[1]
                    control.name = row[2]
                    control.default = row[3]
                    control.variable = row[4]
                    self.transform_kernel.transform_kernel_controls.append(control)
            
            for row in self.design_parameter_controls:
                _, combobox, _, text_name, _, text_value, _, text_variable = row
                # Not all rows contain data. Don't forget to ignore empty ones.
                type_ = combobox.GetStringSelection()
                
                if type_:
                    control = rfp_transform_kernel.TransformKernelControl()
                    control.type = type_
                    control.name = text_name.GetValue().strip()
                    control.default = text_value.GetValue()
                    control.variable = text_variable.GetValue()
                    self.transform_kernel.transform_kernel_controls.append(control)
                                
            return True


    def correct_tab_order(self):
        # wx sets the tab order based on the order in which controls were
        # created. wxGlade writes the control creation statements in an order 
        # that results in logical tabbing. When controls are created after
        # dialog init (as they are in this dialog), that messes up the tab
        # order and it needs to be corrected. 
        
        # Adding to the quirkiness here is the fact that many of the controls
        # in the lower half of this dialog are not always useful (e.g.
        # "show experiment list" when there are no experiments to show). 

        # The correct tab order is Add Param, Remove Param, param controls 
        # left to right, top to bottom, Add Pulse, Remove Pulse, pulse
        # checkboxes top to bottom, List Experiments, leftmost OK/Cancel 
        # button, rightmost OK/Cancel button.
        
        controls = [ ]
        for row in self.design_parameter_controls:
            controls += row
            
        previous = self.ButtonRemoveUserParameter
        for control in controls:
            control.MoveAfterInTabOrder(previous)
            previous = control


    ##### Internal Test Tab helper functions  ##########################################
    #
    # note. these are designed to have functionality equivalent to that in the
    #  tab_pulse_design "inner_notebook" object in the main program
    
    def active_tab(self):
        """ for equivalent functionality as tab_pulse_design """
        return self.tab_transform
    
    def activate_tab(self, tab):
        """ for equivalent functionality as tab_pulse_design """
        pass

    def update_sync_status(self, tab=None):
        """ for equivalent functionality as tab_pulse_design """
        pass
    
    def get_left_neighbor(self, tab):
        """ for equivalent functionality as tab_pulse_design """
        return None

    def get_current_calc_resolution(self):
        # Note. in main program, this value would be taken direct from widget
        # tab = self._inner_notebook.GetPage(0)
        # calc_resolution = tab.TextCalcResolution.GetValue()
        return self.pulse_design.calc_resolution
    
    def run(self, event): 

        #------------------------------------------------------------
        # here is where I inserted code for tab_pulse_design.run()
        # so as to work the same as when in main program.
        # This is a truncated version since I have only one tab here.

        # Set the focus to the run button if possible. This is a fix for 
        # RFPulse bug 19:
        
        
        active_tab = self.tab_transform
        
        run_button = _find_run_button(active_tab)
        if run_button:
            run_button.SetFocus()
   
        tabs = [active_tab,]
        right_neighbors = []

        if tabs:
            # First I have to ensure that each tab's input is valid
            gui_is_valid = True
            for tab in tabs:
                if not tab.validate_gui():
                    gui_is_valid = False
                    break
                
            if gui_is_valid:

                run_approved = True

                if run_approved:
                        
                    # Next I have to move the data from the GUI into the 
                    # objects.
                    for tab in tabs:
                        tab.accept_gui_data()
                
                    # Run 'em.
                    wx.SetCursor(wx.HOURGLASS_CURSOR)
                    for tab in tabs:
                        if not tab.run():
                            break

                wx.SetCursor(wx.NullCursor)
        else:
            common_dialogs.message("These results are already up to date.")            
            
            
        # this is the end of the tab_pulse_design.run() code 
        #----------------------------------------------------------------------
            


    def display_results_text(self):
    
        lines = str(self.experiment)
        lines += "\n\nSimulation Results\n" + "-" * 75 + "\n\n"

        lines += "\n".join([simulation.summary() for simulation 
                                                 in self.experiment.simulations])
            
        common_wx_util.display_text_as_file(lines)


    def say(self, message):
        # Count leading newlines
        i = 0
        for c in message:
            if c != "\n":
                break
            i += 1
            
        # Detach those newlines so I can place them before the timestamp
        newlines = message[:i]
        message = message[i:]
        
        timestamp = str(util_time.now())
        self.TextConsole.AppendText("%s%s: %s" % (newlines, timestamp, message))
        # We force an update in case the next step is a long one, like
        # running a simluation. 
        self.TextConsole.Update()



def hsinc(npts, ncycles, filter='hamming'):
    """
    Returns a sinc function of length npts, with ncycles sinc-cycles. This
    yields a time-bandwidth value of 4 * ncycles
    """    
    t = np.arange(npts) - (npts/2.0)
    t = t / (npts/2.0)
    val = 2*np.pi*ncycles*t + 0.00001
    res = np.sin(val) / val
    if filter == 'hamming':
        res = res * 4 * ncycles * (0.54 + 0.46*np.cos(np.pi*t)) / npts
        
    return res   


if __name__ == "__main__":
    app = wx.PySimpleApp(0)
    wx.InitAllImageHandlers()
    DialogEditor = DialogEditorTransform(None, "", "")
    #DialogEditor = DialogEditorTransform(None, wx.ID_ANY, "")
    app.SetTopWindow(DialogEditor)
    DialogEditor.Show()
    app.MainLoop()