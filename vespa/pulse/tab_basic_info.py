# Python modules


# 3rd party modules
import wx


# Our modules
import vespa.pulse.dialog_machine_specs as dialog_machine_specs
import vespa.pulse.auto_gui.basic_info as basic_info
import vespa.common.util.time_ as util_time
import vespa.common.util.misc as util_misc
import vespa.common.constants as constants
import vespa.common.wx_gravy.common_dialogs as common_dialogs
import vespa.common.rfp_machine_specs as rfp_machine_specs


class TabBasicInfo(basic_info.PanelBasicInfo):
    
    def __init__(self, parent, db, pulse_design):
        
        basic_info.PanelBasicInfo.__init__(self, parent)

        self._inner_notebook = parent
        self._db             = db
        self._pulse_design   = pulse_design
        self._last_run = { }
        self._last_save = { }
        self._current_focus_id = None
        
        self.machine_specs = self._pulse_design.machine_specs
        
        self.initialize_controls()

        # We set _last_run and _last_save unconditionally here, unlike in 
        # transform tabs where they're only set when opening an 
        # existing tab.
        # The difference means that on a new design, the basic info tab is 
        # not out of sync when it's first opened. Nor does it have unsaved 
        # changes, so it can be closed without any complaint from Pulse. 
        self._last_run = self.get_raw_gui_data()
        self._last_save = self.get_raw_gui_data()

        self.Layout()
        self.Fit()
        
        self.Bind(wx.EVT_CHILD_FOCUS, self.on_child_focus)
        self.Bind(wx.EVT_TEXT, self.on_text)


    @property
    def is_saved(self):
        """
        True if the current input is the same as it was when the tab
        was last saved.
        """
        return self._last_save == self.get_raw_gui_data()


    @property
    def is_synced(self):
        """True if the current input is in sync with the results."""
        return self._last_run == self.get_raw_gui_data()


    ##### Event Handlers ######################################################
    def on_text(self, event):
        #
        #                         Hack alert! 
        #
        # This event fires when text is changed in any textbox. If the change
        # happens in one of the controls listed below (name, creator, or 
        # comment), then we do two things. First, if the control's value is
        # present in _last_run then we update it immediately. Second, we 
        # unconditionally update the value of the field in the pulse_design
        # object.
        # We update _last_run because these three fields don't affect the  
        # design results, so we don't want a change in one of these fields 
        # to make the design out of sync. 
        # We update the pulse design because we want these fields to be 
        # save-able even when the design is frozen. Normally, values are 
        # moved from the GUI to the objects only when the user runs the 
        # design. Since frozen designs can't be run, the event that would 
        # normally transfer the name etc. from the GUI into the pulse design 
        # can never happen. This hack sidesteps that problem.
        #
        # So yes, this is a hack, but it's the least ugly option.
        control = event.GetEventObject()
        
        control_text = control.GetValue()
        
        # key and attribute_name tell me which key to look for in _last_run
        # and which attribute to update on the pulse design. At present 
        # they always happen to be the same but I thought it might be 
        # confusing to store them both in the same variable so each gets its
        # own variable.
        key = attribute_name = None
        
        if control.Id == self.TextName.Id:
            key = "name"
            attribute_name = "name"
        elif control.Id == self.TextCreator.Id:
            key = "creator"
            attribute_name = "creator"
        elif control.Id == self.TextComment.Id:
            key = "comment"
            attribute_name = "comment"

        if key:
            if key in self._last_run:
                self._last_run[key] = control_text
            setattr(self._pulse_design, attribute_name, control_text)
        #else:
            # This event happened in a textbox that doesn't participate in 
            # this hack.


    def on_activation(self):
        # This is a faux event handler. wx doesn't call it directly. It's 
        # a notification from my parent (the experiment notebook) to let
        # me know that this tab has become the current one.
        pass


    def on_child_focus(self, event):
        # When the focus changes we take the opportunity to note if the 
        # user has made any changes and, if so, we update the sync status
        # indicator. The design notebook actually does that work; this 
        # is just a trigger.
        
        # The wx doc for wxChildFocusEvent says --
        # "Notice that child window is the direct child of the window 
        # receiving event. Use FindFocus to retreive [sic] the window which  
        # is actually getting focus."
        recipient = wx.Window.FindFocus()

        # I'm not sure why, but we seem to get a surplus of these messages.
        # Here we ensure that the focus really is somewhere else.
        if self._current_focus_id != recipient.Id:
            # Focus changed
            self._current_focus_id = recipient.Id            
            self._inner_notebook.update_sync_status(self)


    def on_edit_machine_specs(self, event):
        dialog = dialog_machine_specs.DialogMachineSpecs( self, 
                                                          self._db,
                                                          self.machine_specs,
                                                          self._pulse_design.is_frozen)
        if dialog.ShowModal() == wx.ID_OK:
            self.machine_specs = dialog.machine_specs
            self.update_machine_specs_summary()
            # Update sync status since the machine specs probably just changed
            self._inner_notebook.update_sync_status(self)


    ##### Helper Functions ####################################################
    
    def accept_gui_data(self, metadata_only=False):
        """ 
        See documentation here:
        http://scion.duhs.duke.edu/vespa/rfpulse/wiki/CommonTabFeatures
        """
        # When metadata_only is True, only the fields that don't alter sync
        # status (currently name & comments) are copied from the GUI.
        d = self.get_cooked_gui_data()
        if metadata_only:
            # Discard everything except name & comment
            d = { "name" : d["name"], "comment" : d["comment"] }

        self._pulse_design.inflate(d)
        
        
    def initialize_controls(self):
        
        self.TextName.ChangeValue(self._pulse_design.name)
        self.StaticTextUUID.SetLabel(self._pulse_design.id)
        self.TextCreator.ChangeValue(self._pulse_design.creator)
        created = self._pulse_design.created.strftime(util_time.DISPLAY_DATE_FORMAT)
        self.StaticTextCreated.SetLabel(created)
        self.TextComment.ChangeValue(self._pulse_design.comment)            

        self.update_machine_specs_summary()
        
        self.TextCalcResolution.ChangeValue(str(self._pulse_design.calc_resolution))    
        
#        # Set values in PulseConvention ComboBox.
#        self.ComboBandwidthType.Clear()        
#        self.ComboBandwidthType.AppendItems( \
#                            [constants.PulseConvention.HALF_HEIGHT['display'], ])
                             # :FIXME: Implement these in SLR so it does not break.
                             # constants.PulseConvention.MINIMUM['display'],
                             # constants.PulseConvention.MAXIMUM['display']])
        
        if self._pulse_design.pulse_bandwidth_type == 'half_height':
            s = "FW at Half Height"
        elif self._pulse_design.pulse_bandwidth_type == 'minimum':
            s = "FW at Min Height"
        elif self._pulse_design.pulse_bandwidth_type == 'maximum':
            s = "FW at Max Height"
        else:
            s = "FW at Half Height"
        self.ChoiceBandwidthType.SetStringSelection(s)
        
        s = str(self._pulse_design.gyromagnetic_nuclei)
        self.ChoiceGyromagneticNuclei.SetStringSelection(s)

        if self._pulse_design.is_frozen:
            for control in (self.TextCreator, 
                            self.TextCalcResolution,
                            self.ChoiceBandwidthType,
                            self.ChoiceGyromagneticNuclei, ):
                control.Disable()
                
            self.ButtonEditMachineSpecs.SetLabel("View...")


    def get_cooked_gui_data(self):
        """ See documentation here:
        http://scion.duhs.duke.edu/vespa/rfpulse/wiki/CommonTabFeatures
        """
        d = self.get_raw_gui_data()
        
        d["machine_specs"]   =  rfp_machine_specs.MachineSpecs(d["machine_specs"])
        d["calc_resolution"] = int(d["calc_resolution"])
        
        if d["pulse_bandwidth_type"] == "FW at Half Height":
            d["pulse_bandwidth_type"] = 'half_height'
        elif d["pulse_bandwidth_type"] == "FW at Min Height":
            d["pulse_bandwidth_type"] = 'minimum'
        elif d["pulse_bandwidth_type"] == "FW at Max Height":
            d["pulse_bandwidth_type"] = 'maximum'
        else:
            d["pulse_bandwidth_type"] = 'half_height'
        
        # d["gyromagnetic_nuclei"] should not need cooking
        
        return d
        
            
    def get_raw_gui_data(self):
        """ See documentation here:
        http://scion.duhs.duke.edu/vespa/rfpulse/wiki/CommonTabFeatures
        """
        d = { }
        
        d["name"] = self.TextName.GetValue().strip()
        d["creator"] = self.TextCreator.GetValue().strip()        
        # Note that we don't call .strip() on the comment. It's totally
        # freeform and whitespace might be important.
        d["comment"] = self.TextComment.GetValue()
        d["calc_resolution"] = self.TextCalcResolution.GetValue().strip()
        #d["pulse_bandwidth_type"] = self.ComboBandwidthType.GetValue().strip()
        d["pulse_bandwidth_type"] = self.ChoiceBandwidthType.GetStringSelection().strip()
        #d["gyromagnetic_nuclei"] = self.ComboGyromagneticNuclei.GetValue().strip()
        d["gyromagnetic_nuclei"] = self.ChoiceGyromagneticNuclei.GetStringSelection().strip()
        d["machine_specs"] = self.machine_specs.deflate(constants.Deflate.DICTIONARY)

        return d


    def on_save(self):
        self._last_save = self.get_raw_gui_data()
        

    def run(self):
        """ See documentation here:
        http://scion.duhs.duke.edu/vespa/rfpulse/wiki/CommonTabFeatures
        """
        self._last_run = self.get_raw_gui_data()
        
        return True
        
    
    def update_uuid(self):
        """A small hack that's called after a 'save as' in order to update
        the displayed UUID (which changes during a 'save as').
        """
        self.StaticTextUUID.SetLabel(self._pulse_design.id)
        

    def update_machine_specs_summary(self):
        # Refreshes the summary info in the machine specs group box.
        summary = ""
        if self.machine_specs.machine_type in constants.MachineType.ALL:
            summary = self.machine_specs.machine_type["display"]
        else:
            summary = self.machine_specs.machine_type
            
        summary += (" %.1fT" % self.machine_specs.field_strength)
        
        self.LabelMachineSpecsSummary.SetLabel(summary)
        
        # Repaint to account for changes to the size of the summary
        self.LabelMachineSpecsSummary.GetContainingSizer().Layout()


    def validate_gui(self, validate_for_save=False):
        """ See documentation here:
        http://scion.duhs.duke.edu/vespa/rfpulse/wiki/CommonTabFeatures
        
        validate_gui() implicitly validates the GUI to see if it can be
        run. This validate_gui() has an extra feature which is that is can 
        validate the GUI for run or for save. On most tabs, they're the same.
        On this tab, the name (which doesn't affect the run) can be invalid
        for save (e.g. if it's blank).
        """
        d = self.get_raw_gui_data()
        
        msg = ""
        
        if validate_for_save:
            if not d["name"]:
                msg = """Please enter a name for the design."""
        else:
            # validate for run
            if not util_misc.is_intable(d["calc_resolution"]):
                msg = """I don't understand the calculation resolution "%s".""" % d["calc_resolution"]

        if msg:
            self._inner_notebook.activate_tab(self)
            common_dialogs.message(msg)
            
        return not bool(msg)
