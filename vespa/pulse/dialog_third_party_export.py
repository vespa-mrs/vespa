# Python modules

import os

# 3rd party modules
import wx
import numpy as np

# Our modules
import vespa.pulse.auto_gui.third_party_export as third_party_export
import vespa.common.pulse_funcs.siemens_export as siemens_export
import vespa.common.util.misc as util_misc
import vespa.common.util.time_ as util_time
import vespa.common.util.config as util_config
import vespa.common.constants as constants
import vespa.common.wx_gravy.common_dialogs as common_dialogs
import vespa.common.wx_gravy.util as wx_util



IDEA_MINSLICE   =  1.0
IDEA_MAXSLICE   = 40.0
VISION_MINSLICE =  1.0
VISION_MAXSLICE = 40.0
DEFAULT_FAMILY  = 'vespa'


def _write_lines(filename, lines):
    """
    Given a filename and a list of lines, joins those lines into a
    newline-delimited string and writes it to filename. If the file 
    exists, it's truncated before writing.
    
    Newlines are always written as Windows newlines (CRLF = 0xd 0xa) because
    that's what Siemens machines expect.
    
    The strings in the list of lines are converted to Unicode before
    writing if they're not Unicode already. If any of the non-Unicode
    strings contain non-ASCII, the conversion will fail. For details,
    see the following:
    http://scion.duhs.duke.edu/vespa/project/wiki/ThePerilsOfStr
    """
    # Note that all_lines will be a Unicode string after this join() 
    # because anything returned from wx will be Unicode, and any 
    # non-Unicode string joined to a Unicode string is "promoted" to
    # Unicode. Atfer the join we explicitly force the Unicode string
    # to UTF-8 so that it will be safe for write().
    lines = "\n".join(lines)
    lines = lines.replace("\n", "\r\n")
    lines = lines.encode("utf-8")
    

    try:
        open(filename, "wb").write(lines)
    except (OSError, IOError) as error_instance:
        msg = """Writing the file "%s" failed. The operating system reports:\n\n"""
        msg %= filename
        msg += error_instance.strerror
        common_dialogs.message(msg, "Pulse - Third Party Export")

#------------------------------------------------------------------------------
# Note. GUI Architecture/Style
#
# Vespa GUI components are designed with WxGlade to speed up development times.
# Only stub functions for event methods are created. The WxGlade files (with
# *.wxg extensions) are stored in the 'wxglade' subdirectory. The ouput of their
# code generation are stored in the 'auto_gui' subdirectory.
#
# Each auto_gui class is inherited into a unique 'vespa' class. Program specific
# initialization and widget event handlers are overloaded to provide program
# specific event handling.
#------------------------------------------------------------------------------


class DialogThirdPartyExport(third_party_export.MyDialog):

    def __init__(self, parent, pulse_design_id, pulse, bandwidth, tip_angle):
        
        if not parent:
            parent = wx.GetApp().GetTopWindow()

        third_party_export.MyDialog.__init__(self, parent)
        
        self.pulse_design_id = pulse_design_id
        self.pulse = pulse
        self.bandwidth = bandwidth
        self.tip_angle = tip_angle
        
        self.format = constants.ThirdPartyExportFormat.IDEA
        
        self.initialize_controls()
        self.populate_controls()

        self.Fit()
        self.Center()



    ##### Event Handlers ######################################################
    
    def on_browse(self, event): 

        default_path = util_config.get_last_export_path()

        fname = self.LabelFilename.GetLabel().strip()
        if fname:
            default_path, fname = os.path.split(fname)
            fname, ext = os.path.splitext(fname)
        else:
            if self.format == constants.ThirdPartyExportFormat.IDEA:
                fname = self.TextIdeaName.GetValue().strip()
                if not fname:
                    fname = 'idea_pulse_export'
            elif self.format == constants.ThirdPartyExportFormat.VISION:
                fname = self.TextVisionName.GetValue().strip()
                if not fname:
                    fname = 'vision_pulse_export'
            elif self.format == constants.ThirdPartyExportFormat.ASCII_MAGN_PHASE:
                fname = 'ascii_magn_phase_pulse_export'
            elif self.format == constants.ThirdPartyExportFormat.ANNOTATED_ASCII_MAGN_PHASE:
                fname = 'annotated_ascii_magn_phase_pulse_export'
            elif self.format == constants.ThirdPartyExportFormat.IDEA_C_HEADER:
                fname = 'idea_pulse_c_header'

        if self.format == constants.ThirdPartyExportFormat.IDEA: 
            fname = fname +'.pta'
            filetype_filter = "Pulse Third Party Export Filename (*.pta)|*.pta"
        elif self.format == constants.ThirdPartyExportFormat.VISION:
            fname = fname +'.pta'
            filetype_filter = "Pulse Third Party Export Filename (*.pta)|*.pta"
        elif self.format == constants.ThirdPartyExportFormat.ASCII_MAGN_PHASE:
            fname = fname +'.txt'
            filetype_filter = "Pulse Third Party Export Filename (*.txt)|*.txt"
        elif self.format == constants.ThirdPartyExportFormat.ANNOTATED_ASCII_MAGN_PHASE:
            fname = fname +'.txt'
            filetype_filter = "Pulse Third Party Export Filename (*.txt)|*.txt"
        if self.format == constants.ThirdPartyExportFormat.IDEA_C_HEADER: 
            fname = fname +'.h'
            filetype_filter = "Pulse Third Party Export Filename (*.h)|*.h"

        
        fname = common_dialogs.save_as(filetype_filter=filetype_filter,
                                       default_path=default_path,
                                       default_filename=fname)
        if fname:
            self.LabelFilename.SetLabel(fname)
            path, _ = os.path.split(fname)
            util_config.set_last_export_path(path)
        
    
    def on_format(self, event):
        index = event.GetEventObject().GetStringSelection().strip()
        self.Freeze()
        if index == 'Siemens-IDEA':
            self.format = constants.ThirdPartyExportFormat.IDEA
            self.PanelSiemensIdea.Show()
            self.PanelSiemensVision.Hide()
            self.PanelAsciiMagnPhase.Hide()
            self.PanelSiemensIdeaCHeader.Hide()
            
        elif index == 'Siemens-Vision':
            self.format = constants.ThirdPartyExportFormat.VISION
            self.PanelSiemensIdea.Hide()
            self.PanelSiemensVision.Show()
            self.PanelAsciiMagnPhase.Hide()
            self.PanelSiemensIdeaCHeader.Hide()

        elif index == 'ASCII - Magn/Phase':
            self.format = constants.ThirdPartyExportFormat.ASCII_MAGN_PHASE
            self.PanelSiemensIdea.Hide()
            self.PanelSiemensVision.Hide()
            self.PanelAsciiMagnPhase.Show()
            self.PanelSiemensIdeaCHeader.Hide()
            
        elif index == 'Annotated ASCII - Magn/Phase':
            self.format = constants.ThirdPartyExportFormat.ANNOTATED_ASCII_MAGN_PHASE
            self.PanelSiemensIdea.Hide()
            self.PanelSiemensVision.Hide()
            self.PanelAsciiMagnPhase.Show()
            self.PanelSiemensIdeaCHeader.Hide()

        if index == 'Siemens-IDEA C Header':
            self.format = constants.ThirdPartyExportFormat.IDEA_C_HEADER
            self.PanelSiemensIdea.Hide()
            self.PanelSiemensVision.Hide()
            self.PanelAsciiMagnPhase.Hide()
            self.PanelSiemensIdeaCHeader.Show()
        
        fname = self.LabelFilename.GetLabel().strip()
        if fname:
            default_path, fname = os.path.splitext(fname)
            if self.format == constants.ThirdPartyExportFormat.IDEA: 
                self.LabelFilename.SetLabel(default_path + '.pta')
            elif self.format == constants.ThirdPartyExportFormat.VISION:
                self.LabelFilename.SetLabel(default_path + '.pta')
            elif self.format == constants.ThirdPartyExportFormat.ASCII_MAGN_PHASE:
                self.LabelFilename.SetLabel(default_path + '.txt')
            elif self.format == constants.ThirdPartyExportFormat.ANNOTATED_ASCII_MAGN_PHASE:
                self.LabelFilename.SetLabel(default_path + '.txt')
            elif self.format == constants.ThirdPartyExportFormat.IDEA_C_HEADER:
                self.LabelFilename.SetLabel(default_path + '.h')

        self.Layout()
        self.Fit()        
        self.Thaw()


    def on_ok(self, event):
        if self.validate_gui():
            
            d = self.get_cooked_gui_data()
            
            if self.format == constants.ThirdPartyExportFormat.IDEA:
                self.do_idea_export(d)
            elif self.format == constants.ThirdPartyExportFormat.VISION:
                self.do_vision_export(d)
            elif self.format == constants.ThirdPartyExportFormat.ASCII_MAGN_PHASE:
                self.do_ascii_magn_phase_export()
            elif self.format == constants.ThirdPartyExportFormat.ANNOTATED_ASCII_MAGN_PHASE:
                self.do_ascii_magn_phase_export(True)
            elif self.format == constants.ThirdPartyExportFormat.IDEA_C_HEADER:
                self.do_idea_c_header_export(d)
            
            # we were successful, close dialog
            self.Close()        


    ##### Internal helper functions  ##########################################

    def initialize_controls(self):
        # performs run-time changes to dialog widgets, but does not
        # put default values into them.
        # We add the OK & Cancel buttons dynamically so that they're in the 
        # right order under OS X, GTK, Windows, etc.
        self.ButtonOk, self.ButtonCancel = \
                    wx_util.add_ok_cancel(self, self.LabelOkCancelPlaceholder,
                                          self.on_ok)


    def populate_controls(self):
        # populates controls with default values using data from the pulse
        txt_bandwidth = str(self.bandwidth)
        txt_tip_angle = str(self.tip_angle)
        
        comment = 'Vespa-Pulse export: Design Name = %s '      \
                  ', UUID = %s, Duration = %s [ms] '            \
                  ', Bandwidth = %s [kHz], and Tip Angle = %s [deg]'
        comment %= (self.pulse.name, self.pulse_design_id, self.pulse.duration, self.bandwidth, self.tip_angle)

        # Note. When Siemens stores an RF pulse in an ASCII file, 
        # the file name is created from the family and pulse name 
        # as default (e.g., family.pulse). For this reason, do not 
        # use periods in the family or pulse names.
        export_name = self.pulse.name.replace(' ', '_')
        export_name = self.pulse.name.replace('.', '_')

        default_path = util_config.get_last_export_path()
        filename = os.path.join(default_path, export_name + ".pta")
        
        self.LabelFilename.SetLabel(filename)

        self.TextIdeaName.ChangeValue(export_name)
        self.TextIdeaComment.SetValue(comment)
        self.TextIdeaMinslice.ChangeValue(str(IDEA_MINSLICE))
        self.TextIdeaMaxslice.ChangeValue(str(IDEA_MAXSLICE))
        if self.tip_angle:
            self.TextIdeaFlipAngle.ChangeValue(str(self.tip_angle))

        self.TextVisionName.ChangeValue(export_name)
        self.TextVisionComment.SetValue(comment)
        self.TextVisionFamilyName.ChangeValue(DEFAULT_FAMILY)
        self.TextVisionMinslice.ChangeValue(str(IDEA_MINSLICE))
        self.TextVisionMaxslice.ChangeValue(str(IDEA_MAXSLICE))
        if self.tip_angle:
            self.TextVisionFlipAngle.ChangeValue(str(self.tip_angle))

        self.TextHeaderVariablesBaseString.ChangeValue("MyCustomRf")
        self.TextHeaderComment.SetValue(comment)
        self.TextHeaderMinslice.ChangeValue(str(IDEA_MINSLICE))
        if self.tip_angle:
            self.TextHeaderFlipAngle.ChangeValue(str(self.tip_angle))

        self.PanelSiemensIdea.Show()
        self.PanelSiemensVision.Hide()
        self.PanelAsciiMagnPhase.Hide()
        self.PanelSiemensIdeaCHeader.Hide()


    def validate_gui(self):
        # Validate controls one by one
        # A non-empty string message indicates failure
        
        d = self.get_raw_gui_data()
        
        msg = ""

        if not d['filename']:
            msg = "Please enter an output location and name for the export file."

        if self.format == constants.ThirdPartyExportFormat.IDEA:

            if not msg:
                if not d["idea_name"]:
                    msg = "Please enter a name for this pulse (do not use '.' or spaces)."
                    
            if not msg:
                if not util_misc.is_floatable(d["idea_minslice"]):
                    msg = """I don't understand the min slice thickness "%s".""" % d["idea_minslice"]

            if not msg:
                if not util_misc.is_floatable(d["idea_maxslice"]):
                    msg = """I don't understand the max slice thickness "%s".""" % d["idea_maxslice"]

            if not msg:
                if not util_misc.is_floatable(d["idea_bandwidth"]):
                    msg = """I don't understand the bandwidth "%s".""" % d["idea_bandwidth"]

            if not msg:
                if not util_misc.is_floatable(d["idea_flip_angle"]):
                    msg = """I don't understand the flip angle "%s".""" % d["idea_flip_angle"]

        elif self.format == constants.ThirdPartyExportFormat.VISION:

            if not msg:
                if not d["vision_name"]:
                    msg = """Please enter a pulse name (do not use '.' or spaces)."""

            if not msg:
                if not d["vision_family_name"]:
                    msg = """Please enter a family name (do not use '.' or spaces). """

            if not msg:
                if not util_misc.is_floatable(d["vision_minslice"]):
                    msg = """I don't understand the min slice thickness "%s".""" % d["vision_minslice"]

            if not msg:
                if not util_misc.is_floatable(d["vision_maxslice"]):
                    msg = """I don't understand the max slice thickness "%s".""" % d["vision_maxslice"]

            if not msg:
                if not util_misc.is_floatable(d["vision_bandwidth"]):
                    msg = """I don't understand the bandwidth "%s".""" % d["vision_bandwidth"]

            if not msg:
                if not util_misc.is_floatable(d["vision_flip_angle"]):
                    msg = """I don't understand the flip angle "%s".""" % d["vision_flip_angle"]

        elif self.format == constants.ThirdPartyExportFormat.IDEA_C_HEADER:

            if not msg:
                if not d["header_base_string"]:
                    msg = "Please enter a variable name base string for this pulse (do not use '.' or spaces)."
                    
            if not msg:
                if not util_misc.is_floatable(d["header_minslice"]):
                    msg = """I don't understand the min slice thickness "%s".""" % d["header_minslice"]

            if not msg:
                if not util_misc.is_floatable(d["header_bandwidth"]):
                    msg = """I don't understand the bandwidth "%s".""" % d["header_bandwidth"]

            if not msg:
                if not util_misc.is_floatable(d["header_flip_angle"]):
                    msg = """I don't understand the flip angle "%s".""" % d["header_flip_angle"]

        if not msg:
            # At this point we know the raw data can be cooked, and all of
            # remaining validation is on the cooked data.
            d = self.get_cooked_gui_data()

        if not msg:
            if self.format == constants.ThirdPartyExportFormat.IDEA:
                if d["idea_bandwidth"] <= 0:
                    msg = "Please enter a bandwidth > 0."
                if d["idea_flip_angle"] <= 0:
                    msg = "Please enter a flip angle >= 0."
            elif self.format == constants.ThirdPartyExportFormat.VISION:
                if d["vision_bandwidth"] <= 0:
                    msg = "Please enter a bandwidth > 0."
                if d["vision_flip_angle"] <= 0:
                    msg = "Please enter a flip angle >= 0."
            elif self.format == constants.ThirdPartyExportFormat.IDEA_C_HEADER:
                if d["header_bandwidth"] <= 0:
                    msg = "Please enter a bandwidth > 0."
                if d["header_flip_angle"] <= 0:
                    msg = "Please enter a flip angle >= 0."

        if msg:
            # validation failed
            common_dialogs.message(msg, "Third Party Export", common_dialogs.I_OK)
            
                                
        return not bool(msg)        


    def get_cooked_gui_data(self):
        """ 
        See documentation here:
        http://scion.duhs.duke.edu/vespa/rfpulse/wiki/CommonTabFeatures
        """
        d = self.get_raw_gui_data()
        
        if self.format == constants.ThirdPartyExportFormat.IDEA:
            for key in ("idea_minslice", "idea_maxslice",
                        "idea_bandwidth", "idea_flip_angle"):
                d[key] = float(d[key])
        elif self.format == constants.ThirdPartyExportFormat.VISION:
            for key in ("vision_minslice", "vision_maxslice", 
                        "vision_bandwidth", "vision_flip_angle"):
                d[key] = float(d[key])
        elif self.format == constants.ThirdPartyExportFormat.IDEA_C_HEADER:
            for key in ("header_minslice", "header_bandwidth", "header_flip_angle"):
                d[key] = float(d[key])
        
        return d


    def get_raw_gui_data(self):
        """ 
        See documentation here:
        http://scion.duhs.duke.edu/vespa/rfpulse/wiki/CommonTabFeatures
        """
        d = { }

        d["filename"] = self.LabelFilename.GetLabel().strip()

        if self.format == constants.ThirdPartyExportFormat.IDEA:
            d["idea_name"] = self.TextIdeaName.GetValue().strip()
            d["idea_comment"] = self.TextIdeaComment.GetValue()
            d["idea_minslice"] = self.TextIdeaMinslice.GetValue().strip()
            d["idea_maxslice"] = self.TextIdeaMaxslice.GetValue().strip()
            d["idea_bandwidth"] = self.TextIdeaBandwidth.GetValue().strip()
            d["idea_flip_angle"] = self.TextIdeaFlipAngle.GetValue().strip()
        
        elif self.format == constants.ThirdPartyExportFormat.VISION:
            d["vision_name"] = self.TextVisionName.GetValue().strip()
            d["vision_comment"] = self.TextVisionComment.GetValue().strip()
            d["vision_family_name"] = self.TextVisionFamilyName.GetValue().strip()
            d["vision_minslice"] = self.TextVisionMinslice.GetValue().strip()
            d["vision_maxslice"] = self.TextVisionMaxslice.GetValue().strip()
            d["vision_bandwidth"] = self.TextVisionBandwidth.GetValue().strip()
            d["vision_flip_angle"] = self.TextVisionFlipAngle.GetValue().strip()

        elif self.format == constants.ThirdPartyExportFormat.IDEA_C_HEADER:
            d["header_base_string"] = self.TextHeaderVariablesBaseString.GetValue().strip()
            d["header_comment"] = self.TextHeaderComment.GetValue()
            d["header_minslice"] = self.TextHeaderMinslice.GetValue().strip()
            d["header_bandwidth"] = self.TextHeaderBandwidth.GetValue().strip()
            d["header_flip_angle"] = self.TextHeaderFlipAngle.GetValue().strip()
                    
        return d

        
    def do_idea_export(self, d):
        waveform  = np.array(self.pulse.rf_waveform)
        
        all_lines = siemens_export.siemens_export( waveform,
                                        constants.ThirdPartyExportFormat.FIELD_B1,
                                        self.pulse.dwell_time,
                                        d["idea_bandwidth"],
                                        d["idea_maxslice"],
                                        d["idea_minslice"],
                                        d["idea_name"],
                                        d["idea_comment"], 
                                        constants.ThirdPartyExportFormat.IDEA,
                                        d["idea_flip_angle"])

        filename = self.LabelFilename.GetLabel()
        _write_lines(filename, all_lines)
        

    def do_vision_export(self, d):
        waveform  = np.array(self.pulse.rf_waveform)

        all_lines = siemens_export.siemens_export(waveform,
                                    constants.ThirdPartyExportFormat.FIELD_B1,
                                    self.pulse.dwell_time,
                                    d["vision_bandwidth"],
                                    d["vision_maxslice"],
                                    d["vision_minslice"],
                                    d["vision_name"],
                                    d["vision_comment"],
                                    constants.ThirdPartyExportFormat.VISION,
                                    d["vision_family_name"],
                                    d["vision_flip_angle"])
                              
        filename = self.LabelFilename.GetLabel()
        _write_lines(filename, all_lines)


    def do_idea_c_header_export(self, d):
        waveform  = np.array(self.pulse.rf_waveform)
        filename = self.LabelFilename.GetLabel()
        fbase = os.path.basename(filename)
        
        all_lines = siemens_export.siemens_export( waveform,
                                        constants.ThirdPartyExportFormat.FIELD_B1,
                                        self.pulse.dwell_time,
                                        d["header_bandwidth"],
                                        200.0,
                                        d["header_minslice"],
                                        '',
                                        d["header_comment"], 
                                        constants.ThirdPartyExportFormat.IDEA_C_HEADER,
                                        '',
                                        d["header_flip_angle"],
                                        d["header_base_string"],
                                        fbase)

        _write_lines(filename, all_lines)

    def do_ascii_magn_phase_export(self, include_annotation=False):
        """
        Outputs rf waveform (and possibly a gradient) to a textfile in columnar
        format,  Amplitude Phase(deg) Gradient, one array entry per line. 
        
        Gradient values are added to a third column only if there is a gradient
        present in the rf_result object, and it is the same length as the 
        rf_waveform array. No error is reported if a gradient is present but 
        not the same length, it is just not included.
        
        """
        all_lines = [ ]
        if include_annotation:
            all_lines.append('# Exported from Vespa-Pulse')
            all_lines.append('# For more information, visit http://scion.duhs.duke.edu/vespa/')
            all_lines.append('# Export Timestamp: ' + util_time.now(util_time.ISO_TIMESTAMP_FORMAT))
            all_lines.append('# Pulse Name: ' + self.pulse.name)
            all_lines.append('# Design UUID: ' + self.pulse_design_id)
            all_lines.append('# Pulse Total Duration: %s [ms]'  % self.pulse.duration)
            all_lines.append('# ')
            all_lines.append('# Note. If a gradient waveform is part of the pulse design, and it is ')
            all_lines.append('#  the same length as the rf_waveform, it will be included as a third')
            all_lines.append('#  column of this output. No error is reported if a different length.')
            all_lines.append('# ----------------------------------------')
        #else:
            # Annotated ASCII is the same as un-annotated ASCII except
            # for the comment lines we added just above.
            
        # Get waveform from results
        waveform = np.array(self.pulse.rf_waveform) 
        
        gflag = False
        if self.pulse.gradient is not None:
            grad_arr = np.array(self.pulse.gradient)         
            gflag = (np.sum(grad_arr) != 0.0)
            if gflag:
                if len(grad_arr) != len(waveform):
                    gflag = False
        
        ampl_arr = np.abs(waveform)
        # Simulation expects angle in degrees so use that as default
        phas_arr = constants.RADIANS_TO_DEGREES * np.angle(waveform)

        # using same float format as Siemens IDEA 
        if gflag:
            format = "%1.9f %1.9f %1.9f" 
            for i in range(len(ampl_arr)):
                line = format % (ampl_arr[i], phas_arr[i], grad_arr[i])
                all_lines.append(line)
        else:
            format = "%1.9f %1.9f" 
            for i in range(len(ampl_arr)):
                line = format % (ampl_arr[i], phas_arr[i])
                all_lines.append(line)

        filename = self.LabelFilename.GetLabel()
        _write_lines(filename, all_lines)

