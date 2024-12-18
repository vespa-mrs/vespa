# Python modules

import os
import math
import tempfile

# 3rd party modules
import wx
#import wx.aui as aui
import wx.lib.agw.aui as aui        # NB. wx.aui version throws odd wxWidgets exception on Close/Exit  ?? Not anymore in wxPython 4.0.6 ??

import wx.stc as wx_stc
import numpy as np

# Our modules
import vespa.simulation.prefs as prefs
import vespa.simulation.plot_panel_editor as plot_panel_editor
import vespa.simulation.util_simulation_config as sim_config
import vespa.simulation.run_experiment_controller as run_experiment_controller
import vespa.simulation.dialog_experiment_list as dialog_experiment_list
import vespa.simulation.mrs_data_basis as mrs_data_basis
import vespa.simulation.build_basis_functions as bbf_module
import vespa.simulation.auto_gui.pulse_sequence_editor as pulse_sequence_editor
import vespa.common.styled_text_control as our_stc
import vespa.common.mrs_pulse_sequence as mrs_pulse_sequence
import vespa.common.mrs_simulation as mrs_simulation
import vespa.common.mrs_experiment as mrs_experiment
import vespa.common.constants as common_constants
import vespa.common.wx_gravy.common_dialogs as common_dialogs
import vespa.common.wx_gravy.util as common_wx_util
import vespa.common.util.ppm as util_ppm
import vespa.common.util.misc as util_misc
import vespa.common.util.time_ as util_time
import vespa.common.util.generic_spectral as util_generic_spectral
import vespa.common.util.misc as util_misc

#import vespa.common.dialog_pulse_project_browser as dialog_pulse_project_browser
import vespa.common.dialog_pulse_design_browser as dialog_pulse_design_browser

PI = math.pi

# This is the row of lines that goes in the console to separate the messages
# of one run from another.
SEPARATOR = "-" * 60

class LoopControlSet(object):
    """A container for a set of loop controls. They are the column label
    and the three textboxes for start, count and size"""
    def __init__(self, label, start):
        self.label = label
        self.start = start

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


class DialogPulseSequenceEditor(pulse_sequence_editor.MyDialog):

    def __init__(self, parent, db, pulse_sequence=None):

        config = sim_config.Config()
        position, size = config.get_window_coordinates("dialog_pulse_sequence_editor")
        self._left, self._top = position
        self._width, self._height = size

        if not parent:
            parent = wx.GetApp().GetTopWindow()

        pulse_sequence_editor.MyDialog.__init__(self, parent)

        self.parent         = parent
        self.db             = db
        self.pulse_sequence = pulse_sequence
        self.experiment     = mrs_experiment.Experiment()

        if self.pulse_sequence:

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
            self.pulse_sequence = mrs_pulse_sequence.PulseSequence()
            self.pulse_sequence.id = util_misc.uuid()
            self.pulse_sequence.experiment_names = []
            title = "New Pulse Sequence"

        self.SetTitle(title)

        #------------------------------
        # Setup plot parameters

        # plot attributes

        self.foreground = "black"
        self.background = "white"
        self._prefs = prefs.PrefsMain()

        # metabolite model attributes

        d = { "frequency" : self.experiment.b0 }
        self.basis = mrs_data_basis.DataBasis(d)
        self.freq1d = mrs_data_basis.DataBasis(d)

        self.maxppm         = self.freq1d.pts2ppm(0)
        self.minppm         = self.freq1d.pts2ppm(self.freq1d.dims[0]-1)
        self.linewidth      = common_constants.DEFAULT_LINEWIDTH
        self.apodization    = None
        self.isotope        = '1H'
        self.phase0         = 0.0
        self.phase1         = 0.0
        self.pivot          = 4.7   # ppm

        #------------------------------
        # Initialize widget controls

        self.initialize_controls()

        self.labelStatus.SetLabel('')

        # Set dialog size & position from the INI file.
        self.SetSize( (self._width, self._height) )
        self.SetPosition( (self._left, self._top) )

        # The sash (a.k.a. grabber) needs to be 550 wide to show the entire
        # left-hand pane under OS X. Under Windows & Ubuntu, it needs less.
        sash_position = 550

        self.SplitHorizontal.SetSashPosition(sash_position, True)

        # self.Center()

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

    def on_add_static_parameter(self, event=None):
        # This event is called both from wx as a proper event and from
        # our code in which case event will be None. Don't rely on it being
        # present!

        # sizer is a GridSizer, thus we add a row to it and fill in the
        # columns with widgets to describe a user parameter

        sizer = self.DesignParameterSizer
        sizer.SetRows(len(self.design_parameter_controls) + 1)

        self.add_parameter_row()
        row = self.design_parameter_controls[-1]

        checkbox, combobox, label_name, text_name, label_value, text_value = row
        sizer.Add(checkbox,    0, wx.ALIGN_CENTER_VERTICAL|wx.RIGHT, 5)
        sizer.Add(combobox,    0, wx.RIGHT, 5)
        sizer.Add(label_name,  0, wx.LEFT|wx.ALIGN_CENTER_VERTICAL, 5)
        sizer.Add(text_name,   1, wx.LEFT|wx.RIGHT|wx.EXPAND, 5)
        sizer.Add(label_value, 0, wx.LEFT|wx.ALIGN_CENTER_VERTICAL, 5)
        sizer.Add(text_value,  1, wx.LEFT|wx.RIGHT|wx.EXPAND, 5)

        self.PanelDesign.Layout()

        self.correct_tab_order()


    def on_remove_static_parameter(self, event):
        # Remove all checked rows from the user static param definition
        # controls.
        # I loop over the rows of controls but I can't delete them right
        # away because it's not kosher to alter a list while one is
        # iterating over it.
        delete_these_rows = [ ]
        for i, row in enumerate(self.design_parameter_controls):
            checkbox, _, _, _, _, _ = row
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
            self.DesignParameterSizer.SetRows(len(self.design_parameter_controls))
            self.PanelDesign.Layout()
        else:
            common_dialogs.message("Please check the parameter you want to remove.", style=common_dialogs.I_OK)


    def on_add_pulse(self, event):
        dialog = dialog_pulse_design_browser.DialogPulseDesignBrowser(self, self.db)
        dialog.ShowModal()
        pulse_design_id = dialog.selected_pulse_design_id
        if pulse_design_id:
            msg = ""

            # We kept PulseSequence attribute named pulse_projects even though
            # Pulse pulse_designs are replacing RFPulse pulse_projects

            ids = [item.id for item in self.pulse_sequence.pulse_projects]
            if pulse_design_id in ids:
                msg = "This pulse is already available in this pulse sequence."

            if not msg:
                pulse_design = self.db.fetch_pulse_design(pulse_design_id)
                if not pulse_design.get_pulse():
                    msg = """The pulse design "%s" doesn't have a final """ \
                          """pulse. Simulation can only use completed """    \
                          """pulse designs.""" % pulse_design.name

            if msg:
                common_dialogs.message(msg)
            else:
                # All is well; add the pulse project

                self.pulse_sequence.pulse_projects.append(pulse_design)

                checkbox, label = self.add_pulse_row(pulse_design)

                sizer = self.DesignPulseSizer
                # sizer is a GridSizer, thus we add a row to it and fill in the
                # columns with widgets to describe a user pulses
                sizer.SetRows(len(self.design_pulse_controls) + 1)
                sizer.Add(checkbox, 0, wx.ALIGN_CENTER_VERTICAL|wx.RIGHT, 5)
                sizer.Add(label,    1, wx.LEFT|wx.ALIGN_CENTER_VERTICAL, 5)

                self.PanelYourRFPulses.Layout()

                self.correct_tab_order()


    def on_remove_pulse(self, event):
        # Remove all checked rows from the user static param definition
        # controls.
        # I loop over the rows of controls but I can't delete them right
        # away because it's not kosher to alter a list while one is
        # iterating over it.
        delete_these_rows = [ ]
        for i, row in enumerate(self.design_pulse_controls):
            checkbox, _ = row
            if checkbox.IsChecked():
                delete_these_rows.append(i)

        # Now I can clean up the list. I do it from the highest index to the
        # lowest. If I did it lowest to highest, then deleting the first
        # item would throw off the indices of the subsequent items.
        delete_these_rows.reverse()

        for row in delete_these_rows:
            for control in self.design_pulse_controls[row]:
                # Been nice knowin' ya!
                control.Destroy()
            del self.design_pulse_controls[row]
            del self.pulse_sequence.pulse_projects[row]

        if delete_these_rows:
            # Refresh the grid.
            self.DesignPulseSizer.SetRows(len(self.design_pulse_controls))
            self.PanelYourRFPulses.Layout()
        else:
            common_dialogs.message("Please check the pulses you want to remove.", style=common_dialogs.I_OK)


    def on_page_changed(self, event):

        index = self.NotebookTools.GetSelection()

        # we check the designer to see if the current set of parameters
        # pass the validation step, and if so, we update the tester widget

        if index == 1:
            if self.validate_pulse_sequence():
                self.initialize_test_controls()
            else:
                self.NotebookTools.SetSelection(0)


    def on_combo_isotope(self, event):
        # update the metabolite comboBox control when an isotope is selected
        val = event.GetEventObject().GetStringSelection()
        if val == "1H":
            # DEFAULT_PROTON_CENTER_PPM should be 4.7
            ppm = common_constants.DEFAULT_PROTON_CENTER_PPM
        else:
            # DEFAULT_XNUCLEI_CENTER_PPM should be 0.0
            ppm = common_constants.DEFAULT_XNUCLEI_CENTER_PPM

        self.freq1d.ctrppm     = ppm
        self.basis.ctrppm      = ppm

        self.populate_metabolites()


    def on_linewidth(self, event):
        # re-calculate the self.apodization array for use in the
        # plot_canvas routine
        self.linewidth = event.GetEventObject().GetValue()
        self.calc_apodization()
        self.plot_canvas()


    def calc_apodization(self):
        xx = np.arange(self.freq1d.dims[0], dtype='float64') / self.freq1d.sw
        type_ = 'Gaussian' if self.CheckGaussian.IsChecked() else 'Lorentzian'
        self.apodization = util_generic_spectral.apodize(xx, self.linewidth, type_)


    def on_spectral_points(self, event):
        # make sure that the points selected are a power of 2 so that the
        # fft routines are happy.  Current min/max are 256 and 32768.
        # Also, basis functions need to be re-calculated here because their
        # array size has changed
        val = event.GetEventObject().GetValue()
        pow2 = 0
        while val > 2**pow2: pow2 = pow2 + 1
        val = 2**pow2
        event.GetEventObject().SetValue(float(val))

        shape = list(self.freq1d.data.shape)
        shape[-1] = val
        self.freq1d.data = np.empty(shape, dtype='complex64')

        shape = list(self.basis.data.shape)
        shape[-1] = val
        self.basis.data = np.empty(shape, dtype='complex64')

        self.build_basis_functions()
        self.calc_apodization()
        self.plot_canvas()


    def on_sweep_width(self, event):
        # Basis functions need to be re-calculated here because their
        # dwell time has changed.
        # Also, sweep width changes the extent of the x-axis in the plot,
        # thus we need to recalculate the max extent of the x-axis and
        # reset_history.  If possible we also try to save the current
        # x-axis window.
        #
        # FIXME - this double plot causes a 'flicker', it would be nice to
        #         somehow just reset the history[0] xlim value rather than
        #         double plotting ...
        val = event.GetEventObject().GetValue()
        self.freq1d.sw     = val
        self.basis.sw      = val
        self.maxppm = self.freq1d.pts2ppm(0)
        self.minppm = self.freq1d.pts2ppm(self.freq1d.dims[0]-1)
        xlim = self.canvas_axes.get_xlim()

        # if x-axis has changed, ensure bounds are appropriate
        xmin = self.minppm
        xmax = self.maxppm
        x0, y0, x1, y1 = self.canvas_axes.dataLim.bounds
        self.canvas_axes.ignore_existing_data_limits = True
        self.canvas_axes.update_datalim([[xmin,y0],[xmax,y1+y0]])

        self.canvas_axes.set_xlim((self.maxppm, self.minppm))
        self.build_basis_functions()
        self.calc_apodization()
        self.plot_canvas()
        self.canvas_axes.set_xlim(xlim)
        self.plot_canvas()


    def on_gaussian(self, event):
        # there are two line shape styles, Gaussian (if checked) and
        # Lorentzian if not checked.
        self.calc_apodization()
        self.plot_canvas()


    def on_magnitude(self, event):
        # there are two modes, Real data plot if not checked and Magnitude
        # data if checked.  These are tested for in plot_canvas()
        self.plot_canvas()


    def on_phase0(self, event):
        val = event.GetEventObject().GetValue()
        self.phase0 = val
        self.plot_canvas()

    def on_phase1(self, event):
        val = event.GetEventObject().GetValue()
        self.phase1 = val
        self.plot_canvas()

    def on_pivot(self, event):
        val = event.GetEventObject().GetValue()
        self.pivot = val
        self.plot_canvas()


    def on_run(self, event):
        # this tests the pulse sequence object for one metabolite and one
        # set of loop values.
        self.say("\nStarting test simulation...\n")

        msg = self.do_run()

        if msg:
            self.say(msg + "\n")
            self.say("Test simulation finished with errors.\n")
        else:
            self.say("Test simulation finished successfully.\n")
            self.build_basis_functions()
            self.plot_canvas()
            self.say("Results plotted to display canvas.\n")

        self.TextConsole.AppendText(SEPARATOR)

        self.ButtonRun.Enable()


    def on_results(self, event):
        # pops up a textual description of the test Experiment and results
        # in a native text editor.  This allows a user to keep a copy of a
        # previous test run to compare to subsequent runs.
        self.display_results_text()


    def on_screen_shot(self, event):
        # pops up a JPEG screen shot of the plot canvas in a native image
        # viewer. This allows a user to keep a copy of a previous test run
        # to compare to subsequent runs.

        figure = self.visualize.canvas.figure
        fd, filename = tempfile.mkstemp(".png")
        os.close(fd)

        figure.savefig( filename,
                        dpi=300,
                        facecolor='w',
                        edgecolor='w',
                        orientation='portrait',
                        #papertype='letter',
                        format=None,
                        transparent=False)

        common_wx_util.display_file(filename)

    def on_experiments(self, event):
        dialog = dialog_experiment_list.DialogExperimentList(self, self.pulse_sequence)
        dialog.ShowModal()

    def on_ok(self, event):
        if self.validate_pulse_sequence():
            # Save my coordinates
            config = sim_config.Config()
            config.set_window_coordinates("dialog_pulse_sequence_editor", self._left, self._top, self._width, self._height)
            config.write()
            # close dialog
            self.EndModal(wx.ID_OK)


    ##### Internal helper functions  ##########################################

    def initialize_controls(self):
        # sets up widgets in designer tab and test tab (including test
        # experiment values) from values in either a new or existing
        # mrs_pulse_sequence object.
        #
        # a few of the FloatSpin controls need to have limits, steps, etc.
        # set dynamically at instantiation.
        #
        # finally if the pulse sequence is 'frozen' most of the widgets
        # are disabled since only the name and comment should be editable in
        # that state

        self.initialize_auinotebook_panes()
        # We add the OK & Cancel buttons dynamically so that they're in the
        # right order under OS X, GTK, Windows, etc.
        self.ButtonOk, self.ButtonCancel = \
                    common_wx_util.add_ok_cancel(self,
                                                 self.LabelOkCancelPlaceholder,
                                                 self.on_ok)
        self.initialize_design_controls()
        self.initialize_experiment_controls()
        self.initialize_test_controls(first_time=True)

        self.FloatSpectralPoints.multiplier = 2
        self.FloatSpectralPoints.SetDigits(0)
        self.FloatSpectralPoints.SetIncrement(1.0)
        self.FloatSpectralPoints.SetRange(256.0,32768.0)
        self.FloatSpectralPoints.SetValue(float(common_constants.DEFAULT_SPECTRAL_POINTS))

        self.FloatSweepWidth.SetDigits(1)
        self.FloatSweepWidth.SetIncrement(500.0)
        self.FloatSweepWidth.SetRange(100.0,100000.0)
        self.FloatSweepWidth.SetValue(float(common_constants.DEFAULT_SWEEP_WIDTH))

        self.FloatLinewidth.SetDigits(1)
        self.FloatLinewidth.SetIncrement(0.5)
        self.FloatLinewidth.SetRange(0.001,1000.0)
        self.FloatLinewidth.SetValue(float(common_constants.DEFAULT_LINEWIDTH))

        self.FloatPhase0.SetDigits(1)
        self.FloatPhase0.SetIncrement(1.0)
        self.FloatPhase0.SetRange(-360.0,360.0)
        self.FloatPhase0.SetValue(0.0)

        self.FloatPhase1.SetDigits(1)
        self.FloatPhase1.SetIncrement(10.0)
        self.FloatPhase1.SetRange(-50000.0,50000.0)
        self.FloatPhase1.SetValue(0.0)

        self.FloatPivot.SetDigits(2)
        self.FloatPivot.SetIncrement(0.1)
        self.FloatPivot.SetRange(-20.0,20.0)
        self.FloatPivot.SetValue(self.pivot)


        if self.pulse_sequence.is_frozen:
            controls = (self.TextCreator,
                        self.TextLoop1Label,
                        self.TextLoop2Label,
                        self.TextLoop3Label,
                        self.ButtonAddStaticParameter,
                        self.ButtonRemoveStaticParameter,
                        self.ButtonAddPulse,
                        self.ButtonRemovePulse,
                        self.PanelTest)

            for control in controls:
                control.Disable()

            for label in self.loop_labels:
                label.Disable()

            msg = "Welcome to the Pulse Sequence Editor.\n" + \
                  "\nThe sequence that you are editing is 'frozen' so only "+ \
                  "the name and comment fields can be edited.\n"
        else:
            msg = "Welcome to the Pulse Sequence Editor.\n" + \
                  "\nStatus messages will appear here as you test your code.\n"

        self.TextConsole.Clear()
        self.TextConsole.AppendText(msg)

        if not self.pulse_sequence.experiment_names:
            # No point in showing this button which will only display an
            # empty list.
            self.ButtonListExperiments.Hide()



    def initialize_auinotebook_panes(self):

        self.AuiNotebookSizer = self.TextAuiPlaceholder.GetContainingSizer()
        self.TextAuiPlaceholder.Destroy()

        style =  aui.AUI_NB_TAB_SPLIT | \
                 aui.AUI_NB_TAB_MOVE | \
                 aui.AUI_NB_TAB_EXTERNAL_MOVE | \
                 wx.NO_BORDER | \
                 aui.AUI_NB_BOTTOM

        self.tab_notebook = aui.AuiNotebook(self.PanelAuiNotebook, style=style)

        self.panel_canvas = wx.Panel(self.tab_notebook, -1)
        self.visualize = plot_panel_editor.PlotPanelEditor(self.panel_canvas,
                                                           self,
                                                        self.labelStatus,
                                                        naxes=1,
                                                        reversex=True,
                                                        zoom='box',
                                                        reference=True,
                                                        middle=True,
                                                        do_zoom_select_event=True,
                                                        do_zoom_motion_event=True,
                                                        do_refs_select_event=True,
                                                        do_refs_motion_event=True,
                                                        do_middle_select_event=True,
                                                        do_middle_motion_event=True,
                                                        xscale_bump = 0.02,
                                                        yscale_bump = 0.06,
                                                        props_zoom=dict(alpha=0.2, facecolor='yellow'),
                                                        props_cursor=dict(alpha=0.1, facecolor='gray') )

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.visualize, 1, wx.LEFT | wx.TOP | wx.EXPAND)
        self.panel_canvas.SetSizer(sizer)
        self.visualize.Fit()

        # NB. the End Of Line (EOL) mode below is set to 2 which results in
        # an LF only, the default mode is 0 = CR,LF and mode = 1 is CR only
        self.tab_sequence_code = our_stc.PythonSTC(self, -1)
        self.tab_sequence_code.EmptyUndoBuffer()
        self.tab_sequence_code.SetEOLMode(2)
        self.tab_sequence_code.Colourise(0, -1)
        # line numbers in the margin
        self.tab_sequence_code.SetMarginType(1, wx_stc.STC_MARGIN_NUMBER)
        self.tab_sequence_code.SetMarginWidth(1, 25)
        if self.pulse_sequence.is_frozen:
            if self.pulse_sequence.sequence_code:
                # need to add code prior to setting 'read only'
                code = self.pulse_sequence.sequence_code
                self.tab_sequence_code.EmptyUndoBuffer()
                self.tab_sequence_code.SetText(code)
            else:
                self.tab_sequence_code.SetText("")
            self.tab_sequence_code.SetReadOnly(True)
        else:
            self.tab_sequence_code.SetText("")

        self.tab_binning_code = our_stc.PythonSTC(self, -1)
        self.tab_binning_code.EmptyUndoBuffer()
        self.tab_binning_code.SetEOLMode(2)
        self.tab_binning_code.Colourise(0, -1)
        # line numbers in the margin
        self.tab_binning_code.SetMarginType(1, wx_stc.STC_MARGIN_NUMBER)
        self.tab_binning_code.SetMarginWidth(1, 25)
        if self.pulse_sequence.is_frozen:
            if self.pulse_sequence.binning_code:
                code = self.pulse_sequence.binning_code
                self.tab_binning_code.EmptyUndoBuffer()
                self.tab_binning_code.SetText(code)
            self.tab_binning_code.SetReadOnly(True)
        else:
            self.tab_binning_code.SetText(common_constants.DEFAULT_BINNING_CODE)

        # The name of the stock tab art changed somewhere between wx 2.8 and 3.0.
        if hasattr(aui, 'AuiDefaultTabArt'):
            tab_art_provider = aui.AuiDefaultTabArt
        else:
            tab_art_provider = aui.AuiGenericTabArt
        self.tab_notebook.SetArtProvider(tab_art_provider())
        self.tab_notebook.AddPage(self.tab_sequence_code," Sequence Code ",False)
        self.tab_notebook.AddPage(self.tab_binning_code," Binning Code ",False)
        self.tab_notebook.AddPage(self.panel_canvas," Visualize ",False)

        sizer = self.AuiNotebookSizer
        sizer.Add(self.tab_notebook, 1, wx.EXPAND, 0)

        # Plot_Panel canvas initializations ----------
        td = 1.0/self.freq1d.sw
        xx = np.arange(self.freq1d.dims[0], dtype='float64')

        self.apodization = np.exp(-(0.6 * PI * self.linewidth * td * xx)**2)
        xx = self.freq1d.pts2ppm(xx)

        self.visualize.figure.set_facecolor(self.background)
        self.visualize.figure.subplots_adjust(left=0.0,right=0.999,bottom=0.001,top=1.0,wspace=0.0,hspace=0.01)
        self.canvas_axes = self.visualize.axes[0]
        self.canvas_axes.set_facecolor(self.background)
        self.canvas_axes.set_xlim(xx.max(),xx.min())


    def initialize_design_controls(self):

        sequence = self.pulse_sequence

        # a persistent reference to this sizer to add/remove controls
        self.DesignParameterSizer = self.InfoLabelParametersPlaceholder.GetContainingSizer()
        self.InfoLabelParametersPlaceholder.Destroy()

        # we changed wxGlade file pulse_sequence_editor.wxg to just declare an empty 0x2 row/col
        # FlexGrid sizer with the name self.DesignPulseSizer stored as an attribute.

        # I force this to a minimum size so that it is at least wide enough
        # to be useable. The size I chose is arbitrary.
        self.TextComment.SetMinSize( (350, -1) )

        # I create loop label aliases that are easier to work with than the
        # names that wxGlade assigns.
        self.loop_labels = [ ]
        for i in range(1, common_constants.RESULTS_SPACE_DIMENSIONS):
            self.loop_labels.append(getattr(self, "TextLoop%dLabel" % i))

        # design_parameter_controls is a list of lists. The inner lists contain
        # one row each of parameter controls.
        self.design_parameter_controls = []

        # design_pulse_controls is a list of lists. The inner lists contain
        # one row each of pulse controls.
        self.design_pulse_controls = []


        self.TextName.ChangeValue(sequence.name)
        self.LabelUuid.SetLabel(sequence.id)
        self.TextCreator.ChangeValue(sequence.creator)
        created = sequence.created.strftime(util_time.DISPLAY_DATE_FORMAT)
        self.LabelCreated.SetLabel(created)
        self.TextComment.ChangeValue(sequence.comment)

        for i, label in enumerate(sequence.loop_labels):
            self.loop_labels[i].ChangeValue(label)

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

        if sequence.user_static_parameters:
            for parameter in sequence.user_static_parameters:
                self.on_add_static_parameter(None)

                row = self.design_parameter_controls[-1]

                _, combobox, _, text_name, _, text_value = row

                combobox.SetStringSelection(parameter.type)
                text_name.ChangeValue(parameter.name)
                text_value.ChangeValue(parameter.default)
        else:
            # Create a single blank row
            self.on_add_static_parameter()

        self.DesignParameterSizer.Layout()

        # Here I hide all existing controls and then figure out
        # how many rows of pulse controls I need for the update.
        # I create one row of pulse controls for each pulse, or one empty
        # row if there are no pulse & we're in edit mode.
        for i, parameter in enumerate(self.design_pulse_controls):
            row = self.design_pulse_controls[i]
            checkbox, text_name = row
            checkbox.SetValue(False)
            text_name.SetLabel("")
            for control in row:
                control.Hide()

        if sequence.pulse_projects:
            sizer = self.DesignPulseSizer
            sizer.SetRows(len(sequence.pulse_projects) + 1)
            for pulse_design in sequence.pulse_projects:
                tmp_name = pulse_design.name+'    [uid = '+pulse_design.id+']'
                checkbox, label = self.add_pulse_row(pulse_design)
                label.SetLabel(tmp_name)
                sizer.Add(checkbox,   0, wx.ALIGN_CENTER_VERTICAL|wx.RIGHT, 5)
                sizer.Add(label, 1, wx.LEFT|wx.ALIGN_CENTER_VERTICAL, 5)

            self.PanelYourRFPulses.Layout()

        self.Layout()
        self.Fit()

        # Sequence & Binning Code tabs

        code = ""
        if sequence.sequence_code:
            code = sequence.sequence_code
        self.tab_sequence_code.EmptyUndoBuffer()
        self.tab_sequence_code.SetText(code)

        code = ""
        if sequence.binning_code:
            code = sequence.binning_code
        self.tab_binning_code.EmptyUndoBuffer()
        self.tab_binning_code.SetText(code)


    def initialize_experiment_controls(self):
        # these are set up separately from the other controls on the test
        # tab to simplify the update of the pulse sequence widget values
        # when the user wants to update them from the designer tab.
        # When that happens, the values for the test Experiment do not
        # need to be changed.
        self.ComboIsotope.AppendItems(self.db.fetch_isotopes())
        self.ComboIsotope.SetStringSelection(self.experiment.isotope)
        self.TextB0.ChangeValue(str(self.experiment.b0))
        self.TextPeakSearchRangeLow.ChangeValue(str(self.experiment.peak_search_ppm_low))
        self.TextPeakSearchRangeHigh.ChangeValue(str(self.experiment.peak_search_ppm_high))
        self.TextBlendTolerancePpm.ChangeValue(str(self.experiment.blend_tolerance_ppm))
        self.TextBlendTolerancePhase.ChangeValue(str(self.experiment.blend_tolerance_phase))
        self.populate_metabolites()


    def initialize_test_controls(self, first_time=False):

        if first_time:
            # a persistent reference to this sizer to add/remove controls
            self.ParameterSizer = self.LabelParametersPlaceholder.GetContainingSizer()
            # Make the column that holds the user's input growable.
            self.ParameterSizer.AddGrowableCol(1, 1)

            self.LabelParametersPlaceholder.Destroy()
            self.LoopSizer = self.LabelLoop1.GetContainingSizer()
            # parameter_controls is a list of lists. The inner lists contain
            # one row each of parameter controls.
            self.parameter_controls = []

            # I create aliases for these that are easier to work with than the
            # names that wxGlade assigns.

            self.loop_controls = []
            for i in range(1, common_constants.RESULTS_SPACE_DIMENSIONS):
                label = getattr(self, "LabelLoop%d" % i)
                start = getattr(self, "TextLoop%dStart" % i)
                self.loop_controls.append(LoopControlSet(label, start))

        self.experiment.pulse_sequence = self.pulse_sequence

        loop_labels = self.pulse_sequence.loop_labels
        parameters = self.pulse_sequence.user_static_parameters

        for i, loop_label in enumerate(loop_labels):
            control_set = self.loop_controls[i]
            control_set.label.SetLabel(loop_label+":")
            control_set.start.ChangeValue("")
            control_set.label.Enable()
            control_set.start.Enable()

        for i, control_set in enumerate(self.loop_controls[len(loop_labels):]):
            control_set.label.SetLabel("Loop %d:" % (len(loop_labels)+i+1))
            control_set.start.ChangeValue("")
            control_set.label.Disable()
            control_set.start.Disable()

        self.LoopSizer.Layout()

        sizer = self.ParameterSizer

        # There are three columns, one for the descriptive label, one for the
        # textbox that holds the default value, and one for another label that
        # describes the type (string, double, etc.)
        sizer.SetCols(3)
        sizer.SetRows(len(parameters))

        # Get rid of any existing controls
        for control_group in self.parameter_controls:
            for control in control_group:
                control.Destroy()
        self.parameter_controls = []

        # Create one row of new controls for each param
        for i, parameter in enumerate(parameters):
            controls = [ ]
            name = parameter.name
            if not name.endswith(":"):
                name += ":"
            name_label = wx.StaticText(self.PanelTest, wx.ID_ANY, name)
            sizer.Add(name_label, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_RIGHT|wx.RIGHT, 5)

            textbox = wx.TextCtrl(self.PanelTest)
            textbox.ChangeValue(parameter.default)
            sizer.Add(textbox, 1, wx.LEFT|wx.RIGHT|wx.EXPAND, 5)

            type_label = wx.StaticText(self.PanelTest, wx.ID_ANY, "(%s)" % parameter.type)
            sizer.Add(type_label, 0, wx.ALIGN_CENTER_VERTICAL)

            self.parameter_controls.append( (name_label, textbox, type_label) )

        # wx sets tab order according to control creation order. Since I
        # just created controls, they'll be *after* the Run button in the
        # tab order which is wrong. Here I correct the tab order.
        if self.parameter_controls:
            last_control = self.parameter_controls[-1][-1]
        else:
            last_control = self.loop_controls[-1].start

        self.PanelTest.Layout()


    def populate_metabolites(self):
        # get list of metabolite names from the DB and set them in ComboBox
        self.ComboMetabolites.Clear()
        isotope = self.ComboIsotope.GetStringSelection()
        metabolites = self.db.fetch_metabolites(isotope, False)
        for metabolite in metabolites:
            self.ComboMetabolites.Append(metabolite.name, metabolite)


    def add_parameter_row(self):
        index = len(self.design_parameter_controls)

        # These controls need to be children of the same window as the other
        # controls. I think that's the notebook, but there's no point to me
        # guessing when I can just ask wx.
        parent = self.TextName.GetParent()

        types = ("Double", "Long", "String")
        row = [ ]
        checkbox = wx.CheckBox(parent)
        row.append(checkbox)
        if self.pulse_sequence.is_frozen:
            # In read-only mode, the checkboxes are useless. I hide them
            # rather than destroying them because the grid sizer still
            # expects a certain number of columns per row.
            checkbox.Hide()

        combobox = wx.ComboBox(parent, style=wx.CB_DROPDOWN|wx.CB_READONLY,
                               size=(100, -1), choices=types)
        # Under OS X, wx2.9 selects the first item in the list automatically.
        combobox.SetSelection(-1)
        if self.pulse_sequence.is_frozen:
            combobox.Disable()
        row.append(combobox)

        row.append(wx.StaticText(parent, wx.ID_ANY, "Name:"))

        textbox = wx.TextCtrl(parent)
        if self.pulse_sequence.is_frozen:
            textbox.Disable()
        row.append(textbox)

        row.append(wx.StaticText(parent, wx.ID_ANY, "Default Value:"))

        textbox = wx.TextCtrl(parent)
        if self.pulse_sequence.is_frozen:
            textbox.Disable()
        row.append(textbox)

        self.design_parameter_controls.append(row)

        return row


    def add_pulse_row(self, pulse_design):
        parent = self.PanelYourRFPulses

        row = [ ]
        checkbox = wx.CheckBox(parent)
        if self.pulse_sequence.is_frozen:
            # In read-only mode, the checkboxes are useless. I hide them
            # rather than destroying them because the grid sizer still
            # expects a certain number of columns per row.
            checkbox.Hide()

        label = wx.StaticText(parent, wx.ID_ANY, pulse_design.name)

        self.design_pulse_controls.append( (checkbox, label) )

        return (checkbox, label)


    def validate_pulse_sequence(self):
        # Validate controls one by one
        # A non-empty string message indicates failure

        msg = ""

        # name
        name = self.TextName.GetValue().strip()
        if not name:
            msg = "Please enter a name for the pulse sequence."

        # creator (optional)
        if not msg:
            creator = self.TextCreator.GetValue().strip()

        # comment
        if not msg:
            # Note that I don't call .strip() on this field. It's totally
            # freeform and whitespace might be important.
            comment = self.TextComment.GetValue()

        # parameters.
        if not msg:
            for i, row in enumerate(self.design_parameter_controls):
                _, combobox, _, text_name, _, text_value = row
                combobox = bool(combobox.GetStringSelection())
                text_name = bool(text_name.GetValue())
                text_value = bool(text_value.GetValue())

                if (combobox != text_name) or (text_name != text_value):
                    msg = "Please complete parameter %d." % (i + 1)

        # RF Pulses need to be checked to make sure they still exist. It's
        # possible for a user to add a reference to one here in Simulation
        # then switch to Pulse and delete the referenced pulse design.
        if not msg:
            for pulse_design in self.pulse_sequence.pulse_projects:
                if not self.db.count_pulse_designs(pulse_design.id):
                    msg = """The pulse design "%s" no longer exists. """   \
                          """Please remove the reference to it."""
                    msg %= pulse_design.name
                    break

        # sequence code.
        if not msg:
            # Since this is someone's code, whitespace might be important
            # and I don't want to strip() it. However, I want to detect
            # the situation where a user has entered only whitespace, so
            # it's appropriate to strip() there.
            sequence_code = self.tab_sequence_code.GetText()

            if not sequence_code.strip():
                msg = "In order to save/test a pulse sequence, there must \n"+\
                      "be PyGamma code in the Sequence Code window. \n\n"+\
                      "Please enter code for the pulse sequence."

        # binning code.
        if not msg:
            binning_code = self.tab_binning_code.GetText()

            if not binning_code.strip():
                msg = "In order to save/test a pulse sequence, there must\n"+\
                      "be PyGamma code in the Binning Code window. \n\n"+\
                      "Please enter code for the binning step."

        if not msg:
            # I still have yet to check whether or not this name is unique.
            # I save it until last because it is a database hit.
            pulse_sequences = self.db.fetch_pulse_sequences_by_name(name)
            ids = [pulse_sequence.id for pulse_sequence in pulse_sequences \
                               if pulse_sequence.id != self.pulse_sequence.id]

            if ids:
                msg = "A pulse sequence with this name already exists.\n"+\
                      "Please change the name to a unique value."

        if msg:
            # validation failed
            common_dialogs.message(msg, None, common_dialogs.X_OK)
            return False
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

            self.pulse_sequence.user_static_parameters = []
            for row in self.design_parameter_controls:
                _, combobox, _, text_name, _, text_value = row
                # Not all rows contain data. Don't forget to ignore empty ones.
                type_ = combobox.GetStringSelection()

                if type_:
                    parameter = mrs_pulse_sequence.UserStaticParameter()
                    parameter.type = type_
                    parameter.name = text_name.GetValue().strip()
                    parameter.default = text_value.GetValue()
                    self.pulse_sequence.user_static_parameters.append(parameter)

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

        previous = self.ButtonRemoveStaticParameter
        for control in controls:
            control.MoveAfterInTabOrder(previous)
            previous = control

        controls = [ ]
        for row in self.design_pulse_controls:
            controls += row

        previous = self.ButtonRemovePulse
        for control in controls:
            control.MoveAfterInTabOrder(previous)
            previous = control


    def do_run(self):
        """

        Runs the single-simulation experiment. Returns "" if the run
        succeeded or a message suitable for display in the console window.
        """
        # Validate controls one by one
        msg = ""

        self.ButtonRun.Disable()

        name         = 'pulse sequence designer test name'
        investigator = 'pulse sequence designer test investigator'
        comment      = 'pulse sequence designer test comment'

        # pulse sequence loops
        if not msg:
            for i, control_set in enumerate(self.loop_controls):
                if control_set.start.IsEnabled():
                    start = control_set.start.GetValue().strip()
                    if start:
                        if not util_misc.is_floatable(start):
                            msg = """I don't understand the start value "%s".""" % start
                    else:
                        msg = "Please enter a start value for loop %d." % (i + 1)
                if msg:
                    # No point in going through the other controls.
                    break

        # pulse sequence static user parameters
        if not msg:
            for parameters in self.parameter_controls:

                key   = parameters[0].GetLabel().strip()
                value = parameters[1].GetValue().strip()
                type_  = parameters[2].GetLabel().strip()
                if type_ == "(Double)":
                    if not util_misc.is_floatable(value):
                        msg = """I don't understand the %s parameter value "%s".""" % (key,value)
                elif type_ == "(Long)":
                    if not util_misc.is_intable(value):
                        msg = """I don't understand the %s parameter value "%s".""" % (key,value)
                elif type_ == "(String)":
                    if value == '':
                        msg = """Please enter a value for the "%s" parameter".""" % key

                if msg:
                    # No point in going through the other controls.
                    break

        # B0
        if not msg:
            b0 = self.TextB0.GetValue().strip()
            if b0:
                if not util_misc.is_floatable(b0):
                    msg = """I don't understand the B0 value "%s".""" % b0
            else:
                msg = "Please enter a B0 value."

        # metabs
        if not msg:
            if self.ComboMetabolites.GetSelection()==wx.NOT_FOUND:
                msg = "Please add a metabolite to this experiment."

        # peak search low
        if not msg:
            peak_search_ppm_low = self.TextPeakSearchRangeLow.GetValue().strip()
            if peak_search_ppm_low:
                if not util_misc.is_floatable(peak_search_ppm_low):
                    msg = """I don't understand the peak search ppm low value "%s".""" % peak_search_ppm_low
            else:
                msg = "Please enter a peak search ppm low value."

        # peak search high
        if not msg:
            peak_search_ppm_high = self.TextPeakSearchRangeHigh.GetValue().strip()
            if peak_search_ppm_high:
                if not util_misc.is_floatable(peak_search_ppm_high):
                    msg = """I don't understand the peak search ppm high value "%s".""" % peak_search_ppm_high
            else:
                msg = "Please enter a peak search ppm high value."

        # blend tolerance PPM
        if not msg:
            blend_tolerance_ppm = self.TextBlendTolerancePpm.GetValue().strip()
            if blend_tolerance_ppm:
                if not util_misc.is_floatable(blend_tolerance_ppm):
                    msg = """I don't understand the blend tolerance value "%s".""" % blend_tolerance_ppm
            else:
                msg = "Please enter both blend tolerance values."

        # blend tolerance phase
        if not msg:
            blend_tolerance_phase = self.TextBlendTolerancePhase.GetValue().strip()
            if blend_tolerance_phase:
                if not util_misc.is_floatable(blend_tolerance_phase):
                    msg = """I don't understand the blend tolerance value "%s".""" % blend_tolerance_phase
            else:
                msg = "Please enter both blend tolerance values."

        if msg:
            common_dialogs.message(msg, None, common_dialogs.I_OK)
            self.ButtonRun.Enable()
            return "\nError in test Experiment parameters setup in Test tab."

        else:

            # All is well, time to Run Simulations ----------------------------
            #
            # this section copies values from widgets into the Experiment
            # object for a number of input parameters

            self.experiment.name                    = name
            self.experiment.investigator            = investigator
            self.experiment.comment                 = comment
            self.experiment.isotope                 = self.ComboIsotope.GetStringSelection()
            self.experiment.b0                      = float(b0)
            self.experiment.peak_search_ppm_low     = float(peak_search_ppm_low)
            self.experiment.peak_search_ppm_high    = float(peak_search_ppm_high)
            self.experiment.blend_tolerance_ppm     = float(blend_tolerance_ppm)
            self.experiment.blend_tolerance_phase   = float(blend_tolerance_phase)
            self.experiment.pulse_sequence          = self.pulse_sequence

            self.experiment.pulse_sequence.sequence_code = self.tab_sequence_code.GetText()
            self.experiment.pulse_sequence.binning_code  = self.tab_binning_code.GetText()

            self.experiment.user_static_parameters = []
            for control_group in self.parameter_controls:
                name_label, textbox, type_label = control_group
                self.experiment.user_static_parameters.append(textbox.GetValue())

            # get test metabolite selected for simulation run
            index = self.ComboMetabolites.GetSelection()
            metabolite = self.ComboMetabolites.GetClientData(index)
            self.experiment.metabolites = [metabolite]

            self.experiment.dims = [ ]
            for i, control_set in enumerate(self.loop_controls):
                if control_set.start.IsEnabled():
                    start  = float(control_set.start.GetValue())
                    step   = 0.0
                    length = 1
                    dim = mrs_experiment.expand_loop(start, step, length)
                else:
                    dim = mrs_experiment.DEFAULT_LOOP

                self.experiment.dims.append(dim)

            results, failed = run_experiment_controller.run_experiment(self.experiment, local_cpu=1)

            self.ButtonRun.Enable()

            if failed:
                message = \
                    run_experiment_controller.exception_to_message(failed[0])
            else:
                message = ""
                # Since we only asked for one simulation, we'll only get one
                # result.
                results = results[0]

                # The dims in the results contain only the metab name; we
                # want the entire metab object. Fortunately we already know
                # what it is.
                results["metabolite"] = metabolite

                self.experiment.simulations = [mrs_simulation.Simulation(results)]

            return message


    def display_results_text(self):

        lines = str(self.experiment)
        lines += "\n\nSimulation Results\n" + "-" * 75 + "\n\n"

        lines += "\n".join([simulation.summary() for simulation
                                                 in self.experiment.simulations])

        common_wx_util.display_text_as_file(lines)


    def plot_canvas(self, reset_history=False):

        # check if any egregious issues with proceeding
        if not self.experiment.metabolites:
            return
        if self.basis.data is None:
            return

        flag_magnitude = self.CheckMagnitude.IsChecked()

        b = self.basis.data[0,0,0,0,:]
        data = b.copy()

        # transform into frequency domain and store
        data[0] = data[0] / 2.0
        data[:] = np.fft.fft(data[:] * self.apodization)/float(self.freq1d.dims[0])
        dim0    = self.freq1d.dims[0]
        pivot   = self.freq1d.ppm2pts(self.pivot)
        phase0  = self.phase0 * common_constants.DEGREES_TO_RADIANS
        phase1  = self.phase1 * common_constants.DEGREES_TO_RADIANS
        phase1  = phase1 * (np.arange(dim0) - pivot) / dim0
        phase   = np.exp(1j * (phase0 + phase1))

        self.freq1d.data = data * phase

        # select data type to plot - real, imaginary, magnitude
        if self.CheckMagnitude.IsChecked():
            data = abs(self.freq1d.data)
        else:
            data = self.freq1d.data.real

        # x-axis creation
        xx = np.arange(self.freq1d.dims[0], dtype='float64')
        xx = self.freq1d.pts2ppm(xx)
        xlim = self.canvas_axes.get_xlim()
        xmin = (min(xlim)<min(xx)).choose(min(xlim),min(xx))
        xmax = (max(xlim)>max(xx)).choose(max(xlim),max(xx))
        xlim = (xmax, xmin)

        # set up dynamic y-scaling based on min/max of data plotted
        dmin = data.min()
        dmax = data.max()
        dr   = dmax - dmin
        dmin = dmin - dr*0.03   # space line collection vertically
        dmax = dmax + dr*0.03
        dr = dmax - dmin
        y0 = dmin
        y1 = dmin + dr

        self.canvas_axes.cla()
        self.visualize.refresh_cursors()

        self.canvas_axes.plot(xx, data,
                              self.foreground,
                              linewidth=self._prefs.line_width)

        self.canvas_axes.autoscale_view(tight=True)

        self.canvas_axes.set_xlim(xlim)

        # show a zero line
        self.canvas_axes.axhline(y = 0*dr,
                                 color=self._prefs.zero_line_color,
                                 linestyle=self._prefs.zero_line_style,
                                 linewidth=self._prefs.line_width)

        # y-lim set needs to stay after zero_line or y-axis gets skewed
        self.canvas_axes.set_ylim(y0, y1)
        bob = self.canvas_axes.dataLim.bounds

        # show the x-axis
        self.visualize.figure.subplots_adjust(left=0.0,right=0.999,bottom=0.05,top=1.0,wspace=0.0,hspace=0.01)

        self.visualize.canvas.draw()


    def build_basis_functions(self):
        if not self.experiment.metabolites:
            return

        if self.experiment.isotope == '1H':
            ppm = common_constants.DEFAULT_PROTON_CENTER_PPM
        else:
            ppm = common_constants.DEFAULT_XNUCLEI_CENTER_PPM

        self.freq1d.resppm     = ppm

        # ensure that datasets and certain boundaries
        # reflect the experiment settings

        self.freq1d.frequency = self.experiment.b0
        self.basis.frequency = self.experiment.b0
        self.maxppm = self.freq1d.pts2ppm(0)
        self.minppm = self.freq1d.pts2ppm(self.freq1d.dims[0]-1)

        self.basis.data = bbf_module.build_basis_functions(self.experiment,
                                                           self.freq1d.dims[0],
                                                           self.freq1d.sw,
                                                           self.freq1d.resppm,
                                                           local_cpu=1)

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

