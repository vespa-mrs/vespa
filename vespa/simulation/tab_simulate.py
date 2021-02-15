# Python modules
import math
import os

# 3rd party modules
import wx
import numpy as np

# Our modules
import vespa.simulation.prefs as prefs
import vespa.simulation.constants as constants
import vespa.simulation.mrs_data_basis as mrs_data_basis
import vespa.simulation.build_basis_functions as bbf_module
import vespa.common.util.ppm as util_ppm


import vespa.simulation.util_simulation_config as util_simulation_config
import vespa.simulation.auto_gui.simulate as simulate
import vespa.common.mrs_experiment as mrs_experiment
import vespa.common.mrs_simulation as mrs_simulation
import vespa.common.wx_gravy.common_dialogs as common_dialogs
import vespa.common.wx_gravy.util as common_wx_util
import vespa.common.util.time_ as util_time
import vespa.common.util.misc as util_misc
import vespa.common.constants as common_constants
import vespa.simulation.run_experiment_controller as run_experiment_controller

PI = math.pi

#------------------------------------------------------------------------------
# Note. GUI Architecture/Style
#
# Most Vespa GUI components are designed using WxGlade for rapid developement.
# WxGlade is interactive and users can preview the resultant window/panel/dialog.
# However, event functions can be specified but only stub functions with those
# names are created. WxGlade (*.wxg) files are in the 'wxglade' subdirectory.
# WxGlade generated modules are stored in 'auto_gui' subdirectory.
#
# Each WxGlade created GUI class is inherited into a unique 'vespa' class, where
# program specific initialization and other functionality are written. The stub
# event functions are overloaded here to provide program specific event handling.
#------------------------------------------------------------------------------


def force_to_float(a_potential_float):
    """Given arbitrary data, returns a float, always. If a_potential_float
    cannot be made into a float, returns NaN (i.e. float("nan")) which is
    a float that's not equal to anything.
    """
    try:
        the_float = float(a_potential_float)
    except ValueError:
        the_float = float("NaN")

    return the_float


class LoopControlSet(object):
    """A container for a set of loop controls. They are the column label
    and the three textboxes for start, count and size"""
    def __init__(self, label, start, count, size):
        self.label = label
        self.start = start
        self.count = count
        self.size = size


class TabSimulate(simulate.PanelSimulateUI):

    def __init__(self, tab_experiment, top, db, experiment, is_new):

        simulate.PanelSimulateUI.__init__(self, tab_experiment)

        self._top = wx.GetApp().GetTopWindow()
        self._tab_experiment = tab_experiment
        self._prefs = prefs

        self.db         = db
        self.top        = top
        self.experiment = experiment
        self.statusbar  = top.statusbar

        self.pulse_sequences = [ ]

        # last_save & last_run contain the state of the GUI (represented as a
        # dict) was last save/run. This is almost the same as the state of
        # the experiment, but without the baggage of all the simulations.
        self.last_save      = None
        self.last_run       = None
        self.last_exception = None
        self.is_new         = is_new

        # I create aliases for the loop contols that are easier to work with
        # than the names that wxGlade assigns.
        self.loop_controls = [ ]
        for i in range(1, common_constants.RESULTS_SPACE_DIMENSIONS):
            label = getattr(self, "LabelLoop%d" % i)
            start = getattr(self, "TextLoop%dStart" % i)
            count = getattr(self, "TextLoop%dStepCount" % i)
            size  = getattr(self, "TextLoop%dStepSize" % i)
            self.loop_controls.append(LoopControlSet(label, start, count, size))

        # I need a persistent reference to this sizer for adding and removing controls
        self.ParameterSizer = self.LabelParametersPlaceholder.GetContainingSizer()
        self.LabelParametersPlaceholder.Destroy()

        self.parameter_controls = []

        self.initialize_controls()

        # Save the state of the experiment as described by the GUI.

        if is_new and experiment.simulations:
            self.last_save = { }
        else:
            self.last_save = self.get_raw_gui_data()

        if not self.is_new or experiment.simulations:
            self.last_run = self.get_raw_gui_data()

        self.Layout()
        self.Fit()

        # Grab key events so I can look for Ctrl/Cmd+A ==> select all
        self.ListMetabolitesInUse.Bind(wx.EVT_KEY_DOWN, self.on_key_down)
        self.ListMetabolitesAvailable.Bind(wx.EVT_KEY_DOWN, self.on_key_down)


    # =======================================================
    #
    #           GUI Setup Calls
    #
    # =======================================================

    def initialize_controls(self):
        """
        Initializes the controls to be the right size or have the right
        range or number of decimal places. It typically does not set the
        default value (that's for populate_controls method to do). This
        method does the one-time setup bits.

        """
        self.TextName.ChangeValue(self.experiment.name)
        self.LabelUuid.SetLabel(self.experiment.id)
        self.TextInvestigator.ChangeValue(self.experiment.investigator)
        created = self.experiment.created.strftime(util_time.DISPLAY_DATE_FORMAT)
        self.LabelCreated.SetLabel(created)
        self.TextComment.ChangeValue(self.experiment.comment)

        self.ComboIsotope.AppendItems(self.db.fetch_isotopes())
        self.ComboIsotope.SetStringSelection(self.experiment.isotope)

        self.populate_metabolites(True)

        self.TextB0.ChangeValue(str(self.experiment.b0))
        self.TextPeakSearchRangeLow.ChangeValue(str(self.experiment.peak_search_ppm_low))
        self.TextPeakSearchRangeHigh.ChangeValue(str(self.experiment.peak_search_ppm_high))
        self.TextBlendTolerancePpm.ChangeValue(str(self.experiment.blend_tolerance_ppm))
        self.TextBlendTolerancePhase.ChangeValue(str(self.experiment.blend_tolerance_phase))

        self.populate_pulse_sequences()

        textboxes = [parameter_control_set[1] for parameter_control_set \
                                              in self.parameter_controls]
        for parameter, textbox in \
                    zip(self.experiment.user_static_parameters, textboxes):
            textbox.ChangeValue(parameter)

        self.set_control_enabled()


    # =======================================================
    #
    #           Event Handlers
    #
    # =======================================================

    def on_combo_isotope(self, event):
        self.populate_metabolites()

    def on_button_left(self, event):
        self.move_metabolites(self.ListMetabolitesAvailable,
                              self.ListMetabolitesInUse)

    def on_button_right(self, event):
        self.move_metabolites(self.ListMetabolitesInUse,
                              self.ListMetabolitesAvailable)

    def on_key_down(self, event):
        if common_wx_util.is_select_all(event):
            focus = self.FindFocus()
            if focus == self.ListMetabolitesInUse:
                listbox = self.ListMetabolitesInUse
            elif focus == self.ListMetabolitesAvailable:
                listbox = self.ListMetabolitesAvailable
            else:
                listbox = None

            if listbox:
                for i in range(listbox.GetCount()):
                    listbox.Select(i)
        event.Skip()


    def on_pulse_sequence_selected(self, event):
        self.pulse_sequence_setup()


    def on_run(self, event):

        self.statusbar.SetStatusText('',3)
        self.ButtonRun.Disable()

        # Validate controls
        msg = self.is_gui_valid()

        if not msg:
            if self.last_run:
                # This has been run once before. Has it changed since then?
                current = self.get_raw_gui_data()

                # For purposes of determining whether or not the results
                # are current, these fields don't matter so I ensure that
                # they're always equal in the comparison below.
                for key in ("name", "comment", "investigator"):
                    current[key] = self.last_run[key]

                if current == self.last_run:
                    msg = "This experiment's results are already up to date."
        if msg:
            common_dialogs.message(msg, None, common_dialogs.I_OK)
            self.ButtonRun.Enable()
        else:
            # All is well, time to run simulations -------------------------

            self.statusbar.SetStatusText('Running experiment, please wait.',0)

            # Save these as default/last chosen values.
            wx.GetApp().vespa.isotope = self.experiment.isotope
            wx.GetApp().vespa.b0      = self.experiment.b0

            # Save the state of the experiment as described by the GUI.
            self.last_run = self.get_raw_gui_data()

            if self.is_new:
                # New experiments start from scratch every time.
                self.experiment.dims = [ [ ] ] * len(self.experiment.dims)
                self.experiment.metabolites = [ ]
                self.experiment.simulations = [ ]

            # When running an existing experiment (one that's been saved,
            # closed and re-opened), we know that the significant values
            # haven't changed. Therefore, we don't need to re-run all metabs,
            # just the ones that have been added.
            # In order to figure out which metabs are new, we have to note
            # here which ones are currently in the experiment.
            current_metabolite_ids = [metabolite.id for metabolite
                                                in self.experiment.metabolites]

            self.populate_experiment_from_gui()

            # Figure out which metabs were just added to the experiment.
            new_metabolites = [metabolite for metabolite
                                          in self.experiment.metabolites
                                          if metabolite.id
                                          not in current_metabolite_ids]

            results, exception = \
                run_experiment_controller.run_experiment(self.experiment,
                                                         new_metabolites)

            # Strange but true -- we reduce all exception information to
            # a boolean failed/did not fail.
            #
            # bjs - 2018-08 added a bit more info on return, BUT need to keep msg short!
            if exception:
                self.last_exception = exception
                self.last_run = {}
                self.statusbar.SetStatusText('Run returned errors!',2)
                self.statusbar.SetStatusText(str(exception[0][1]),3)
            else:
                self.last_exception = None

                # Turn the results (a list of dicts) into Simulation objects.

                # The result dicts contain only metabolite names, not UUIDs
                # or metabolite objects. Here we create a dictionary that
                # maps those names to metabolite objects.
                metabolite_map = [(metabolite.name, metabolite) for metabolite in new_metabolites]
                metabolite_map = dict(metabolite_map)

                for i, result in enumerate(results):
                    # Replace result["metabolite"] with the metab object
                    # rather than just the metab name
                    result["metabolite"] = metabolite_map[result["metabolite"]]

                    simulation = mrs_simulation.Simulation(result)
                    self.experiment.simulations.append(simulation)

                self.experiment.simulations.sort()

                # refresh Visualize and Simulate panels to reflect new states
                #  - turn on visualize widgets
                #  - turn off widgets in simulate tab that should not be altered
                self._tab_experiment.refresh_visualize()

                self.set_control_enabled()

            self.statusbar.SetStatusText(' ',0)

            # The run button is only disabled if the Experiment is "public"
            self.ButtonRun.Enable()


    #=======================================================
    #
    #           Internal methods start here
    #
    #=======================================================

    def close(self):
        """
        Asks the user if it is OK to close the experiment. Returns False
        if the user rejects the close (i.e. hits "cancel" on the "OK to
        close?" message box), True otherwise. If there are no changes,
        the user is not prompted at all (and True is returned).

        """
        gui = self.get_raw_gui_data()
        if self.last_save != gui:
            name = gui["name"]
            msg = """The experiment "%s" """ % name
            msg += "has unsaved changes. Simulation is about to discard those changes."
            answer = common_dialogs.message(msg, None,
                                    wx.ICON_INFORMATION | wx.OK | wx.CANCEL)

            return (answer == wx.OK)
        else:
            # The GUI hasn't changed since last save; it's OK to close
            return True


    def get_cooked_gui_data(self):
        """Returns a dict containing the data from the GUI in cooked form.
        They keys to the dict are strings representing the field names.

        "Cooked" is in contrast to the raw form returned by
        get_raw_gui_data() (q.v.). Raw data is cooked by (a) forcing strings
        to floats as appropriate. If a string that should be a float isn't a
        valid float, it's represented by the float value NaN as returned
        by force_to_float().

        The metabs and the pulse seq are untouched. The former is a list of
        Metabolite objects, the latter a PulseSequence object.
        """
        d = self.get_raw_gui_data()

        for key in ("b0", "peak_search_ppm_low", "peak_search_ppm_high",
                    "blend_tolerance_ppm", "blend_tolerance_phase"):
            d[key] = force_to_float(d[key])

        for loop in d["pulse_sequence_loops"]:
            loop["start"] = float(loop["start"])
            loop["step"] = float(loop["step"])
            loop["length"] = int(loop["length"])

        return d


    def get_raw_gui_data(self):
        """Returns a dict containing the data from the GUI in almost-raw
        form. They keys to the dict are strings representing the field names.
        In most cases they match the experiment attribute names which makes
        it possible to pass this dictionary to experiment.inflate() (after
        some verification and modification).

        "Almost-raw" means that strip() has been called on all of the
        text fields but nothing else has been done to the values received
        from the GUI. Objects (metabs & pulse seqs) are represented as
        objects, they're not deflated to dicts.

        Note: this method performs no type translation, so the value for
        e.g. d["b0"] is a string, not a float. Furthermore, this method
        doesn't guarantee validity, so d["b0"] might well be "My brain hurts!"
        instead of something valid.
        """
        d = { }

        d["name"] = self.TextName.GetValue().strip()
        d["investigator"] = self.TextInvestigator.GetValue().strip()
        # Note that we don't call .strip() on the comment. It's totally
        # freeform and whitespace might be important.
        d["comment"] = self.TextComment.GetValue()
        d["isotope"] = self.ComboIsotope.GetStringSelection()
        d["metabolites"] = [self.ListMetabolitesInUse.GetClientData(i) \
                        for i in range(self.ListMetabolitesInUse.GetCount())]
        d["b0"] = self.TextB0.GetValue().strip()
        d["peak_search_ppm_low"] = self.TextPeakSearchRangeLow.GetValue().strip()
        d["peak_search_ppm_high"] = self.TextPeakSearchRangeHigh.GetValue().strip()
        d["blend_tolerance_ppm"] = self.TextBlendTolerancePpm.GetValue().strip()
        d["blend_tolerance_phase"] = self.TextBlendTolerancePhase.GetValue().strip()
        d["pulse_sequence"] = None
        i = self.ComboPulseSequence.GetCurrentSelection()
        if i != wx.NOT_FOUND:
            d["pulse_sequence"] = self.pulse_sequences[i]

        d["pulse_sequence_loops"] = [ ]
        if d["pulse_sequence"]:
            # Figure out how many loop values there should be.
            loop_count = len(d["pulse_sequence"].loop_labels)
            loop_controls = self.loop_controls[:loop_count]
            # Go get 'em
            for control_set in loop_controls:
                loop_d = { }
                loop_d["start"] = control_set.start.GetValue().strip()
                loop_d["length"] = control_set.count.GetValue().strip()
                loop_d["step"] = control_set.size.GetValue().strip()
                d["pulse_sequence_loops"].append(loop_d)
        #else:
            # There's no pulse seq selected so there can't possibly be any
            # loop values.

        d["user_static_parameters"] = [ ]
        for parameters in self.parameter_controls:
            parameter = { }
            parameter["name"] = parameters[0].GetLabel().strip()
            # When we get the value for a user param, we don't strip() it.
            # Whitespace might be significant!
            parameter["value"] = parameters[1].GetValue()
            parameter["type"] = parameters[2].GetLabel().strip()

            d["user_static_parameters"].append(parameter)

        return d


    def is_gui_valid(self, check_for_unique_name=False):
        """
        Examines the contents of the GUI and determines whether or not
        what the user has entered is valid. If this method finds something
        it doesn't like, it returns a message string indicating the problem
        that's appropriate for displaying in a message box. Otherwise it
        returns None.
        """
        msg = None

        d = self.get_raw_gui_data()

        # name
        name = d["name"]
        if not name:
            msg = "Please enter a name for the experiment."

        if check_for_unique_name:
            if self.db.count_experiments(name=name):
                msg = "Please enter a unique name for the experiment."

        # metabs
        if not msg:
            if not d["metabolites"]:
                msg = "Please add at least one metabolite to this experiment."

        # B0
        if not msg:
            b0 = d["b0"]

            if b0:
               if not util_misc.is_floatable(b0):
                   msg = """I don't understand the main field value "%s".""" % b0
            else:
                msg = "Please enter a main field value."

        # peak search low
        if not msg:
            peak_search_ppm_low = d["peak_search_ppm_low"]

            if peak_search_ppm_low:
               if not util_misc.is_floatable(peak_search_ppm_low):
                   msg = """I don't understand the peak search range low value "%s".""" % peak_search_ppm_low
            else:
                msg = "Please enter a peak search range low value."

        # peak search high
        if not msg:
            peak_search_ppm_high = d["peak_search_ppm_high"]

            if peak_search_ppm_high:
               if not util_misc.is_floatable(peak_search_ppm_high):
                   msg = """I don't understand the peak search range high value "%s".""" % peak_search_ppm_high
            else:
                msg = "Please enter a peak search range high value."

        # blend tolerance PPM
        if not msg:
            blend_tolerance_ppm = d["blend_tolerance_ppm"]

            if blend_tolerance_ppm:
               if not util_misc.is_floatable(blend_tolerance_ppm):
                   msg = """I don't understand the blend tolerance value "%s".""" % blend_tolerance_ppm
            else:
                msg = "Please enter both blend tolerance values."

        # blend tolerance phase
        if not msg:
            blend_tolerance_phase = d["blend_tolerance_phase"]

            if blend_tolerance_phase:
               if not util_misc.is_floatable(blend_tolerance_phase):
                   msg = """I don't understand the blend tolerance value "%s".""" % blend_tolerance_phase
            else:
                msg = "Please enter both blend tolerance values."

        # pulse sequence
        if not msg:
            if not d["pulse_sequence"]:
                msg = "Please select a pulse sequence for this experiment."

        # pulse sequence loops
        if not msg:
            for i, loop in enumerate(d["pulse_sequence_loops"]):
                start = loop["start"]
                length = loop["length"]
                step = loop["step"]
                if start:
                   if not util_misc.is_floatable(start):
                       msg = """I don't understand the start value "%s".""" % start
                else:
                    msg = "Please enter a start value for loop %d." % (i + 1)

                if not msg:
                    if length:
                       if util_misc.is_intable(length):
                           length = int(length)
                           if length <= 0:
                               msg = "Please enter a step count > 0."
                       else:
                           msg = """I don't understand the step count "%s".""" % length
                    else:
                        msg = "Please enter a step count for loop %d." % (i + 1)

                if not msg:
                    if step:
                       if not util_misc.is_floatable(step):
                           msg = """I don't understand the step size "%s".""" % step
                    else:
                        msg = "Please enter a step size for loop %d." % (i + 1)

                if msg:
                    # No point in going through the other controls.
                    break

        # pulse sequence static user parameters
        if not msg:
            for parameter in d["user_static_parameters"]:
                name = parameter["name"]
                value = parameter["value"]
                type_ = parameter["type"]

                if type_ == "(Double)":
                    if not util_misc.is_floatable(value):
                        msg = """I don't understand the %s parameter value "%s".""" % (name, value)
                elif type_ == "(Long)":
                    if not util_misc.is_intable(value):
                        msg = """I don't understand the %s parameter value "%s".""" % (name, value)
                elif type_ == "(String)":
                    if not value:
                        msg = """Please enter a value for the "%s" parameter".""" % name

                if msg:
                    # No point in going through the other controls.
                    break

        return msg


    def move_metabolites(self, source, target):
        """Moves a metab from the source listbox to the target listbox"""
        indices = source.GetSelections()

        if indices:
            # It's very simple to remove metabs from one list and move them
            # to another BUT if I remove items from lowest to highest, the
            # higher indices will change. (e.g. If the user has selected
            # items 0 and 4, when I remove 0, 4 will move to 3.)
            # Therefore I have to ensure that my list of indices which
            # GetSelections() returns in arbitrary order is reverse sorted.
            for index in sorted(indices, reverse=True):
                metabolite = source.GetClientData(index)
                target.Append(metabolite.name, metabolite)
                source.Delete(index)

            self.populate_metabolites()
        else:
            msg = "Please select the metabolite(s) you want to move."
            common_dialogs.message(msg, None, common_dialogs.I_OK)


    def populate_experiment_from_gui(self):
        """Copies the editable elements of the GUI into self.experiment.
        The editable elements of the GUI include all fields that the user can
        edit on a new experiment. Notably, this does NOT include UUID and
        creation date.

        This function does no validation. It assumes that the values in the
        GUI are valid (e.g. ready to be passed to float() and int()).
        """
        cooked = self.get_cooked_gui_data()

        # The cooked GUI data is almost perfect as-is for passing to
        # experiment.inflate(). I just need to make a few modifications.
        # First, I simplify the the pulse seq params a little.
        cooked["user_static_parameters"] = \
                    [parameter["value"] for parameter
                                        in cooked["user_static_parameters"]]

        # Next, I assign the metabs and pulse seq directly since they
        # haven't been deflated to dicts (and in fact don't need to be).
        self.experiment.pulse_sequence = cooked["pulse_sequence"]
        self.experiment.metabolites = cooked["metabolites"]

        # Last but not least, I delete them from the dict so that
        # experiment.inflate() doesn't try to use them.
        del cooked["pulse_sequence"]
        del cooked["metabolites"]

        self.experiment.inflate(cooked)

        # Now the experiment is populated except for the dims which are not
        # so straightforward.

        # The values the user entered into the pulse seq loop controls (if
        # any) determine how many simulations the expt has. Each expt. has
        # 4 dimensions of simulations: 3 numeric plus a list of metabs. An
        # experiment's simulations attribute is a flattened list of the 4D
        # space.
        # We start with default values ([0]) for each dimension.
        dims = [mrs_experiment.DEFAULT_LOOP] * \
                              (common_constants.RESULTS_SPACE_DIMENSIONS - 1)

        # Translate the populated loops into dims
        for i, loop in enumerate(cooked["pulse_sequence_loops"]):
            dims[i] = mrs_experiment.expand_loop(loop["start"], loop["step"],
                                                 loop["length"])

        self.experiment.dims = dims


    def populate_metabolites(self, first_time=False):
        """Clears & repopulates the list of available metabs. The flag
        first_time should be True during init, False otherwise.
        """
        isotope = self.ComboIsotope.GetStringSelection()

        if first_time:
            # During init I need to seed the "in use" listbox
            in_use = self.experiment.metabolites
        else:
            # There's no guarantee that the metabs the user has currently
            # added to the experiment will still be valid after I fetch a new
            # list from the database. I preserve as many as possible.
            in_use = [self.ListMetabolitesInUse.GetClientData(i) for i in \
                                range(self.ListMetabolitesInUse.GetCount())]

        self.ListMetabolitesInUse.Clear()
        self.ListMetabolitesAvailable.Clear()

        metabolites = self.db.fetch_metabolites(isotope, False)

        in_use_id = [obj.id for obj in in_use]

        # Now metabolites holds the set of metabs appropriate for the
        # current isotope. Call this set C. The list in_use holds the set
        # of metabs selected by the user. Call this set U. At this point,
        # U is not necessarily a subset of C, but it should (and soon will)
        # be.

        # I'll call the set of available metabs A. A = C - U
        available = [metabolite for metabolite in metabolites
                                                if metabolite.id not in in_use_id]

        # U = U intersection C. I could use Python sets to perform these
        # set operations, but sets don't preserve order and I want to
        # retain the database's ordering.
        in_use = [metabolite for metabolite in metabolites
                                            if metabolite.id in in_use_id]

        for metabolite in in_use:
            self.ListMetabolitesInUse.Append(metabolite.name, metabolite)

        for metabolite in available:
            self.ListMetabolitesAvailable.Append(metabolite.name, metabolite)

        # Ensure first item in each list is showing.
        for metab_list in (self.ListMetabolitesAvailable, self.ListMetabolitesInUse):
            if metab_list.GetCount():
                metab_list.SetFirstItem(0)


    def populate_pulse_sequences(self):
        """Populates the list of available pulse sequences and selects
        the one indicated in the experiment (if any).
        """
        self.pulse_sequences = self.db.fetch_pulse_sequences()

        names = [pulse_sequence.name for pulse_sequence in self.pulse_sequences]

        self.ComboPulseSequence.AppendItems(names)

        if self.experiment.pulse_sequence:
            for i, pulse_sequence in enumerate(self.pulse_sequences):
                if self.experiment.pulse_sequence.id == pulse_sequence.id:
                    self.ComboPulseSequence.SetSelection(i)
                    break

        self.pulse_sequence_setup()


    def pulse_sequence_setup(self):
        """
        Sets up the panel appropriately based on which pulse sequence
        is selected

        """
        i = self.ComboPulseSequence.GetCurrentSelection()
        if i == wx.NOT_FOUND:
            # No selection
            loop_labels = [ ]
            parameters = [ ]
        else:
            loop_labels = self.pulse_sequences[i].loop_labels
            parameters = self.pulse_sequences[i].user_static_parameters

        for i, label in enumerate(loop_labels):
            control_set = self.loop_controls[i]
            control_set.label.SetLabel(label)
            control_set.start.ChangeValue("")
            control_set.count.ChangeValue("")
            control_set.size.ChangeValue("")
            control_set.start.Enable()
            control_set.count.Enable()
            control_set.size.Enable()

        for control_set in self.loop_controls[len(loop_labels):]:
            control_set.label.SetLabel("")
            control_set.start.ChangeValue("")
            control_set.count.ChangeValue("")
            control_set.size.ChangeValue("")
            control_set.start.Disable()
            control_set.count.Disable()
            control_set.size.Disable()

        sizer = self.ParameterSizer

        # There are three columns, one for the descriptive label, one for the
        # textbox that holds the default value, and one for another label that
        # describes the type (string, double, etc.)
        sizer.SetCols(3)
        sizer.SetRows(len(parameters))
        if len(parameters):
            sizer.AddGrowableCol(1, 1)

        # Get rid of any existing controls
        for control_group in self.parameter_controls:
            for control in control_group:
                control.Destroy()
        self.parameter_controls = [ ]

        # Create one row of new controls for each param
        for i, parameter in enumerate(parameters):
            controls = []
            name = parameter.name
            if not name.endswith(":"):
                name += ":"
            name_label = wx.StaticText(self, wx.ID_ANY, name)
            sizer.Add(name_label, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_RIGHT|wx.RIGHT, 5)

            textbox = wx.TextCtrl(self)
            textbox.ChangeValue(parameter.default)
            sizer.Add(textbox, 1, wx.LEFT|wx.RIGHT|wx.EXPAND, 5)

            type_label = wx.StaticText(self, wx.ID_ANY, "(%s)" % parameter.type)
            sizer.Add(type_label, 0, wx.ALIGN_CENTER_VERTICAL)
            sizer.Layout()

            self.parameter_controls.append( (name_label, textbox, type_label) )

        # wx sets tab order according to control creation order. Since I
        # just created controls, they'll be *after* the Run button in the
        # tab order which is wrong. Here I correct the tab order.
        if self.parameter_controls:
            last_control = self.parameter_controls[-1][-1]
        else:
            last_control = self.loop_controls[-1].size

        self.ButtonRun.MoveAfterInTabOrder(last_control)

        # Note that this function can get called during init but also when
        # pulse seq info changes during editing. In the latter case, I only
        # want to populate the pulse seq loop controls if the selected pulse
        # seq was already present in the experiment.
        i = self.ComboPulseSequence.GetCurrentSelection()
        if (i != wx.NOT_FOUND)             and                              \
           self.experiment.pulse_sequence  and                              \
           (self.experiment.pulse_sequence.id == self.pulse_sequences[i].id):
            for i, dim in enumerate(self.experiment.dims):
                control_set = self.loop_controls[i]
                if control_set.start.IsEnabled():
                    start, step, length = mrs_experiment.deconstruct_loop(dim)
                    control_set.start.ChangeValue(str(start))
                    control_set.size.ChangeValue(str(step))
                    control_set.count.ChangeValue(str(length))

        self.Layout()
        self.Refresh()


    def save_experiment(self):
        """Saves the experiment, if possible, and returns True. Also returns
        True if no save is needed.
        Returns False if the experiment can't be saved (because it needs to
        be run or contains invalid data).
        """
        msg = ""

        if self.is_new:
            # (self.is_new == True) ==> nothing important in self.last_save
            if self.last_run:
                # OK, the expt has been run at least once which is a necessary
                # precondition for saving.
                # Has it changed since it was last run?
                last_run = self.last_run
                gui = self.get_raw_gui_data()
                # When deciding whether or not this needs to be re-run,
                # name & comment changes are not significant so I force them
                # to always match.
                last_run["name"] = gui["name"]
                last_run["comment"] = gui["comment"]

                if last_run != gui:
                    # Significant changes have been made since the last run.
                    msg = "Please re-run this experiment before saving it."
            else:
                # Can't save an experiment that hasn't been run
                msg = "Please run this experiment before saving it."
        else:
            # This is an experiment that's been saved once already. Life is
            # simpler in this case because only incidentals or metabs can
            # be changed. In addition, the latter can only be added to.
            # So here we check to see if the last run metabs are the same as
            # the metabs currently in the GUI. The safe way to do this is
            # by comparing UUIDs (as opposed to comparing objects which may
            # be different instances of the Metabolite class that represent
            # the same metabolite).
            last_run = sorted([metabolite.id for metabolite
                                             in self.last_run["metabolites"]])
            raw = sorted([metabolite.id for metabolite
                                        in self.get_raw_gui_data()["metabolites"]])
            if last_run != raw:
                # Metabs have changed since the last run.
                msg = "Please re-run this experiment before saving it."

        if not msg:
            # so far so good
            msg = self.is_gui_valid(self.is_new)

        if msg:
            common_dialogs.message(msg, None, common_dialogs.I_OK)
            rc = False
        else:
            rc = True
            self.populate_experiment_from_gui()

            if self.is_new:
                self.db.insert_experiment(self.experiment)
            else:
                if self.db.check_experiment_id(self.experiment.id):
                    # This is already in the database
                    self.db.replace_experiment(self.experiment)
                else:
                    # Not in database AND not 'is_new' so maybe is a derived
                    # experiment that has been locked to prevent editing
                    self.db.insert_experiment(self.experiment)

            self.is_new = False
            self.last_save = self.get_raw_gui_data()

            self.set_control_enabled()


    def set_control_enabled(self):
        """
        Disables controls in the GUI as appropriate. Assumes that
        everything is currently enabled.

        """
        # For new experiments, evething is left enabled. For existing
        # experiments, the GUI is much more locked down.
        if not self.is_new:
            self.TextInvestigator.Disable()
            self.ComboIsotope.Disable()
            self.ComboIsotope.Disable()

            self.TextB0.Disable()
            self.TextPeakSearchRangeLow.Disable()
            self.TextPeakSearchRangeHigh.Disable()
            self.TextBlendTolerancePpm.Disable()
            self.TextBlendTolerancePhase.Disable()
            self.ComboPulseSequence.Disable()

            for i, dim in enumerate(self.experiment.dims):
                control_set = self.loop_controls[i]
                control_set.start.Disable()
                control_set.size.Disable()
                control_set.count.Disable()

            for control in self.parameter_controls:
                # user parameters saved in tuple = (name_label, textbox, type_label)
                control[1].Disable()

            self.ButtonMoveRight.Disable()

            self.ButtonRun.SetLabel("Run Additional Metabolite(s)")

            # if the Experiment has been "made public" ie. exported to others,
            # then it's frozen which means almost none of it can be changed.
            if self.experiment.is_public == True:
                self.ButtonMoveLeft.Disable()
                self.ButtonMoveRight.Disable()
                self.ButtonRun.Disable()

            self.Layout()


