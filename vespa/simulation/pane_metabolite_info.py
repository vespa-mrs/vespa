# Python modules


# 3rd party modules
import wx

# Our modules
import vespa.simulation.auto_gui.metabolite_info as gui_metabolite_info
import vespa.simulation.dialog_experiment_list as dialog_experiment_list
import vespa.common.mrs_metabolite as mrs_metabolite
import vespa.common.mrs_spin as mrs_spin
import vespa.common.mrs_j_coupling as mrs_j_coupling
import vespa.common.util.time_ as util_time
import vespa.common.util.misc as util_misc
import vespa.common.constants as common_constants
import vespa.common.wx_gravy.common_dialogs as common_dialogs
import vespa.common.wx_gravy.util as wx_util

# max spins of 14 is a somewhat arbitrary number, mostly just based on
# limits of the real world. The dialog is coded so that you can change this
# number and the dialog will adjust. That's the intent, anyway!
MAX_SPINS = 14

# This tell me how many J coupling values there are for a given number
# of spins. e.g with MAX_SPINS = 14 --
# {  1:  0,  2:  1,  3:  3,  4:  6, 5: 10, 6: 15, 7: 21, 8: 28, 9: 36, 10: 45,
#   11: 55, 12: 66, 13: 78, 14: 91}
J_COUPLINGS_PER_SPIN = dict([(i, sum(range(i))) for i in range(1, MAX_SPINS + 1)])

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


class PaneMetaboliteInfo(gui_metabolite_info.MetaboliteInfoUI):

    def __init__(self, parent, db, metabolite=None, read_only=False):
        gui_metabolite_info.MetaboliteInfoUI.__init__(self, parent)

        # This db handle is only used to (a) get a list of isotopes and
        # (b) check that a metab name is unique. It's not for writing!
        self.db = db

        # This dialog has three states of "editabililty" --
        # - When creating a new metab, all fields are editable
        # - When editing a frozen metab, only name and comments are editable.
        #   A metab is frozen if its is_frozen flag is set OR if it's
        #   used in an experiment.
        # - When viewing an existing metab, no controls are editable. In
        #   this case, the caller will set the read_only param = True.
        
        if metabolite:
            self.read_only = read_only
            self.metabolite = metabolite
        else:
            # This is a new metab
            # When creating a new metab, read_only can't be True.
            self.read_only = False
            # Create a new metab object
            self.metabolite = mrs_metabolite.Metabolite()
            # I give brand new metab a spin because the init
            # code is simpler if I can assume that all metabs have at least
            # one spin.
            # Note that a metab with just one spin has zero J couplings so
            # this is a logically complete minimalist metabolite.
            spin = mrs_spin.Spin()
            spin.isotope = common_constants.DEFAULT_ISOTOPE
            spin.chemical_shift = 0.0
            self.metabolite.spins = [ spin ]
            # This non-standard attribute must exist for all metabs
            # in this dialog. The "List Experiments" button expects it.
            self.metabolite.experiment_names = [ ]
            
        self.metabolite_in_use = bool(self.metabolite.experiment_names)

        # These are references to the controls that I'll create in
        # initialize_controls().
        self.shift_textboxes = [ ]
        self.isotope_comboboxes = [ ]
        self.j_textboxes = [ ]

        self.initialize_controls(db.fetch_isotopes())
        
        self.InitDialog()


    # EVENT HANDLER Redefinitions ----------------------------------------
    def on_spin_count_change(self, event):
        spin_count = int(self.ComboNumberOfSpins.GetStringSelection())

        # All of the textboxes with indices <= the spin count are enabled,
        # the others are disabled and emptied.
        for control in self.shift_textboxes[:spin_count]:
            if not control.IsEnabled():
                control.Enable()
                control.ChangeValue(str(0.0))
            #else:
                # It's already enabled; leave it alone.

        for control in self.isotope_comboboxes[:spin_count]:
            if not control.IsEnabled():
                control.Enable()
                if not control.GetValue():
                    control.SetStringSelection(common_constants.DEFAULT_ISOTOPE)
            #else:
                # It's already enabled; leave it alone.

        for control in (self.shift_textboxes[spin_count:] +
                        self.isotope_comboboxes[spin_count:]):
            control.Disable()
            control.SetValue("")

        # En/disable j coupling textboxes
        values_per_row = list(range(spin_count - 1, 0, -1))
        values_per_row += [ 0 ] * (MAX_SPINS - len(values_per_row))
        
        for row, textboxes in enumerate(self.j_textboxes):
            for column, textbox in enumerate(textboxes):
                if column < values_per_row[row]:
                    if not textbox.IsEnabled():
                        textbox.Enable()
                        textbox.ChangeValue(str(0.0))
                    #else:
                        # It's already enabled; leave it alone.
                else:
                    textbox.Disable()
                    textbox.ChangeValue("")
                    

    def on_ok(self, event):
        # Validate controls one by one
        propogate_event = True
        msg = ""

        # name
        name = self.TextName.GetValue().strip()        
        if not name:
            msg = "Please enter a name for the metabolite."

        # creator
        if not msg:
            creator = self.TextCreator.GetValue().strip()        

        # comment
        if not msg:
            # Note that I don't call .strip() on this field. It's totally
            # freeform and whitespace might be important.
            comment = self.TextComment.GetValue()
            
        # spins
        if not msg:
            for i, textbox in enumerate(self.shift_textboxes):
                if textbox.IsEnabled():
                    shift = textbox.GetValue().strip()
                    if shift:
                        if not util_misc.is_floatable(shift):
                            msg = """I don't understand the shift "%s".""" % shift
                    else:
                        msg = "Please enter a shift value for spin %d." % (i + 1)
                        
                if msg:
                    # Bail out of the loop; no point in checking further.
                    break

        # J couplings
        if not msg:
            # Note: the call to sum() here flattens my list of lists.
            j_textboxes = sum(self.j_textboxes, [ ])
            for textbox in j_textboxes:
                if textbox.IsEnabled():
                    j_coupling = textbox.GetValue().strip()
                    if j_coupling:
                        if not util_misc.is_floatable(j_coupling):
                            msg = """I don't understand the J coupling "%s".""" % j_coupling
                    else:
                        msg = "Please enter a value for each J coupling."
                        
                if msg:
                    # Bail out of the loop; no point in checking further.
                    break
                    
        if not msg:
            # I still have yet to check whether or not this name is unique.
            # I save it until last because it requires a database hit.
            name = self.TextName.GetValue().strip()
            metabolites = self.db.fetch_metabolites_by_name(name)
            ids = [metabolite.id for metabolite in metabolites \
                                     if metabolite.id != self.metabolite.id]

            if ids:
                msg = "A metabolite with this name already exists."
                
        if msg:
            common_dialogs.message(msg, None, common_dialogs.X_OK)
            propogate_event = False
        else:
            # All is well
            self.metabolite.name = name    
            self.metabolite.creator = creator
            self.metabolite.comment = comment

            if self.metabolite.is_frozen:
                # Nothing to do; these controls can't be edited
                pass
            else:
                spin_count = int(self.ComboNumberOfSpins.GetStringSelection())
                self.metabolite.spins = \
                                [mrs_spin.Spin() for i in range(spin_count)]
                            
                for spin, textbox, combobox in \
                            zip(self.metabolite.spins, self.shift_textboxes, 
                                self.isotope_comboboxes):
                    spin.chemical_shift = float(textbox.GetValue())
                    spin.isotope = combobox.GetStringSelection()
    
                # As I (re)create the J couplings, I need to associate them
                # with the appropriate spins.
                # Note: the call to sum() here flattens my list of lists.
                j_textboxes = sum(self.j_textboxes, [ ])
                j_textboxes = [textbox for textbox in j_textboxes \
                                                       if textbox.IsEnabled()]

                self.metabolite.j_couplings = [ ]
                j_coupling_count = J_COUPLINGS_PER_SPIN[spin_count]
                spin1_index = 0
                spin2_index = 1                
                for i in range(j_coupling_count):
                    j_coupling = mrs_j_coupling.JCoupling()
                    j_coupling.value = float(j_textboxes[i].GetValue())
                    j_coupling.spin1 = self.metabolite.spins[spin1_index]
                    j_coupling.spin2 = self.metabolite.spins[spin2_index]
                
                    self.metabolite.j_couplings.append(j_coupling)
                
                    spin2_index += 1
                    if spin2_index == spin_count:
                        spin1_index += 1
                        spin2_index = spin1_index + 1
            #else:
                # Nothing more to do in this case; the rest of the controls
                # on the dialog are disabled.
     

        # Since we're in a pane, we have to call event.Skip() to
        # get the OK event passed to the dialog.
        if propogate_event:
            event.Skip()


    def on_experiments(self, event):
        dialog = dialog_experiment_list.DialogExperimentList(self, self.metabolite)
        dialog.ShowModal()



    #########################       Internal Helper functions

    def initialize_controls(self, all_isotopes):
        # Lots of work done here. Controls are disabled and enabled as
        # necessary, created, hidden, populated, etc.
        
        # Note that when self.read_only is True, everything is disabled 
        # except for the "Done" button. When metab.is_frozen is true, only
        # some controls are disabled.
        
        self.TextName.ChangeValue(self.metabolite.name)
        self.LabelUuid.SetLabel(self.metabolite.id)
        self.TextCreator.ChangeValue(self.metabolite.creator)
        if not self.metabolite.created:
            # Set to current time
            self.metabolite.created = util_time.now()
        created = self.metabolite.created.strftime(util_time.DISPLAY_DATE_FORMAT)
        self.LabelCreated.SetLabel(created)

        self.TextComment.ChangeValue(self.metabolite.comment)
                
        # populate ComboNumberOfSpins and select the appropriate item
        self.ComboNumberOfSpins.SetItems(["%d" % i for i in range(1, MAX_SPINS + 1)])
        self.ComboNumberOfSpins.SetSelection(len(self.metabolite.spins) - 1)
        self.ComboNumberOfSpins.Enable(not self.read_only)

        # Disable controls as appropriate.
        controls = [ ]
        if self.metabolite.is_frozen:
            controls = (self.ComboNumberOfSpins, self.TextCreator)
            
        if self.read_only:
            controls = (self.TextName, self.TextCreator, self.TextComment,
                        self.ComboNumberOfSpins)
            
        for control in controls:
            control.Disable()
                
        # Next, build the list of shift/isotope text/combobox pairs. These
        # live in a flex grid control that I find at runtime by marking it
        # with a placeholder.
        # Populating the grid sizer at runtime makes it easier to ensure
        # that all of the controls get consistent styles, sizes, etc. It
        # also makes it easy for me to get a handle to each of the textboxes.
        sizer = self.LabelSpinListPlaceholder.GetContainingSizer()
        self.LabelSpinListPlaceholder.Destroy()

        # The columns are a numeric label, the shift textbox and the isotope
        # textbox.
        sizer.SetCols(3)
        # There's one row for the column headers plus one row for each spin.
        sizer.SetRows(MAX_SPINS + 1)

        # Upper left grid corner is just empty space
        sizer.AddSpacer( 1 )
        # These are the column headers
        sizer.Add(wx.StaticText(self, label="Shift"))
        sizer.Add(wx.StaticText(self, label="Isotope"))

        # Build the grid contents
        label_style = wx.ALIGN_CENTER_VERTICAL | wx.RIGHT | wx.ALIGN_RIGHT
        label_border = 5

        # A width of 65 pixels is big enough to hold 000.0000 (or one fewer
        # diigts plus a negative sign) on OS X, GTK and Win XP.
        textbox_width = 65
        textbox_size = wx.Size(textbox_width, -1)
        for i in range(MAX_SPINS):
            sizer.Add(wx.StaticText(self, label=str(i + 1)), 0,
                                    label_style, label_border)

            textbox = wx.TextCtrl(self, size=textbox_size)
            sizer.Add(textbox)
            # Disable by default and enable only if necessary
            textbox.Disable()
            self.shift_textboxes.append(textbox)

            combobox = wx.ComboBox(self, choices=all_isotopes,
                                    size=textbox_size, style=wx.CB_DROPDOWN)

            sizer.Add(combobox)
            # Disable by default and enable only if necessary
            combobox.Disable()
            self.isotope_comboboxes.append(combobox)
        
        # Populate & enable the controls I just created.
        for spin, textbox, combobox in \
                            zip(self.metabolite.spins, self.shift_textboxes, 
                                self.isotope_comboboxes):
            textbox.ChangeValue(str(spin.chemical_shift))
            combobox.SetStringSelection(spin.isotope)
            if self.metabolite.is_frozen or self.read_only:
                # Leave these disabled
                pass
            else:
                textbox.Enable()
                combobox.Enable()

        sizer.Layout()

        # Here I build the J coupling textboxes using a technique similar to
        # the one I used for the shift & isotope controls.
        containing_sizer = self.LabelJCouplingsPlaceholder.GetContainingSizer()
        self.LabelJCouplingsPlaceholder.Destroy()
        
        # There're MAX_SPINS-1 data rows/columns plus one row/column for the
        # headers.
        sizer = wx.GridSizer(MAX_SPINS, MAX_SPINS, 2, 2)
        containing_sizer.Add(sizer)
        
        # Getting the label text centered under OS X and WinXP is just a
        # matter of creating the label the same size as the textbox and
        # applying the wx.ALIGN_CENTRE style. That doesn't work under GTK.
        # Instead I allow the label to size itself and tell its containing
        # sizer to center it horizontally.
        # This trick works under OSX, WinXP & GTK but relies (AFAICT) on the
        # undocumented fact that each grid cell is a unique sizer. I'd rather
        # not rely on that but that's the only way I've found centering will
        # work under GTK.
        label_style = wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL

        # The first thing I create in the grid spacer is the upper left corner
        # which is just empty space. As of wxPython 3.0.0 this is just set to 
        # a width of 1 and it seems to work ok - bjs
        sizer.AddSpacer(1)
        
        # Column headers
        for column in range(1, MAX_SPINS):
            sizer.Add(wx.StaticText(self, wx.ID_ANY, "%d" % (column + 1)),
                      flag=label_style)
        
        # Textboxes, spacers and labels.
        for row in range(1, MAX_SPINS):
            self.j_textboxes.append([ ])
            for column in range(1, MAX_SPINS + 1):
                if column < row:
                    # spacer to indent this row
                    sizer.AddSpacer( 1 )
                    
                if column == row:
                    # row label
                    sizer.Add(wx.StaticText(self, wx.ID_ANY, "%d" % row),
                              flag=label_style)
                              
                if column > row:
                    # textbox
                    textbox = wx.TextCtrl(self, size=textbox_size)
                    # I disable the textbox by default and enable it later
                    # if necessary.
                    textbox.Disable()
                    self.j_textboxes[-1].append(textbox)
                    sizer.Add(textbox)

        sizer.Layout()

        # Now I populate the J coupling textboxes. 
        # Basically, I insert values in textboxes until I run out of space in
        # that row, then I move to the next row. The number of elements
        # per row is calculated per row. e.g. with 6 spins: [5, 4, 3, 2, 1]
        elements_per_row = list(range(len(self.metabolite.spins) - 1, 0, -1))

        row = 0
        column = 0
        for i, j_coupling in enumerate(self.metabolite.j_couplings):
            textbox = self.j_textboxes[row][column]
            if self.metabolite.is_frozen or self.read_only:
                # Leave this disabled
                pass
            else:
                textbox.Enable()
            textbox.ChangeValue(str(j_coupling.value))

            # If my index has passed an end-of-row marker, it's time to move
            # to the next row.
            if (i + 1) == sum(elements_per_row[:row + 1]):
                # end of a row
                row += 1
                column = 0
            else:
                column += 1
        
        if not self.metabolite.experiment_names:
            # No point showing this button if it will show an empty list
            self.ButtonListExperiments.Destroy()

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

        

class DialogMetaboliteInfo(wx.Dialog):

    def __init__(self, parent, db, metabolite=None, read_only=False):
        if not parent:
            parent = wx.GetApp().GetTopWindow()

        if metabolite:
            title = "Metabolite " + metabolite.name
            
            states = [ ]
            if metabolite.is_frozen:
                states.append("frozen")
            else:
                states.append("not frozen")
            if metabolite.experiment_names:
                states.append("in use")
            else:
                states.append("not in use")
            if metabolite.is_public:
                states.append("public")
            else:
                states.append("private")
            if not metabolite.deactivated:
                states.append("active")
            if metabolite.deactivated:
                states.append("not active")
                
            if states:
                title += " (" + ", ".join(states) + ")"
        else:
            title = "New Metabolite"

        wx.Dialog.__init__(self, parent, wx.ID_ANY, title)

        self.pane = PaneMetaboliteInfo(self, db, metabolite, read_only)

        # Layout with sizers
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.pane, 0, wx.EXPAND|wx.ALL, 10)
        self.SetSizer(sizer)
        sizer.Fit(self)

        self.Center()

