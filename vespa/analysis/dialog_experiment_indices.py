# Python modules

# 3rd party modules
import wx

# Our modules
import vespa.common.mrs_experiment as mrs_experiment
import vespa.common.wx_gravy.util as wx_util
import vespa.common.wx_gravy.common_dialogs as common_dialogs
import vespa.analysis.auto_gui.dialog_experiment_indices as dialog_experiment_indices

# Number of places right of decimal when displaying a float number
_SIGNIFICANT_DIGITS = 6


def _find_best_column_size(listctrl, column):
    # ListCtrls can be sized according to the header or the largest content
    # item. Sometimes sizing by the header means content gets clipped, and
    # sometimes the reverse happens. This function figures out which of the
    # two is larger and returns that.
    listctrl.SetColumnWidth(column, wx.LIST_AUTOSIZE)
    value_width = listctrl.GetColumnWidth(column)

    listctrl.SetColumnWidth(column, wx.LIST_AUTOSIZE_USEHEADER)
    header_width = listctrl.GetColumnWidth(column)

    return max(value_width, header_width)


class DialogExperimentIndices(dialog_experiment_indices.MyDialog):
    """
    Displays a dialog box that contains a list of the dimensions in the
    experiment. The user is required to pick one or hit cancel.

    When the dialog closes, dialog.selected_dim_id is set to
    the id of the selected dimension, or 0 if the user didn't choose a dim.
    
    The dim_ids parameter to this dialog should be a dict keyed by dim id
    mapping to the dim tuple, e.g. { 42: {0.1, 3.4, 0.0) }
    """
    def __init__(self, parent, db, experiment_id, dim_ids):
        dialog_experiment_indices.MyDialog.__init__(self, parent)

        self.db = db

        # Find the preview for this experiment so we can get the pulse seq.
        preview = db.fetch_experiment_preview(experiment_id)

        pulse_sequence = db.fetch_pulse_sequence(preview.pulse_sequence_id)

        # All we need from the pulse seq is loop labels.
        loop_labels = pulse_sequence.loop_labels

        self.selected_dim_id = 0

        # Build a dict keyed by the dimension tuple (e.g. (0, 1.3, 2.9) )
        # that maps dimensions to the unique id assigned by the database.
        # We'll need this to map the tuple selected by the user to the value
        # that this dialog returns.
        self.dim_ids = dict((value, key) for key, value in dim_ids.items())

        # Reduce the dims to a list of tuples containing the unique values
        # for each dim.
        self.dims = list(self.dim_ids.keys())
        self.dims = list(zip(*self.dims))
        # Filter out duplicates & sort what's left
        self.dims = [sorted(list(set(dim))) for dim in self.dims]
        
        # Make a shorthand reference for the labels above the lists
        self.heading_labels = [getattr(self, "LabelLoop%d" % i) for i 
                                                                in range(1, 4)] 

        self._initialize_controls(loop_labels)

        # We add the Open & Cancel buttons dynamically so that they're in the
        # right order under OS X, GTK, Windows, etc.
        self.ButtonOpen, self.ButtonCancel = \
            wx_util.add_ok_cancel(self, self.LabelOKCancelPlaceholder,
                                  self.on_ok)

        width, height = self.GetSize()
        self.SetSize( (width, 300) )
        self.Layout()
        self.Center()

        wx.CallAfter(self._autosize_listctrls)

        # Under Linux, the focus is elsewhere (?) unless I manually
        # set it to the list
        self.ListLoop1.SetFocus()

        self.Bind(wx.EVT_SIZE,  self.on_size)


    #################     Event Handlers   ##################
    
    def on_list_select(self, event):
        # This is called when the user makes a selection in any of the lists.
        # Figure out which listctrl was selected
        for i in range(1, 4):
            listctrl = getattr(self, "ListLoop%d" % i)
            if event.GetId() == listctrl.GetId():
                break
                
        # Update the heading to reflect what was selected
        self._set_list_heading(listctrl, i - 1)
    
    
    def on_size(self, event):
        # Correct label wrapping if necessary.
        self._wrap_instructions()


    def on_ok(self, event):
        msg = ""

        dims = [ ]

        if not msg:
            for i in range(len(self.dims)):
                listctrl = getattr(self, "ListLoop%d" % (i + 1))
                if listctrl.IsShown():
                    selections = wx_util.get_selected_item_indices(listctrl)
                    if not selections:
                        msg = "Please select a value for loop %d." % (i + 1)
                        break
                    else:
                        dims.append(self.dims[i][selections[0]])
                else:
                    # This is an empty dim
                    dims.append(0)

        if msg:
            common_dialogs.message(msg)
        else:
            self.selected_dim_id = self.dim_ids[tuple(dims)]

            self.Close()



    #################     Internal Helpers    ##################

    def _autosize_listctrls(self):
        # Here's some good old fashioned fun! There will be 1-3 lists showing,
        # and the optimal size for each list depends on its content. (Optimal
        # size means the minimum width required to avoid the need for a
        # horizontal scrollbar.) Since the lists are in a grid sizer that 
        # gives equal size to each list based on the size of the dialog, 
        # there's a lot of inter-related factors at work in determining the
        # size of each lists.
        # By the time this code executes, the columns inside the lists have
        # been sized to fit their content. Therefore we know that if the
        # width of the two columns plus a constant (here called FUDGE_FACTOR)
        # to account for the vertical scrollbar, mistakes, etc. is greater
        # than the width of the list control, then the list is not optimally
        # sized (wx will add a horizontal scrollbar). So, if a list isn't
        # optimally sized, we increase the width of the dialog by a fixed 
        # percentage, allow the dialog to resize itself which in turn resizes
        # (widens) the list, and then we measure again. 
        # We're done when no list control needs to be resized.
        FUDGE_FACTOR = 25
        resized = False

        for i, dim in enumerate(self.dims):
            listctrl = getattr(self, "ListLoop%d" % (i + 1))

            if listctrl.IsShown():
                column_widths = listctrl.GetColumnWidth(0) + listctrl.GetColumnWidth(1)
                listctrl_width, _ = listctrl.GetSize()
                if (column_widths + FUDGE_FACTOR) > listctrl_width:
                    # List is not optimally sized. Increase dialog size by 10%.
                    width, height = self.GetSize()
                    self.SetSize( (int(width * 1.1), height) )
                    self.Layout()

                    resized = True
                    break
            # else:
                # We don't bother with hidden list controls

        if resized:
            # The dialog was resized. Allow the resize to occur and then 
            # call this method again. 
            wx.CallAfter(self._autosize_listctrls)
        else:
            # All the lists are optimally sized -- we're done!
            self._wrap_instructions()

            self.Center()


    def _initialize_controls(self, loop_labels):
        for i, dim in enumerate(self.dims):
            listctrl = getattr(self, "ListLoop%d" % (i + 1))
            label = getattr(self, "LabelLoop%d" % (i + 1))

            if dim == mrs_experiment.DEFAULT_LOOP:
                # This is an empty loop.
                listctrl.Hide()
                label.Hide()
            else:
                # Build the ListCtrl columns, set label & font
                listctrl.InsertColumn(0, "Index", wx.LIST_FORMAT_RIGHT)
                listctrl.InsertColumn(1, loop_labels[i], wx.LIST_FORMAT_RIGHT)
                label.SetLabel("Loop %d" % (i + 1))
                
                self._set_list_heading(listctrl, i)

                # We use monospace font so that the padding we use (spaces)
                # and the numbers will line up correctly.
                wx_util.set_font_to_monospace(listctrl)

                # Figure out the width (in digits) of the max index value and
                # the max dim. This allows us to right justify the numbers
                # in the lists which makes them much easier to read.
                index_digits = len(str(len(dim)))
                value_digits = len("%d" % max(dim))

                is_int = all([int(value) == value for value in dim])

                if is_int:
                    # All dim values are ints
                    formatter = "%d"
                else:
                    # Some values are floats. Format them as such and account
                    # for the width of the digits to the right of the decimal.
                    formatter = "%." + str(_SIGNIFICANT_DIGITS) + "f"
                    value_digits += _SIGNIFICANT_DIGITS

                # Add some padding
                formatter = "   " + formatter

                # Populate the list
                for i, value in enumerate(dim):
                    listctrl.InsertItem(i, str(i + 1).rjust(index_digits))

                    listctrl.SetItem(i, 1, (formatter % value).rjust(value_digits))

                # Size the columns optimally
                listctrl.SetColumnWidth(0, _find_best_column_size(listctrl, 0))
                listctrl.SetColumnWidth(1, _find_best_column_size(listctrl, 1))


    def _set_list_heading(self, listctrl, i):
        # Figure out which item (if any) is selected and update the heading
        # to contain that info. We do this because the background color for
        # a selected item is very hard to see on Windows 7, so it's hard for
        # users to see what they've selected. This is only a problem under
        # Win7. It's fine under Win XP, OS X, and Linux/GTK.
        selections = wx_util.get_selected_item_indices(listctrl)
        label = "Loop %d" % (i + 1)
        if selections:
            label += ( " = %s" % listctrl.GetItem(selections[0], 1).GetText().strip())
        
        self.heading_labels[i].SetLabel(label)
        
        # This recenters the text
        self.heading_labels[i].GetContainingSizer().Layout()
        

    def _wrap_instructions(self):
        wx_util.wrap_label(self.LabelInstructions, self)

