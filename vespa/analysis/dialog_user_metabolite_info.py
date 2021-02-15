# Python modules
import os

# 3rd party modules
import wx
import wx.grid as gridlib

# Our modules
import vespa.analysis.util_import as util_import
import vespa.analysis.util_analysis_config as util_analysis_config
import vespa.common.wx_gravy.util as wx_util

from vespa.common.wx_gravy.common_dialogs import pickfile, message, E_OK
from vespa.analysis.auto_gui.dialog_user_metabolite_info import MyDialog


#------------------------------------------------------------------------------
# Note. GUI Architecture/Style
#
# Vespa GUI components are designed using WxGlade. WxGlade files (with *.wxg
# extensions) are stored in the 'wxglade' subdirectory under each application.
# Ouput of their code generation are stored in the 'auto_gui' subdirectory.
#
# Each GUI Class is inherited into a unique 'vespa' Class, where program specific
# initialization and other functionality are written.  Original stub functions
# for event handlers are overloaded to provide program specific event handling.
#------------------------------------------------------------------------------


class DialogUserMetaboliteInfo(MyDialog):
    """
    Used to organize user defined information about metabolites. Things like
    alternative metabolite name spellings and literature concentration values.

    """
    def __init__(self, parent, dataset):

        if not parent:
            parent = wx.GetApp().GetTopWindow()

        MyDialog.__init__(self, parent)

        self._dataset  = dataset
        self.parent    = parent
        self.metinfo   = self._dataset.user_prior.metinfo
        self.new_index = 0
        self.limits    = None

        self._initialize_controls()

        # Add OK & Cancel buttons dynamically for proper order under OSX, GTK, Windows
        r = wx_util.add_ok_cancel(self, self.LabelOkCancelPlaceholder, self.on_ok)
        self.ButtonOk, self.ButtonCancel = r

        self.Layout()
        self.Fit()
        self.SetSize( (580, 700) )

        self._improve_height()

        # Under Linux, focus is elsewhere (?) unless manually set to the Grid
        self.ButtonAddMetabolite.SetFocus()


    ##### Event Handlers ######################################################

    def on_reset(self, event):
        self.GridData.DeleteCols(pos=0,numCols=6)
        self._initialize_controls()
        self.Layout()
        self.Fit()
        self.SetSize((580, 700))
        self._improve_height()

    def on_load_from_preset(self, event):

        default_path = util_analysis_config.get_path("load_metinfo")
        filename = pickfile(filetype_filter="Spectra Preset (*.xml,*.xml.gz)|*.xml;*.xml.gz",
                            multiple=False, default_path=default_path)
        if filename:
            msg = ""
            try:
                importer = util_import.MetinfoImporter(filename)
            except IOError:
                msg = """I can't read the Metinfo XML file "%s".""" % filename
            except SyntaxError:
                msg = """The Metinfo file "%s" isn't valid Vespa Interchange File Format.""" % filename

            if msg:
                message(msg, "Dialog User Metabolite Info - Open Metinfo File", E_OK)
            else:
                wx.BeginBusyCursor()
                results = importer.go()        # Time to rock and roll!
                wx.EndBusyCursor()

                metinfo = results[0]

                self.GridData.DeleteCols(pos=0, numCols=self.GridData.GetNumberCols())
                self.GridData.DeleteRows(pos=0, numRows=self.GridData.GetNumberRows())
                self._initialize_controls(metinfo=metinfo)
                self.Layout()
                self.Fit()
                self.SetSize((580, 700))
                self._improve_height()

                path, _ = os.path.split(filename)
                util_analysis_config.set_path("load_metinfo", path)

    def on_select_all(self, event):
        for i in range(self.GridData.GetNumberRows()):
            self.GridData.SetCellValue(i, 0, '1')
        self.GridData.ForceRefresh()

    def on_deselect_all(self, event):
        for i in range(self.GridData.GetNumberRows()):
            self.GridData.SetCellValue(i, 0, '')
        self.GridData.ForceRefresh()

    def on_add_metabolite(self, event):
        self.new_index += 1
        if self.GridData.AppendRows():
            nrow = self.GridData.GetNumberRows() - 1
            self.GridData.SetCellValue(nrow, 0, '')
            self.GridData.SetCellValue(nrow, 1, 'new_entry'+str(self.new_index))
            self.GridData.SetCellAlignment(wx.ALIGN_CENTER, nrow, 1)
            self.GridData.SetCellValue(nrow, 2, 'new_abbr'+str(self.new_index))
            self.GridData.SetCellAlignment(wx.ALIGN_CENTER, nrow, 2)
            self.GridData.SetCellValue(nrow, 3, '1')
            self.GridData.SetCellAlignment(wx.ALIGN_CENTER, nrow, 3)
            self.GridData.SetCellValue(nrow, 4, '1.0')
            self.GridData.SetCellAlignment(wx.ALIGN_CENTER, nrow, 4)
            self.GridData.SetCellValue(nrow, 5, '250.0')
            self.GridData.SetCellAlignment(wx.ALIGN_CENTER, nrow, 5)
        wx.CallAfter(self._improve_height)

    def on_remove_selected(self, event):
        nrow = self.GridData.GetNumberRows()
        self.GridData.BeginBatch()
        remove = [i for i in range(nrow) if self.GridData.GetCellValue(i, 0) == '1']
        for i in remove[::-1]:
            self.GridData.DeleteRows(i, 1)
        self.GridData.EndBatch()
        wx.CallAfter(self._improve_height)
        self.Layout()
        self.Refresh()

    def on_ok(self,event):
        # No changes made until user selects "OK" rather than "Cancel" so
        # here we update the metinfo object to reflect values in dialog

        # TODO, bjs, likely we need to do some error checking here for the
        # string widgets to ensure that there are no blank lines

        names, abbrs, spins, concs, t2 = self.get_values()

        self.metinfo.full_names     = names
        self.metinfo.abbreviations  = abbrs
        self.metinfo.spins          = spins
        self.metinfo.concentrations = concs
        self.metinfo.t2decays       = t2

        self.EndModal(True)

    def on_cancel(self, event):
        # Nothing to do here since changes are only written on "OK"
        self.EndModal(False)


    ##### Internal helper functions  ##########################################

    def _improve_height(self):
        # No matter what format the user chooses, we use a scrolled window to
        # contain the list of metabs because with a fixed-sized (non-scrolled)
        # panel, experiments that contain a lot of metabs (20+) can cause
        # this dialog to exceed the display height. However, the scrolled
        # window has its own problem: it doesn't size itself intelligently.
        # It always has the same (small) default height.
        #
        # So here we check to see if the scrolled window is actually scrolled;
        # i.e. it contains content that it can't completely display in its
        # current area. If so, we increase the dialog's height exactly
        # enough to allow the scrolled window to expand so that the user
        # doesn't have to scroll to see the contents.
        _, display_height = wx.GetDisplaySize()
        dialog_width, dialog_height = self.GetSize()

        # Compare virtual height with real height.
        # delta is how much bigger it needs to be to display all of its
        # content without scrolling.
        _, v_height = self.ScrolledWindowGridData.GetVirtualSize()
        _, r_height = self.ScrolledWindowGridData.GetClientSize()
        delta = v_height - r_height

        # max_delta is the max we can increase the dialog height before it
        # exceeds the display area. Note that wx reports the raw display
        # area without accounting for things like the Windows taskbar or the
        # OS X dock. Actual space available for use by applications may be
        # less. To account for this we subtract a fudge factor of 132 pixels.
        # This is pretty arbitrary, although on my Mac laptop the OS X dock
        # and the top menu occupy 132 pixels and the dock is huge relative
        # to the Windows taskbar and Gnome's similar thingies, so
        # hopefully this will be sufficient everywhere.
        max_delta = (display_height - dialog_height) - 132

        delta = min(delta, max_delta)

        if delta > 0:
            self.SetSize( (dialog_width, dialog_height + delta) )

        self.Center()

    def _initialize_controls(self, metinfo=None):

        if metinfo is None:
            metinfo = self.metinfo

        #----------------------------------------------------------------------
        # sets up other widgets in metabolite info dialog from values in the
        # metinfo object sent into the dialog initialization

        headings = ['','Metabolite\nFull Name', 'Abbreviation', 'Number\nof Spins',
                    'Literature\nConc [mM]', 'Literature\nT2 decay [ms]']

        self.limits = [[],[],[],[1,20], [0.01,50000.0], [0.01,1000.0]]

        self.GridData.ClearGrid()
        self.GridData.AppendCols(6)
        self.GridData.AppendRows(len(metinfo.full_names))

        self.GridData.SetColFormatBool(0)
        self.GridData.SetColFormatNumber(3)     # Cols 1,2 strings by default
        self.GridData.SetColFormatFloat(4, precision=2)
        self.GridData.SetColFormatFloat(5, precision=2)

        for col, labl in enumerate(headings):
            self.GridData.SetColLabelValue(col, labl)

        self.GridData.BeginBatch()
        for row, line in enumerate(metinfo.to_lines()):
            self.GridData.SetCellValue( row, 0, '')
            self.GridData.SetCellEditor(row, 0, gridlib.GridCellBoolEditor())
            self.GridData.SetCellRenderer(row, 0, gridlib.GridCellBoolRenderer())
            for col, item in enumerate(line):
                if col+1 in [1,2]:
                    self.GridData.SetCellValue(row, col + 1, str(item))
                elif col+1==3:
                    self.GridData.SetCellValue( row, col + 1, str(int(item)))
                    self.GridData.SetCellEditor(row, col + 1, gridlib.GridCellNumberEditor(min=self.limits[3][0], max=self.limits[3][1]))
                elif col+1 in [4,5]:
                    self.GridData.SetCellValue( row, col + 1, str(float(item)))
                    self.GridData.SetCellEditor(row, col + 1, gridlib.GridCellFloatEditor())
                self.GridData.SetCellAlignment(wx.ALIGN_CENTER, row, col + 1)
        self.GridData.EndBatch()

        self.GridData.SetColLabelAlignment(wx.ALIGN_CENTER, wx.ALIGN_CENTER)
        self.GridData.HideRowLabels()
        self.GridData.SetMargins(0, 0)
        self.GridData.AutoSizeColumns(True)

        self.Bind(gridlib.EVT_GRID_CELL_CHANGED, self.on_cell_change)

    def on_cell_change(self, evt):

        lim = self.limits
        row, col = evt.GetRow(), evt.GetCol()
        value = self.GridData.GetCellValue(row,col)

        #print("on_cell_change: (%d,%d) %s\n" % (row, col, evt.GetPosition()))
        # Show how to stay in a cell that has bad data.  We can't just
        # call SetGridCursor here since we are nested inside one so it
        # won't have any effect.  Instead, set coordinates to move to in
        # idle time.
        #print("on_cell_change - Value = %s" % value)

        if col in [3,4,5]:
            if float(value) < lim[col][0]:
                self.GridData.SetCellValue(row, col, str(lim[col][0]))
                self.GridData.SetCellTextColour(row, col, wx.RED)
            elif float(value) > lim[col][1]:
                self.GridData.SetCellValue(row, col, str(lim[col][1]))
                self.GridData.SetCellTextColour(row, col, wx.RED)
            else:
                self.GridData.SetCellTextColour(row, col, wx.BLACK)

    def get_values(self):

        grid = self.GridData
        nrow = grid.GetNumberRows()
        names = [grid.GetCellValue(i,1) for i in range(nrow)]
        abbrs = [grid.GetCellValue(i,2) for i in range(nrow)]
        spins = [int(grid.GetCellValue(i,3)) for i in range(nrow)]
        concs = [float(grid.GetCellValue(i,4)) for i in range(nrow)]
        t2s   = [float(grid.GetCellValue(i,5)) for i in range(nrow)]

        return names, abbrs, spins, concs, t2s


########################################################################

class MyForm(wx.Frame):

    def __init__(self, dataset):

        self.dataset = dataset
        wx.Frame.__init__(self, None, wx.ID_ANY, "Testing the Dialog")
        panel = wx.Panel(self, wx.ID_ANY)
        buttn = wx.Button(panel, label="Show Dialog")
        buttn.Bind(wx.EVT_BUTTON, self.on_dialog)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(buttn, 0, wx.ALL | wx.CENTER, 5)
        panel.SetSizer(sizer)

    def on_dialog(self, event):

        dialog = DialogUserMetaboliteInfo(self, self.dataset)
        code = dialog.ShowModal()
        if code == wx.ID_OK:
            print("DialogUserMetaboliteInfo  User hit OK")
        elif code == wx.ID_CANCEL:
            print("DialogUserMetaboliteInfo  User Cancelled")
        dialog.Destroy()


# ----------------------------------------------------------------------
# Run the program

if __name__ == "__main__":

    import vespa.analysis.mrs_dataset as mrs_dataset

    app = wx.App(False)
    frame = MyForm(mrs_dataset.Dataset())
    frame.Show()
    app.MainLoop()

