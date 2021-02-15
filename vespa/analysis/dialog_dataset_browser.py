# Python modules

# 3rd party modules
import wx

# Our modules
import vespa.analysis.mrs_dataset as mrs_dataset
import vespa.analysis.auto_gui.dataset_browser as dataset_browser
import vespa.common.wx_gravy.util as wx_util

from pathlib import Path

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


class DialogDatasetBrowser(dataset_browser.MyDialog):
    """
    Opens a dialog that allows the user to select a dataset from all of the
    datasets open in tabs in Analysis.

    After the dialog closes, the dialog's dataset attribute contains the 
    selected dataset, or None if the user hits cancel.
    
    The datasets param passed to __init__() should be a dict of dataset
    names (tab names) mapping to datasets.

    """
    def __init__(self, datasets, parent=None, omit_self=None):
        if not parent:
            parent = wx.GetApp().GetTopWindow()

        dataset_browser.MyDialog.__init__(self, parent)

        _sort_datasets = [(name,datasets[name]) for name in sorted(datasets.keys())]
        
        if omit_self and isinstance(omit_self, mrs_dataset.Dataset):
            _all_datasets  = [item for item in _sort_datasets if item[1].id != omit_self.id]
        else:
            _all_datasets = _sort_datasets
        
        # self._all_datasets is where I stash the dict of all available 
        # datasets from which the user can choose.
        # self.dataset contains the dataset that the user choose and
        # is only relevant as this dialog closes.
        self._all_datasets = _all_datasets
        self.dataset = None

        # Add Open & Cancel buttons dynamically so in right order under OS X, GTK, Windows
        self.ButtonOK, _ = wx_util.add_ok_cancel(self, self.LabelOkCancelPlaceholder, self.on_ok)
        self._populate_list()

        wx_util.set_font_to_monospace(self.Listbox)
        
        self.SetSize( (700, 350) )
        self.Layout()
        self.Center()
        self.Listbox.SetFocus()


    ########################    Event Handlers    #################

    def on_list_double_click(self, event):
        self.on_ok()
        

    def on_list_selection(self, event):
        # Enable OK if something selected.
        self.ButtonOK.Enable( (len(self.Listbox.GetSelections()) > 0) )


    def on_ok(self, event=None):
        # This can be called programmatically, in which case event is None.
        i = self.Listbox.GetSelection()
        if i != wx.NOT_FOUND:
            # Just in case OK clicked without something selected
            self.dataset = self.Listbox.GetClientData(i)
        self.Close()


    ########################    "Private" Methods    #################
        
    def _populate_list(self):
        # FIXME PS - the dataset sorting code doesn't work as intended. The
        # problem is that "Dataset10" sorts *before* "Dataset9". That can
        # be fixed but alphabetical probably isn't the best sort order
        # anyway. A better choice would be to mimic the order of the tabs
        # as they're open in the notebook, but that will require some
        # help from NotebookDatasets.

        # Each list item is [tab name]: [data filename] - [data source] where the
        # data source is the first listed in the dataset. When constructing the
        # listbox entries, I right justify each tab and filename so that they line
        # up regardless of how many characters are in the name, e.g. --
        #     Dataset1:     bob.xml - blah blah blah
        #    Dataset52:  thomas.xml - blah blah blah
        #   Dataset999:   chuck.xml - blah blah blah

        wx.SetCursor(wx.HOURGLASS_CURSOR)

        # Reset controls
        self.ButtonOK.Enable(False)
        self.Listbox.Clear()

        tags = []
        fnames = []
        sources = []
        dsets = []

        ntags = 0
        nfnames = 0

        for item in self._all_datasets:
            source = item[1].blocks["raw"].data_source
            fname = Path(source).name
            tags.append(item[0])
            fnames.append(fname)
            sources.append(source)
            dsets.append(item[1])
            if ntags < len(item[0]) + 2: ntags = len(item[0]) + 2
            if nfnames < len(fname) + 2: nfnames = len(fname)

        for tag, fname, source, ds in zip(tags,fnames,sources,dsets):
            s = "%s: %s - %s" % (tag.rjust(ntags), fname.rjust(nfnames), source)
            self.Listbox.Append(s)
            self.Listbox.SetClientData((self.Listbox.GetCount() - 1), ds)

        wx.SetCursor(wx.NullCursor)

########################################################################

class MyForm(wx.Frame):

    def __init__(self, datasets):

        self.datasets = datasets
        wx.Frame.__init__(self, None, wx.ID_ANY, "Testing the Dialog")
        panel = wx.Panel(self, wx.ID_ANY)
        buttn = wx.Button(panel, label="Show Dialog")
        buttn.Bind(wx.EVT_BUTTON, self.on_dialog)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(buttn, 0, wx.ALL | wx.CENTER, 5)
        panel.SetSizer(sizer)

    def on_dialog(self, event):

        dialog = DialogDatasetBrowser(self.datasets, parent=self)
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

    ds1 = mrs_dataset.Dataset()
    ds2 = mrs_dataset.Dataset()
    ds3 = mrs_dataset.Dataset()

    ds1.blocks["raw"].data_sources = ["C:/users/bsoher/chuck1.py",]
    ds2.blocks["raw"].data_sources = ["D:/data/home/bob2.py",]
    ds3.blocks["raw"].data_sources = ["C:/temp/dreck/thomasina3.py",]

    datasets = {'bob1':ds1, 'bob2':ds2, 'bob3':ds3}

    app = wx.App(False)
    frame = MyForm(datasets)
    frame.Show()
    app.MainLoop()