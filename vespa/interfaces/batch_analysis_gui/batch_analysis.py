# Python modules

import argparse
import os
import sys
import platform
import collections


# 3rd party modules
import wx
import wx.html

import vespa.analysis.auto_gui.gui_batch as gui_batch
import vespa.common.wx_gravy.common_dialogs as common_dialogs
import vespa.common.util.time_ as util_time


# Our modules
import vespa.interfaces.batch_analysis_gui.analysis_cli as cli



def _configure_combo(control, choices, selection=''):
        lines = [choices[key] for key in list(choices.keys())]
        control.Clear()
        control.AppendItems( lines )
        if selection in lines:
            control.SetStringSelection(selection)
        else:
            control.SetStringSelection(lines[0])


class DataTypes(object):
    """ Supported data types for batch processing """
    # These constants are arbitrary and may change.
    WBNAA = 'wbnaa'

    # Items for the data type combo pull-down
    choices = collections.OrderedDict(( (WBNAA , "WBNAA"),  ))


class RedirectText(object):
    """
    This is used to redirect std.out and if wanted, std.err text output to
    a WxTextCtrl window in the widget.  Sort of slick.

    """
    def __init__(self,aWxTextCtrl):
        self.out=aWxTextCtrl

    def write(self,string):
        self.out.WriteText(string)


class BatchAnalysis(gui_batch.GuiBatch):

    def __init__(self, parent, id=wx.ID_ANY):
        gui_batch.GuiBatch.__init__(self, parent, id)

        self.data_type  = 'wbnaa'
        self.top_level_directory = ''
        self.extension  = '.dat'
        self.presetfile = ''
        self.csvfile    = ''
        self.savexml    = False
        self.twodir     = False
        self.verbose    = False
        self.debug      = False

        _configure_combo(self.ComboDataType, DataTypes.choices)

        # redirect text here
        redir = RedirectText(self.TextConsole)
        sys.stdout = redir


    def on_data_type(self, event):
        index = self.ComboDataType.GetSelection()
        val = list(DataTypes.choices.keys())[index]
        self.data_type = val

    def on_extension(self, event):
        val = event.GetEventObject().GetValue()
        # ensure extension starts with '.' and no extra whitespace
        val = val.strip()
        if val[0] != '.':
            val = '.'+val
        event.GetEventObject().SetValue(val)
        self.extension = val

    def on_browse_top_level_directory(self, event):
        dirname = common_dialogs.pickdir()
        if dirname:
            self.top_level_directory = dirname
            self.StaticTopLevelDirectory.SetLabel(dirname)

    def on_browse_preset_file(self, event):
        filename = common_dialogs.pickfile()
        if filename:
            self.presetfile = filename
            self.StaticPresetFile.SetLabel(filename)

    def on_browse_csv_file(self, event):
        filename = common_dialogs.pickfile(new_file=True)
        if filename:
            self.csvfile = filename
            self.StaticCsvFile.SetLabel(filename)

    def on_save_xml(self, event):
        self.savexml = event.GetEventObject().GetValue()

    def on_two_directory(self, event):
        self.twodir = event.GetEventObject().GetValue()

    def on_verbose(self, event):
        self.verbose = event.GetEventObject().GetValue()

    def on_debug(self, event):
        self.debug = event.GetEventObject().GetValue()

    def on_quit(self, event):
        self.Destroy()

    def on_do_batch(self, event):

        msg = ''
        if not self.top_level_directory:
            msg = "Error - Select a top level directory for MRS files."
        if not self.presetfile:
            msg = "Error - Select a preset file for processing."
        if not self.csvfile:
            msg = "Error - Select a CSV output file."

        if msg:
            self.say(msg + "\n")
            self.say("Batch run finished with errors.\n")
            return

        # this gets all files *.dat in all subdirectories

        datafiles = []
        for dirpath, dirnames, filenames in os.walk(self.top_level_directory):
            for filename in [f for f in filenames if f.endswith(self.extension)]:
                datafiles.append(os.path.join(dirpath, filename))

        if self.debug:
            msg = "Debug: Listing data file names only - no processing."
            self.say(msg + "\n")
            for datafile in datafiles:
                self.say(datafile + "\n")
            return

        msg = "Start Processing - "
        self.say(msg + "\n")
        presetfile = self.presetfile
        csvfile    = self.csvfile

        for datafile in datafiles:
            cli.analysis_cli(datafile,
                             self.data_type,
                             self.presetfile,
                             self.csvfile,
                             savexml=self.savexml,
                             twodir=self.twodir,
                             verbose=self.verbose,
                             debug=self.debug)

        msg = "End Processing - "
        self.say(msg + "\n")

    def on_help(self, event):
        page = self.get_help_page()
        helpDlg = HelpDlg(self, page)
        helpDlg.Show()

    def get_help_page(self):
        if self.data_type == 'wbnaa':
            return page_wbnaa
        else:
            return page_null

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
        # We force an update in case the next step is a long one
        self.TextConsole.Update()


class HelpDlg(wx.Frame):

    def __init__(self, parent, page, title="Batch Analysis - Help", size=(600,400)):

        wx.Frame.__init__(self, parent, wx.ID_ANY, title=title, size=size)
        html = wx.html.HtmlWindow(self)
        html.SetPage( page )


page_wbnaa = r"""
<HTML>
<BODY>

<h2>Help - Fitting WBNAA</h2>

<h3>Overview</h3>
<p>Fill in values for all the widgets from top to bottom. After all values are set, click on the Do Batch button. If you do not want to perform the fit, or just want to exit the program, hit the Quit button. </p>

<h3>Type of Data to be Fitted - Pull Down</h3>
<p>At the moment the only selection in this pull down is for fitting WBNAA data from Siemens VB data. </p>

<h3>Top Level Directory to search for Data Files</h3>
<p>Use the Browse button to select a directory that contains all data file in one or more subdirectories. You can type in an extension string to further restrict the files returned. You have to hit return after typing in the Extension string. </p>

<h3>Analysis Preset File</h3>
<p>Use the Browse button to select an XML file containing the settings to be used to fit the WBNAA data. Note. this file will be used to fit all files. </p>

<h3>CSV Output File</h3>
<p>Use the Browse button to select a file into which to output fitting results. If the file does not exist it will be created. If it exists, the results will be appended to that file. Results are saved as comma separated text column output. </p>

<h3>Flag Settings</h3>
<p><b>Save Data to XML</b> - Fitted Analysis dataset with provenance will be saved to and XML file at the level of the original data file. A unique filename will be created from concatenating the original file name and one or two parent directory names. </p>
<p><b>Use Two Level Naming</b> - Used by the Save to XML command to determine if one or two parent directory names are used to create output filename. </p>
<p><b>Verbose</b> - Progress messages will be printed to the Console Output window. </p>
<p><b>Debug</b> - If this flag is set, NO PROCESSING will be performed. Only a list of data file names will be printed to the Console Output window. Use this to check if all data files are being found by your search directory/extension settings. </p>

<h3>Quit and Do Batch Buttons</h3>
<p><b>Quit</b> - select to leave the program without performing fitting or when fitting is done. </p>
<p><b>Do Batch</b> - select to perform fitting. Program is not closed when finished. </p>

<h3>Console Output - Text Window</h3>
<p>Information about the processing or errors encountered are printed to this window. It is read only. </p>

</BODY>
</HTML>
"""

page_null = r"""
<HTML>
<BODY>

<h2>Help - No Page Found</h2>

<p>No help page found for the data type selected. </p>

</BODY>
</HTML>
"""


#--------------------------------------------------------------------

if __name__ == "__main__":
    app = wx.App(0)
    app.SetAppName("Batch Analysis")
    frame1 = BatchAnalysis(None, -1)
    frame1.Show()
    app.MainLoop()

