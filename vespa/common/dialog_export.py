# Python modules


# 3rd party modules
import wx
import os

# Our modules
import vespa.common.auto_gui.export as gui_export
import vespa.common.wx_gravy.common_dialogs as common_dialogs
import vespa.common.wx_gravy.util as wx_util
import vespa.common.util.config as util_config

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


class DialogExport(gui_export.MyDialog):
    """Displays a dialog that gathers the information needed for an
    export (filename, comment & compression on/off). When the dialog is 
    closed, the user's input is recorded in the following variables:
        dialog.export   - True if the user clicked Export, False if he
                          clicked Cancel.
        dialog.comment  - Comment text (string, may be blank).
        dialog.filename - Selected filename, fully qualified.
        dialog.compress - True if compression selected.
    
    show_warning controls whether or not this dialog displays the warning 
    about exported objects becoming public & frozen. Analysis sets this to 
    False since it doesn't store objects in the database.
    """
    def __init__(self, parent=None, show_warning=True, filename=None):
        if not parent:
            parent = wx.GetApp().GetTopWindow()

        gui_export.MyDialog.__init__(self, parent)
        
        # These variables contain the user's selections.
        self.export = False
        self.comment = ""
        self.filename = ""
        self.compress = False
        
        # We add the Export & Cancel buttons dynamically so that they're in  
        # the right order under OS X, GTK, Windows, etc.
        wx_util.add_ok_cancel(self, self.LabelButtonPlaceholder,
                              self.on_export, ok_text="E&xport")
        
        # Set the suggested filename
        if filename != None:
            self.default_filename = filename
        else:
            self.default_filename = "vespa_%s_export.xml" % wx.GetApp().GetAppName()
        self.default_filename = self.default_filename.lower()        
        path = util_config.get_last_export_path()
        self.LabelFilename.SetLabel(os.path.join(path, self.default_filename))
        
        if show_warning:
            font = self.LabelWarning.GetFont()
            font.SetWeight(wx.FONTWEIGHT_BOLD)
            self.LabelWarning.SetFont(font)
        else:
            self.PanelWarning.Hide()
        
        # It's a good idea to make this dialog pretty wide so that the 
        # entire filename is visible even if it's long.
        self.SetSize( (700, 350) )
        
        self.Layout()

        if show_warning:
            wx_util.wrap_label(self.LabelWarning, self)
        
        self.Center()

        
    def on_browse(self, event):
        # Show the file picker dialog
        filename = self.default_filename
                
        default_path = util_config.get_last_export_path()
        
        if self.CheckboxCompress.GetValue():
            filename += ".gz"
            filetype_filter = "Vespa Export Files (*.gz)|*.gz"
        else:
            filetype_filter = "Vespa Export Files (*.xml)|*.xml"
        
        filename = common_dialogs.save_as(filetype_filter=filetype_filter,
                                          default_path=default_path,
                                          default_filename=filename)
                                          
        if filename:
            self.LabelFilename.SetLabel(filename)
            path, _ = os.path.split(filename)
            util_config.set_last_export_path(path)


    def on_check_compress(self, event):
        # We try to keep the filename in sync with the compression setting.
        # This isn't necessary since the importer doesn't take its clues
        # from the filename, but it hints to the user to choose a sensible
        # name.
        filename = self.LabelFilename.GetLabel()
                
        if filename:
            if self.CheckboxCompress.GetValue():
                if not filename.endswith(".gz"):
                    self.LabelFilename.SetLabel(filename + ".gz")
            else:
                if filename.endswith(".gz"):
                    self.LabelFilename.SetLabel(filename[:-3])

            
    def on_export(self, event):
        self.filename = self.LabelFilename.GetLabel()
        self.comment = self.TextComment.GetValue()
        self.compress = self.CheckboxCompress.GetValue()

        path, _ = os.path.split(self.filename)
        
        if not os.path.exists(path) or not os.path.isdir(path):
            msg = """The path "%s" doesn't exist or is not a directory.""" % path
            common_dialogs.message(msg, "Vespa Export", common_dialogs.E_OK)
        else:
            self.export = True
            self.Close()

