# Python modules

# 3rd party modules
import wx

# Our modules
import vespa.analysis.tab_base as tab_base
import vespa.analysis.util_menu as util_menu
import vespa.analysis.auto_gui.raw as raw
import vespa.common.wx_gravy.util as wx_util
import vespa.common.util.misc as util_misc


#------------------------------------------------------------------------------
#
#  Tab RAW 
#
#------------------------------------------------------------------------------

class TabRaw(tab_base.Tab, raw.PanelRawUI):
    
    # self-identify tab to notebook, value does not matter, its presence is sufficient.
    IS_RAW = True
    
    def __init__(self, tab_dataset, top, block):       
        
        raw.PanelRawUI.__init__(self, tab_dataset.NotebookDataset)
        
        tab_base.Tab.__init__(self, tab_dataset, top, None)
        
        self.top    = top               # application frame
        self.block  = block              # processing object
        
        self.initialize_controls()
        self.populate_controls()

        wx_util.set_font_to_monospace(self.TextHeaderInformationDisplay)
        
        self.TextHeaderInformationDisplay.Bind(wx.EVT_LEFT_DCLICK, 
                                               self.on_header_double_click)


    #=======================================================
    #
    #           GUI Setup Handlers 
    #
    #=======================================================

    def initialize_controls(self):
        """ 
        Initializes the controls to be the right size or have the right
        range or number of decimal places. It typically does not set the
        default value (that's for populate_controls method to do). This
        method does the one-time setup bits.
        """
        # Under OS X, the call to SetSashPosition() doesn't work unless
        # we use wx.CallAfter().
        wx.CallAfter(self.WindowSplitter.SetSashPosition, 400, True)

 
    def populate_controls(self, preset=False):
        """ 
        Populates the raw data tab with values from the dataset.raw 
        object. It's meant to be called when a new data object is loaded.

        """
        raw = self.dataset.blocks["raw"]

        self.ListDataSources.Clear()
        self.ListDataSources.SetItems(raw.data_sources)

        self.TextHeaderInformationDisplay.Clear()
        
        if self.ListDataSources.GetCount():
            # Select the first item. Since programmatic selection doesn't
            # trigger an event, we have to call the event handler manually.
            self.ListDataSources.SetSelection(0)
            self.on_list_item_select()
        

    
    #=======================================================
    #
    #           Menu Event Handlers 
    #
    #=======================================================

    def on_menu_view_output(self, event):
        event_id = event.GetId()
        
        if event_id == util_menu.ViewIdsRaw.OUTPUT_HEADER_TEXT:
            self._display_header_text()


    #=======================================================
    #
    #           wx Event Handlers  
    #
    #=======================================================
    
    def on_header_double_click(self, event):
        # If the user double clicks on the displayed header, we open it in
        # the default text editor so thatthe user can copy, print, etc.
        header = self.TextHeaderInformationDisplay.GetValue()
        wx_util.display_text_as_file(header)


    def on_list_double_click(self, event):
        # Display the header for the selected item in a text editor so that
        # the user can copy, print, etc.
        index = self.ListDataSources.GetSelection()
        if index != wx.NOT_FOUND:
            header = self.dataset.blocks["raw"].headers[index]
            wx_util.display_text_as_file(header)


    def on_list_item_select(self, event=None):
        # Display the header for the selected item
        # Since this is sometimes pretty cryptic, I also precede it with
        # a printout of the raw data block parameters set from the header.
        # Sometimes this is called from this code (i.e. by us, not wx) and
        # in that case the event object isn't present.
        header0 = str(self.dataset.blocks["raw"])
        
        index = self.ListDataSources.GetSelection()
        if index != wx.NOT_FOUND:
            header = self.dataset.blocks["raw"].headers[index]
            header = util_misc.normalize_newlines(header)
        else:
            header = ""
        
        hdr = header0 + '\n\n----------- Original Header Data -----------\n' + header
        
        self.TextHeaderInformationDisplay.SetValue(hdr)




    #=======================================================
    #
    #           Public Methods
    #
    #=======================================================
    def on_voxel_change(self, voxel):
        index = voxel[0]
        if index < self.ListDataSources.GetCount() - 1:
            # Select the item in the list that corresponds to the selected voxel
            self.ListDataSources.SetSelection(index)
            self.on_list_item_select()


    #
    # Note. the process(), plot() and process_and_plot() methods are are
    #       inherited from the base class and are basically no-ops
    #
        

    #=======================================================
    #
    #           Internal Helper Functions  
    #
    #=======================================================

    def _display_header_text(self):
    
        lines = "\nCurrent Header\n" + "-" * 75 + "\n\n"
        lines += str(self.TextHeaderInformationDisplay.GetValue())
            
        wx_util.display_text_as_file(lines)



