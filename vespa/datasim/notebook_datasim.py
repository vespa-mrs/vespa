# Python modules

import os

# 3rd party modules
import wx
import wx.html
#import wx.aui as aui
import wx.lib.agw.aui as aui        # NB. wx.aui version throws odd wxWidgets exception on Close/Exit ?? Not anymore in wxPython 4.0.6 ??

# Our modules
import vespa.datasim.tab_datasim as tab_datasim
import vespa.common.wx_gravy.util as wx_util
import vespa.common.wx_gravy.notebooks as vespa_notebooks
import vespa.common.util.misc as util_misc



class NotebookDatasim(vespa_notebooks.VespaAuiNotebook):
    # I need the path to the welcome tab image which is in vespa/common.
    _path = util_misc.get_vespa_install_directory()
    _path = os.path.join(_path, "common", "resources", "prior_welcome.png")

    WELCOME_TAB_TEXT = """
    <html><body>
    <h1>Vespa - DataSim</h1>
    <img src="%s" alt="Time-Freq Plots" />
    <p><b>Currently there are no datasims loaded.</b></p>
    <p>Use the Datasim menu to Open a saved DataSim or browse for a Vespa-Simulation to create a New datasim.</p>
    </body></html>
    """ % _path
    # I tidy up my namespace by deleting this temporary variable.
    del _path
    
    
    def __init__(self, top):

        vespa_notebooks.VespaAuiNotebook.__init__(self, top)

        self.top    = top
        self.count  = 0
        
        self.show_welcome_tab()

        self.Bind(aui.EVT_AUINOTEBOOK_PAGE_CLOSE, self.on_tab_close)
        self.Bind(aui.EVT_AUINOTEBOOK_PAGE_CLOSED, self.on_tab_closed)
        self.Bind(aui.EVT_AUINOTEBOOK_PAGE_CHANGED, self.on_tab_changed)
        

    #=======================================================
    #
    #           Global and Menu Event Handlers 
    #
    #=======================================================

    def on_menu_view_option(self, event):
        if self.active_tab:
            self.active_tab.on_menu_view_option(event)

    def on_menu_view_output(self, event):
        if self.active_tab:
            self.active_tab.on_menu_view_output(event)

    def on_tab_changed(self, event):

        self._set_title()
            
        if self.active_tab:
            self.active_tab.on_activation()
            
            
    def on_tab_close(self, event):
        """
        This is a two step event. Here we give the user a chance to cancel 
        the Close action. If user selects to continue, then the on_tab_closed()
        event will also fire.  
        
        """
        msg = "Are you sure you want to close this datasim?"
        if wx.MessageBox(msg, "Close Datasim", wx.YES_NO, self) != wx.YES:
            event.Veto()



    def on_tab_closed(self, event):        
        """
        At this point the tab is already closed and the datasim removed from
        memory.        
        """
        if not self.tabs:
            self.show_welcome_tab()

        self._set_title()


    #=======================================================
    #
    #           Public methods shown below
    #             in alphabetical order 
    #
    #=======================================================

    def add_datasim_tab(self, datasim=None):

        # If the welcome tab is open, close it.
        if self.is_welcome_tab_open:
            self.remove_tab(index=0)

        self.count += 1
        name = "Datasim%d" % self.count

        # create new notebook tab with process controls 
        tab = tab_datasim.TabDatasim(self, self.top, datasim)
        self.AddPage(tab, name, True)


    def close_datasim(self):
        if self.active_tab:
            wx_util.send_close_to_active_tab(self)


    #=======================================================
    #
    #           Internal methods shown below
    #             in alphabetical order 
    #
    #=======================================================

    def _set_title(self):
        title = "Datasim"

        if self.active_tab:
            tab = self.active_tab

            if tab.datasim:
                title += " - " + tab.datasim.experiment.name

        wx.GetApp().GetTopWindow().SetTitle(title)



    ###############    Overloaded event handlers re. wx.lib.aui.auibook.py 
    ###############    
    
    
    def OnTabButton(self, event):
        """
        Handles the ``EVT_AUINOTEBOOK_BUTTON`` event for :class:`AuiNotebook`.

        :param `event`: a :class:`AuiNotebookEvent` event to be processed.
        
        This is overloaded here because the original calls Destroy() through
        a CallAfter() method. This messes up a number of code steps downstream
        that need to know how many tabs remain in the notebook (e.g. putting
        the StartTab back in place if the last Experiment tab is closed).
        
        """

        tabs = event.GetEventObject()
        button_id = event.GetInt()

        if button_id == aui.AUI_BUTTON_CLOSE:

            selection = event.GetSelection()

            if selection == -1:

                # if the close button is to the right, use the active
                # page selection to determine which page to close
                selection = tabs.GetActivePage()

            if selection == -1 or not tabs.GetEnabled(selection):
                return

            if selection != -1:

                close_wnd = tabs.GetWindowFromIdx(selection)

                if close_wnd.GetName() == "__fake__page__":
                    # This is a notebook preview
                    previous_active, page_status = close_wnd.__previousStatus
                    for page, status in zip(tabs.GetPages(), page_status):
                        page.enabled = status

                    main_idx = self._tabs.GetIdxFromWindow(close_wnd)
                    self.DeletePage(main_idx)

                    if previous_active >= 0:
                        tabs.SetActivePage(previous_active)
                        page_count = tabs.GetPageCount()
                        selection = -1

                        for page in range(page_count):
                            # remove the page from the source tabs
                            page_info = tabs.GetPage(page)
                            if page_info.active:
                                selection = page
                                break

                        tabs.DoShowHide()
                        self.DoSizing()
                        tabs.Refresh()

                        if selection >= 0:
                            wx.CallAfter(tabs.MakeTabVisible, selection, self)

                    # Don't fire the event
                    return

                # ask owner if it's ok to close the tab
                e = aui.AuiNotebookEvent(aui.wxEVT_COMMAND_AUINOTEBOOK_PAGE_CLOSE, self.GetId())
                idx = self._tabs.GetIdxFromWindow(close_wnd)
                e.SetSelection(idx)
                e.SetOldSelection(event.GetSelection())
                e.SetEventObject(self)
                self.GetEventHandler().ProcessEvent(e)
                if not e.IsAllowed():
                    return

                if repr(close_wnd.__class__).find("AuiMDIChildFrame") >= 0:
                    close_wnd.Close()

                else:
                    main_idx = self._tabs.GetIdxFromWindow(close_wnd)
                    self.DeletePage(main_idx)
                    #wx.CallAfter(self.DeletePage, main_idx)

                # notify owner that the tab has been closed
                e2 = aui.AuiNotebookEvent(aui.wxEVT_COMMAND_AUINOTEBOOK_PAGE_CLOSED, self.GetId())
                e2.SetSelection(idx)
                e2.SetEventObject(self)
                self.GetEventHandler().ProcessEvent(e2)

                if self.GetPageCount() == 0:
                    mgr = self.GetAuiManager()
                    win = mgr.GetManagedWindow()
                    win.SendSizeEvent()    
