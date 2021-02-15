# Python modules

import os

# 3rd party modules
import wx
import wx.html
#import wx.aui as aui
import wx.lib.agw.aui as aui        # NB. wx.aui version throws odd wxWidgets exception on Close/Exit ?? Not anymore in wxPython 4.0.6 ??
import numpy as np

# Our modules
import vespa.simulation.tab_experiment as tab_experiment
import vespa.simulation.dialog_calculate_add_sub as dialog_calculate_add_sub
import vespa.common.util.misc as util_misc
import vespa.common.wx_gravy.common_dialogs as common_dialogs
import vespa.common.wx_gravy.util as wx_util
import vespa.common.wx_gravy.notebooks as vespa_notebooks

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


class NotebookExperiments(vespa_notebooks.VespaAuiNotebook):
    """This is the notebook that contains all of the experiment tabs (or the
    welcome tab if there aren't any experiments open). The app instantiates
    exactly one instance of this class so it's used as a singleton, although
    the class itself doesn't enforce singleton-ness.
    """
    # I need the path to the welcome page image which is in vespa/common.
    _path = util_misc.get_vespa_install_directory()
    _path = os.path.join(_path, "common", "resources", "simulation_welcome.png")

    WELCOME_TAB_TEXT = """
    <html><body>
    <h1>Vespa - Simulation</h1>
    <img src="%s" alt="Time-Freq Plots" />
    <p><b>Currently there are no experiments loaded.</b></p>
    <p>Use the Experiment menu to Open an existing Experiment or 
    create a New experiment.</p>
    </body></html>
    """ % _path
    # I tidy up my namespace by deleting this temporary variable.
    del _path


    def __init__(self, top, db):
        vespa_notebooks.VespaAuiNotebook.__init__(self, top)

        self.db    = db
        self.top   = top
        self.count = 0
        
        # I make my get_open_experiments() function available to other parts
        # of the app.
        wx.GetApp().vespa.get_open_experiments = self.get_open_experiments

        self.show_welcome_tab()
        
        self.Bind(aui.EVT_AUINOTEBOOK_PAGE_CLOSE, self.on_tab_close)
        self.Bind(aui.EVT_AUINOTEBOOK_PAGE_CLOSED, self.on_tab_closed)
        self.Bind(aui.EVT_AUINOTEBOOK_PAGE_CHANGED, self.on_tab_changed)
        

    #####################    Event handlers start here 

    def on_menu_view_option(self, event):
        if self.active_tab:
            self.active_tab.on_menu_view_option(event)


    def on_menu_view_output(self, event):
        if self.active_tab:
            self.active_tab.on_menu_view_output(event)


    def on_tab_changed(self, event):
        # Clear the status bar
        statusbar = self.top.statusbar
        for i in range(statusbar.GetFieldsCount()):
            statusbar.SetStatusText("", i)

        name = ""
        
        # Create an appropriate name for whatever is selected.
        tab = self.active_tab
        if tab:
            name = tab.experiment.name
            if not name:
                name = "New Experiment"
                
            name = " - " + name
                
            tab.on_activation()

        wx.GetApp().GetTopWindow().SetTitle("Simulation" + name)
            

    def on_tab_close(self, event):
        veto = False

        if self.active_tab:
            if not self.active_tab.close():
                veto = True
        else:
            veto = True
            
        if veto:
            event.Veto()


    def on_tab_closed(self, event):
        if not self.GetPageCount():
            self.show_welcome_tab()
        
    ###############    Internal & shared helper methods start here 
    ###############    They're in alphabetical order

    def activate_experiment_tab(self, experiment_id):
        """Activates the tab containing the experiment identified by the 
        experiment_id param. Returns True if it found a tab matching the id,
        False otherwise.
        """
        activated = False
        
        app = wx.GetApp()
        
        for i, experiment in enumerate(self.get_open_experiments()):
            if experiment.id == experiment_id:
                activated = True
                self.SetSelection(i)
                break

        return activated
            

    def add_experiment_tab(self, experiment, is_new):
        """Creates a new experiment tab (which becomes the active tab) and 
        populates it with the experiment's data if appropriate.
        """        
        # If the welcome tab is open, close it.
        if self.is_welcome_tab_open:
            self.DeletePage(0) 
        
        # Create new notebook with Experiment controls and display
        self.count += 1
        etab = tab_experiment.TabExperiment(self, self.top, self.db,
                                            experiment, is_new)
       
        self.AddPage(etab, "Experiment%d" % self.count, True)
        if is_new:
            etab.SetSelection(0)
            etab.simulate.SetFocus()
            

    def copy_tab(self):
        if self.active_tab:
            msg = "The copy will be made from the last time this " \
                  "experiment was saved or run."
            common_dialogs.message(msg, "Copy Tab to New Experiment", 
                                   common_dialogs.I_OK)
            
            experiment = self.active_tab.experiment.clone()
            experiment.id = util_misc.uuid()
            # Add "copy" to the name. This is just a convenience for the 
            # user so that when they save a copied tab we usually won't have
            # to bug them to change the name first. It's not guaranteed to be
            # unique, but it will be most of the time and that's good enough.
            experiment.name += " copy"
            
            self.add_experiment_tab(experiment, True)
            

    def add_sub_tab(self):
        if self.active_tab:
            
            experiment = self.active_tab.experiment
            
            # check if there is a dimension that has at least 2 entries in it
            # to which to apply the add/sub operation
            
            if (len(experiment.dims[0]) < 2 and \
                len(experiment.dims[1]) < 2 and \
                len(experiment.dims[2]) < 2) :
                msg = "The active experiment contains no loops with " \
                      "dimensionality of two or more. Can not perform an " \
                      "add/subtract operation. Returning."
                common_dialogs.message(msg, "Error Calculating Add/Sub Experiment", 
                                       common_dialogs.I_OK)
                return
                
            dialog = dialog_calculate_add_sub.DialogCalculateAddSub(self, experiment)

            if dialog.ShowModal() == wx.ID_OK:

                idim = dialog.index
                            
                original   = self.active_tab.experiment
                experiment = self.active_tab.experiment.clone()
                experiment.id = util_misc.uuid()
                # Add "Add-Sub" to the name. This is just a convenience for the 
                # user so that when they save a copied tab we usually won't have
                # to bug them to change the name first. It's not guaranteed to be
                # unique, but it will be most of the time and that's good enough.
                experiment.name += " Add-Sub"
                 
                # do the add/sub calculations on the dimension selected len(experiment.simulations)
                shape = [len(experiment.metabolites)]
                shape += [len(dim) for dim in experiment.dims]
                                
                esim = np.array(experiment.simulations)
                osim = np.array(original.simulations)
                esim.shape = shape[::-1]
                osim.shape = shape[::-1]
                
                if dialog.index == 0:
                    for imet in range(shape[0]):
                        for i in range(shape[2]):
                            for j in range(shape[3]):
                                esim[j,i,0,imet].add_simulation(     osim[j,i,1,imet])
                                esim[j,i,1,imet].subtract_simulation(osim[j,i,0,imet])
                elif dialog.index == 1:
                    for imet in range(shape[0]):
                        for i in range(shape[1]):
                            for j in range(shape[3]):
                                esim[j,0,i,imet].add_simulation(     osim[j,1,i,imet])
                                esim[j,1,i,imet].subtract_simulation(osim[j,0,i,imet])
                elif dialog.index == 2:
                    for imet in range(shape[0]):
                        for i in range(shape[1]):
                            for j in range(shape[3]):
                                esim[0,j,i,imet].add_simulation(     osim[1,j,i,imet])
                                esim[1,j,i,imet].subtract_simulation(osim[0,j,i,imet])
                
                # Add tab to the notebook with the Add-Sub results in it. 
                # We set 'is_public' to True and do not add this as a 'new' 
                # experiment to force the Simulate tab to gray out the Run
                # button and not allow users to add any metabolite. This is
                # because it would not make sense to add more On/Off results
                # to Add/Sub values.
                experiment.is_public = True
                experiment.comment += '\n'+dialog.comment
                self.add_experiment_tab(experiment, False)    
            #else:
                # User hit cancel, don't calculate a new tab.
        


    def close_all_experiment_tabs(self):
        """Attemps to close all open experiment tabs. If a tab contains 
        unsaved changes, the user will be prompted and can cancel the close.
        
        Returns True if all tabs were closed or False if any were cancelled.
        """
        close_accepted = True
        # While there are experiments open and the user hasn't cancelled,
        # close the current tab.
        i = len(self.get_open_experiments())
        while i and close_accepted:
            self.close_experiment_tab()
            # There should be one less experiment tab now. If not, the
            # user cancelled the close.
            i -= 1
            close_accepted = (i == len(self.get_open_experiments()))
            
        return close_accepted
        

    def close_experiment_tab(self):
        """Attempts to close the active experiment tab. If the tab contains 
        unsaved changes, the user will be prompted and can cancel the close.
        """
        if self.active_tab:
            wx_util.send_close_to_active_tab(self)

        # It would be nice to return True/False indicating whether or not the
        # user cancelled the close, but that's difficult
        # to do. The event sent above simulates a click which sets off a 
        # cascade of other events, including an EVT_AUINOTEBOOK_PAGE_CLOSE
        # which is what we trap to interrupt the close. It's nice that we
        # don't have to simulate the cascade "by hand", but it also means 
        # we're disconnected from the outcome.


    def get_open_experiments(self):
        return [tab.experiment for tab in self.tabs]


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
        
