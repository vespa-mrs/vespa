#!/usr/bin/env python

# Python modules

import os
import webbrowser
import struct

# 3rd party modules
import wx
#import wx.aui as aui
import wx.adv as wx_adv
import wx.lib.agw.aui as aui        # NB. wx.aui version throws odd wxWidgets exception on Close/Exit  ?? Not anymore in wxPython 4.0.6 ??


# Our modules
# We import vespa.common.wx_gravy.plot_panel here even though we don't use it
# in this module. It needs to be the first Vespa import to avoid matplotlib 
# sensitivity to import order. 
# See http://scion.duhs.duke.edu/vespa/project/ticket/26
import vespa.common.wx_gravy.plot_panel
import vespa.pulse.util_db as util_db
import vespa.pulse.util_menu as util_menu
import vespa.pulse.notebook_pulse_designs as notebook_pulse_designs
import vespa.pulse.dialog_third_party_export as dialog_third_party_export
import vespa.pulse.dialog_manage_pulse_designs as dialog_manage_pulse_designs
import vespa.pulse.dialog_manage_transform_kernels as dialog_manage_transform_kernels
import vespa.pulse.dialog_manage_machine_specs as dialog_manage_machine_specs
import vespa.pulse.util_pulse_config as rf_config
import vespa.common.images as images
import vespa.common.util.init as util_init
import vespa.common.util.misc as util_misc
import vespa.common.wx_gravy.common_dialogs as common_dialogs
import vespa.common.wx_gravy.util as wx_util
import vespa.common.dialog_pulse_design_browser as dialog_pulse_design_browser


#----------------------------------------------------------------------


class Main(wx.Frame):

    def __init__(self, db, position, size):
        
        self._left, self._top = position
        self._width, self._height = size
        
        style = wx.CAPTION | wx.CLOSE_BOX | wx.MINIMIZE_BOX |           \
                wx.MAXIMIZE_BOX | wx.SYSTEM_MENU | wx.RESIZE_BORDER |   \
                wx.CLIP_CHILDREN

        wx.Frame.__init__(self, None, wx.ID_ANY, "Pulse",
                          (self._left, self._top),
                          (self._width, self._height), style)

        self.db = db
        self.pulse_design_notebook = None
        
        # GUI Creation ----------------------------------------------
        
        aui_manager = aui.AuiManager()
        aui_manager.SetManagedWindow(self)
        # We don't use this reference to the AUI manager outside of this 
        # code, but if we don't keep a ref to it then Python will garbage
        # collect it and the AUI manager will die which usually results in 
        # a segfault. So here we create a reference that we don't use. 
        self.__aui_manager = aui_manager

        self.SetIcon(images.mondrian_like_icon3_pix32.GetIcon())
        #self.SetIcon(images.Mondrian.GetIcon())

        self.statusbar = self.CreateStatusBar(4, 0)
        self.statusbar.SetStatusText("Ready")
        
        # I make the status bar globally available because multiple places
        # in the app want to use it.
        wx.GetApp().vespa.statusbar = self.statusbar

        util_menu.bar = util_menu.PulseMenuBar(self)
        self.SetMenuBar(util_menu.bar)
        
        create0 = self.db.get_transform_kernel_create_menu_items()
        modify0 = self.db.get_transform_kernel_modify_menu_items()
        
        # search create list for any we want hidden
        #
        # using this to remove buggy kernels whose presence would be missed if 
        # user opened a pulse project which had used this kernel previously. Point
        # here is to make it hard for buggy kernels to be used going forward.
        #
        # FIXME - bjs, May want to move this list to constants in future 
        
        hide = ['3cb7a9e2-bb23-4934-ae17-96c01b8d93c7',     # Import from RFPulse - only needed for upgrade.
                '814e1d0e-7068-4583-9110-72f0ccf4914e',]    # orig Import from File with int() dwell bug
        
        create = [item for item in create0 if item[0] not in hide]
        modify = [item for item in modify0 if item[0] not in hide]
                
        util_menu.bar.update_transforms(create,modify)
        

        self.build_panes(aui_manager)
        self.bind_events()



    ##############                                    ############
    ##############     Internal helpers are below     ############
    ##############       in alphabetical order        ############
    ##############                                    ############

    def bind_events(self):
        self.Bind(wx.EVT_CLOSE, self.on_self_close)
        self.Bind(wx.EVT_SIZE, self.on_self_coordinate_change)
        self.Bind(wx.EVT_MOVE, self.on_self_coordinate_change)
        self.Bind(wx.EVT_ERASE_BACKGROUND, self.on_erase_background)
        
        
    def build_panes(self, aui_manager):

        self.pulse_design_notebook = notebook_pulse_designs.NotebookPulseDesigns(self, self.db)

        # create center pane
        aui_manager.AddPane(self.pulse_design_notebook, 
                            aui.AuiPaneInfo().
                            Name("pulse designs").
                            CenterPane().
                            PaneBorder(False))

        # "commit" all changes made to AuiManager
        aui_manager.Update()


    ##############                                    ############
    ##############      Event handlers are below      ############
    ##############       in alphabetical order        ############
    ##############                                    ############

    def on_erase_background(self, event):
        event.Skip()


    def on_self_close(self, event):
        if self.pulse_design_notebook.close_all():
            # Save my coordinates
            config = rf_config.Config()

            config.set_window_coordinates("main", self._left, self._top, 
                                          self._width, self._height)
            config.write()
            
            self.__aui_manager.UnInit()      # needed to avoid wx._core.wxAssertionError: C++ assertion "GetEventHandler() == this" failed at ..\..\src\common\wincmn.cpp 
            self.Destroy()


    def on_self_coordinate_change(self, event):
        # This is invoked for move & size events
        if self.IsMaximized() or self.IsIconized():
            # Bah, forget about this. Recording coordinates doesn't make sense
            # when the window is maximized or minimized. This is only a
            # concern on Windows; GTK and OS X don't produce move or size
            # events when a window is minimized or maximized.
            pass
        else:
            if event.GetEventType() == wx.wxEVT_MOVE:
                self._left, self._top = self.GetPosition()
            else:
                # This is a size event
               self._width, self._height = self.GetSize()
        

    ##############                                    ############
    ##############       Menu handlers are below      ############
    ##############       in the order they appear     ############
    ##############             on the menu            ############
    ##############                                    ############

    ############    Pulse Design menu

    def on_new_pulse_design(self, event):
        wx.SetCursor(wx.HOURGLASS_CURSOR)
        self.pulse_design_notebook.add_pulse_design_tab()
        wx.SetCursor(wx.NullCursor)


    def on_open_pulse_design(self, event):

        dialog = dialog_pulse_design_browser.DialogPulseDesignBrowser(self, self.db)
        dialog.ShowModal()
        if dialog.selected_pulse_design_id:
            wx.SetCursor(wx.HOURGLASS_CURSOR)
            pulse_design = self.db.fetch_pulse_design(dialog.selected_pulse_design_id)
            # At this point we have retrieved waveform/gradient results, but we
            # do not have Bloch simulated the mx,my,mz profile results. We will
            # do this step within each tab, since there are display dependent 
            # variables that have to be retrieved at the tab_transform level.
            self.pulse_design_notebook.add_pulse_design_tab(pulse_design)
            wx.SetCursor(wx.NullCursor)


    def on_run_pulse_design(self, event):
        if self.pulse_design_notebook.active_tab:
            self.pulse_design_notebook.active_tab.run()


    def on_save_pulse_design(self, event):
        self.pulse_design_notebook.save_pulse_design()


    def on_save_as_pulse_design(self, event):
        self.pulse_design_notebook.save_as_pulse_design()


    def on_close_pulse_design(self, event):
        self.pulse_design_notebook.close_pulse_design()


    def on_copy_pulse_design_tab(self, event):
        wx.SetCursor(wx.HOURGLASS_CURSOR)
        self.Freeze()
        self.pulse_design_notebook.copy_tab()
        self.Thaw()
        wx.SetCursor(wx.NullCursor)


    def on_third_party_export(self, event):
        tab = self.pulse_design_notebook.active_tab
        
        if tab:
            if tab.is_synced:
                # All of the tabs in the pulse design are in sync.
                export = True
            else:
                # At least one tab is not in sync. Give the user the option
                # to halt the export. 
                msg = "You have edited this design since it was "      \
                      "last opened or run. Pulse can only export "    \
                      "the design as it was when it was last opened "  \
                      "or run.\n\n"                                     \
                      "Do you want to continue with the export?"
          
                export = (wx.YES == common_dialogs.message(msg, "Continue Export", 
                                                           common_dialogs.Q_YES_NO))

            if export:       
                pulse     = tab.pulse_design.get_pulse()       # returns MinimalistPulse obj
                bandwidth = tab.pulse_design.get_bandwidth()   # in kHz
                tip_angle = tab.pulse_design.get_tip_angle()   # in degrees
                
                if pulse:
                    dialog = dialog_third_party_export.DialogThirdPartyExport(self, 
                                                              tab.pulse_design.id,
                                                              pulse, bandwidth, tip_angle)
                    dialog.ShowModal()
                else:
                    msg = """This pulse design doesn't have a final """    \
                          """pulse. You can only use export completed """   \
                          """pulse designs to third party formats."""
                    common_dialogs.message(msg)
        #else:
            # Ignore the menu click -- Only the welcome tab is open
        


    ############    Management  menu   ###############

    def on_manage_pulse_designs(self, event):
        dialog = dialog_manage_pulse_designs.DialogManagePulseDesigns(self, self.db)
        dialog.ShowModal()

    def on_manage_transform_kernels(self, event):
        dialog = dialog_manage_transform_kernels.DialogManageTransformKernels(self, self.db)
        dialog.ShowModal()
        # update items in the Add Transforms menu
        create = self.db.get_transform_kernel_create_menu_items()
        modify = self.db.get_transform_kernel_modify_menu_items()
        util_menu.bar.update_transforms(create,modify)        
        
    def on_manage_machine_specs_templates(self, event):
        dialog = dialog_manage_machine_specs.DialogManageMachineSpecs(self, self.db)
        dialog.ShowModal()


    ############    Transforms menu    ############

    def on_transform_item(self, event):
        # This is triggered for all transform items created from entries in the
        # Vespa database. They hang off of "Add Transforms".

        # Get the id of the menu item that triggered this event
        id_ = event.GetId()

        # Find the transforms that are in the database
        kernel_id, kernel_name = util_menu.bar.get_transform_menu_item_info(id_)
        
        #print "Add Transform, kernel_name = "+kernel_name+"  kernel_id = "+str(kernel_id)
        
        kernel = self.db.fetch_transform_kernel(kernel_id)
         
        self.Freeze()
        self.pulse_design_notebook.add_transform(kernel)
        self.Thaw()



    ############    View  menu
    # View options affect only the experiment tab and so it's up to the
    # experiment notebook to react to them.
            
    def on_menu_view_option(self, event):    
        self.pulse_design_notebook.on_menu_view_option(event)

    def on_menu_view_output(self, event):
        self.pulse_design_notebook.on_menu_view_output(event)

    ############    Help menu

    def on_user_manual(self, event):
        path = util_misc.get_vespa_install_directory()
        path = os.path.join(path, "docs", "pulse_user_manual.pdf")
        wx_util.display_file(path)

    def on_pulse_help_online(self, event):
        webbrowser.open("http://scion.duhs.duke.edu/vespa/pulse", 1)

    def on_vespa_help_online(self, event):
        webbrowser.open("http://scion.duhs.duke.edu/vespa", 1)

    def on_about(self, event):
        bit = str(8 * struct.calcsize('P')) + '-bit Python'
        info = wx_adv.AboutDialogInfo()
        info.SetVersion(util_misc.get_vespa_version())
        info.SetCopyright("Copyright 2010, Brian J. Soher. All rights reserved.")
        info.SetDescription("Pulse is an advanced toolkit for generating RF pulses for MR simulations. Running on "+bit)
        info.SetWebSite("http://scion.duhs.duke.edu/vespa/")
        wx_adv.AboutBox(info)


    def on_show_inspection_tool(self, event):
        wx_util.show_wx_inspector(self)





#--------------------------------------------------------------------

if __name__ == "__main__":
    app, db_filename = util_init.init_app("Pulse")

    db = util_db.Database(db_filename)

    if util_misc.get_platform() == "windows":
        # code to allow Windows to set Windows taskbar icon correctly
        # - without this, Win was using default icon for python.exe
        # - bugfix from, https://stackoverflow.com/questions/1551605/how-to-set-applications-taskbar-icon-in-windows-7/1552105#1552105

        import ctypes
        # myappid is an arbitrary string - but different for each Vespa App
        myappid = u'mycompany.myproduct.subproduct.pulse' 
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)    

    # All my other settings are in pulse.ini
    config = rf_config.Config()
    
    position, size = config.get_window_coordinates("main")
    frame = Main(db, position, size)
    app.SetTopWindow(frame)
    frame.Show()
    app.MainLoop()

