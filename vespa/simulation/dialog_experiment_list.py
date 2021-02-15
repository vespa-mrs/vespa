# Python modules


# 3rd party modules
import wx

# Our modules
import vespa.simulation.auto_gui.experiment_list as gui_experiment_list
import vespa.common.wx_gravy.common_dialogs as common_dialogs
import vespa.common.wx_gravy.util as wx_util

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


class DialogExperimentList(gui_experiment_list.MyDialog):
    def __init__(self, parent, metabolite):
        if not parent:
            parent = wx.GetApp().GetTopWindow()

        gui_experiment_list.MyDialog.__init__(self, parent)

        self.metabolite = metabolite
        
        self.SetTitle("Experiments Using %s" % metabolite.name)

        names = [name for name in metabolite.experiment_names]

        self.ListExperiments.SetItems(names)

        self.Layout()

        self.Center()

        self.ListExperiments.SetFocus()


    def on_copy(self, event):
        s = "Simulation %s:\n" % self.GetTitle()
        s += "\n".join([name for name in self.metabolite.experiment_names])
        wx_util.copy_to_clipboard(s)


    def on_close(self, event):
        self.Close()


