# Python modules


# 3rd party modules
import wx

# Our modules
import vespa.pulse.auto_gui.referrers as gui_referrers
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


class DialogReferrers(gui_referrers.MyDialog):
    def __init__(self, parent, pulse_design):
        if not parent:
            parent = wx.GetApp().GetTopWindow()

        gui_referrers.MyDialog.__init__(self, parent)

        self.pulse_design = pulse_design
        
        self.SetTitle("Pulse Sequences Using %s" % pulse_design.name)

        names = [referrer[1] for referrer in pulse_design.referrers]

        self.ListPulseSequences.SetItems(names)

        self.Layout()

        self.Center()

        self.ButtonClose.SetFocus()


    def on_copy(self, event):
        s = "Vespa %s:\n" % self.GetTitle()
        s += "\n".join([referrer[1] for referrer in self.pulse_design.referrers])
        wx_util.copy_to_clipboard(s)


    def on_close(self, event):
        self.Close()


