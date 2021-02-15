# Python modules


# 3rd party modules
import wx

# Our modules
import vespa.pulse.auto_gui.view_pulse_design as view_pulse_design
import vespa.pulse.dialog_referrers as dialog_referrers


class DialogViewPulseDesign(view_pulse_design.MyDialog):
    """Displays the dialog for viewing a pulse_design in HTML"""
    def __init__(self, parent, pulse_design):
        if not parent:
            parent = wx.GetApp().GetTopWindow()

        view_pulse_design.MyDialog.__init__(self, parent) 
        
        self.pulse_design = pulse_design

        html_sizer = self.LabelHtmlPlaceholder.GetContainingSizer()
        parent = self.LabelHtmlPlaceholder.GetParent()
        self.LabelHtmlPlaceholder.Destroy()
        
        self.html_ctrl = wx.html.HtmlWindow(parent)
        
        html_sizer.Add(self.html_ctrl, 1, wx.EXPAND|wx.ALIGN_TOP)

        self.html_ctrl.SetPage(pulse_design.as_html())
        
        self.SetTitle("View Pulse Design %s" % pulse_design.name)

        self.SetSize( (500, 600) )

        self.Layout()

        self.Center()


    def on_show_referrers(self, event):
        dialog = dialog_referrers.DialogReferrers(self, self.pulse_design)
        dialog.ShowModal()
