# Python modules


# 3rd party modules
import wx

# Our modules
import vespa.pulse.auto_gui.view_transform_kernel as view_transform_kernel
import vespa.pulse.dialog_referrers as dialog_referrers


class DialogViewTransformKernel(view_transform_kernel.MyDialog):
    """
    Displays the dialog for viewing a transform_kernel in HTML
    
    """
    def __init__(self, parent, transform_kernel):
        
        if not parent:
            parent = wx.GetApp().GetTopWindow()

        view_transform_kernel.MyDialog.__init__(self, parent) 
        
        self.transform_kernel = transform_kernel

        html_sizer = self.LabelHtmlPlaceholder.GetContainingSizer()
        parent = self.LabelHtmlPlaceholder.GetParent()
        self.LabelHtmlPlaceholder.Destroy()
        
        self.html_ctrl = wx.html.HtmlWindow(parent)
        
        html_sizer.Add(self.html_ctrl, 1, wx.EXPAND|wx.ALIGN_TOP)

        self.html_ctrl.SetPage(transform_kernel.as_html())
        
        self.SetTitle("View Transform Kernel %s" % transform_kernel.name)

        self.SetSize( (500, 600) )

        self.Layout()

        self.Center()


    def on_show_referrers(self, event):
        dialog = dialog_referrers.DialogReferrers(self, self.transform_kernel)
        dialog.ShowModal()
