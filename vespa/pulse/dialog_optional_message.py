# Python modules


# 3rd party modules
import wx

# Our modules
import vespa.pulse.auto_gui.optional_message as optional_message
import vespa.common.wx_gravy.util as wx_util


class DialogOptionalMessage(optional_message.MyDialog):
    """Displays a messgae box with one unusual feature, which is a 
    checkbox that says "Don't show this message again". If the user checks
    that before closing the dialog, the dialog attribute dont_show_again
    is True when the dialog returns.
    
    The message param to the constructor is the message that will be 
    displayed to the user.
    
    Note that this "message box" doesn't offer a lot of the features of 
    a regular message box. For instance, it has only OK & Cancel buttons; 
    there's no way to display Yes & No buttons instead, for instance. Those 
    features could be implemented, but at present they're not.
    """
    def __init__(self, parent, message, title=""):
        if not parent:
            parent = wx.GetApp().GetTopWindow()
            
        self.dont_show_again = False

        optional_message.MyDialog.__init__(self, parent)
        
        # Default return code is cancel
        self.SetReturnCode(wx.CANCEL)

        self.LabelMessage.SetLabel(message)
        
        if not title:
            title = wx.GetApp().GetAppName()
        self.SetTitle(title)
        
        # We add the OK & Cancel buttons dynamically so that they're in the 
        # right order under OS X, GTK, Windows, etc.
        self.ButtonOk, self.ButtonCancel = \
                    wx_util.add_ok_cancel(self, self.LabelOkCancelPlaceholder,
                                          self.on_ok)

        self.ButtonOk.SetFocus()
        self.ButtonOk.SetDefault()

        self.SetSize( (500, 200) )

        self.Layout()

        # There's one last thing left to do. In wx, labels (wxStaticText) 
        # don't wrap automatically, and the message in self.LabelWarning is
        # often longer than can display in the label once the dialog is
        # sized (above). 
        # wxStaticText objects have a .Wrap() method that requires an explicit
        # pixel width. The correct width to use is that of the sizer, but 
        # the sizer's width seems to include it's left & right borders which
        # are not meant to be drawn on. So I get the sizer's width, shrink
        # it a little, pass that to .Wrap() and then lay out the dialog 
        # again in case LabelMessage is several lines tall. 
        width = self.LabelMessage.GetContainingSizer().GetSize().width
        
        self.LabelMessage.Wrap(int(width * .97))
        
        self.Layout()
        
        self.Fit()
        
        self.Center()
        

    def on_ok(self, event):
        if self.CheckDontShowAgain.IsChecked():
            self.dont_show_again = True

        self.EndModal(wx.OK)

