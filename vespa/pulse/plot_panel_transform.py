# Python modules


# 3rd party modules
import wx

# Our modules
import vespa.common.wx_gravy.plot_panel as plot_panel
import vespa.common.wx_gravy.plot_panel_points as plot_panel_points


class PlotPanelTransform(plot_panel_points.PlotPanelPoints):

    def __init__(self, parent, **kwargs):

        plot_panel_points.PlotPanelPoints.__init__(self, parent, **kwargs)

        # The status bar (which this code uses a lot) is available on the
        # global vespa object which hangs off of the app object.
        # Here we create an alias for the call to SetStatusText(). In 
        # Python, each "." in a reference takes a function call or two 
        # to resolve (hasattr() and getattr()), so fewer dots == fewer calls
        # == more efficient. This rarely matters, but in this code which is
        # called on every mouse move, a little extra efficiency can't hurt.
        try:
            self.set_status_text = wx.GetApp().vespa.statusbar.SetStatusText
        except:
            self.top = parent.GetTopLevelParent()
            self.set_status_text = self.redirect_status

        self.set_color((255, 255, 255))

    def redirect_status(self, msg, loc):
        if loc == 0:
            self.top.LabelStatus0.SetLabel(msg)
        elif loc == 1:
            self.top.LabelStatus1.SetLabel(msg)
        elif loc == 2:
            self.top.LabelStatus2.SetLabel(msg)
        elif loc == 3:
            self.top.LabelStatus3.SetLabel(msg)

    # EVENT FUNCTIONS -----------------------------------------------

    def on_motion(self, xdata, ydata, val, bounds, iaxis):
        self.set_status_text(" X = %.3f" % xdata, 0)
        self.set_status_text(" Y = %.5f" % ydata, 1)
        self.set_status_text(" Value = %.5f" % val[0], 2)
        self.set_status_text(" ", 3)

    def on_zoom_select(self, xmin, xmax, val, ymin, ymax, reset=False, iplot=None):
        pass

    def on_zoom_motion(self, xmin, xmax, val, ymin, ymax, iplot=None):
        self.set_status_text(" X Range = %.3f to %.3f" % (xmin, xmax), 0)
        self.set_status_text(" Y Range = %.5f to %.5f" % (ymin, ymax), 1)
        self.set_status_text(" Value = %.5f" % val[0], 2)
        self.set_status_text(" dX = %.2f  dY = %.1f" % (xmax - xmin, ymax - ymin), 3)

    def on_refs_select(self, xmin, xmax, val, iplot=None):
        pass

    def on_refs_motion(self, xmin, xmax, val, iplot=None):
        self.set_status_text(" X Range = %.3f to %.3f" % (xmin, xmax), 0)
        self.set_status_text(" ", 1)
        self.set_status_text(" Value = %.5f" % val[0], 2)
        self.set_status_text(" dX = %.2f " % (xmax - xmin), 3)
