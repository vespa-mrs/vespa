# Python modules



# 3rd party modules
import wx
import numpy as np

# Our modules
import vespa.common.wx_gravy.plot_panel as plot_panel




def find_nearest(a, a0):
    "Element in nd array `a` closest to the scalar value `a0`"
    return np.abs(a - a0).argmin()
    

class PlotPanelTransformContour(plot_panel.PlotPanel):
    """
    This version of PlotPanel allows you to show and interact with a contour
    map. Most data is not saved in a contour map, only the paths along the 
    isocontour lines.  I want to be able to report the z values at a given 
    x,y location. So here we save the x,y,z values every time we set data into
    the contour map using the plot_contour() method. This way we can track 
    the z values as the cursor is moved.
    """
    
    def __init__(self, parent, **kwargs):
        
        plot_panel.PlotPanel.__init__( self, parent, **kwargs )
        
        self.xvals = None
        self.yvals = None
        self.zvals = None

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
        
        self.set_color( (255,255,255) )


    def redirect_status(self, msg, loc):
        if loc == 0:
            self.top.LabelStatus0.SetLabel(msg)
        elif loc == 1:
            self.top.LabelStatus1.SetLabel(msg)
        elif loc == 2:
            self.top.LabelStatus2.SetLabel(msg)
        elif loc == 3:
            self.top.LabelStatus3.SetLabel(msg)


    def get_values(self, event):
        """
        Generic utility function that polls the axes that the mouse is within
        to return a list of data values at the x location of the cursor.
        
        """
        value = []
        x0, y0, x1, y1 = event.inaxes.dataLim.bounds
    
        if self.zvals is None:
            return [0.0]
    
        # need this specific test because of sporadic event calls where
        # event.inaxes == None and the tests below throw exceptions
        if event.inaxes:
            x  = event.xdata
            y  = event.ydata
            xx = find_nearest(self.xvals, x)
            yy = find_nearest(self.yvals, y)
            
            return [self.zvals[yy,xx]]
    
        return [0.0]


    def _on_move(self, event):
        """
        This is the internal method that organizes the data that is sent to the
        external user defined event handler for motion events. In here we 
        gather data values from either collection or line plots, determine 
        which axis we are in, then call the (hopefully) overloaded on_motion()
        method
        
        """
        if event.inaxes == None or not self.do_motion_event: return
        value = []
        x0, y0, x1, y1 = bounds = event.inaxes.dataLim.bounds
    
        # need this specific test because of sporadic event calls where
        # event.inaxes == None and the tests below throw exceptions
        if event.inaxes:
            value = self.get_values(event)

        if value == []: 
            value = [0.0]
        
        iaxis = None
        for i,axis in enumerate(self.axes):
            if axis == event.inaxes:
                iaxis = i

        self.on_motion(event.xdata, event.ydata, value, bounds, iaxis)    
    

    # EVENT FUNCTIONS -----------------------------------------------
    
    def on_motion(self, xdata, ydata, val, bounds, iaxis):
        self.set_status_text( " X = %.3f" % xdata, 0)
        self.set_status_text( " Y = %.5f" % ydata, 1)
        self.set_status_text( " Value = %.5f" % val[0], 2)
        self.set_status_text( " " , 3)
            
    def on_zoom_select(self, xmin, xmax, val, ymin, ymax, iplot=None):
        pass

    def on_zoom_motion(self, xmin, xmax, val, ymin, ymax, iplot=None):
        self.set_status_text( " X Range = %.3f to %.3f" % (xmin, xmax), 0)
        self.set_status_text( " Y Range = %.5f to %.5f" % (ymin, ymax), 1)
        self.set_status_text( " Value = %.5f" % val[0], 2)
        self.set_status_text( " dX = %.2f  dY = %.1f" % (xmax-xmin, ymax-ymin), 3)

    def on_refs_select(self, xmin, xmax, val, iplot=None):
        pass

    def on_refs_motion(self, xmin, xmax, val, iplot=None):
        self.set_status_text( " X Range = %.3f to %.3f" % (xmin, xmax), 0)
        self.set_status_text( " " , 1)
        self.set_status_text( " Value = %.5f" % val[0], 2)
        self.set_status_text( " dX = %.2f " % (xmax-xmin), 3)


    
