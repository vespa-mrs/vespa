# Python modules



# 3rd party modules
import numpy as np


# Our modules
import vespa.common.wx_gravy.plot_panel as plot_panel


        

class PlotPanelIntegral(plot_panel.PlotPanel):
    
    def __init__(self, parent, tab, **kwargs):
        
        plot_panel.PlotPanel.__init__( self, parent, **kwargs )

        # parent is the containing panel for this plot_panel, it is used
        # in resize events, the tab attribute is the AUI Notebook tab
        # that contains this plot_panel

        self.parent = parent
        self.tab    = tab
        
        self.set_color( (255,255,255) )


    # EVENT FUNCTIONS -----------------------------------------------
    
    def on_motion(self, xdata, ydata, value, bounds, iaxis):

        if self.tab.integral_data is None: return
        npts = self.tab.integral_data.shape[0]
        xindx = int(npts * (-(xdata - bounds[0]) / (-bounds[2])))
        xindx = np.clip(xindx, 0, npts-1)+1
        self.tab.statusbar.SetStatusText(( "X-Index = "+str(xindx)), 0)
        self.tab.statusbar.SetStatusText(( "X-Axis  = "+str(xdata)), 1)
        self.tab.statusbar.SetStatusText(( "Value = "  +str(value)), 2)   

            
    def on_zoom_select(self, xmin, xmax, val, ymin, ymax, reset=False):
        pass

    def on_zoom_motion(self, xmin, xmax, val, ymin, ymax):
        pass

    def on_refs_select(self, xmin, xmax, val):
        pass

    def on_refs_motion(self, xmin, xmax, val):
       pass
   
    def on_middle_press(self, xloc, yloc, indx, bounds=None, xdata=None, ydata=None):
        
        if self.tab.integral_data is None: return
        
        if bounds is not None and xdata is not None:
            npts = self.tab.integral_data.shape[0]
            xindx = int(npts * (-(xdata - bounds[0]) / (-bounds[2])))
            xindx = np.clip(xindx, 0, npts-1)+1

            # check for display mode and fill x-axis
            mode = self.tab.ChoiceDisplayMode.GetCurrentSelection()
            if mode == 0:
                return
            elif mode == 1:
                self.tab.SpinIndex1.SetValue(xindx)
                self.tab.LabelIndex1Value.SetLabel(str(self.tab.experiment.dims[0][xindx-1]))
            elif mode == 2:
                self.tab.SpinIndex2.SetValue(xindx)
                self.tab.LabelIndex2Value.SetLabel(str(self.tab.experiment.dims[1][xindx-1]))
            elif mode == 3:
                self.tab.SpinIndex3.SetValue(xindx)
                self.tab.LabelIndex3Value.SetLabel(str(self.tab.experiment.dims[2][xindx-1]))
            else:
                return
            self.tab.plot_canvas()
            self.tab.plot_contour()




    
