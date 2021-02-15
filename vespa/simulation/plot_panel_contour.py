# Python modules



# 3rd party modules
import numpy as np

# Our modules
import vespa.common.wx_gravy.plot_panel as plot_panel


        

class PlotPanelContour(plot_panel.PlotPanel):
    
    def __init__(self, parent, tab, **kwargs):
        
        plot_panel.PlotPanel.__init__( self, parent, **kwargs )

        # parent is the containing panel for this plot_panel, it is used
        # in resize events, the tab attribute is the AUI Notebook tab
        # that contains this plot_panel
        
        self.parent = parent
        self.tab    = tab
        
        self.set_color( (255,255,255) )


    # EVENT FUNCTIONS -----------------------------------------------
    
    def on_motion(self, xdata, ydata, val, bounds, iaxis):
        
        mode = self.tab.ChoiceContourMode.GetCurrentSelection()
        if mode == 0:
            xlen = len(self.tab.experiment.dims[0])
            ylen = len(self.tab.experiment.dims[1])
            label1 = "X,Y Locs = "
        elif mode == 1:
            xlen = len(self.tab.experiment.dims[1])
            ylen = len(self.tab.experiment.dims[2])
            label1 = "X,Z Locs = "
        elif mode == 2:
            xlen = len(self.tab.experiment.dims[1])
            ylen = len(self.tab.experiment.dims[2])
            label1 = "Y,Z Locs = "
        if xlen == 1 and ylen ==1:
            return
        if bounds[2] == 0 or bounds[3] == 0:    # before data loaded
            return

        xindx = int(xlen * (-(xdata - bounds[0]) / (-bounds[2])))
        xindx = np.clip(xindx, 0, xlen-1)+1
        yindx = int(ylen * (-(ydata - bounds[1]) / (-bounds[3])))
        yindx = np.clip(yindx, 0, ylen-1)+1
        self.tab.statusbar.SetStatusText(( label1+str(xindx)+","+str(yindx)), 0)
        self.tab.statusbar.SetStatusText("", 1)
        self.tab.statusbar.SetStatusText(( "Value = "  +str(val)), 2)   

        
            
    def on_zoom_select(self, xmin, xmax, val, ymin, ymax, reset=False):
        pass

    def on_zoom_motion(self, xmin, xmax, val, ymin, ymax):
        pass

    def on_refs_select(self, xmin, xmax, val):
        pass

    def on_refs_motion(self, xmin, xmax, val):
        pass

    def on_middle_press(self, xloc, yloc, indx, bounds=None, xdata=None, ydata=None):
        
        if bounds is not None and xdata is not None and ydata is not None:

            mode = self.tab.ChoiceContourMode.GetCurrentSelection()
            if mode == 0:
                xlen = len(self.tab.experiment.dims[0])
                ylen = len(self.tab.experiment.dims[1])
            elif mode == 1:
                xlen = len(self.tab.experiment.dims[1])
                ylen = len(self.tab.experiment.dims[2])
            elif mode == 2:
                xlen = len(self.tab.experiment.dims[1])
                ylen = len(self.tab.experiment.dims[2])
            if xlen == 1 and ylen ==1:
                return
            if bounds[2] == 0 or bounds[3] == 0:    # before data loaded
                return

            xindx = int(xlen * (-(xdata - bounds[0]) / (-bounds[2])))
            xindx = np.clip(xindx, 0, xlen-1)+1

            yindx = int(ylen * (-(ydata - bounds[1]) / (-bounds[3])))
            yindx = np.clip(yindx, 0, ylen-1)+1

            # check for display mode and fill x-axis
            if mode == 0:
                self.tab.SpinIndex1.SetValue(xindx)
                self.tab.SpinIndex2.SetValue(yindx)
                self.tab.LabelIndex1Value.SetLabel(str(self.tab.experiment.dims[0][xindx-1]))
                self.tab.LabelIndex2Value.SetLabel(str(self.tab.experiment.dims[1][yindx-1]))
            elif mode == 1:
                self.tab.SpinIndex1.SetValue(xindx)
                self.tab.SpinIndex3.SetValue(yindx)
                self.tab.LabelIndex1Value.SetLabel(str(self.tab.experiment.dims[0][xindx-1]))
                self.tab.LabelIndex3Value.SetLabel(str(self.tab.experiment.dims[2][yindx-1]))
            elif mode == 2:
                self.tab.SpinIndex2.SetValue(xindx)
                self.tab.SpinIndex3.SetValue(yindx)
                self.tab.LabelIndex2Value.SetLabel(str(self.tab.experiment.dims[1][xindx-1]))
                self.tab.LabelIndex3Value.SetLabel(str(self.tab.experiment.dims[2][yindx-1]))
            else:
                return
            self.tab.plot_canvas()
            self.tab.plot_integral()
    
