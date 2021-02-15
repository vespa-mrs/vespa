# Python modules


# 3rd party modules

# Our modules
import vespa.common.wx_gravy.plot_panel as plot_panel


        

class PlotPanelEditor(plot_panel.PlotPanel):
    
    def __init__(self, parent, tab, status, **kwargs):
        
        plot_panel.PlotPanel.__init__( self, parent, **kwargs )

        self.parent = parent
        self.tab    = tab
        self.status = status
        
        self.set_color( (255,255,255) )


    # EVENT FUNCTIONS -----------------------------------------------
    
    def on_motion(self, xdata, ydata, val, bounds, iaxis):
        txt = 'X,Y=%1.4f,%1.4f    Value=%s' % (xdata,ydata,str(val))
        self.status.SetLabel(txt)
            
    def on_zoom_select(self, xmin, xmax, val, ymin, ymax, reset=False):
        txt = 'Xmin=%1.4f  Xmax=%1.4f  Width=%1.4f    Value=%s' % (xmin,xmax, xmax-xmin, str(val))
        self.status.SetLabel(txt)

    def on_zoom_motion(self, xmin, xmax, val, ymin, ymax):
        if ymin is not None:
            txt = 'Xmin=%1.4f Xmax=%1.4f   Ymin=%1.4f Ymax=%1.4f   Value=%s' % (xmin,xmax,ymin,ymax,str(val))
        else:
            txt = 'Xmin=%1.4f  Xmax=%1.4f  Width=%1.4f    Value=%s' % (xmin,xmax, xmax-xmin, str(val))
        self.status.SetLabel(txt)

    def on_refs_select(self, xmin, xmax, val):
        txt = 'Xmin=%1.4f  Xmax=%1.4f  Width=%1.4f    Value=%s' % (xmin,xmax, xmax-xmin, str(val))
        self.status.SetLabel(txt)

    def on_refs_motion(self, xmin, xmax, val):
        txt = 'Xmin=%1.4f  Xmax=%1.4f  Width=%1.4f    Value=%s' % (xmin,xmax, xmax-xmin, str(val))
        self.status.SetLabel(txt)


    def on_middle_select(self, xstr, ystr, xend, yend, iplot):
        pass


    def on_middle_motion(self, xcur, ycur, xprev, yprev, iplot):
        # The mouse probably moved in both the X and Y directions, but to
        # make phasing easier for the user to do accurately, we only pay
        # attention to the larger of the two movements.

        dy = ycur - yprev
        dx = xcur - xprev
        delta  = 0
        deltab = 0

        if abs(dy) > abs(dx):
            # 0 order phase
            if dy != 0.0:
                self.tab.phase0 = (self.tab.phase0 + dy) % 360
                self.tab.FloatPhase0.SetValue(self.tab.phase0)
        else:
            # first order phase, x10 makes Phase1 changes more noticeable
            if dx != 0.0:
                self.tab.phase1 = self.tab.phase1 + dx*10
                self.tab.FloatPhase1.SetValue(self.tab.phase1)
        self.tab.plot_canvas()

    
