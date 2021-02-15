# Python modules


# 3rd party modules
import numpy as np

# Our modules
import vespa.common.util.ppm as util_ppm
import vespa.common.wx_gravy.plot_panel as plot_panel
import vespa.common.wx_gravy.util as wx_util
from wx.lib.agw.floatspin import FloatSpin, EVT_FLOATSPIN, FS_LEFT, FS_RIGHT, FS_CENTRE, FS_READONLY,  FixedPoint

        

class PlotPanelPlot1d(plot_panel.PlotPanel):
    
    def __init__(self, parent, tab, prefs, **kwargs):
        
        plot_panel.PlotPanel.__init__( self, parent, **kwargs )

        # parent is the containing panel for this plot_panel, it is used
        # in resize events, the tab attribute is the AUI Notebook tab
        # that contains this plot_panel

        self.parent = parent
        self.tab    = tab
        self._prefs = prefs
        self.count  = 0
        
        self.set_color( (255,255,255) )


    # EVENT FUNCTIONS -----------------------------------------------
    
    def on_motion(self, xdata, ydata, val, bounds, iaxis):
        if self._prefs.xaxis_hertz:
            hz  = xdata
            ppm = self.tab.freq1d.pts2ppm(self.tab.freq1d.hz2pts(xdata))
        else:
            ppm = xdata
            hz  = self.tab.freq1d.pts2hz(self.tab.freq1d.ppm2pts(xdata))
        
        # if in a stack plot, determine which plot we are in so we can
        # display the appropriate value from the val list
        ydiff = bounds[3]
        if ydiff == 0: ydiff = 1
        iplot = int(len(val) * (ydata / ydiff))
        iplot = np.clip(iplot, 0, len(val)-1)
        
        # show metabolite name label
        metabolite_label = ''
        if self.tab.plot_labels:
            mode = self.tab.ChoiceDisplayMode.GetCurrentSelection()
            if mode == 0:
                metabolite_label = self.tab.plot_labels[iplot]
            else:
                metabolite_label = self.tab.plot_labels[-1]
            if self.tab.CheckSumPlots.IsChecked():
                metabolite_label = "Summed Metabolites"
                
        self.tab.statusbar.SetStatusText(( "PPM = "  +str(ppm)  ), 0)
        self.tab.statusbar.SetStatusText(( "Hz = "   +str(hz)   ), 1)
        self.tab.statusbar.SetStatusText(( "Value = "+str(val[iplot])), 2)   
        self.tab.statusbar.SetStatusText(( metabolite_label  ), 3)
        
            
    def on_zoom_select(self, xmin, xmax, val, ymin, ymax, reset=False):
        if self._prefs.xaxis_hertz:
            xmax = self.tab.freq1d.pts2ppm(self.tab.freq1d.hz2pts(xmax))
            xmin = self.tab.freq1d.pts2ppm(self.tab.freq1d.hz2pts(xmin))
            xmax = np.clip(xmax, self.tab.minppm, self.tab.maxppm)
            xmin = np.clip(xmin, self.tab.minppm, self.tab.maxppm)
        self.tab.FloatXaxisMin.SetValue(xmin)
        self.tab.FloatXaxisMax.SetValue(xmax)


    def on_zoom_motion(self, xmin, xmax, val, ymin, ymax):
        if self._prefs.xaxis_hertz:
            hz_str = xmin
            hz_end = xmax
            if hz_str > hz_end: hz_str, hz_end = hz_end, hz_str
            ppm_str = hz_str / self.tab.freq1d.frequency
            ppm_end = hz_end / self.tab.freq1d.frequency
        else:
            ppm_str = np.clip(xmin, self.tab.minppm, self.tab.maxppm)
            ppm_end = np.clip(xmax, self.tab.minppm, self.tab.maxppm)
            if ppm_str > ppm_end: ppm_str, ppm_end = ppm_end, ppm_str
            hz_str = ppm_str * self.tab.freq1d.frequency
            hz_end = ppm_end * self.tab.freq1d.frequency

        self.tab.FloatXaxisMax.SetValue(FixedPoint(str(xmax),10))
        self.tab.FloatXaxisMin.SetValue(FixedPoint(str(xmin),10))
        
        delta_ppm = -1*(ppm_str - ppm_end)  # keeps delta positive
        delta_hz  = delta_ppm * self.tab.freq1d.frequency
        self.tab.statusbar.SetStatusText((" PPM Range = %.2f to %.2f" % (ppm_end, ppm_str)), 0)
        self.tab.statusbar.SetStatusText((" Hz Range = %.1f to %.1f" % (hz_str, hz_end)), 1)
        self.tab.statusbar.SetStatusText((" dPPM = %.2f  dHz = %.1f" % (delta_ppm, delta_hz)), 2)


    def on_refs_select(self, xmin, xmax, val):
        
        if self._prefs.xaxis_hertz:
            xmax = self.tab.freq1d.pts2ppm(self.tab.freq1d.hz2pts(xmax))
            xmin = self.tab.freq1d.pts2ppm(self.tab.freq1d.hz2pts(xmin))
            xmax = np.clip(xmax, self.tab.minppm, self.tab.maxppm)
            xmin = np.clip(xmin, self.tab.minppm, self.tab.maxppm)
        self.tab.FloatCursorMax.SetValue(FixedPoint(str(xmax),10))
        self.tab.FloatCursorMin.SetValue(FixedPoint(str(xmin),10))
        self.tab.plot_integral()
        self.tab.plot_contour()
        self.count = 0

    def on_refs_motion(self, xmin, xmax, val):

        if self._prefs.xaxis_hertz:
            hz_str = xmin
            hz_end = xmax
            if hz_str > hz_end: hz_str, hz_end = hz_end, hz_str
            ppm_str = hz_str / self.tab.freq1d.frequency
            ppm_end = hz_end / self.tab.freq1d.frequency
        else:
            ppm_str = xmin
            ppm_end = xmax
            if ppm_str > ppm_end: ppm_str, ppm_end = ppm_end, ppm_str
            hz_str = ppm_str * self.tab.freq1d.frequency
            hz_end = ppm_end * self.tab.freq1d.frequency

        self.tab.FloatCursorMax.SetValue(FixedPoint(str(xmax),10))
        self.tab.FloatCursorMin.SetValue(FixedPoint(str(xmin),10))

        delta_ppm = -1 * (ppm_str - ppm_end)
        delta_hz = delta_ppm * self.tab.freq1d.frequency
        self.tab.statusbar.SetStatusText((" PPM Range = %.2f to %.2f" % (ppm_end, ppm_str)), 0)
        self.tab.statusbar.SetStatusText((" Hz Range = %.1f to %.1f" % (hz_str, hz_end)), 1)
        self.tab.statusbar.SetStatusText((" dPPM = %.2f  dHz = %.1f" % (delta_ppm, delta_hz)), 2)
        self.count += 1
        if (self.count % 2) == 0:
            self.tab.plot_integral()

