# Python modules


# 3rd party modules
import wx
import numpy as np

# Our modules
import vespa.common.wx_gravy.plot_panel_spectrum as plot_panel_spectrum
        

class PlotPanelDatasim(plot_panel_spectrum.PlotPanelSpectrum):
    
    def __init__(self, parent, tab, tab_dataset, **kwargs):

        tab.SizerSplitterWindow.Fit(tab)         # bugfix wxGlade 0.9.6 to 1.0.0

        plot_panel_spectrum.PlotPanelSpectrum.__init__( self, parent, **kwargs )

        # tab is the containing widget for this plot_panel, it is used
        # in resize events, the tab attribute is the AUI Notebook tab
        # that contains this plot_panel

        self.tab   = tab  
        self.tab_dataset = tab_dataset
        self.top   = wx.GetApp().GetTopWindow()
        
        self.ref_locations = 0,self.tab.datasim.spectral_dims[0]-1     # reference span locations in ppm
        
        self.set_color( (255,255,255) )


    # EVENT FUNCTIONS -----------------------------------------------
    
    def on_motion(self, xdata, ydata, val, bounds, iaxis):
        
        if self.prefs.xaxis_hertz:
            hz  = xdata
            ppm = xdata / self.tab.datasim.frequency
        else:
            ppm = xdata
            hz  = xdata * self.tab.datasim.frequency
        
        value = 0.0
        if iaxis in [0,1,2]:
            value = val[0]
            
        self.top.statusbar.SetStatusText( " PPM = %.3f" % (ppm, ), 0)
        self.top.statusbar.SetStatusText( " Hz = %.3f"  % (hz,  ), 1)
        self.top.statusbar.SetStatusText(( " Value = "+str(value)), 2)   

        
    
    def on_scroll(self, button, step, iaxis):
        
        self.set_vertical_scale(step)
        self.tab.FloatScale.SetValue(self.vertical_scale)
            
            
    def on_zoom_select(self, xmin, xmax, val, ymin, ymax, reset=False, iplot=None):
        if reset:
            # we only need to bother with setting the vertical scale in here
            # if we did a reset of the x,y axes.
            self.vertical_scale = self.dataymax
            self.tab.FloatScale.SetValue(self.dataymax)
        

    def on_zoom_motion(self, xmin, xmax, val, ymin, ymax, iplot=None):
        
        if self.prefs.xaxis_hertz:
            hz_str  = xmin
            hz_end  = xmax
            if hz_str > hz_end: hz_str, hz_end = hz_end, hz_str
            ppm_str = hz_str / self.tab.datasim.frequency
            ppm_end = hz_end / self.tab.datasim.frequency
        else: 
            ppm_str  = xmin
            ppm_end  = xmax
            if ppm_str > ppm_end: ppm_str, ppm_end = ppm_end, ppm_str
            hz_str = ppm_str * self.tab.datasim.frequency
            hz_end = ppm_end * self.tab.datasim.frequency
        delta_ppm = -1*(ppm_str - ppm_end)  # keeps delta positive
        delta_hz  = -1*(hz_str - hz_end)
        self.top.statusbar.SetStatusText(( " PPM Range = %.2f to %.2f" % (ppm_str, ppm_end)), 0)
        self.top.statusbar.SetStatusText(( " Hz Range = %.1f to %.1f" % (hz_str, hz_end)), 1)
        self.top.statusbar.SetStatusText(( " dPPM = %.2f  dHz = %.1f" % (delta_ppm, delta_hz)), 2)
            
        
    def on_refs_select(self, xmin, xmax, val, reset=False, iplot=None):

        ppm_str = xmax
        ppm_end = xmin
        if self.prefs.xaxis_hertz:
            ppm_str = ppm_str / self.tab.datasim.frequency
            ppm_end = ppm_end / self.tab.datasim.frequency

        if ppm_str < ppm_end: ppm_str, ppm_end = ppm_end, ppm_str

        # Calculate area of span
        all_areas, all_rms = self.calculate_area()
        area = all_areas[0]
        rms  = all_rms[0]
        self.top.statusbar.SetStatusText(' Area = %1.5g  RMS = %1.5g' % (area,rms), 3)


    def on_refs_motion(self, xmin, xmax, val, iplot=None):

        if self.prefs.xaxis_hertz:
            hz_str  = xmin
            hz_end  = xmax
            if hz_str > hz_end: hz_str, hz_end = hz_end, hz_str
            ppm_str = hz_str / self.frequency
            ppm_end = hz_end / self.frequency
        else:
            ppm_str  = xmin
            ppm_end  = xmax
            if ppm_str > ppm_end: ppm_str, ppm_end = ppm_end, ppm_str
            hz_str = ppm_str * self.frequency
            hz_end = ppm_end * self.frequency
        delta_ppm = -1*(ppm_str - ppm_end)  # keeps delta positive
        delta_hz  = -1*(hz_str - hz_end)
        self.top.statusbar.SetStatusText(( " PPM Range = %.2f to %.2f" % (ppm_str, ppm_end)), 0)
        self.top.statusbar.SetStatusText(( " Hz Range = %.1f to %.1f" % (hz_str, hz_end)), 1)
        self.top.statusbar.SetStatusText(( " dPPM = %.2f  dHz = %.1f" % (delta_ppm, delta_hz)), 2)

        all_areas, all_rms = self.calculate_area()
        area = all_areas[0]
        rms  = all_rms[0]
        self.top.statusbar.SetStatusText(' Area = %1.5g  RMS = %1.6g' % (area,rms), 3)
        

    def on_middle_select(self, xstr, ystr, xend, yend, iplot):
        pass


    def on_middle_press(self, xloc, yloc, iplot):
        pass

        
    def on_middle_motion(self, xcur, ycur, xprev, yprev, iplot):
        pass
        
    
        
  
    
