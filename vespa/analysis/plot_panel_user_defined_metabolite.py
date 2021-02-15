# Python modules

# 3rd party modules
import wx
import numpy as np

# Our modules
import vespa.common.wx_gravy.plot_panel_spectrum as plot_panel_spectrum

class PlotPanelUserDefinedMetabolite(plot_panel_spectrum.PlotPanelSpectrum):
    
    def __init__(self, parent, **kwargs):
        '''
        This is a customization of the PlotPanel object for use in a specific
        location in our application. The basic matplotlib functionality is
        contained in the base class, all we are doing here is overwriting the
        various canvas Event functions to mesh with our application.
        
        parent - the widget to which the PlotPanel is directly attached

        dataset - the currently active dataset
        
        '''
        super().__init__(parent, **kwargs)

        self._dataset = kwargs["dataset"]
        self._statusbar = wx.GetApp().GetTopWindow().statusbar

        # these are in points to facilitate area calculations
        self.ref_locations = [0,self._dataset.spectral_dims[0]-1 ]  
        
        bob = 10  
        

    # EVENT FUNCTIONS -----------------------------------------------
    
    def on_motion(self, xdata, ydata, val, bounds, iaxis):
        if self.prefs.xaxis_hertz:
            hz  = xdata
            ppm = xdata / self._dataset.frequency
        else:
            ppm = xdata
            hz  = xdata * self._dataset.frequency

        self._statusbar.SetStatusText(" PPM = %.3f" % (ppm, ), 0)
        self._statusbar.SetStatusText(" Hz = %.3f"  % (hz,  ), 1)
        
        # this may give wrong values in an axes that has multiple lines
        # displayed since there is no sure way to determine which line
        # the cursor is closest to
        self._statusbar.SetStatusText((" Value = "+str(val[0])), 2)   
    
    def on_scroll(self, button, step, iaxis):
        self.set_vertical_scale(step)

            
    def on_zoom_motion(self, xmin, xmax, val, ymin, ymax, iplot=None):
        if self.prefs.xaxis_hertz:
            hz_str  = xmin
            hz_end  = xmax
            if hz_str > hz_end: hz_str, hz_end = hz_end, hz_str
            ppm_str = hz_str / self._dataset.frequency
            ppm_end = hz_end / self._dataset.frequency
        else:
            ppm_str  = xmin
            ppm_end  = xmax
            if ppm_str > ppm_end: ppm_str, ppm_end = ppm_end, ppm_str
            hz_str = ppm_str * self._dataset.frequency
            hz_end = ppm_end * self._dataset.frequency
    
        delta_ppm = -1*(ppm_str - ppm_end)  # keeps delta positive
        delta_hz  = delta_ppm * self._dataset.frequency
        self._statusbar.SetStatusText((" PPM Range = %.2f to %.2f" % (ppm_str, ppm_end)), 0)
        self._statusbar.SetStatusText((" Hz Range = %.1f to %.1f" % (hz_str, hz_end)), 1)
        self._statusbar.SetStatusText((" dPPM = %.2f  dHz = %.1f" % (delta_ppm, delta_hz)), 2)


    def on_refs_select(self, xmin, xmax, val, reset=False, iplot=None):
        # Calculate area of span
        area, rms = self.calculate_area()
        txt = ' Area = %1.5g  RMS = %1.5g' % (area[1], rms[1])  # summed values
        self._statusbar.SetStatusText(txt, 3)


    def on_refs_motion(self, xmin, xmax, val, iplot=None):
        if self.prefs.xaxis_hertz:
            hz_str  = xmin
            hz_end  = xmax
            if hz_str > hz_end: hz_str, hz_end = hz_end, hz_str
            ppm_str = hz_str / self._dataset.frequency
            ppm_end = hz_end / self._dataset.frequency
        else:
            ppm_str  = xmin
            ppm_end  = xmax
            if ppm_str > ppm_end: ppm_str, ppm_end = ppm_end, ppm_str
            hz_str = ppm_str * self._dataset.frequency
            hz_end = ppm_end * self._dataset.frequency
                    
        delta_ppm = -1*(ppm_str - ppm_end) # keeps delta positive
        delta_hz  = delta_ppm * self._dataset.frequency
        self._statusbar.SetStatusText((" PPM Range = %.2f to %.2f" % (ppm_str, ppm_end)), 0)
        self._statusbar.SetStatusText((" Hz Range = %.1f to %.1f" % (hz_str, hz_end)), 1)
        self._statusbar.SetStatusText((" dPPM = %.2f  dHz = %.1f" % (delta_ppm, delta_hz)), 2)

        # convert ppm to points
        dataset = self._dataset
        pts_end = (dataset.spectral_dims[0] / 2) - (dataset.frequency * (ppm_str - dataset.resppm) / dataset.spectral_hpp)
        pts_end = np.where(pts_end > 0, pts_end, 0)        
        pts_str = (dataset.spectral_dims[0] / 2) - (dataset.frequency * (ppm_end - dataset.resppm) / dataset.spectral_hpp)
        pts_str = np.where(pts_str > 0, pts_str, 0)        

        # store reference span start/end values as PPM
        self.ref_locations = pts_str, pts_end

        area, rms = self.calculate_area()
        txt = ' Area = %1.5g  RMS = %1.5g' % (area[1], rms[1])  # summed values
        self._statusbar.SetStatusText(txt, 3)

