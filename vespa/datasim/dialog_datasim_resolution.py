# Python modules


# 3rd party modules
import wx

# Our modules
import vespa.datasim.auto_gui.datasim_resolution as datasim_resolution
import vespa.common.wx_gravy.util as wx_util



#------------------------------------------------------------------------------
# Note. GUI Architecture/Style
#
# Many of the GUI components in Vespa are designed using the WxGlade 
# application to speed up development times. The GUI components are designed
# interactively and users can preview the resultant window/panel/dialog, but
# while event functions can be specified, only stub functions with those 
# names are created. The WxGlade files (with *.wxg extensions) are stored in
# the 'wxglade' subdirectory. The ouput of their code generation are stored
# in the 'auto_gui' subdirectory. 
#
# To used these GUI classes, each one is inherited into a unique 'vespa'
# class, where program specific initialization and other functionality are
# written.  Also, the original stub functions for widget event handlers are 
# overloaded to provide program specific event handling.
#------------------------------------------------------------------------------

FREQUENCY_MIN = 1.0
FREQUENCY_MAX = 1000.0
SPECTRAL_POINTS_MIN = 32
SPECTRAL_POINTS_MAX = 32768
SWEEP_WIDTH_MIN = 100.0
SWEEP_WIDTH_MAX = 20000.0

class DialogDatasimResolution(datasim_resolution.MyDialog):

    def __init__(self, parent, sw, pts, freq, ppm_start, ppm_end):
        
        if not parent:
            parent = wx.GetApp().GetTopWindow()
        
        datasim_resolution.MyDialog.__init__(self, parent, wx.ID_ANY)

        # Add OK & Cancel here so they're in right order for OS X, GTK, Windows, etc.
        r = wx_util.add_ok_cancel(self, self.LabelOkCancelPlaceholder, self.on_ok)
        self.ButtonOk, self.ButtonCancel = r

        self.points          = pts
        self.sweep_width     = sw
        self.frequency       = freq
        self.ppm_range_start = ppm_start        # FIXME bjs - these need!!! to match up with the sw/pts above
        self.ppm_range_end   = ppm_end

        self.FloatFrequency.SetDigits(3)
        self.FloatFrequency.SetIncrement(1.0)
        self.FloatFrequency.SetRange(FREQUENCY_MIN,FREQUENCY_MAX)
        self.FloatFrequency.SetValue(float(freq))
        
        self.FloatSpectralPoints.multiplier = 2
        self.FloatSpectralPoints.SetDigits(0)
        self.FloatSpectralPoints.SetIncrement(1.0)
        self.FloatSpectralPoints.SetRange(SPECTRAL_POINTS_MIN,SPECTRAL_POINTS_MAX)
        self.FloatSpectralPoints.SetValue(float(pts))
    
        self.FloatSweepWidth.SetDigits(1)
        self.FloatSweepWidth.SetIncrement(500.0)
        self.FloatSweepWidth.SetRange(SWEEP_WIDTH_MIN,SWEEP_WIDTH_MAX)
        self.FloatSweepWidth.SetValue(float(sw))

        ppm_min = 4.7 - (self.sweep_width/self.frequency)/2.0
        ppm_max = 4.7 + (self.sweep_width/self.frequency)/2.0

        self.ppm_range_start = ppm_start if ppm_start > ppm_min else ppm_min
        self.ppm_range_end   = ppm_end   if ppm_end   < ppm_max else ppm_max

        self.FloatPpmRangeStart.SetDigits(3)
        self.FloatPpmRangeStart.SetIncrement(0.1) 
        self.FloatPpmRangeStart.SetRange(ppm_min,ppm_max)
        self.FloatPpmRangeStart.SetValue(float(self.ppm_range_start))

        self.FloatPpmRangeEnd.SetDigits(3)
        self.FloatPpmRangeEnd.SetIncrement(0.1)
        self.FloatPpmRangeEnd.SetRange(ppm_min,ppm_max)
        self.FloatPpmRangeEnd.SetValue(float(self.ppm_range_end))
        
        self.SetTitle("Spectral Resolution and Metabolite Inclusion")
        
        self.Fit()
        self.Layout() 
        self.Center()

        # Under Ubuntu, the focus is nowhere unless set explicitly.
        self.FloatSpectralPoints.SetFocus()



    ##### Event Handlers ######################################################

    def on_frequency(self, event): 
        self.frequency = event.GetEventObject().GetValue()

    
    def on_sweep_width(self, event): 
        self.sweep_width = event.GetEventObject().GetValue()


    def on_spectral_points(self, event): 
        # ensure points selected are a power of 2 for fft routine.  Also,
        # re-calculate basis functions here because array size may have changed
        obj = event.GetEventObject()
        val = obj.GetValue()
        pow2 = 0
        while val > 2**pow2: pow2 = pow2 + 1
        val = 2**pow2
        obj.SetValue(float(val))
        self.points = obj.GetValue()


    def on_ppm_range_start(self, event): 
        
        ppmstr = event.GetEventObject().GetValue()
        ppmend = self.ppm_range_end
        
        if ppmstr > ppmend: ppmstr, ppmend = ppmend, ppmstr
        
        self.ppm_range_start = ppmstr
        self.ppm_range_end   = ppmend
        self.FloatPpmRangeStart.SetValue(ppmstr)
        self.FloatPpmRangeEnd.SetValue(ppmend)

    def on_ppm_range_end(self, event): 
        
        ppmend = event.GetEventObject().GetValue()
        ppmstr = self.ppm_range_start
        
        if ppmstr > ppmend: ppmstr, ppmend = ppmend, ppmstr
        
        self.ppm_range_start = ppmstr
        self.ppm_range_end   = ppmend
        self.FloatPpmRangeStart.SetValue(ppmstr)
        self.FloatPpmRangeEnd.SetValue(ppmend)


    def on_ok(self, event):
        self.EndModal(wx.ID_OK)        


    ##### Internal helper functions  ##########################################



