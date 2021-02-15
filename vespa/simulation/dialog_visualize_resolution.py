# Python modules


# 3rd party modules
import wx

# Our modules
import vespa.simulation.constants as constants
import vespa.simulation.auto_gui.visualize_resolution as visualize_resolution
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


class DialogVisualizeResolution(visualize_resolution.MyDialog):

    def __init__(self, parent, sw, pts):
        
        if not parent:
            parent = wx.GetApp().GetTopWindow()
        
        visualize_resolution.MyDialog.__init__(self, parent, wx.ID_ANY)

        # We add the OK & Cancel buttons dynamically so that they're in the 
        # right order under OS X, GTK, Windows, etc.
        self.ButtonOk, self.ButtonCancel = \
                    wx_util.add_ok_cancel(self, self.LabelOkCancelPlaceholder,
                                          self.on_ok)
        
        self.points      = pts
        self.sweep_width = sw
        
        self.FloatSpectralPoints.multiplier = 2
        self.FloatSpectralPoints.SetDigits(0)
        self.FloatSpectralPoints.SetIncrement(1.0)
        self.FloatSpectralPoints.SetRange(constants.SPECTRAL_POINTS_MIN,constants.SPECTRAL_POINTS_MAX)
        self.FloatSpectralPoints.SetValue(float(pts))
    
        self.FloatSweepWidth.SetDigits(1)
        self.FloatSweepWidth.SetIncrement(500.0)
        self.FloatSweepWidth.SetRange(constants.SWEEP_WIDTH_MIN,constants.SWEEP_WIDTH_MAX)
        self.FloatSweepWidth.SetValue(float(sw))
        
        self.SetTitle("Basis Spectra Resolution")
        
        self.Fit()
        self.Layout()
        self.Center()

        # Under Ubuntu, the focus is nowhere to be found unless I set it
        # explicitly.
        self.FloatSpectralPoints.SetFocus()



    ##### Event Handlers ######################################################
    
    def on_sweep_width(self, event): 
        
        self.sweep_width = event.GetEventObject().GetValue()


    def on_spectral_points(self, event): 
        # make sure that the points selected are a power of 2 so that the
        # fft routines are happy.  Also, basis functions need to be 
        # re-calculated here because their array size has changed
        val = event.GetEventObject().GetValue()
        pow2 = 0
        while val > 2**pow2: pow2 = pow2 + 1
        val = 2**pow2
        event.GetEventObject().SetValue(float(val))

        self.points = event.GetEventObject().GetValue()

    
    def on_ok(self, event):
        self.EndModal(wx.ID_OK)        


    ##### Internal helper functions  ##########################################



