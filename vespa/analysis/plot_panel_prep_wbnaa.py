# Python modules

# 3rd party modules
import wx

# Our modules
from vespa.common.wx_gravy.plot_panel_spectrum import PlotPanelSpectrum
from vespa.common.wx_gravy.plot_panel_points   import PlotPanelPoints

        

class PlotPanelPrepWbnaa(PlotPanelSpectrum):
    
    def __init__(self, parent, tab, tab_dataset, **kwargs):
        '''
        This is a customization of the PlotPanelSpectrum object for use in  
        our application. The basic matplotlib functionality is contained
        in the base class, all we are doing here is overwriting the
        various canvas Event functions to mesh with our application.
        
        parent - the widget to which the PlotPanelSpectrum is directly attached
        
        tab    - the Notebook tab in which the PlotPanel resides, this gives
                 access to tab.block, the processing block object for events
        
        '''
        tab.SizerSplitterWindow.Fit(tab)  # bugfix wxGlade 0.9.6 to 1.0.0

        super().__init__(parent, **kwargs)

        self.tab = tab  
        self.tab_dataset = tab_dataset
        self.top = wx.GetApp().GetTopWindow()
        
        # these are in points to facilitate area calculations
        self.ref_locations = 0,self.tab.block.dims[0]-1     
        
        self.set_color( (255,255,255) )


    # EVENT FUNCTIONS -----------------------------------------------
    
    def on_motion(self, xdata, ydata, val, bounds, iaxis):

        if iaxis in (0,1,2):
            if self.prefs.xaxis_hertz:
                hz  = xdata
                ppm = xdata / self.tab.dataset.frequency
            else:
                ppm = xdata
                hz  = xdata * self.tab.dataset.frequency
    
            self.top.statusbar.SetStatusText( " PPM = %.3f" % (ppm, ), 0)
            self.top.statusbar.SetStatusText( " Hz = %.3f"  % (hz,  ), 1)
            self.top.statusbar.SetStatusText(( " Value = "+str(val[0])), 2)   
        else: 
            self.top.statusbar.SetStatusText( " " , 0)
            self.top.statusbar.SetStatusText( " ", 1)
            self.top.statusbar.SetStatusText(( " Value = "+str(round(val[0]))), 2)   

        
    
    def on_scroll(self, button, step, iaxis):
        if iaxis in (0,1,2):
            self.set_vertical_scale(step)              
            self.tab_dataset.FloatScale.SetValue(self.vertical_scale)
            
            
    def on_zoom_select(self, xmin, xmax, val, ymin, ymax, reset=False, iplot=None):
        if iplot in (0,1,2):
            self.tab_dataset.FloatScale.SetValue(self.dataymax)


    def on_zoom_motion(self, xmin, xmax, val, ymin, ymax, iplot=None):
        
        if iplot in (0,1,2):
            if self.prefs.xaxis_hertz:
                hz_str  = xmin
                hz_end  = xmax
                if hz_str > hz_end: hz_str, hz_end = hz_end, hz_str
                ppm_str = hz_str / self.tab.dataset.frequency
                ppm_end = hz_end / self.tab.dataset.frequency
            else:
                ppm_str  = xmin
                ppm_end  = xmax
                if ppm_str > ppm_end: ppm_str, ppm_end = ppm_end, ppm_str
                hz_str = ppm_str * self.tab.dataset.frequency
                hz_end = ppm_end * self.tab.dataset.frequency
    
            delta_ppm = -1*(ppm_str - ppm_end)  # keeps delta positive
            delta_hz  = delta_ppm * self.tab.dataset.frequency
            self.top.statusbar.SetStatusText(( " PPM Range = %.2f to %.2f" % (ppm_str, ppm_end)), 0)
            self.top.statusbar.SetStatusText(( " Hz Range = %.1f to %.1f" % (hz_str, hz_end)), 1)
            self.top.statusbar.SetStatusText(( " dPPM = %.2f  dHz = %.1f" % (delta_ppm, delta_hz)), 2)
            
        
    def on_refs_select(self, xmin, xmax, val, reset=False, iplot=None):
        if iplot in (0,1,2):
            # Calculate area of span
            area, rms = self.calculate_area()
            if self.prefs.area_calc_plot_a:
                index = 0
            else:
                index = 1
            self.top.statusbar.SetStatusText(self.tab.build_area_text(area[index], rms[index]), 3)


    def on_refs_motion(self, xmin, xmax, val, iplot=None):

        if iplot in (0,1,2):
            if self.prefs.xaxis_hertz:
                hz_str  = xmin
                hz_end  = xmax
                if hz_str > hz_end: hz_str, hz_end = hz_end, hz_str
                ppm_str = hz_str / self.tab.dataset.frequency
                ppm_end = hz_end / self.tab.dataset.frequency
            else:
                ppm_str  = xmin
                ppm_end  = xmax
                if ppm_str > ppm_end: ppm_str, ppm_end = ppm_end, ppm_str
                hz_str = ppm_str * self.tab.dataset.frequency
                hz_end = ppm_end * self.tab.dataset.frequency
            delta_ppm = -1*(ppm_str - ppm_end)
            delta_hz  = delta_ppm * self.tab.dataset.frequency
            self.top.statusbar.SetStatusText(( " PPM Range = %.2f to %.2f" % (ppm_str, ppm_end)), 0)
            self.top.statusbar.SetStatusText(( " Hz Range = %.1f to %.1f" % (hz_str, hz_end)), 1)
            self.top.statusbar.SetStatusText(( " dPPM = %.2f  dHz = %.1f" % (delta_ppm, delta_hz)), 2)
    
            area, rms = self.calculate_area()
            if self.prefs.area_calc_plot_a:
                index = 0
            else:
                index = 1
            
            self.top.statusbar.SetStatusText(self.tab.build_area_text(area[index], rms[index]), 3)

        

    def on_middle_select(self, xstr, ystr, xend, yend, iplot):
        if iplot in (0,1,2):
            self.tab.process_and_plot()
        else:
            print('boing 11')


    def on_middle_motion(self, xcur, ycur, xprev, yprev, iplot):

        if iplot in (0,1,2):
            # note. the use of phase in FIDSUM is different than in other parts of
            #  the Analysis program. Here it is a part of the processing, thus we
            #  do NOT use the inherent phasing in the plot_panel_spectrum. Rather
            #  phase is processed into each fid and summed to apply add/subtract
            #  effects into the summed spectrum (plot B). Thus we do not call the
            #  built in phasing method in this event handler.

            # The mouse probably moved in both the X and Y directions, but to
            # make phasing easier for the user to do accurately, we only pay
            # attention to the larger of the two movements.
    
            dy = ycur - yprev
            #dx = xcur - xprev

            index = self.tab.fid_index

            if iplot == 0:
                self.tab.block.phase_0[index] += dy
#                self.tab.block.phase_0[index] = self.tab.block.phase_0[index] % 360
            elif iplot == 1 or iplot == 2:
                self.tab.block.phase_0 += dy
#                self.tab.block.phase_0 = self.tab.block.phase_0 % 360

# not going to do Phase1 interactive for now, bjs, bust save for future. - not working completely.
#
#             if abs(dy) > abs(dx):
#                 # 0 order phase
#                 if iplot == 0:
#                     self.tab.block.phase_0[index] += dy
#                     #self.tab.block.phase_0[index] = self.tab.block.phase_0[index] % 360
#                 elif iplot == 1 or iplot == 2:
#                     self.tab.block.phase_0 += dy
#                     #self.tab.block.phase_0 = self.tab.block.phase_0 % 360
#             else:
#                 delta = dx*10
#                 self.tab.block.set.global_phase1 += dx
    
            val = self.tab.block.phase_0[index]
            self.tab.FloatPhase0.SetValue(val)
            #val = self.tab.block.set.global_phase1
            #self.tab.FloatPhase1.SetValue(val)
     
            self.tab.process_and_plot()
            
            # Calculate the new area after phasing
            area, rms = self.calculate_area()
            if self.prefs.area_calc_plot_a:
                area = area[0]
                rms = rms[0]
            else:
                area = area[1]
                rms = rms[1]
            
            self.top.statusbar.SetStatusText(self.tab.build_area_text(area, rms), 3)



class PlotPanelPrepWbnaaSeries(PlotPanelPoints):
    
    def __init__(self, parent, tab, tab_dataset, **kwargs):
        '''
        This is a customization of the PlotPanelPoints object for use in  
        our application. The basic matplotlib functionality is contained
        in the base class, all we are doing here is overwriting the
        various canvas Event functions to mesh with our application.
        
        parent - the widget to which the PlotPanelSpectrum is directly attached
        
        tab    - the Notebook tab in which the PlotPanel resides, this gives
                 access to tab.block, the processing block object for events
        
        '''
        
        PlotPanelPoints.__init__( self, parent, **kwargs )

        self.tab = tab  
        self.tab_dataset = tab_dataset
        self.top = wx.GetApp().GetTopWindow()
        
        # these are in points to facilitate area calculations
        self.ref_locations = 0,self.tab.block.dims[0]-1     
        
        self.set_color( (255,255,255) )


    # EVENT FUNCTIONS -----------------------------------------------
    
    def on_motion(self, xdata, ydata, val, bounds, iaxis):  

        self.top.statusbar.SetStatusText( " X Loc = %.3f" % (xdata, ), 0)
        self.top.statusbar.SetStatusText( " " , 1)
        self.top.statusbar.SetStatusText(( " Value = "+str(val[0])), 2)   
       
    
    def on_scroll(self, button, step, iaxis):
        self.set_vertical_scale(step)              
        self.tab_dataset.FloatScale.SetValue(self.vertical_scale)

            
    def on_zoom_select(self, xmin, xmax, val, ymin, ymax, reset=False, iplot=None):
        pass


    def on_zoom_motion(self, xmin, xmax, val, ymin, ymax, iplot=None):

        delta_x  = xmax - xmin
        self.top.statusbar.SetStatusText(( " X Range = %.2f to %.2f" % (xmin, xmax)), 0)
        self.top.statusbar.SetStatusText( " " , 1)
        self.top.statusbar.SetStatusText(( " delta X = %.2f" % (delta_x,)), 2)
            
        
    def on_refs_select(self, xmin, xmax, val, reset=False, iplot=None, rindex=None):
        pass


    def on_refs_motion(self, xmin, xmax, val, iplot=None):

        delta_x = xmax - xmin
        self.top.statusbar.SetStatusText(( " X Range = %.2f to %.2f" % (xmin, xmax)), 0)
        self.top.statusbar.SetStatusText( " " , 1)
        self.top.statusbar.SetStatusText(( " delta X = %.2f" % (delta_x,)), 2)
        

    def on_middle_select(self, xstr, ystr, xend, yend, iplot):
        pass


    def on_middle_motion(self, xcur, ycur, xprev, yprev, iplot):
        pass
        

    def on_middle_press(self, xloc, yloc, iplot, bounds=None, xdata=None, ydata=None):

        if xdata is not None:
            xindex = int(round(xdata))
            self.tab.block.toggle_exclude_index(self.tab.dataset, xindex)
            val = ','.join(str(x) for x in self.tab.block.exclude_indices)
            self.tab.TextDataExclusion.SetValue(val)
            self.tab.TextDataExclusion.SetValue(val)
            
            self.tab.fid_index = xindex 
            self.tab.SpinFidIndex.SetValue(self.tab.fid_index)
            shft = self.tab.block.frequency_shift[self.tab.fid_index]
            phas = self.tab.block.phase_0[self.tab.fid_index]
            self.tab.FloatPeakShiftValue.SetValue(shft)        
            self.tab.FloatPhase0.SetValue(phas)
            self.tab.FloatPhase1.SetValue(self.tab.block.set.global_phase1)
            
            self.tab.process_and_plot()   
