# Python modules


# 3rd party modules
import wx

# Our modules
import vespa.common.wx_gravy.plot_panel_spectrum as plot_panel_spectrum



class PlotPanelWatref(plot_panel_spectrum.PlotPanelSpectrum):

    def __init__(self, parent, tab, tab_dataset, **kwargs):
        '''
        This is a customization of the PlotPanel object for use in a specific
        location in our application. The basic matplotlib functionality is
        contained in the base class, all we are doing here is overwriting the
        various canvas Event functions to mesh with our application.

        parent      - the widget to which the PlotPanel is directly attached
        tab         - the Spectral tab in which the PlotPanel resides
        tab._dataset - the Dataset tab in which the Spectral tab resides

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
        if self.tab.block.set.fft:
            if self.prefs.xaxis_hertz:
                hz  = xdata
                ppm = xdata / self.tab.dataset.frequency
            else:
                ppm = xdata
                hz  = xdata * self.tab.dataset.frequency

            self.top.statusbar.SetStatusText( " PPM = %.3f" % (ppm, ), 0)
            self.top.statusbar.SetStatusText( " Hz = %.3f"  % (hz,  ), 1)
        else:
            self.top.statusbar.SetStatusText( " Time [ms] = %.2f" % (xdata, ), 0)
            self.top.statusbar.SetStatusText( " " , 1)

        self.top.statusbar.SetStatusText(( " Value = "+str(val[0])), 2)


    def on_scroll(self, button, step, iaxis):
        self.set_vertical_scale(step)
        self.tab_dataset.FloatScale.SetValue(self.vertical_scale)


    def on_zoom_select(self, xmin, xmax, val, ymin, ymax, reset=False, iplot=None):
        self.tab_dataset.FloatScale.SetValue(self.dataymax)


    def on_zoom_motion(self, xmin, xmax, val, ymin, ymax, iplot=None):
        if self.tab.block.set.fft:
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
            self.top.statusbar.SetStatusText(( " Hz Range = %.1f to %.1f"  % (hz_str, hz_end)), 1)
            self.top.statusbar.SetStatusText(( " dPPM = %.2f  dHz = %.1f"  % (delta_ppm, delta_hz)), 2)
        else:
            self.top.statusbar.SetStatusText(( " Time Range [ms] = %.2f to %.2f" % (xmin, xmax)), 0)
            self.top.statusbar.SetStatusText(  " " , 1)
            self.top.statusbar.SetStatusText(( " dTime [ms] = %.2f" % (xmax-xmin, )), 2)


    def on_refs_select(self, xmin, xmax, val, reset=False, iplot=None):
        # Calculate area of span
        area, rms = self.calculate_area()
        if self.prefs.area_calc_plot_a:
            index = 0
            labl = 'A'
        elif self.prefs.area_calc_plot_b:
            index = 1
            labl = 'B'
        elif self.prefs.area_calc_plot_c:
            index = 2
            labl = 'C'

        self.top.statusbar.SetStatusText(self.tab.build_area_text(area[index], rms[index], plot_label=labl), 3)


    def on_refs_motion(self, xmin, xmax, val, iplot=None):
        if self.tab.block.set.fft:
            if self.prefs.xaxis_hertz:
                hz_str  = xmin
                hz_end  = xmax
                if hz_str > hz_end: hz_str, hz_end = hz_end, hz_str
                ppm_str = hz_str / self.tab.dataset.frequency
                ppm_end = hz_end /self.tab.dataset.frequency
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
        else:
            self.top.statusbar.SetStatusText(( " Time Range [ms] = %.2f to %.2f" % (xmin, xmax)), 0)
            self.top.statusbar.SetStatusText(  " " , 1)
            self.top.statusbar.SetStatusText(( " dTime [ms] = %.2f" % (xmax-xmin, )), 2)

        # Calculate area of span
        area, rms = self.calculate_area()
        if self.prefs.area_calc_plot_a:
            index = 0
            labl = 'A'
        elif self.prefs.area_calc_plot_b:
            index = 1
            labl = 'B'
        elif self.prefs.area_calc_plot_c:
            index = 2
            labl = 'C'
        self.top.statusbar.SetStatusText(self.tab.build_area_text(area[index], rms[index], plot_label=labl), 3)


    def on_middle_select(self, xstr, ystr, xend, yend, iplot):
        pass


    def on_middle_motion(self, xcur, ycur, xprev, yprev, iplot):
        '''
        Interactive phasing occurs when the user holds down the middle button
        and moves around.  Up/down for phase0 and left/right for phase1.

        The actual calls to plot_panel_spectrum.set_phase_0() or phase_1()
        occur at the dataset level.  That is because phase0/1 can be changed
        from multiple locations (main spectral plot interactive, svd filter
        plot interactive, voigt tab plot interactive, main spectral phase 0
        and phase 1 widgets, voigt initial value widgets. (Note. this is
        true for B0 shift as well by not dealt with here). In tab._dataset.py,
        the set_phase_0() and set_phase_1() methods can be called from any
        tab without sibling tabs talking to each other.

        Due to the 'sync' check box on the Spectral tab, phase changes on
        one dataset tab can also cause phase changes on another dataset tab.
        To accomodate this rule, phase events actually call an event handler
        at the notebook_dataset level. This polls whether sync is on, whether
        there is a datasetB displayed that needs updating too, and iterates
        through all dataset tabs to find out which tabs are affected by this
        phase event. When it has a list of affected datasets, it then calls
        one by one the set_phase_x() calls in each dataset to update the
        block values and view plots.

        An important differentiation to understand is that now that phase is
        inherently taken care of in the plot_panel_spectral object, the phase
        value stored in the block AND the new phase displayed in various plots
        have to be explicitly applied.  It is no longer sufficient to just set
        the phase0/1 in the block and then call process_and_plot to get the
        plots to update.  We go to this trouble in order to gain flexibility
        and speed in the interactive phasing and general display of data in
        each plot_panel_spectrum derived class.

        '''
        if iplot not in (0,1): return
        voxel = self.tab_dataset.voxel

        # The mouse probably moved in both the X and Y directions, but to
        # make phasing easier for the user to do accurately, we only pay
        # attention to the larger of the two movements.

        dy = ycur - yprev
        dx = xcur - xprev

        # determine label(s) of dataset(s) being phased
        if self.tab.do_sync:
            poll_labels = [self.tab_dataset.indexAB[0],self.tab_dataset.indexAB[1]]
        else:
            if iplot==0:
                poll_labels = [self.tab_dataset.indexAB[0]]
            else:
                poll_labels = [self.tab_dataset.indexAB[1]]

        index = [iplot]
        if self.tab.do_sync:
            index = [0,1]

        if abs(dy) > abs(dx):
            # 0 order phase
            delta = dy
            self.tab.top.notebook_datasets.global_poll_phase(poll_labels, delta, voxel, do_zero=True)
        else:
            pass
            # first order phase, x10 makes Phase1 changes more noticeable
            delta = dx*10
            self.tab.top.notebook_datasets.global_poll_phase(poll_labels, delta, voxel, do_zero=False)

        # Calculate the new area after phasing
        area, rms = self.calculate_area()
        if self.prefs.area_calc_plot_a:
            index = 0
            labl = 'A'
        elif self.prefs.area_calc_plot_b:
            index = 1
            labl = 'B'
        elif self.prefs.area_calc_plot_c:
            index = 2
            labl = 'C'
        self.top.statusbar.SetStatusText(self.tab.build_area_text(area[index], rms[index], plot_label=labl), 3)


