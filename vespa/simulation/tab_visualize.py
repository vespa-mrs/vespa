# Python modules

import copy
import math
import os
import tempfile

# 3rd party modules
import wx
import numpy as np
import wx.lib.agw.aui as aui        # NB. wx.aui version throws odd wxWidgets exception on Close/Exit ?? Not anymore in wxPython 4.0.6 ??
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm


# Our modules
import vespa.simulation.prefs as prefs
import vespa.simulation.constants as constants
import vespa.simulation.util_menu as util_menu
import vespa.simulation.tab_visualize as tab_visualize
import vespa.simulation.tab_simulate as tab_simulate
import vespa.simulation.plot_panel_plot1d as plot_panel_plot1d
import vespa.simulation.plot_panel_integral as plot_panel_integral
import vespa.simulation.plot_panel_contour as plot_panel_contour
import vespa.simulation.mrs_data_basis as mrs_data_basis
import vespa.simulation.build_basis_functions as bbf_module
import vespa.common.util.ppm as util_ppm
import vespa.common.util.generic_spectral as util_generic_spectral
import vespa.common.wx_gravy.util as wx_util
import vespa.common.constants as common_constants
import vespa.common.mrs_experiment as mrs_experiment
import vespa.common.wx_gravy.common_dialogs as common_dialogs
import vespa.common.wx_gravy.notebooks as vespa_notebooks

import vespa.simulation.constants as constants
import vespa.simulation.dialog_visualize_resolution as dialog_visualize_resolution
import vespa.simulation.auto_gui.visualize as visualize

PI = math.pi


#------------------------------------------------------------------------------
# Note. GUI Architecture/Style
#
# Most Vespa GUI components are designed using WxGlade for rapid developement.
# WxGlade is interactive and users can preview the resultant window/panel/dialog.
# However, event functions can be specified but only stub functions with those
# names are created. WxGlade (*.wxg) files are in the 'wxglade' subdirectory.
# WxGlade generated modules are stored in 'auto_gui' subdirectory.
#
# Each WxGlade created GUI class is inherited into a unique 'vespa' class, where
# program specific initialization and other functionality are written. The stub
# event functions are overloaded here to provide program specific event handling.
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------


def _configure_combo(control, choices, selection=''):
        lines = list(choices.values())
        control.SetItems(lines)
        if selection in lines:
            control.SetStringSelection(selection)
        else:
            control.SetStringSelection(lines[0])


def _paired_event(obj_min, obj_max):
        val_min = obj_min.GetValue()
        val_max = obj_max.GetValue()
        pmin = min(val_min, val_max)
        pmax = max(val_min, val_max)
        obj_min.SetValue(pmin)
        obj_max.SetValue(pmax)
        return pmin, pmax


class TabVisualize(visualize.PanelVisualizeUI):

    def __init__(self, tab_experiment, top, db, experiment, prefs):

        visualize.PanelVisualizeUI.__init__(self, tab_experiment)

        self._top = wx.GetApp().GetTopWindow()
        self._tab_experiment = tab_experiment
        self._prefs = prefs

        self.db = db
        self.top = top
        self.experiment = experiment
        self.statusbar = top.statusbar

        # metabolite model attributes

        d = { "frequency" : self.experiment.b0 }
        self.basis  = mrs_data_basis.DataBasis(d)
        self.freq1d = mrs_data_basis.DataBasis(d)
        self.freq2d = mrs_data_basis.DataBasis(d)

        self.maxppm         = self.freq1d.pts2ppm(0)
        self.minppm         = self.freq1d.pts2ppm(self.freq1d.dims[0]-1)
        self.index1         = 0
        self.index2         = 0
        self.index3         = 0

        self.linewidth      = common_constants.DEFAULT_LINEWIDTH
        self.apodization    = None
        self.integral_data  = None
        self.plot_labels    = None

        self.initialize_controls()
        self.populate_controls()
        self.build_basis_functions()
        self.plot_canvas()

        # Grab key events so I can look for Ctrl/Cmd+A ==> select all
        self.ListPlotMetabolites.Bind(wx.EVT_KEY_DOWN, self.on_key_down)

        # If the sash position isn't recorded in the INI file, we use the
        # arbitrary-ish value of 400.
        if not self._prefs.sash_position:
            self._prefs.sash_position = 400

        # Under OS X, wx sets the sash position to 10 (why 10?) *after*
        # this method is done. So setting the sash position here does no
        # good. We use wx.CallAfter() to (a) set the sash position and
        # (b) fake an EVT_SPLITTER_SASH_POS_CHANGED.
        wx.CallAfter(self.SplitterWindow.SetSashPosition, self._prefs.sash_position, True)
        wx.CallAfter(self.on_splitter)


    # =======================================================
    #
    #           GUI Setup Calls
    #
    # =======================================================

    def initialize_controls(self):
        """
        Initializes the controls to be the right size or have the right
        range or number of decimal places. It typically does not set the
        default value (that's for populate_controls method to do). This
        method does the one-time setup bits.

        """

        #-------------------------------------------------------------
        # Set up the view tabs

        self.AuiNotebookSizer = self.TextAuiPlaceholder.GetContainingSizer()
        self.TextAuiPlaceholder.Destroy()

        style =  aui.AUI_NB_TAB_SPLIT | \
                 aui.AUI_NB_TAB_MOVE | \
                 aui.AUI_NB_TAB_EXTERNAL_MOVE | \
                 wx.NO_BORDER | \
                 aui.AUI_NB_BOTTOM

        self.tab_notebook = aui.AuiNotebook(self.PanelAuiNotebook, agwStyle=style)

        self.panel_plot1d = wx.Panel(self.tab_notebook, -1)
        self.plot1d = plot_panel_plot1d.PlotPanelPlot1d(self.panel_plot1d,
                                                        self,
                                                        self._prefs,
                                                        naxes=1,
                                                        reversex=True,
                                                        zoom='span',
                                                        reference=True,
                                                        do_zoom_select_event=True,
                                                        do_zoom_motion_event=True,
                                                        do_refs_select_event=True,
                                                        do_refs_motion_event=True,
                                                        uses_collections=True,
                                                        props_zoom=dict(alpha=0.2, facecolor='yellow'),
                                                        props_cursor=dict(alpha=0.1, facecolor='gray'))

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.plot1d, 1, wx.LEFT | wx.TOP | wx.EXPAND)
        self.panel_plot1d.SetSizer(sizer)
        self.plot1d.Fit()

        self.panel_integral = wx.Panel(self.tab_notebook, -1)
        self.integral = plot_panel_integral.PlotPanelIntegral(self.panel_integral,
                                                              self,
                                                              naxes=1,
                                                              reversex=False,
                                                              zoom='span',
                                                              reference=True,
                                                              middle=True,
                                                              do_zoom_select_event=False,
                                                              do_zoom_motion_event=True,
                                                              do_refs_select_event=False,
                                                              do_refs_motion_event=True,
                                                              do_middle_press_event=True,
                                                              props_zoom=dict(alpha=0.2, facecolor='yellow'),
                                                              props_cursor=dict(alpha=0.1, facecolor='gray'))

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.integral, 1, wx.LEFT | wx.TOP | wx.EXPAND)
        self.panel_integral.SetSizer(sizer)
        self.integral.Fit()

        self.panel_contour = wx.Panel(self.tab_notebook, -1)
        self.contour = plot_panel_contour.PlotPanelContour(self.panel_contour,
                                                           self,
                                                           naxes=1,
                                                           reversex=False,
                                                           zoom='box',
                                                           reference=False,
                                                           middle=True,
                                                           do_zoom_select_event=False,
                                                           do_zoom_motion_event=True,
                                                           do_refs_select_event=False,
                                                           do_refs_motion_event=False,
                                                           do_middle_press_event=True,
                                                           props_zoom=dict(alpha=0.2, facecolor='yellow'),
                                                           props_cursor=dict(alpha=0.1, facecolor='gray'))

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.contour, 1, wx.LEFT | wx.TOP | wx.EXPAND)
        self.panel_contour.SetSizer(sizer)
        self.contour.Fit()

        #------------------------------------------------------------------
        # Insert all panels into AUI Notebook and decorate as needed

        # The name of the stock tab art changed somewhere between wx 2.8 and 3.0.
        if hasattr(aui, 'AuiDefaultTabArt'):
            tab_art_provider = aui.AuiDefaultTabArt
        else:
            tab_art_provider = aui.AuiGenericTabArt
        self.tab_notebook.SetArtProvider(tab_art_provider())
        self.tab_notebook.AddPage(self.panel_plot1d," 1D Plot ",False)
        self.tab_notebook.AddPage(self.panel_integral," Integral Plot ",False)
        self.tab_notebook.AddPage(self.panel_contour," Contour Plot ",False)

        self.AuiNotebookSizer.Add(self.tab_notebook, 1, wx.EXPAND, 0)

        # ---------------------------------------------------------------------
        # plot initializations

        xx = np.arange(self.freq1d.dims[0], dtype='float64')
        tt = xx / self.freq1d.sw
        self.apodization = util_generic_spectral.apodize(tt, self.linewidth, 'Gaussian')

        xx = self.freq1d.pts2ppm(xx)

        self.plot1d.figure.set_facecolor('white')
        self.plot1d_axes = self.plot1d.axes[0]
        self.plot1d_axes.set_facecolor(self._prefs.bgcolor)
        self.plot1d_axes.set_xlim(xx.max(),xx.min())
        self.plot1d_axes.axes.get_yaxis().set_visible(False)
        if not self._prefs.xaxis_show:
            self.plot1d_axes.axes.set_xticklabels([])

        self.integral.figure.set_facecolor('white')
        self.integral_axes = self.integral.axes[0]
        self.integral_axes.set_facecolor(self._prefs.bgcolor)
        self.integral_axes.set_xlim(0,1)
        if not self._prefs.integral_xaxis_show:
            self.integral_axes.axes.set_xticklabels([])
        if not self._prefs.integral_yaxis_show:
            self.integral_axes.axes.set_yticklabels([])

        self.contour.figure.set_facecolor('white')
        self.contour_axes = self.contour.axes[0]
        self.contour_axes.set_facecolor(self._prefs.bgcolor)
        self.contour_axes.set_xlim(0,1)
        self.contour_axes.set_ylim(0,1)
        if not self._prefs.contour_axes_show:
            self.contour_axes.axes.set_xticklabels([])
            self.contour_axes.axes.set_yticklabels([])


    def update_plot_controls(self):
        """ Legacy method name ... deprecated. """
        self.populate_controls()


    def populate_controls(self):

        # This update routine is located on the Experiment_Tab level because
        # changes to the Simulate Tab can affect the Visualize Tab, but not
        # the other way around. So all Simulate functions are in the Simulate
        # module, but not all Visualize functions are in Visualize module.

        # update the Visualize tab controls

        self.PanelVisualize.Enable(True)
        expt = self.experiment

        # a few calculations

        npts = self.freq1d.dims[0]
        sw = self.freq1d.sw
        field = expt.b0
        hpp = sw / npts

        self.maxppm = self.freq1d.pts2ppm(0)
        self.minppm = self.freq1d.pts2ppm(npts - 1)

        # update Plot control

        # There are N+1 display modes, where N = the number of non-empty
        # dims in the experiment. "1D plot" is always an option.
        # Before repopulating this list of choices, we grab the current
        # selection so we can restore it if possible.
        selection = self.ChoiceDisplayMode.GetStringSelection()
        choices = ["1D Plot"]
        # The other options are a stack plot for each non-empty dim.
        choices += ["Stack Plot Index %d" % (i + 1) for i, dim
                    in enumerate(self.experiment.dims)
                    if dim != mrs_experiment.DEFAULT_LOOP]
        self.ChoiceDisplayMode.SetItems(choices)

        # Restore the selection if it's still in the combobox.
        selection = self.ChoiceDisplayMode.FindString(selection)
        if selection == wx.NOT_FOUND:
            selection = 0
        self.ChoiceDisplayMode.SetSelection(selection)

        # setting the Xaxis and Cursor ranges depends on the initial setup
        # of the Visualize panel to have urealistically large numbers in
        # these fields (ie. -999 to 999)

        xaxmax = self.FloatXaxisMax.GetValue()
        xaxmin = self.FloatXaxisMin.GetValue()
        curmax = self.FloatCursorMax.GetValue()
        curmin = self.FloatCursorMin.GetValue()

        if xaxmax > self.maxppm: xaxmax = self.maxppm
        if xaxmin > self.minppm: xaxmin = self.minppm
        if curmax > self.maxppm: curmax = self.maxppm
        if curmin > self.minppm: curmin = self.minppm

        #        if self.plot_options.x_axis_type != PlotOptions.X_AXIS_PPM:
        #
        #            #FIXME - BJS found that the following holds
        #            # ppm2hz != pts2hz(ppm2pts())
        #            # hz2ppm != pts2ppm(hz2pts())
        #            # Need to check original IDL code
        #
        #            xaxmax = self.freq1d.pts2hz(self.freq1d.ppm2pts(xaxmax))
        #            xaxmin = self.freq1d.pts2hz(self.freq1d.ppm2pts(xaxmin))
        #            curmax = self.freq1d.pts2hz(self.freq1d.ppm2pts(curmax))
        #            curmin = self.freq1d.pts2hz(self.freq1d.ppm2pts(curmin))
        #            maxppm = self.freq1d.pts2hz(self.freq1d.ppm2pts(self.maxppm))
        #            minppm = self.freq1d.pts2hz(self.freq1d.ppm2pts(self.minppm))
        #        else:
        #            maxppm = self.maxppm
        #            minppm = self.minppm
        maxppm = self.maxppm
        minppm = self.minppm

        self.FloatXaxisMax.SetValue(xaxmax)
        self.FloatXaxisMin.SetValue(xaxmin)
        self.FloatCursorMax.SetValue(curmax)
        self.FloatCursorMin.SetValue(curmin)
        self.FloatXaxisMax.SetRange(minppm, maxppm)
        self.FloatXaxisMin.SetRange(minppm, maxppm)
        self.FloatCursorMax.SetRange(minppm, maxppm)
        self.FloatCursorMin.SetRange(minppm, maxppm)

        for i in range(common_constants.RESULTS_SPACE_DIMENSIONS - 1):
            spin_control = getattr(self, "SpinIndex%d" % (i + 1))
            spin_control.SetRange(1, len(self.experiment.dims[i]))
            # Undex OS X (not Win or GTK), it's necessary to explicitly set
            # the value to get it to display. Setting just the range is not
            # sufficient.
            spin_control.SetValue(1)

        # By default I hide all of the loop control panels and then show
        # them as needed. Rather than hardcoding their names here, I just
        # use getattr until it fails.
        i = 1
        panel = getattr(self, "panel_loop%d_controls" % i)
        while panel:
            panel.Hide()
            i += 1
            try:
                panel = getattr(self, "panel_loop%d_controls" % i)
            except AttributeError:
                panel = None

        # Now show & populate the enabled loop panels (if any).
        if expt.pulse_sequence:
            for i, loop_label in enumerate(expt.pulse_sequence.loop_labels):
                j = i + 1
                panel = getattr(self, "panel_loop%d_controls" % j)
                panel.Show()

                # Show the labels for these loops and their values.
                control = getattr(self, "LabelLoop%dLabel" % j)
                control.SetLabel(loop_label + " = ")
                control = getattr(self, "LabelIndex%dValue" % j)
                control.SetLabel(str(self.experiment.dims[i][0]))
                # panel.Layout()

        # Deselect everything in the PlotMetabs listbox before resetting its
        # contents. This is a fix/workaround 
        i = self.ListPlotMetabolites.GetCount()
        while i:
            i -= 1
            self.ListPlotMetabolites.Deselect(i)

        names = [metabolite.name for metabolite in expt.metabolites]
        self.ListPlotMetabolites.SetItems(names)

        self.CheckGrayscale.SetValue(self._prefs.contour_grayscale)
        self.SpinContourLevels.SetValue(self._prefs.contour_levels)

        if self.linewidth < constants.LINEWIDTH_MIN: self.linewidth = constants.LINEWIDTH_MIN
        if self.linewidth > constants.LINEWIDTH_MAX: self.linewidth = constants.LINEWIDTH_MAX
        if self.freq1d.sw < constants.SWEEP_WIDTH_MIN: self.freq1d.sw = constants.SWEEP_WIDTH_MIN
        if self.freq1d.sw > constants.SWEEP_WIDTH_MAX: self.freq1d.sw = constants.SWEEP_WIDTH_MAX
        if self.freq1d.dims[0] < constants.SPECTRAL_POINTS_MIN: self.freq1d.dims[0] = constants.SPECTRAL_POINTS_MIN
        if self.freq1d.dims[0] > constants.SPECTRAL_POINTS_MAX: self.freq1d.dims[0] = constants.SPECTRAL_POINTS_MAX

        self.FloatLinewidth.SetValue(float(self.linewidth))
        self.FloatLinewidth.SetDigits(1)
        self.FloatLinewidth.SetIncrement(0.5)
        self.FloatLinewidth.SetRange(constants.LINEWIDTH_MIN, constants.LINEWIDTH_MAX)

        self.PanelVisualize.Layout()
        # self.PanelVisualize.Refresh()
        # self.PanelVisualize.Update()

        # Oddly enough, .Refresh() and .Update() don't make visible the
        # controls hidden/shown above (or if it makes them visible, it
        # positions them incorrectly). So we resort to a trick to force
        # the pane to reconsider the layout of its contents.
        # Resizing is a good way to do that. Under Windows & OS X it's
        # sufficient to subtract then add one to the pane's width. The net
        # change is zero, but it has the desired effect.
        # Under GTK, however, the window size must actually change. If we
        # always used the same delta (+1 or -1), then the pane would gradually
        # get larger or smaller every time the user selects a new water
        # filter.
        # To avoid this, we vary the delta based on whether or not the
        # current pane width is odd or even. No matter how many times the
        # user selects from the water filter list, the size will vary at
        # most one pixel from the original width.
        window = self.top._mgr.GetPane("experiments").window

        size = window.GetSize()
        size.width += (-1 if (size.width % 2) else 1)
        window.SetSize(size)





    # =======================================================
    #
    #           Event Handlers
    #
    # =======================================================

    def on_key_down(self, event):
        if wx_util.is_select_all(event):
            names = [metabolite.name for metabolite in self.experiment.metabolites]
            for item in names:
                self.ListPlotMetabolites.SetStringSelection(item)
            self.plot_canvas(reset_history=True)
            self.plot_integral()
        else:
            # Don't eat the keystroke! If we eat it, users can't use the
            # keyboard to move around in the list or select items.
            event.Skip()

    def on_display_mode(self, event):
        self.plot_canvas(reset_history=True)
        self.plot_integral()

    def on_xaxis_max(self, event):
        self.update_widget_ranges()
        xmax = self.FloatXaxisMax.GetValue()
        xmin = self.FloatXaxisMin.GetValue()
        if self._prefs.xaxis_hertz:
            freq1d = self.freq1d
            xmax = freq1d.pts2hz(freq1d.ppm2pts(xmax))
            xmin = freq1d.pts2hz(freq1d.ppm2pts(xmin))
        if xmax == xmin:
            xmin *= 0.999
        self.plot1d_axes.set_xlim((xmax,xmin))
        self.plot1d.canvas.draw()

    def on_xaxis_min(self, event):
        self.update_widget_ranges()
        xmax = self.FloatXaxisMax.GetValue()
        xmin = self.FloatXaxisMin.GetValue()
        if self._prefs.xaxis_hertz:
            freq1d = self.freq1d
            xmax = freq1d.pts2hz(freq1d.ppm2pts(xmax))
            xmin = freq1d.pts2hz(freq1d.ppm2pts(xmin))
        if xmax == xmin:
            xmin *= 0.999
        self.plot1d_axes.set_xlim((xmax,xmin))
        self.plot1d.canvas.draw()

    def on_cursor_max(self, event):
        self.update_widget_ranges()
        xmax = self.FloatCursorMax.GetValue()
        xmin = self.FloatCursorMin.GetValue()
        if self._prefs.xaxis_hertz:
            freq1d = self.freq1d
            xmax = freq1d.pts2hz(freq1d.ppm2pts(xmax))
            xmin = freq1d.pts2hz(freq1d.ppm2pts(xmin))
        if xmax == xmin:
            xmin *= 0.999
        if self.plot1d.refs != None:
            self.plot1d.refs.set_span(xmin,xmax)
            self.plot1d.canvas.draw()
            self.plot_integral()
            self.plot_contour()

    def on_cursor_min(self, event):
        self.update_widget_ranges()
        xmax = self.FloatCursorMax.GetValue()
        xmin = self.FloatCursorMin.GetValue()
        if self._prefs.xaxis_hertz:
            freq1d = self.freq1d
            xmax = freq1d.pts2hz(freq1d.ppm2pts(xmax))
            xmin = freq1d.pts2hz(freq1d.ppm2pts(xmin))
        if xmax == xmin:
            xmin *= 0.999
        if self.plot1d.refs != None:
            self.plot1d.refs.set_span(xmin,xmax)
            self.plot1d.canvas.draw()
            self.plot_integral()
            self.plot_contour()

    def on_index(self, event, dim):
         # This pseudo event handler is called by the on_indexN() event
         # handlers (below)
        index = event.GetEventObject().GetValue() - 1
        control = getattr(self, "LabelIndex%dValue" % dim)
        # Under Windows/wx 3.x, SetRange() has a side effect of firing one of these events with
        # index == -1 which raises an IndexError if we don't discard/ignore it.
        if index != -1:
            control.SetLabel(str(self.experiment.dims[dim - 1][index]))
        self.plot_canvas()
        self.plot_integral()

    def on_index1(self, event):
        self.on_index(event, 1)

    def on_index2(self, event):
        self.on_index(event, 2)

    def on_index3(self, event):
        self.on_index(event, 3)

    def on_metabolites(self, event):
        self.plot_canvas(reset_history=True)
        self.plot_integral()
        self.plot_contour()

    def on_sum_plots(self, event):
        self.plot_canvas(reset_history=True)

    def on_grayscale(self, event):
        self._prefs.contour_grayscale = event.GetEventObject().IsChecked()
        self.plot_contour(recalc=False)

    def on_contour_levels(self, event):
        self._prefs.contour_levels = event.GetEventObject().GetValue()
        self.plot_contour(recalc=False)

    def on_contour_mode(self, event):
        index = event.GetEventObject().GetCurrentSelection()

        # Contour plots involve two dimensions, and at least one of the
        # dims must have a length > 1. The code below verifies that that's
        # true. If it's not, we attempt to find and select a contour plot
        # that's valid.
        dim_lengths = [len(dim) for dim in self.experiment.dims]
        if index == 0:
            if (dim_lengths[0] == 1) and (dim_lengths[1] == 1):
                if (dim_lengths[2] != 1):
                    index = 1
                else:
                    self.plot_contour(reset=True)
                    return
        elif index == 1:
            if (dim_lengths[0] == 1) and (dim_lengths[2] == 1):
                if (dim_lengths[1] != 1):
                    index = 0
                else:
                    self.plot_contour(reset=True)
                    return
        elif index == 2:
            if (dim_lengths[1] == 1) and (dim_lengths[2] == 1):
                if (dim_lengths[0] != 1):
                    index = 0
                else:
                    self.plot_contour(reset=True)
                    return
        self._prefs.contour_mode = index
        event.GetEventObject().SetSelection(index)
        self.plot_contour()

    def on_linewidth(self, event):
        self.linewidth = event.GetEventObject().GetValue()
        self.set_apodization()
        self.plot_canvas()
        self.plot_integral()

    def on_resolution(self, event):

        pts = self.freq1d.dims[0]
        sw  = self.freq1d.sw
        dialog = dialog_visualize_resolution.DialogVisualizeResolution(self, sw, pts)
        if dialog.ShowModal() == wx.ID_OK:

            sw  = dialog.sweep_width
            pts = int(dialog.points)

            self.freq1d.sw = sw
            self.freq2d.sw = sw
            self.basis.sw  = sw

            self.build_basis_functions(dim0=pts)
            self.update_widget_ranges()
            self.set_apodization(dim0=pts)
            self.plot_canvas(reset_history=True)
            self.plot_integral()

        dialog.Destroy()

    def on_ascii_display(self, event):
        self.display_results_text()

    def on_splitter(self, event=None):
        # This is sometimes called programmatically, in which case event is None
        self._prefs.sash_position = self.SplitterWindow.GetSashPosition()


    #=======================================================
    #
    #           Internal methods start here
    #
    #=======================================================

    def populate_metabolite_list(self):
        """Clears & repopulates the list of metabs."""
        self.ListPlotMetabolites.Clear()
        metabolites = [self.db.fetch_metabolite(metabolite.id) \
                                for metabolite in self.experiment.metabolites]

        names = [metabolite.name for metabolite in metabolites]
        self.ListPlotMetabolites.SetItems(names)


    def update_widget_ranges(self):
        self.maxppm = self.freq1d.pts2ppm(0)
        self.minppm = self.freq1d.pts2ppm(self.freq1d.dims[0]-1)
        xaxmax = self.FloatXaxisMax.GetValue()
        xaxmin = self.FloatXaxisMin.GetValue()
        curmax = self.FloatCursorMax.GetValue()
        curmin = self.FloatCursorMin.GetValue()
        if xaxmax > self.maxppm: xaxmax = self.maxppm
        if xaxmin < self.minppm: xaxmin = self.minppm
        if curmax > self.maxppm: curmax = self.maxppm
        if curmin < self.minppm: curmin = self.minppm
        self.FloatXaxisMax.SetRange(self.minppm, self.maxppm)
        self.FloatXaxisMin.SetRange(self.minppm, self.maxppm)
        self.FloatCursorMax.SetRange(self.minppm, self.maxppm)
        self.FloatCursorMin.SetRange(self.minppm, self.maxppm)
        self.FloatXaxisMax.SetValue(xaxmax)
        self.FloatXaxisMin.SetValue(xaxmin)
        self.FloatCursorMax.SetValue(curmax)
        self.FloatCursorMin.SetValue(curmin)


    def build_basis_functions(self, dim0=None):
        """
        Here we construct FID signals from the prior information from the
        experiment. These are rendered at the spectral resolution (number
        of points and sweep width) specified by the user. But, no line
        broadening is applied.

        By performing this recalculation for all metabolites under all
        conditions (ie. the FOR loops in the experiment) we are able to
        browse through the displayed spectra more quickly. However, for
        large numbers of simulations in an experiment, this recalculation
        can take a long time.

        We added the dim0 keyword here because when the user sets a new
        spectral resolution using the spectral dialog, we need to know
        what the new spectral dimension (aka. freq1d.dims[0]) value is,
        but since the basis functions have not yet been recalculated, we
        can not use the "dims" property to return a correct value from the
        data.shape attribute.  This keyword gets us around this conundrum.

        """
        if not self.experiment.metabolites:
            return

        if self.experiment.isotope == '1H':
            self.freq1d.resppm     = common_constants.DEFAULT_PROTON_CENTER_PPM
        else:
            self.freq1d.resppm     = common_constants.DEFAULT_XNUCLEI_CENTER_PPM

        # ensure that datasets and certain boundaries
        # reflect the experiment settings

        self.basis.frequency      = self.experiment.b0
        self.freq1d.frequency     = self.experiment.b0
        self.freq2d.frequency     = self.experiment.b0

        # can not use util_ppm methods because dims not properly set at this point
        if not dim0:
            dim0 = self.freq1d.dims[0]

        sw     = self.freq1d.sw
        freq   = self.freq1d.frequency
        resppm = self.freq1d.resppm
        self.maxppm = ( ((dim0/2) - (0)     ) * ((sw/dim0) / freq) ) + resppm
        self.minppm = ( ((dim0/2) - (dim0-1)) * ((sw/dim0) / freq) ) + resppm

        # if x-axis has changed, ensure bounds are appropriate
        axes = self.plot1d.axes[0]
        xmin = self.minppm
        xmax = self.maxppm
        x0, y0, x1, y1 = axes.dataLim.bounds
        # bjs - this was returning 'inf' values on init prior to any data being plotted
        #  not sure why it did not crash before. May need to put dummy data in on
        #  creation of the plot_panel
        axes.ignore_existing_data_limits = True
        if not all(np.isfinite([x0,y0,x1,y1])):
            axes.update_datalim([[xmin,0.0],[xmax,1.0]])
        else:
            axes.update_datalim([[xmin, y0],[xmax,y1+y0]])

        bx0, by0, bx1, by1 = axes.dataLim.bounds


        # We use a callback function to provide a progress indicator. However,
        # it only works under Linux/GTK.
        callback = self.build_basis_progress if ("__WXGTK__" in wx.PlatformInfo) else None

        self.top.statusbar.SetStatusText('Calculating spectra for Results', 0)
        self.basis.data = bbf_module.build_basis_functions(self.experiment, dim0, sw, resppm, callback)
        self.set_apodization()
        self.top.statusbar.SetStatusText(' ',0)


    def build_basis_progress(self):
        """
        Updates the status bar with build basis progress.
        Passed to and called by bbf_module.build_basis_functions().
        """
        _, progress = bbf_module.global_progress

        status = ("Basis functions %s complete...") % progress
        self.top.statusbar.SetStatusText(status)


    def display_results_text(self):

        lines = str(self.experiment)
        lines += "\n\nSimulation Results\n" + "-" * 75 + "\n\n"
        lines += "\n".join([simulation.summary() for simulation in self.experiment.simulations])
        wx_util.display_text_as_file(lines)


    def plot_canvas(self, reset_history=False):

        # check if any egregious issues with proceeding
        if not self.experiment.metabolites:
            return
        if self.basis.data is None:
            return

        # get indices of data to include
        imets = self.ListPlotMetabolites.GetSelections()
        # Under Windows & GTK, the list returned by GetSelections() is always
        # sorted smallest to largest. Under OS X, it's not sorted and might
        # be e.g. (2, 1, 0) or even (1, 0, 2). We sort it in order to ensure
        # consistent behavior.
        imets = sorted(imets)
        nmet  = len(imets)
        if nmet == 0:
            return

        self.plot_labels = [self.ListPlotMetabolites.GetString(i) for i in imets]
        self.plot_labels.reverse()
        index1 = self.SpinIndex1.GetValue()-1
        index2 = self.SpinIndex2.GetValue()-1
        index3 = self.SpinIndex3.GetValue()-1

        # check for display mode and fill data array
        mode = self.ChoiceDisplayMode.GetCurrentSelection()
        if mode == 0:
            b = self.basis.data[index3,index2,index1,imets,:]
        elif mode == 1:
            b = self.basis.data[index3,index2,:,imets[0],:]
        elif mode == 2:
            b = self.basis.data[index3,:,index1,imets[0],:]
        elif mode == 3:
            b = self.basis.data[:,index2,index1,imets[0],:]

        data = b.copy()
        data = data[::-1,:]

        if self.CheckSumPlots.IsChecked():
            data = data.sum(axis=0)
            data.shape = 1,data.shape[0]
            x00, y00, x11, y11 = self.plot1d_axes.dataLim.bounds

        npts = self.freq1d.dims[0]
        # transform into frequency domain and store
        if data.shape[0] == 1:
            data[:] = np.fft.fft(data[:] * self.apodization) / npts
        else:
            for i in range(data.shape[0]):
                data[i,:] = np.fft.fft(data[i,:] * self.apodization) / npts

        self.freq1d.data = data

        # select data type to plot - real, imaginary, magnitude
        if self._prefs.data_type_real:
            data = self.freq1d.data.real
            color = self._prefs.line_color_real
        elif self._prefs.data_type_imaginary:
            data = self.freq1d.data.imag
            color = self._prefs.line_color_imaginary
        elif self._prefs.data_type_magnitude:
            data = abs(self.freq1d.data)
            color = self._prefs.line_color_magnitude

        # check widget ranges if data has changed
        self.update_widget_ranges()

        # x-axis creation, zoom and cursor checks
        xx   = np.arange(self.freq1d.dims[0], dtype='float64')
        xmax = self.FloatXaxisMax.GetValue()
        xmin = self.FloatXaxisMin.GetValue()
        if self._prefs.xaxis_ppm:
            xx = self.freq1d.pts2ppm(xx)
        else:
            xx = self.freq1d.pts2hz(xx)
            xmax = self.freq1d.pts2hz(self.freq1d.ppm2pts(xmax))
            xmin = self.freq1d.pts2hz(self.freq1d.ppm2pts(xmin))
        xlim = (xmax, xmin)

        xmax = self.FloatCursorMax.GetValue()
        xmin = self.FloatCursorMin.GetValue()
        if self._prefs.xaxis_hertz:
            xmax = self.freq1d.pts2hz(self.freq1d.ppm2pts(xmax))
            xmin = self.freq1d.pts2hz(self.freq1d.ppm2pts(xmin))
        if self.plot1d.refs != None:
            self.plot1d.refs.set_span(xmin,xmax)

        # set up dynamic y-scaling based on min/max of data plotted
        ticklocs = []
        dmin = data.min()
        dmax = data.max()
        dr   = abs(dmax - dmin)
        dmin = dmin - dr*0.03   # space line collection vertically
        dmax = dmax + dr*0.03
        dr = dmax - dmin
        y0 = dmin
        y1 = dmin + data.shape[0] * dr

        # set up line collection requirements
        segs = []
        for i in range(data.shape[0]):
            segs.append(np.hstack((xx[:,np.newaxis], data[i,:,np.newaxis])))
            ticklocs.append(i*dr)
        offsets = np.zeros((data.shape[0],2), dtype=float)
        offsets[:,1] = ticklocs

        lines = mpl.collections.LineCollection(segs, offsets=offsets,
                                               linewidth=self._prefs.line_width)
        lines.set_color(color)

        self.plot1d_axes.cla()
        self.plot1d_axes.set_xlim(xlim)     # has to happen after cla()
        self.plot1d.refresh_cursors()
        self.plot1d_axes.add_collection(lines)
        self.plot1d_axes.set_yticks(ticklocs)

        if self._prefs.zero_line_show:
            for i in range(data.shape[0]):
                self.plot1d_axes.axhline(y = i*dr,
                                         color=self._prefs.zero_line_color,
                                         linestyle=self._prefs.zero_line_style,
                                         linewidth=self._prefs.line_width)

        # y-lim set needs to stay after zero_line or y-axis gets skewed
        self.plot1d_axes.set_ylim(y0, y1)
        x0, boundy0, x1, boundy1 = self.plot1d_axes.dataLim.bounds
        self.plot1d_axes.ignore_existing_data_limits = True
        self.plot1d_axes.update_datalim([[x0,y0],[x1+x0,y1]])

        self.plot1d_axes.axes.get_yaxis().set_visible(False)
        if not self._prefs.xaxis_show:
            self.plot1d_axes.axes.set_xticklabels([])

        self.plot1d.canvas.draw()


    def plot_contour(self, reset=False, recalc=True):

        # check if any egregious issues with proceeding
        if not self._prefs.contour_plot_show:
            reset = True
        if reset:
            self.contour_axes.cla()
            self.contour.canvas.draw()
            return
        if self.basis.data is None:
            return

        if recalc:
            # get indices of data to include
            imets = self.ListPlotMetabolites.GetSelections()
            # Under Windows & GTK, the list returned by GetSelections() is
            # always sorted smallest to largest. Under OS X, it's not sorted
            # and might be e.g. (2, 1, 0) or even (1, 0, 2). We sort it in
            # order to ensure consistent behavior.
            imets = sorted(imets)
            nmet  = len(imets)
            if nmet == 0:
                return

            index1 = self.SpinIndex1.GetValue()-1
            index2 = self.SpinIndex2.GetValue()-1
            index3 = self.SpinIndex3.GetValue()-1

            # check for display mode and fill data array
            mode = self.ChoiceContourMode.GetCurrentSelection()
            if mode == 0:
                view = self.basis.data[index3,:,:,imets[0],:]   # data dims are X,Y
                xlen = len(self.experiment.dims[0])
                ylen = len(self.experiment.dims[1])
            elif mode == 1:
                view = self.basis.data[:,index2,:,imets[0],:]   # data dims are X,Z
                xlen = len(self.experiment.dims[0])
                ylen = len(self.experiment.dims[2])
            elif mode == 2:
                view = self.basis.data[:,:,index1,imets[0],:]   # data dims are Y,Z
                xlen = len(self.experiment.dims[1])
                ylen = len(self.experiment.dims[2])
            if xlen == 1 and ylen ==1:
                return
            view = view.copy()
            view.shape = xlen*ylen, self.basis.dims[0]

            # transform into frequency domain and store
            for i in range(view.shape[0]):
                view[i,:] = np.fft.fft(view[i,:] * self.apodization) / self.basis.dims[0]
        else:
            view = self.freq2d.data
            # check for display mode and fill data array
            mode = self.ChoiceContourMode.GetCurrentSelection()
            if mode == 0:
                xlen = len(self.experiment.dims[0])
                ylen = len(self.experiment.dims[1])
            elif mode == 1:
                xlen = len(self.experiment.dims[1])
                ylen = len(self.experiment.dims[2])
            elif mode == 2:
                xlen = len(self.experiment.dims[1])
                ylen = len(self.experiment.dims[2])
            if xlen == 1 and ylen ==1:
                return

        if reset:
            view = view * 0

        # save data so updates don't need recalculation
        self.freq2d.data = view

        # get locations for and plot cursors
        x0 = self.freq1d.ppm2pts(self.FloatCursorMax.GetValue())
        x1 = self.freq1d.ppm2pts(self.FloatCursorMin.GetValue())
        x0 = int(np.clip(x0, 0, self.freq1d.dims[0]-1))
        x1 = int(np.clip(x1, 1, self.freq1d.dims[0]-1))
        view = view[:,x0:x1]
        view = view.sum(axis=1)
        view.shape = xlen,ylen

        # select data type to plot - real, imaginary, magnitude
        if self._prefs.data_type_real:
            view = view.real
        elif self._prefs.data_type_imaginary:
            view = view.imag
        elif self._prefs.data_type_magnitude:
            view = abs(view)

        self.contour_axes.clear()

        if not self._prefs.contour_axes_show:
            self.contour_axes.axes.set_xticklabels([])
            self.contour_axes.axes.set_yticklabels([])

        # if data is all the same, the contour plot fails, I think it
        # is because it can not decide where to put multiple level lines,
        # in this case, set an alert message and return
        self.statusbar.SetStatusText((''), 2)
        if abs((view.max()-view.min())/view.mean())<0.001:
             txt = "Contour can not plot constant areas"
             self.statusbar.SetStatusText((txt), 2)
             self.contour.canvas.draw()
             return

        # in some cases it is desirable to expand a 1D array out to pretend
        # that it is a 2d array, this code checks which dimension is len = 1
        # and expands the array into a square matrix for display
        if xlen == 1:
            y = np.ndarray(ylen*ylen,dtype='float')
            for i in np.arange(ylen):
                for j in np.arange(ylen):
                    y[j+i*ylen] = view[0,j]
            view = y
            view.shape = ylen,ylen
            xlen = ylen

        if ylen == 1:
            y = np.ndarray(xlen*xlen,dtype='float')
            for i in np.arange(xlen):
                for j in np.arange(xlen):
                    y[j+i*xlen] = view[j,0]
            view = y
            view.shape = xlen,xlen
            ylen = xlen

        # FIXME bjs - magic number multiplier to not get underflow exception
        # the 10000 multiplier below was necessary to avoid an odd
        # FloatPointError exception that was being throw for various
        # levels of contour.  With bigger values in view, it seems
        # to not occur. I should normalize this scaling using a
        # constant so as to be able to control it more easily, but
        # not today.

        levels = self.SpinContourLevels.GetValue()
        if self.CheckGrayscale.IsChecked():
            self.contour_axes.imshow(view, cmap=cm.gray, extent=(1,xlen,1,ylen), origin='lower')
        self.contour_axes.contour(view*10000, levels, extent=(1,xlen,1,ylen))
        self.contour.canvas.draw()

        self.contour.canvas.draw()  # np.max(view)


    def plot_integral(self, reset=False):

        # check if any egregious issues with proceeding
        if not self._prefs.integral_plot_show:
            return
        if self.freq1d.data is None:
            return

        # check that there are metabolites selected
        imets = self.ListPlotMetabolites.GetSelections()
        if len(imets) == 0:
            return

        # select data type to plot - real, imaginary, magnitude
        if self._prefs.data_type_real:
            view = self.freq1d.data.real
        elif self._prefs.data_type_imaginary:
            view = self.freq1d.data.imag
        elif self._prefs.data_type_magnitude:
            view = abs(self.freq1d.data)
        if len(view.shape)<2:
            return

        # check for display mode and fill x-axis
        mode = self.ChoiceDisplayMode.GetCurrentSelection()
        if mode == 0:
            return
        else:
            xx = np.array(self.experiment.dims[mode - 1])
        if len(xx)<2:
            return

        if self._prefs.xaxis_hertz:
            xx = xx * self.freq1d.frequency   # convert generally to Hz

        x0 = self.freq1d.ppm2pts(self.FloatCursorMax.GetValue())
        x1 = self.freq1d.ppm2pts(self.FloatCursorMin.GetValue())
        x0 = int(np.clip(x0, 0, self.freq1d.dims[0]-1))
        x1 = int(np.clip(x1, 1, self.freq1d.dims[0]-1))
        view = view[:,x0:x1]
        y = view.sum(axis=1)
        y = y[::-1]

        if reset:
            y = y * 0
            ymin = -1
            ymax = 1
        else:
            ymin = y.min() - (abs(y.min())*0.03)
            ymax = y.max() + (abs(y.max())*0.03)

        self.integral_data = y              # used during motion events
        self.integral_axes.clear()
        self.integral.refresh_cursors()
        self.integral_axes.plot(xx, y, linewidth=self._prefs.line_width)
        self.integral_axes.set_xlim(xx.min(), xx.max())
        self.integral_axes.set_ylim(ymin, ymax)
        if not self._prefs.integral_xaxis_show:
            self.integral_axes.axes.set_xticklabels([])
        if not self._prefs.integral_yaxis_show:
            self.integral_axes.axes.set_yticklabels([])

        self.integral.canvas.draw()


    def set_apodization(self, dim0=None):
        """
        We added the dim0 keyword here because when the user sets a new
        spectral resolution using the spectral dialog, we need to know
        what the new spectral dimension (aka. freq1d.dims[0]) value is,
        but since the basis functions have not yet been recalculated, we
        can not use the "dims" property to return a correct value from the
        data.shape attribute.  This keyword gets us around this conundrum.

        """
        if not dim0:
            dim0 = self.freq1d.dims[0]

        xx = np.arange(dim0, dtype='float64') / self.freq1d.sw
        type_ = 'Gaussian' if self._prefs.line_shape_gaussian else 'Lorentzian'
        self.apodization = util_generic_spectral.apodize(xx, self.linewidth, type_)




