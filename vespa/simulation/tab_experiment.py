# Python modules

import math
import os
import tempfile


# 3rd party modules
import wx
#import wx.aui as aui
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


class TabExperiment(vespa_notebooks.VespaAuiNotebook):

    def __init__(self, parent, top, db, experiment, is_new):
        style = wx.BORDER_NONE
        
        agw_style = aui.AUI_NB_TAB_SPLIT            |   \
                    aui.AUI_NB_TAB_MOVE             |   \
                    aui.AUI_NB_TAB_EXTERNAL_MOVE    |   \
                    aui.AUI_NB_BOTTOM 

        vespa_notebooks.VespaAuiNotebook.__init__(self, parent, style, agw_style)

        # global attributes

        self.top        = top
        self.parent     = parent
        self.db         = db
        self.experiment = experiment
        self.statusbar  = top.statusbar
        self.is_new     = is_new

        # plot attributes
        self._prefs = prefs.PrefsMain()

        self.simulate  = None
        self.visualize = None

        self.initialize_controls()

        self.Bind(wx.EVT_WINDOW_DESTROY, self.on_destroy, self)


    #=======================================================
    #
    #           Event handlers start here
    #
    #=======================================================


    def on_activation(self):
        # This is a faux event handler. wx doesn't call it directly. It's
        # a notification from my parent (the experiment notebook) to let
        # me know that this tab has become the current one.

        # Force the View menu to match the current plot options.
        util_menu.bar.set_menu_from_state(self._prefs.menu_state)


    def on_destroy(self, event):
        self._prefs.save()


    def on_menu_view_option(self, event):

        reset_history = False
        plot_canvas   = True
        plot_integral = False
        plot_contour  = False

        event_id = event.GetId()

        if self._prefs.handle_event(event_id):
            if event_id in (util_menu.ViewIds.INTEGRAL_XAXIS_SHOW,
                            util_menu.ViewIds.INTEGRAL_YAXIS_SHOW,
                            util_menu.ViewIds.CONTOUR_AXES_SHOW,
                            util_menu.ViewIds.INTEGRAL_PLOT_SHOW,
                            util_menu.ViewIds.CONTOUR_PLOT_SHOW,
                           ):
                # For these events, we don't need to bother replotting the
                # main plot.
                plot_canvas = False

            # if event_id in (util_menu.ViewIds.XAXIS_PPM,
            #                 util_menu.ViewIds.XAXIS_HERTZ,
            #                ):
            #     reset_history = True
            #     self.visualize.plot1d.reversex = (event_id==util_menu.ViewIds.XAXIS_PPM)

            if event_id in (util_menu.ViewIds.DATA_TYPE_REAL,
                            util_menu.ViewIds.DATA_TYPE_IMAGINARY,
                            util_menu.ViewIds.DATA_TYPE_MAGNITUDE,
                            util_menu.ViewIds.LINE_SHAPE_GAUSSIAN,
                            util_menu.ViewIds.LINE_SHAPE_LORENTZIAN,
                            util_menu.ViewIds.INTEGRAL_PLOT_SHOW,
                            util_menu.ViewIds.INTEGRAL_XAXIS_SHOW,
                            util_menu.ViewIds.INTEGRAL_YAXIS_SHOW,
                           ):
                plot_integral = True

            if event_id in (util_menu.ViewIds.CONTOUR_PLOT_SHOW,
                            util_menu.ViewIds.CONTOUR_AXES_SHOW,
                           ):
                plot_contour = True

            if event_id in (util_menu.ViewIds.LINE_SHAPE_GAUSSIAN,
                            util_menu.ViewIds.LINE_SHAPE_LORENTZIAN,
                           ):
                self.visualize.set_apodization()

            if plot_canvas:
                self.visualize.plot_canvas(reset_history=reset_history)
            if plot_integral:
                self.visualize.plot_integral()
            if plot_contour:
                self.visualize.plot_contour()


    def on_menu_view_output(self, event):
        event_id = event.GetId()

        if event_id == util_menu.ViewIds.OUTPUT_TEXT_RESULTS:
            self.visualize.display_results_text()
        else:
            # "fands" is short for "formats and sources"
            fands = {
                util_menu.ViewIds.OUTPUT_1D_STACK_TO_PNG : ("PNG", self.visualize.plot1d),
                util_menu.ViewIds.OUTPUT_INTEGRAL_TO_PNG : ("PNG", self.visualize.integral),
                util_menu.ViewIds.OUTPUT_CONTOUR_TO_PNG  : ("PNG", self.visualize.contour),
                util_menu.ViewIds.OUTPUT_1D_STACK_TO_SVG : ("SVG", self.visualize.plot1d),
                util_menu.ViewIds.OUTPUT_INTEGRAL_TO_SVG : ("SVG", self.visualize.integral),
                util_menu.ViewIds.OUTPUT_CONTOUR_TO_SVG  : ("SVG", self.visualize.contour),
                util_menu.ViewIds.OUTPUT_1D_STACK_TO_EPS : ("EPS", self.visualize.plot1d),
                util_menu.ViewIds.OUTPUT_INTEGRAL_TO_EPS : ("EPS", self.visualize.integral),
                util_menu.ViewIds.OUTPUT_CONTOUR_TO_EPS  : ("EPS", self.visualize.contour),
                util_menu.ViewIds.OUTPUT_1D_STACK_TO_PDF : ("PDF", self.visualize.plot1d),
                util_menu.ViewIds.OUTPUT_INTEGRAL_TO_PDF : ("PDF", self.visualize.integral),
                util_menu.ViewIds.OUTPUT_CONTOUR_TO_PDF  : ("PDF", self.visualize.contour),
                      }

            format, source = fands[event_id]
            lformat = format.lower()
            filter_ = "%s files (*.%s)|*.%s" % (format, lformat, lformat)

            filename = common_dialogs.save_as("", filter_)

            if filename:
                msg = ''
                try:
                    source.figure.savefig(filename,
                                          dpi=300,
                                          facecolor='w',
                                          edgecolor='w',
                                          orientation='portrait',
                                          #papertype='letter',
                                          format=None,
                                          transparent=False)
                except IOError:
                    msg = """I can't write the file "%s".""" % filename

                if msg:
                    common_dialogs.message(msg, style=common_dialogs.E_OK)


    #=======================================================
    #
    #           Internal methods start here
    #
    #=======================================================

    def initialize_controls(self):
        """
        This is called once during init to build the visualize and
        simulate tabs. This code could just as well be inline in __init__().

        """

        self.visualize = tab_visualize.TabVisualize(self,
                                                    self.top,
                                                    self.db,
                                                    self.experiment,
                                                    self._prefs)

        self.simulate   = tab_simulate.TabSimulate(self,
                                                   self.top,
                                                   self.db,
                                                   self.experiment,
                                                   self.is_new)

        self.AddPage(self.simulate, "  Simulate  ", select=False)
        self.AddPage(self.visualize,"  Visualize  ",select=(not self.is_new))


    #==================================================================
    #
    #       Shared/public methods start here
    #
    # In most (all?) cases these are only shared between this tab and
    # its two child tabs (Simulate & Visualize)
    #==================================================================

    def close(self):
        return self.simulate.close()


    def save_experiment(self):
        """
        Saves the experiment, if possible, and returns True. If the
        experiment can't be saved (because it contains invalid data, e.g.),
        returns False.

        This is quite similar to an event handler since it is called
        from main where it's invoked by a menu item.

        """
        return self.simulate.save_experiment()


    def refresh_visualize(self):
        """ Typically called from Simulate tab when an experiment completes. """
        self.visualize.build_basis_functions()
        self.visualize.update_plot_controls()


    def update_plot_controls(self):
        """
        This update routine is located on the Experiment_Tab level because
        changes to the Simulate Tab can affect the Visualize Tab

        update the Visualize tab controls

        """
        self.visualize.update_plot_controls()


