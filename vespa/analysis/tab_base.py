# Python modules

# 3rd party modules
import wx

# Our modules
import vespa.analysis.util_menu as util_menu


#------------------------------------------------------------------------------
#
#  Tab base class
#
#------------------------------------------------------------------------------

class Tab(object):
    """
    This is an abstract base class for Analysis tabs. It provides services
    common to all the tabs (like the dataset property) and describes a common
    subset of features.

    """

    def __init__(self, tab_dataset, top, PrefsClass):

        self._top = wx.GetApp().GetTopWindow()
        self._tab_dataset = tab_dataset
        self._prefs = (PrefsClass() if PrefsClass else None)

        # Most tabs have a plot panel which needs to init the scale to a sane
        # value after data has been passed to the view. _scale_intialized
        # tracks whether or not this has happened. It's intialized to False,
        # set to True once the scale is intialized and never changes thereafter.
        self._scale_intialized = False

        # Every tab has a view attr. The ones that don't have a view (like the
        # raw tab) can just leave this set to None.
        self.view = None

        self.Bind(wx.EVT_WINDOW_DESTROY, self.on_destroy, self)


    @property
    def is_raw(self):
        return hasattr(self, "IS_RAW")

    @property
    def is_prep_fidsum(self):
        return hasattr(self, "IS_PREP_FIDSUM")

    @property
    def is_prep_timeseries(self):
        return hasattr(self, "IS_PREP_TIMESERIES")

    @property
    def is_voigt(self):
        return hasattr(self, "IS_VOIGT")

    @property
    def is_giso(self):
        return hasattr(self, "IS_GISO")

    @property
    def is_spectral(self):
        return hasattr(self, "IS_SPECTRAL")

    @property
    def is_watref(self):
        return hasattr(self, "IS_WATREF")


    @property
    def dataset(self):
        return self._top.datasets.get(self._tab_dataset.indexAB[0], None)

    #=======================================================
    #
    #           wx Event Handlers
    #
    #=======================================================

    def on_destroy(self, event):
        if self._prefs:
            self._prefs.save()



    #=======================================================
    #
    #           Public Methods
    #
    #=======================================================

    def on_activation(self):
        """
        Sets the dataset tab's scale floatspin to the appropriate value for
        this tab. Also forces the view menu to match this tab's prefs.

        This should be called from the dataset notebook whenever a new tab
        becomes the current tab.
        """
        if self.view:
            self._tab_dataset.FloatScale.SetValue(self.view.vertical_scale)

        if self._prefs:
            # Force the View menu to match the current plot options.
            util_menu.bar.set_menu_from_state(self._prefs.menu_state)


    def on_voxel_change(self, voxel):
        pass


    def process(self, entry='all'):
        """
        The process(), plot() and process_and_plot() methods are standard in
        all processing tabs. They are called to update the data in the chain
        object, the plot_panel in the View side of the tab or both.

        The steps that are actually run for each call to process() can be
        controlled using the entry keyword.
        """
        pass


    def process_and_plot(self, entry='all'):
        """
        The process(), plot() and process_and_plot() methods are standard in
        all processing tabs. They are called to update the data in the chain
        object, the plot_panel in the View side of the tab or both.

        The steps that are actually run for each call to process() can be
        controlled using the entry keyword.
        """
        pass


    def plot(self):
        """
        The process(), plot() and process_and_plot() methods are standard in
        all processing tabs. They are called to update the data in the chain
        object, the plot_panel in the View side of the tab or both.
        """
        pass


    def build_area_text(self, area, rms, plot_label=None):
        txt = ' Area = %1.5g  RMS = %1.5g' % (area, rms)
        if plot_label:
            # used to indicate if area value from Plot A, B, C, D
            txt = plot_label+': '+txt
        return txt

