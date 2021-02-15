# Python modules
import abc

# 3rd party modules

# Our modules
import vespa.analysis.util_menu as util_menu
import vespa.analysis.util_analysis_config as util_analysis_config
import vespa.common.prefs as prefs
import vespa.common.util.xml_ as util_xml


"""See common/prefs.py for info on the classes below."""

class AnalysisPrefs(prefs.Prefs, metaclass=abc.ABCMeta):
    def __init__(self, id_class):
        prefs.Prefs.__init__(self, util_menu.bar, id_class)

    @property
    def _ConfigClass(self):
        """Returns the appropriate ConfigObj class for this app, specifically
        util_analysis_config.Config.
        """
        return util_analysis_config.Config


class PrefsPrepFidsum(AnalysisPrefs):
    def __init__(self):
        AnalysisPrefs.__init__(self, util_menu.ViewIdsPrepFidsum)


    @property
    def _ini_section_name(self):
        return "prep_fidsum_prefs"


    def deflate(self):
        # Call my base class deflate
        d = AnalysisPrefs.deflate(self)

        # Add my custom stuff
        for name in ("sash_position",
                     "line_color_imaginary", 
                     "line_color_magnitude", 
                     "line_color_real", 
                     "zero_line_color", 
                     "zero_line_style",
                     "zero_line_plot_color",
                     "zero_line_plot_style",
                     "line_width",
                    ):
            d[name] = getattr(self, name)

        return d


    def inflate(self, source):
        # Call my base class inflate
        AnalysisPrefs.inflate(self, source)

        # Add my custom stuff
        for name in ("sash_position", ):
            setattr(self, name, int(source[name]))

        for name in ("line_color_imaginary", 
                     "line_color_magnitude", 
                     "line_color_real",
                     "zero_line_color", 
                     "zero_line_style",  
                     "zero_line_plot_color",
                     "zero_line_plot_style",
                    ):
            setattr(self, name, source[name])

        for name in ("line_width", ):
            setattr(self, name, float(source[name]))


class PrefsSpectral(AnalysisPrefs):
    def __init__(self):
        AnalysisPrefs.__init__(self, util_menu.ViewIdsSpectral)


    @property
    def _ini_section_name(self):
        return "spectral_prefs"


    def deflate(self):
        # Call my base class deflate
        d = AnalysisPrefs.deflate(self)

        # Add my custom stuff
        for name in ("sash_position_main", 
                     "sash_position_svd", 
                     "line_color_imaginary", 
                     "line_color_magnitude", 
                     "line_color_real", 
                     "line_color_svd", 
                     "zero_line_color", 
                     "zero_line_style",
                     "line_width", 
                    ):
            d[name] = getattr(self, name)

        return d


    def inflate(self, source):
        # Call my base class inflate
        AnalysisPrefs.inflate(self, source)

        # Add my custom stuff
        for name in ("sash_position_main", 
                     "sash_position_svd", ):
            setattr(self, name, int(source[name]))

        for name in ("line_color_imaginary", 
                     "line_color_magnitude", 
                     "line_color_real",
                     "line_color_svd", 
                     "zero_line_color", 
                     "zero_line_style", 
                    ):
            setattr(self, name, source[name])

        for name in ("line_width", ):
            setattr(self, name, float(source[name]))


class PrefsVoigt(AnalysisPrefs):
    def __init__(self):
        AnalysisPrefs.__init__(self, util_menu.ViewIdsVoigt)

        # FIXME PS - magic number
        self.plotx = [PrefsVoigtPlotX(i) for i in range(4)]


    @property
    def _ini_section_name(self):
        return "voigt_prefs"


    @property
    def menu_state(self):
        state = { }
        for menu_item_id in self._id_name_map:
            name = self._id_name_map[menu_item_id]
            state[menu_item_id] = getattr(self, name)

        for plotx in self.plotx:
            state.update(plotx.menu_state)

        return state


    @property
    def n_plots(self):
        """Returns the # of plots that should be visible. Read-only."""
        if self.n_plots_1:
            return 1
        elif self.n_plots_2:
            return 2
        elif self.n_plots_3:
            return 3
        elif self.n_plots_4:
            return 4


    def handle_event(self, menu_item_id):
        if menu_item_id in self._id_name_map:
            # This is an ordinary boolean. My base class can handle it.
            return AnalysisPrefs.handle_event(self, menu_item_id)
        else:
            for i, plotx in enumerate(self.plotx):
                if plotx.handle_event(menu_item_id):
                    return True

            return False


    def deflate(self):
        # Call my base class deflate
        d = AnalysisPrefs.deflate(self)

        # Get rid of the plotx items
        del d["plotx"]

        # Add my custom stuff
        for name in ("sash_position", 
                     "line_color_raw", 
                     "line_color_fit", 
                     "line_color_base",
                     "line_color_init", 
                     "line_color_weight", 
                     "line_color_imaginary", 
                     "line_color_magnitude", 
                     "line_color_real",
                     "zero_line_color", 
                     "zero_line_style",
                     "line_width", 
                     "csv_qa_metab_labels",
                    ):
            d[name] = getattr(self, name)

        return d


    def inflate(self, source):
        # Due to a bug, for a while we wrote "n_plots" to the INI file, so it
        # might be present in source. If it is, we need to get rid of it, 
        # otherwise we'll try to overwrite the read-only "n_plots" property
        # during inflate.
        if "n_plots" in source:
            del source["n_plots"]

        # Call my base class inflate
        AnalysisPrefs.inflate(self, source)

        # Add my custom stuff
        self.sash_position = int(source.get("sash_position", 0))

        for name in ("line_color_raw", 
                     "line_color_fit", 
                     "line_color_base",
                     "line_color_init", 
                     "line_color_weight", 
                     "line_color_imaginary", 
                     "line_color_magnitude", 
                     "line_color_real",
                     "zero_line_color", 
                     "zero_line_style", 
                     ):
            setattr(self, name, source[name])

        for name in ("line_width", 
                     ):
            setattr(self, name, float(source[name]))

        self.csv_qa_metab_labels = util_xml.BOOLEANS[source.get("csv_qa_metab_labels", False)]


    def save(self):
        AnalysisPrefs.save(self)

        for plotx in self.plotx:
            plotx.save()


class PrefsVoigtPlotX(AnalysisPrefs):
    def __init__(self, index):
        self._index = index
        AnalysisPrefs.__init__(self, util_menu.PlotXIds)

        # At this point (after my base class init), the keys in _id_name_map 
        # are 4-tuples. I need to reduce them to single ints to make this 
        # _id_name_map look mimic the behavior of _id_name_map in all of the
        # other Prefs classes.
        keys = list(self._id_name_map.keys())

        self._id_name_map = \
                    dict( (key[self.index], value) for key, value 
                                                   in self._id_name_map.items()) 


    @property
    def index(self):
        return self._index


    @property
    def _ini_section_name(self):
        # The ini section name is 'voigt_plot_x' where x = 'a', 'b', 'c', or 'd'.
        return "voigt_plot_%s" % chr(ord('a') + self.index)


    def handle_event(self, menu_item_id):
        # This version of handle_event differs from the base class version 
        # because I can't guarantee that the id passed to me belongs to my
        # menu at all.
        if menu_item_id in self._id_name_map:
            return AnalysisPrefs.handle_event(self, menu_item_id)
        else:
            # Not my menu
            return False


class PrefsGiso(AnalysisPrefs):
    def __init__(self):
        AnalysisPrefs.__init__(self, util_menu.ViewIdsGiso)

        # FIXME PS - magic number
        self.plotx = [PrefsGisoPlotX(i) for i in range(4)]


    @property
    def _ini_section_name(self):
        return "giso_prefs"


    @property
    def menu_state(self):
        state = { }
        for menu_item_id in self._id_name_map:
            name = self._id_name_map[menu_item_id]
            state[menu_item_id] = getattr(self, name)

        for plotx in self.plotx:
            state.update(plotx.menu_state)

        return state


    @property
    def n_plots(self):
        """Returns the # of plots that should be visible. Read-only."""
        if self.n_plots_1:
            return 1
        elif self.n_plots_2:
            return 2
        elif self.n_plots_3:
            return 3
        elif self.n_plots_4:
            return 4


    def handle_event(self, menu_item_id):
        if menu_item_id in self._id_name_map:
            # This is an ordinary boolean. My base class can handle it.
            return AnalysisPrefs.handle_event(self, menu_item_id)
        else:
            for i, plotx in enumerate(self.plotx):
                if plotx.handle_event(menu_item_id):
                    return True

            return False


    def deflate(self):
        # Call my base class deflate
        d = AnalysisPrefs.deflate(self)

        # Get rid of the plotx items
        del d["plotx"]

        # Add my custom stuff
        for name in ("sash_position", 
                     "line_color_raw", 
                     "line_color_fit", 
                     "line_color_base",
                     "line_color_init", 
                     "line_color_weight", 
                     "line_color_imaginary", 
                     "line_color_magnitude", 
                     "line_color_real",
                     "zero_line_color", 
                     "zero_line_style",
                     "line_width", 
                    ):
            d[name] = getattr(self, name)

        return d


    def inflate(self, source):
        # Due to a bug, for a while we wrote "n_plots" to the INI file, so it
        # might be present in source. If it is, we need to get rid of it, 
        # otherwise we'll try to overwrite the read-only "n_plots" property
        # during inflate.
        if "n_plots" in source:
            del source["n_plots"]

        # Call my base class inflate
        AnalysisPrefs.inflate(self, source)

        # Add my custom stuff
        self.sash_position = int(source.get("sash_position", 0))

        for name in ("line_color_raw", 
                     "line_color_fit", 
                     "line_color_base",
                     "line_color_init", 
                     "line_color_weight", 
                     "line_color_imaginary", 
                     "line_color_magnitude", 
                     "line_color_real",
                     "zero_line_color", 
                     "zero_line_style", ):
            setattr(self, name, source[name])

        for name in ("line_width", ):
            setattr(self, name, float(source[name]))


    def save(self):
        AnalysisPrefs.save(self)

        for plotx in self.plotx:
            plotx.save()



class PrefsGisoPlotX(AnalysisPrefs):
    def __init__(self, index):
        self._index = index
        AnalysisPrefs.__init__(self, util_menu.PlotXIds)

        # At this point (after my base class init), the keys in _id_name_map 
        # are 4-tuples. I need to reduce them to single ints to make this 
        # _id_name_map look mimic the behavior of _id_name_map in all of the
        # other Prefs classes.
        keys = list(self._id_name_map.keys())

        self._id_name_map = \
                    dict( (key[self.index], value) for key, value 
                                                   in self._id_name_map.items()) 


    @property
    def index(self):
        return self._index


    @property
    def _ini_section_name(self):
        # The ini section name is 'voigt_plot_x' where x = 'a', 'b', 'c', or 'd'.
        return "giso_plot_%s" % chr(ord('a') + self.index)


    def handle_event(self, menu_item_id):
        # This version of handle_event differs from the base class version 
        # because I can't guarantee that the id passed to me belongs to my
        # menu at all.
        if menu_item_id in self._id_name_map:
            return AnalysisPrefs.handle_event(self, menu_item_id)
        else:
            # Not my menu
            return False



class PrefsWatref(AnalysisPrefs):
    def __init__(self):
        AnalysisPrefs.__init__(self, util_menu.ViewIdsWatref)


    @property
    def _ini_section_name(self):
        return "watref_prefs"


    def deflate(self):
        # Call my base class deflate
        d = AnalysisPrefs.deflate(self)

        # Add my custom stuff
        for name in ("sash_position_main",
                     "csv_qa_metab_labels", 
                    ):
            d[name] = getattr(self, name)

        return d


    def inflate(self, source):
        # Call my base class inflate
        AnalysisPrefs.inflate(self, source)

        # Add my custom stuff
        for name in ("sash_position_main",):
            setattr(self, name, int(source[name]))
            
        self.csv_qa_metab_labels = util_xml.BOOLEANS[source.get("csv_qa_metab_labels", False)]
            

