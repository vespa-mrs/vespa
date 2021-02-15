# Python modules

import abc

# 3rd party modules

# Our modules
import vespa.datasim.util_menu as util_menu
import vespa.datasim.util_datasim_config as util_datasim_config
import vespa.common.prefs as prefs

"""See common/prefs.py for info on the classes below."""

class DatasimPrefs(prefs.Prefs, metaclass=abc.ABCMeta):
    def __init__(self, id_class):
        prefs.Prefs.__init__(self, util_menu.bar, id_class)

    @property
    def _ConfigClass(self):
        """Returns the appropriate ConfigObj class for this app."""
        return util_datasim_config.Config


class PrefsMain(DatasimPrefs):
    def __init__(self):
        DatasimPrefs.__init__(self, util_menu.ViewIds)


    @property
    def _ini_section_name(self):
        return "main_prefs"


    def deflate(self):
        # Call my base class deflate
        d = DatasimPrefs.deflate(self)

        # Add my custom stuff
        for name in ("sash_position",
                     "line_color_baseline",
                     "line_color_imaginary", 
                     "line_color_magnitude",
                     "line_color_metabolite", 
                     "line_color_real", 
                     "zero_line_color", 
                     "zero_line_style", 
                     "line_width", 
                    ):
            d[name] = getattr(self, name)

        return d


    def inflate(self, source):
        # Call my base class inflate
        DatasimPrefs.inflate(self, source)

        # Add my custom stuff
        for name in ("sash_position", ):
            setattr(self, name, int(source[name]))

        for name in ("line_color_baseline",
                     "line_color_imaginary", 
                     "line_color_magnitude",
                     "line_color_metabolite", 
                     "line_color_real", 
                     "zero_line_color", 
                     "zero_line_style",
                    ):
            setattr(self, name, source[name])

        for name in ("line_width", ):
            setattr(self, name, float(source[name]))


