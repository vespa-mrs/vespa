# Python modules

import abc

# 3rd party modules

# Our modules
import vespa.simulation.util_menu as util_menu
import vespa.simulation.util_simulation_config as util_simulation_config
import vespa.common.prefs as prefs

"""See common/prefs.py for info on the classes below."""

class SimulationPrefs(prefs.Prefs, metaclass=abc.ABCMeta):
    def __init__(self, id_class):
        prefs.Prefs.__init__(self, util_menu.bar, id_class)

    @property
    def _ConfigClass(self):
        """Returns the appropriate ConfigObj class for this app."""
        return util_simulation_config.Config


class PrefsMain(SimulationPrefs):
    def __init__(self):
        SimulationPrefs.__init__(self, util_menu.ViewIds)


    @property
    def _ini_section_name(self):
        return "main_prefs"


    def deflate(self):
        # Call my base class deflate
        d = SimulationPrefs.deflate(self)

        # Add my custom stuff
        for name in ("contour_levels", "sash_position",
                     "line_color_real", "line_color_imaginary",
                     "line_color_magnitude", "zero_line_color", 
                     "zero_line_style", "line_width", 
                     "contour_grayscale",                     
                    ):
            d[name] = getattr(self, name)

        return d


    def inflate(self, source):
        # Call my base class inflate
        SimulationPrefs.inflate(self, source)

        # Add my custom stuff
        for name in ("contour_levels",
                     "sash_position"):
            setattr(self, name, int(source[name]))

        for name in ("line_color_real", "line_color_imaginary", 
                     "line_color_magnitude", "zero_line_color", 
                     "zero_line_style",
                    ):
            setattr(self, name, source[name])

        for name in ("line_width", ):
            setattr(self, name, float(source[name]))

        for name in ("contour_grayscale", ):
            # FIXME PS - assuming source is a ConfigObj Section
            setattr(self, name, source.as_bool(name))



