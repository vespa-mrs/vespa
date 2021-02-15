# Python modules

import abc

# 3rd party modules

# Our modules
import vespa.pulse.util_menu as util_menu
import vespa.pulse.util_pulse_config as util_pulse_config
import vespa.common.prefs as prefs

"""See common/prefs.py for info on the classes below."""

class PulsePrefs(prefs.Prefs, metaclass=abc.ABCMeta):
    def __init__(self, id_class):
        prefs.Prefs.__init__(self, util_menu.bar, id_class)

    @property
    def _ConfigClass(self):
        """Returns the appropriate ConfigObj class for this app."""
        return util_pulse_config.Config


class PrefsMain(PulsePrefs):
    def __init__(self):
        PulsePrefs.__init__(self, util_menu.ViewIds)


    @property
    def _ini_section_name(self):
        return "main_prefs"


    def deflate(self):
        # Call my base class deflate
        d = PulsePrefs.deflate(self)

        # Add my custom stuff
        for name in ("sash_position",
                     "data_type_real",
                     "data_type_real_imaginary",
                     "line_color_real", 
                     "line_color_imaginary", 
                     "line_color_magnitude", 
                     "line_color_phase_degrees",
                     "zero_line_color", 
                     "zero_line_style", 
                     "line_width", 
                    ):
            d[name] = getattr(self, name)

        return d


    def inflate(self, source):
        # Call my base class inflate
        PulsePrefs.inflate(self, source)

        # Add my custom stuff
        for name in ("sash_position", ):
            setattr(self, name, int(source[name]))

        for name in ("line_color_real", 
                     "line_color_imaginary", 
                     "line_color_magnitude", 
                     "line_color_phase_degrees",
                     "zero_line_color", 
                     "zero_line_style",
                    ):
            setattr(self, name, source[name])

        for name in ("data_type_real",
                     "data_type_real_imaginary",
                     "xaxis_show",
                     "zero_line_show",
                    ):
            setattr(self, name, bool(source[name]))

        for name in ("line_width", ):
            setattr(self, name, float(source[name]))


