# Python modules


# 3rd party modules

# Our modules
import vespa.common.util.config as util_config

def get_window_coordinates(window_name):
    """A shortcut for the Config object method of the same name."""
    config = Config()
    return config.get_window_coordinates(window_name)


class Config(util_config.BaseConfig):

    # This base class initializes itself with the Default values from
    # 'default_ini_file_contents.py', then checks these values against the
    # actual values in the (maybe) existent INI file. If some values are
    # not yet in the INI file, the default values are used automagically.

    def __init__(self):
        util_config.BaseConfig.__init__(self, "simulation.ini")  
        

