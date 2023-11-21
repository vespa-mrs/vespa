# Python modules

import os

# 3rd party modules

# need for inline processing - no wx
try:
    import wx
except:
    wx = None


# Our modules
import vespa.common.configobj as configobj
import vespa.common.util.misc as util_misc
import vespa.common.default_ini_file_content as default_ini_file_content

def get_last_export_path():
    """A shortcut for the Config object method of the same name."""
    config = VespaConfig()
    return config.get_last_export_path()

def set_last_export_path(path):
    """A shortcut for the Config object method of the same name."""
    config = VespaConfig()
    config.set_last_export_path(path)
    config.write()


class BaseConfig(configobj.ConfigObj):
    """This abstract base class is a special purpose version of a 
    configobj.ConfigObj which is a 3rd party INI file parser that we've 
    adopted. It's a lot nicer than the Python standard library's ConfigParser.
    
    There's documentation for ConfigObj here:
    http://www.voidspace.org.uk/python/configobj.html
    
    This class focuses on saving config values (like window size & position). 
    Methods like get_window_coordinates() return wxWidgets-friendly tuples.
    
    Note that this class doesn't automatically save changes to disk. It's
    up to the caller to call the .write() method. Since the entire INI file
    is rewritten when .write() is called, it's not a good idea to keep one
    of these objects around for a long time before calling .write() because
    you might overwrite changes written by another part of the app.

    As mentioned above, this is meant to be an abstract base class. In truth,
    one can instantiate this class directly, but it's not that useful to
    do so. 
    
    There's subclasses of this class for reading specific Vespa INI files,
    such as the VespaConfig class. Each subclass first reads the default
    INI file content and then updates (overrides) that with whatever settings 
    the user provides.
    """
    def __init__(self, filename):
        """ Constructor. The filename should be something like "analysis.ini".
        It will be prefixed with an app-specific version of the user's data 
        directory.
        """
        # I build my initial/base config from the the default INI content
        ini_name = os.path.basename(filename)
        ini_name, _ = os.path.splitext(ini_name)
        default = default_ini_file_content.DEFAULT_INI_FILE_CONTENT[ini_name]        
        default = default.split("\n")
        configobj.ConfigObj.__init__(self, infile=default, encoding="utf-8")
        
        # Since I initialized my base class with a list instead of providing
        # a filename, self.filename is still None. That needs to be corrected.
        filename = os.path.join(util_misc.get_data_dir(), filename)
        self.filename = filename
        
        # Any settings the user provides trump the defaults
        user_config = configobj.ConfigObj(infile=filename, encoding="utf-8")
        
        self.merge(user_config)

    
    ##########    Public functions follow in get/set pairs    ##########

    def get_last_export_path(self):
        """Returns the last path from which the user exported a file via this
        application. The path is always fully-qualified, is never blank, and
        always points to an existing path. 
        """
        path = ""
        
        if "general" in self:
            path = self["general"].get("last_export_path", "")

        if not os.path.exists(path):
            path = util_misc.get_documents_dir()

        return path


    def set_last_export_path(self, path):
        """Sets the last path to which the user exported a file via this
        application.
        """
        path = os.path.abspath(path)
        if "general" not in self:
            self["general"] = { }
        self["general"]["last_export_path"] = path

    
    def get_last_file_open_path(self):
        """Returns the last path from which the user opened a file via this
        application. The path is always fully-qualified, is never blank, and
        always points to an existing path.
        """
        path = ""
        
        if "general" in self:
            path = self["general"].get("last_file_open_path", "")

        if not os.path.exists(path):
            path = util_misc.get_documents_dir()

        return path

    
    def set_last_file_open_path(self, path):
        """Sets the last path from which the user opened a file via this
        application.
        """
        path = os.path.abspath(path)
        if "general" not in self:
            self["general"] = { }
        self["general"]["last_file_open_path"] = path

    
    def get_window_coordinates(self, window_name):
        """A convenience function that returns a tuple of tuples representing
        (position, size) like so:
           (left, top), (width, height)
        These can be passed directly to a wxWidgets window's __init__()
        function.
        
        This function sanity checks the values in the INI file and adjusts
        them if they're out of range.
        """
        # I can't figure out why, but once in a while these numbers get 
        # written as non-whole numbers so I need to be prepared to handle 
        # not just ints but also floats.
        left = self[window_name].as_float("left")
        top = self[window_name].as_float("top")
        width = self[window_name].as_float("width")
        height = self[window_name].as_float("height")

        if wx is None:
            screen_width  = 1024
            screen_height = 720
        else:
            # I sanity check these numbers. No part of the window should be
            # off the screen
            screen_width = wx.SystemSettings.GetMetric(wx.SYS_SCREEN_X)
            screen_height = wx.SystemSettings.GetMetric(wx.SYS_SCREEN_Y)
                
        # Ensure left & top are non-negative
        left = max(left, 0)
        top = max(top, 0)

        # If the right edge of the window would be off the screen, decrease
        # left & width by 10% until that's no longer the case.
        while (left + width >= screen_width):
            left  *= 0.9
            width *= 0.9
            
        # Do the same for the bottom edge of the screen.
        while (top + height >= screen_height):
            top    *= 0.9
            height *= 0.9
        
        return ( (int(left), int(top)), (int(width), int(height)) )
    
    
    def set_window_coordinates(self, window_name, left, top, width, height):
        """A convenience function that sets window position and size."""
        self[window_name]["left"] = left
        self[window_name]["top"] = top
        self[window_name]["width"] = width
        self[window_name]["height"] = height
    
    
    def get_window_maximized(self, window_name):
        """Returns true if the window's maximize flag is set, False otherwise."""
        maximized = False
        
        # "maximize" might not be present, or it might contain a non-boolean
        # value. Neither of these are fatal errors.
        try:
            maximized = self[window_name].as_bool("maximized")
        except (KeyError, ValueError):
            pass
            
        return maximized


    def set_window_maximized(self, window_name, maximized):
        self[window_name]["maximized"] = bool(maximized)


    def write(self, outfile=None, section=None):
        """This is a thin wrapper on ConfigObj.write(). See that documentation
        for details. This adds sorting of keys within a section.
        """
        if not section:
            for section_ in list(self.values()):
                # Make a copy of the section, clear it, and then re-add the
                # keys in sorted order.
                d = section_.copy()
                section_.clear()
                for key in sorted(d.keys()):
                    section_[key] = d[key]

        return configobj.ConfigObj.write(self, outfile, section)


        

class VespaConfig(BaseConfig):
    def __init__(self):
        BaseConfig.__init__(self, "vespa.ini")


    @property
    def show_wx_inspector(self):
        return ("debug" in self) and                        \
               ("show_wx_inspector" in self["debug"]) and   \
               self["debug"].as_bool("show_wx_inspector")
