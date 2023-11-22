# Python modules

import abc

# 3rd party modules

# Our modules


"""There are three layers of Prefs classes. They are --

       Prefs
         |
         |
  App-specific Prefs (e.g. AnalysisPrefs, SimulationPrefs, etc.)
         |
         |
  Tab-specific Prefs (e.g. PrefsVoigt)         


You'll never need to instantiate classes from the first two layers. In fact, 
you can't -- they're abstract base classes. They define interfaces and provide
some common functions, but you don't need them directly.

Most of the interesting stuff happens at the bottom (Prefs) layer, which 
speaks to the fact that most of this code is the same. 

The app-specific Prefs classes contain very little. In fact, they exist only
to provide a property that tells the class which INI file to use.

The tab-specific prefs mostly exist to provide the _ini_section_name property
and inflate/deflate specifics.

BJS - additional comments 
----------------------------------

So I admit that I do not completely understand how these Prefs work.  So what
I write here is mostly empirical and some trial/error practicality.

Most of the  Prefs attributes are filled by turning the ViewIdsXxxx values in 
util_menu.py module from UPPER_CASE to lower_case and using these strings as
argument to setattr().  The ViewIdsXxxx values specific to the Tab being
instantiated are used (as listed in the __init__ call at the top of the tab's
__init__ code).  Yes, double inits.  

A second way to set a Prefs attribute is in the local prefs.py module in the
inflate() method for whatever tab (e.g. PrefsSpectral class). The attributes
listed in this location are NOT modified by a menu item, thus I consider them
to be 'static' prefs. Modifiable in the INI file, but not dynamically in the 
program.  But, can also be values other than Boolean, which is nice.

So, when should we add something to the common/default_ini_file_content?  By
trial and error, I deduced that you HAVE to add a value here if you are going
to put the prefs attribute in the inflate() method.  If you are setting the
value in the util_menu module, you can set a default here, or just let the 
menu handling code default on its own.   

BJS - additional comments 
----------------------------------

So I admit that I do not completely understand how these Prefs work.  So what
I write here is mostly empirical and some trial/error practicality.

Most of the  Prefs attributes are filled by turning the ViewIdsXxxx values in 
util_menu.py module from UPPER_CASE to lower_case and using these strings as
argument to setattr().  The ViewIdsXxxx values specific to the Tab being
instantiated are used (as listed in the __init__ call at the top of the tab's
__init__ code).  Yes, double inits.  

A second way to set a Prefs attribute is in the local prefs.py module in the
inflate() method for whatever tab (e.g. PrefsSpectral class). The attributes
listed in this location are NOT modified by a menu item, thus I consider them
to be 'static' prefs. Modifiable in the INI file, but not dynamically in the 
program.  But, can also be values other than Boolean, which is nice.

So, when should we add something to the common/default_ini_file_content?  By
trial and error, I deduced that you HAVE to add a value here if you are going
to put the prefs attribute in the inflate() method.  If you are setting the
value in the util_menu module, you can set a default here, or just let the 
menu handling code default on its own.   

"""


class Prefs(object, metaclass=abc.ABCMeta):
    """Prefs is an abstract base class for tracking user preferences. It 
    implements some work-saving magic that allows its subclasses to ignore
    most boolean & radio-type preferences reflected on the menu (e.g. show zero
    line, data type real/imag/magnitude, etc.).

    In particular, each instance of a Prefs object has a number of 
    "automatic" attributes. These attrs are not explicitly enumerated in the
    class definition. Rather, they are created by inference at class 
    instantiation based on the names of the constants in one of the 
    util_menu.ViewIdXxxx classes.
    """
    __metaclass__ = abc.ABCMeta


    def __init__(self, menu_bar, id_class):
        self._menu_bar = menu_bar

        # Determine the names of my auto attributes.
        self._auto_names = [name.lower() for name in id_class.BOOLEANS]
        # Create each auto name as an attr and initialize it to False.
        f = lambda name: setattr(self, name, False)
        list(map(f, self._auto_names))

        # In a few of my methods I need to map a menu item id to the name of
        # one of my automatic attributes. Here's where I build that map.
        self._id_name_map = { }

        f = lambda name: getattr(id_class, name.upper())
        self._id_name_map = dict( (f(name), name) for name 
                                                  in self._auto_names)

        # (Almost?) all of our INI sections have an entry for background
        # color, so we define it here so that it will get picked up during
        # inflate(). I initialize it to a poor default so that if this
        # default ever gets used by accident, it will be immediately clear.
        self.bgcolor = "black"

        config = self._ConfigClass()

        section = config[self._ini_section_name]

        self.inflate(section)


    @property
    @abc.abstractmethod
    def _ConfigClass(self):
        """Returns the appropriate ConfigObj class for this app (e.g.
        util_analysis_config.ConfigObj). Must be implemented by the AppPrefs
        class layer.
        """
        raise NotImplementedError


    @property
    @abc.abstractmethod
    def _ini_section_name(self):
        """Returns the INI file section name for these preferences, e.g.
        'spectral_prefs'. Read only. Subclasses *must* implement this property.
        """
        raise NotImplementedError


    @property
    def _state(self):
        """Returns a dict containing the true/false state of each of the
        automatic attributes. The keys are attribute names.
        """
        return dict((name, getattr(self, name)) for name in self._auto_names)


    @property
    def menu_state(self):
        """Returns a dict containing the true/false state of each of the
        automatic attributes. The keys are menu item ids.
        """
        state = { }
        for menu_item_id in self._id_name_map:
            name = self._id_name_map[menu_item_id]
            state[menu_item_id] = getattr(self, name)

        return state


    def __str__(self):
        return str(self).encode("utf-8")


    def __unicode__(self):
        """Returns state of all attrs (not just automatic ones) as 'attr: value'
        sorted by attribute name."""
        d = self.deflate()

        pairs = [ "%s: %s" % (key, d[key]) for key in sorted(d.keys()) ]

        return '\n'.join(pairs)


    def deflate(self):
        """Returns a dict containing all the attributes on this object
        (not just automatic ones). Attributes that start with underscore
        are ignored.

        Unlike many other objects in Vespa, there is no option for deflating
        prefs to XML.
        """
        return dict( (key, value) for key, value in self.__dict__.items()
                                                 if not key.startswith('_') )


    def inflate(self, source):
        """Given a dictionary, populates this object's attributes. If the 
        dict is a ConfigObj instance, automatic attributes will be 
        converted to boolean during inflation.

        Unlike many other objects in Vespa, there is no option for inflating
        prefs from XML.
        """
        is_configobj = (hasattr(source, "as_bool"))

        for key in source:
            if hasattr(self, key):
                if is_configobj and (key in self._auto_names):
                    value = source.as_bool(key)
                else:
                    value = source[key]

                setattr(self, key, value)


    def handle_event(self, menu_item_id):
        """Given the id of a menu item that has changed state, updates this
        object's attributes to reflect that of the menu. Returns True if 
        any attributes changed, False otherwise.
        """
        # Since many menu items are radio items, one click can often cause
        # two state changes as one item is turned on and another is turned off.
        # Rather than trying to figure out what changed, I use the brute force
        # approach. Every time this is called, I check the state of every 
        # boolean item on the menu. 

        # First, I determine which top level menu contains this item.
        menu = self._menu_bar.find_top_level_parent(menu_item_id)
        if menu:
            state = self._menu_bar.get_top_level_menu_state(menu)

            # state maps menu item ids to their boolean values. I want attribute 
            # names as keys, not menu item ids. Remap!
            state = dict( (self._id_name_map[key], value) for key, value 
                                                          in  state.items() )

            if state != self._state:
                # State change ==> one of the menu items that I track was selected.
                # Update my attributes.
                for name, value in state.items():
                    setattr(self, name, value)
                return True
            else:
                # Nothing changed.
                return False
        else:
            # menu_item_id refers to a menu item that I don't care about.
            return False


    def save(self):
        """Saves these preferences to the INI file."""
        config = self._ConfigClass()
        config[self._ini_section_name] = self.deflate()
        config.write()


