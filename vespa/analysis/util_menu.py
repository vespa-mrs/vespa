# Python modules

# 3rd party modules
import wx

# Our modules
import vespa.common.menu as common_menu
import vespa.common.util.config as util_common_config
import vespa.common.util.misc as util_misc



########################################################################
# Collection of menu-related constants, functions and utilities.
# The menu bar build function lives here, as does the menu definition.
########################################################################


class ViewIdsRaw(common_menu.IdContainer):
    """ Container for ids of menu items to which we need explicit references. """
    OUTPUT_HEADER_TEXT = "replace me"
    
    
class ViewIdsPrepFidsum(common_menu.IdContainer):
    """
    A container for the ids of all of the menu items to which we need
    explicit references.
    """
    ZERO_LINE_SHOW = "replace me"
    ZERO_LINE_TOP = "replace me"
    ZERO_LINE_MIDDLE = "replace me"
    ZERO_LINE_BOTTOM = "replace me"

    XAXIS_SHOW = "replace me"
    XAXIS_PPM = "replace me"
    XAXIS_HERTZ = "replace me"

    DATA_TYPE_REAL = "replace me"
    DATA_TYPE_IMAGINARY = "replace me"
    DATA_TYPE_MAGNITUDE = "replace me"

    AREA_CALC_PLOT_A = "replace me"
    AREA_CALC_PLOT_B = "replace me"
    
    ZERO_LINE_PLOT_SHOW = "replace me"
    ZERO_LINE_PLOT_TOP = "replace me"
    ZERO_LINE_PLOT_MIDDLE = "replace me"
    ZERO_LINE_PLOT_BOTTOM = "replace me"

    VIEW_TO_PNG  = "replace me"
    VIEW_TO_SVG  = "replace me"
    VIEW_TO_PDF  = "replace me"
    VIEW_TO_EPS  = "replace me"


class ViewIdsSpectral(common_menu.IdContainer):
    """
    A container for the ids of all of the menu items to which we need
    explicit references.
    """
    ZERO_LINE_SHOW = "replace me"
    ZERO_LINE_TOP = "replace me"
    ZERO_LINE_MIDDLE = "replace me"
    ZERO_LINE_BOTTOM = "replace me"

    XAXIS_SHOW = "replace me"
    XAXIS_PPM = "replace me"
    XAXIS_HERTZ = "replace me"

    DATA_TYPE_REAL = "replace me"
    DATA_TYPE_IMAGINARY = "replace me"
    DATA_TYPE_MAGNITUDE = "replace me"
    DATA_TYPE_SUMMED = "replace me"
    
    PLOT_C_FUNCTION_NONE = "replace me"
    PLOT_C_FUNCTION_A_MINUS_B = "replace me"
    PLOT_C_FUNCTION_B_MINUS_A = "replace me"
    PLOT_C_FUNCTION_A_PLUS_B = "replace me"
    
    AREA_CALC_PLOT_A = "replace me"
    AREA_CALC_PLOT_B = "replace me"
    AREA_CALC_PLOT_C = "replace me"

    USER_BUTTON_PHASING = "replace me"
    USER_BUTTON_AREA = "replace me"

    VIEW_TO_PNG  = "replace me"
    VIEW_TO_SVG  = "replace me"
    VIEW_TO_PDF  = "replace me"
    VIEW_TO_EPS  = "replace me"
    
    SVD_DATA_TO_VIFF = "replace me"
    SVD_PEAKS_CHECKED_SUM_TO_VIFF = "replace me"
    SVD_FIDS_CHECKED_SUM_TO_VIFF = "replace me"
    SVD_PEAKS_CHECKED_DIFF_TO_VIFF = "replace me"
    SVD_TABLE_VALUES_TO_XML = "replace me"
    SVD_TABLE_VALUES_TO_CSV = "replace me"

    PLOTS_TO_BINARY = "replace me"
    PLOTS_TO_ASCII = "replace me"


class ViewIdsPrepTimeseries(common_menu.IdContainer):
    """
    A container for the ids of all of the menu items to which we need
    explicit references.
    """
    ZERO_LINE_SHOW = "replace me"
    ZERO_LINE_TOP = "replace me"
    ZERO_LINE_MIDDLE = "replace me"
    ZERO_LINE_BOTTOM = "replace me"

    XAXIS_SHOW  = "replace me"
    XAXIS_PPM   = "replace me"
    XAXIS_HERTZ = "replace me"

    DATA_TYPE_REAL      = "replace me"
    DATA_TYPE_IMAGINARY = "replace me"
    DATA_TYPE_MAGNITUDE = "replace me"

    AREA_CALC_PLOT_A = "replace me"
    AREA_CALC_PLOT_B = "replace me"
    

class ViewIdsVoigt(common_menu.IdContainer):
    ZERO_LINE_SHOW   = "replace me"
    ZERO_LINE_TOP    = "replace me"
    ZERO_LINE_MIDDLE = "replace me"
    ZERO_LINE_BOTTOM = "replace me"

    XAXIS_SHOW  = "replace me"
    XAXIS_PPM   = "replace me"
    XAXIS_HERTZ = "replace me"

    N_PLOTS_1 = "replace me"
    N_PLOTS_2 = "replace me"
    N_PLOTS_3 = "replace me"
    N_PLOTS_4 = "replace me"

    AREA_CALC_PLOT_A = "replace me"
    AREA_CALC_PLOT_B = "replace me"
    AREA_CALC_PLOT_C = "replace me"
    AREA_CALC_PLOT_D = "replace me"

    VIEW_TO_PNG = "replace me"
    VIEW_TO_SVG = "replace me"
    VIEW_TO_PDF = "replace me"
    VIEW_TO_EPS = "replace me"
    
    RESULTS_CSV_CURRENT   = "replace me"
    RESULTS_CSV_ALL       = "replace me"
    RESULTS_LCM_PDF       = "replace me"
    RESULTS_LCM_PDF_MULTI = "replace me"
    RESULTS_LCM_PNG       = "replace me"
    RESULTS_BRP512_PNG    = "replace me"
    RESULTS_BRP1024_PNG   = "replace me"
    RESULTS_2PLOT_PDF     = "replace me"
    RESULTS_2PLOT_PNG     = "replace me"
    RESULTS_4PLOT_PDF     = "replace me"
    RESULTS_4PLOT_PNG     = "replace me"

    DEBUG_POPUP_CURRENT = "replace me"
    DEBUG_CSV_CURRENT   = "replace me"


class ViewIdsGiso(common_menu.IdContainer):
    ZERO_LINE_SHOW   = "replace me"
    ZERO_LINE_TOP    = "replace me"
    ZERO_LINE_MIDDLE = "replace me"
    ZERO_LINE_BOTTOM = "replace me"

    XAXIS_SHOW  = "replace me"
    XAXIS_PPM   = "replace me"
    XAXIS_HERTZ = "replace me"

    N_PLOTS_1 = "replace me"
    N_PLOTS_2 = "replace me"
    N_PLOTS_3 = "replace me"
    N_PLOTS_4 = "replace me"

    AREA_CALC_PLOT_A = "replace me"
    AREA_CALC_PLOT_B = "replace me"
    AREA_CALC_PLOT_C = "replace me"
    AREA_CALC_PLOT_D = "replace me"

    VIEW_TO_PNG = "replace me"
    VIEW_TO_SVG = "replace me"
    VIEW_TO_PDF = "replace me"
    VIEW_TO_EPS = "replace me"
    
    RESULTS_LCM_PDF   = "replace me"
    RESULTS_LCM_PNG   = "replace me"
    RESULTS_2PLOT_PDF = "replace me"
    RESULTS_2PLOT_PNG = "replace me"
    RESULTS_4PLOT_PDF = "replace me"
    RESULTS_4PLOT_PNG = "replace me"

class ViewIdsWatref(common_menu.IdContainer):
    """A container for the ids of all of the menu items to which we need
    explicit references.
    """
    #OUTPUT_RESULTS_TEXT     = "replace me"
    RESULTS_CSV_CURRENT       = "replace me"
    RESULTS_CSV_ALL           = "replace me"
    RESULTS_LCM_PDF           = "replace me"
    RESULTS_LCM_PDF_MULTI     = "replace me"
    RESULTS_LCM_PNG           = "replace me"
        
    


class PlotXIds(common_menu.IdContainer):
    @classmethod
    def enumerate_booleans(klass, menu):
        # Same docstring as base class
        __doc__ == common_menu.IdContainer.enumerate_booleans.__doc__

        # This code is the same as the base class implementation except that
        # it's aware that the constants in this class contain a list of 
        # ids rather than just a single id.
        klass.BOOLEANS = [ ]

        # as of wxPython 3.0.0 (?) wx ids have their own class wx.WindowIDRef
        # which has 'value' and 'Id' attributes, so we can scan the whole 
        # util_menu object to see what items it has that are in ALL_CAPS

        for attribute_name, value in klass.__dict__.items():
            # Make sure this is one of my attributes.
            if attribute_name.isupper() and util_misc.is_iterable(value) and \
               isinstance(value[0], wx.WindowIDRef):
                # OK, it's one of mine. See if it is boolean.
                item = menu.FindItemById(value[0])
                if item.IsCheckable():
                    klass.BOOLEANS.append(attribute_name)

        # I make them a tuple just to emphasize that they're static.
        klass.BOOLEANS = tuple(klass.BOOLEANS)


    @classmethod
    def init_ids(klass):
        # Same docstring as base class
        __doc__ == common_menu.IdContainer.init_ids.__doc__

        # This code is the same as the base class implementation except that
        # it's aware that the constants in this class contain a list of 
        # ids rather than just a single id.
        for key, value in klass.__dict__.items():
            if key.isupper() and (value == "replace me"):
                # It's one of ours
                # PLOT_X menu items get a tuple of ids, one for each plot
                # visible on Voigt. Note that these must be tuples, not
                # lists. That's because the Prefs class uses them as
                # dictionary keys, and tuples are hashable (because they
                # can't be modified) while lists are not.
                value = tuple([wx.NewIdRef() for i in range(4)])

                setattr(klass, key, value)


    RAW_DATA                        = "replace me"
    FITTED_DATA                     = "replace me"
    BASELINE                        = "replace me"
    COMBO_RAW_MINUS_FIT             = "replace me"
    COMBO_RAW_MINUS_BASE            = "replace me"
    COMBO_RAW_MINUS_FIT_MINUS_BASE  = "replace me"
    COMBO_RAW_AND_FIT               = "replace me"
    COMBO_RAW_AND_BASE              = "replace me"
    COMBO_RAW_AND_FIT_PLUS_BASE     = "replace me"
    COMBO_RAW_MINUS_BASE_AND_FIT    = "replace me"
    COMBO_FIT_PLUS_BASE             = "replace me"
    COMBO_FIT_AND_BASE              = "replace me"
    COMBO_RAW_AND_INIT_MODEL        = "replace me"
    COMBO_RAW_AND_WT_ARR            = "replace me"
    DATA_TYPE_REAL                  = "replace me"
    DATA_TYPE_IMAGINARY             = "replace me"
    DATA_TYPE_MAGNITUDE             = "replace me"
    DATA_TYPE_SUMMED                = "replace me"



# When main creates an instance of AnalysisMenuBar(), it sets the variable
# below to that instance. It's a convenience. It's the same as 
# wx.GetApp().GetTopWindow().GetMenuBar(), but much easier to type.
bar = None

# _PLOT_MENU_INSERT_INDEX tells us where to start inserting the 'Plot X' menu
# items when voigt or giso are shown. 
_PLOT_MENU_INSERT_INDEX = 3


class _AnalysisMenu(wx.Menu):
    """
    A menu item that knows its own label (standard wx.Menus don't).
    
    Code outside of this module shouldn't need to create instances of this.
    """
    def __init__(self, label=""):
        wx.Menu.__init__(self)
        
        self.label = label


class AnalysisMenuBar(common_menu.VespaMenuBar):
    """
    A subclass of wx.MenuBar that adds some Analysis-specific functions
    and constants.
    
    The TYPE_XXX constants are for the show_menus() method.

    Analysis only needs one instance of this class per invocation. It's a 
    singleton class.
    """
    DEVELOP_FLAG            = True
    
    TYPE_NONE               = ""
    TYPE_START              = "start"
    TYPE_RAW                = "raw"
    TYPE_PREP_FIDSUM        = "prep_fidsum"
    TYPE_PREP_TIMESERIES    = "prep_timeseries"
    TYPE_SPECTRAL           = "spectral"
    TYPE_VOIGT              = "voigt"
    TYPE_GISO               = "giso"
    TYPE_WATREF             = "watref"
    
    def __init__(self, main, import_items):
        common_menu.VespaMenuBar.__init__(self, main)
        
        ViewIdsRaw.init_ids()
        ViewIdsPrepFidsum.init_ids()
        ViewIdsPrepTimeseries.init_ids()
        ViewIdsSpectral.init_ids()
        ViewIdsVoigt.init_ids()
        ViewIdsGiso.init_ids()
        ViewIdsWatref.init_ids()

        PlotXIds.init_ids()

        # self._current tracks which menu type is currently visible.
        self._current = None
        self._main = main

        # _get_menu_data() is called just once, right here. 
        vals = _get_menu_data(main)
        
        file_, processing, view_start, view_raw, view_prep_fidsum,  \
        view_prep_timeseries, view_spectral, view_voigt, view_giso, \
        view_watref, plot_menus, help = vals

        # We save some of the menu data so we can re-create the view menus
        # on the fly.
        self._menu_data = { }

        self._menu_data[self.TYPE_START]            = view_start
        self._menu_data[self.TYPE_RAW]              = view_raw
        self._menu_data[self.TYPE_PREP_FIDSUM]      = view_prep_fidsum
        self._menu_data[self.TYPE_PREP_TIMESERIES]  = view_prep_timeseries
        self._menu_data[self.TYPE_SPECTRAL]         = view_spectral
        self._menu_data[self.TYPE_VOIGT]            = view_voigt
        self._menu_data[self.TYPE_GISO]             = view_giso
        self._menu_data[self.TYPE_WATREF]           = view_watref

        # Here we construct the 'Plot X' menus and cache them in 
        # self._plot_menus. At this point we just build them, we don't show
        # them yet.
        self._plot_menus = [ ]
        ASCII_A = ord('A')
        for i, menu in enumerate(plot_menus):
            self._plot_menus.append(_create_menu(main, "&Plot %s" % chr(ASCII_A + i), menu))

        # Build the four top-level menus that are always present. Note that 
        # view is created empty.
        file_       = _create_menu(main, "&File", file_)
        processing  = _create_menu(main, "&Processing", processing)
        view        = _create_menu(main, "&View", { })
        help        = _create_menu(main, "&Help", help)

        for menu in (file_, processing, view, help):
            self.Append(menu, menu.label)

        # I need to enumerate the booleans on the constant classes. I can't
        # do this until the menus are created. The view menus are created only
        # as needed. In practice, I need to populate the BOOLEANS attribute
        # before that, so I create a set of dummy menus here just so I can
        # enumerate the booleans.
        for menu_definition, IdClass in ( (view_raw,             ViewIdsRaw),
                                          (view_prep_fidsum,     ViewIdsPrepFidsum),
                                          (view_prep_timeseries, ViewIdsPrepTimeseries ),
                                          (view_spectral,        ViewIdsSpectral),
                                          (view_voigt,           ViewIdsVoigt),
                                          (view_giso,            ViewIdsGiso),
                                          (view_watref,          ViewIdsWatref),
                                        ):
            dummy = _create_menu(main, "Dummy", menu_definition)
            IdClass.enumerate_booleans(dummy)

        # The Plot X menus are already created, so I don't have to create 
        # dummies for them.
        PlotXIds.enumerate_booleans(self._plot_menus[0])

        self._user_menu_items = { }

        for module_name in list(import_items.keys()):

            # A top level keyword with the string 'separator' in it will cause
            # a menu separator to be added to the menu            
            if 'separator' in module_name.lower():
                self.file_import_menu.AppendSeparator()
            else:

                section = import_items[module_name]
    
                # Here's what a section looks like:
                #    [acme]
                #    path=/Users/fred/acme.py
                #    class_name=RawReaderAcme
                #    menu_item_text=Acme
                #    ini_file_name=import_acme

                path = section["path"]
                    
                # Replace '\t' in the text with an actual tab (ASCII 9). This
                # allows users to define accelerators so that they can access 
                # their menu item with a single keystroke.
                menu_item_text = section["menu_item_text"]
                menu_item_text = menu_item_text.replace("\\t", "\t")

                id_ = wx.NewIdRef()
                menu_item = wx.MenuItem(self.file_import_menu, id_, text=menu_item_text)
                self.file_import_menu.Append(menu_item)     # bjs_ccx 

                # Save the module_name & INI file name associated with this menu item.
                self._user_menu_items[id_] = (module_name, section["ini_file_name"])

                main.Bind(wx.EVT_MENU, main.on_import_user_item, menu_item)


        self.show_menus(self.TYPE_START)
        


    @property
    def file_menu(self):
        """A convenience property that returns the File menu"""
        return self.GetMenu(self.FindMenu("File"))


    @property
    def file_import_menu(self):
        """A convenience property that returns the File/Import submenu"""
        id_ = self.file_menu.FindItem("Import")
        return self.file_menu.FindItemById(id_).GetSubMenu()


    def get_user_menu_item_info(self, menu_id):
        """Given the wx id of a menu item defined via 
        analysis_import_menu_additions.ini, returns a 2-tuple of 
        (class, ini_name) where class is the RawReader class to be 
        instantiated, and ini_name is the analysis.ini file key for 
        reading/saving the last path used when opening this particular format
        """
        return self._user_menu_items[menu_id]


    def show_menus(self, type_):
        """
        Given one of the type constants (TYPE_RAW, TYPE_SPECTRAL, 
        TYPE_VOIGT, TYPE_GISO), shows/hides the appropriate menus for the mode.
        """
        if type_ == self._current:
            # do nothing
            pass
        else:
            if self._current == self.TYPE_VOIGT:
                # Plot menus are visible; hide them.
                plot_menu_labels = [menu.label for menu in self._plot_menus]

                for menu in self.top_level_menus:
                    if menu.label in plot_menu_labels:
                        self.Remove(self.FindMenu(menu.label))
            elif self._current == self.TYPE_GISO:
                # Plot menus are visible; hide them.
                plot_menu_labels = [menu.label for menu in self._plot_menus]

                for menu in self.top_level_menus:
                    if menu.label in plot_menu_labels:
                        self.Remove(self.FindMenu(menu.label))

            # Rebuild the view menu by deleting everything from it and then 
            # reappending the appropriate items.
            while self.view_menu.GetMenuItemCount():
                #self.view_menu.DeleteItem(self.view_menu.FindItemByPosition(0))
                self.view_menu.Delete(self.view_menu.FindItemByPosition(0))

            _append_items(self._main, self.view_menu, self._menu_data[type_])

            if type_ == self.TYPE_VOIGT:
                # add plot menus
                for menu in self._plot_menus[::-1]:
                    self.Insert(_PLOT_MENU_INSERT_INDEX, menu, menu.label)
                    # Under wxPython 2.9, the menus I add with this call to 
                    # Insert() don't have their label set. I think it's a bug,
                    # but I can't recreate it outside of this app. Manually
                    # setting the label here is a workaround.
                    self.SetMenuLabel(_PLOT_MENU_INSERT_INDEX, menu.label)
            elif type_ == self.TYPE_GISO:
                # add plot menus
                for menu in self._plot_menus[::-1]:
                    self.Insert(_PLOT_MENU_INSERT_INDEX, menu, menu.label)
                    # Under wxPython 2.9, the menus I add with this call to 
                    # Insert() don't have their label set. I think it's a bug,
                    # but I can't recreate it outside of this app. Manually
                    # setting the label here is a workaround.
                    self.SetMenuLabel(_PLOT_MENU_INSERT_INDEX, menu.label)


        self._current = type_
        


# ================    Module Internal Use Only    =======================


def _append_items(main, menu, items):
    """Given a reference to main, a top-level menu, and dict of items, turns
    the latter into menu items and appends them to the top-level menu."""
    for item in items:
        # Item can be a separator, a tuple representing a simple menu item
        # (like the Experiment/Open menu item) or a tuple representing a
        # menu item with a submenu (like View/Plot Color/...).

        # In the second case, the tuple is 2-4 items long, where the items
        # are (label, handler, kind, id).

        # In the third case, the tuple is always 2 items long. The items
        # are a label and a list of submenu items.
        if item == common_menu.SEPARATOR:
            menu.AppendSeparator()
        else:
            if callable(item[1]):
                # This is a regular menu item because it has an event 
                # handler (method) as its second element.
                common_menu.create_menu_item(main, menu, *item)
            else:
                # This is a menu item with subitems
                label, subitems = item
                submenu = _create_menu(main, "", subitems)
                menu.Append(wx.ID_ANY, label, submenu)       


def _create_menu(main, label, items):
    menu = _AnalysisMenu(label)

    _append_items(main, menu, items)

    return menu
       

def _get_menu_data(main):
    # Note that wx treats the ids wx.ID_EXIT and wx.ID_ABOUT specially by 
    # moving them to their proper location on the Mac. wx will also change
    # the text of the ID_EXIT item to "Quit" as is standard under OS X. 
    # Quit is also the standard under Gnome but unfortunately wx doesn't seem
    # to change Exit --> Quit there, so our menu looks a little funny under
    # Gnome.
    file_ = (
            ("&Open\tCTRL+O",       main.on_open_viff),        
            ("&Import", (
                ("VIFF Raw Data (*.xml)",     main.on_import_mrs_data_raw),
            )),
            common_menu.SEPARATOR,
            ("S&ave\tCTRL+S",       main.on_save_viff),       
            ("Save As...",          main.on_save_as_viff),       
            ("Close\tCTRL+W",       main.on_close_dataset),
            ("Close All",           main.on_close_all),
            common_menu.SEPARATOR,
            ("Presets...", (
                ("Load from File ",  main.on_load_preset_from_file),       
#                ("Load from Tab",    main.on_load_preset_from_tab),       
                ("Save to File",     main.on_save_preset_to_file),
            )),
            common_menu.SEPARATOR,
            ("Exit",       main.on_close_window, wx.ITEM_NORMAL, wx.ID_EXIT),
            )

    processing = (
            ("Add Fitting Tab - Voigt", main.on_add_voigt_tab, wx.ITEM_NORMAL),
#            ("Add Fitting Tab - Giso",  main.on_add_giso_tab, wx.ITEM_NORMAL),
            ("Add WatRef Quant Tab",    main.on_add_watref_tab, wx.ITEM_NORMAL),
            common_menu.SEPARATOR,
            ("Edit User Spectral Information...",   main.on_user_prior,           wx.ITEM_NORMAL),
            ("Edit User Metabolite Information...", main.on_user_metabolite_info, wx.ITEM_NORMAL),
                  )

    view_start = tuple()

    view_raw = (
            ("Output", (
                ("View Selected Header -> Text", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIdsRaw.OUTPUT_HEADER_TEXT),
                       )
            ),
           )
    

    view_prep_fidsum = (
            ("Zero Line - Spectral", (
                ("Show",   main.on_menu_view_option, wx.ITEM_CHECK, ViewIdsPrepFidsum.ZERO_LINE_SHOW),
                common_menu.SEPARATOR,
                ("Top",    main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsPrepFidsum.ZERO_LINE_TOP),
                ("Middle", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsPrepFidsum.ZERO_LINE_MIDDLE),
                ("Bottom", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsPrepFidsum.ZERO_LINE_BOTTOM))),
            ("X-Axis - Spectral", (
                ("Show", main.on_menu_view_option,   wx.ITEM_CHECK, ViewIdsPrepFidsum.XAXIS_SHOW),
                common_menu.SEPARATOR,
                ("PPM",   main.on_menu_view_option,  wx.ITEM_RADIO, ViewIdsPrepFidsum.XAXIS_PPM),
                ("Hertz", main.on_menu_view_option,  wx.ITEM_RADIO, ViewIdsPrepFidsum.XAXIS_HERTZ))),
            common_menu.SEPARATOR,
            ("Data Type - Spectral", (
                ("Real",      main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsPrepFidsum.DATA_TYPE_REAL),
                ("Imaginary", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsPrepFidsum.DATA_TYPE_IMAGINARY),
                ("Magnitude", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsPrepFidsum.DATA_TYPE_MAGNITUDE),
                                     )
            ),
            common_menu.SEPARATOR,
            ("Area Calculation - Spectral", (
                ("From Plot A", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsPrepFidsum.AREA_CALC_PLOT_A),
                ("From Plot B", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsPrepFidsum.AREA_CALC_PLOT_B),
                                            )
            ),
            common_menu.SEPARATOR,
            ("Zero Line - Series", (
                ("Show",   main.on_menu_view_option, wx.ITEM_CHECK, ViewIdsPrepFidsum.ZERO_LINE_PLOT_SHOW),
                common_menu.SEPARATOR,
                ("Top",    main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsPrepFidsum.ZERO_LINE_PLOT_TOP),
                ("Middle", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsPrepFidsum.ZERO_LINE_PLOT_MIDDLE),
                ("Bottom", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsPrepFidsum.ZERO_LINE_PLOT_BOTTOM),
                                   )
            ),
            common_menu.SEPARATOR,
            ("Output", (
                ("View -> PNG", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIdsPrepFidsum.VIEW_TO_PNG),
                ("View -> SVG", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIdsPrepFidsum.VIEW_TO_SVG),
                ("View -> PDF", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIdsPrepFidsum.VIEW_TO_PDF),
                ("View -> EPS", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIdsPrepFidsum.VIEW_TO_EPS),
            )),


           )      
           
    view_prep_timeseries = (
            ("Zero Line", (
                ("Show",   main.on_menu_view_option, wx.ITEM_CHECK, ViewIdsPrepTimeseries.ZERO_LINE_SHOW),
                common_menu.SEPARATOR,
                ("Top",    main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsPrepTimeseries.ZERO_LINE_TOP),
                ("Middle", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsPrepTimeseries.ZERO_LINE_MIDDLE),
                ("Bottom", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsPrepTimeseries.ZERO_LINE_BOTTOM))),
            ("X-Axis", (
                ("Show", main.on_menu_view_option,   wx.ITEM_CHECK, ViewIdsPrepTimeseries.XAXIS_SHOW),
                common_menu.SEPARATOR,
                ("PPM",   main.on_menu_view_option,  wx.ITEM_RADIO, ViewIdsPrepTimeseries.XAXIS_PPM),
                ("Hertz", main.on_menu_view_option,  wx.ITEM_RADIO, ViewIdsPrepTimeseries.XAXIS_HERTZ))),
            common_menu.SEPARATOR,
            ("Data Type", (
                ("Real",      main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsPrepTimeseries.DATA_TYPE_REAL),
                ("Imaginary", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsPrepTimeseries.DATA_TYPE_IMAGINARY),
                ("Magnitude", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsPrepTimeseries.DATA_TYPE_MAGNITUDE),
                          )
            ),
            common_menu.SEPARATOR,
            ("Area Calculation", (
                ("From Plot A", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsPrepTimeseries.AREA_CALC_PLOT_A),
                ("From Plot B", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsPrepTimeseries.AREA_CALC_PLOT_B),
                                 )
            )
           )                 
        
    view_spectral = (
            ("Zero Line", (
                ("Show",   main.on_menu_view_option, wx.ITEM_CHECK, ViewIdsSpectral.ZERO_LINE_SHOW),
                common_menu.SEPARATOR,
                ("Top",    main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsSpectral.ZERO_LINE_TOP),
                ("Middle", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsSpectral.ZERO_LINE_MIDDLE),
                ("Bottom", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsSpectral.ZERO_LINE_BOTTOM))),
            ("X-Axis", (
                ("Show", main.on_menu_view_option,   wx.ITEM_CHECK, ViewIdsSpectral.XAXIS_SHOW),
                common_menu.SEPARATOR,
                ("PPM",   main.on_menu_view_option,  wx.ITEM_RADIO, ViewIdsSpectral.XAXIS_PPM),
                ("Hertz", main.on_menu_view_option,  wx.ITEM_RADIO, ViewIdsSpectral.XAXIS_HERTZ))),
            common_menu.SEPARATOR,
            ("Data Type", (
                ("Real",      main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsSpectral.DATA_TYPE_REAL),
                ("Imaginary", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsSpectral.DATA_TYPE_IMAGINARY),
                ("Magnitude", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsSpectral.DATA_TYPE_MAGNITUDE),
                common_menu.SEPARATOR,
                ("Summed",    main.on_menu_view_option, wx.ITEM_CHECK, ViewIdsSpectral.DATA_TYPE_SUMMED),
                          )
            ),
            common_menu.SEPARATOR,
            ("Area Calculation", (
                ("From Plot A", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsSpectral.AREA_CALC_PLOT_A),
                ("From Plot B", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsSpectral.AREA_CALC_PLOT_B),
                ("From Plot C", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsSpectral.AREA_CALC_PLOT_C))),
            ("Plot C Function", (
                ("None",         main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsSpectral.PLOT_C_FUNCTION_NONE),
                ("Residual A-B", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsSpectral.PLOT_C_FUNCTION_A_MINUS_B),
                ("Residual B-A", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsSpectral.PLOT_C_FUNCTION_B_MINUS_A),
                ("Sum A+B",      main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsSpectral.PLOT_C_FUNCTION_A_PLUS_B))),
            ("User Button Function", (
                ("Automatic Phasing", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsSpectral.USER_BUTTON_PHASING),
                ("Output Area Value", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsSpectral.USER_BUTTON_AREA))),
            common_menu.SEPARATOR,
            ("Output", (
                ("View -> PNG", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIdsSpectral.VIEW_TO_PNG),
                ("View -> SVG", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIdsSpectral.VIEW_TO_SVG),
                ("View -> PDF", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIdsSpectral.VIEW_TO_PDF),
                ("View -> EPS", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIdsSpectral.VIEW_TO_EPS),
                common_menu.SEPARATOR,
                ("Plot A", (
                    ("to Binary", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIdsSpectral.PLOTS_TO_BINARY),
                    ("to ASCII",  main.on_menu_view_output, wx.ITEM_NORMAL, ViewIdsSpectral.PLOTS_TO_ASCII))),
                common_menu.SEPARATOR,
                ("SVD Exports", (
                    ("Plot A Spectrum -> VIFF Raw Data", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIdsSpectral.SVD_DATA_TO_VIFF),
                    ("Plot B Spectrum Sum -> VIFF Raw Data",  main.on_menu_view_output, wx.ITEM_NORMAL, ViewIdsSpectral.SVD_PEAKS_CHECKED_SUM_TO_VIFF),
                    ("Plot B FID Data Sum -> VIFF Raw Data", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIdsSpectral.SVD_FIDS_CHECKED_SUM_TO_VIFF),
                    ("Plot C Spectrum Diff -> VIFF Raw Data", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIdsSpectral.SVD_PEAKS_CHECKED_DIFF_TO_VIFF),
                    ("HLSVD Table Values -> VIFF XML", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIdsSpectral.SVD_TABLE_VALUES_TO_XML),
                    ("HLSVD Table Values -> CSV File", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIdsSpectral.SVD_TABLE_VALUES_TO_CSV),
                                )
                 ),
            )))

    view_voigt = (
            ("Zero Line", (
                ("Show",   main.on_menu_view_option, wx.ITEM_CHECK, ViewIdsVoigt.ZERO_LINE_SHOW),
                common_menu.SEPARATOR,
                ("Top",    main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsVoigt.ZERO_LINE_TOP),
                ("Middle", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsVoigt.ZERO_LINE_MIDDLE),
                ("Bottom", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsVoigt.ZERO_LINE_BOTTOM))),
            ("X-Axis", (
                ("Show", main.on_menu_view_option,   wx.ITEM_CHECK, ViewIdsVoigt.XAXIS_SHOW),
                common_menu.SEPARATOR,
                ("PPM",   main.on_menu_view_option,  wx.ITEM_RADIO, ViewIdsVoigt.XAXIS_PPM),
                ("Hertz", main.on_menu_view_option,  wx.ITEM_RADIO, ViewIdsVoigt.XAXIS_HERTZ))),
            ("Number of Plots", (
                ("1", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsVoigt.N_PLOTS_1),
                ("2", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsVoigt.N_PLOTS_2),
                ("3", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsVoigt.N_PLOTS_3),
                ("4", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsVoigt.N_PLOTS_4))),
            common_menu.SEPARATOR,
            ("Area Calculation", (
                ("From Plot A", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsVoigt.AREA_CALC_PLOT_A),
                ("From Plot B", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsVoigt.AREA_CALC_PLOT_B),
                ("From Plot C", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsVoigt.AREA_CALC_PLOT_C),
                ("From Plot D", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsVoigt.AREA_CALC_PLOT_D))),
            ("Output", (
                ("View -> PNG", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIdsVoigt.VIEW_TO_PNG),
                ("View -> SVG", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIdsVoigt.VIEW_TO_SVG),
                ("View -> PDF", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIdsVoigt.VIEW_TO_PDF),
                ("View -> EPS", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIdsVoigt.VIEW_TO_EPS))),
            ("Results to File", (
                ("Text-CSV Layout", (
                  ("Current Voxel", main.on_menu_view_results, wx.ITEM_NORMAL, ViewIdsVoigt.RESULTS_CSV_CURRENT),
                  ("All Voxels",    main.on_menu_view_results, wx.ITEM_NORMAL, ViewIdsVoigt.RESULTS_CSV_ALL))),
                ("LCM Layout", (
                  ("PDF Format (Ph0/1 corr)",            main.on_menu_view_results, wx.ITEM_NORMAL, ViewIdsVoigt.RESULTS_LCM_PDF),
                  ("PDF Format Multi-page (Ph0/1 corr)", main.on_menu_view_results, wx.ITEM_NORMAL, ViewIdsVoigt.RESULTS_LCM_PDF_MULTI),
                  ("PNG Format (Ph0/1 corr)",            main.on_menu_view_results, wx.ITEM_NORMAL, ViewIdsVoigt.RESULTS_LCM_PNG))),
                ("Analysis BRP Layout", (
                  ("PNG 512 Format (Ph0/1 corr)",  main.on_menu_view_results, wx.ITEM_NORMAL, ViewIdsVoigt.RESULTS_BRP512_PNG),
                  ("PNG 1024 Format (Ph0/1 corr)", main.on_menu_view_results, wx.ITEM_NORMAL, ViewIdsVoigt.RESULTS_BRP1024_PNG))),
                ("Analysis 2-Plot Layout", (
                  ("PDF Format (Ph0/1 corr)", main.on_menu_view_results, wx.ITEM_NORMAL, ViewIdsVoigt.RESULTS_2PLOT_PDF),
                  ("PNG Format (Ph0/1 corr)", main.on_menu_view_results, wx.ITEM_NORMAL, ViewIdsVoigt.RESULTS_2PLOT_PNG))),
                ("Analysis 4-Plot Layout", (
                  ("PDF Format (Ph0/1 corr)", main.on_menu_view_results, wx.ITEM_NORMAL, ViewIdsVoigt.RESULTS_4PLOT_PDF),
                  ("PNG Format (Ph0/1 corr)", main.on_menu_view_results, wx.ITEM_NORMAL, ViewIdsVoigt.RESULTS_4PLOT_PNG)))
                )
            ),
            ("Debug", (
                ("Current Popup",     main.on_menu_view_debug, wx.ITEM_NORMAL, ViewIdsVoigt.DEBUG_POPUP_CURRENT),
                ("Current->Text-CSV", main.on_menu_view_debug, wx.ITEM_NORMAL, ViewIdsVoigt.DEBUG_CSV_CURRENT))),

           )

    view_giso = (
            ("Zero Line", (
                ("Show",   main.on_menu_view_option, wx.ITEM_CHECK, ViewIdsGiso.ZERO_LINE_SHOW),
                common_menu.SEPARATOR,
                ("Top",    main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsGiso.ZERO_LINE_TOP),
                ("Middle", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsGiso.ZERO_LINE_MIDDLE),
                ("Bottom", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsGiso.ZERO_LINE_BOTTOM))),
            ("X-Axis", (
                ("Show", main.on_menu_view_option,   wx.ITEM_CHECK, ViewIdsGiso.XAXIS_SHOW),
                common_menu.SEPARATOR,
                ("PPM",   main.on_menu_view_option,  wx.ITEM_RADIO, ViewIdsGiso.XAXIS_PPM),
                ("Hertz", main.on_menu_view_option,  wx.ITEM_RADIO, ViewIdsGiso.XAXIS_HERTZ))),
            ("Number of Plots", (
                ("1", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsGiso.N_PLOTS_1),
                ("2", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsGiso.N_PLOTS_2),
                ("3", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsGiso.N_PLOTS_3),
                ("4", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsGiso.N_PLOTS_4))),
            common_menu.SEPARATOR,
            ("Area Calculation", (
                ("From Plot A", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsGiso.AREA_CALC_PLOT_A),
                ("From Plot B", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsGiso.AREA_CALC_PLOT_B),
                ("From Plot C", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsGiso.AREA_CALC_PLOT_C),
                ("From Plot D", main.on_menu_view_option, wx.ITEM_RADIO, ViewIdsGiso.AREA_CALC_PLOT_D))),
            ("Output", (
                ("View -> PNG", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIdsGiso.VIEW_TO_PNG),
                ("View -> SVG", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIdsGiso.VIEW_TO_SVG),
                ("View -> PDF", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIdsGiso.VIEW_TO_PDF),
                ("View -> EPS", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIdsGiso.VIEW_TO_EPS))),
            ("Results to File", (
                ("LCM Layout", (
                  ("PDF Format", main.on_menu_view_results, wx.ITEM_NORMAL, ViewIdsGiso.RESULTS_LCM_PDF),
                  ("PNG Format", main.on_menu_view_results, wx.ITEM_NORMAL, ViewIdsGiso.RESULTS_LCM_PNG))),
                ("Analysis 2-Plot Layout", (
                  ("PDF Format", main.on_menu_view_results, wx.ITEM_NORMAL, ViewIdsGiso.RESULTS_2PLOT_PDF),
                  ("PNG Format", main.on_menu_view_results, wx.ITEM_NORMAL, ViewIdsGiso.RESULTS_2PLOT_PNG))),
                ("Analysis 4-Plot Layout", (
                  ("PDF Format", main.on_menu_view_results, wx.ITEM_NORMAL, ViewIdsGiso.RESULTS_4PLOT_PDF),
                  ("PNG Format", main.on_menu_view_results, wx.ITEM_NORMAL, ViewIdsGiso.RESULTS_4PLOT_PNG))))),

           )
           
    view_watref = (
            ("Results to File", (
                # ("Text-CSV Layout -> Current Voxel", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIdsWatref.OUTPUT_CURRENT_VOXEL_CSV),
                # ("Text-CSV Layout -> All Voxels",    main.on_menu_view_output, wx.ITEM_NORMAL, ViewIdsWatref.OUTPUT_ALL_VOXELS_CSV),
                ("Text-CSV Layout", (
                    ("Current Voxel", main.on_menu_view_results, wx.ITEM_NORMAL, ViewIdsWatref.RESULTS_CSV_CURRENT),
                    ("All Voxels", main.on_menu_view_results, wx.ITEM_NORMAL, ViewIdsWatref.RESULTS_CSV_ALL))),
                ("LCM Layout", (
                    ("PDF Format (Ph0/1 corr)", main.on_menu_view_results, wx.ITEM_NORMAL, ViewIdsWatref.RESULTS_LCM_PDF),
                    ("PDF Format Multi-page (Ph0/1 corr)", main.on_menu_view_results, wx.ITEM_NORMAL, ViewIdsWatref.RESULTS_LCM_PDF_MULTI),
                    ("PNG Format (Ph0/1 corr)", main.on_menu_view_results, wx.ITEM_NORMAL, ViewIdsWatref.RESULTS_LCM_PNG))),
            )
            ),
           )
    

    # There are four plot menus. They all look the same, but the menu items 
    # have unique ids.
    plot_menus = [ ]
    for i in range(4):
        plot_menus.append(
        ( 
            ("Plot Type", (
                ("Raw Data",            main.on_menu_plot_x, wx.ITEM_RADIO, PlotXIds.RAW_DATA[i]),
                ("Fitted Data",         main.on_menu_plot_x, wx.ITEM_RADIO, PlotXIds.FITTED_DATA[i]),
                ("Baseline Data",       main.on_menu_plot_x, wx.ITEM_RADIO, PlotXIds.BASELINE[i]),
                ("Raw-Fit",             main.on_menu_plot_x, wx.ITEM_RADIO, PlotXIds.COMBO_RAW_MINUS_FIT[i]),
                ("Raw-Base",            main.on_menu_plot_x, wx.ITEM_RADIO, PlotXIds.COMBO_RAW_MINUS_BASE[i]),
                ("Raw-Fit-Base",        main.on_menu_plot_x, wx.ITEM_RADIO, PlotXIds.COMBO_RAW_MINUS_FIT_MINUS_BASE[i]),
                ("Raw and Fit",         main.on_menu_plot_x, wx.ITEM_RADIO, PlotXIds.COMBO_RAW_AND_FIT[i]),
                ("Raw and Base",        main.on_menu_plot_x, wx.ITEM_RADIO, PlotXIds.COMBO_RAW_AND_BASE[i]),
                ("Raw and (Fit+Base)",  main.on_menu_plot_x, wx.ITEM_RADIO, PlotXIds.COMBO_RAW_AND_FIT_PLUS_BASE[i]),
                ("(Raw-Base) and Fit",  main.on_menu_plot_x, wx.ITEM_RADIO, PlotXIds.COMBO_RAW_MINUS_BASE_AND_FIT[i]),
                ("Fit+Base",            main.on_menu_plot_x, wx.ITEM_RADIO, PlotXIds.COMBO_FIT_PLUS_BASE[i]),
                ("Fit and Base",        main.on_menu_plot_x, wx.ITEM_RADIO, PlotXIds.COMBO_FIT_AND_BASE[i]),
                ("Raw and InitialModel",main.on_menu_plot_x, wx.ITEM_RADIO, PlotXIds.COMBO_RAW_AND_INIT_MODEL[i]),
                ("Raw and WeightArray", main.on_menu_plot_x, wx.ITEM_RADIO, PlotXIds.COMBO_RAW_AND_WT_ARR[i]),
                          )
            ),
            common_menu.SEPARATOR,
            ("Data Type",   (
                ("Real",      main.on_menu_plot_x, wx.ITEM_RADIO, PlotXIds.DATA_TYPE_REAL[i]),
                ("Imaginary", main.on_menu_plot_x, wx.ITEM_RADIO, PlotXIds.DATA_TYPE_IMAGINARY[i]),
                ("Magnitude", main.on_menu_plot_x, wx.ITEM_RADIO, PlotXIds.DATA_TYPE_MAGNITUDE[i]),
                common_menu.SEPARATOR,
                ("Summed Lines",  main.on_menu_plot_x, wx.ITEM_CHECK, PlotXIds.DATA_TYPE_SUMMED[i]),
                            )
            )
        )
        )

    help = [
#            ("&User Manual",            main.on_user_manual),
            ("&Analysis Online User Manual",   main.on_analysis_online_user_manual),
            ("&Vespa Help Online",      main.on_vespa_help_online),
            ("&About", main.on_about, wx.ITEM_NORMAL, wx.ID_ABOUT),
           ]
    
    if util_common_config.VespaConfig().show_wx_inspector:  
        help.append(common_menu.SEPARATOR)
        help.append( ("Show Inspection Tool", main.on_show_inspection_tool) )

    return file_, processing, view_start, view_raw, view_prep_fidsum, view_prep_timeseries, \
           view_spectral, view_voigt, view_giso, view_watref, plot_menus, help



# def import_from_dotted_path(dotted_names, path=None):
#     """
#     import_from_dotted_path('foo.bar') -> from foo import bar; return bar
#
#     """
#     next_module, remaining_names = dotted_names.split('.', 1)
#     fp, pathname, description = imp.find_module(next_module, path)
#     module = imp.load_module(next_module, fp, pathname, description)
#
#     if hasattr(module, remaining_names):
#         return getattr(module, remaining_names)
#     if '.' not in remaining_names:
#         return module
#     return import_from_dotted_path(remaining_names, path=module.__path__)




