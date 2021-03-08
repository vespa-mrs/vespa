# Python modules


# 3rd party modules
import wx

# Our modules
import vespa.common.menu as common_menu
import vespa.common.util.config as util_common_config


########################################################################
# This is a collection of menu-related constants, functions and utilities. 
# The function that builds the menu bar lives here, as does the menu 
# definition.
########################################################################


class ViewIds(common_menu.IdContainer):
    """A container for the ids of all of the menu items to which we need 
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

    PLOT_VIEW_FINAL = "replace me"
    PLOT_VIEW_ALL = "replace me"

    EXPERIMENT_TO_TEXT = "replace me"
    
    VIEW_TO_PNG = "replace me"
    VIEW_TO_SVG = "replace me"
    VIEW_TO_PDF = "replace me"
    VIEW_TO_EPS = "replace me"


# When main creates an instance of DatasimMenuBar(), it sets the variable
# below to that instance. It's a convenience. It's the same as 
# wx.GetApp().GetTopWindow().GetMenuBar(), but much easier to type.
bar = None


class DatasimMenuBar(common_menu.VespaMenuBar):
    """A subclass of wx.MenuBar that adds some app-specific functions
    and constants.
    
    There should be only one instance of this class per invocation of the 
    app. It's a singleton class.
    """
    
    def __init__(self, main):
        common_menu.VespaMenuBar.__init__(self, main)
        
        ViewIds.init_ids()

        # _get_menu_data() is called just once, right here. 
        datasim, view, help = _get_menu_data(main)

        # Build the top-level menus that are always present. 
        datasim = common_menu.create_menu(main, "Datasim", datasim)
        view = common_menu.create_menu(main, "&View", view)
        help = common_menu.create_menu(main, "&Help", help)

        for menu in (datasim, view, help):
            self.Append(menu, menu.label)

        ViewIds.enumerate_booleans(self.view_menu)


# ================    Module Internal Use Only    =======================


def _get_menu_data(main):
    # Note that wx treats the ids wx.ID_EXIT and wx.ID_ABOUT specially by 
    # moving them to their proper location on the Mac. wx will also change
    # the text of the ID_EXIT item to "Quit" as is standard under OS X. 
    # Quit is also the standard under Gnome but unfortunately wx doesn't seem
    # to change Exit --> Quit there, so our menu looks a little funny under
    # Gnome.

    prior = (
                ("N&ew Datasim from Experiment...\tCTRL+N",    main.on_new),
                ("O&pen Datasim...\tCTRL+O",   main.on_open),
                common_menu.SEPARATOR,
                ("S&ave\tCTRL+S",      main.on_save_viff),
                ("S&ave As...",        main.on_save_as_viff),
                common_menu.SEPARATOR,
                ("Close\tCTRL+W",       main.on_close_datasim),
                common_menu.SEPARATOR,
                ("E&xport Spectrum", (
                    ("to VIFF Raw Data...",   main.on_export_spectrum_viff),
                    ("to Siemens *.rda...",   main.on_export_spectrum_siemens_rda))),
                ("E&xport Monte Carlo", (
                    ("to VIFF Raw Data...",   main.on_export_monte_carlo_viff), )),
                common_menu.SEPARATOR,
                ("&Exit",                   main.on_self_close))
    view = (
                ("Zero Line", (
                    ("Show",   main.on_menu_view_option, wx.ITEM_CHECK, ViewIds.ZERO_LINE_SHOW),
                    common_menu.SEPARATOR,
                    ("Top",    main.on_menu_view_option, wx.ITEM_RADIO, ViewIds.ZERO_LINE_TOP),
                    ("Middle", main.on_menu_view_option, wx.ITEM_RADIO, ViewIds.ZERO_LINE_MIDDLE),
                    ("Bottom", main.on_menu_view_option, wx.ITEM_RADIO, ViewIds.ZERO_LINE_BOTTOM))),
                ("X-Axis", (
                    ("Show",   main.on_menu_view_option, wx.ITEM_CHECK, ViewIds.XAXIS_SHOW),
                    common_menu.SEPARATOR,
                    ("PPM",    main.on_menu_view_option,  wx.ITEM_RADIO, ViewIds.XAXIS_PPM),
                    ("Hertz",  main.on_menu_view_option,  wx.ITEM_RADIO, ViewIds.XAXIS_HERTZ))),
                common_menu.SEPARATOR,
                ("Data Type", (
                    ("Real",      main.on_menu_view_option, wx.ITEM_RADIO, ViewIds.DATA_TYPE_REAL),
                    ("Imaginary", main.on_menu_view_option, wx.ITEM_RADIO, ViewIds.DATA_TYPE_IMAGINARY),
                    ("Magnitude", main.on_menu_view_option, wx.ITEM_RADIO, ViewIds.DATA_TYPE_MAGNITUDE),
                    common_menu.SEPARATOR,
                    ("Summed",    main.on_menu_view_option, wx.ITEM_CHECK, ViewIds.DATA_TYPE_SUMMED))),
                common_menu.SEPARATOR,
                ("Plot Views", (
                    ("Final Only", main.on_menu_view_option, wx.ITEM_RADIO, ViewIds.PLOT_VIEW_FINAL),
                    ("All Three",  main.on_menu_view_option, wx.ITEM_RADIO, ViewIds.PLOT_VIEW_ALL))),
                common_menu.SEPARATOR,
                ("Output Experiment Text", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIds.EXPERIMENT_TO_TEXT),
                ("Output Plots", (
                    ("View to PNG", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIds.VIEW_TO_PNG),
                    ("View to SVG", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIds.VIEW_TO_SVG),
                    ("View to EPS", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIds.VIEW_TO_EPS),
                    ("View to PDF", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIds.VIEW_TO_PDF)  )))
    help = (
                ("&User Manual",          main.on_user_manual),
                ("&DataSim Help Online",  main.on_datasim_help_online),
                ("&Vespa Help Online",    main.on_vespa_help_online),
                ("&About", main.on_about, wx.ITEM_NORMAL, wx.ID_ABOUT),
           )

    if util_common_config.VespaConfig().show_wx_inspector:
        help = list(help)
        help.append(common_menu.SEPARATOR)
        help.append( ("Show Inspection Tool", main.on_show_inspection_tool) )
    
    return (prior, view, help)          

