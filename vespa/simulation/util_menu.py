# Python modules


# 3rd party modules
import wx

# Our modules
import vespa.common.menu as common_menu
import vespa.common.util.config as util_common_config


class ViewIds(common_menu.IdContainer):
    # This one-off class is a container for the ids of all of the menu items
    # to which I need explicit references.
    ZERO_LINE_SHOW = common_menu.IdContainer.PLACEHOLDER

    XAXIS_SHOW = common_menu.IdContainer.PLACEHOLDER
    XAXIS_PPM = common_menu.IdContainer.PLACEHOLDER
    XAXIS_HERTZ = common_menu.IdContainer.PLACEHOLDER

    INTEGRAL_PLOT_SHOW = common_menu.IdContainer.PLACEHOLDER
    INTEGRAL_XAXIS_SHOW = common_menu.IdContainer.PLACEHOLDER
    INTEGRAL_YAXIS_SHOW = common_menu.IdContainer.PLACEHOLDER
    CONTOUR_PLOT_SHOW = common_menu.IdContainer.PLACEHOLDER
    CONTOUR_AXES_SHOW = common_menu.IdContainer.PLACEHOLDER

    DATA_TYPE_REAL = common_menu.IdContainer.PLACEHOLDER
    DATA_TYPE_IMAGINARY = common_menu.IdContainer.PLACEHOLDER
    DATA_TYPE_MAGNITUDE = common_menu.IdContainer.PLACEHOLDER

    LINE_SHAPE_GAUSSIAN = common_menu.IdContainer.PLACEHOLDER
    LINE_SHAPE_LORENTZIAN = common_menu.IdContainer.PLACEHOLDER

    OUTPUT_1D_STACK_TO_PNG = common_menu.IdContainer.PLACEHOLDER
    OUTPUT_1D_STACK_TO_SVG = common_menu.IdContainer.PLACEHOLDER
    OUTPUT_1D_STACK_TO_EPS = common_menu.IdContainer.PLACEHOLDER
    OUTPUT_1D_STACK_TO_PDF = common_menu.IdContainer.PLACEHOLDER
    OUTPUT_INTEGRAL_TO_PNG = common_menu.IdContainer.PLACEHOLDER
    OUTPUT_INTEGRAL_TO_SVG = common_menu.IdContainer.PLACEHOLDER
    OUTPUT_INTEGRAL_TO_EPS = common_menu.IdContainer.PLACEHOLDER
    OUTPUT_INTEGRAL_TO_PDF = common_menu.IdContainer.PLACEHOLDER
    OUTPUT_CONTOUR_TO_PNG = common_menu.IdContainer.PLACEHOLDER
    OUTPUT_CONTOUR_TO_SVG = common_menu.IdContainer.PLACEHOLDER
    OUTPUT_CONTOUR_TO_EPS = common_menu.IdContainer.PLACEHOLDER
    OUTPUT_CONTOUR_TO_PDF = common_menu.IdContainer.PLACEHOLDER
    OUTPUT_TEXT_RESULTS = common_menu.IdContainer.PLACEHOLDER


# When main creates the menu bar it sets the variable below to that instance. 
# It's a convenience. It's the same as wx.GetApp().GetTopWindow().GetMenuBar(),
# but much easier to type.
bar = None


class SimulationMenuBar(common_menu.VespaMenuBar):
    """A subclass of wx.MenuBar that adds some app-specific functions
    and constants.

    Simulation only needs one instance of this class per invocation. It's a 
    singleton class.
    """
    def __init__(self, main):
        common_menu.VespaMenuBar.__init__(self, main)
        
        ViewIds.init_ids()

        # _get_menu_data() is called just once, right here. 
        experiment, manage, view, help = _get_menu_data(main)

        # Build the top-level menus that are always present. 
        experiment = common_menu.create_menu(main, "&Experiment", experiment)
        manage = common_menu.create_menu(main, "&Management", manage)
        view = common_menu.create_menu(main, "&View", view)
        help = common_menu.create_menu(main, "&Help", help)

        for menu in (experiment, manage, view, help):
            self.Append(menu, menu.label)

        ViewIds.enumerate_booleans(self.view_menu)



def _get_menu_data(main):
    # Note that wx treats the ids wx.ID_EXIT and wx.ID_ABOUT specially by 
    # moving them to their proper location on the Mac. wx will also change
    # the text of the ID_EXIT item to "Quit" as is standard under OS X. 
    # Quit is also the standard under Gnome but unfortunately wx doesn't seem
    # to change Exit --> Quit there, so our menu looks a little funny under
    # Gnome.
    experiment = (
        ("N&ew\tCTRL+N",  main.on_new_experiment),
        ("O&pen\tCTRL+O", main.on_open_experiment),
        common_menu.SEPARATOR,
#        ("Copy &Tab to New\tCTRL+T", main.on_copy_tab),
        ("Derivations...", (
            ("Copy &Tab to New\tCTRL+T",     main.on_copy_tab),
            ("Calculate Add/Sub Tab in New", main.on_add_sub_tab))),
        common_menu.SEPARATOR,
        ("S&ave\tCTRL+S", main.on_save_experiment),
        ("Close\tCTRL+W", main.on_close_experiment),
        common_menu.SEPARATOR,
        ("Third Party Export To", (
            ("Analysis Prior XML...",   main.on_menu_third_analysis_prior),
            ("MIDAS Prior XML...",      main.on_menu_third_midas_prior),
            ("LCModel...",              main.on_menu_third_lcmodel),
            ("jMRUI Text...",           main.on_menu_third_jmrui_text),
            ("GAVA Text...",            main.on_menu_third_gava_text))),
        common_menu.SEPARATOR,
        ("Exit", main.on_self_close, wx.ITEM_NORMAL, wx.ID_EXIT),
     )

    manage = (
                ("Manage &Experiments",      main.on_manage_experiments),
                ("Manage &Metabolites",      main.on_manage_metabolites),
                ("Manage &Pulse Sequences",  main.on_manage_pulse_sequences)
             )

    view = (
        ("Show Zero Line", main.on_menu_view_option, wx.ITEM_CHECK, ViewIds.ZERO_LINE_SHOW),
        ("X-Axis", (
            ("Show", main.on_menu_view_option, wx.ITEM_CHECK, ViewIds.XAXIS_SHOW),
            common_menu.SEPARATOR,
            ("PPM", main.on_menu_view_option, wx.ITEM_RADIO, ViewIds.XAXIS_PPM),
            ("Hertz", main.on_menu_view_option, wx.ITEM_RADIO, ViewIds.XAXIS_HERTZ))),
        common_menu.SEPARATOR,
        ("Data Type", (
            ("Real",      main.on_menu_view_option, wx.ITEM_RADIO, ViewIds.DATA_TYPE_REAL),
            ("Imaginary", main.on_menu_view_option, wx.ITEM_RADIO, ViewIds.DATA_TYPE_IMAGINARY),
            ("Magnitude", main.on_menu_view_option, wx.ITEM_RADIO, ViewIds.DATA_TYPE_MAGNITUDE))),
        ("Line Shape", (
            ("Gaussian",   main.on_menu_view_option, wx.ITEM_RADIO, ViewIds.LINE_SHAPE_GAUSSIAN),
            ("Lorentzian", main.on_menu_view_option, wx.ITEM_RADIO, ViewIds.LINE_SHAPE_LORENTZIAN))),
        common_menu.SEPARATOR,
        ("Contour Plot", (
            ("Show", main.on_menu_view_option, wx.ITEM_CHECK, ViewIds.CONTOUR_PLOT_SHOW),
            ("Show Axes", main.on_menu_view_option, wx.ITEM_CHECK, ViewIds.CONTOUR_AXES_SHOW),
                         )
        ),
        common_menu.SEPARATOR,
        ("Integral Plot", (
            ("Show", main.on_menu_view_option, wx.ITEM_CHECK, ViewIds.INTEGRAL_PLOT_SHOW),
            ("Show X-axis", main.on_menu_view_option, wx.ITEM_CHECK, ViewIds.INTEGRAL_XAXIS_SHOW),
            ("Show Y-axis", main.on_menu_view_option, wx.ITEM_CHECK, ViewIds.INTEGRAL_YAXIS_SHOW),
                         )
        ),
        common_menu.SEPARATOR,
        ("Output", (
            ("1D/Stackplot", (
                ("Plot to PNG", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIds.OUTPUT_1D_STACK_TO_PNG),
                ("Plot to SVG", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIds.OUTPUT_1D_STACK_TO_SVG),
                ("Plot to EPS", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIds.OUTPUT_1D_STACK_TO_EPS),
                ("Plot to PDF", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIds.OUTPUT_1D_STACK_TO_PDF))),
            ("Integral Plot", (
                ("Plot to PNG", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIds.OUTPUT_INTEGRAL_TO_PNG),
                ("Plot to SVG", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIds.OUTPUT_INTEGRAL_TO_SVG),
                ("Plot to EPS", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIds.OUTPUT_INTEGRAL_TO_EPS),
                ("Plot to PDF", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIds.OUTPUT_INTEGRAL_TO_PDF))),
            ("Contour Plot", (
                ("Plot to PNG", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIds.OUTPUT_CONTOUR_TO_PNG),
                ("Plot to SVG", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIds.OUTPUT_CONTOUR_TO_SVG),
                ("Plot to EPS", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIds.OUTPUT_CONTOUR_TO_EPS),
                ("Plot to PDF", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIds.OUTPUT_CONTOUR_TO_PDF))),
            ("Text Results",    main.on_menu_view_output, wx.ITEM_NORMAL, ViewIds.OUTPUT_TEXT_RESULTS)
                    )
        )
    )
    
    help = [
                ("&User Manual",             main.on_user_manual),
                ("&Simulation Help Online",  main.on_simulation_help_online),
                ("&Vespa Help Online",       main.on_vespa_help_online),
                ("&About", main.on_about, wx.ITEM_NORMAL, wx.ID_ABOUT),
           ]

    if util_common_config.VespaConfig().show_wx_inspector:
        help.append(common_menu.SEPARATOR)
        help.append( ("Show Inspection Tool", main.on_show_inspection_tool) )

    return (experiment, manage, view, help)


