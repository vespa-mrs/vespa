# Python modules


# 3rd party modules
import wx

# Our modules
import vespa.common.menu as common_menu
import vespa.common.util.config as util_common_config

from vespa.common.menu import _VespaMenu


class ViewIds(common_menu.IdContainer):
    # This one-off class is a container for the ids of all of the menu items
    # to which I need explicit references.
    ZERO_LINE_SHOW = common_menu.IdContainer.PLACEHOLDER

    XAXIS_SHOW = common_menu.IdContainer.PLACEHOLDER
    # FIXME PS - temporarily disabling PPM/kHz 
    # XAXIS_PPM = common_menu.IdContainer.PLACEHOLDER
    # XAXIS_KILOHERTZ = common_menu.IdContainer.PLACEHOLDER

    DATA_TYPE_REAL = common_menu.IdContainer.PLACEHOLDER
    DATA_TYPE_REAL_IMAGINARY = common_menu.IdContainer.PLACEHOLDER

    VIEW_TO_PNG = common_menu.IdContainer.PLACEHOLDER
    VIEW_TO_SVG = common_menu.IdContainer.PLACEHOLDER
    VIEW_TO_EPS = common_menu.IdContainer.PLACEHOLDER
    VIEW_TO_PDF = common_menu.IdContainer.PLACEHOLDER


# When main creates the menu bar it sets the variable below to that instance. 
# It's a convenience. It's the same as wx.GetApp().GetTopWindow().GetMenuBar(),
# but much easier to type.
bar = None


class PulseMenuBar(common_menu.VespaMenuBar):
    def __init__(self, main):
        common_menu.VespaMenuBar.__init__(self, main)
        
        ViewIds.init_ids()

        # _get_menu_data() is called just once, right here. 
        pulse_design, transforms, manage, view, help = _get_menu_data(main)

        # Build the top-level menus that are always present. 
        pulse_design = common_menu.create_menu(main, "&Pulse Design", pulse_design)
        transforms = common_menu.create_menu(main, "Add &Transforms", transforms)
        manage = common_menu.create_menu(main, "&Management", manage)
        view = common_menu.create_menu(main, "&View", view)
        help = common_menu.create_menu(main, "&Help", help)

        for menu in (pulse_design, transforms, manage, view, help):
            self.Append(menu, menu.label)

        ViewIds.enumerate_booleans(self.view_menu)

        # Last but not least, the self._transform_menu_items dict here to keep
        # track of the transforms that exist in the database. The dict is 
        # keyed by the id of the menu item (which we generate dynamically) and 
        # the values are 2-tuples of (transform_id, transform_name). main.py 
        # uses this dictionary to add transforms to a design notebook.
        self._transform_menu_items = { }
        self.main = main


    @property
    def transforms_menu(self):
        """A convenience property that returns the File menu"""
        return self.GetMenu(self.FindMenu("Add Transforms"))

#     @property
#     def create_menu(self):
#         """A convenience property that returns the File menu"""
#         menu = self.transforms_menu
#         return self.GetMenu(menu.FindItem("Create"))


    def get_transform_menu_item_info(self, menu_id):
        """
        Given the wx id of a menu item, returns a 2-tuple of 
        (transform_id, transform_name) 
        """
        return self._transform_menu_items[menu_id]
    
    
    def update_transforms(self, create_items, modify_items=[]):      
    
        old_items = self.transforms_menu.GetMenuItems()
        for item in old_items:
            self.transforms_menu.Remove(item.GetId())
        
        self._transform_menu_items = { }
        
        create_menu = _VespaMenu('')
        
        for item in create_items:
            id_ = wx.NewIdRef()
            menu_item = wx.MenuItem(create_menu, id_, text=item[1])
            create_menu.Append(menu_item)
            # Save the data associated with this menu item.
            self._transform_menu_items[id_] = (item[0], item[1]+'_name')
            self.main.Bind(wx.EVT_MENU, self.main.on_transform_item, menu_item)
        
        self.transforms_menu.Append(wx.ID_ANY, "Create", create_menu)
        
        for item in modify_items:
            id_ = wx.NewIdRef()
            menu_item = wx.MenuItem(self.transforms_menu, id_,text=item[1])
            self.transforms_menu.Append(menu_item)
            # Save the data associated with this menu item.
            self._transform_menu_items[id_] = (item[0], item[1]+'_name')
            self.main.Bind(wx.EVT_MENU, self.main.on_transform_item, menu_item)
    
        

def _get_menu_data(main):
    # Note that wx treats the ids wx.ID_EXIT and wx.ID_ABOUT specially by 
    # moving them to their proper location on the Mac. wx will also change
    # the text of the ID_EXIT item to "Quit" as is standard under OS X. 
    # Quit is also the standard under Gnome but unfortunately wx doesn't seem
    # to change Exit --> Quit there, so our menu looks a little funny under
    # Gnome.
    
    pulseproj = (
                    ("N&ew\tCTRL+N",  main.on_new_pulse_design),
                    ("Copy &Tab to New\tCTRL+T", main.on_copy_pulse_design_tab),                     
                    ("O&pen\tCTRL+O", main.on_open_pulse_design),
                    common_menu.SEPARATOR,
                    ("R&un\tCTRL+R", main.on_run_pulse_design),
                    ("S&ave\tCTRL+S", main.on_save_pulse_design),
                    # Ctrl+Shift+S is the preferred accelerator for Save As
                    # under OS X and GTK. Microsoft doesn't provide one at
                    # all. For consistency, we'll leave off the accelerator.
                    ("Save A&s", main.on_save_as_pulse_design),
                    ("Close\tCTRL+W", main.on_close_pulse_design),
                    common_menu.SEPARATOR,
                    ("Third Party Export...", main.on_third_party_export),
                    common_menu.SEPARATOR,
                    ("Exit",          main.on_self_close, wx.ITEM_NORMAL, wx.ID_EXIT)
                )

    manage = (
                ("Manage &Pulse Designs",      main.on_manage_pulse_designs),
                ("Manage &Transform Kernels",   main.on_manage_transform_kernels),
                ("Manage &Machine Specs Templates", main.on_manage_machine_specs_templates)
             )    

    transforms = (
                    ("Create", (
                        
                        # FIXME PS - temporarily disabling PPM/kHz 
                         common_menu.SEPARATOR,
                        # ("[mm]", main.on_menu_view_option, wx.ITEM_RADIO, ViewIds.XAXIS_PPM),
                        # ("[kHz]", main.on_menu_view_option, wx.ITEM_RADIO, ViewIds.XAXIS_KILOHERTZ)
                        )
                    ),
                  )
                  

    view = (
            ("Show Zero Line", main.on_menu_view_option, wx.ITEM_CHECK, ViewIds.ZERO_LINE_SHOW),
            ("X-Axis", (
                ("Show", main.on_menu_view_option, wx.ITEM_CHECK, ViewIds.XAXIS_SHOW),
                # FIXME PS - temporarily disabling PPM/kHz 
                # common_menu.SEPARATOR,
                # ("[mm]", main.on_menu_view_option, wx.ITEM_RADIO, ViewIds.XAXIS_PPM),
                # ("[kHz]", main.on_menu_view_option, wx.ITEM_RADIO, ViewIds.XAXIS_KILOHERTZ)
                )
            ),
            ("Data Type", (
                ("Real",            main.on_menu_view_option, wx.ITEM_RADIO, ViewIds.DATA_TYPE_REAL),
                ("Real+Imaginary",  main.on_menu_view_option, wx.ITEM_RADIO, ViewIds.DATA_TYPE_REAL_IMAGINARY))),
            common_menu.SEPARATOR,
            ("Output", (
                ("Plot to PNG", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIds.VIEW_TO_PNG),
                ("Plot to SVG", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIds.VIEW_TO_SVG),
                ("Plot to EPS", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIds.VIEW_TO_EPS),
                ("Plot to PDF", main.on_menu_view_output, wx.ITEM_NORMAL, ViewIds.VIEW_TO_PDF)),
            )
           )

    help = [
#                ("&User Manual",          main.on_user_manual),
                ("&Pulse Online User Manual",    main.on_pulse_online_user_manual),
                ("&Vespa Help Online",    main.on_vespa_help_online),
                ("&About", main.on_about, wx.ITEM_NORMAL, wx.ID_ABOUT),
           ]

    if util_common_config.VespaConfig().show_wx_inspector:
        help.append(common_menu.SEPARATOR)
        help.append( ("Show Inspection Tool", main.on_show_inspection_tool) )

    return (pulseproj, transforms, manage, view, help)


