# Python modules

import abc

# 3rd party modules
import wx

# Our modules
import vespa.common.util.config as util_common_config


##########################################################################
# This is a collection of menu-related constants, functions and utilities. 
# The function that builds the menu bar lives here, as does the menu 
# definition.
##########################################################################


# SEPARATOR defines a menu separator for use in append_items().
SEPARATOR = (None, )

class _VespaMenu(wx.Menu):
    """A menu item that knows its own label (standard wx.Menus don't)."""
    def __init__(self, label=""):
        wx.Menu.__init__(self)
        
        self.label = label


def append_items(main, menu, items):
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
        if item == SEPARATOR:
            menu.AppendSeparator()
        else:
            if callable(item[1]):
                # This is a regular menu item because it has an event 
                # handler (method) as its second element.
                create_menu_item(main, menu, *item)
            else:
                # This is a menu item with subitems
                label, subitems = item
                submenu = create_menu(main, "", subitems)
                menu.Append(wx.ID_ANY, label, submenu)    



def create_menu(main, label, items):
    """Given a reference to the main event handler for an app (the one that
    handles all of the menu events), a label (e.g. '&View', and a dict of 
    items, creates a menu, populates it with the items and returns it.

    See append_items() for details on the dict format.
    """
    menu = _VespaMenu(label)

    append_items(main, menu, items)

    return menu


def create_menu_item(main, menu, label, handler=None, kind=wx.ITEM_NORMAL,
                      id_=wx.ID_ANY):
    """Given a reference to the main frame, a menu instance, and an item 
    label, adds a new menu item to the menu with that label text. 

    If an event handler (method) is passed, the menu item is bound to it.

    The param kind can be wx.ITEM_CHECK, wx.ITEM_RADIO, etc.

    The param id_ can be a wx special ID like wx.ID_EXIT or wx.ID_ABOUT.     
    """
    menu_item = menu.Append(id_, label, "", kind)
    if handler:
        main.Bind(wx.EVT_MENU, handler, menu_item)


class IdContainer(object, metaclass=abc.ABCMeta):
    """IdContainer is an abstract base class for the menu item constant classes
    It exists to contain the method definitions. Subclasses will add on a 
    bunch of app-specific constants that are used for menu ids.

    These ids allow us to handle menu item events and track preferences through
    a lot of meta-code; i.e. code that can figure out on the fly what item was
    selected and the preference to which it relates. 
 
    Usually we assign values to the constants in the class definition. However, 
    these constants represent menu items and must get their values from 
    wx.NewIdRef().

    That complicates things. This code is part of a class definition, 
    so it is executed as soon as the file is imported. At that point, the 
    wx.App object has not yet been created, so calling most wx functions 
    raises an error. So code that uses this base class must initialize its
    menu constants to PLACEHOLDER. They'll be replaced when init_ids() is
    called in the subclass __init__().

    The method enumerate_booleans() should also be called in subclass 
    __init__(). Like init_ids(), this method is only called once. It figures out
    which class attribute names (e.g. 'ZERO_LINE_TOP') are constants, and 
    which of those relate to boolean menu items. Again, due to the fact that 
    class init happens before wx is initialized means that we have to create
    the list of booleans via an explicit call in menu bar init.    
    """

    PLACEHOLDER = "replace me"


    # BOOLEANS is a tuple containing the name of this class's constant 
    # attributes that relate to booleans. e.g. 'ZERO_LINE_TOP', 
    # 'ZERO_LINE_MIDDLE', etc., but *not* 'VIEW_TO_PNG' (since it's not 
    # boolean).
    BOOLEANS = tuple()

    @classmethod
    def enumerate_booleans(klass, menu):
        """Given a menu that contains the items with the ids defined in this
        constant class, populates the BOOLEANS attribute on this class. This 
        only needs to be called once since the booleans don't change once 
        they're enumerated.
        """
        klass.BOOLEANS = [ ]

        # as of wxPython 3.0.0 (?) wx ids have their own class wx.WindowIDRef
        # which has 'value' and 'Id' attributes, so we can scan the whole 
        # util_menu object to see what items it has that are in ALL_CAPS

        for attribute_name, value in klass.__dict__.items():
            # Make sure this is one of my attributes.
            if attribute_name.isupper() and isinstance(value, wx.WindowIDRef):
                # OK, it's one of mine. See if it is boolean.
                item = menu.FindItemById(value)
                if not item:
                    raise ValueError("Unused menu constant: '%s'" % attribute_name)
                if item.IsCheckable():
                    klass.BOOLEANS.append(attribute_name)

        # I make them a tuple just to emphasize that they're static.
        klass.BOOLEANS = tuple(klass.BOOLEANS)


    @classmethod
    def init_ids(klass):
        """A single-use method that assigns values to the constants. See the 
        lengthy comment at the top of the class as to why this exists.
        """
        # Rather than hardcoding the constant names here, I use introspection
        # to enumerate them. I rely on the fact that all of our constants
        # are UPPER_CASE and are initialized to PLACEHOLDER.
        for key, value in klass.__dict__.items():
            if key.isupper() and (value == klass.PLACEHOLDER):
                # It's one of ours
                value = wx.NewIdRef()

                setattr(klass, key, value)
            #else:
                # It's not something we created; don't touch it!


class VespaMenuBar(wx.MenuBar):
    """A subclass of wx.MenuBar that adds some Vespa-specific goodies."""
    def __init__(self, main):
        wx.MenuBar.__init__(self)


    @property
    def view_menu(self):
        """A convenience property that returns the View menu"""
        return self.GetMenu(self.FindMenu("View"))


    @property
    def top_level_menus(self):
        """Returns a list of the top-level menus"""
        return [self.GetMenu(i) for i in range(self.GetMenuCount())]


    def is_checked(self, id_): 
        """True if the item identified by id_ is checked, False otherwise.""" 
        return self.FindItemById(id_).IsChecked()


    def find_top_level_parent(self, menu_item_id):
        """Given the id of a menu item, returns the top level menu that 
        contains the item, or None if no such menu is found.

        There are two caveats. First, the menu_item_id must be associated 
        with a boolean menu item (checkable or radio style). Second, the 
        parent menu must of one of the *currently visible* top level menus.
        """
        for top_level_menu in self.top_level_menus:
            state = self.get_top_level_menu_state(top_level_menu)
            if menu_item_id in state:
                return top_level_menu
        return None


    def get_top_level_menu_state(self, top_level_menu):
        """Given a top-level menu, returns a dict where the keys are the ids 
        of all of the boolean menu items on that menu. The values of the dict 
        reflect the True/False state of the menu item.

        If the menu has submenus, the menu items on the submenus are included
        in the dict, but the dict is a flattened representation of the menu.
        """
        state = { }
        
        for item in top_level_menu.GetMenuItems():
            if item.IsCheckable():
                state[item.Id] = item.IsChecked()
            elif item.IsSubMenu():
                # Traverse this submenu's items
                for subitem in item.GetSubMenu().GetMenuItems():
                    if subitem.IsCheckable():
                        state[subitem.Id] = subitem.IsChecked()

        return state


    def set_menu_from_state(self, state):
        """Given a dict created by get_top_level_menu_state(), checks and 
        unchecks menu items as appropriate.
        """
        for key, value in state.items():
            # Under normal circumstances we should find a menu item for each 
            # key. However, assuming that will *always* be true makes it 
            # difficult to temporarily remove (comment out) a menu item. So 
            # here we just quietly ignore references to menu items that don't 
            # exist.
            menu_item = self.FindItemById(key)
            if menu_item:
                menu_item.Check(value)



    # ===================    Internal Use Only    ==========================

    def _print_menu(self, identifier):
        # Code for debugging. Given an int, string or wx.Menu instance,
        # prints the labels of all the items on the menu.
        # When identifier is an int, it's treated as an index into the list
        # of top-level menus. 
        # When identifier is a string, this code finds and prints the menu 
        # with that name.
        # When idenfitier is an wx.Menu instance, it prints that instance.
        
        if isinstance(identifier, wx.Menu):
            menu = identifier
        else:
            try:
                int(identifier)
            except ValueError:
                # it's a string
                identifier = self.FindMenu(identifier)

            menu = self.GetMenu(identifier)
        
        for item in menu.GetMenuItems():
            if item.IsSubMenu():
                # Traverse this submenu's items
                for subitem in item.GetSubMenu().GetMenuItems():
                    print("%d: %s:" % (subitem.Id, subitem.GetItemLabelText()))


