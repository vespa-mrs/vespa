
"""\
FloatSpinMultiplier objects

@copyright: 2016-2019     Brian J. Soher
@license: MIT (see LICENSE.txt) - THIS PROGRAM COMES WITH NO WARRANTY
"""

import wx

from edit_windows import ManagedBase, EditStylesMixin
from tree import Node
import common, compat, config, misc
import new_properties as np
from collections import OrderedDict

from wx.lib.agw.floatspin import FloatSpin, EVT_FLOATSPIN, FS_LEFT, FS_RIGHT, FS_CENTRE, FS_READONLY
from .floatspin_multiplier_base import *



class EditFloatSpinMultiplier(ManagedBase, EditStylesMixin):
    """\
    Class to handle FloatSpinMultiplier objects
    
    """
    _PROPERTIES = ["Widget", "range", "value", "increment", "digits", "extrastyle", "style"]
    PROPERTIES = ManagedBase.PROPERTIES + _PROPERTIES + ManagedBase.EXTRA_PROPERTIES

    def __init__(self, name, parent, id, sizer, pos):

        # Initialise parent classes
        ManagedBase.__init__(self, name, 'FloatSpinMultiplier', parent, id, sizer, pos)
        EditStylesMixin.__init__(self)

        # initialise instance properties
        self.range      = np.FloatRangePropertyA( "0.0, 100.0" )
        self.value      = np.SpinDoublePropertyA(0.0, val_range=(0.0,100.0), immediate=True, default_value=0.0)
        self.multiplier = np.SpinDoublePropertyA(1.1, val_range=(0.0,100.0), immediate=True, default_value=1.0)
        self.digits     = np.SpinPropertyA(3, val_range=(0,20), immediate=True, default_value=3)
        self.extrastyle = WidgetExtraStyleProperty()

        if config.preferences.default_border:
            self.border = config.preferences.default_border_size
            self.flag = wx.ALL


    def create_widget(self):
        mi,ma  = self.properties["range"].get_tuple()
        value  = self.properties["value"].get()
        mult   = self.properties["multiplier"].get()
        digits = self.properties["digits"].get()
        
        self.widget = FloatSpinMultiplier(self.parent.widget, self.id, min_val=mi, max_val=ma,
                                          value=value, multiplier=mult, digits=digits)


    def finish_widget_creation(self, sel_marker_parent=None, re_add=True):
        ManagedBase.finish_widget_creation(self, sel_marker_parent, re_add)
        self.widget.Bind(wx.EVT_CHILD_FOCUS, self._on_set_focus)
        self.widget.Bind(wx.EVT_SET_FOCUS, self._on_set_focus)
        self.widget.Bind(wx.EVT_SPIN, self.on_set_focus)


    def _on_set_focus(self, event):
        # don't set focused_widget during event, as this may cause crashes
        if not misc.focused_widget is self:
            misc.set_focused_widget(self, delayed=True)
        event.Skip()


    def properties_changed(self, modified):  # from SpinCtrlDouble

        if not modified or "range" in modified and self.widget:
            mi,ma = self.properties["range"].get_tuple()
            self.widget.SetRange(mi, ma)
            self.properties["value"].set_range(mi,ma)
            self.properties["multiplier"].set_range(mi,ma)

        if not modified or "multiplier" in modified and self.widget:
            self.widget.SetIncrement(self.multiplier)

        if not modified or "digits" in modified and self.widget:
            self.widget.SetDigits(self.digits)

        if not modified or "value" in modified or "range" in modified:
            # check that value is inside range
            value_p = self.properties["value"]
            if value_p.is_active():
                mi,ma = self.properties["range"].get_tuple()
                value = value_p.get()
                if value<mi:
                    value_p.set(mi)
                    value = mi
                elif value>ma:
                    value_p.set(ma)
                    value = ma
                if self.widget:
                    self.widget.SetValue(value)

        # FS style changes do not need any action here
        
        EditStylesMixin.properties_changed(self, modified)
        ManagedBase.properties_changed(self, modified)

# end of class EditFloatSpinMultiplier



class WidgetExtraStyleProperty(np.WidgetStyleProperty):
    # for Extra-Styles widget/settings for FloatSpin; e.g. FS_LEFT, FS_RIGHT, FS_CENTER,FS_READONLY

    FLAG_DESCRIPTION = OrderedDict()
    FLAG_DESCRIPTION['Extra-Styles'] = ['FS_LEFT', 'FS_CENTRE', 'FS_RIGHT', 'FS_READONLY']
    FLAG_NAMES  = sum( list(FLAG_DESCRIPTION.values()), [] )

    EXCLUDES = {'FS_LEFT'   : set(['FS_CENTRE', 'FS_RIGHT']),
                'FS_CENTRE' : set(['FS_LEFT',   'FS_RIGHT']),
                'FS_RIGHT'  : set(['FS_LEFT',   'FS_CENTRE']) }

    STYLE_DEFS = {
        'FS_LEFT':     { 'desc': _('Same as wxTE_LEFT for wxTextCtrl: the numerical text is left aligned.'),
                         'exclude': 'FS_CENTRE|FS_RIGHT'  },
        'FS_CENTRE':   { 'desc': _('Same as wxTE_CENTER for wxTextCtrl: the numerical text is center aligned.'),
                         'exclude': 'FS_LEFT|FS_RIGHT'  },
        'FS_RIGHT':    { 'desc': _('Same as wxTE_RIGHT for wxTextCtrl: the numerical text is right aligned.'),
                         'exclude': 'FS_LEFT|FS_CENTRE'  },
        'FS_READONLY': { 'desc': _('The numericalS text is read only.') },
        }

    def _ensure_values(self):
        """ from the base classs - but needs some rethink for this way of using it """
        self._values = []  # the associated flag values
        widget_writer = self.owner.widget_writer

        for name in self._names:
            wx_name = widget_writer.cn_f(name)
            if not wx_name:  # cn_f() returns an empty string if the given styles are not supported
                self._values.append(None)
                continue
            try:
                self._values.append( self.owner.wxname2attr(wx_name) )  # FIXME bjs - this always fails
            except:
                self._values.append(None)                               # but this is not objectionable as a fallback
                continue

    def set_owner(self, owner, attname):
        """ style information is taken from class constants above """
        np._CheckListProperty.set_owner(self, owner, attname) 
        self.style_defs = self.STYLE_DEFS
        self.styles     = self.FLAG_DESCRIPTION
        self._names     = sum( list(self.styles.values()), [] )
        self._values    = None
        self.set('FS_LEFT')
        self.default_value = set(set(['FS_LEFT']))
        self.set("")
        self.modified = False

    def create_editor(self, panel, sizer):
        """ similar to base class but simpler """
        self._choices = [] # the checkboxes
        self._ensure_values()

        tooltips = self._create_tooltip_text()

        static_box = wx.StaticBox(panel, -1, _("Extra-Style"), style=wx.FULL_REPAINT_ON_RESIZE)
        box_sizer = wx.StaticBoxSizer(static_box, wx.VERTICAL)
        for name, flag_value in zip(self._names, self._values):
            if name in self.STYLE_DEFS:
                style_def = self.STYLE_DEFS[name]
            else:
                style_def = {}
            if "obsolete" in style_def or "rename_to" in style_def:
                self._choices.append(None)
                continue
            checkbox = wx.CheckBox(panel, -1, name)

            if name in tooltips:
                compat.SetToolTip( checkbox, tooltips[name] )

            self._choices.append(checkbox)
            box_sizer.Add(checkbox)

        sizer.Add(box_sizer, 0, wx.ALL | wx.EXPAND, 5)

        self.update_display(True)
        for checkbox in self._choices:
            if checkbox is not None:
                checkbox.Bind(wx.EVT_CHECKBOX, self.on_checkbox)     


def builder(parent, sizer, pos, number=[1]):
    """\
    factory function for EditFloatSpin objects.
    """
    name = 'float_spin_multiplier_%d' % number[0]
    while common.app_tree.has_name(name):
        number[0] += 1
        name = 'float_spin_multiplier_%d' % number[0]
        
    with parent.frozen():
        spin = EditFloatSpinMultiplier(name, parent, wx.ID_ANY, sizer, pos)
        spin.properties["style"].set_to_default()
        spin.check_defaults()
        node = Node(spin)
        spin.node = node
        if parent.widget: spin.create()
    common.app_tree.insert(node, sizer.node, pos - 1)


def xml_builder(attrs, parent, sizer, sizeritem, pos=None):
    """ factory function to build EditFloatSpinMultiplier objects from a XML file """
    from xml_parse import XmlParsingError
    try:
        name = attrs['name']
    except KeyError:
        raise XmlParsingError(_("'name' attribute missing"))
    if sizer is None or sizeritem is None:
        raise XmlParsingError(_("sizer or sizeritem object cannot be None"))

    spin = EditFloatSpinMultiplier( name, parent, wx.ID_ANY, sizer, pos )
    spin.properties["value"].set_active(False)
    node = Node(spin)
    spin.node = node
    if pos is None:
        common.app_tree.add(node, sizer.node)
    else:
        common.app_tree.insert(node, sizer.node, pos-1)
    return spin


def initialize():
    """\
    initialization function for the module: returns a wxBitmapButton to be
    added to the main palette.
    """
    import os
    common.widgets['EditFloatSpinMultiplier'] = builder
    common.widgets_from_xml['EditFloatSpinMultiplier'] = xml_builder

    return common.make_object_button('EditFloatSpinMultiplier', os.path.join(os.path.dirname(__file__), 'floatspinmultiplier.xpm'))


